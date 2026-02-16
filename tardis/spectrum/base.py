import warnings

import numpy as np
from astropy import units as u

from tardis.io.hdf_writer_mixin import HDFWriterMixin
from tardis.spectrum.formal_integral.base import IntegrationError
from tardis.spectrum.spectrum import TARDISSpectrum
from tardis.spectrum.virtual_packet_solver import VirtualPacketSolver
from tardis.util.base import (
    quantity_linspace,
)


class SpectrumSolver(HDFWriterMixin):
    hdf_properties = [
        "montecarlo_virtual_luminosity",
        "spectrum_real_packets",
        "spectrum_virtual_packets",
        "spectrum_real_packets_reabsorbed",
        "spectrum_integrated",
    ]

    hdf_name = "spectrum"

    def __init__(
        self, transport_state, spectrum_frequency_grid, integrator_settings
    ):
        self.transport_state = transport_state
        self.spectrum_frequency_grid = spectrum_frequency_grid
        self._montecarlo_virtual_luminosity = u.Quantity(
            np.zeros_like(self.spectrum_frequency_grid.value), "erg / s"
        )  # should be init with v_packets_energy_hist
        self._integrator = None
        self.integrator_settings = integrator_settings
        self._spectrum_integrated = None
        self.virtual_packet_state = None
        self.virtual_packet_solver = None

    def setup_optional_spectra(
        self,
        transport_state,
        virtual_packet_luminosity=None,
        integrator=None,
        simulation_state=None,
        transport=None,
        plasma=None,
        opacity_state=None,
        macro_atom_state=None,
    ):
        """Set up the solver to handle real and virtual spectra

        Parameters
        ----------
        transport_state : MonteCarloTransportState
            The transport state to be used to compute the spectra
        v_packets_energy_hist : np.ndarray
            Virtual packets energy histogram, unnormalized
        integrator : FormalIntegratorSolver, optional
            Integrator to compute the integrated spectrum with
        """
        self.transport_state = transport_state
        if virtual_packet_luminosity is not None:
            self._montecarlo_virtual_luminosity.value[:] = (
                virtual_packet_luminosity
            )
        self._integrator = integrator
        self.simulation_state = simulation_state
        self.opacity_state = opacity_state
        self.transport = transport
        self.plasma = plasma
        self.macro_atom_state = macro_atom_state

    @property
    def spectrum_real_packets(self):
        return TARDISSpectrum(
            self.spectrum_frequency_grid, self.montecarlo_emitted_luminosity
        )

    @property
    def spectrum_real_packets_reabsorbed(self):
        return TARDISSpectrum(
            self.spectrum_frequency_grid, self.montecarlo_reabsorbed_luminosity
        )

    @property
    def spectrum_virtual_packets(self):
        # Use VirtualPacketState if available (postprocessing), otherwise fall back to old method
        if self.virtual_packet_state is not None:
            luminosity = (
                u.Quantity(self.virtual_packet_state.luminosity, "erg / s")[:-1]
                / self.transport_state.time_of_simulation.value
            )
            return TARDISSpectrum(self.spectrum_frequency_grid, luminosity)

        if np.all(self.montecarlo_virtual_luminosity == 0):
            warnings.warn(
                "SpectrumSolver.spectrum_virtual_packets "
                "is zero. Virtual packets are now generated via postprocessing. "
                "Please call sim.generate_virtual_spectrum() after running the simulation "
                "with no_of_virtual_packets > 0",
                UserWarning,
            )

        return TARDISSpectrum(
            self.spectrum_frequency_grid, self.montecarlo_virtual_luminosity
        )

    @property
    def spectrum_integrated(self):
        if self._spectrum_integrated is None and self.integrator is not None:
            # This was changed from unpacking to specific attributes as compute
            # is not used in calculate_spectrum
            try:
                self._spectrum_integrated = self.integrator.solve(
                    self.spectrum_frequency_grid[:-1],
                    self.simulation_state,
                    self.transport,
                    self.opacity_state,
                    self.plasma.atomic_data,
                    self.plasma.electron_densities,
                    self.macro_atom_state,
                )
            except IntegrationError:
                # if integration is impossible or fails, return an empty spectrum
                warnings.warn(
                    "The FormalIntegrator is not yet implemented for the full "
                    "relativity mode or continuum processes. "
                    "Please run with config option enable_full_relativity: "
                    "False and continuum_processes_enabled: False "
                    "This RETURNS AN EMPTY SPECTRUM!",
                    UserWarning,
                )
                self._spectrum_integrated = TARDISSpectrum(
                    np.array([np.nan, np.nan]) * u.Hz,
                    np.array([np.nan]) * u.erg / u.s,
                )
        return self._spectrum_integrated

    @property
    def integrator(self):
        if self._integrator is None:
            warnings.warn(
                "SpectrumSolver.integrator: "
                "The FormalIntegrator is not yet available."
                "Please run the montecarlo simulation at least once.",
                UserWarning,
            )
        if self.transport_state.enable_full_relativity:
            raise NotImplementedError(
                "The FormalIntegrator is not yet implemented for the full "
                "relativity mode. "
                "Please run with config option enable_full_relativity: "
                "False."
            )
        return self._integrator

    @property
    def montecarlo_reabsorbed_luminosity(self):
        return u.Quantity(
            np.histogram(
                self.transport_state.reabsorbed_packet_nu,
                weights=self.transport_state.reabsorbed_packet_luminosity,
                bins=self.spectrum_frequency_grid,
            )[0],
            "erg / s",
        )

    @property
    def montecarlo_emitted_luminosity(self):
        return u.Quantity(
            np.histogram(
                self.transport_state.emitted_packet_nu,
                weights=self.transport_state.emitted_packet_luminosity,
                bins=self.spectrum_frequency_grid,
            )[0],
            "erg / s",
        )

    @property
    def montecarlo_virtual_luminosity(self):
        return (
            self._montecarlo_virtual_luminosity[:-1]
            / self.transport_state.time_of_simulation.value
        )

    def generate_virtual_spectrum(
        self,
        number_of_vpackets,
        enable_full_relativity=False,
        tau_russian=10.0,
        survival_probability=0.0,
        v_packet_spawn_start_frequency=0.0,
        v_packet_spawn_end_frequency=np.inf,
        temporary_v_packet_bins=10000,
    ):
        """
        Generate virtual packets via postprocessing from tracker data.

        This method creates a VirtualPacketSolver, extracts spawn events from
        the tracker data, and generates virtual packets in postprocessing.

        Parameters
        ----------
        number_of_vpackets : int
            Number of virtual packets to generate per spawn event
        enable_full_relativity : bool, optional
            Enable full relativistic effects, by default False
        tau_russian : float, optional
            Russian roulette optical depth threshold, by default 10.0
        survival_probability : float, optional
            Survival probability for Russian roulette, by default 0.0
        v_packet_spawn_start_frequency : float, optional
            Start frequency for vpacket spawning [Hz], by default 0.0
        v_packet_spawn_end_frequency : float, optional
            End frequency for vpacket spawning [Hz], by default np.inf
        temporary_v_packet_bins : int, optional
            Initial size of temporary storage arrays, by default 10000

        Returns
        -------
        VirtualPacketState
            State containing virtual packet data and spectrum

        Notes
        -----
        Requires that the simulation was run with tracking enabled.
        The tracker_full_df must be available in transport_state.
        """
        if self.transport_state.tracker_full_df is None:
            raise ValueError(
                "Cannot generate virtual packets: tracker_full_df is None. "
                "Please run the simulation with tracking enabled "
                "(set rpacket_tracking to True in the montecarlo configuration)."
            )

        if number_of_vpackets == 0:
            warnings.warn(
                "number_of_vpackets is 0. No virtual packets will be generated.",
                UserWarning,
            )

        # Create solver if not already created or if parameters changed
        self.virtual_packet_solver = VirtualPacketSolver(
            enable_full_relativity=enable_full_relativity,
            tau_russian=tau_russian,
            survival_probability=survival_probability,
            v_packet_spawn_start_frequency=v_packet_spawn_start_frequency,
            v_packet_spawn_end_frequency=v_packet_spawn_end_frequency,
            number_of_vpackets=number_of_vpackets,
            temporary_v_packet_bins=temporary_v_packet_bins,
        )

        # Generate virtual packets
        self.virtual_packet_state = (
            self.virtual_packet_solver.generate_virtual_packets(
                self.transport_state.tracker_full_df,
                self.transport_state.geometry_state,
                self.transport_state.opacity_state,
                self.transport_state.time_explosion.cgs.value,
                self.spectrum_frequency_grid.value,
                self.transport_state.time_of_simulation.cgs.value,
            )
        )

        return self.virtual_packet_state

    def solve(self, transport_state):
        """Solve the spectra

        Parameters
        ----------
        transport_state: MonteCarloTransportState
            The transport state to be used to compute the spectra

        Returns
        -------
        tuple(TARDISSpectrum)
            Real, virtual and integrated spectra, if available
        """
        self.transport_state = transport_state

        return (
            self.spectrum_real_packets,
            self.spectrum_virtual_packets,
            self.spectrum_integrated,
        )

    @classmethod
    def from_config(cls, config):
        spectrum_frequency_grid = quantity_linspace(
            config.spectrum.stop.to("Hz", u.spectral()),
            config.spectrum.start.to("Hz", u.spectral()),
            num=config.spectrum.num + 1,
        )

        return cls(
            transport_state=None,
            spectrum_frequency_grid=spectrum_frequency_grid,
            integrator_settings=config.spectrum.integrated,
        )
