import warnings

import numpy as np
from astropy import units as u

from tardis.io.util import HDFWriterMixin
from tardis.spectrum.formal_integral import IntegrationError
from tardis.spectrum.spectrum import TARDISSpectrum
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

    def __init__(self, transport_state, spectrum_frequency_grid):
        self.transport_state = transport_state
        self.spectrum_frequency_grid = spectrum_frequency_grid
        self._montecarlo_virtual_luminosity = u.Quantity(
            np.zeros_like(self.spectrum_frequency_grid.value), "erg / s"
        )  # should be init with v_packets_energy_hist
        self._integrator = None
        self.integrator_settings = None
        self._spectrum_integrated = None

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
        if np.all(self.montecarlo_virtual_luminosity == 0):
            warnings.warn(
                "SpectrumSolver.spectrum_virtual_packets "
                "is zero. Please run the montecarlo simulation with "
                "no_of_virtual_packets > 0",
                UserWarning,
            )

        return TARDISSpectrum(
            self.spectrum_frequency_grid, self.montecarlo_virtual_luminosity
        )

    @property
    def spectrum_integrated(self):
        if self._spectrum_integrated is None:
            # This was changed from unpacking to specific attributes as compute
            # is not used in calculate_spectrum
            try:
                self._spectrum_integrated = self.integrator.calculate_spectrum(
                    self.spectrum_frequency_grid[:-1],
                    points=self.integrator_settings.points,
                    interpolate_shells=self.integrator_settings.interpolate_shells,
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
                return TARDISSpectrum(
                    np.array([np.nan, np.nan]) * u.Hz,
                    np.array([np.nan]) * u.erg / u.s,
                )
        return self._spectrum_integrated

    @property
    def integrator(self):
        if self._integrator is None:
            warnings.warn(
                "MontecarloTransport.integrator: "
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

    def calculate_emitted_luminosity(
        self, luminosity_nu_start, luminosity_nu_end
    ):
        """
        Calculate emitted luminosity.

        Parameters
        ----------
        luminosity_nu_start : astropy.units.Quantity
        luminosity_nu_end : astropy.units.Quantity

        Returns
        -------
        astropy.units.Quantity
        """
        luminosity_wavelength_filter = (
            self.transport_state.emitted_packet_nu > luminosity_nu_start
        ) & (self.transport_state.emitted_packet_nu < luminosity_nu_end)

        return self.transport_state.emitted_packet_luminosity[
            luminosity_wavelength_filter
        ].sum()

    def calculate_reabsorbed_luminosity(
        self, luminosity_nu_start, luminosity_nu_end
    ):
        """
        Calculate reabsorbed luminosity.

        Parameters
        ----------
        luminosity_nu_start : astropy.units.Quantity
        luminosity_nu_end : astropy.units.Quantity

        Returns
        -------
        astropy.units.Quantity
        """
        luminosity_wavelength_filter = (
            self.transport_state.reabsorbed_packet_nu > luminosity_nu_start
        ) & (self.transport_state.reabsorbed_packet_nu < luminosity_nu_end)

        return self.transport_state.reabsorbed_packet_luminosity[
            luminosity_wavelength_filter
        ].sum()

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
        )
