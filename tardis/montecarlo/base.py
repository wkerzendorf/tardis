import logging
import warnings

from astropy import units as u
from tardis import constants as const
from numba import set_num_threads
from numba import cuda

from tardis.montecarlo.spectrum import TARDISSpectrum

from tardis.util.base import quantity_linspace
from tardis.io.util import HDFWriterMixin
from tardis.montecarlo import packet_source as source
from tardis.montecarlo.montecarlo_numba.formal_integral import FormalIntegrator
from tardis.montecarlo.montecarlo_numba.estimators import initialize_estimators
from tardis.montecarlo import montecarlo_configuration as mc_config_module
from tardis.montecarlo.montecarlo_state import MonteCarloTransportState

from tardis.montecarlo.montecarlo_numba import montecarlo_radial1d
from tardis.montecarlo.montecarlo_numba.numba_interface import (
    configuration_initialize,
)
from tardis.montecarlo.montecarlo_numba import numba_config
from tardis.io.logger import montecarlo_tracking as mc_tracker

import numpy as np

logger = logging.getLogger(__name__)


# TODO: refactor this into more parts
class MontecarloTransportSolver(HDFWriterMixin):
    """
    This class is designed as an interface between the Python part and the
    montecarlo C-part
    """

    hdf_properties = [
        "output_nu",
        "output_energy",
        "nu_bar_estimator",
        "j_estimator",
        "montecarlo_virtual_luminosity",
        "last_interaction_in_nu",
        "last_interaction_type",
        "last_line_interaction_in_id",
        "last_line_interaction_out_id",
        "last_line_interaction_shell_id",
        "packet_luminosity",
        "spectrum",
        "spectrum_virtual",
        "spectrum_reabsorbed",
        "time_of_simulation",
        "emitted_packet_mask",
    ]

    vpacket_hdf_properties = [
        "virt_packet_nus",
        "virt_packet_energies",
        "virt_packet_initial_rs",
        "virt_packet_initial_mus",
        "virt_packet_last_interaction_in_nu",
        "virt_packet_last_interaction_type",
        "virt_packet_last_line_interaction_in_id",
        "virt_packet_last_line_interaction_out_id",
        "virt_packet_last_line_interaction_shell_id",
    ]

    hdf_name = "transport"

    def __init__(
        self,
        spectrum_frequency,
        virtual_spectrum_spawn_range,
        disable_electron_scattering,
        enable_reflective_inner_boundary,
        enable_full_relativity,
        inner_boundary_albedo,
        line_interaction_type,
        integrator_settings,
        v_packet_settings,
        spectrum_method,
        packet_source,
        enable_virtual_packet_logging,
        nthreads=1,
        debug_packets=False,
        logger_buffer=1,
        tracking_rpacket=False,
        use_gpu=False,
    ):
        # inject different packets
        self.disable_electron_scattering = disable_electron_scattering
        self.spectrum_frequency = spectrum_frequency
        self.virtual_spectrum_spawn_range = virtual_spectrum_spawn_range
        self.enable_reflective_inner_boundary = enable_reflective_inner_boundary
        self.inner_boundary_albedo = inner_boundary_albedo
        self.enable_full_relativity = enable_full_relativity
        numba_config.ENABLE_FULL_RELATIVITY = enable_full_relativity
        self.line_interaction_type = line_interaction_type
        self.integrator_settings = integrator_settings
        self.v_packet_settings = v_packet_settings
        self.spectrum_method = spectrum_method
        self._integrator = None
        self._spectrum_integrated = None
        self.use_gpu = use_gpu

        self.virt_logging = enable_virtual_packet_logging
        self.virt_packet_last_interaction_type = np.ones(2) * -1
        self.virt_packet_last_interaction_in_nu = np.ones(2) * -1.0
        self.virt_packet_last_line_interaction_in_id = np.ones(2) * -1
        self.virt_packet_last_line_interaction_out_id = np.ones(2) * -1
        self.virt_packet_last_line_interaction_shell_id = np.ones(2) * -1
        self.virt_packet_nus = np.ones(2) * -1.0
        self.virt_packet_energies = np.ones(2) * -1.0
        self.virt_packet_initial_rs = np.ones(2) * -1.0
        self.virt_packet_initial_mus = np.ones(2) * -1.0

        self.packet_source = packet_source

        # Setting up the Tracking array for storing all the RPacketTracker instances
        self.rpacket_tracker = None

        # Set number of threads
        self.nthreads = nthreads

        # set up logger based on config
        mc_tracker.DEBUG_MODE = debug_packets
        mc_tracker.BUFFER = logger_buffer

        mc_config_module.RPACKET_TRACKING = tracking_rpacket

        if self.spectrum_method == "integrated":
            self.optional_hdf_properties.append("spectrum_integrated")

    def _initialize_packets(self, simulation_state, no_of_packets, iteration):
        # the iteration (passed as seed_offset) is added each time to preserve randomness
        # across different simulations with the same temperature,
        # for example.
        mc_config_module.packet_seeds = self.packet_source.create_packet_seeds(
            no_of_packets, iteration
        )

        # Create packets
        self.packet_collection = self.packet_source.create_packets(
            no_of_packets
        )

        self._output_nu = self.packet_collection.output_nus
        self._output_energy = self.packet_collection.output_energies

        self.last_line_interaction_in_id = -1 * np.ones(
            no_of_packets, dtype=np.int64
        )
        self.last_line_interaction_out_id = -1 * np.ones(
            no_of_packets, dtype=np.int64
        )
        self.last_line_interaction_shell_id = -1 * np.ones(
            no_of_packets, dtype=np.int64
        )
        self.last_interaction_type = -1 * np.ones(no_of_packets, dtype=np.int64)
        self.last_interaction_in_nu = np.zeros(no_of_packets, dtype=np.float64)

        self._montecarlo_virtual_luminosity = u.Quantity(
            np.zeros_like(self.spectrum_frequency.value), "erg / s"
        )

    @property
    def spectrum(self):
        return TARDISSpectrum(
            self.spectrum_frequency, self.montecarlo_emitted_luminosity
        )

    @property
    def spectrum_reabsorbed(self):
        return TARDISSpectrum(
            self.spectrum_frequency, self.montecarlo_reabsorbed_luminosity
        )

    @property
    def spectrum_virtual(self):
        if np.all(self.montecarlo_virtual_luminosity == 0):
            warnings.warn(
                "MontecarloTransport.spectrum_virtual"
                "is zero. Please run the montecarlo simulation with"
                "no_of_virtual_packets > 0",
                UserWarning,
            )

        return TARDISSpectrum(
            self.spectrum_frequency, self.montecarlo_virtual_luminosity
        )

    @property
    def spectrum_integrated(self):
        if self._spectrum_integrated is None:
            # This was changed from unpacking to specific attributes as compute
            # is not used in calculate_spectrum
            self._spectrum_integrated = self.integrator.calculate_spectrum(
                self.spectrum_frequency[:-1],
                points=self.integrator_settings.points,
                interpolate_shells=self.integrator_settings.interpolate_shells,
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
        if self.enable_full_relativity:
            raise NotImplementedError(
                "The FormalIntegrator is not yet implemented for the full "
                "relativity mode. "
                "Please run with config option enable_full_relativity: "
                "False."
            )
        return self._integrator

    def run(
        self,
        simulation_state,
        plasma,
        no_of_packets,
        no_of_virtual_packets=0,
        iteration=0,
        total_iterations=0,
        show_progress_bars=True,
    ):
        """
        Run the montecarlo calculation

        Parameters
        ----------
        model : tardis.model.SimulationState
        plasma : tardis.plasma.BasePlasma
        no_of_packets : int
        no_of_virtual_packets : int
        total_iterations : int
            The total number of iterations in the simulation.

        Returns
        -------
        None
        """

        set_num_threads(self.nthreads)

        if not plasma.continuum_interaction_species.empty:
            gamma_shape = plasma.gamma.shape
        else:
            gamma_shape = (0, 0)

        # Initializing estimator array
        estimators = initialize_estimators(
            plasma.tau_sobolevs.shape, gamma_shape
        )

        self._initialize_packets(
            simulation_state,
            no_of_packets,
            iteration,
        )

        self.mc_state = MonteCarloTransportState(
            self.packet_collection,
            estimators,
            simulation_state.volume.cgs.copy(),
        )
        configuration_initialize(self, no_of_virtual_packets)
        montecarlo_radial1d(
            simulation_state,
            plasma,
            iteration,
            self.packet_collection,
            self.mc_state.estimators,
            total_iterations,
            show_progress_bars,
            self,
        )
        self._integrator = FormalIntegrator(simulation_state, plasma, self)

    def legacy_return(self):
        return (
            self.mc_state.packet_collection.output_nus,
            self.mc_state.packet_collection.output_energies,
            self.mc_state.estimators.j_estimator,
            self.mc_state.estimators.nu_bar_estimator,
            self.last_line_interaction_in_id,
            self.last_line_interaction_out_id,
            self.last_interaction_type,
            self.last_line_interaction_shell_id,
        )

    def get_line_interaction_id(self, line_interaction_type):
        return ["scatter", "downbranch", "macroatom"].index(
            line_interaction_type
        )

    @property
    def output_nu(self):
        return u.Quantity(self._output_nu, u.Hz)

    @property
    def output_energy(self):
        return u.Quantity(self._output_energy, u.erg)

    @property
    def virtual_packet_nu(self):
        try:
            return u.Quantity(self.virt_packet_nus, u.Hz)
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virtual_packet_nu:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virtual_packet_energy(self):
        try:
            return u.Quantity(self.virt_packet_energies, u.erg)
        except AttributeError:
            warnings.warn(
                "MontecarloTransport.virtual_packet_energy:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def virtual_packet_luminosity(self):
        try:
            return self.virtual_packet_energy / self.time_of_simulation
        except TypeError:
            warnings.warn(
                "MontecarloTransport.virtual_packet_luminosity:"
                "Set 'virtual_packet_logging: True' in the configuration file"
                "to access this property"
                "It should be added under 'virtual' property of 'spectrum' property",
                UserWarning,
            )
            return None

    @property
    def montecarlo_reabsorbed_luminosity(self):
        return u.Quantity(
            np.histogram(
                self.reabsorbed_packet_nu,
                weights=self.reabsorbed_packet_luminosity,
                bins=self.spectrum_frequency,
            )[0],
            "erg / s",
        )

    @property
    def montecarlo_emitted_luminosity(self):
        return u.Quantity(
            np.histogram(
                self.emitted_packet_nu,
                weights=self.emitted_packet_luminosity,
                bins=self.spectrum_frequency,
            )[0],
            "erg / s",
        )

    @property
    def montecarlo_virtual_luminosity(self):
        return (
            self._montecarlo_virtual_luminosity[:-1]
            / self.time_of_simulation.value
        )

    @classmethod
    def from_config(
        cls, config, packet_source, enable_virtual_packet_logging=False
    ):
        """
        Create a new MontecarloTransport instance from a Configuration object.

        Parameters
        ----------
        config : tardis.io.config_reader.Configuration
        virtual_packet_logging : bool

        Returns
        -------
        MontecarloTransport
        """
        if config.plasma.disable_electron_scattering:
            logger.warn(
                "Disabling electron scattering - this is not physical."
                "Likely bug in formal integral - "
                "will not give same results."
            )
            numba_config.SIGMA_THOMSON = 1e-200
            # mc_config_module.disable_electron_scattering = True
        else:
            logger.debug("Electron scattering switched on")
            numba_config.SIGMA_THOMSON = const.sigma_T.to("cm^2").value
            # mc_config_module.disable_electron_scattering = False

        spectrum_frequency = quantity_linspace(
            config.spectrum.stop.to("Hz", u.spectral()),
            config.spectrum.start.to("Hz", u.spectral()),
            num=config.spectrum.num + 1,
        )
        running_mode = config.spectrum.integrated.compute.upper()

        if running_mode == "GPU":
            if cuda.is_available():
                use_gpu = True
            else:
                raise ValueError(
                    """The GPU option was selected for the formal_integral,
                    but no CUDA GPU is available."""
                )
        elif running_mode == "AUTOMATIC":
            if cuda.is_available():
                use_gpu = True
            else:
                use_gpu = False
        elif running_mode == "CPU":
            use_gpu = False
        else:
            raise ValueError(
                """An invalid option for compute was passed. The three
                valid values are 'GPU', 'CPU', and 'Automatic'."""
            )

        mc_config_module.disable_line_scattering = (
            config.plasma.disable_line_scattering
        )

        mc_config_module.INITIAL_TRACKING_ARRAY_LENGTH = (
            config.montecarlo.tracking.initial_array_length
        )

        return cls(
            spectrum_frequency=spectrum_frequency,
            virtual_spectrum_spawn_range=config.montecarlo.virtual_spectrum_spawn_range,
            enable_reflective_inner_boundary=config.montecarlo.enable_reflective_inner_boundary,
            inner_boundary_albedo=config.montecarlo.inner_boundary_albedo,
            enable_full_relativity=config.montecarlo.enable_full_relativity,
            line_interaction_type=config.plasma.line_interaction_type,
            integrator_settings=config.spectrum.integrated,
            v_packet_settings=config.spectrum.virtual,
            spectrum_method=config.spectrum.method,
            disable_electron_scattering=config.plasma.disable_electron_scattering,
            packet_source=packet_source,
            debug_packets=config.montecarlo.debug_packets,
            logger_buffer=config.montecarlo.logger_buffer,
            enable_virtual_packet_logging=(
                config.spectrum.virtual.virtual_packet_logging
                | enable_virtual_packet_logging
            ),
            nthreads=config.montecarlo.nthreads,
            tracking_rpacket=config.montecarlo.tracking.track_rpacket,
            use_gpu=use_gpu,
        )
