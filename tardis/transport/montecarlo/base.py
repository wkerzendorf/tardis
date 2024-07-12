import logging

from astropy import units as u
from numba import cuda, set_num_threads

from tardis import constants as const
from tardis.io.logger import montecarlo_tracking as mc_tracker
from tardis.io.util import HDFWriterMixin
from tardis.transport.montecarlo import (
    montecarlo_main_loop,
)
from tardis.transport.montecarlo.configuration import (
    constants,
    montecarlo_globals,
)
from tardis.transport.montecarlo.configuration.base import (
    configuration_initialize,
    MonteCarloConfiguration,
)
from tardis.transport.montecarlo.estimators.radfield_mc_estimators import (
    initialize_estimator_statistics,
)
from tardis.transport.montecarlo.formal_integral import FormalIntegrator
from tardis.transport.montecarlo.montecarlo_transport_state import (
    MonteCarloTransportState,
)
from tardis.transport.montecarlo.numba_interface import (
    opacity_state_initialize,
)
from tardis.transport.montecarlo.packet_trackers import (
    rpacket_trackers_to_dataframe,
)
from tardis.util.base import (
    quantity_linspace,
    refresh_packet_pbar,
    update_iterations_pbar,
)

logger = logging.getLogger(__name__)


# TODO: refactor this into more parts
class MonteCarloTransportSolver(HDFWriterMixin):
    """
    This class modifies the MonteCarloTransportState to solve the radiative
    transfer problem.
    """

    hdf_properties = ["transport_state"]

    hdf_name = "transport"

    def __init__(
        self,
        spectrum_frequency,
        virtual_spectrum_spawn_range,
        enable_full_relativity,
        line_interaction_type,
        integrator_settings,
        spectrum_method,
        packet_source,
        enable_virtual_packet_logging=False,
        nthreads=1,
        debug_packets=False,
        logger_buffer=1,
        use_gpu=False,
        montecarlo_configuration=None,
    ):
        # inject different packets
        self.spectrum_frequency = spectrum_frequency
        self.virtual_spectrum_spawn_range = virtual_spectrum_spawn_range
        self.enable_full_relativity = enable_full_relativity
        self.line_interaction_type = line_interaction_type
        self.integrator_settings = integrator_settings
        self.spectrum_method = spectrum_method
        self._integrator = None

        self.use_gpu = use_gpu

        self.enable_vpacket_tracking = enable_virtual_packet_logging
        self.montecarlo_configuration = montecarlo_configuration

        self.packet_source = packet_source

        # Setting up the Tracking array for storing all the RPacketTracker instances
        self.rpacket_tracker = None

        # Set number of threads
        self.nthreads = nthreads

        # set up logger based on config
        mc_tracker.DEBUG_MODE = debug_packets
        mc_tracker.BUFFER = logger_buffer

    def initialize_transport_state(
        self,
        simulation_state,
        plasma,
        no_of_packets,
        no_of_virtual_packets=0,
        iteration=0,
    ):
        if not plasma.continuum_interaction_species.empty:
            gamma_shape = plasma.gamma.shape
        else:
            gamma_shape = (0, 0)

        packet_collection = self.packet_source.create_packets(
            no_of_packets, seed_offset=iteration
        )
        estimators = initialize_estimator_statistics(
            plasma.tau_sobolevs.shape, gamma_shape
        )

        geometry_state = simulation_state.geometry.to_numba()
        opacity_state = opacity_state_initialize(
            plasma,
            self.line_interaction_type,
            montecarlo_globals.DISABLE_LINE_SCATTERING,
            montecarlo_globals.CONTINUUM_PROCESSES_ENABLED,
        )
        transport_state = MonteCarloTransportState(
            packet_collection,
            estimators,
            spectrum_frequency=self.spectrum_frequency,
            geometry_state=geometry_state,
            opacity_state=opacity_state,
        )

        transport_state.integrator_settings = self.integrator_settings
        transport_state._integrator = FormalIntegrator(
            simulation_state, plasma, self
        )

        self.montecarlo_configuration.NUMBER_OF_VPACKETS = no_of_virtual_packets
        self.montecarlo_configuration.TEMPORARY_V_PACKET_BINS = (
            no_of_virtual_packets
        )

        return transport_state

    def run(
        self,
        transport_state,
        time_explosion,
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
        self.transport_state = transport_state

        number_of_vpackets = self.montecarlo_configuration.NUMBER_OF_VPACKETS

        (
            v_packets_energy_hist,
            last_interaction_tracker,
            vpacket_tracker,
            rpacket_trackers,
        ) = montecarlo_main_loop(
            transport_state.packet_collection,
            transport_state.geometry_state,
            time_explosion.cgs.value,
            transport_state.opacity_state,
            transport_state.radfield_mc_estimators,
            transport_state.spectrum_frequency.value,
            number_of_vpackets,
            iteration=iteration,
            show_progress_bars=show_progress_bars,
            total_iterations=total_iterations,
            montecarlo_configuration=self.montecarlo_configuration,
        )

        transport_state._montecarlo_virtual_luminosity.value[:] = (
            v_packets_energy_hist
        )
        transport_state.last_interaction_type = last_interaction_tracker.types
        transport_state.last_interaction_in_nu = last_interaction_tracker.in_nus
        transport_state.last_line_interaction_in_id = (
            last_interaction_tracker.in_ids
        )
        transport_state.last_line_interaction_out_id = (
            last_interaction_tracker.out_ids
        )
        transport_state.last_line_interaction_shell_id = (
            last_interaction_tracker.shell_ids
        )

        if montecarlo_globals.ENABLE_VPACKET_TRACKING and (
            number_of_vpackets > 0
        ):
            transport_state.vpacket_tracker = vpacket_tracker

        update_iterations_pbar(1)
        refresh_packet_pbar()
        # Condition for Checking if RPacket Tracking is enabled
        if montecarlo_globals.ENABLE_RPACKET_TRACKING:
            transport_state.rpacket_tracker = rpacket_trackers

        if self.transport_state.rpacket_tracker is not None:
            self.transport_state.rpacket_tracker_df = (
                rpacket_trackers_to_dataframe(
                    self.transport_state.rpacket_tracker
                )
            )
        transport_state.virt_logging = (
            montecarlo_globals.ENABLE_VPACKET_TRACKING
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
            logger.warning(
                "Disabling electron scattering - this is not physical."
                "Likely bug in formal integral - "
                "will not give same results."
            )
            constants.SIGMA_THOMSON = 1e-200
        else:
            logger.debug("Electron scattering switched on")
            constants.SIGMA_THOMSON = const.sigma_T.to("cm^2").value

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
            use_gpu = bool(cuda.is_available())
        elif running_mode == "CPU":
            use_gpu = False
        else:
            raise ValueError(
                """An invalid option for compute was passed. The three
                valid values are 'GPU', 'CPU', and 'Automatic'."""
            )

        montecarlo_configuration = MonteCarloConfiguration()
        configuration_initialize(
            montecarlo_configuration, config, enable_virtual_packet_logging
        )

        return cls(
            spectrum_frequency=spectrum_frequency,
            virtual_spectrum_spawn_range=config.montecarlo.virtual_spectrum_spawn_range,
            enable_full_relativity=config.montecarlo.enable_full_relativity,
            line_interaction_type=config.plasma.line_interaction_type,
            integrator_settings=config.spectrum.integrated,
            spectrum_method=config.spectrum.method,
            packet_source=packet_source,
            debug_packets=config.montecarlo.debug_packets,
            logger_buffer=config.montecarlo.logger_buffer,
            enable_virtual_packet_logging=(
                config.spectrum.virtual.virtual_packet_logging
                | enable_virtual_packet_logging
            ),
            nthreads=config.montecarlo.nthreads,
            use_gpu=use_gpu,
            montecarlo_configuration=montecarlo_configuration,
        )
