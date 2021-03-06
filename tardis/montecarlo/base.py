from astropy import units as u, constants as const

from scipy.special import zeta

from tardis.montecarlo import montecarlo, packet_source


import numpy as np

class MontecarloRunner(object):
    """
    This class is designed as an interface between the Python part and the
    montecarlo C-part
    """

    w_estimator_constant = ((const.c ** 2 / (2 * const.h)) *
                            (15 / np.pi ** 4) * (const.h / const.k_B) ** 4 /
                            (4 * np.pi)).cgs.value

    t_rad_estimator_constant = ((np.pi**4 / (15 * 24 * zeta(5, 1))) *
                                (const.h / const.k_B)).cgs.value


    def __init__(self, seed):
        self.packet_source = packet_source.BlackBodySimpleSource(seed)



    def _initialize_montecarlo_arrays(self, model):
        """
        Initialize the output arrays of the montecarlo simulation.

        Parameters
        ----------

        model: ~Radial1DModel
        """

        no_of_packets = model.current_no_of_packets
        no_of_shells = model.tardis_config.structure.no_of_shells
        self._output_nu = np.ones(no_of_packets, dtype=np.float64) * -99.0
        self._output_energy = np.ones(no_of_packets, dtype=np.float64) * -99.0
        self.last_line_interaction_in_id = -1 * np.ones(no_of_packets, dtype=np.int64)
        self.last_line_interaction_out_id = -1 * np.ones(no_of_packets, dtype=np.int64)
        self.last_line_interaction_shell_id = -1 * np.ones(no_of_packets, dtype=np.int64)
        self.last_interaction_type = -1 * np.ones(no_of_packets, dtype=np.int64)
        self.last_interaction_in_nu = np.zeros(no_of_packets, dtype=np.float64)

        #Estimators
        self.j_estimator = np.zeros(no_of_shells, dtype=np.float64)
        self.nu_bar_estimator = np.zeros(no_of_shells, dtype=np.float64)
        self.j_blue_estimator = np.zeros_like(
            model.plasma_array.tau_sobolevs.values)


    def _initialize_geometry_arrays(self, structure):
        """
        Generate the cgs like geometry arrays for the montecarlo part

        Parameters
        ----------

        structure: ~ConfigurationNameSpace
        """
        self.r_inner_cgs = structure.r_inner.to('cm').value
        self.r_outer_cgs = structure.r_outer.to('cm').value
        self.v_inner_cgs = structure.v_inner.to('cm/s').value

    def _initialize_packets(self, T, no_of_packets):
        nus, mus, energies = self.packet_source.create_packets(T, no_of_packets)
        self.input_nu = nus
        self.input_mu = mus
        self.input_energy = energies


    def run(self, model, no_of_virtual_packets, nthreads=1):
        """
        Running the TARDIS simulation

        Parameters
        ----------

        :param model:
        :param no_of_virtual_packets:
        :param nthreads:
        :return:
        """
        self.time_of_simulation = model.time_of_simulation
        self.volume = model.tardis_config.structure.volumes
        self._initialize_montecarlo_arrays(model)
        self._initialize_geometry_arrays(model.tardis_config.structure)

        self._initialize_packets(model.t_inner.value,
                                 model.current_no_of_packets)
        montecarlo.montecarlo_radial1d(
            model, self, virtual_packet_flag=no_of_virtual_packets,
            nthreads=nthreads)

    def legacy_return(self):
        return (self.output_nu, self.output_energy,
                self.j_estimator, self.nu_bar_estimator,
                self.last_line_interaction_in_id,
                self.last_line_interaction_out_id,
                self.last_interaction_type,
                self.last_line_interaction_shell_id)

    def get_line_interaction_id(self, line_interaction_type):
        return ['scatter', 'downbranch', 'macroatom'].index(
            line_interaction_type)


    @property
    def output_nu(self):
        return u.Quantity(self._output_nu, u.Hz)

    @property
    def output_energy(self):
        return u.Quantity(self._output_energy, u.erg)

    @property
    def packet_luminosity(self):
        return self.output_energy / self.time_of_simulation

    @property
    def emitted_packet_mask(self):
        return self.output_energy >=0

    @property
    def emitted_packet_nu(self):
        return self.output_nu[self.emitted_packet_mask]

    @property
    def reabsorbed_packet_nu(self):
        return self.output_nu[~self.emitted_packet_mask]

    @property
    def reabsorbed_packet_luminosity(self):
        return -self.packet_luminosity[~self.emitted_packet_mask]


    @property
    def emitted_packet_luminosity(self):
        return self.packet_luminosity[self.emitted_packet_mask]

    @property
    def reabsorbed_packet_luminosity(self):
        return -self.packet_luminosity[~self.emitted_packet_mask]

    def calculate_radiationfield_properties(self):
        """
        Calculate an updated radiation field from the :math:`\\bar{nu}_\\textrm{estimator}` and :math:`\\J_\\textrm{estimator}`
        calculated in the montecarlo simulation. The details of the calculation can be found in the documentation.

        Parameters
        ----------

        nubar_estimator : ~np.ndarray (float)

        j_estimator : ~np.ndarray (float)

        Returns
        -------

        t_rad : ~astropy.units.Quantity (float)

        w : ~numpy.ndarray (float)

        """


        t_rad = (self.t_rad_estimator_constant * self.nu_bar_estimator
                / self.j_estimator)
        w = self.j_estimator / (4 * const.sigma_sb.cgs.value * t_rad ** 4
                                * self.time_of_simulation.value
                                * self.volume.value)

        return t_rad * u.K, w
