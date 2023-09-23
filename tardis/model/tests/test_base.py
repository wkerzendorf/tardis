import os
from pathlib import Path
import pytest
import pandas as pd
from astropy import units as u
from numpy.testing import assert_almost_equal, assert_array_almost_equal

from tardis.io.configuration.config_reader import Configuration
from tardis.model import SimulationState
from tardis.io.decay import IsotopeAbundances


class TestModelFromPaper1Config:
    @pytest.fixture(autouse=True)
    def setup(self, example_configuration_dir):
        self.config = Configuration.from_yaml(
            example_configuration_dir / "paper1_tardis_configv1.yml"
        )
        self.simulation_state = SimulationState.from_config(self.config)

    def test_abundances(self):
        oxygen_abundance = self.config.model.abundances.O
        assert_array_almost_equal(
            oxygen_abundance, self.simulation_state.abundance.loc[8].values
        )

    def test_velocities(self):
        velocity = self.config.model.structure.velocity
        assert_almost_equal(
            velocity.start.cgs.value, self.simulation_state.v_inner[0].cgs.value
        )
        assert_almost_equal(
            velocity.stop.cgs.value, self.simulation_state.v_outer[-1].cgs.value
        )
        assert len(self.simulation_state.v_outer) == velocity.num

    def test_densities(self):
        assert_almost_equal(
            self.simulation_state.density[0].cgs.value,
            (7.542803599143591e-14 * u.Unit("g/cm^3")).value,
        )
        assert_almost_equal(
            self.simulation_state.density[-1].cgs.value,
            (1.432259798833509e-15 * u.Unit("g/cm^3")).value,
        )

    def test_time_explosion(self):
        assert_almost_equal(
            self.simulation_state.time_explosion.to(u.day).value, 13.0
        )


class TestModelFromASCIIDensity:
    @pytest.fixture(autouse=True)
    def setup(self, example_model_file_dir):
        self.config = Configuration.from_yaml(
            example_model_file_dir / "tardis_configv1_ascii_density.yml"
        )
        self.simulation_state = SimulationState.from_config(self.config)

    def test_velocities(self):
        assert self.simulation_state.v_inner.unit == u.Unit("cm/s")
        assert_almost_equal(self.simulation_state.v_inner[0].value, 1e4 * 1e5)

    def test_abundances(self):
        oxygen_abundance = self.config.model.abundances.O
        assert_array_almost_equal(
            oxygen_abundance, self.simulation_state.abundance.loc[8].values
        )


class TestModelFromArtisDensity:
    @pytest.fixture(autouse=True)
    def setup(self, example_model_file_dir):
        self.config = Configuration.from_yaml(
            example_model_file_dir / "tardis_configv1_artis_density.yml"
        )
        self.simulation_state = SimulationState.from_config(self.config)

    def test_velocities(self):
        assert self.simulation_state.v_inner.unit == u.Unit("cm/s")
        assert_almost_equal(
            self.simulation_state.v_inner[0].value, 1.259375e03 * 1e5
        )

    def test_abundances(self):
        oxygen_abundance = self.config.model.abundances.O
        assert_array_almost_equal(
            oxygen_abundance, self.simulation_state.abundance.loc[8].values
        )


class TestModelFromArtisDensityAbundances:
    @pytest.fixture(autouse=True)
    def setup(self, example_model_file_dir):
        self.config = Configuration.from_yaml(
            example_model_file_dir / "tardis_configv1_artis_density.yml"
        )
        self.config.model.abundances.type = "file"
        self.config.model.abundances.filename = "artis_abundances.dat"
        self.config.model.abundances.filetype = "artis"
        self.simulation_state = SimulationState.from_config(self.config)

    def test_velocities(self):
        assert self.simulation_state.v_inner.unit == u.Unit("cm/s")
        assert_almost_equal(
            self.simulation_state.v_inner[0].value, 1.259375e03 * 1e5
        )

    def test_abundances(self):
        assert_almost_equal(
            self.simulation_state.abundance.loc[14, 54], 0.21864420000000001
        )


class TestModelFromArtisDensityAbundancesVSlice:
    @pytest.fixture(autouse=True)
    def setup(self, example_model_file_dir):

        self.config = Configuration.from_yaml(
            example_model_file_dir / "tardis_configv1_artis_density_v_slice.yml"
        )
        self.config.model.abundances.type = "file"
        self.config.model.abundances.filename = "artis_abundances.dat"
        self.config.model.abundances.filetype = "artis"
        self.simulation_state = SimulationState.from_config(self.config)

    def test_velocities(self):
        assert self.simulation_state.v_inner.unit == u.Unit("cm/s")
        assert_almost_equal(
            self.simulation_state.v_inner[0].to(u.km / u.s).value, 9000
        )

    def test_abundances(self):
        assert_almost_equal(
            self.simulation_state.abundance.loc[14, 31], 2.156751e-01
        )


class TestModelFromUniformDensity:
    @pytest.fixture(autouse=True)
    def setup(self, example_configuration_dir):
        self.config = Configuration.from_yaml(
            example_configuration_dir / "tardis_configv1_uniform_density.yml"
        )
        self.simulation_state = SimulationState.from_config(self.config)

    def test_density(self):
        assert_array_almost_equal(
            self.simulation_state.density.to(u.Unit("g / cm3")).value, 1.0e-14
        )


class TestModelFromInitialTinner:
    @pytest.fixture(autouse=True)
    def setup(self, example_configuration_dir):
        self.config = Configuration.from_yaml(
            example_configuration_dir / "tardis_configv1_uniform_density.yml"
        )
        self.config.plasma.initial_t_inner = 2508 * u.K
        self.simulation_state = SimulationState.from_config(self.config)

    def test_initial_temperature(self):
        assert_almost_equal(self.simulation_state.t_inner.value, 2508)


class TestModelFromArtisDensityAbundancesAllAscii:
    @pytest.fixture(autouse=True)
    def setup(self, example_model_file_dir):
        self.config = Configuration.from_yaml(
            example_model_file_dir / "tardis_configv1_ascii_density_abund.yml"
        )
        self.config.model.structure.filename = "density.dat"
        self.config.model.abundances.filename = "abund.dat"
        self.simulation_state = SimulationState.from_config(self.config)

    def test_velocities(self):
        assert self.simulation_state.v_inner.unit == u.Unit("cm/s")
        assert_almost_equal(
            self.simulation_state.v_inner[0].to(u.km / u.s).value, 11000
        )

    def test_abundances(self):
        assert_almost_equal(self.simulation_state.abundance.loc[14, 0], 0.1)
        assert_almost_equal(self.simulation_state.abundance.loc[14, 1], 0.2)
        assert_almost_equal(self.simulation_state.abundance.loc[14, 2], 0.2)
        assert_almost_equal(self.simulation_state.abundance.loc[14, 3], 0.2)
        assert_almost_equal(self.simulation_state.abundance.loc[14, 4], 0.2)
        assert_almost_equal(self.simulation_state.abundance.loc[14, 5], 0.2)
        assert_almost_equal(self.simulation_state.abundance.loc[14, 6], 0.0)
        assert_almost_equal(self.simulation_state.abundance.loc[6, 0], 0.0)
        assert_almost_equal(self.simulation_state.abundance.loc[6, 1], 0.0)
        assert_almost_equal(self.simulation_state.abundance.loc[6, 2], 0.0)
        assert_almost_equal(self.simulation_state.abundance.loc[6, 3], 0.0)
        assert_almost_equal(self.simulation_state.abundance.loc[6, 4], 0.0)
        assert_almost_equal(self.simulation_state.abundance.loc[6, 5], 0.0)
        assert_almost_equal(self.simulation_state.abundance.loc[6, 6], 0.5)

    def test_densities(self):
        assert_almost_equal(
            self.simulation_state.density[0].to(u.Unit("g/cm3")).value,
            9.7656229e-11 / 13.0**3,
        )
        assert_almost_equal(
            self.simulation_state.density[1].to(u.Unit("g/cm3")).value,
            4.8170911e-11 / 13.0**3,
        )
        assert_almost_equal(
            self.simulation_state.density[2].to(u.Unit("g/cm3")).value,
            2.5600000e-11 / 13.0**3,
        )
        assert_almost_equal(
            self.simulation_state.density[3].to(u.Unit("g/cm3")).value,
            1.4450533e-11 / 13.0**3,
        )
        assert_almost_equal(
            self.simulation_state.density[4].to(u.Unit("g/cm3")).value,
            8.5733893e-11 / 13.0**3,
        )
        assert_almost_equal(
            self.simulation_state.density[5].to(u.Unit("g/cm3")).value,
            5.3037103e-11 / 13.0**3,
        )
        assert_almost_equal(
            self.simulation_state.density[6].to(u.Unit("g/cm3")).value,
            3.3999447e-11 / 13.0**3,
        )


def test_ascii_reader_power_law(example_configuration_dir):
    config = Configuration.from_yaml(
        example_configuration_dir / "tardis_configv1_density_power_law_test.yml"
    )
    simulation_state = SimulationState.from_config(config)

    expected_densites = [
        3.29072513e-14,
        2.70357804e-14,
        2.23776573e-14,
        1.86501954e-14,
        1.56435277e-14,
        1.32001689e-14,
        1.12007560e-14,
        9.55397475e-15,
        8.18935779e-15,
        7.05208050e-15,
        6.09916083e-15,
        5.29665772e-15,
        4.61758699e-15,
        4.04035750e-15,
        3.54758837e-15,
        3.12520752e-15,
        2.76175961e-15,
        2.44787115e-15,
        2.17583442e-15,
        1.93928168e-15,
    ]

    assert simulation_state.no_of_shells == 20
    for i, mdens in enumerate(expected_densites):
        assert_almost_equal(
            simulation_state.density[i].to(u.Unit("g / (cm3)")).value, mdens
        )


def test_ascii_reader_exponential_law(example_configuration_dir):
    config = Configuration.from_yaml(
        example_configuration_dir
        / "tardis_configv1_density_exponential_test.yml"
    )
    simulation_state = SimulationState.from_config(config)

    expected_densites = [
        5.18114795e-14,
        4.45945537e-14,
        3.83828881e-14,
        3.30364579e-14,
        2.84347428e-14,
        2.44740100e-14,
        2.10649756e-14,
        1.81307925e-14,
        1.56053177e-14,
        1.34316215e-14,
        1.15607037e-14,
        9.95038990e-15,
        8.56437996e-15,
        7.37143014e-15,
        6.34464872e-15,
        5.46088976e-15,
        4.70023138e-15,
        4.04552664e-15,
        3.48201705e-15,
        2.99699985e-15,
    ]
    expected_unit = "g / (cm3)"

    assert simulation_state.no_of_shells == 20
    for i, mdens in enumerate(expected_densites):
        assert_almost_equal(simulation_state.density[i].value, mdens)
        assert simulation_state.density[i].unit == u.Unit(expected_unit)


@pytest.fixture
def simple_isotope_abundance():
    index = pd.MultiIndex.from_tuples(
        [(6, 14), (12, 28)], names=["atomic_number", "mass_number"]
    )
    abundance = [[0.2] * 20] * 2
    return IsotopeAbundances(abundance, index=index)


def test_model_decay(simple_isotope_abundance, example_configuration_dir):
    config = Configuration.from_yaml(
        example_configuration_dir / "tardis_configv1_verysimple.yml"
    )
    simulation_state = SimulationState.from_config(config)

    simulation_state.raw_isotope_abundance = simple_isotope_abundance
    decayed = simple_isotope_abundance.decay(
        simulation_state.time_explosion
    ).as_atoms()
    norm_factor = 1.4

    assert_almost_equal(
        simulation_state.abundance.loc[8][0],
        simulation_state.raw_abundance.loc[8][0] / norm_factor,
        decimal=4,
    )
    assert_almost_equal(
        simulation_state.abundance.loc[14][0],
        (simulation_state.raw_abundance.loc[14][0] + decayed.loc[14][0])
        / norm_factor,
        decimal=4,
    )
    assert_almost_equal(
        simulation_state._abundance.loc[12][5],
        (simulation_state.raw_abundance.loc[12][5] + decayed.loc[12][5])
        / norm_factor,
        decimal=4,
    )
    assert_almost_equal(
        simulation_state.abundance.loc[6][12],
        (decayed.loc[6][12]) / norm_factor,
        decimal=4,
    )


@pytest.mark.parametrize(
    ("index", "expected"),
    [
        (0, 1.00977478e45),
        (10, 1.98154804e45),
        (19, 3.13361319e45),
    ],
)
def test_radial_1D_geometry_volume(simulation_verysimple, index, expected):
    geometry = simulation_verysimple.simulation_state.model_state.geometry
    volume = geometry.volume

    assert volume.unit == u.Unit("cm3")
    assert_almost_equal(volume[index].value, expected, decimal=-40)


@pytest.mark.parametrize(
    ("index", "expected"),
    [
        ((8, 0), 539428198),
        ((8, 1), 409675383),
        ((8, 2), 314387928),
        ((12, 0), 56066111),
        ((12, 1), 42580098),
        ((12, 2), 32676283),
        ((14, 0), 841032262),
        ((14, 1), 638732300),
        ((14, 2), 490167906),
        ((16, 0), 269136275),
        ((16, 1), 204398856),
        ((16, 2), 156857199),
        ((18, 0), 45482957),
        ((18, 1), 34542591),
        ((18, 2), 26508241),
        ((20, 0), 34001569),
        ((20, 1), 25822910),
        ((20, 2), 19816693),
    ],
)
def test_composition_elemental_number_density(
    simulation_verysimple, index, expected
):
    comp = simulation_verysimple.simulation_state.model_state.composition

    assert_almost_equal(
        comp.elemental_number_density.loc[index], expected, decimal=-2
    )


@pytest.mark.parametrize(
    ("index", "expected"),
    [
        ((8, 0), 1.4471412e31),
        ((16, 10), 2.6820129e30),
        ((20, 19), 1.3464444e29),
    ],
)
def test_model_state_mass(simulation_verysimple, index, expected):
    model_state = simulation_verysimple.simulation_state.model_state

    assert_almost_equal((model_state.mass).loc[index], expected, decimal=-27)


@pytest.mark.parametrize(
    ("index", "expected"),
    [
        ((8, 0), 5.4470099e53),
        ((16, 10), 5.0367073e52),
        ((20, 19), 2.0231745e51),
    ],
)
def test_model_state_number(simulation_verysimple, index, expected):
    model_state = simulation_verysimple.simulation_state.model_state

    assert_almost_equal((model_state.number).loc[index], expected, decimal=-47)


@pytest.fixture
def non_uniform_model_state(atomic_dataset, example_model_file_dir):
    config = Configuration.from_yaml(
        example_model_file_dir / "tardis_configv1_isotope_iabund.yml"
    )
    atom_data = atomic_dataset
    model = SimulationState.from_config(config, atom_data=atom_data)
    return model.model_state


@pytest.mark.parametrize(
    ("index", "expected"),
    [
        ((1, 0), 1.67378172e-24),
        ((28, 0), 9.51707707e-23),
        ((28, 1), 9.54725917e-23),
    ],
)
def test_radial_1d_model_atomic_mass(non_uniform_model_state, index, expected):
    atomic_mass = non_uniform_model_state.composition.atomic_mass

    assert_almost_equal(
        atomic_mass.loc[index],
        expected,
        decimal=30,
    )


class TestModelStateFromNonUniformAbundances:
    @pytest.fixture
    def model_state(self, non_uniform_model_state):
        return non_uniform_model_state

    def test_atomic_mass(self, model_state):
        atomic_mass = model_state.composition.atomic_mass
        assert_almost_equal(atomic_mass.loc[(1, 0)], 1.67378172e-24, decimal=30)
        assert_almost_equal(
            atomic_mass.loc[(28, 0)], 9.51707707e-23, decimal=30
        )
        assert_almost_equal(
            atomic_mass.loc[(28, 1)], 9.54725917e-23, decimal=30
        )

    def test_elemental_number_density(self, model_state):
        number = model_state.composition.elemental_number_density
        assert_almost_equal(number.loc[(1, 0)], 0)
        assert_almost_equal(number.loc[(28, 0)], 10825427.035, decimal=2)
        assert_almost_equal(number.loc[(28, 1)], 1640838.763, decimal=2)

    def test_number(self, model_state):
        number = model_state.number
        assert_almost_equal(number.loc[(1, 0)], 0)
        assert_almost_equal(number.loc[(28, 0)], 1.53753476e53, decimal=-47)
        assert_almost_equal(number.loc[(28, 1)], 4.16462779e52, decimal=-47)


###
# Save and Load
###


@pytest.fixture(scope="module", autouse=True)
def to_hdf_buffer(hdf_file_path, simulation_verysimple):
    simulation_verysimple.simulation_state.to_hdf(hdf_file_path, overwrite=True)


model_scalar_attrs = ["t_inner"]


@pytest.mark.parametrize("attr", model_scalar_attrs)
def test_hdf_model_scalars(hdf_file_path, simulation_verysimple, attr):
    path = os.path.join("model", "scalars")
    expected = pd.read_hdf(hdf_file_path, path)[attr]
    actual = getattr(simulation_verysimple.simulation_state, attr)
    if hasattr(actual, "cgs"):
        actual = actual.cgs.value
    assert_almost_equal(actual, expected)


model_nparray_attrs = ["w", "v_inner", "v_outer"]


@pytest.mark.parametrize("attr", model_nparray_attrs)
def test_hdf_model_nparray(hdf_file_path, simulation_verysimple, attr):
    path = os.path.join("model", attr)
    expected = pd.read_hdf(hdf_file_path, path)
    actual = getattr(simulation_verysimple.simulation_state, attr)
    if hasattr(actual, "cgs"):
        actual = actual.cgs.value
    assert_almost_equal(actual, expected.values)
