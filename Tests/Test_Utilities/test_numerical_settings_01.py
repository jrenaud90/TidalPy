"""Tests that the configured ``[numerical]`` settings reach the C++ code that uses them after ``update_constants``."""
import math

import pytest

import TidalPy
import TidalPy.constants
from TidalPy.constants import update_constants
from TidalPy.Rheology import Maxwell
from TidalPy.Cooling.cooling import ConductiveCooling
from TidalPy.Radiogenics.radiogenics import FixedRadiogenics
from TidalPy.Structures.layers.base import BaseLayer
from TidalPy.Structures.worlds.layered import LayeredWorld
from TidalPy.Rheology import Elastic
from TidalPy.RadialSolver.solver import radial_solver

_DEFAULT_FLOOR = 1.0e-100
# Large enough that a guarded result moves by many orders of magnitude, small enough to stay unphysical.
_TEST_FLOOR = 1.0e-5
_DEFAULT_CONTINUITY_RTOL = 1.0e-6
_DEFAULT_START_RADIUS_FRACTION = 0.90
_PLANET_RADIUS = 1.0e6


@pytest.fixture
def numerical_setter():
    """Set any ``[numerical]`` key for one test and restore every change afterwards."""
    originals = {}

    def set_value(key, value):
        originals.setdefault(key, TidalPy.config["numerical"][key])
        TidalPy.config["numerical"][key] = value
        update_constants()

    yield set_value
    for key, value in originals.items():
        TidalPy.config["numerical"][key] = value
    update_constants()


@pytest.mark.parametrize(
    "key, expected",
    [
        ("numerical_floor", _DEFAULT_FLOOR),
        ("layer_continuity_rtol", _DEFAULT_CONTINUITY_RTOL),
        ("max_start_radius_fraction", _DEFAULT_START_RADIUS_FRACTION),
        ("frequency_match_rtol", 1.0e-9),
        ("tides_3d_latitude_nodes", 16),
        ("tides_3d_longitude_nodes", 64),
        ("tides_3d_radial_slices", 16),
        ("minimum_nusselt", 1.0),
        ("eos_invert_rtol", 1.0e-13),
        ("eos_invert_max_iters", 60),
    ])
def test_default_is_wired_through(key, expected):
    """The packaged default is in the config and in TidalPy.constants."""
    assert TidalPy.config["numerical"][key] == expected
    assert getattr(TidalPy.constants, key) == expected


def test_rheology_guard_uses_the_configured_floor(numerical_setter):
    """At zero forcing frequency the Maxwell modulus uses the configured frequency floor."""
    maxwell = Maxwell()
    numerical_setter("numerical_floor", _DEFAULT_FLOOR)
    at_default = maxwell.calc_complex_modulus(5.0e10, 1.0e19, 0.0)
    numerical_setter("numerical_floor", _TEST_FLOOR)
    at_test = maxwell.calc_complex_modulus(5.0e10, 1.0e19, 0.0)

    # The imaginary part is viscosity * guarded_frequency, so it is the floor itself scaled.
    assert at_default.imag == pytest.approx(_DEFAULT_FLOOR, rel=1e-9)
    assert at_test.imag == pytest.approx(_TEST_FLOOR, rel=1e-9)
    assert at_test.real > at_default.real


def test_radiogenics_guard_uses_the_configured_floor(numerical_setter):
    """The average half life is raised to the configured floor before it becomes a decay constant."""
    seconds = 1.0e6
    mass = 1.0e22
    heat_production = 1.0e-11
    model = FixedRadiogenics(fixed_heat_production=heat_production,
                             average_half_life=1.0, ref_time=0.0)

    numerical_setter("numerical_floor", _DEFAULT_FLOOR)
    assert model.calc_heating(seconds, mass) == 0.0

    # A huge floor, equal to the elapsed time, is needed to see the guard: a smaller one still underflows to zero.
    numerical_setter("numerical_floor", seconds)
    assert model.calc_heating(seconds, mass) == pytest.approx(0.5 * mass * heat_production, rel=1e-9)

    # A zero half life short-circuits to the undecayed rate, independent of the floor.
    undecayed = FixedRadiogenics(fixed_heat_production=heat_production,
                                 average_half_life=0.0, ref_time=0.0)
    assert undecayed.calc_heating(seconds, mass) == pytest.approx(mass * heat_production, rel=1e-12)


def test_cooling_guard_uses_the_configured_floor(numerical_setter):
    """A zero thickness layer's conductive flux divides by the configured floor."""
    model = ConductiveCooling()
    delta_temp, thermal_conductivity = 100.0, 3.0
    arguments = (delta_temp, 0.0, 9.8, 3500.0, 1.0e21, thermal_conductivity, 1.0e-6, 3.0e-5)

    numerical_setter("numerical_floor", _DEFAULT_FLOOR)
    at_default = model.calc_cooling(*arguments).cooling_flux
    thickness_floor = 1.0e3
    numerical_setter("numerical_floor", thickness_floor)
    at_test = model.calc_cooling(*arguments).cooling_flux

    assert at_default == pytest.approx(thermal_conductivity * delta_temp / _DEFAULT_FLOOR, rel=1e-9)
    assert at_test == pytest.approx(thermal_conductivity * delta_temp / thickness_floor, rel=1e-9)


def test_restoring_the_floor_restores_the_result(numerical_setter):
    """Setting the floor back restores the original result."""
    maxwell = Maxwell()
    numerical_setter("numerical_floor", _DEFAULT_FLOOR)
    before = maxwell.calc_complex_modulus(5.0e10, 1.0e19, 0.0)
    numerical_setter("numerical_floor", _TEST_FLOOR)
    numerical_setter("numerical_floor", _DEFAULT_FLOOR)
    assert maxwell.calc_complex_modulus(5.0e10, 1.0e19, 0.0) == before


def _two_layer_world(gap):
    """A world whose outer layer starts `gap` meters above the inner layer's 1e6 m outer radius."""
    world = LayeredWorld("continuity", 2.0e6, 1.0e23)
    world.add_layer(BaseLayer("inner", 0, 0.0, 1.0e6, 5.0e22))
    return world, BaseLayer("outer", 1, 1.0e6 + gap, 2.0e6, 5.0e22)


def test_continuity_rtol_decides_whether_a_gap_is_accepted(numerical_setter):
    """A 10 m gap on a 1e6 m boundary (1e-5 relative) is refused at the default and accepted at 1e-4."""
    numerical_setter("layer_continuity_rtol", _DEFAULT_CONTINUITY_RTOL)
    world, outer = _two_layer_world(10.0)
    with pytest.raises(ValueError, match="not continuous"):
        world.add_layer(outer)

    numerical_setter("layer_continuity_rtol", 1.0e-4)
    world, outer = _two_layer_world(10.0)
    world.add_layer(outer)
    assert world.num_layers == 2


def test_continuous_geometry_is_accepted_at_the_default(numerical_setter):
    """Continuous layers are accepted at the default tolerance."""
    numerical_setter("layer_continuity_rtol", _DEFAULT_CONTINUITY_RTOL)
    world, outer = _two_layer_world(0.0)
    world.add_layer(outer)
    assert world.num_layers == 2


def _homogeneous_solve(starting_radius):
    """One layer supplied-moduli solve, which is where the starting radius rules are applied."""
    import numpy as np

    slices = 30
    frequency = 2.0 * math.pi / 86400.0
    radius = np.linspace(0.0, _PLANET_RADIUS, slices)
    density = np.full(slices, 5000.0)
    viscosity = np.full(slices, 1.0e19)
    complex_shear = Maxwell().calc_complex_modulus_vectorize_modulus(
        np.full(slices, 5.0e10), viscosity, frequency)
    complex_bulk = Elastic().calc_complex_modulus_vectorize_modulus(
        np.full(slices, 1.0e11), viscosity, frequency)
    return radial_solver(
        radius.copy(),
        density.copy(),
        complex_bulk.copy(),
        complex_shear.copy(),
        frequency,
        5000.0,
        ("solid",),
        (False,),
        (False,),
        np.asarray((_PLANET_RADIUS,)),
        degree_l=2,
        solve_for=("tidal",),
        starting_radius=starting_radius,
        nondimensionalize=True,
        integration_method="DOP853",
        integration_rtol=1e-8,
        integration_atol=1e-10,
        max_num_steps=5_000_000,
        raise_on_fail=False)


def test_start_radius_just_below_the_fraction_is_accepted(numerical_setter):
    """A starting radius just below the default fraction solves."""
    numerical_setter("max_start_radius_fraction", _DEFAULT_START_RADIUS_FRACTION)
    assert _homogeneous_solve(0.89 * _PLANET_RADIUS).success


def test_start_radius_above_the_fraction_is_rejected(numerical_setter):
    """A starting radius above the fraction raises, reporting the configured fraction."""
    numerical_setter("max_start_radius_fraction", _DEFAULT_START_RADIUS_FRACTION)
    with pytest.raises(ValueError, match=r"above 90% of the planet radius"):
        _homogeneous_solve(0.91 * _PLANET_RADIUS)


def test_start_radius_fraction_is_configurable(numerical_setter):
    """A loosened fraction accepts a starting radius the default refuses."""
    numerical_setter("max_start_radius_fraction", 0.95)
    assert _homogeneous_solve(0.91 * _PLANET_RADIUS).success
    with pytest.raises(ValueError, match=r"above 95% of the planet radius"):
        _homogeneous_solve(0.96 * _PLANET_RADIUS)


def test_frequency_match_rtol_decides_which_modes_share_a_frequency(numerical_setter):
    """Tidal modes differing by a part in 1e-6 stay distinct at the default and merge at a looser tolerance."""
    from TidalPy.Tides.potential import global_potential
    n = 2.0e-5
    # Spin at n (1 + 1e-6): modes such as (2, 2, 0, 1) at 3n - 2 spin and (2, 0, 1, 1) at n differ by a part in 1e-6.
    args = (1.0e6, n, n * (1.0 + 1.0e-6), 0.05, 0.0, 4.0e8, 1.0e27, 6.674e-11, 2, 2, 2, "off")

    numerical_setter("frequency_match_rtol", 1.0e-9)
    unique_at_default = len(global_potential(*args)[2])
    numerical_setter("frequency_match_rtol", 1.0e-4)
    unique_when_loose = len(global_potential(*args)[2])
    assert unique_when_loose < unique_at_default


def test_3d_quadrature_setting_reaches_the_heating_integral(numerical_setter):
    """A coarse configured radial quadrature changes the 3D heating total exactly as the same argument does."""
    from TidalPy.constants import G
    from TidalPy.Structures import build_world
    world = build_world("io")
    world.solve_eos(G_to_use=G)
    orbit = dict(orbital_frequency=4.11e-5, spin_frequency=4.11e-5, eccentricity=0.0041, obliquity=0.0,
                 semi_major_axis=4.217e8, host_mass=1.898e27)

    def total(**kwargs):
        return float(world.calc_3d_tides(
            radial_summed=True, latitude_summed=True, longitude_summed=True, **orbit, **kwargs)["total"])

    fine = total()
    numerical_setter("tides_3d_radial_slices", 2)
    coarse_by_config = total()
    numerical_setter("tides_3d_radial_slices", 16)
    coarse_by_argument = total(radial_slices=2)
    assert math.isclose(coarse_by_config, coarse_by_argument, rel_tol=1.0e-12)
    assert not math.isclose(coarse_by_config, fine, rel_tol=1.0e-6)
    assert math.isclose(total(), fine, rel_tol=1.0e-12)


def test_minimum_nusselt_floors_the_convection_model(numerical_setter):
    """A sub-critical convective layer sits on the configured Nusselt floor."""
    from TidalPy.Cooling.cooling import ConvectiveCooling
    # A stiff, thin layer: sub-critical, so the Nusselt number sits on the floor.
    inputs = (100.0, 1.0e4, 9.8, 3300.0, 1.0e24, 4.0, 1.0e-6, 3.0e-5)
    numerical_setter("minimum_nusselt", 2.0)
    at_default = ConvectiveCooling().calc_cooling(*inputs)
    assert at_default.nusselt == 2.0
    numerical_setter("minimum_nusselt", 3.0)
    at_three = ConvectiveCooling().calc_cooling(*inputs)
    assert at_three.nusselt == 3.0
    assert at_three.boundary_layer_thickness == pytest.approx(1.0e4 / 3.0)
    assert at_three.cooling_flux == pytest.approx(1.5 * at_default.cooling_flux)


@pytest.mark.parametrize("model_name", ("birch_murnaghan", "vinet"))
def test_eos_inversion_default_comes_from_the_config(numerical_setter, model_name):
    """Every way of building a compressible EOS takes and stores the configured inversion settings."""
    from TidalPy.Material.eos.material_eos import BirchMurnaghanEOS, VinetEOS, make_material_eos
    eos_class = BirchMurnaghanEOS if model_name == "birch_murnaghan" else VinetEOS
    numerical_setter("eos_invert_rtol", 1.0e-9)
    numerical_setter("eos_invert_max_iters", 25)
    built = (eos_class(3500.0, 1.3e11, 4.5),
             make_material_eos(model_name, {"reference_density_kg_m3": 3500.0}),
             make_material_eos(model_name))
    for eos in built:
        assert (eos.invert_rtol, eos.invert_max_iters) == (1.0e-9, 25)
        config = eos.get_config_dict()
        assert (config["invert_rtol"], config["invert_max_iters"]) == (1.0e-9, 25)


def test_eos_inversion_is_fixed_when_the_model_is_built(numerical_setter):
    """A later config change reaches new EOS models only, and an explicit value always wins."""
    from TidalPy.Material.eos.material_eos import BirchMurnaghanEOS
    numerical_setter("eos_invert_rtol", 1.0e-9)
    existing = BirchMurnaghanEOS(3500.0, 1.3e11, 4.5)
    explicit = BirchMurnaghanEOS(3500.0, 1.3e11, 4.5, invert_rtol=1.0e-11, invert_max_iters=80)
    numerical_setter("eos_invert_rtol", 1.0e-7)
    assert existing.invert_rtol == 1.0e-9
    assert BirchMurnaghanEOS(3500.0, 1.3e11, 4.5).invert_rtol == 1.0e-7
    assert (explicit.invert_rtol, explicit.invert_max_iters) == (1.0e-11, 80)


def test_eos_inversion_tolerance_changes_the_density(numerical_setter):
    """A one step, loose inversion returns a different density than a tight one."""
    from TidalPy.Material.eos.material_eos import BirchMurnaghanEOS
    pressure = 5.0e10
    numerical_setter("eos_invert_rtol", 1.0e-13)
    tight = BirchMurnaghanEOS(3500.0, 1.3e11, 4.5)
    numerical_setter("eos_invert_rtol", 1.0e-1)
    numerical_setter("eos_invert_max_iters", 1)
    loose = BirchMurnaghanEOS(3500.0, 1.3e11, 4.5)
    assert tight.calc_density(pressure, 300.0) != loose.calc_density(pressure, 300.0)
