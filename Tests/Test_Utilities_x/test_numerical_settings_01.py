"""The configured ``[numerical]`` settings reach the C++ code that uses them.

Values that were `constexpr` literals inside the C++ headers now live in the configuration and reach C++
through the shared config singleton. These tests pin that wiring: changing a value and calling
``update_constants_x`` has to change observable behavior, otherwise the C++ side is still reading a literal.

- ``numerical_floor``: the magnitude a guarded denominator is raised to, in rheology_x, cooling_x, and
  radiogenics_x.
- ``layer_continuity_rtol``: how far a layer's inner radius may sit from the previous layer's outer radius.
- ``max_start_radius_fraction``: how close to the surface a radial-solver integration may start.
"""
import math

import pytest

import TidalPy
from TidalPy.constants import update_constants_x
from TidalPy.rheology_x import Maxwell
from TidalPy.cooling_x.cooling import ConductiveCooling
from TidalPy.radiogenics_x.radiogenics import FixedRadiogenics
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.rheology_x import Elastic
from TidalPy.RadialSolver_x.solver import radial_solver

_DEFAULT_FLOOR = 1.0e-100
# Large enough that a guarded result moves by many orders of magnitude, small enough to stay unphysical.
_TEST_FLOOR = 1.0e-5


@pytest.fixture
def numerical_setter():
    """Set any ``[numerical]`` key for one test and restore every change afterwards."""
    originals = {}

    def set_value(key, value):
        originals.setdefault(key, TidalPy.config_x["numerical"][key])
        TidalPy.config_x["numerical"][key] = value
        update_constants_x()

    yield set_value
    for key, value in originals.items():
        TidalPy.config_x["numerical"][key] = value
    update_constants_x()


def test_default_floor_is_wired_through():
    from TidalPy.constants import numerical_floor
    assert TidalPy.config_x["numerical"]["numerical_floor"] == _DEFAULT_FLOOR
    assert numerical_floor == _DEFAULT_FLOOR


def test_rheology_guard_uses_the_configured_floor(numerical_setter):
    """At zero forcing frequency the Maxwell compliance divides by the guarded frequency."""
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
    """The average half life is guarded before it becomes a decay constant.

    The guard's effect is only observable where the guarded half life is comparable to the elapsed
    time: below that the decay constant is so large that ``c_safe_exp`` underflows to zero whatever
    the floor is. So this uses a deliberately huge floor (restored by the fixture) with a one-second
    half life and one megasecond elapsed, where guarding turns a total underflow into exactly half the
    undecayed rate.
    """
    seconds = 1.0e6
    mass = 1.0e22
    heat_production = 1.0e-11
    model = FixedRadiogenics(fixed_heat_production=heat_production,
                             average_half_life=1.0, ref_time=0.0)

    numerical_setter("numerical_floor", _DEFAULT_FLOOR)
    assert model.calc_heating(seconds, mass) == 0.0   # 1e6 half lives elapsed: fully decayed

    numerical_setter("numerical_floor", seconds)                             # the half life is lifted to the elapsed time
    assert model.calc_heating(seconds, mass) == pytest.approx(0.5 * mass * heat_production, rel=1e-9)

    # A zero half life still short-circuits to the undecayed rate, independent of the floor.
    undecayed = FixedRadiogenics(fixed_heat_production=heat_production,
                                 average_half_life=0.0, ref_time=0.0)
    assert undecayed.calc_heating(seconds, mass) == pytest.approx(mass * heat_production, rel=1e-12)


def test_cooling_guard_uses_the_configured_floor(numerical_setter):
    """Conductive flux divides by the layer thickness, so a zero-thickness layer hits the guard."""
    model = ConductiveCooling()
    delta_temp, thermal_conductivity = 100.0, 3.0
    arguments = (delta_temp, 0.0, 9.8, 3500.0, 1.0e21, thermal_conductivity, 1.0e-6, 3.0e-5)

    numerical_setter("numerical_floor", _DEFAULT_FLOOR)
    at_default = model.calc_cooling(*arguments).cooling_flux
    thickness_floor = 1.0e3
    numerical_setter("numerical_floor", thickness_floor)
    at_test = model.calc_cooling(*arguments).cooling_flux

    # flux = k * dT / guarded_thickness, so the guarded thickness is the floor exactly.
    assert at_default == pytest.approx(thermal_conductivity * delta_temp / _DEFAULT_FLOOR, rel=1e-9)
    assert at_test == pytest.approx(thermal_conductivity * delta_temp / thickness_floor, rel=1e-9)


def test_restoring_the_floor_restores_the_result(numerical_setter):
    """update_constants_x is repeatable, so a session can change the floor back."""
    maxwell = Maxwell()
    numerical_setter("numerical_floor", _DEFAULT_FLOOR)
    before = maxwell.calc_complex_modulus(5.0e10, 1.0e19, 0.0)
    numerical_setter("numerical_floor", _TEST_FLOOR)
    numerical_setter("numerical_floor", _DEFAULT_FLOOR)
    assert maxwell.calc_complex_modulus(5.0e10, 1.0e19, 0.0) == before


# =====================================================================================================================
# layer_continuity_rtol
# =====================================================================================================================
_DEFAULT_CONTINUITY_RTOL = 1.0e-6


def _two_layer_world(gap):
    """A world whose outer layer starts `gap` metres above the inner layer's 1e6 m outer radius."""
    world = LayeredWorld("continuity", 2.0e6, 1.0e23)
    world.add_layer(PhysicsLayer("inner", 0, 0.0, 1.0e6, 5.0e22))
    return world, PhysicsLayer("outer", 1, 1.0e6 + gap, 2.0e6, 5.0e22)


def test_default_continuity_rtol_is_wired_through():
    from TidalPy.constants import layer_continuity_rtol
    assert TidalPy.config_x["numerical"]["layer_continuity_rtol"] == _DEFAULT_CONTINUITY_RTOL
    assert layer_continuity_rtol == _DEFAULT_CONTINUITY_RTOL


def test_continuity_rtol_decides_whether_a_gap_is_accepted(numerical_setter):
    """A 10 m gap on a 1e6 m boundary is 1e-5 relative: outside the default, inside a loosened one."""
    numerical_setter("layer_continuity_rtol", _DEFAULT_CONTINUITY_RTOL)
    world, outer = _two_layer_world(10.0)
    with pytest.raises(ValueError, match="not continuous"):
        world.add_layer(outer)

    numerical_setter("layer_continuity_rtol", 1.0e-4)
    world, outer = _two_layer_world(10.0)
    world.add_layer(outer)
    assert world.num_layers == 2


def test_continuous_geometry_is_accepted_at_the_default(numerical_setter):
    numerical_setter("layer_continuity_rtol", _DEFAULT_CONTINUITY_RTOL)
    world, outer = _two_layer_world(0.0)
    world.add_layer(outer)
    assert world.num_layers == 2


# =====================================================================================================================
# max_start_radius_fraction
# =====================================================================================================================
_DEFAULT_START_RADIUS_FRACTION = 0.90
_PLANET_RADIUS = 1.0e6


def _homogeneous_solve(starting_radius):
    """One-layer supplied-moduli solve, which is where the starting-radius rules are applied."""
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
        radius.copy(), density.copy(), complex_bulk.copy(), complex_shear.copy(),
        frequency, 5000.0, ("solid",), (False,), (False,), np.asarray((_PLANET_RADIUS,)),
        degree_l=2, solve_for=("tidal",), starting_radius=starting_radius,
        nondimensionalize=True, integration_method="DOP853", integration_rtol=1e-8,
        integration_atol=1e-10, max_num_steps=5_000_000, raise_on_fail=False)


def test_default_start_radius_fraction_is_wired_through():
    from TidalPy.constants import max_start_radius_fraction
    assert TidalPy.config_x["numerical"]["max_start_radius_fraction"] == _DEFAULT_START_RADIUS_FRACTION
    assert max_start_radius_fraction == _DEFAULT_START_RADIUS_FRACTION


def test_start_radius_just_below_the_fraction_is_accepted(numerical_setter):
    numerical_setter("max_start_radius_fraction", _DEFAULT_START_RADIUS_FRACTION)
    assert _homogeneous_solve(0.89 * _PLANET_RADIUS).success


def test_start_radius_above_the_fraction_is_rejected(numerical_setter):
    """The message reports the configured fraction rather than a hard-coded 90%."""
    numerical_setter("max_start_radius_fraction", _DEFAULT_START_RADIUS_FRACTION)
    with pytest.raises(ValueError, match=r"above 90% of the planet radius"):
        _homogeneous_solve(0.91 * _PLANET_RADIUS)


def test_start_radius_fraction_is_configurable(numerical_setter):
    """Loosening the fraction accepts a starting radius the default refuses."""
    numerical_setter("max_start_radius_fraction", 0.95)
    assert _homogeneous_solve(0.91 * _PLANET_RADIUS).success
    with pytest.raises(ValueError, match=r"above 95% of the planet radius"):
        _homogeneous_solve(0.96 * _PLANET_RADIUS)


# =====================================================================================================================
# frequency_match_rtol: which tidal modes share one frequency
# =====================================================================================================================
def test_default_frequency_match_rtol_is_wired_through():
    from TidalPy.constants import frequency_match_rtol
    assert TidalPy.config_x["numerical"]["frequency_match_rtol"] == 1.0e-9
    assert frequency_match_rtol == 1.0e-9


def test_frequency_match_rtol_decides_which_modes_share_a_frequency(numerical_setter):
    """Two modes whose frequencies differ by a part in 1e-6 are distinct at the default and one at a looser value."""
    from TidalPy.Tides_x.potential import global_potential
    n = 2.0e-5
    # Spin at n (1 + 1e-6): modes such as (2, 2, 0, 1) at 3n - 2 spin and (2, 0, 1, 1) at n differ by a part in
    # 1e-6, which the default tolerance keeps apart and a looser one merges.
    args = (1.0e6, n, n * (1.0 + 1.0e-6), 0.05, 0.0, 4.0e8, 1.0e27, 6.674e-11, 2, 2, 2, "off")

    numerical_setter("frequency_match_rtol", 1.0e-9)
    unique_at_default = len(global_potential(*args)[2])
    numerical_setter("frequency_match_rtol", 1.0e-4)
    unique_when_loose = len(global_potential(*args)[2])
    assert unique_when_loose < unique_at_default


# =====================================================================================================================
# tides_3d_*: the quadrature resolutions of calc_3d_tides
# =====================================================================================================================
def test_default_3d_quadrature_is_wired_through():
    from TidalPy.constants import tides_3d_latitude_nodes, tides_3d_longitude_nodes, tides_3d_radial_slices
    numerical = TidalPy.config_x["numerical"]
    assert (numerical["tides_3d_latitude_nodes"], numerical["tides_3d_longitude_nodes"],
            numerical["tides_3d_radial_slices"]) == (16, 64, 16)
    assert (tides_3d_latitude_nodes, tides_3d_longitude_nodes, tides_3d_radial_slices) == (16, 64, 16)


def test_3d_quadrature_setting_reaches_the_heating_integral(numerical_setter):
    """A coarse radial quadrature moves the collapsed total; the same value passed as an argument moves it the same."""
    import math
    from TidalPy.constants import G
    from TidalPy.structures_x import build_world
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


# =====================================================================================================================
# minimum_nusselt: the floor of the convection model
# =====================================================================================================================
def test_default_minimum_nusselt_is_wired_through():
    from TidalPy.constants import minimum_nusselt
    assert TidalPy.config_x["numerical"]["minimum_nusselt"] == 2.0
    assert minimum_nusselt == 2.0


def test_minimum_nusselt_floors_the_convection_model(numerical_setter):
    from TidalPy.cooling_x.cooling import ConvectiveCooling
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
