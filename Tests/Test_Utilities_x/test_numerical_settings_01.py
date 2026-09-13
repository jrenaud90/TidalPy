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
