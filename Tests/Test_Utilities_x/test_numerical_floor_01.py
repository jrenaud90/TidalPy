"""The shared numerical floor (``[numerical].numerical_floor``) reaches the C++ guards.

rheology_x, cooling_x, and radiogenics_x each guard a division against a denominator that can reach
zero. The floor used to be an ``inline constexpr`` literal repeated in all three headers; it is now one
configured value on the shared C++ config singleton. These tests pin that wiring: changing the value and
calling ``update_constants_x`` has to change what the guards substitute, otherwise a module is still
reading a literal.
"""
import pytest

import TidalPy
from TidalPy.constants import update_constants_x
from TidalPy.rheology_x import Maxwell
from TidalPy.cooling_x.cooling import ConductiveCooling
from TidalPy.radiogenics_x.radiogenics import FixedRadiogenics

_DEFAULT_FLOOR = 1.0e-100
# Large enough that a guarded result moves by many orders of magnitude, small enough to stay unphysical.
_TEST_FLOOR = 1.0e-5


@pytest.fixture
def floor_setter():
    """Set ``[numerical].numerical_floor`` for one test and restore it afterwards."""
    original = TidalPy.config_x["numerical"]["numerical_floor"]

    def set_floor(value):
        TidalPy.config_x["numerical"]["numerical_floor"] = value
        update_constants_x()

    yield set_floor
    TidalPy.config_x["numerical"]["numerical_floor"] = original
    update_constants_x()


def test_default_floor_is_wired_through():
    from TidalPy.constants import numerical_floor
    assert TidalPy.config_x["numerical"]["numerical_floor"] == _DEFAULT_FLOOR
    assert numerical_floor == _DEFAULT_FLOOR


def test_rheology_guard_uses_the_configured_floor(floor_setter):
    """At zero forcing frequency the Maxwell compliance divides by the guarded frequency."""
    maxwell = Maxwell()
    floor_setter(_DEFAULT_FLOOR)
    at_default = maxwell.calc_complex_modulus(5.0e10, 1.0e19, 0.0)
    floor_setter(_TEST_FLOOR)
    at_test = maxwell.calc_complex_modulus(5.0e10, 1.0e19, 0.0)

    # The imaginary part is viscosity * guarded_frequency, so it is the floor itself scaled.
    assert at_default.imag == pytest.approx(_DEFAULT_FLOOR, rel=1e-9)
    assert at_test.imag == pytest.approx(_TEST_FLOOR, rel=1e-9)
    assert at_test.real > at_default.real


def test_radiogenics_guard_uses_the_configured_floor(floor_setter):
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

    floor_setter(_DEFAULT_FLOOR)
    assert model.calc_heating(seconds, mass) == 0.0   # 1e6 half lives elapsed: fully decayed

    floor_setter(seconds)                             # the half life is lifted to the elapsed time
    assert model.calc_heating(seconds, mass) == pytest.approx(0.5 * mass * heat_production, rel=1e-9)

    # A zero half life still short-circuits to the undecayed rate, independent of the floor.
    undecayed = FixedRadiogenics(fixed_heat_production=heat_production,
                                 average_half_life=0.0, ref_time=0.0)
    assert undecayed.calc_heating(seconds, mass) == pytest.approx(mass * heat_production, rel=1e-12)


def test_cooling_guard_uses_the_configured_floor(floor_setter):
    """Conductive flux divides by the layer thickness, so a zero-thickness layer hits the guard."""
    model = ConductiveCooling()
    delta_temp, thermal_conductivity = 100.0, 3.0
    arguments = (delta_temp, 0.0, 9.8, 3500.0, 1.0e21, thermal_conductivity, 1.0e-6, 3.0e-5)

    floor_setter(_DEFAULT_FLOOR)
    at_default = model.calc_cooling(*arguments).cooling_flux
    thickness_floor = 1.0e3
    floor_setter(thickness_floor)
    at_test = model.calc_cooling(*arguments).cooling_flux

    # flux = k * dT / guarded_thickness, so the guarded thickness is the floor exactly.
    assert at_default == pytest.approx(thermal_conductivity * delta_temp / _DEFAULT_FLOOR, rel=1e-9)
    assert at_test == pytest.approx(thermal_conductivity * delta_temp / thickness_floor, rel=1e-9)


def test_restoring_the_floor_restores_the_result(floor_setter):
    """update_constants_x is repeatable, so a session can change the floor back."""
    maxwell = Maxwell()
    floor_setter(_DEFAULT_FLOOR)
    before = maxwell.calc_complex_modulus(5.0e10, 1.0e19, 0.0)
    floor_setter(_TEST_FLOOR)
    floor_setter(_DEFAULT_FLOOR)
    assert maxwell.calc_complex_modulus(5.0e10, 1.0e19, 0.0) == before
