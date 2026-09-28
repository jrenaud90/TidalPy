"""A new solve_eos retires Love numbers, tidal results, and layer heating until they are solved again."""
import cmath
import math

import pytest

from TidalPy.Structures import build_world


# Io-Jupiter-like orbital state.
_STATE = dict(
    orbital_frequency=4.11e-5, spin_frequency=4.11e-5, eccentricity=0.0041, obliquity=0.0,
    semi_major_axis=4.217e8, host_mass=1.898e27)


def _is_nan(value: complex) -> bool:
    return cmath.isnan(complex(value))


@pytest.fixture()
def solved_io():
    world = build_world("io")
    world.solve_eos()
    return world


def test_love_numbers_are_retired_by_an_eos_resolve(solved_io):
    """Radial Love results read NaN after an EOS re-solve and return unchanged with the next Love solve."""
    world = solved_io
    world.solve_love_numbers(frequency=4.11e-5)
    assert world.love_solved and world.love_success
    k2 = world.love_number_k
    assert math.isfinite(k2.real)
    assert math.isfinite(world.get_love_radial_y(0.9 * world.radius).real)

    world.solve_eos()
    assert not world.love_solved
    assert not world.love_success
    for value in (world.love_number_k, world.love_number_h, world.love_number_l,
                  world.get_love_number_k(0), world.get_love_surface_y(0, 0),
                  world.get_love_radial_y(0.9 * world.radius)):
        assert _is_nan(value)

    world.solve_love_numbers(frequency=4.11e-5)
    assert world.love_solved
    assert world.love_number_k == pytest.approx(k2, rel=1.0e-9)


def test_analytic_love_numbers_are_retired_too(solved_io):
    """Homogeneous Love results are retired by an EOS re-solve."""
    world = solved_io
    world.solve_love_numbers(frequency=4.11e-5, love_method="homogeneous")
    assert world.love_solved and math.isfinite(world.love_number_k.real)
    world.solve_eos()
    assert not world.love_solved
    assert not world.love_success
    assert _is_nan(world.love_number_k)


def test_tides_are_retired_by_an_eos_resolve(solved_io):
    """Tidal results and layer heating read NaN after an EOS re-solve and return with the next calc_tides."""
    world = solved_io
    world.calc_tides(**_STATE)
    assert world.tides_solved
    heating = world.get_tidal_heating()
    assert heating > 0.0
    layer_heating = [layer.get_tidal_heating() for layer in world]
    assert any(value > 0.0 for value in layer_heating)

    world.solve_eos()
    assert not world.tides_solved
    assert math.isnan(world.get_tidal_heating())
    assert all(math.isnan(value) for value in world.get_tidal_potential_derivatives())
    assert world.get_num_tidal_modes() == 0
    assert _is_nan(world.get_tidal_love_k(2, 2, 0, 1))
    for index, layer in enumerate(world):
        assert math.isnan(world.get_layer_tidal_heating(index))
        assert math.isnan(layer.get_tidal_heating())

    world.calc_tides(**_STATE)
    assert world.tides_solved
    assert world.get_tidal_heating() == pytest.approx(heating, rel=1.0e-9)


def test_spin_and_obliquity_setters_leave_a_tidal_result_alone(solved_io):
    """Setting spin or obliquity does not retire a tidal result (calc_tides takes its own)."""
    world = solved_io
    world.calc_tides(**_STATE)
    heating = world.get_tidal_heating()
    world.set_spin_frequency(2.0 * _STATE["spin_frequency"])
    world.set_obliquity(0.1)
    assert world.tides_solved
    assert world.get_tidal_heating() == heating
