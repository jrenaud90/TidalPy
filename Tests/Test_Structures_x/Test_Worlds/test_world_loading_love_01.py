"""Load Love numbers through the tidal machinery.

Tidal and load Love numbers differ only in the surface boundary condition the radial solver applies, so
every tide entry point takes a ``loading`` flag that is passed straight down to that stage; heating,
potential derivatives, and the 3D fields are then computed from whichever set was requested by the same
code. These tests pin the three things that could rot: that the flag actually reaches the solver and
changes the answer, that it agrees with the pre-existing ``solve_love_numbers(solve_for='loading')``
route, and that every path which cannot produce load Love numbers refuses instead of silently ignoring
the request.
"""
import math

import numpy as np
import pytest

from TidalPy.structures_x.configs import build_world
from TidalPy.Tides_x.classes import make_tide


_HOST_MASS = 1.898e27
_SEMI_MAJOR_AXIS = 4.217e8
_ECCENTRICITY = 0.0041
_ORBITAL_FREQUENCY = math.sqrt(6.674e-11 * _HOST_MASS / _SEMI_MAJOR_AXIS ** 3)
_ORBIT = (_ORBITAL_FREQUENCY, _ORBITAL_FREQUENCY, _ECCENTRICITY, 0.0, _SEMI_MAJOR_AXIS, _HOST_MASS)


def _layered_world():
    """Io, whose bundled config already selects the rheology tide model and the radial solver."""
    world = build_world("io")
    world.solve_eos()
    world.set_spin_frequency(_ORBITAL_FREQUENCY)
    return world


def test_loading_changes_the_love_numbers():
    """The flag has to reach the radial solver; if it were dropped the two solves would agree."""
    world = _layered_world()
    world.calc_tides(*_ORBIT)
    tidal_k = world.love_number_k
    world.calc_tides(*_ORBIT, loading=True)
    loading_k = world.love_number_k

    assert tidal_k != loading_k
    # Tidal k2 is positive with a negative imaginary part; the load response reverses both.
    assert tidal_k.real > 0.0 and tidal_k.imag < 0.0
    assert loading_k.real < 0.0 and loading_k.imag > 0.0


def test_loading_matches_the_standalone_love_solve():
    """calc_tides(loading=True) must agree with the route that existed before it."""
    world = _layered_world()
    world.calc_tides(*_ORBIT, loading=True)
    through_tides = world.love_number_k

    world.solve_love_numbers(_ORBITAL_FREQUENCY, solve_for="loading")
    assert world.love_number_k == pytest.approx(through_tides, rel=1e-12)


def test_calc_tides_loading_alias_matches_the_keyword():
    world = _layered_world()
    world.calc_tides(*_ORBIT, loading=True)
    keyword_heating = world.get_tidal_heating()

    world.calc_tides_loading(*_ORBIT)
    assert world.get_tidal_heating() == pytest.approx(keyword_heating, rel=1e-12)


def test_loading_heating_is_negative_with_this_sign_convention():
    """Documented behavior, pinned so it cannot change quietly.

    The heating expression carries -Im(k). Load Love numbers have the opposite imaginary sign to tidal
    ones, so substituting them returns a negative number whose magnitude is the dissipation rate. The
    flag deliberately reuses the tidal machinery rather than deriving a separate loading heating.
    """
    world = _layered_world()
    world.calc_tides(*_ORBIT)
    tidal_heating = world.get_tidal_heating()
    world.calc_tides(*_ORBIT, loading=True)
    loading_heating = world.get_tidal_heating()

    assert tidal_heating > 0.0
    assert loading_heating < 0.0


@pytest.mark.parametrize("method", ("homogeneous",))
def test_loading_rejects_love_methods_that_never_run_the_radial_solver(method):
    world = _layered_world()
    world.set_tide_config(love_method=method)
    with pytest.raises(RuntimeError, match="radial solver"):
        world.calc_tides(*_ORBIT, loading=True)
    # The same world still solves tidally, so the rejection is specific to the loading request.
    world.calc_tides(*_ORBIT)
    assert world.tides_solved


@pytest.mark.parametrize("model", ("cpl", "ctl", "ctl_q"))
def test_loading_rejects_analytic_tide_models(model):
    """An analytic tide model never runs a radial solve, so it would otherwise ignore the flag."""
    world = _layered_world()
    world.set_tide_model(make_tide(model))
    with pytest.raises(RuntimeError, match="analytic"):
        world.calc_tides(*_ORBIT, loading=True)
    world.calc_tides(*_ORBIT)
    assert world.tides_solved


def test_loading_rejects_a_world_without_a_radial_solver():
    """A gas giant carries an analytic tide model, so the flag has nothing to act on."""
    world = build_world("jupiter_simple")
    with pytest.raises(RuntimeError, match="analytic"):
        world.calc_tides(*_ORBIT, loading=True)
    with pytest.raises(RuntimeError, match="analytic"):
        world.calc_tides_loading(*_ORBIT)
    world.calc_tides(*_ORBIT)
    assert world.tides_solved


def test_3d_tides_accepts_loading_and_changes_the_field():
    world = _layered_world()
    radii = np.linspace(2.0e5, 1.80e6, 6)
    colatitudes = np.linspace(0.2, math.pi - 0.2, 4)
    longitudes = np.linspace(0.0, 2.0 * math.pi, 5)

    tidal = np.asarray(world.calc_3d_tides(
        *_ORBIT, radii=radii, colatitudes=colatitudes, longitudes=longitudes)["heating"])
    loading = np.asarray(world.calc_3d_tides(
        *_ORBIT, radii=radii, colatitudes=colatitudes, longitudes=longitudes, loading=True)["heating"])

    assert tidal.shape == loading.shape
    assert not np.allclose(np.nan_to_num(tidal), np.nan_to_num(loading))


def test_calc_3d_tides_loading_alias_matches_the_keyword():
    world = _layered_world()
    radii = np.linspace(2.0e5, 1.80e6, 5)
    colatitudes = np.linspace(0.3, math.pi - 0.3, 3)
    longitudes = np.linspace(0.0, 2.0 * math.pi, 4)
    grid = dict(radii=radii, colatitudes=colatitudes, longitudes=longitudes)

    keyword = np.asarray(world.calc_3d_tides(*_ORBIT, loading=True, **grid)["heating"])
    alias = np.asarray(world.calc_3d_tides_loading(*_ORBIT, **grid)["heating"])
    np.testing.assert_allclose(np.nan_to_num(alias), np.nan_to_num(keyword), rtol=1e-12)


def test_every_3d_entry_point_takes_the_loading_keyword():
    """The flag is on all of them, not just the two with alias wrappers."""
    world = _layered_world()
    radii = np.asarray([1.0e6, 1.5e6])
    colatitudes = np.asarray([0.7, 1.9])

    assert math.isfinite(world.get_3d_tidal_heating(*_ORBIT, 1.5e6, 1.0, loading=True))
    batch = np.asarray(world.get_3d_tidal_heating_array(*_ORBIT, radii, colatitudes, loading=True))
    assert batch.shape == radii.shape
    displacements = world.calc_3d_displacements(
        *_ORBIT, radii, colatitudes, np.asarray([0.0, 3.0]), np.asarray([0.0]), loading=True)
    assert "radial" in displacements or len(displacements) > 0
