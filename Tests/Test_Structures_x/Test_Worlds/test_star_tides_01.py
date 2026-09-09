"""Tests for tides on a StarWorld (a layerless world).

The analytic tide pipeline (cpl/ctl/ctl_q) is common to every world type: it lives on
``c_BaseWorld``, so even a star (which has no layers and no EOS) can dissipate tidally. The
rheology model needs the radial solver and a layered interior, so it must be rejected on a star.
"""
import math

import pytest

from TidalPy.constants import G
from TidalPy.structures_x.worlds.stellar import StarWorld
from TidalPy.Tides_x.classes.tide import make_tide


_STAR_RADIUS = 6.957e8     # [m] (solar)
_STAR_MASS   = 1.989e30    # [kg] (solar)

# A close, eccentric companion forcing the star.
_HOST_MASS = 1.898e27
_SMA       = 1.0e10
_N         = 3.0e-6
_ECC       = 0.05


def _star(model="cpl", config=None):
    star = StarWorld("test_star", _STAR_RADIUS, _STAR_MASS)
    star.set_tide_model(make_tide(model, config if config is not None else {"fixed_k": [0.03], "fixed_q": [1.0e6]}))
    star.set_tide_config(min_degree_l=2, max_degree_l=2,
                         eccentricity_truncation=2, obliquity_truncation=0)
    return star


def _solve(world):
    world.calc_tides(orbital_frequency=_N, spin_frequency=_N, eccentricity=_ECC,
                     obliquity=0.0, semi_major_axis=_SMA, host_mass=_HOST_MASS)


def test_star_has_tide_api():
    star = StarWorld("s", _STAR_RADIUS, _STAR_MASS)
    assert star.tide_model_set is False
    star.set_tide_model(make_tide("cpl", {"fixed_k": [0.03], "fixed_q": [1.0e6]}))
    assert star.tide_model_set is True


def test_star_cpl_tides_positive_heating():
    star = _star("cpl")
    assert not star.tides_solved
    _solve(star)
    assert star.tides_solved
    assert star.get_num_tidal_modes() > 0
    assert star.get_tidal_heating() > 0.0


def test_star_cpl_matches_analytic_rate():
    """A synchronous, low-e star reproduces the analytic CPL heating rate."""
    k2, q2 = 0.03, 1.0e6
    star = _star("cpl", {"fixed_k": [k2], "fixed_q": [q2]})
    _solve(star)
    expected = (21.0 / 2.0) * (k2 / q2) * G * _HOST_MASS ** 2 * _STAR_RADIUS ** 5 * _N * _ECC ** 2 / _SMA ** 6
    assert math.isclose(star.get_tidal_heating(), expected, rel_tol=5.0e-3)


def test_star_potential_derivatives_present():
    star = _star("cpl")
    _solve(star)
    dUdM, dUdw, dUdO = star.get_tidal_potential_derivatives()
    assert abs(dUdM) > 0.0


def test_star_ctl_positive_heating():
    star = _star("ctl", {"fixed_k": [0.03], "fixed_dt": [10.0]})
    _solve(star)
    assert star.get_tidal_heating() > 0.0


def test_star_rejects_rheology_model():
    """The rheology model needs the radial solver; a layerless star must reject it."""
    star = _star("rheology")
    with pytest.raises(RuntimeError):
        _solve(star)


def test_star_calc_tides_without_model_raises():
    star = StarWorld("bare", _STAR_RADIUS, _STAR_MASS)
    assert not star.tide_model_set
    with pytest.raises(RuntimeError):
        _solve(star)
