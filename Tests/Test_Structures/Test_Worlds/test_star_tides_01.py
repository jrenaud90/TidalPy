"""Analytic tides (cpl, ctl) on a layerless StarWorld, and rejection of the rheology tide model there."""
import math

import pytest

from TidalPy.constants import G
from TidalPy.Structures.worlds.stellar import StarWorld
from TidalPy.Tides.classes.tide import make_tide


_STAR_RADIUS = 6.957e8     # [m] (solar)
_STAR_MASS   = 1.989e30    # [kg] (solar)

_HOST_MASS = 1.898e27
_SMA       = 1.0e10
_N         = 3.0e-6
_ECC       = 0.05


def _star(model="cpl", config=None):
    star = StarWorld("test_star", _STAR_RADIUS, _STAR_MASS)
    star.set_tide_model(make_tide(model, config if config is not None else {"fixed_k": [0.03], "fixed_q": [1.0e6]}))
    # Truncation 2 keeps the leading-order e^2 dissipation of the analytic rate below.
    star.set_tide_config(min_degree_l=2, max_degree_l=2,
                         eccentricity_truncation=2, obliquity_truncation=0)
    return star


def _solve(world):
    world.calc_tides(orbital_frequency=_N, spin_frequency=_N, eccentricity=_ECC,
                     obliquity=0.0, semi_major_axis=_SMA, host_mass=_HOST_MASS)


def test_star_tide_model_flag_and_missing_model_raises():
    star = StarWorld("s", _STAR_RADIUS, _STAR_MASS)
    assert star.tide_model_set is False
    with pytest.raises(RuntimeError):
        _solve(star)
    star.set_tide_model(make_tide("cpl", {"fixed_k": [0.03], "fixed_q": [1.0e6]}))
    assert star.tide_model_set is True


@pytest.mark.parametrize("model, config", [
    ("cpl", {"fixed_k": [0.03], "fixed_q": [1.0e6]}),
    ("ctl", {"fixed_k": [0.03], "fixed_dt_s": [10.0]}),
], ids=["cpl", "ctl"])
def test_star_analytic_tides_positive_heating(model, config):
    star = _star(model, config)
    assert not star.tides_solved
    _solve(star)
    assert star.tides_solved
    assert star.get_num_tidal_modes() > 0
    assert star.get_tidal_heating() > 0.0
    dUdM, dUdw, dUdO = star.get_tidal_potential_derivatives().values()
    assert abs(dUdM) > 0.0


def test_star_cpl_matches_analytic_rate():
    """A synchronous star at truncation 2 reproduces the analytic CPL heating rate."""
    k2, q2 = 0.03, 1.0e6
    star = _star("cpl", {"fixed_k": [k2], "fixed_q": [q2]})
    _solve(star)
    expected = (21.0 / 2.0) * (k2 / q2) * G * _HOST_MASS ** 2 * _STAR_RADIUS ** 5 * _N * _ECC ** 2 / _SMA ** 6
    assert math.isclose(star.get_tidal_heating(), expected, rel_tol=1.0e-12)


def test_star_rejects_rheology_model():
    """The rheology model needs the radial solver, which a layerless star lacks."""
    star = _star("rheology", {})
    with pytest.raises(RuntimeError):
        _solve(star)
