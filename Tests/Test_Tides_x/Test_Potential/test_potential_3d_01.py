"""Tests for the dynamic Kaula 3D tidal-potential engine (TidalPy.Tides_x.potential.potential_3d).

The engine builds every active tidal mode (l, m, p, q) from the eccentricity/obliquity functions plus
the associated Legendre functions, returning each mode's degree, signed frequency, and potential
angular factor (U + 5 derivatives) at a point. These tests check the returned shapes, the expected
mode frequencies for a known configuration, and the finiteness/consistency of the output.
"""
import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Tides_x.potential.potential_3d import tidal_potential_3d_modes


# A representative non-synchronous, eccentric, no-obliquity state.
_N = 2.0 * np.pi / 86400.0
_SPIN = 1.37 * _N            # non-commensurate spin: no mode lands exactly at zero frequency
_ECC = 0.05
_HOST = 1.0e26
_SMA = 1.0e9
_R = 6.0e6
_COLAT, _LON, _T = 1.1, 0.7, 1234.0


def _modes(**kw):
    args = dict(
        orbital_frequency=_N, spin_frequency=_SPIN, eccentricity=_ECC, obliquity=0.0,
        host_mass=_HOST, semi_major_axis=_SMA, planet_radius=_R,
        colatitude=_COLAT, longitude=_LON, time=_T, G_to_use=G,
    )
    args.update(kw)
    return tidal_potential_3d_modes(**args)


def test_shapes_consistent():
    degrees, freqs, pots = _modes(max_degree_l=2, eccentricity_truncation=3, obliquity_truncation=0)
    assert degrees.ndim == 1 and freqs.ndim == 1 and pots.ndim == 2
    assert degrees.shape[0] == freqs.shape[0] == pots.shape[0]
    assert pots.shape[1] == 6
    assert degrees.shape[0] > 0
    assert np.all(np.isfinite(freqs)) and np.all(np.isfinite(pots))
    assert np.all(degrees == 2)


def test_no_obliquity_orders_are_even():
    """At zero obliquity only m = 0 and m = 2 harmonics survive (F_2mp = 0 for odd m at i = 0), so the
    dU/dphi terms vanish for m = 0 and are nonzero for m = 2 — the mode set carries no m = 1 content."""
    degrees, freqs, pots = _modes(max_degree_l=2, eccentricity_truncation=3, obliquity_truncation=0)
    # Frequencies must be integer-ish combinations q*n and (2+q)*n - 2*spin.
    for f in freqs:
        # f = a*n + b*spin for small integers a, b -> check it matches one such combo within tol.
        found = False
        for a in range(-6, 7):
            for b in (0, -2, 2):
                if abs(f - (a * _N + b * _SPIN)) < 1e-18 + 1e-9 * abs(f):
                    found = True
        assert found, f"unexpected mode frequency {f}"


def test_obliquity_activates_m1_modes():
    """A nonzero obliquity turns on the m = 1 (P_21) harmonics, adding modes vs the no-obliquity set."""
    _, freqs_no_obl, _ = _modes(max_degree_l=2, eccentricity_truncation=3, obliquity_truncation=0)
    _, freqs_obl, _ = _modes(max_degree_l=2, eccentricity_truncation=3, obliquity_truncation=2,
                             obliquity=0.2)
    assert freqs_obl.shape[0] > freqs_no_obl.shape[0]


def test_higher_degree_adds_modes():
    _, freqs_l2, degrees_l2 = _modes(max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)
    degrees_l3, freqs_l3, _ = _modes(max_degree_l=3, eccentricity_truncation=2, obliquity_truncation=0)
    assert 3 in set(degrees_l3.tolist())
    assert freqs_l3.shape[0] >= degrees_l2.shape[0]


def test_bad_truncation_raises():
    with pytest.raises(ValueError):
        _modes(max_degree_l=2, eccentricity_truncation=99, obliquity_truncation=0)
