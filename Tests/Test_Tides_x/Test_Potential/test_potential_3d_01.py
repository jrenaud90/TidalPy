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
        colatitude=_COLAT, longitude=_LON, G_to_use=G,
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
    """Raising max_degree_l to 3 adds degree-3 modes on top of the degree-2 set."""
    degrees_l2, freqs_l2, _ = _modes(max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)
    degrees_l3, freqs_l3, _ = _modes(max_degree_l=3, eccentricity_truncation=2, obliquity_truncation=0)
    assert set(degrees_l2.tolist()) == {2}
    assert 3 in set(degrees_l3.tolist())
    assert freqs_l3.shape[0] > freqs_l2.shape[0]


def test_bad_truncation_raises():
    with pytest.raises(NotImplementedError):
        _modes(max_degree_l=2, eccentricity_truncation=99, obliquity_truncation=0)


# =====================================================================================================================
# Numeric amplitude validation vs an independent reference (Kaula eccentricity functions)
# =====================================================================================================================
# At zero obliquity each (m, q) frequency hosts exactly one active p (F_220(0) = 3 with F_221 = F_222 = 0;
# F_201(0) = -1/2 with F_200 = F_202 = 0), so a mode is identified uniquely by its signed frequency. The
# potential factor of mode (2, 2, 0, q) is linear in G_20q(e), so amplitude ratios between q values at a
# fixed point equal ratios of the eccentricity functions with every normalization and phase convention
# cancelling. Reference series (Kaula 1964; Murray & Dermott 1999):
#     G_200 = 1 - 5/2 e^2 + 13/16 e^4 - 35/288 e^6
#     G_201 = 7/2 e - 123/16 e^3 + 489/128 e^5
#     G_20-1 = -1/2 e + 1/16 e^3 - 5/384 e^5
#     G_202 = 17/2 e^2 - 115/6 e^4
_ECC_RATIO = 0.1


def _amp_at(freqs, pots, target_freq):
    """Return the unique mode potential factor at the target signed frequency."""
    hits = [i for i, f in enumerate(freqs) if abs(f - target_freq) < 1e-9 * max(abs(target_freq), _N)]
    assert len(hits) == 1, f"expected one mode at frequency {target_freq}, found {len(hits)}"
    return pots[hits[0], 0]


@pytest.mark.parametrize('q, g_ratio_tol', ((1, 1.0e-5), (-1, 1.0e-5), (2, 1.0e-3)))
def test_amplitude_ratios_match_kaula_eccentricity_functions(q, g_ratio_tol):
    """The (2, 2, 0, q)/(2, 2, 0, 0) amplitude ratios equal G_20q/G_200 from the literature series.

    The q = 2 reference series stops at e^4, so its comparison tolerance is looser.
    """
    e = _ECC_RATIO
    g_200 = 1.0 - 5.0 * e**2 / 2.0 + 13.0 * e**4 / 16.0 - 35.0 * e**6 / 288.0
    g_by_q = {
        1: 7.0 * e / 2.0 - 123.0 * e**3 / 16.0 + 489.0 * e**5 / 128.0,
        -1: -e / 2.0 + e**3 / 16.0 - 5.0 * e**5 / 384.0,
        2: 17.0 * e**2 / 2.0 - 115.0 * e**4 / 6.0,
    }
    _, freqs, pots = _modes(max_degree_l=2, eccentricity_truncation=10, obliquity_truncation=0,
                            eccentricity=e)
    amp_q0 = _amp_at(freqs, pots, 2.0 * _N - 2.0 * _SPIN)
    amp_q = _amp_at(freqs, pots, (2.0 + q) * _N - 2.0 * _SPIN)
    ratio = amp_q / amp_q0
    # The common (complex) phase factor cancels in the ratio, leaving the real G ratio.
    assert abs(ratio.imag) < 1.0e-10
    np.testing.assert_allclose(ratio.real, g_by_q[q] / g_200, rtol=g_ratio_tol)


def test_amplitude_colatitude_dependence_matches_legendre():
    """The m = 2 amplitude scales as P_22(cos theta) = 3 sin^2(theta) between two colatitudes."""
    theta_1, theta_2 = 1.1, 0.6
    _, freqs_1, pots_1 = _modes(max_degree_l=2, eccentricity_truncation=3, obliquity_truncation=0,
                                colatitude=theta_1)
    _, freqs_2, pots_2 = _modes(max_degree_l=2, eccentricity_truncation=3, obliquity_truncation=0,
                                colatitude=theta_2)
    target = 2.0 * _N - 2.0 * _SPIN
    ratio = _amp_at(freqs_1, pots_1, target) / _amp_at(freqs_2, pots_2, target)
    expected = np.sin(theta_1)**2 / np.sin(theta_2)**2
    assert abs(ratio.imag) < 1.0e-10
    np.testing.assert_allclose(ratio.real, expected, rtol=1.0e-10)


def test_amplitude_truncation_convergence():
    """The (2, 2, 0, 1) amplitude converges with rising eccentricity truncation (e = 0.1)."""
    target = 3.0 * _N - 2.0 * _SPIN

    def amp(truncation):
        _, freqs, pots = _modes(max_degree_l=2, eccentricity_truncation=truncation,
                                obliquity_truncation=0, eccentricity=_ECC_RATIO)
        return _amp_at(freqs, pots, target)

    reference = amp(20)
    errors = {trunc: abs(amp(trunc) - reference) / abs(reference) for trunc in (2, 4, 10)}
    assert errors[4] < errors[2]
    assert errors[10] < errors[4]
    assert errors[10] < 1.0e-10
