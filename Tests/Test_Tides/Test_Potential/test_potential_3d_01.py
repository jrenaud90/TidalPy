"""The 3D tidal-potential mode engine ``tidal_potential_3d_modes``: shapes, frequencies, and amplitudes."""
import numpy as np
import pytest
from scipy.special import eval_legendre

from TidalPy.constants import G
from TidalPy.Tides.potential.potential_3d import tidal_potential_3d_modes

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
    """Degrees, frequencies, and six-column potential rows share a length and are finite."""
    degrees, freqs, pots = _modes(max_degree_l=2, eccentricity_truncation=6, obliquity_truncation=0)
    assert degrees.ndim == 1 and freqs.ndim == 1 and pots.ndim == 2
    assert degrees.shape[0] == freqs.shape[0] == pots.shape[0]
    assert pots.shape[1] == 6
    assert degrees.shape[0] > 0
    assert np.all(np.isfinite(freqs)) and np.all(np.isfinite(pots))
    assert np.all(degrees == 2)


def test_no_obliquity_orders_are_even():
    """At zero obliquity every frequency is a n + b spin with b in (0, +-2): no m = 1 modes."""
    degrees, freqs, pots = _modes(max_degree_l=2, eccentricity_truncation=6, obliquity_truncation=0)
    for f in freqs:
        found = False
        for a in range(-8, 9):
            for b in (0, -2, 2):
                if abs(f - (a * _N + b * _SPIN)) < 1e-18 + 1e-9 * abs(f):
                    found = True
        assert found, f"unexpected mode frequency {f}"


def test_obliquity_activates_m1_modes():
    """A nonzero obliquity adds the m = 1 modes."""
    _, freqs_no_obl, _ = _modes(max_degree_l=2, eccentricity_truncation=6, obliquity_truncation=0)
    _, freqs_obl, _ = _modes(max_degree_l=2, eccentricity_truncation=6, obliquity_truncation=2,
                             obliquity=0.2)
    assert freqs_obl.shape[0] > freqs_no_obl.shape[0]


def test_higher_degree_adds_modes():
    """Raising max_degree_l to 3 adds degree-3 modes to the degree-2 set."""
    degrees_l2, freqs_l2, _ = _modes(max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)
    degrees_l3, freqs_l3, _ = _modes(max_degree_l=3, eccentricity_truncation=2, obliquity_truncation=0)
    assert set(degrees_l2.tolist()) == {2}
    assert 3 in set(degrees_l3.tolist())
    assert freqs_l3.shape[0] > freqs_l2.shape[0]


def test_bad_truncation_raises():
    with pytest.raises(NotImplementedError):
        _modes(max_degree_l=2, eccentricity_truncation=99, obliquity_truncation=0)


# At zero obliquity each (m, q) frequency hosts one active p, so a signed frequency identifies a mode, and the
# (2, 2, 0, q) amplitude is linear in G_20q(e): amplitude ratios are G ratios with every convention cancelling.
_ECC_RATIO = 0.1


def _amp_at(freqs, pots, target_freq):
    """The unique mode potential factor at the target signed frequency."""
    hits = [i for i, f in enumerate(freqs) if abs(f - target_freq) < 1e-9 * max(abs(target_freq), _N)]
    assert len(hits) == 1, f"expected one mode at frequency {target_freq}, found {len(hits)}"
    return pots[hits[0], 0]


# The q = 2 reference series stops at e^4, so its tolerance is looser.
@pytest.mark.parametrize('q, g_ratio_tol', ((1, 1.0e-5), (-1, 1.0e-5), (2, 1.0e-3)))
def test_amplitude_ratios_match_kaula_eccentricity_functions(q, g_ratio_tol):
    """(2, 2, 0, q) / (2, 2, 0, 0) amplitude ratios equal G_20q / G_200 (Kaula 1964; Murray & Dermott 1999)."""
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
    assert abs(ratio.imag) < 1.0e-10
    np.testing.assert_allclose(ratio.real, g_by_q[q] / g_200, rtol=g_ratio_tol)


def test_amplitude_colatitude_dependence_matches_legendre():
    """The m = 2 amplitude scales as P_22(cos theta) = 3 sin^2(theta)."""
    theta_1, theta_2 = 1.1, 0.6
    _, freqs_1, pots_1 = _modes(max_degree_l=2, eccentricity_truncation=6, obliquity_truncation=0,
                                colatitude=theta_1)
    _, freqs_2, pots_2 = _modes(max_degree_l=2, eccentricity_truncation=6, obliquity_truncation=0,
                                colatitude=theta_2)
    target = 2.0 * _N - 2.0 * _SPIN
    ratio = _amp_at(freqs_1, pots_1, target) / _amp_at(freqs_2, pots_2, target)
    expected = np.sin(theta_1)**2 / np.sin(theta_2)**2
    assert abs(ratio.imag) < 1.0e-10
    np.testing.assert_allclose(ratio.real, expected, rtol=1.0e-10)


def test_amplitude_truncation_convergence():
    """The (2, 2, 0, 1) amplitude converges with rising eccentricity truncation."""
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


def _direct_potential(
        degree_l,
        colatitude,
        longitude,
        time,
        eccentricity,
        obliquity):
    """(G M / r) (R / r)^l P_l(cos psi), psi the angle between the point and the host."""
    # Orbit: ascending node on the body's x axis, periapse at the node, tilted by the obliquity, periapse at t = 0.
    # The body spins prograde about z, and longitude runs with the rotation from the x axis at t = 0.
    mean_anomaly = _N * time
    ecc_anomaly = mean_anomaly
    for _ in range(60):
        ecc_anomaly = mean_anomaly + eccentricity * np.sin(ecc_anomaly)
    distance = _SMA * (1.0 - eccentricity * np.cos(ecc_anomaly))
    true_anomaly = 2.0 * np.arctan2(np.sqrt(1.0 + eccentricity) * np.sin(ecc_anomaly / 2.0),
                                    np.sqrt(1.0 - eccentricity) * np.cos(ecc_anomaly / 2.0))
    host = np.array([np.cos(true_anomaly),
                     np.sin(true_anomaly) * np.cos(obliquity),
                     np.sin(true_anomaly) * np.sin(obliquity)])
    rotation = _SPIN * time
    host_body = np.array([host[0] * np.cos(rotation) + host[1] * np.sin(rotation),
                          -host[0] * np.sin(rotation) + host[1] * np.cos(rotation),
                          host[2]])
    point = np.array([np.sin(colatitude) * np.cos(longitude),
                      np.sin(colatitude) * np.sin(longitude),
                      np.cos(colatitude)])
    cos_psi = float(point @ host_body)
    return G * _HOST / distance * (_R / distance)**degree_l * eval_legendre(degree_l, cos_psi)


@pytest.mark.parametrize('degree_l', (2, 3, 4, 5))
def test_summed_modes_match_direct_potential(degree_l):
    """Re[sum U_c e^{i omega t}] over one degree's modes equals the host's point-mass potential of that degree."""
    rng = np.random.default_rng(degree_l)
    scale = G * _HOST / _SMA * (_R / _SMA)**degree_l
    for _ in range(10):
        colatitude = rng.uniform(0.1, np.pi - 0.1)
        longitude = rng.uniform(0.0, 2.0 * np.pi)
        time = rng.uniform(0.0, 3.0e5)
        eccentricity = rng.uniform(0.0, 0.05)
        obliquity = rng.uniform(0.0, 1.0)
        _, freqs, pots = _modes(
            colatitude=colatitude,
            longitude=longitude,
            eccentricity=eccentricity,
            obliquity=obliquity,
            min_degree_l=degree_l,
            max_degree_l=degree_l,
            eccentricity_truncation=20,
            obliquity_truncation="gen")
        summed = float(np.sum(np.real(pots[:, 0] * np.exp(1j * freqs * time))))
        direct = _direct_potential(degree_l, colatitude, longitude, time, eccentricity, obliquity)
        assert abs(summed - direct) < 1.0e-12 * scale
