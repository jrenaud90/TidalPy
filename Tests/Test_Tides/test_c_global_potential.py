"""The global (1D) tidal potential: mode maps, frequencies, strengths, and heating for an Io-like body."""
from math import isclose, isfinite

import pytest
import numpy as np

from TidalPy.constants import G
from TidalPy.Tides.potential import global_potential
from TidalPy.Tides.potential.potential_common import ModeMap

HOST_MASS = 1.8982e27
PLANET_RADIUS = 1.8216e6
ORBITAL_FREQUENCY = 2.0 * np.pi / (1.769 * 86400.0)
SEMI_MAJOR_AXIS = 4.217e8


def _global_potential(**kwargs):
    """Synchronous degree-2 `global_potential`, called by keyword so an argument reorder cannot shuffle inputs."""
    arguments = dict(
        planet_radius=PLANET_RADIUS,
        orbital_frequency=ORBITAL_FREQUENCY,
        spin_frequency=ORBITAL_FREQUENCY,
        eccentricity=0.0,
        obliquity=0.0,
        semi_major_axis=SEMI_MAJOR_AXIS,
        host_mass=HOST_MASS,
        G_to_use=G,
        min_degree_l=2,
        max_degree_l=2,
    )
    arguments.update(kwargs)
    return global_potential(**arguments)


@pytest.mark.parametrize('degree_l', (2, 3, 4))
@pytest.mark.parametrize('obliquity_truncation', ('gen', 2, 4, 'off'))
@pytest.mark.parametrize('eccentricity_truncation', (2, 4, 6, 8, 10, 20))
def test_global_potential_basic(degree_l, obliquity_truncation, eccentricity_truncation):
    """Return types, mode coefficients and frequencies, and finite nonnegative heating for every truncation."""
    mode_map, unique_freq_index_map, unique_freq_list, potential_dict = _global_potential(
        obliquity=0.01,
        eccentricity=0.0041,
        min_degree_l=degree_l,
        max_degree_l=degree_l,
        obliquity_truncation=obliquity_truncation,
        eccentricity_truncation=eccentricity_truncation)

    assert isinstance(mode_map, ModeMap)
    assert isinstance(unique_freq_list, list)
    assert isinstance(potential_dict, dict)

    # Every truncation keeps the e^2 heating terms, so there are always active modes.
    assert len(mode_map) > 0
    assert len(unique_freq_list) > 0
    assert len(potential_dict) > 0

    mode_keys = set()
    for (l, m, p, q), mode_data in mode_map:
        mode_keys.add((l, m, p, q))
        assert l == degree_l
        assert 0 <= m <= l
        assert 0 <= p <= l
        assert len(mode_data) == 4
        mode_val, mode_strength, n_coeff, o_coeff = mode_data
        assert isinstance(mode_val, float)
        assert isinstance(mode_strength, float)
        assert isinstance(n_coeff, int)
        assert isinstance(o_coeff, int)
        assert n_coeff == l - 2 * p + q
        assert o_coeff == -m
        assert isclose(mode_val, n_coeff * ORBITAL_FREQUENCY + o_coeff * ORBITAL_FREQUENCY, rel_tol=1e-12)
        assert isfinite(mode_strength)

    assert mode_keys == set(potential_dict.keys())

    for key, (dU_dM, dU_dw, dU_dO, E_dot) in potential_dict.items():
        assert isinstance(dU_dM, float)
        assert isinstance(dU_dw, float)
        assert isinstance(dU_dO, float)
        assert isinstance(E_dot, float)
        assert isfinite(dU_dM) and isfinite(dU_dw) and isfinite(dU_dO) and isfinite(E_dot)
        assert E_dot >= 0.0

    for freq, num_instances in unique_freq_list:
        assert isinstance(freq, float)
        assert freq > 0.0
        assert isinstance(num_instances, int)
        assert num_instances >= 1


@pytest.mark.parametrize('obliquity_truncation', ('gen', 2, 'off'))
def test_global_potential_zero_obliquity(obliquity_truncation):
    """Zero obliquity gives degree-2 modes, and only m = 0 and m = 2 under the 'off' truncation."""
    mode_map, unique_freq_index_map, unique_freq_list, potential_dict = _global_potential(
        obliquity=0.0,
        eccentricity=0.1,
        obliquity_truncation=obliquity_truncation,
        eccentricity_truncation=10)

    assert len(mode_map) > 0
    assert len(potential_dict) > 0
    for (l, m, p, q), mode_data in mode_map:
        assert l == 2
        if obliquity_truncation == 'off':
            assert m in (0, 2)


@pytest.mark.parametrize('eccentricity_truncation', (2, 4, 6, 8, 10, 20))
def test_global_potential_zero_eccentricity(eccentricity_truncation):
    """Zero eccentricity keeps only q = 0 modes (obliquity alone forces the synchronous m = 1 modes)."""
    mode_map, unique_freq_index_map, unique_freq_list, potential_dict = _global_potential(
        obliquity=0.3,
        eccentricity=0.0,
        obliquity_truncation='gen',
        eccentricity_truncation=eccentricity_truncation)

    assert len(mode_map) > 0
    for (l, m, p, q), mode_data in mode_map:
        assert q == 0


def test_global_potential_synchronous_zero_obliquity():
    """A circular, zero-obliquity, synchronous orbit has only zero-frequency modes, which are skipped."""
    mode_map, unique_freq_index_map, unique_freq_list, potential_dict = _global_potential(
        obliquity=0.0,
        eccentricity=0.0,
        obliquity_truncation='off',
        eccentricity_truncation=2)

    assert len(mode_map) == 0
    assert len(potential_dict) == 0


def test_global_potential_multi_degree():
    """A degree range 2 to 3 returns modes from both degrees."""
    mode_map, unique_freq_index_map, unique_freq_list, potential_dict = _global_potential(
        spin_frequency=ORBITAL_FREQUENCY * 1.5,
        obliquity=0.1,
        eccentricity=0.1,
        min_degree_l=2,
        max_degree_l=3,
        obliquity_truncation='gen',
        eccentricity_truncation=4)

    assert len(mode_map) > 0
    degree_ls_found = set()
    for (l, m, p, q), mode_data in mode_map:
        degree_ls_found.add(l)
        assert l in (2, 3)
        assert 0 <= m <= l
        assert 0 <= p <= l
    assert 2 in degree_ls_found
    assert 3 in degree_ls_found


def test_global_potential_mode_strength_normalization():
    """The largest absolute mode strength is normalized to 1."""
    mode_map, unique_freq_index_map, unique_freq_list, potential_dict = _global_potential(
        spin_frequency=ORBITAL_FREQUENCY * 1.5,
        obliquity=0.2,
        eccentricity=0.1,
        obliquity_truncation='gen',
        eccentricity_truncation=4)

    assert len(mode_map) > 0
    max_abs_strength = 0.0
    for (l, m, p, q), mode_data in mode_map:
        max_abs_strength = max(max_abs_strength, abs(mode_data[1]))
    assert isclose(max_abs_strength, 1.0, rel_tol=1e-12)


def test_global_potential_nonsynchronous():
    """At 2:1 spin-orbit each mode frequency is n_coeff n + o_coeff spin."""
    spin_frequency = 2.0 * ORBITAL_FREQUENCY
    mode_map, unique_freq_index_map, unique_freq_list, potential_dict = _global_potential(
        spin_frequency=spin_frequency,
        obliquity=0.0,
        eccentricity=0.05,
        obliquity_truncation='off',
        eccentricity_truncation=4)

    assert len(mode_map) > 0
    for (l, m, p, q), mode_data in mode_map:
        mode_val, mode_strength, n_coeff, o_coeff = mode_data
        assert isclose(mode_val, n_coeff * ORBITAL_FREQUENCY + o_coeff * spin_frequency, rel_tol=1e-12)


@pytest.mark.parametrize('eccentricity_truncation', (2, 6))
def test_global_potential_synchronous_heating_matches_analytic(eccentricity_truncation):
    """Summed heating per unit k2/Q matches (21/2) G M_host^2 R^5 n e^2 / a^6 at low eccentricity."""
    eccentricity = 0.01
    mode_map, unique_freq_index_map, unique_freq_list, potential_dict = _global_potential(
        eccentricity=eccentricity,
        obliquity_truncation='off',
        eccentricity_truncation=eccentricity_truncation)

    heating_per_k_over_q = sum(E_dot for (dU_dM, dU_dw, dU_dO, E_dot) in potential_dict.values())
    expected = (21.0 / 2.0) * G * HOST_MASS**2 * PLANET_RADIUS**5 * ORBITAL_FREQUENCY * eccentricity**2 \
        / SEMI_MAJOR_AXIS**6
    assert isclose(heating_per_k_over_q, expected, rel_tol=5.0e-3)


@pytest.mark.parametrize('obliquity_truncation, eccentricity_truncation', [(5, 4), (2, 7)],
                         ids=['obliquity', 'eccentricity'])
def test_global_potential_unsupported_truncation(obliquity_truncation, eccentricity_truncation):
    """An untabulated obliquity or eccentricity truncation raises."""
    with pytest.raises(NotImplementedError):
        _global_potential(
            obliquity=0.1,
            eccentricity=0.1,
            obliquity_truncation=obliquity_truncation,
            eccentricity_truncation=eccentricity_truncation)


def test_global_potential_spot_check_l2_off_sync():
    """Synchronous, 'off' obliquity, e > 0: only m = 0 and m = 2 modes with q != 0 and nonzero frequency."""
    mode_map, unique_freq_index_map, unique_freq_list, potential_dict = _global_potential(
        obliquity=0.0,
        eccentricity=0.1,
        obliquity_truncation='off',
        eccentricity_truncation=2)

    assert len(mode_map) > 0
    # Synchronous: both (2, 0, 1, q) and (2, 2, 0, q) have frequency q n, so q = 0 is static and skipped.
    for (l, m, p, q), mode_data in mode_map:
        assert l == 2
        assert m in (0, 2)
        assert q != 0
        assert abs(mode_data[0]) > 0.0
