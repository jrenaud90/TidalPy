from math import isclose, isfinite

import pytest
import numpy as np

from TidalPy.constants import G
from TidalPy.Tides_x.potential import global_potential
from TidalPy.Tides_x.potential.potential_common import ModeMap


# Io-Jupiter-like test parameters
host_mass = 1.8982e27        # Jupiter mass [kg]
planet_radius = 1.8216e6     # Io radius [m]
orbital_period = 1.769 * 86400.0  # Io orbital period [s]
orbital_frequency = 2.0 * np.pi / orbital_period
semi_major_axis = 4.217e8    # Io semi-major axis [m]


def _global_potential(**kwargs):
    """Call `global_potential` by keyword so an argument reorder cannot silently shuffle the inputs."""
    arguments = dict(
        planet_radius=planet_radius,
        orbital_frequency=orbital_frequency,
        spin_frequency=orbital_frequency,  # synchronous
        eccentricity=0.0,
        obliquity=0.0,
        semi_major_axis=semi_major_axis,
        host_mass=host_mass,
        G_to_use=G,
        min_degree_l=2,
        max_degree_l=2,
    )
    arguments.update(kwargs)
    return global_potential(**arguments)


@pytest.mark.parametrize('degree_l', (2, 3, 4))
@pytest.mark.parametrize('obliquity_truncation', ('gen', 1, 2, 'off'))
@pytest.mark.parametrize('eccentricity_truncation', (2, 4, 6, 8, 10, 20))
def test_global_potential_basic(degree_l, obliquity_truncation, eccentricity_truncation):
    """Tests that global_potential runs without error and returns correct types for various parameters."""

    mode_map, unique_freq_index_map, unique_freq_list, potential_dict = _global_potential(
        obliquity=0.01,
        eccentricity=0.0041,
        min_degree_l=degree_l,
        max_degree_l=degree_l,
        obliquity_truncation=obliquity_truncation,
        eccentricity_truncation=eccentricity_truncation
    )

    # Check return types
    assert isinstance(mode_map, ModeMap)
    assert isinstance(unique_freq_list, list)
    assert isinstance(potential_dict, dict)

    # Every truncation keeps the e^2 heating terms, so a synchronous eccentric orbit always has active modes.
    assert len(mode_map) > 0
    assert len(unique_freq_list) > 0
    assert len(potential_dict) > 0

    # Check that modes and potentials have the same keys
    mode_keys = set()
    for (l, m, p, q), mode_data in mode_map:
        mode_keys.add((l, m, p, q))
        assert l == degree_l
        assert 0 <= m <= l
        assert 0 <= p <= l
        # mode_data is (mode, mode_strength, n_coeff, o_coeff)
        assert len(mode_data) == 4
        mode_val, mode_strength, n_coeff, o_coeff = mode_data
        assert isinstance(mode_val, float)
        assert isinstance(mode_strength, float)
        assert isinstance(n_coeff, int)
        assert isinstance(o_coeff, int)
        # n_coeff should be l - 2p + q
        assert n_coeff == l - 2 * p + q
        # o_coeff should be -m
        assert o_coeff == -m
        assert isclose(mode_val, n_coeff * orbital_frequency + o_coeff * orbital_frequency, rel_tol=1e-12)
        assert isfinite(mode_strength)

    potential_keys = set(potential_dict.keys())
    assert mode_keys == potential_keys

    # Check potential dict values
    for key, (dU_dM, dU_dw, dU_dO, E_dot) in potential_dict.items():
        assert isinstance(dU_dM, float)
        assert isinstance(dU_dw, float)
        assert isinstance(dU_dO, float)
        assert isinstance(E_dot, float)
        assert isfinite(dU_dM) and isfinite(dU_dw) and isfinite(dU_dO) and isfinite(E_dot)
        # Heating should always be non-negative (it's |mode| * host_mass * common_coeff)
        assert E_dot >= 0.0

    # Check unique frequency list entries
    for freq, num_instances in unique_freq_list:
        assert isinstance(freq, float)
        assert freq > 0.0
        assert isinstance(num_instances, int)
        assert num_instances >= 1


@pytest.mark.parametrize('obliquity_truncation', ('gen', 2, 'off'))
def test_global_potential_zero_obliquity(obliquity_truncation):
    """Tests global_potential with zero obliquity (common simplification)."""

    mode_map, unique_freq_index_map, unique_freq_list, potential_dict = _global_potential(
        obliquity=0.0,
        eccentricity=0.1,
        obliquity_truncation=obliquity_truncation,
        eccentricity_truncation=10
    )

    assert len(mode_map) > 0
    assert len(potential_dict) > 0

    # With zero obliquity, only m=0 and m=2 modes should be present for l=2
    # (m=1 modes vanish when obliquity is zero in the 'off' truncation).
    for (l, m, p, q), mode_data in mode_map:
        assert l == 2
        if obliquity_truncation == 'off':
            assert m in (0, 2)


@pytest.mark.parametrize('eccentricity_truncation', (2, 4, 6, 8, 10, 20))
def test_global_potential_zero_eccentricity(eccentricity_truncation):
    """Tests global_potential with zero eccentricity."""

    mode_map, unique_freq_index_map, unique_freq_list, potential_dict = _global_potential(
        obliquity=0.3,
        eccentricity=0.0,
        obliquity_truncation='gen',
        eccentricity_truncation=eccentricity_truncation
    )

    # Obliquity alone forces the m = 1 modes at synchronous rotation.
    assert len(mode_map) > 0

    # With zero eccentricity, only q=0 modes should survive (G_lpq is only nonzero for q=0 at e=0).
    for (l_, m_, p_, q_), mode_data in mode_map:
        assert q_ == 0


def test_global_potential_synchronous_zero_obliquity():
    """Tests that synchronous rotation with zero obliquity yields no nonzero-frequency modes
    for the dominant (l=2, off obliquity, low eccentricity) case, except eccentricity-driven ones."""

    mode_map, unique_freq_index_map, unique_freq_list, potential_dict = _global_potential(
        obliquity=0.0,
        eccentricity=0.0,
        obliquity_truncation='off',
        eccentricity_truncation=2
    )

    # With zero obliquity and zero eccentricity at synchronous rotation:
    # The only nonzero modes come from eccentricity. At e=0, very few modes survive.
    # All surviving modes should have zero frequency (synchronous + zero ecc).
    # Since zero-frequency modes are skipped, the mode map may be empty.
    # This is a valid physical result.
    assert len(mode_map) == 0
    assert len(potential_dict) == 0


def test_global_potential_multi_degree():
    """Tests global_potential spanning multiple degree l values."""

    mode_map, unique_freq_index_map, unique_freq_list, potential_dict = _global_potential(
        spin_frequency=orbital_frequency * 1.5,  # non-synchronous spin
        obliquity=0.1,
        eccentricity=0.1,
        min_degree_l=2,
        max_degree_l=3,
        obliquity_truncation='gen',
        eccentricity_truncation=4
    )

    assert len(mode_map) > 0

    # Check that we get modes from both l=2 and l=3
    degree_ls_found = set()
    for (l, m, p, q), mode_data in mode_map:
        degree_ls_found.add(l)
        assert l in (2, 3)
        assert 0 <= m <= l
        assert 0 <= p <= l

    assert 2 in degree_ls_found
    assert 3 in degree_ls_found


def test_global_potential_mode_strength_normalization():
    """Tests that the relative mode strengths are normalized (max abs value == 1)."""

    mode_map, unique_freq_index_map, unique_freq_list, potential_dict = _global_potential(
        spin_frequency=orbital_frequency * 1.5,
        obliquity=0.2,
        eccentricity=0.1,
        obliquity_truncation='gen',
        eccentricity_truncation=4
    )

    assert len(mode_map) > 0

    max_abs_strength = 0.0
    for (l, m, p, q), mode_data in mode_map:
        mode_val, mode_strength, n_coeff, o_coeff = mode_data
        max_abs_strength = max(max_abs_strength, abs(mode_strength))

    # The maximum absolute mode strength should be 1.0 after normalization.
    assert isclose(max_abs_strength, 1.0, rel_tol=1e-12)


def test_global_potential_nonsynchronous():
    """Tests global_potential with non-synchronous spin to verify mode frequency calculation."""

    spin_frequency = 2.0 * orbital_frequency  # 2:1 spin-orbit

    mode_map, unique_freq_index_map, unique_freq_list, potential_dict = _global_potential(
        spin_frequency=spin_frequency,
        obliquity=0.0,
        eccentricity=0.05,
        obliquity_truncation='off',
        eccentricity_truncation=4
    )

    assert len(mode_map) > 0

    # Check mode frequency calculation: mode = n_coeff * n + o_coeff * Omega
    for (l, m, p, q), mode_data in mode_map:
        mode_val, mode_strength, n_coeff, o_coeff = mode_data
        expected_mode = n_coeff * orbital_frequency + o_coeff * spin_frequency
        assert isclose(mode_val, expected_mode, rel_tol=1e-12)


@pytest.mark.parametrize('eccentricity_truncation', (2, 6))
def test_global_potential_synchronous_heating_matches_analytic(eccentricity_truncation):
    """At synchronous rotation, zero obliquity, and low eccentricity, the summed mode heating per unit k2/Q
    equals the classic (21/2) G M_host^2 R^5 n e^2 / a^6, from the lowest truncation up (exactly at level 2)."""

    eccentricity = 0.01

    mode_map, unique_freq_index_map, unique_freq_list, potential_dict = _global_potential(
        eccentricity=eccentricity,
        obliquity_truncation='off',
        eccentricity_truncation=eccentricity_truncation
    )

    heating_per_k_over_q = sum(E_dot for (dU_dM, dU_dw, dU_dO, E_dot) in potential_dict.values())
    expected = (21.0 / 2.0) * G * host_mass**2 * planet_radius**5 * orbital_frequency * eccentricity**2 \
        / semi_major_axis**6
    assert isclose(heating_per_k_over_q, expected, rel_tol=5.0e-3)


def test_global_potential_unsupported_obliquity_truncation():
    """Tests that an unsupported obliquity truncation raises NotImplementedError."""

    with pytest.raises(NotImplementedError):
        _global_potential(
            obliquity=0.1,
            eccentricity=0.1,
            obliquity_truncation=5,
            eccentricity_truncation=4
        )

def test_global_potential_unsupported_eccentricity_truncation():
    """Tests that an unsupported eccentricity truncation raises NotImplementedError."""

    with pytest.raises(NotImplementedError):
        _global_potential(
            obliquity=0.1,
            eccentricity=0.1,
            obliquity_truncation=2,
            eccentricity_truncation=7
        )

def test_global_potential_spot_check_l2_off_sync():
    """Spot check: l=2, off obliquity, synchronous rotation.
    With e > 0, the dominant mode should be (2, 0, 1, -1) and (2, 2, 0, 1) etc."""

    mode_map, unique_freq_index_map, unique_freq_list, potential_dict = _global_potential(
        obliquity=0.0,
        eccentricity=0.1,
        obliquity_truncation='off',
        eccentricity_truncation=2
    )

    assert len(mode_map) > 0

    # With zero obliquity + 'off' truncation, only m=0 (p=1) and m=2 (p=0) contribute.
    # For synchronous rotation, mode = (l-2p+q)*n + (-m)*n = (2-2p+q-m)*n
    # So for (2,0,1,q): mode = q*n ; for (2,2,0,q): mode = (2+q-2)*n = q*n
    # The q=0 modes have zero frequency and are skipped.
    for (l, m, p, q), mode_data in mode_map:
        assert l == 2
        assert m in (0, 2)
        assert q != 0  # zero-freq modes are skipped
        mode_val = mode_data[0]
        # Verify mode is nonzero
        assert abs(mode_val) > 0.0
