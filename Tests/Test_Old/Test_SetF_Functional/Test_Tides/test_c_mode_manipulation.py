import numba
import numpy as np

import TidalPy
from TidalPy.utilities.performance.numba import use_numba


def test_eccentricity_multi_l_calc():
    from TidalPy.tides.modes.mode_calc_helper import eccentricity_functions_lookup

    # Test a few truncation levels; these tests can take a long time to run so only doing a spot check on e^2 and e^10
    eccen_trunc_funcset_2 = eccentricity_functions_lookup[2]
    eccen_trunc_funcset_10 = eccentricity_functions_lookup[10]

    for eccen_funcset in [eccen_trunc_funcset_2, eccen_trunc_funcset_10]:
        for tidal_order_l, eccen_func in eccen_funcset.items():

            # Perform float calculation
            e_result_float = eccen_func(0.3)
            for order_l in range(2, tidal_order_l + 1):
                for p, p_result in e_result_float[order_l].items():
                    assert type(p) in [int, np.int32]
                    for q, q_result in p_result.items():
                        assert type(q) in [int, np.int32]
                        assert type(q_result) in [float, np.float64]

            # Perform array calculation
            e_result_array = eccen_func(np.linspace(0.1, 0.9, 4))
            for order_l in range(2, tidal_order_l + 1):
                for p, p_result in e_result_array[order_l].items():
                    for q, q_result in p_result.items():
                        assert type(q_result) is np.ndarray


def test_inclination_multi_l_calc():
    from TidalPy.tides.modes.mode_calc_helper import inclination_functions_lookup

    # Test a few truncation levels
    inclination_funcset_on = inclination_functions_lookup[True]
    inclination_funcset_off = inclination_functions_lookup[False]

    for inclin_funcset in [inclination_funcset_on, inclination_funcset_off]:
        for tidal_order_l, inclin_func in inclin_funcset.items():

            # Perform float calculation
            i_result_float = inclin_func(0.3)
            for order_l in range(2, tidal_order_l + 1):
                for (p, m), result in i_result_float[order_l].items():
                    assert type(p) in [int, np.int32]
                    assert type(m) in [int, np.int32]
                    assert type(result) in [float, np.float64]

            # Perform array calculation
            i_result_array = inclin_func(np.linspace(0.1, 0.9, 4))
            for order_l in range(2, tidal_order_l + 1):
                for (p, m), result in i_result_array[order_l].items():
                    assert type(p) in [int, np.int32]
                    assert type(m) in [int, np.int32]
                    assert type(result) is np.ndarray


# Float testing
spin_frequency, orbital_frequency = 0.001, 0.002
semi_major_axis = 100.
eccentricity, obliquity = 0.1, 0.2
radius = 10.
gravity = 9.
density = 1000.
shear_modulus = 1.e9
static_complex_compliance = shear_modulus**(-1) + 1.0j * (shear_modulus * 0.1)**(-1)
tidal_scale = 0.9
tidal_host_mass = 1.e12
tidal_susceptibility = 1.e22

# Array testing
spin_frequency_array = np.linspace(0.001, 0.005, 5)
orbital_frequency_array = 0.5 * spin_frequency_array
semi_major_axis_array = 100. * spin_frequency_array
eccentricity_array = np.linspace(0.1, 0.5, 5)
obliquity_array = np.linspace(0.1, 0.5, 5)
tidal_susceptibility_array = tidal_susceptibility * np.ones_like(orbital_frequency_array)


def test_calculate_and_collapse_modes():
    from TidalPy.utilities.performance import njit
    from TidalPy.tides.modes.mode_manipulation import calculate_terms, collapse_modes
    from TidalPy.tides.modes.mode_calc_helper import inclination_functions_lookup, eccentricity_functions_lookup

    # These tests can take a long time to run, so only doing a spot check on l=2,3 and e^2, e^10
    for order_l in [2, 3]:
        for truncation_level in [2, 10]:

            # Calculate eccentricity results
            eccentricity_func = eccentricity_functions_lookup[truncation_level][order_l]
            eccentricity_results_float = eccentricity_func(eccentricity)
            eccentricity_results_array = eccentricity_func(eccentricity_array)

            # Calculate obliquity results
            inclination_func = inclination_functions_lookup[True][order_l]
            obliquity_results_float = inclination_func(obliquity)
            obliquity_results_array = inclination_func(obliquity_array)

            if use_numba:
                assert isinstance(eccentricity_results_float, numba.typed.typeddict.Dict)
                assert isinstance(eccentricity_results_array, numba.typed.typeddict.Dict)
                assert isinstance(obliquity_results_float, numba.typed.typeddict.Dict)
                assert isinstance(obliquity_results_array, numba.typed.typeddict.Dict)
            else:
                assert type(eccentricity_results_float) is dict
                assert type(eccentricity_results_array) is dict
                assert type(obliquity_results_float) is dict
                assert type(obliquity_results_array) is dict

            # Calculate tidal frequencies and modes
            unique_freq_float, tidal_results_float = \
                calculate_terms(
                    spin_frequency, orbital_frequency, semi_major_axis, radius,
                    eccentricity_results_float, obliquity_results_float
                    )
            unique_freq_array, tidal_results_array = \
                calculate_terms(
                    spin_frequency_array, orbital_frequency_array, semi_major_axis_array, radius,
                    eccentricity_results_array, obliquity_results_array
                    )

            @njit()
            def build_numba_dict(freq_dict, static_comp):

                fake_index = list(freq_dict.keys())[0]
                fake_freq = freq_dict[fake_index]
                complex_comp_dict = {(-100, -100): fake_freq * (1. + 1.j)}

                for freq_sig in freq_dict:
                    complex_comp_dict[freq_sig] = static_comp

                del complex_comp_dict[(-100, -100)]
                return complex_comp_dict

            complex_comp_float = build_numba_dict(unique_freq_float, static_complex_compliance)
            complex_comp_array = build_numba_dict(
                unique_freq_array,
                static_complex_compliance * np.ones_like(spin_frequency_array)
                )

            if use_numba:
                assert isinstance(unique_freq_float, numba.typed.typeddict.Dict)
                assert isinstance(tidal_results_float, numba.typed.typeddict.Dict)
            else:
                assert type(unique_freq_float) is dict
                assert type(tidal_results_float) is dict

            # Collapse modes
            result_float = \
                collapse_modes(
                    gravity, radius, density, shear_modulus, tidal_scale, tidal_host_mass,
                    tidal_susceptibility, complex_comp_float, tidal_results_float, order_l,
                    cpl_ctl_method=False
                    )

            tidal_heating, dUdM, dUdw, dUdO, love_number_by_orderl, negative_imk_by_orderl, effective_q_by_orderl = \
                result_float

            assert type(tidal_heating) in [float, np.float64]
            assert type(dUdM) in [float, np.float64]
            assert type(dUdw) in [float, np.float64]
            assert type(dUdO) in [float, np.float64]
            assert type(love_number_by_orderl) in [dict, numba.typed.typeddict.Dict]
            assert len(love_number_by_orderl) == order_l - 2 + 1
            assert type(love_number_by_orderl[2]) in [complex, np.complex64]
            assert type(negative_imk_by_orderl) in [dict, numba.typed.typeddict.Dict]
            assert len(negative_imk_by_orderl) == order_l - 2 + 1
            assert type(negative_imk_by_orderl[2]) in [float, np.float64]
            assert type(effective_q_by_orderl) in [dict, numba.typed.typeddict.Dict]
            assert len(effective_q_by_orderl) == order_l - 2 + 1
            assert type(effective_q_by_orderl[2]) in [float, np.float64]

            result_array = \
                collapse_modes(
                    gravity, radius, density, shear_modulus, tidal_scale, tidal_host_mass,
                    tidal_susceptibility_array, complex_comp_array, tidal_results_array, order_l,
                    cpl_ctl_method=False
                    )

            tidal_heating, dUdM, dUdw, dUdO, love_number_by_orderl, negative_imk_by_orderl, effective_q_by_orderl = \
                result_array

            assert type(tidal_heating) is np.ndarray
            assert type(dUdM) is np.ndarray
            assert type(dUdw) is np.ndarray
            assert type(dUdO) is np.ndarray
            assert type(love_number_by_orderl) in [dict, numba.typed.typeddict.Dict]
            assert len(love_number_by_orderl) == order_l - 2 + 1
            assert type(love_number_by_orderl[2]) is np.ndarray
            assert type(negative_imk_by_orderl) in [dict, numba.typed.typeddict.Dict]
            assert len(negative_imk_by_orderl) == order_l - 2 + 1
            assert type(negative_imk_by_orderl[2]) is np.ndarray
            assert type(effective_q_by_orderl) in [dict, numba.typed.typeddict.Dict]
            assert len(effective_q_by_orderl) == order_l - 2 + 1
            assert type(effective_q_by_orderl[2]) is np.ndarray
