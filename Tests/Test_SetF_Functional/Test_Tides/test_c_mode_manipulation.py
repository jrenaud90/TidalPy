import numpy as np
import numba

import TidalPy
TidalPy.use_disk = False

use_numba = TidalPy.configurations['use_numba']


def test_eccentricity_multi_l_calc():

    from TidalPy.tides.modeCalcHelper import eccentricity_functions_lookup

    # Test a few truncation levels
    eccen_trunc_funcset_2 = eccentricity_functions_lookup[2]
    eccen_trunc_funcset_10 = eccentricity_functions_lookup[10]
    eccen_trunc_funcset_20 = eccentricity_functions_lookup[20]

    for eccen_funcset in [eccen_trunc_funcset_2, eccen_trunc_funcset_10, eccen_trunc_funcset_20]:
        for tidal_order_l, eccen_func in eccen_funcset.items():

            # Perform float calculation
            e_result_float = eccen_func(0.3)
            for order_l in range(2, tidal_order_l+1):
                for p, p_result in e_result_float[order_l].items():
                    assert type(p) is int
                    for q, q_result in p_result.items():
                        assert type(q) is int
                        assert type(q_result) is float

            # Perform array calculation
            e_result_array = eccen_func(np.linspace(0.1, 0.9, 4))
            for order_l in range(2, tidal_order_l+1):
                for p, p_result in e_result_array[order_l].items():
                    for q, q_result in p_result.items():
                        assert type(q_result) is np.ndarray


def test_inclination_multi_l_calc():

    from TidalPy.tides.modeCalcHelper import inclination_functions_lookup

    # Test a few truncation levels
    inclination_funcset_on = inclination_functions_lookup[True]
    inclination_funcset_off = inclination_functions_lookup[False]

    for inclin_funcset in [inclination_funcset_on, inclination_funcset_off]:
        for tidal_order_l, inclin_func in inclin_funcset.items():

            # Perform float calculation
            i_result_float = inclin_func(0.3)
            for order_l in range(2, tidal_order_l+1):
                for (p, m), result in i_result_float[order_l].items():
                    assert type(p) is int
                    assert type(m) is int
                    assert type(result) is float

            # Perform array calculation
            i_result_array = inclin_func(np.linspace(0.1, 0.9, 4))
            for order_l in range(2, tidal_order_l+1):
                for (p, m), result in i_result_array[order_l].items():
                    assert type(p) is int
                    assert type(m) is int
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

    from TidalPy.tides.mode_manipulation import calculate_terms, collapse_modes
    from TidalPy.tides.modeCalcHelper import inclination_functions_lookup, eccentricity_functions_lookup

    for order_l in [2, 3, 5]:
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
                calculate_terms(spin_frequency, orbital_frequency, semi_major_axis, radius,
                                eccentricity_results_float, obliquity_results_float)
            unique_freq_array, tidal_results_array = \
                calculate_terms(spin_frequency_array, orbital_frequency_array, semi_major_axis_array, radius,
                                eccentricity_results_array, obliquity_results_array)

            complex_comp_float = \
                tuple([static_complex_compliance for _ in unique_freq_float])
            complex_comp_array = \
                tuple([static_complex_compliance * np.ones_like(spin_frequency_array) for _ in unique_freq_array])

            if use_numba:
                assert isinstance(unique_freq_float, numba.typed.typeddict.Dict)
                assert isinstance(tidal_results_float, numba.typed.typeddict.Dict)
            else:
                assert type(unique_freq_float) is dict
                assert type(tidal_results_float) is dict

            # Collapse modes
            result_float = \
                collapse_modes(gravity, radius, density, shear_modulus, tidal_scale, tidal_host_mass,
                               tidal_susceptibility, complex_comp_float, tidal_results_float, order_l,
                               cpl_ctl_method=False)

            tidal_heating, dUdM, dUdw, dUdO, love_number, negative_imk = result_float

            assert type(tidal_heating) is float
            assert type(dUdM) is float
            assert type(dUdw) is float
            assert type(dUdO) is float
            assert type(love_number) is complex
            assert type(negative_imk) is float

            result_array = \
                collapse_modes(gravity, radius, density, shear_modulus, tidal_scale, tidal_host_mass,
                               tidal_susceptibility_array, complex_comp_array, tidal_results_array, order_l,
                               cpl_ctl_method=False)

            tidal_heating, dUdM, dUdw, dUdO, love_number, negative_imk = result_array

            assert type(tidal_heating) is np.ndarray
            assert type(dUdM) is np.ndarray
            assert type(dUdw) is np.ndarray
            assert type(dUdO) is np.ndarray
            assert type(love_number) is np.ndarray
            assert type(negative_imk) is np.ndarray
