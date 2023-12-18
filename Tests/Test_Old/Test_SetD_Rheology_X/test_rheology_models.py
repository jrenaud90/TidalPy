""" Test TidalPy's rheology models used to find complex shear or bulk. """
import math

import numpy as np
import pytest

import TidalPy



from TidalPy.rheology.models import find_rheology

viscosity = 1.0e18
shear_mod = 50.0e9
frequency = 1.0e-6

expected_results = {
    'off'           : ('ElasticRheology', shear_mod + 0.0j),
    'elastic'       : ('ElasticRheology', shear_mod + 0.0j),
    'newton'        : ('NewtonRheology', 0.0 + 1.0j * (viscosity * frequency)),
    'viscous'       : ('NewtonRheology', 0.0 + 1.0j * (viscosity * frequency)),
    'maxwell'       : ('MaxwellRheology', 49875311720.69826 + 2493765586.034913j),
    'voigt'         : ('VoigtRheology', 250000000000 + 19999999999.999996j),
    'voigtkelvin'   : ('VoigtRheology', 250000000000 + 19999999999.999996j),
    'burgers'       : ('BurgersRheology', 41585201404.18996 + 2286083021.4902835j),
    'andrade'       : ('AndradeRheology', 36746190752.118454 + 5984215754.873647j),
    'sundberg'      : ('SundbergCooperRheology', 32061580147.854797 + 4874982764.688812j),
    'sundbergcooper': ('SundbergCooperRheology', 32061580147.854797 + 4874982764.688812j)
    }


@pytest.mark.parametrize('model_name', expected_results.keys())
def test_build_rheology_models(model_name):
    # Find the correct rheology class
    RheologyClass = find_rheology(model_name)

    # Instantiate it.
    rheology_instance = RheologyClass()
    assert rheology_instance.class_name == expected_results[model_name][0]


@pytest.mark.parametrize('model_name', expected_results.keys())
def test_rheology_models_accuracy(model_name):
    # Find the correct rheology class
    RheologyClass = find_rheology(model_name)

    # Instantiate it.
    rheology_instance = RheologyClass()

    # Get result
    complex_shear = rheology_instance(frequency, shear_mod, viscosity)

    assert math.isfinite(complex_shear.real)
    assert math.isfinite(complex_shear.imag)
    assert type(complex_shear) == complex
    assert math.isclose(complex_shear.real, expected_results[model_name][1].real)
    assert math.isclose(complex_shear.imag, expected_results[model_name][1].imag)


@pytest.mark.parametrize('model_name', expected_results.keys())
def test_rheology_models_optional_args(model_name):
    if model_name in ('off', 'elastic', 'viscous', 'newton', 'maxwell'):
        # These models do not have other inputs
        assert True
    else:
        # Find the correct rheology class
        RheologyClass = find_rheology(model_name)

        if model_name == 'andrade':
            # Alpha, Zeta
            rheology_instance = RheologyClass((0.1, 10.))
            complex_shear = rheology_instance(frequency, shear_mod, viscosity)
            assert math.isfinite(complex_shear.real)
            assert math.isfinite(complex_shear.imag)
            assert math.isclose(complex_shear.real, 31941504195.267456)
            assert math.isclose(complex_shear.imag, 2830070091.949614)
        elif model_name in ('voigt', 'voigtkelvin'):
            # Shear Scale, Viscosity Scale
            rheology_instance = RheologyClass((0.1, 100.))
            complex_shear = rheology_instance(frequency, shear_mod, viscosity)
            assert math.isfinite(complex_shear.real)
            assert math.isfinite(complex_shear.imag)
            assert math.isclose(complex_shear.real, 5000000000.000003)
            assert math.isclose(complex_shear.imag, 100000000000000.03)
        elif model_name == 'burgers':
            # Shear Scale, Viscosity Scale
            rheology_instance = RheologyClass((0.1, 100.))
            complex_shear = rheology_instance(frequency, shear_mod, viscosity)
            assert math.isfinite(complex_shear.real)
            assert math.isfinite(complex_shear.imag)
            assert math.isclose(complex_shear.real, 49872810621.07933)
            assert math.isclose(complex_shear.imag, 2518576873.3377433)
        elif model_name in ('sundberg', 'sundbergcooper'):
            # Shear Scale, Viscosity Scale, alpha, zeta
            rheology_instance = RheologyClass((0.1, 100., 0.1, 10.))
            complex_shear = rheology_instance(frequency, shear_mod, viscosity)
            assert math.isfinite(complex_shear.real)
            assert math.isfinite(complex_shear.imag)
            assert math.isclose(complex_shear.real, 31939692573.57115)
            assert math.isclose(complex_shear.imag, 2840191640.197709)


@pytest.mark.parametrize('model_name', expected_results.keys())
def test_rheology_models_changing_args(model_name):
    if model_name in ('off', 'elastic', 'viscous', 'newton', 'maxwell'):
        # These models do not have other inputs
        assert True
    else:
        # Find the correct rheology class
        RheologyClass = find_rheology(model_name)

        # Instantiate with no args
        rheology_instance = RheologyClass()

        # Run once
        complex_shear = rheology_instance(frequency, shear_mod, viscosity)

        # Change args
        if model_name in ('sundberg', 'sundbergcooper'):
            # S-C has 4 args.
            rheology_instance.change_args((0.9, 0.9, 0.9, 0.9))
        else:
            # Others have 2.
            rheology_instance.change_args((0.9, 0.9))

        # Run and make sure results are different.
        complex_shear_2 = rheology_instance(frequency, shear_mod, viscosity)
        assert not math.isclose(complex_shear.real, complex_shear_2.real)
        assert not math.isclose(complex_shear.imag, complex_shear_2.imag)


@pytest.mark.parametrize('model_name', expected_results.keys())
@pytest.mark.parametrize('vectorized_freq', (True, False))
def test_rheology_models_arrays(model_name, vectorized_freq):
    # Find the correct rheology class
    RheologyClass = find_rheology(model_name)

    # Instantiate with no args
    rheology_instance = RheologyClass()

    N = 100
    output_array = np.empty(N, dtype=np.complex128, order='C')
    if vectorized_freq:
        # With a vectorized frequency
        freq_array = frequency * np.ones(N, dtype=np.float64, order='C')
        rheology_instance.vectorize_frequency(freq_array, shear_mod, viscosity, output_array)
    else:
        # With vectorized thermals (shear and viscosity)
        shear_array = shear_mod * np.ones(N, dtype=np.float64, order='C')
        visc_array = viscosity * np.ones(N, dtype=np.float64, order='C')
        rheology_instance.vectorize_modulus_viscosity(frequency, shear_array, visc_array, output_array)

    # Make an expected array and test.
    expected_array = expected_results[model_name][1] * np.ones(N, dtype=np.complex128, order='C')
    assert np.allclose(output_array, expected_array)
