import numpy as np
from math import isclose, isnan

from TidalPy.rheology import (
    Elastic,
    Newton,
    Maxwell,
    Voigt,
    Burgers,
    Andrade,
    SundbergCooper
    )
from TidalPy.rheology import find_rheology

rheology_classes = {
    'Elastic': Elastic,
    'Newton': Newton,
    'Maxwell': Maxwell,
    'Voigt': Voigt,
    'Burgers': Burgers,
    'Andrade': Andrade,
    'SundbergCooper': SundbergCooper
}

def test_rheology_realistic():
    """ Tests rheologies with default parameters and realistic inputs. """

    viscosity = 1.0e18
    shear = 50.0e9
    frequency = 1.0e-6

    expected_results = {
        'Elastic': (5.0000000000e+10, 0.0000000000e+00),
        'Newton': (0.0000000000e+00, 1.0000000000e+12),
        'Maxwell': (4.9875311721e+10, 2.4937655860e+09),
        'Voigt': (2.5000000000e+11, 2.0000000000e+10),
        'Burgers': (4.1585201404e+10, 2.2860830215e+09),
        'Andrade': (3.6746190752e+10, 5.9842157549e+09),
        'SundbergCooper': (3.2061580148e+10, 4.8749827647e+09)
    }

    for rheology_name, rheology_class in rheology_classes.items():
        rheology_instance = rheology_class()
        expected_real, expected_imag = expected_results[rheology_name]
        result = rheology_instance(frequency, shear, viscosity)
        assert type(result) == complex
        assert not isnan(np.real(result))
        assert not isnan(np.imag(result))

        assert isclose(np.real(result), expected_real, rel_tol=1.0e-6, abs_tol=1.0e-12)
        assert isclose(np.imag(result), expected_imag, rel_tol=1.0e-6, abs_tol=1.0e-12)

        # Results should not change if frequency is negative
        result = rheology_instance(-frequency, shear, viscosity)
        assert isclose(np.real(result), expected_real, rel_tol=1.0e-6, abs_tol=1.0e-12)
        assert isclose(np.imag(result), expected_imag, rel_tol=1.0e-6, abs_tol=1.0e-12)


def test_rheology_zero_frequency():
    """ Tests rheologies with default parameters and zero frequency """

    viscosity = 1.0e18
    shear = 50.0e9
    frequency = 0.

    expected_results = {
        'Elastic': (shear, 0.0000000000e+00),
        'Newton': (0.0000000000e+00, 0.0),
        'Maxwell': (0.0, 0.0),
        'Voigt': (5. * shear, 0.0),  # voigt shear offset * shear
        'Burgers': (0.0, 0.0),
        'Andrade': (0.0, 0.0),
        'SundbergCooper': (0.0, 0.0)
    }

    for rheology_name, rheology_class in rheology_classes.items():
        rheology_instance = rheology_class()
        expected_real, expected_imag = expected_results[rheology_name]
        result = rheology_instance(frequency, shear, viscosity)
        assert type(result) == complex
        assert not isnan(np.real(result))
        assert not isnan(np.imag(result))

        assert isclose(np.real(result), expected_real, rel_tol=1.0e-6, abs_tol=1.0e-12)
        assert isclose(np.imag(result), expected_imag, rel_tol=1.0e-6, abs_tol=1.0e-12)

        # Results should not change if frequency is negative
        result = rheology_instance(-frequency, shear, viscosity)
        assert isclose(np.real(result), expected_real, rel_tol=1.0e-6, abs_tol=1.0e-12)
        assert isclose(np.imag(result), expected_imag, rel_tol=1.0e-6, abs_tol=1.0e-12)


def test_rheology_inf_frequency():
    """ Tests rheologies with default parameters and infinite frequency. """

    viscosity = 1.0e18
    shear = 50.0e9
    frequency = np.inf

    expected_results = {
        'Elastic': (shear, 0.0000000000e+00),
        'Newton': (0.0, np.inf),
        'Maxwell': (shear, 0.0),
        'Voigt': (0.0, np.inf),  # voigt shear offset * shear
        'Burgers': (shear, 0.0),
        'Andrade': (shear, 0.0),
        'SundbergCooper': (shear, 0.0)
    }

    for rheology_name, rheology_class in rheology_classes.items():
        rheology_instance = rheology_class()
        expected_real, expected_imag = expected_results[rheology_name]
        result = rheology_instance(frequency, shear, viscosity)
        assert type(result) == complex

        assert not isnan(np.real(result))
        assert not isnan(np.imag(result))

        assert isclose(np.real(result), expected_real, rel_tol=1.0e-6, abs_tol=1.0e-12)
        assert isclose(np.imag(result), expected_imag, rel_tol=1.0e-6, abs_tol=1.0e-12)

        # Results should not change if frequency is negative
        result = rheology_instance(-frequency, shear, viscosity)
        assert isclose(np.real(result), expected_real, rel_tol=1.0e-6, abs_tol=1.0e-12)
        assert isclose(np.imag(result), expected_imag, rel_tol=1.0e-6, abs_tol=1.0e-12)


def test_rheology_vectorize_frequency():
    """ Tests rheologies vectorize frequency method. """

    viscosity = 1.0e18
    shear = 50.0e9
    frequency = np.logspace(-6, -3, 10)

    for rheology_name, rheology_class in rheology_classes.items():
        rheology_instance = rheology_class()
        result = np.empty(frequency.size, dtype=np.complex128)
        rheology_instance.vectorize_frequency(frequency, shear, viscosity, result)
        assert type(result) == np.ndarray
        assert result.dtype == np.complex128


def test_rheology_vectorize_modulus_viscosity():
    """ Tests rheologies vectorize modulus and viscosity method. """

    viscosity = 1.0e18 * np.linspace(0.1, 10, 10)
    shear = 50.0e9 * np.linspace(1., 10, 10)
    frequency = 1.0e-6

    for rheology_name, rheology_class in rheology_classes.items():
        rheology_instance = rheology_class()
        result = np.empty(shear.size, dtype=np.complex128)
        rheology_instance.vectorize_modulus_viscosity(frequency, shear, viscosity, result)
        assert type(result) == np.ndarray
        assert result.dtype == np.complex128


def test_find_rheology():
    """ Test the find rheology helper function. """

    for rheology_name, rheology_class in rheology_classes.items():
        test_rheo_class = find_rheology(rheology_name)
        assert test_rheo_class is rheology_class
