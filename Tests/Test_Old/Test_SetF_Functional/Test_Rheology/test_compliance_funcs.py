import numba
import numpy as np

import TidalPy


from TidalPy.rheology.complex_compliance import known_models, known_model_const_args

comp_defaults = TidalPy.config['layers']['rock']['rheology']


def test_all_models():
    print('Complex Compliance Model Testing...')

    for model_name, model_func in known_models.items():

        print(f'\t Working on {model_name}')

        # Test Frequency Array
        compliance, viscosity = np.asarray(1.), np.asarray(1.)
        frequency = np.linspace(0., 10., 10)

        res = model_func(frequency, compliance, viscosity)
        assert res.shape == frequency.shape

        # Test Compliance/Viscosity Array
        compliance, viscosity = np.linspace(0., 10., 10), np.linspace(0., 10., 10)

        # TODO: The way frequency=0 checks are implemented requires that frequency always have the same shape
        #    as other input arrays.
        # frequency = np.asarray(1.)
        frequency = np.ones_like(viscosity)
        res = model_func(frequency, compliance, viscosity)
        assert res.shape == compliance.shape

        # Test Mix Array
        xx, yy = np.meshgrid(np.linspace(0., 10., 10), np.linspace(0., 10., 20))
        compliance = xx
        viscosity = np.copy(xx)
        frequency = yy
        try:
            res = model_func(frequency, compliance, viscosity)
            assert res.shape == compliance.shape
        except numba.errors.TypingError:
            # TODO: This is an expected (but unfortunate) error wherein numba does not support 2D+ arrays.
            #    Current workaround is to use flatten() on the 2d array but save the shape.
            shape = xx.shape
            compliance = compliance.flatten()
            viscosity = viscosity.flatten()
            frequency = frequency.flatten()
            res = model_func(frequency, compliance, viscosity)
            # Return to original shape
            res = res.reshape(shape)
            compliance = compliance.reshape(shape)
            assert res.shape == compliance.shape

        # Non array test
        res = model_func(0.001, 1.0, 2.0)
        assert type(res) in [complex, np.complex128]

def test_complex_compliance_zero_freq():
    """ Test that all complex compliance functions return the expected result for when frequency is zero"""

    # Test floats
    frequency = 0.
    compliance = (50.e9)**(-1)
    viscosity = 1.0e18
    for model_name, model_func in known_models.items():

        # Get the additional inputs for this compliance function
        comp_inputs = tuple([comp_defaults[comp_input] for comp_input in known_model_const_args[model_name]])

        if model_name == 'fixed_q':
            complex_compliance = model_func(frequency, compliance, viscosity, planet_beta=10., quality_factor=100.)
            test_value_real = -19. / (2. * 10.)
        elif model_name in ('voigt',):
            complex_compliance = model_func(frequency, compliance, viscosity, *comp_inputs)
            voigt_comp_offset = comp_inputs[0]
            test_value_real = voigt_comp_offset * compliance
        elif model_name in ('burgers',):
            complex_compliance = model_func(frequency, compliance, viscosity, *comp_inputs)
            voigt_comp_offset = comp_inputs[0]
            test_value_real = (1. + voigt_comp_offset) * compliance
        elif model_name in ('andrade', 'sundberg', 'andrade_freq', 'sundberg_freq'):
            complex_compliance = model_func(frequency, compliance, viscosity, *comp_inputs)
            test_value_real = 1.0e100
        elif model_name in ('newton',):
            complex_compliance = model_func(frequency, compliance, viscosity, *comp_inputs)
            test_value_real = 0.
        elif model_name in ('elastic',):
            complex_compliance = model_func(frequency, compliance, viscosity, *comp_inputs)
            test_value_real = compliance
        else:
            complex_compliance = model_func(frequency, compliance, viscosity, *comp_inputs)
            test_value_real = compliance

        # For zero frequency the real part of the complex compliance should equal the regular compliance. Imag = 0.
        assert np.real(complex_compliance) == test_value_real
        assert np.imag(complex_compliance) == 0.

    # Test arrays - comp / visc array
    frequency = 0.
    compliance = np.linspace(30.e9, 50.e9, 4, dtype=np.float64)
    viscosity = np.logspace(22., 18., 4, dtype=np.float64)
    for model_name, model_func in known_models.items():

        # Get the additional inputs for this compliance function
        comp_inputs = tuple([comp_defaults[comp_input] for comp_input in known_model_const_args[model_name]])

        if model_name == 'fixed_q':
            complex_compliance = model_func(frequency, compliance, viscosity, planet_beta=10., quality_factor=100.)
            test_value_real = -19. / (2. * 10.) * np.ones_like(compliance)
        elif model_name in ('voigt',):
            complex_compliance = model_func(frequency, compliance, viscosity, *comp_inputs)
            voigt_comp_offset = comp_inputs[0]
            test_value_real = voigt_comp_offset * compliance
        elif model_name in ('burgers',):
            complex_compliance = model_func(frequency, compliance, viscosity, *comp_inputs)
            voigt_comp_offset = comp_inputs[0]
            test_value_real = (1. + voigt_comp_offset) * compliance
        elif model_name in ('andrade', 'sundberg', 'andrade_freq', 'sundberg_freq',):
            complex_compliance = model_func(frequency, compliance, viscosity, *comp_inputs)
            test_value_real = 1.0e100 * np.ones_like(compliance)
        elif model_name in ('newton',):
            complex_compliance = model_func(frequency, compliance, viscosity, *comp_inputs)
            test_value_real = 0.e100 * np.ones_like(compliance)
        elif model_name in ('elastic',):
            complex_compliance = model_func(frequency, compliance, viscosity, *comp_inputs)
            test_value_real = compliance
        else:
            complex_compliance = model_func(frequency, compliance, viscosity, *comp_inputs)
            test_value_real = compliance

        # For zero frequency the real part of the complex compliance should equal the regular compliance. Imag = 0.
        assert np.allclose(np.real(complex_compliance), test_value_real)
        assert np.all(np.imag(complex_compliance) == 0.)

    # Test arrays - freq array
    frequency = np.asarray((1.e-3, 1.e-6, 0., -1.e-6, -1.e-3))
    compliance = 50.e9
    viscosity = 1.e18
    for model_name, model_func in known_models.items():

        # Get the additional inputs for this compliance function
        comp_inputs = tuple([comp_defaults[comp_input] for comp_input in known_model_const_args[model_name]])

        if model_name == 'fixed_q':
            complex_compliance = model_func(frequency, compliance, viscosity, planet_beta=10., quality_factor=100.)
            test_value_real = -19. / (2. * 10.)
        elif model_name in ('voigt',):
            complex_compliance = model_func(frequency, compliance, viscosity, *comp_inputs)
            voigt_comp_offset = comp_inputs[0]
            test_value_real = voigt_comp_offset * compliance
        elif model_name in ('burgers',):
            complex_compliance = model_func(frequency, compliance, viscosity, *comp_inputs)
            voigt_comp_offset = comp_inputs[0]
            test_value_real = (1. + voigt_comp_offset) * compliance
        elif model_name in ('andrade', 'sundberg', 'andrade_freq', 'sundberg_freq'):
            complex_compliance = model_func(frequency, compliance, viscosity, *comp_inputs)
            test_value_real = 1.0e100
        elif model_name in ('newton'):
            complex_compliance = model_func(frequency, compliance, viscosity, *comp_inputs)
            test_value_real = 0.0
        elif model_name in ('elastic'):
            complex_compliance = model_func(frequency, compliance, viscosity, *comp_inputs)
            test_value_real = compliance
        else:
            complex_compliance = model_func(frequency, compliance, viscosity, *comp_inputs)
            test_value_real = compliance

        # For zero frequency the real part of the complex compliance should equal the regular compliance. Imag = 0.
        assert np.real(complex_compliance[frequency==0.]) == test_value_real
        assert np.imag(complex_compliance[frequency==0.]) == 0.


