import numpy as np
import numba

import TidalPy
TidalPy.config['stream_level'] = 'WARNING'
TidalPy.reinit()

from TidalPy.rheology.complex_compliance import known_models

def test_all_models():

    print('Complex Compliance Model Testing...')

    for model_name, model_func in known_models.items():

        print(f'\t Working on {model_name}')

        if '_array' in model_name:
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

        else:
            # Non array test
            res = model_func(0.001, 1.0, 2.0)
            assert type(res) in [complex, np.complex]




