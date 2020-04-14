import numpy as np

from TidalPy.rheology.compliance import known_models

def test_all_models():

    print('Complex Compliance Model Testing...')

    for model_name, model_func in known_models.items():

        print(f'\t Working on {model_name}')

        # Test Frequency Array
        compliance, viscosity = np.asarray(1.), np.asarray(1.)
        frequency = np.linspace(0., 10., 10)
        res = model_func(compliance, viscosity, frequency)
        assert res.shape == frequency.shape

        # Test Compliance/Viscosity Array
        compliance, viscosity = np.linspace(0., 10., 10), np.linspace(0., 10., 10)
        frequency = np.asarray(1.)
        res = model_func(compliance, viscosity, frequency)
        assert res.shape == compliance.shape

        # Test Mix Array
        xx, yy = np.meshgrid(np.linspace(0., 10., 10), np.linspace(0., 10., 20))
        compliance = xx
        viscosity = np.copy(xx)
        frequency = yy
        x = model_func(compliance, viscosity, frequency)
        assert res.shape == xx.shape




