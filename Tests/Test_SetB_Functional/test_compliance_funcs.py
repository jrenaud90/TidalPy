import numpy as np

from TidalPy.rheology import complex_compliances

def test_all_models():

    for model_func in complex_compliances:

        print(f'Testing {model_func.__name__}')

        # Test Frequency Array
        compliance, viscosity = 1., 1.
        frequency = np.linspace(0., 10., 10)
        res = model_func(compliance, viscosity, frequency)
        assert res.shape == frequency.shape

        # Test Compliance/Viscosity Array
        compliance, viscosity = np.linspace(0., 10., 10), np.linspace(0., 10., 10)
        frequency = 1.
        res = model_func(compliance, viscosity, frequency)
        assert res.shape == compliance.shape

        # Test Mix Array
        xx, yy = np.meshgrid(np.linspace(0., 10., 10), np.linspace(0., 10., 20))
        compliance = xx
        viscosity = np.copy(xx)
        frequency = yy
        x = model_func(compliance, viscosity, frequency)
        assert res.shape == xx.shape




