import numpy as np

from TidalPy.rheology.partialMelt import known_models


def test_all_models():

    print('Partial Melting Model Testing...')

    for model_name, model_func in known_models.items():

        print(f'\t Working on {model_name}')

        # Temperature, compliance, and viscosity should always have the same shape.
        #     So we can just do a single array testing

        # Single Value Test
        premelt_shear = 1.
        melt_fraction = np.asarray(1.)
        temperature, premelt_viscosity, liquid_viscosity = np.asarray(1.), np.asarray(1.), np.asarray(1.)

        postmelt_visc, postmelt_comp = \
            model_func(temperature, melt_fraction, premelt_viscosity, liquid_viscosity, premelt_shear)
        assert postmelt_visc.shape == melt_fraction.shape
        assert postmelt_comp.shape == melt_fraction.shape

        # Array Value Test
        premelt_shear = 1.
        melt_fraction = np.linspace(0., 1., 11)
        temperature, premelt_viscosity, liquid_viscosity = \
            np.linspace(10., 11., 11), np.linspace(10., 11., 11), np.linspace(10., 11., 11)

        postmelt_visc, postmelt_comp = \
            model_func(temperature, melt_fraction, premelt_viscosity, liquid_viscosity, premelt_shear)
        assert postmelt_visc.shape == melt_fraction.shape
        assert postmelt_comp.shape == melt_fraction.shape