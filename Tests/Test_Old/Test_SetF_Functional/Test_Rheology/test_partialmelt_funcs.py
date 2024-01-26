import numpy as np

import TidalPy


solidus = 1600.
liquidus = 2000.
liquid_shear = 1.0e-5
fs_visc_power_slope = 27000.0
fs_visc_power_phase = 1.0
fs_shear_power_slope = 82000.0
fs_shear_power_phase = 40.6
crit_melt_frac = 0.5
crit_melt_frac_width = 0.05
hn_visc_slope_1 = 13.5
hn_visc_slope_2 = 370.0
hn_shear_param_1 = 40000.0
hn_shear_param_2 = 25.0
hn_shear_falloff_slope = 700.0

henning_input = (liquid_shear, crit_melt_frac, crit_melt_frac_width, hn_visc_slope_1, hn_visc_slope_2, hn_shear_param_1,
                 hn_shear_param_2, hn_shear_falloff_slope)
spohn_input = (liquid_shear, fs_visc_power_slope, fs_visc_power_phase, fs_shear_power_slope, fs_shear_power_phase)


def test_melt_fraction():
    from TidalPy.rheology.partial_melt import calculate_melt_fraction, calculate_melt_fraction_array

    # Test float version
    mf = calculate_melt_fraction(1600., 1600., 2000.)
    assert mf == 0.
    mf = calculate_melt_fraction(2000., 1600., 2000.)
    assert mf == 1.
    mf = calculate_melt_fraction(1800., 1600., 2000.)
    assert mf == .5

    # Test array version
    mf = calculate_melt_fraction_array(1600. * np.ones(10), 1600., 2000.)
    np.testing.assert_allclose(mf, 0.)
    mf = calculate_melt_fraction_array(2000. * np.ones(10), 1600., 2000.)
    np.testing.assert_allclose(mf, 1.)
    mf = calculate_melt_fraction_array(1800. * np.ones(10), 1600., 2000.)
    np.testing.assert_allclose(mf, 0.5)


def test_reverse_melt_fraction():
    from TidalPy.rheology.partial_melt import calculate_temperature_frommelt_array, calculate_temperature_frommelt

    # Test float version
    temp = calculate_temperature_frommelt(0., 1600., 2000.)
    assert temp == 1600.
    temp = calculate_temperature_frommelt(1., 1600., 2000.)
    assert temp == 2000.
    temp = calculate_temperature_frommelt(.5, 1600., 2000.)
    assert temp == 1800.

    # Test array version
    temp = calculate_temperature_frommelt_array(0. * np.ones(10), 1600., 2000.)
    np.testing.assert_allclose(temp, 1600.)
    temp = calculate_temperature_frommelt_array(1. * np.ones(10), 1600., 2000.)
    np.testing.assert_allclose(temp, 2000.)
    temp = calculate_temperature_frommelt_array(.5 * np.ones(10), 1600., 2000.)
    np.testing.assert_allclose(temp, 1800.)


def test_all_models():
    # TODO: Currently this does not have a nice programic-way to test `all` models (pull in their inputs, etc.)
    #    This should be possible using the partial_melt module's get_partial_melt_model_default_inputs()

    print('Partial Melting Model Testing...')

    # Test float versions
    from TidalPy.rheology.partial_melt.melting_models import off, spohn, henning
    melt_fraction = 0.5
    temperature = 1800.
    premelt_viscosity = 1.e20
    liquid_viscosity = 1.0
    premelt_shear = 50.e9

    #    Testing `off` model
    off_visco, off_shear = off(melt_fraction, premelt_viscosity, premelt_shear)
    assert type(off_visco) in [float, np.float64]
    assert type(off_shear) in [float, np.float64]
    assert off_visco == premelt_viscosity
    assert off_shear == premelt_shear

    #    Testing `spohn` model
    spohn_visco, spohn_shear = spohn(melt_fraction, temperature, liquid_viscosity, *spohn_input)
    assert type(spohn_visco) in [float, np.float64]
    assert type(spohn_shear) in [float, np.float64]

    #    Testing `henning` model
    henn_visco, henn_shear = henning(
        melt_fraction, temperature, premelt_viscosity, liquid_viscosity, premelt_shear,
        *henning_input
        )
    assert type(henn_visco) in [float, np.float64]
    assert type(henn_shear) in [float, np.float64]

    # Test Array versions
    melt_fraction = 0.5 * np.ones(10)
    temperature = 1800. * np.ones(10)
    premelt_viscosity = 1.e20 * np.ones(10)
    liquid_viscosity = 1.0 * np.ones(10)
    premelt_shear = 50.e9 * np.ones(10)

    #    Testing `off` model
    off_visco, off_shear = off(melt_fraction, premelt_viscosity, premelt_shear)
    assert off_visco.shape == melt_fraction.shape
    assert off_shear.shape == melt_fraction.shape
    np.testing.assert_allclose(off_visco, premelt_viscosity)
    np.testing.assert_allclose(off_shear, premelt_shear)

    #    Testing `spohn` model
    spohn_visco, spohn_shear = spohn(melt_fraction, temperature, liquid_viscosity, *spohn_input)
    assert spohn_visco.shape == melt_fraction.shape
    assert spohn_shear.shape == melt_fraction.shape

    #    Testing `henning` model
    henn_visco, henn_shear = henning(
        melt_fraction, temperature, premelt_viscosity, liquid_viscosity,
        premelt_shear, *henning_input
        )
    assert henn_visco.shape == melt_fraction.shape
    assert henn_shear.shape == melt_fraction.shape
