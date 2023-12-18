import numpy as np

import TidalPy


from TidalPy.structures import build_world, build_from_world

# Build basic layered io
io_base = build_world('io_simple')


def test_radiogenic_calc_fixed_in_layered_world():
    """ This test will build a LayeredWorld and see if radiogenic calculations are performed as expected """

    # Test a model where the radiogenics are off
    radio_dict = {'layers': {'Mantle': {'is_tidal': False, 'radiogenics': {'model': 'off'}}}}
    io_off = build_from_world(io_base, radio_dict)

    # Check that the model was loaded correctly
    assert io_off.Mantle.radiogenics.model == 'off'

    # Set the time and make sure that the radiogenic heating is set as expected
    io_off.set_state(time=0.)
    assert io_off.Mantle.radiogenics.heating == 0.

    # Test a model where the radiogenics are on but fixed
    radio_dict = {
        'layers': {
            'Mantle': {
                'is_tidal': False, 'radiogenics': {
                    'model'                : 'fixed',
                    'fixed_heat_production': 100.,
                    'average_half_life'    : 1.,
                    'ref_time'             : 10.
                    }
                }
            }
        }
    io_fixed = build_from_world(io_base, radio_dict)

    # Set time to reference time. The heat production should equal that in the dict.
    io_fixed.set_state(time=10.)
    expected_heating = 100. * io_fixed.Mantle.mass
    np.testing.assert_approx_equal(io_fixed.Mantle.radiogenics.heating, expected_heating)

    # Change the time to an early point
    io_fixed.set_state(time=1.)
    gamma = np.log(0.5) / 1.
    expected_heating = io_fixed.Mantle.mass * 100. * np.exp(gamma * (1. - 10.))
    np.testing.assert_approx_equal(io_fixed.Mantle.radiogenics.heating, expected_heating)

    # Change the time to a later point
    io_fixed.set_state(time=100.)
    gamma = np.log(0.5) / 1.
    expected_heating = io_fixed.Mantle.mass * 100. * np.exp(gamma * (100. - 10.))
    np.testing.assert_approx_equal(io_fixed.Mantle.radiogenics.heating, expected_heating)

    # Check if arrays work
    time = np.linspace(1., 100., 10)
    io_fixed.set_state(time=time)
    gamma = np.log(0.5) / 1.
    expected_heating = io_fixed.Mantle.mass * 100. * np.exp(gamma * (time - 10.))
    np.testing.assert_allclose(io_fixed.Mantle.radiogenics.heating, expected_heating)


def test_radiogenic_calc_isotope_in_layered_world():
    """ This test will build a LayeredWorld and see if radiogenic calculations are performed as expected """

    # Test a model where there are some fake isotopes
    radio_dict = {
        'layers': {
            'Mantle': {
                'is_tidal': False, 'radiogenics': {
                    'model'   : 'isotope',
                    'ref_time': 10.,
                    'isotopes': {
                        'iso1': {
                            'iso_mass_fraction'    : 0.001,
                            'hpr'                  : 0.5,
                            'half_life'            : 1.,
                            'element_concentration': 0.002
                            },
                        'iso2': {
                            'iso_mass_fraction'    : 0.005,
                            'hpr'                  : 0.75,
                            'half_life'            : .1,
                            'element_concentration': 0.002
                            }
                        }
                    }
                }
            }
        }
    io_iso = build_from_world(io_base, radio_dict)
    assert io_iso.Mantle.radiogenics.model == 'isotope'

    def expected_heating_func(time):
        specific_heating = np.zeros_like(time)
        specific_heating += .001 * .002 * .5 * np.exp((np.log(0.5) / 1.) * (time - 10.))
        specific_heating += .005 * .002 * .75 * np.exp((np.log(0.5) / .1) * (time - 10.))
        return specific_heating * io_iso.Mantle.mass

    # Set time to reference time. The heat production should equal that in the dict.
    io_iso.set_state(time=10.)
    expected_heating = expected_heating_func(10.)
    np.testing.assert_approx_equal(io_iso.Mantle.radiogenics.heating, expected_heating)

    # Change the time to an early point
    io_iso.set_state(time=1.)
    expected_heating = expected_heating_func(1.)
    np.testing.assert_approx_equal(io_iso.Mantle.radiogenics.heating, expected_heating)

    # Change the time to a later point
    io_iso.set_state(time=100.)
    expected_heating = expected_heating_func(100.)
    np.testing.assert_approx_equal(io_iso.Mantle.radiogenics.heating, expected_heating)

    # Check if arrays work
    time = np.linspace(1., 100., 10)
    io_iso.set_state(time=time)
    expected_heating = expected_heating_func(time)
    np.testing.assert_allclose(io_iso.Mantle.radiogenics.heating, expected_heating)
