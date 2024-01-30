import numpy as np

import TidalPy


from TidalPy.structures import build_from_world, build_world
from TidalPy.structures.orbit import PhysicsOrbit

star_config_ = {
    'tides_on': False,
    }
star_base = build_world('55cnc')

config_ = {
    'force_spin_sync': False,
    'type'           : 'layered',
    "tides"          : {
        "model"                      : "layered",
        'eccentricity_truncation_lvl': 2,
        'max_tidal_order_l'          : 2,
        'obliquity_tides_on'         : True
        },
    "layers"         : {
        "Core"  : {
            "is_tidal": False,
            "rheology": {
                'complex_compliance': {
                    'model': 'maxwell'
                    }
                }
            },
        "Mantle": {
            "is_tidal": True
            }
        }
    }
io_base = build_world('io_simple')

def test_single_body_derivative_calc():
    """ This test will load a tidally active body and a tidally inactive host, set the temperatures and orbit then
        check that single body tidal calculations are occurring as expected. """

    star_to_use = build_from_world(star_base, new_config=star_config_)
    world_to_use = build_from_world(io_base, new_config=config_, new_name='io_world')

    orbit = PhysicsOrbit(star_to_use, tidal_host=star_to_use, tidal_bodies=world_to_use)

    # Test loading in layer temperature and performing necessary steps to eventually calculate the complex compliance
    #    once a frequency is set.
    world_to_use.mantle.temperature = 1200.
    assert type(world_to_use.mantle.viscosity) in [float, np.float64]
    assert type(world_to_use.mantle.shear_modulus) in [float, np.float64]

    # Test Floats
    # Set orbital frequency and eccentricity - check that spin locking worked.
    world_to_use.set_state(orbital_period=50., eccentricity=0.2, obliquity=np.radians(10.), spin_period=10.)
    # For a world in synchronous rotation the spin period should equal the new orbital period.
    assert not world_to_use.is_spin_sync
    assert world_to_use.tides_on
    assert world_to_use.eccentricity_truncation_lvl == 2
    assert world_to_use.max_tidal_order_lvl == 2
    assert type(world_to_use.tidal_susceptibility) in [float, np.float64]
    assert type(world_to_use.global_love_by_orderl[2]) in [complex, np.complex128]
    assert type(world_to_use.global_negative_imk_by_orderl[2]) in [float, np.float64]
    assert type(world_to_use.effective_q_by_orderl[2]) in [float, np.float64]
    assert type(world_to_use.dUdM) in [float, np.float64]
    assert type(world_to_use.dUdO) in [float, np.float64]
    assert type(world_to_use.dUdw) in [float, np.float64]
    assert type(world_to_use.tidal_heating_global) in [float, np.float64]
    # Check that dual body calculations were done.
    assert orbit._last_calc_used_dual_body == False
    assert orbit.get_orbital_motion_time_derivative(world_to_use) is \
           orbit.get_orbital_motion_time_derivative(star_to_use)
    assert orbit.get_eccentricity_time_derivative(world_to_use) is \
           orbit.get_eccentricity_time_derivative(star_to_use)
    assert orbit.get_semi_major_axis_time_derivative(world_to_use) is \
           orbit.get_semi_major_axis_time_derivative(star_to_use)
    # Check that types are as expected
    assert type(orbit.get_orbital_motion_time_derivative(world_to_use)) in [float, np.float64]
    assert type(orbit.get_eccentricity_time_derivative(world_to_use)) in [float, np.float64]
    assert type(orbit.get_semi_major_axis_time_derivative(world_to_use)) in [float, np.float64]

    # Test arrays
    # Set orbital frequency and eccentricity - check that spin locking worked.
    orb_period = np.linspace(10., 50., 10)
    world_to_use.set_state(
        orbital_period=orb_period, eccentricity=0.2, obliquity=np.radians(10.),
        spin_period=0.5 * orb_period
        )
    # For a world in synchronous rotation the spin period should equal the new orbital period.
    assert not world_to_use.is_spin_sync
    assert world_to_use.tides_on
    assert world_to_use.eccentricity_truncation_lvl == 2
    assert world_to_use.max_tidal_order_lvl == 2
    assert type(world_to_use.tidal_susceptibility) == np.ndarray
    assert type(world_to_use.global_love_by_orderl[2]) == np.ndarray
    assert type(world_to_use.global_negative_imk_by_orderl[2]) == np.ndarray
    assert type(world_to_use.effective_q_by_orderl[2]) == np.ndarray
    assert type(world_to_use.dUdM) == np.ndarray
    assert type(world_to_use.dUdO) == np.ndarray
    assert type(world_to_use.dUdw) == np.ndarray
    assert type(world_to_use.tidal_heating_global) == np.ndarray
    # Check that dual body calculations were done.
    assert orbit._last_calc_used_dual_body == False
    # Check that types are as expected
    assert type(orbit.get_orbital_motion_time_derivative(world_to_use)) == np.ndarray
    assert type(orbit.get_eccentricity_time_derivative(world_to_use)) == np.ndarray
    assert type(orbit.get_semi_major_axis_time_derivative(world_to_use)) == np.ndarray


def test_dual_body_derivative_calc():
    """ This test will load a tidally active host and body, set the temperatures and orbits then check that dual body
        tidal calculations are occurring as expected. """

    star_to_use = build_from_world(star_base, new_config=star_config_)
    host_to_use = build_from_world(io_base, new_config=config_, new_name='io_host')
    world_to_use = build_from_world(io_base, new_config=config_, new_name='io_world')
    orbit = PhysicsOrbit(star_to_use, tidal_host=host_to_use, tidal_bodies=world_to_use)

    # Test loading in layer temperature and performing necessary steps to eventually calculate the complex compliance
    #    once a frequency is set.
    host_to_use.mantle.temperature = 1400.
    world_to_use.mantle.temperature = 1200.
    assert type(world_to_use.mantle.viscosity) in [float, np.float64]
    assert type(world_to_use.mantle.shear_modulus) in [float, np.float64]
    assert type(host_to_use.mantle.viscosity) in [float, np.float64]
    assert type(host_to_use.mantle.shear_modulus) in [float, np.float64]

    # Test Floats
    # Set orbital frequency and eccentricity - check that spin locking worked.
    world_to_use.set_state(orbital_period=50., eccentricity=0.2, obliquity=np.radians(10.), spin_period=10.)
    host_to_use.set_state(obliquity=np.radians(5.), spin_period=8.)
    # For a world in synchronous rotation the spin period should equal the new orbital period.
    for world_to_check in [host_to_use, world_to_use]:
        assert not world_to_check.is_spin_sync
        assert world_to_check.tides_on
        assert world_to_check.eccentricity_truncation_lvl == 2
        assert world_to_check.max_tidal_order_lvl == 2
        assert type(world_to_check.tidal_susceptibility) in [float, np.float64]
        assert type(world_to_check.global_love_by_orderl[2]) in [complex, np.complex128]
        assert type(world_to_check.global_negative_imk_by_orderl[2]) in [float, np.float64]
        assert type(world_to_check.effective_q_by_orderl[2]) in [float, np.float64]
        assert type(world_to_check.dUdM) in [float, np.float64]
        assert type(world_to_check.dUdO) in [float, np.float64]
        assert type(world_to_check.dUdw) in [float, np.float64]
        assert type(world_to_check.tidal_heating_global) in [float, np.float64]
    # Check that dual body calculations were done.
    assert orbit._last_calc_used_dual_body == True
    assert orbit.get_orbital_motion_time_derivative(world_to_use) is \
           orbit.get_orbital_motion_time_derivative(host_to_use)
    assert orbit.get_eccentricity_time_derivative(world_to_use) is \
           orbit.get_eccentricity_time_derivative(host_to_use)
    assert orbit.get_semi_major_axis_time_derivative(world_to_use) is \
           orbit.get_semi_major_axis_time_derivative(host_to_use)
    # Check that types are as expected
    assert type(orbit.get_orbital_motion_time_derivative(world_to_use)) in [float, np.float64]
    assert type(orbit.get_eccentricity_time_derivative(world_to_use)) in [float, np.float64]
    assert type(orbit.get_semi_major_axis_time_derivative(world_to_use)) in [float, np.float64]

    # Test arrays
    # Set orbital frequency and eccentricity - check that spin locking worked.
    orb_period = np.linspace(10., 50., 10)
    world_to_use.set_state(
        orbital_period=orb_period, eccentricity=0.2, obliquity=np.radians(10.),
        spin_period=0.5 * orb_period
        )
    host_to_use.set_state(obliquity=np.radians(5.), spin_period=0.2 * orb_period)
    # For a world in synchronous rotation the spin period should equal the new orbital period.
    for world_to_check in [host_to_use, world_to_use]:
        assert not world_to_check.is_spin_sync
        assert world_to_check.tides_on
        assert world_to_check.eccentricity_truncation_lvl == 2
        assert world_to_check.max_tidal_order_lvl == 2
        assert type(world_to_check.tidal_susceptibility) == np.ndarray
        assert type(world_to_check.global_love_by_orderl[2]) == np.ndarray
        assert type(world_to_check.global_negative_imk_by_orderl[2]) == np.ndarray
        assert type(world_to_check.effective_q_by_orderl[2]) == np.ndarray
        assert type(world_to_check.dUdM) == np.ndarray
        assert type(world_to_check.dUdO) == np.ndarray
        assert type(world_to_check.dUdw) == np.ndarray
        assert type(world_to_check.tidal_heating_global) == np.ndarray
    # Check that dual body calculations were done.
    assert orbit._last_calc_used_dual_body == True
    # Check that types are as expected
    assert type(orbit.get_orbital_motion_time_derivative(world_to_use)) == np.ndarray
    assert type(orbit.get_eccentricity_time_derivative(world_to_use)) == np.ndarray
    assert type(orbit.get_semi_major_axis_time_derivative(world_to_use)) == np.ndarray


def test_dual_body_derivative_calc_temperature_array():
    """ This test will load a tidally active host and body, set the temperatures and orbits then check that dual body
        tidal calculations are occurring as expected. Temperature property will be set as an array in the array test"""

    star_to_use = build_from_world(star_base, new_config=star_config_)
    host_to_use = build_from_world(io_base, new_config=config_, new_name='io_host')
    world_to_use = build_from_world(io_base, new_config=config_, new_name='io_world')
    orbit = PhysicsOrbit(star_to_use, tidal_host=host_to_use, tidal_bodies=world_to_use)

    # Test loading in layer temperature and performing necessary steps to eventually calculate the complex compliance
    #    once a frequency is set.
    host_to_use.mantle.temperature = 1400.
    world_to_use.mantle.temperature = 1200.
    assert type(world_to_use.mantle.viscosity) in [float, np.float64]
    assert type(world_to_use.mantle.shear_modulus) in [float, np.float64]
    assert type(host_to_use.mantle.viscosity) in [float, np.float64]
    assert type(host_to_use.mantle.shear_modulus) in [float, np.float64]

    # Test Floats
    # Set orbital frequency and eccentricity - check that spin locking worked.
    world_to_use.set_state(orbital_period=50., eccentricity=0.2, obliquity=np.radians(10.), spin_period=10.)
    host_to_use.set_state(obliquity=np.radians(5.), spin_period=8.)
    # For a world in synchronous rotation the spin period should equal the new orbital period.
    for world_to_check in [host_to_use, world_to_use]:
        assert not world_to_check.is_spin_sync
        assert world_to_check.tides_on
        assert world_to_check.eccentricity_truncation_lvl == 2
        assert world_to_check.max_tidal_order_lvl == 2
        assert type(world_to_check.tidal_susceptibility) in [float, np.float64]
        assert type(world_to_check.global_love_by_orderl[2]) in [complex, np.complex128]
        assert type(world_to_check.global_negative_imk_by_orderl[2]) in [float, np.float64]
        assert type(world_to_check.effective_q_by_orderl[2]) in [float, np.float64]
        assert type(world_to_check.dUdM) in [float, np.float64]
        assert type(world_to_check.dUdO) in [float, np.float64]
        assert type(world_to_check.dUdw) in [float, np.float64]
        assert type(world_to_check.tidal_heating_global) in [float, np.float64]
    # Check that dual body calculations were done.
    assert orbit._last_calc_used_dual_body == True
    assert orbit.get_orbital_motion_time_derivative(world_to_use) is \
           orbit.get_orbital_motion_time_derivative(host_to_use)
    assert orbit.get_eccentricity_time_derivative(world_to_use) is \
           orbit.get_eccentricity_time_derivative(host_to_use)
    assert orbit.get_semi_major_axis_time_derivative(world_to_use) is \
           orbit.get_semi_major_axis_time_derivative(host_to_use)
    # Check that types are as expected
    assert type(orbit.get_orbital_motion_time_derivative(world_to_use)) in [float, np.float64]
    assert type(orbit.get_eccentricity_time_derivative(world_to_use)) in [float, np.float64]
    assert type(orbit.get_semi_major_axis_time_derivative(world_to_use)) in [float, np.float64]

    # Test arrays
    # Set orbital frequency and eccentricity - check that spin locking worked.
    host_to_use.mantle.temperature = np.linspace(1100., 1500., 10)
    world_to_use.mantle.temperature = np.linspace(1000., 1500., 10)
    # For a world in synchronous rotation the spin period should equal the new orbital period.
    for world_to_check in [host_to_use, world_to_use]:
        assert not world_to_check.is_spin_sync
        assert world_to_check.tides_on
        assert world_to_check.eccentricity_truncation_lvl == 2
        assert world_to_check.max_tidal_order_lvl == 2
        assert type(world_to_check.tidal_susceptibility) in [float, np.float64]
        assert type(world_to_check.global_love_by_orderl[2]) == np.ndarray
        assert type(world_to_check.global_negative_imk_by_orderl[2]) == np.ndarray
        assert type(world_to_check.effective_q_by_orderl[2]) == np.ndarray
        assert type(world_to_check.dUdM) == np.ndarray
        assert type(world_to_check.dUdO) == np.ndarray
        assert type(world_to_check.dUdw) == np.ndarray
        assert type(world_to_check.tidal_heating_global) == np.ndarray
    # Check that dual body calculations were done.
    assert orbit._last_calc_used_dual_body == True
    # Check that types are as expected
    assert type(orbit.get_orbital_motion_time_derivative(world_to_use)) == np.ndarray
    assert type(orbit.get_eccentricity_time_derivative(world_to_use)) == np.ndarray
    assert type(orbit.get_semi_major_axis_time_derivative(world_to_use)) == np.ndarray
