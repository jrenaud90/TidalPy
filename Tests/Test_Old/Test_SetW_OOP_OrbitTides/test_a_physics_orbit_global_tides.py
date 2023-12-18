import numpy as np

import TidalPy


from TidalPy.structures import build_world, build_from_world
from TidalPy.structures.orbit import PhysicsOrbit
from TidalPy.utilities.conversions import rads2days, days2rads, semi_a2orbital_motion

star = build_world('55cnc')
world = build_world('earth_simple')


def test_apply_physics_orbit():
    """ This will test applying a Physical orbit to the system """

    orbit = PhysicsOrbit(star, tidal_host=star, tidal_bodies=world)

    assert orbit.tidal_host is star
    assert orbit.star is star
    assert orbit.tidal_objects[0] is star
    assert orbit.tidal_objects[1] is world


def test_set_orbital_frequency_and_eccentricity():
    """ This will test applying a Physical orbit to the system """

    world_tidal = build_from_world(world, new_config={'tides_on': True})
    orbit = PhysicsOrbit(star, tidal_host=star, tidal_bodies=world_tidal)

    # Test orbital frequency - Float
    orbit.set_state(world_tidal, orbital_period=50.)
    # Check it was copied over everywhere
    assert orbit.get_orbital_period(world_tidal) == 50.
    assert world_tidal.orbital_period == 50.
    np.testing.assert_approx_equal(rads2days(world_tidal.orbital_frequency), 50.)
    np.testing.assert_allclose(
        rads2days(semi_a2orbital_motion(world_tidal.semi_major_axis, orbit.tidal_host.mass, world_tidal.mass)),
        50.
        )

    # Test orbital frequency - Array
    orb_period = np.linspace(10., 50., 10)
    orbit.set_state(world_tidal, orbital_period=orb_period)
    # Check it was copied over everywhere
    np.testing.assert_allclose(orbit.get_orbital_period(world_tidal), orb_period)
    np.testing.assert_allclose(world_tidal.orbital_period, orb_period)
    np.testing.assert_allclose(rads2days(world_tidal.orbital_frequency), orb_period)
    np.testing.assert_allclose(
        rads2days(semi_a2orbital_motion(world_tidal.semi_major_axis, orbit.tidal_host.mass, world_tidal.mass)),
        orb_period
        )

    # Test eccentricity - Float
    orbit.set_state(world_tidal, eccentricity=0.1)
    # Check it was copied over everywhere
    assert orbit.get_eccentricity(world_tidal) == 0.1
    assert world_tidal.eccentricity == 0.1

    # Test eccentricity - Array
    eccentricities = np.linspace(0.1, 0.3, 10)
    orbit.set_state(world_tidal, eccentricity=eccentricities)
    # Check it was copied over everywhere
    np.testing.assert_allclose(orbit.get_eccentricity(world_tidal), eccentricities)
    np.testing.assert_allclose(world_tidal.eccentricity, eccentricities)


def test_global_tidal_calculation_cpl_synchronous_rotation_no_obliquity():
    """ This will test the global tidal heating calculation assuming a cpl model assuming synchronous rotation
        and no obliquity. """

    config_ = {
        'force_spin_sync': True,
        'type'           : 'simple_tidal',
        'mass'           : 5.972e24,
        # TODO: Bug when you remove the below. For some reason the old world's slice number gets set to None.
        'slices'         : 100,
        "tides"          : {
            "model"                      : "global_approx",
            "fixed_q"                    : 125.0,
            "use_ctl"                    : False,
            'eccentricity_truncation_lvl': 2,
            'max_tidal_order_l'          : 2,
            'obliquity_tides_on'         : False
            }
        }
    world_cpl = build_from_world(world, new_config=config_)
    orbit = PhysicsOrbit(star, tidal_host=star, tidal_bodies=world_cpl)

    # Test Floats
    # Set orbital frequency and eccentricity - check that spin locking worked.
    orbit.set_state(world_cpl, orbital_period=50., eccentricity=0.2)
    # For a world in synchronous rotation the spin period should equal the new orbital period.
    assert world_cpl.is_spin_sync
    assert world_cpl.tides_on
    assert world_cpl.fixed_q == 125.
    assert world_cpl.eccentricity_truncation_lvl == 2
    assert world_cpl.max_tidal_order_lvl == 2
    assert world_cpl.spin_period == 50.
    assert type(world_cpl.tidal_susceptibility) in [float, np.float64]
    assert type(world_cpl.global_love_by_orderl[2]) in [complex, np.complex128]
    assert type(world_cpl.global_negative_imk_by_orderl[2]) in [float, np.float64]
    assert type(world_cpl.effective_q_by_orderl[2]) in [float, np.float64]
    assert type(world_cpl.dUdM) in [float, np.float64]
    assert type(world_cpl.dUdO) in [float, np.float64]
    assert type(world_cpl.dUdw) in [float, np.float64]
    assert type(world_cpl.tidal_heating_global) in [float, np.float64]
    # Check the orbital derivatives
    assert orbit._last_calc_used_dual_body == False
    assert type(orbit.get_eccentricity_time_derivative(world_cpl)) in [float, np.float64]
    assert type(orbit.get_semi_major_axis_time_derivative(world_cpl)) in [float, np.float64]
    assert type(orbit.get_orbital_motion_time_derivative(world_cpl)) in [float, np.float64]
    # Check the tidal body's orbital derivatives match the hosts
    assert orbit.get_eccentricity_time_derivative(world_cpl) is orbit.get_eccentricity_time_derivative(star)
    assert orbit.get_semi_major_axis_time_derivative(world_cpl) is orbit.get_semi_major_axis_time_derivative(star)
    assert orbit.get_orbital_motion_time_derivative(world_cpl) is orbit.get_orbital_motion_time_derivative(star)
    # Check that the derivatives are being stored in the tidal world correctly
    assert orbit.get_eccentricity_time_derivative(world_cpl) is world_cpl.eccentricity_time_derivative
    assert orbit.get_semi_major_axis_time_derivative(world_cpl) is world_cpl.semi_major_axis_time_derivative
    assert orbit.get_orbital_motion_time_derivative(world_cpl) is world_cpl.orbital_motion_time_derivative

    # Test arrays
    # Set orbital frequency and eccentricity - check that spin locking worked.
    orb_period = np.linspace(10., 50., 10)
    orbit.set_state(world_cpl, orbital_period=orb_period, eccentricity=0.2)
    # For a world in synchronous rotation the spin period should equal the new orbital period.
    assert world_cpl.is_spin_sync
    assert world_cpl.tides_on
    assert world_cpl.fixed_q == 125.
    assert world_cpl.eccentricity_truncation_lvl == 2
    assert world_cpl.max_tidal_order_lvl == 2
    np.testing.assert_allclose(world_cpl.spin_period, orb_period)
    assert type(world_cpl.tidal_susceptibility) == np.ndarray
    assert type(world_cpl.global_love_by_orderl[2]) == np.ndarray
    assert type(world_cpl.global_negative_imk_by_orderl[2]) == np.ndarray
    assert type(world_cpl.effective_q_by_orderl[2]) == np.ndarray
    assert type(world_cpl.dUdM) == np.ndarray
    assert type(world_cpl.dUdO) == np.ndarray
    assert type(world_cpl.dUdw) == np.ndarray
    assert type(world_cpl.tidal_heating_global) == np.ndarray
    # Check the orbital derivatives
    assert type(orbit.get_eccentricity_time_derivative(world_cpl)) == np.ndarray
    assert type(orbit.get_semi_major_axis_time_derivative(world_cpl)) == np.ndarray
    assert type(orbit.get_orbital_motion_time_derivative(world_cpl)) == np.ndarray


def test_global_tidal_calculation_cpl_synchronous_rotation_with_obliquity():
    """ This will test the global tidal heating calculation assuming a cpl model assuming synchronous rotation
        with an obliquity """

    config_ = {
        'force_spin_sync': True,
        'type'           : 'simple_tidal',
        'mass'           : 5.972e24,
        'slices'         : 100,
        "tides"          : {
            "model"                      : "global_approx",
            "fixed_q"                    : 125.0,
            "use_ctl"                    : False,
            'eccentricity_truncation_lvl': 2,
            'max_tidal_order_l'          : 2,
            'obliquity_tides_on'         : True
            }
        }

    world_cpl = build_from_world(world, new_config=config_)
    orbit = PhysicsOrbit(star, tidal_host=star, tidal_bodies=world_cpl)

    # Test Floats
    # Set orbital frequency and eccentricity - check that spin locking worked.
    world_cpl.set_state(orbital_period=50., eccentricity=0.2, obliquity=np.radians(10.))
    # For a world in synchronous rotation the spin period should equal the new orbital period.
    assert world_cpl.is_spin_sync
    assert world_cpl.tides_on
    assert world_cpl.fixed_q == 125.
    assert world_cpl.eccentricity_truncation_lvl == 2
    assert world_cpl.max_tidal_order_lvl == 2
    assert world_cpl.spin_period == 50.
    assert type(world_cpl.tidal_susceptibility) in [float, np.float64]
    assert type(world_cpl.global_love_by_orderl[2]) in [complex, np.complex128]
    assert type(world_cpl.global_negative_imk_by_orderl[2]) in [float, np.float64]
    assert type(world_cpl.effective_q_by_orderl[2]) in [float, np.float64]
    assert type(world_cpl.dUdM) in [float, np.float64]
    assert type(world_cpl.dUdO) in [float, np.float64]
    assert type(world_cpl.dUdw) in [float, np.float64]
    assert type(world_cpl.tidal_heating_global) in [float, np.float64]
    # Check the orbital derivatives
    assert type(orbit.get_eccentricity_time_derivative(world_cpl)) in [float, np.float64]
    assert type(orbit.get_semi_major_axis_time_derivative(world_cpl)) in [float, np.float64]
    assert type(orbit.get_orbital_motion_time_derivative(world_cpl)) in [float, np.float64]

    # Test arrays
    # Set orbital frequency and eccentricity - check that spin locking worked.
    orb_period = np.linspace(10., 50., 10)
    world_cpl.set_state(orbital_period=orb_period, eccentricity=0.2, obliquity=np.radians(10.))
    # For a world in synchronous rotation the spin period should equal the new orbital period.
    assert world_cpl.is_spin_sync
    assert world_cpl.tides_on
    assert world_cpl.fixed_q == 125.
    assert world_cpl.eccentricity_truncation_lvl == 2
    assert world_cpl.max_tidal_order_lvl == 2
    np.testing.assert_allclose(world_cpl.spin_period, orb_period)
    assert type(world_cpl.tidal_susceptibility) == np.ndarray
    assert type(world_cpl.global_love_by_orderl[2]) == np.ndarray
    assert type(world_cpl.global_negative_imk_by_orderl[2]) == np.ndarray
    assert type(world_cpl.effective_q_by_orderl[2]) == np.ndarray
    assert type(world_cpl.dUdM) == np.ndarray
    assert type(world_cpl.dUdO) == np.ndarray
    assert type(world_cpl.dUdw) == np.ndarray
    assert type(world_cpl.tidal_heating_global) == np.ndarray
    # Check the orbital derivatives
    assert type(orbit.get_eccentricity_time_derivative(world_cpl)) == np.ndarray
    assert type(orbit.get_semi_major_axis_time_derivative(world_cpl)) == np.ndarray
    assert type(orbit.get_orbital_motion_time_derivative(world_cpl)) == np.ndarray


def test_global_tidal_calculation_cpl_synchronous_rotation_with_eccen_oblique_array():
    """ This will test the global tidal heating calculation assuming a cpl model assuming synchronous rotation
        with an array version of eccentricity and obliquity """

    config_ = {
        'force_spin_sync': True,
        'type'           : 'simple_tidal',
        'mass'           : 5.972e24,
        'slices'         : 100,
        "tides"          : {
            "model"                      : "global_approx",
            "fixed_q"                    : 125.0,
            "use_ctl"                    : False,
            'eccentricity_truncation_lvl': 2,
            'max_tidal_order_l'          : 2,
            'obliquity_tides_on'         : True
            }
        }

    world_cpl = build_from_world(world, new_config=config_)
    orbit = PhysicsOrbit(star, tidal_host=star, tidal_bodies=world_cpl)

    # Test Floats
    # Set orbital frequency and eccentricity - check that spin locking worked.
    eccen = np.linspace(0.1, 0.3, 10)
    obliqu = np.radians(np.linspace(0., 10., 10))
    world_cpl.set_state(orbital_period=50., eccentricity=eccen, obliquity=obliqu)
    # For a world in synchronous rotation the spin period should equal the new orbital period.
    assert world_cpl.is_spin_sync
    assert world_cpl.tides_on
    assert world_cpl.fixed_q == 125.
    assert world_cpl.eccentricity_truncation_lvl == 2
    assert world_cpl.max_tidal_order_lvl == 2
    assert world_cpl.spin_period == 50.
    assert type(world_cpl.tidal_susceptibility) in [float, np.float64]
    assert type(world_cpl.global_love_by_orderl[2]) in [complex, np.complex128]
    assert type(world_cpl.global_negative_imk_by_orderl[2]) in [float, np.float64]
    assert type(world_cpl.effective_q_by_orderl[2]) in [float, np.float64]
    assert type(world_cpl.dUdM) == np.ndarray
    assert type(world_cpl.dUdO) == np.ndarray
    assert type(world_cpl.dUdw) == np.ndarray
    assert type(world_cpl.tidal_heating_global) == np.ndarray
    # Check the orbital derivatives
    assert type(orbit.get_eccentricity_time_derivative(world_cpl)) == np.ndarray
    assert type(orbit.get_semi_major_axis_time_derivative(world_cpl)) == np.ndarray
    assert type(orbit.get_orbital_motion_time_derivative(world_cpl)) == np.ndarray


def test_global_tidal_calculation_ctl_synchronous_rotation_with_obliquity():
    """ This will test the global tidal heating calculation assuming a ctl model assuming synchronous rotation
        with an obliquity """

    config_ = {
        'force_spin_sync': True,
        'type'           : 'simple_tidal',
        'mass'           : 5.972e24,
        'slices'         : 100,
        "tides"          : {
            "model"                      : "global_approx",
            "use_ctl"                    : True,
            'eccentricity_truncation_lvl': 2,
            'max_tidal_order_l'          : 2,
            'obliquity_tides_on'         : True,
            'ctl_calc_method'            : 'linear_simple',
            'fixed_dt'                   : (1. / 125.) * days2rads(50.)**(-1),
            }
        }

    world_ctl = build_from_world(world, new_config=config_)
    orbit = PhysicsOrbit(star, tidal_host=star, tidal_bodies=world_ctl)

    # Test Floats
    # Set orbital frequency and eccentricity - check that spin locking worked.
    world_ctl.set_state(orbital_period=50., eccentricity=0.2, obliquity=np.radians(10.))
    # For a world in synchronous rotation the spin period should equal the new orbital period.
    assert world_ctl.is_spin_sync
    assert world_ctl.tides_on
    assert world_ctl.eccentricity_truncation_lvl == 2
    assert world_ctl.max_tidal_order_lvl == 2
    assert world_ctl.spin_period == 50.
    assert type(world_ctl.tidal_susceptibility) in [float, np.float64]
    assert type(world_ctl.global_love_by_orderl[2]) in [complex, np.complex128]
    assert type(world_ctl.global_negative_imk_by_orderl[2]) in [float, np.float64]
    assert type(world_ctl.effective_q_by_orderl[2]) in [float, np.float64]
    assert type(world_ctl.dUdM) in [float, np.float64]
    assert type(world_ctl.dUdO) in [float, np.float64]
    assert type(world_ctl.dUdw) in [float, np.float64]
    assert type(world_ctl.tidal_heating_global) in [float, np.float64]
    # Check the orbital derivatives
    assert type(orbit.get_eccentricity_time_derivative(world_ctl)) in [float, np.float64]
    assert type(orbit.get_semi_major_axis_time_derivative(world_ctl)) in [float, np.float64]
    assert type(orbit.get_orbital_motion_time_derivative(world_ctl)) in [float, np.float64]

    # Test arrays
    # Set orbital frequency and eccentricity - check that spin locking worked.
    orb_period = np.linspace(10., 50., 10)
    world_ctl.set_state(orbital_period=orb_period, eccentricity=0.2, obliquity=np.radians(10.))
    # For a world in synchronous rotation the spin period should equal the new orbital period.
    assert world_ctl.is_spin_sync
    assert world_ctl.tides_on
    assert world_ctl.eccentricity_truncation_lvl == 2
    assert world_ctl.max_tidal_order_lvl == 2
    np.testing.assert_allclose(world_ctl.spin_period, orb_period)
    assert type(world_ctl.tidal_susceptibility) == np.ndarray
    assert type(world_ctl.global_love_by_orderl[2]) == np.ndarray
    assert type(world_ctl.global_negative_imk_by_orderl[2]) == np.ndarray
    assert type(world_ctl.effective_q_by_orderl[2]) == np.ndarray
    assert type(world_ctl.dUdM) == np.ndarray
    assert type(world_ctl.dUdO) == np.ndarray
    assert type(world_ctl.dUdw) == np.ndarray
    assert type(world_ctl.tidal_heating_global) == np.ndarray
    # Check the orbital derivatives
    assert type(orbit.get_eccentricity_time_derivative(world_ctl)) == np.ndarray
    assert type(orbit.get_semi_major_axis_time_derivative(world_ctl)) == np.ndarray
    assert type(orbit.get_orbital_motion_time_derivative(world_ctl)) == np.ndarray


def test_global_tidal_calculation_ctl_nsr_with_obliquity():
    """ This will test the global tidal heating calculation assuming a ctl model assuming NSR
        with an obliquity """

    config_ = {
        'force_spin_sync': False,
        'type'           : 'simple_tidal',
        'mass'           : 5.972e24,
        'slices'         : 100,
        "tides"          : {
            "model"                      : "global_approx",
            "use_ctl"                    : True,
            'eccentricity_truncation_lvl': 2,
            'max_tidal_order_l'          : 2,
            'obliquity_tides_on'         : True,
            'ctl_calc_method'            : 'linear_simple',
            'fixed_dt'                   : (1. / 125.) * days2rads(50.)**(-1),
            }
        }

    world_ctl = build_from_world(world, new_config=config_)
    orbit = PhysicsOrbit(star, tidal_host=star, tidal_bodies=world_ctl)

    # Test Floats
    # Set orbital frequency and eccentricity - check that spin locking worked.
    world_ctl.set_state(
        orbital_period=50., eccentricity=0.2, obliquity=np.radians(10.),
        spin_period=10.
        )
    # For a world in synchronous rotation the spin period should equal the new orbital period.
    assert not world_ctl.is_spin_sync
    assert world_ctl.tides_on
    assert world_ctl.eccentricity_truncation_lvl == 2
    assert world_ctl.max_tidal_order_lvl == 2
    assert world_ctl.spin_period == 10.
    assert type(world_ctl.tidal_susceptibility) in [float, np.float64]
    assert type(world_ctl.global_love_by_orderl[2]) in [complex, np.complex128]
    assert type(world_ctl.global_negative_imk_by_orderl[2]) in [float, np.float64]
    assert type(world_ctl.effective_q_by_orderl[2]) in [float, np.float64]
    assert type(world_ctl.dUdM) in [float, np.float64]
    assert type(world_ctl.dUdO) in [float, np.float64]
    assert type(world_ctl.dUdw) in [float, np.float64]
    assert type(world_ctl.tidal_heating_global) in [float, np.float64]
    # Check the orbital derivatives
    assert type(orbit.get_eccentricity_time_derivative(world_ctl)) in [float, np.float64]
    assert type(orbit.get_semi_major_axis_time_derivative(world_ctl)) in [float, np.float64]
    assert type(orbit.get_orbital_motion_time_derivative(world_ctl)) in [float, np.float64]

    # Test arrays
    # Set orbital frequency and eccentricity - check that spin locking worked.
    orb_period = np.linspace(10., 50., 10)
    world_ctl.set_state(
        orbital_period=orb_period, eccentricity=0.2, obliquity=np.radians(10.),
        spin_period=0.5 * orb_period
        )
    # For a world in synchronous rotation the spin period should equal the new orbital period.
    assert not world_ctl.is_spin_sync
    assert world_ctl.tides_on
    assert world_ctl.eccentricity_truncation_lvl == 2
    assert world_ctl.max_tidal_order_lvl == 2
    np.testing.assert_allclose(world_ctl.spin_period, 0.5 * orb_period)
    assert type(world_ctl.tidal_susceptibility) == np.ndarray
    assert type(world_ctl.global_love_by_orderl[2]) == np.ndarray
    assert type(world_ctl.global_negative_imk_by_orderl[2]) == np.ndarray
    assert type(world_ctl.effective_q_by_orderl[2]) == np.ndarray
    assert type(world_ctl.dUdM) == np.ndarray
    assert type(world_ctl.dUdO) == np.ndarray
    assert type(world_ctl.dUdw) == np.ndarray
    assert type(world_ctl.tidal_heating_global) == np.ndarray
    # Check the orbital derivatives
    assert type(orbit.get_eccentricity_time_derivative(world_ctl)) == np.ndarray
    assert type(orbit.get_semi_major_axis_time_derivative(world_ctl)) == np.ndarray
    assert type(orbit.get_orbital_motion_time_derivative(world_ctl)) == np.ndarray
