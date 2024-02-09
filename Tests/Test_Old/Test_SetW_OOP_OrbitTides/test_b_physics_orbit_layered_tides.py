import numpy as np

import TidalPy


from TidalPy.structures import build_world, build_from_world
from TidalPy.structures.orbit import PhysicsOrbit

star = build_world('55cnc')
world = build_world('io_simple')


def test_layered_tide_model_maxwell_rheology():
    """ This test will load a maxwell model into a layered Io and test the tidal calculations """

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

    world_tidal = build_from_world(world, new_config=config_)
    orbit = PhysicsOrbit(star, tidal_host=star, tidal_bodies=world_tidal)

    # Test loading in layer temperature and performing necessary steps to eventually calculate the complex compliance
    #    once a frequency is set.
    world_tidal.mantle.temperature = 1400.
    assert type(world_tidal.mantle.viscosity) in [float, np.float64]
    assert type(world_tidal.mantle.shear_modulus) in [float, np.float64]

    # Test Floats
    # Set orbital frequency and eccentricity - check that spin locking worked.
    world_tidal.set_state(
        orbital_period=50., eccentricity=0.2, obliquity=np.radians(10.),
        spin_period=10.
        )
    # For a world in synchronous rotation the spin period should equal the new orbital period.
    assert not world_tidal.is_spin_sync
    assert world_tidal.tides_on
    assert world_tidal.eccentricity_truncation_lvl == 2
    assert world_tidal.max_tidal_order_lvl == 2
    assert world_tidal.spin_period == 10.
    assert type(world_tidal.tidal_susceptibility) in [float, np.float64]
    assert type(world_tidal.global_love_by_orderl[2]) in [complex, np.complex128]
    assert type(world_tidal.global_negative_imk_by_orderl[2]) in [float, np.float64]
    assert type(world_tidal.effective_q_by_orderl[2]) in [float, np.float64]
    assert type(world_tidal.dUdM) in [float, np.float64]
    assert type(world_tidal.dUdO) in [float, np.float64]
    assert type(world_tidal.dUdw) in [float, np.float64]
    assert type(world_tidal.tidal_heating_global) in [float, np.float64]

    # Test arrays
    # Set orbital frequency and eccentricity - check that spin locking worked.
    orb_period = np.linspace(10., 50., 10)
    world_tidal.set_state(
        orbital_period=orb_period, eccentricity=0.2, obliquity=np.radians(10.),
        spin_period=0.5 * orb_period
        )
    # For a world in synchronous rotation the spin period should equal the new orbital period.
    assert not world_tidal.is_spin_sync
    assert world_tidal.tides_on
    assert world_tidal.eccentricity_truncation_lvl == 2
    assert world_tidal.max_tidal_order_lvl == 2
    np.testing.assert_allclose(world_tidal.spin_period, 0.5 * orb_period)
    assert type(world_tidal.tidal_susceptibility) == np.ndarray
    assert type(world_tidal.global_love_by_orderl[2]) == np.ndarray
    assert type(world_tidal.global_negative_imk_by_orderl[2]) == np.ndarray
    assert type(world_tidal.effective_q_by_orderl[2]) == np.ndarray
    assert type(world_tidal.dUdM) == np.ndarray
    assert type(world_tidal.dUdO) == np.ndarray
    assert type(world_tidal.dUdw) == np.ndarray
    assert type(world_tidal.tidal_heating_global) == np.ndarray
