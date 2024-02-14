import numpy as np
from numba.typed.typeddict import Dict

import TidalPy


from TidalPy.toolbox.quick_tides import (dual_dissipation_from_dict_or_world_instance,
                                         quick_dual_body_tidal_dissipation)


def test_quick_tidal_dissipation_fixedq_float():
    """ This will test the quick tidal dissipation calculator - using all floats for fixed-q rheology """

    orbital_period = 45.
    host_mass = 1.e27
    host_radius = 1.e6
    host_gravity = 10.
    host_density = 5000.
    host_spin_period = 1.2
    host_moi = 1.e5
    host_obliquity = 0.1

    target_radius = 1.e6
    target_mass = 1.e24
    target_gravity = 10.
    target_density = 5000.
    target_moi = 1.e5
    target_spin_period = 1.2
    target_obliquity = 0.1
    eccentricity = 0.1
    rheology = 'cpl'
    dissipation_results = \
        quick_dual_body_tidal_dissipation(
            radii=(host_radius, target_radius),
            masses=(host_mass, target_mass),
            gravities=(host_gravity, target_gravity),
            densities=(host_density, target_density),
            mois=(host_moi, target_moi),
            spin_periods=(host_spin_period, target_spin_period),
            obliquities=(host_obliquity, target_obliquity),
            tidal_scales=(1., 1.),
            rheologies=rheology,
            eccentricity=eccentricity, orbital_period=orbital_period,
            max_tidal_order_l=2, eccentricity_truncation_lvl=2,
            use_obliquity=True, fixed_qs=(10000., 120.)
            )

    assert type(dissipation_results) == dict
    for obj_name in ['host', 'secondary']:
        world_results = dissipation_results[obj_name]
        assert type(world_results['tidal_heating']) in [float, np.float64]
        assert type(world_results['dUdM']) in [float, np.float64]
        assert type(world_results['dUdw']) in [float, np.float64]
        assert type(world_results['dUdO']) in [float, np.float64]
        assert type(world_results['love_number_by_orderl']) in [Dict, dict]
        assert type(world_results['negative_imk_by_orderl']) in [Dict, dict]
        assert type(world_results['effective_q_by_orderl']) in [Dict, dict]
        assert len(world_results['effective_q_by_orderl']) == 1
        assert type(world_results['love_number_by_orderl'][2]) in [complex, np.complex128]
        assert type(world_results['negative_imk_by_orderl'][2]) in [float, np.float64]
        assert type(world_results['effective_q_by_orderl'][2]) in [float, np.float64]
    assert type(dissipation_results['eccentricity_derivative']) in [float, np.float64]
    assert type(dissipation_results['semi_major_axis_derivative']) in [float, np.float64]


def test_quick_tidal_dissipation_fixedq_orbital_array():
    """ This will test the quick tidal dissipation calculator - using orbital array with the fixed-q rheology """

    orbital_period = np.linspace(10., 40., 10)
    host_mass = 1.e27
    host_radius = 1.e6
    host_gravity = 10.
    host_density = 5000.
    host_spin_period = 1.2
    host_moi = 1.e5
    host_obliquity = 0.1

    target_radius = 1.e6
    target_mass = 1.e24
    target_gravity = 10.
    target_density = 5000.
    target_moi = 1.e5
    target_spin_period = 1.2
    target_obliquity = 0.1
    eccentricity = 0.1
    rheology = 'cpl'
    dissipation_results = \
        quick_dual_body_tidal_dissipation(
            radii=(host_radius, target_radius),
            masses=(host_mass, target_mass),
            gravities=(host_gravity, target_gravity),
            densities=(host_density, target_density),
            mois=(host_moi, target_moi),
            spin_periods=(host_spin_period, target_spin_period),
            obliquities=(host_obliquity, target_obliquity),
            tidal_scales=(1., 1.),
            rheologies=rheology,
            eccentricity=eccentricity, orbital_period=orbital_period,
            max_tidal_order_l=2, eccentricity_truncation_lvl=2,
            use_obliquity=True, fixed_qs=(10000., 120.)
            )

    assert type(dissipation_results) == dict
    for obj_name in ['host', 'secondary']:
        world_results = dissipation_results[obj_name]

        assert type(world_results['tidal_heating']) == np.ndarray
        assert type(world_results['dUdM']) == np.ndarray
        assert type(world_results['dUdw']) == np.ndarray
        assert type(world_results['dUdO']) == np.ndarray
        assert type(world_results['love_number_by_orderl']) in [Dict, dict]
        assert type(world_results['negative_imk_by_orderl']) in [Dict, dict]
        assert type(world_results['effective_q_by_orderl']) in [Dict, dict]
        assert len(world_results['effective_q_by_orderl']) == 1
        assert type(world_results['love_number_by_orderl'][2]) == np.ndarray
        assert type(world_results['negative_imk_by_orderl'][2]) == np.ndarray
        assert type(world_results['effective_q_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['eccentricity_derivative']) == np.ndarray
    assert type(dissipation_results['semi_major_axis_derivative']) == np.ndarray


def test_quick_tidal_dissipation_fixedq_orbital_array_NSR():
    """ This will test the quick tidal dissipation calculator - using orbital array with the fixed-q rheology - NSR """

    orbital_period = np.linspace(10., 40., 10)
    host_mass = 1.e27
    host_radius = 1.e6
    host_gravity = 10.
    host_density = 5000.
    host_spin_period = np.linspace(5., 10., 10)
    host_moi = 1.e5
    host_obliquity = 0.1

    target_radius = 1.e6
    target_mass = 1.e24
    target_gravity = 10.
    target_density = 5000.
    target_moi = 1.e5
    target_spin_period = np.linspace(5., 10., 10)
    target_obliquity = 0.1
    eccentricity = 0.1
    rheology = 'cpl'
    dissipation_results = \
        quick_dual_body_tidal_dissipation(
            radii=(host_radius, target_radius),
            masses=(host_mass, target_mass),
            gravities=(host_gravity, target_gravity),
            densities=(host_density, target_density),
            mois=(host_moi, target_moi),
            spin_periods=(host_spin_period, target_spin_period),
            obliquities=(host_obliquity, target_obliquity),
            tidal_scales=(1., 1.),
            rheologies=rheology,
            eccentricity=eccentricity, orbital_period=orbital_period,
            max_tidal_order_l=2, eccentricity_truncation_lvl=2,
            use_obliquity=True, fixed_qs=(10000., 120.)
            )

    assert type(dissipation_results) == dict
    for obj_name in ['host', 'secondary']:
        world_results = dissipation_results[obj_name]

        assert type(world_results['tidal_heating']) == np.ndarray
        assert type(world_results['dUdM']) == np.ndarray
        assert type(world_results['dUdw']) == np.ndarray
        assert type(world_results['dUdO']) == np.ndarray
        assert type(world_results['love_number_by_orderl']) in [Dict, dict]
        assert type(world_results['negative_imk_by_orderl']) in [Dict, dict]
        assert type(world_results['effective_q_by_orderl']) in [Dict, dict]
        assert len(world_results['effective_q_by_orderl']) == 1
        assert type(world_results['love_number_by_orderl'][2]) == np.ndarray
        assert type(world_results['negative_imk_by_orderl'][2]) == np.ndarray
        assert type(world_results['effective_q_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['eccentricity_derivative']) == np.ndarray
    assert type(dissipation_results['semi_major_axis_derivative']) == np.ndarray


def test_quick_tidal_dissipation_maxwell_orbital_array_NSR():
    """ This will test the quick tidal dissipation calculator - using orbital array with the Maxwell rheology - NSR """

    viscosity = 1.e22
    shear = 1.e10

    orbital_period = np.linspace(10., 40., 10)
    host_mass = 1.e27
    host_radius = 1.e6
    host_gravity = 10.
    host_density = 5000.
    host_spin_period = np.linspace(5., 10., 10)
    host_moi = 1.e5
    host_obliquity = 0.1

    target_radius = 1.e6
    target_mass = 1.e24
    target_gravity = 10.
    target_density = 5000.
    target_moi = 1.e5
    target_spin_period = np.linspace(5., 10., 10)
    target_obliquity = 0.1
    eccentricity = 0.1
    rheology = 'maxwell'
    dissipation_results = \
        quick_dual_body_tidal_dissipation(
            radii=(host_radius, target_radius),
            masses=(host_mass, target_mass),
            gravities=(host_gravity, target_gravity),
            densities=(host_density, target_density),
            mois=(host_moi, target_moi),
            spin_periods=(host_spin_period, target_spin_period),
            obliquities=(host_obliquity, target_obliquity),
            tidal_scales=(1., 1.),
            viscosities=(viscosity, viscosity), shear_moduli=(shear, shear),
            rheologies=rheology,
            eccentricity=eccentricity, orbital_period=orbital_period,
            max_tidal_order_l=2, eccentricity_truncation_lvl=2,
            use_obliquity=True, fixed_qs=(10000., 120.)
            )

    assert type(dissipation_results) == dict
    for obj_name in ['host', 'secondary']:
        world_results = dissipation_results[obj_name]

        assert type(world_results['tidal_heating']) == np.ndarray
        assert type(world_results['dUdM']) == np.ndarray
        assert type(world_results['dUdw']) == np.ndarray
        assert type(world_results['dUdO']) == np.ndarray
        assert type(world_results['love_number_by_orderl']) in [Dict, dict]
        assert type(world_results['negative_imk_by_orderl']) in [Dict, dict]
        assert type(world_results['effective_q_by_orderl']) in [Dict, dict]
        assert len(world_results['effective_q_by_orderl']) == 1
        assert type(world_results['love_number_by_orderl'][2]) == np.ndarray
        assert type(world_results['negative_imk_by_orderl'][2]) == np.ndarray
        assert type(world_results['effective_q_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['eccentricity_derivative']) == np.ndarray
    assert type(dissipation_results['semi_major_axis_derivative']) == np.ndarray


def test_quick_tidal_dissipation_andrade_orbital_array_NSR():
    """ This will test the quick tidal dissipation calculator - using orbital array with the Andrade rheology - NSR """

    viscosity = 1.e22
    shear = 1.e10

    orbital_period = np.linspace(10., 40., 10)
    host_mass = 1.e27
    host_radius = 1.e6
    host_gravity = 10.
    host_density = 5000.
    host_spin_period = np.linspace(5., 10., 10)
    host_moi = 1.e5
    host_obliquity = 0.1

    target_radius = 1.e6
    target_mass = 1.e24
    target_gravity = 10.
    target_density = 5000.
    target_moi = 1.e5
    target_spin_period = np.linspace(5., 10., 10)
    target_obliquity = 0.1
    eccentricity = 0.1
    rheology = 'andrade'
    andrade_inputs = (0.2, 1.)
    dissipation_results = \
        quick_dual_body_tidal_dissipation(
            radii=(host_radius, target_radius),
            masses=(host_mass, target_mass),
            gravities=(host_gravity, target_gravity),
            densities=(host_density, target_density),
            mois=(host_moi, target_moi),
            spin_periods=(host_spin_period, target_spin_period),
            obliquities=(host_obliquity, target_obliquity),
            tidal_scales=(1., 1.),
            viscosities=(viscosity, viscosity), shear_moduli=(shear, shear),
            rheologies=rheology, complex_compliance_inputs=andrade_inputs,
            eccentricity=eccentricity, orbital_period=orbital_period,
            max_tidal_order_l=2, eccentricity_truncation_lvl=2,
            use_obliquity=True, fixed_qs=(10000., 120.)
            )

    assert type(dissipation_results) == dict
    for obj_name in ['host', 'secondary']:
        world_results = dissipation_results[obj_name]

        assert type(world_results['tidal_heating']) == np.ndarray
        assert type(world_results['dUdM']) == np.ndarray
        assert type(world_results['dUdw']) == np.ndarray
        assert type(world_results['dUdO']) == np.ndarray
        assert type(world_results['love_number_by_orderl']) in [Dict, dict]
        assert type(world_results['negative_imk_by_orderl']) in [Dict, dict]
        assert type(world_results['effective_q_by_orderl']) in [Dict, dict]
        assert len(world_results['effective_q_by_orderl']) == 1
        assert type(world_results['love_number_by_orderl'][2]) == np.ndarray
        assert type(world_results['negative_imk_by_orderl'][2]) == np.ndarray
        assert type(world_results['effective_q_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['eccentricity_derivative']) == np.ndarray
    assert type(dissipation_results['semi_major_axis_derivative']) == np.ndarray


def test_quick_tidal_dissipation_andrade_visco_array_NSR():
    """ This will test the quick tidal dissipation calculator - using viscosity array with the Andrade rheology - NSR """

    viscosity = np.logspace(10., 40., 10)
    shear = 1.e10

    orbital_period = 40.
    host_mass = 1.e27
    host_radius = 1.e6
    host_gravity = 10.
    host_density = 5000.
    host_spin_period = 10.
    host_moi = 1.e5
    host_obliquity = 0.1

    target_radius = 1.e6
    target_mass = 1.e24
    target_gravity = 10.
    target_density = 5000.
    target_moi = 1.e5
    target_spin_period = 10.
    target_obliquity = 0.1
    eccentricity = 0.1
    rheology = 'andrade'
    andrade_inputs = (0.2, 1.)
    dissipation_results = \
        quick_dual_body_tidal_dissipation(
            radii=(host_radius, target_radius),
            masses=(host_mass, target_mass),
            gravities=(host_gravity, target_gravity),
            densities=(host_density, target_density),
            mois=(host_moi, target_moi),
            spin_periods=(host_spin_period, target_spin_period),
            obliquities=(host_obliquity, target_obliquity),
            tidal_scales=(1., 1.),
            viscosities=(viscosity, viscosity), shear_moduli=(shear, shear),
            rheologies=rheology, complex_compliance_inputs=andrade_inputs,
            eccentricity=eccentricity, orbital_period=orbital_period,
            max_tidal_order_l=2, eccentricity_truncation_lvl=2,
            use_obliquity=True, fixed_qs=(10000., 120.)
            )

    assert type(dissipation_results) == dict
    for obj_name in ['host', 'secondary']:
        world_results = dissipation_results[obj_name]

        assert type(world_results['tidal_heating']) == np.ndarray
        assert type(world_results['dUdM']) == np.ndarray
        assert type(world_results['dUdw']) == np.ndarray
        assert type(world_results['dUdO']) == np.ndarray
        assert type(world_results['love_number_by_orderl']) in [Dict, dict]
        assert type(world_results['negative_imk_by_orderl']) in [Dict, dict]
        assert type(world_results['effective_q_by_orderl']) in [Dict, dict]
        assert len(world_results['effective_q_by_orderl']) == 1
        assert type(world_results['love_number_by_orderl'][2]) == np.ndarray
        assert type(world_results['negative_imk_by_orderl'][2]) == np.ndarray
        assert type(world_results['effective_q_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['eccentricity_derivative']) == np.ndarray
    assert type(dissipation_results['semi_major_axis_derivative']) == np.ndarray


def test_quick_tidal_dissipation_andrade_visco_array_higher_l():
    """ This will test the quick tidal dissipation calculator - using viscosity array with the Andrade rheology
        Max l = 3, eccentricity trunc = 4 """

    viscosity = np.logspace(10., 40., 10)
    shear = 1.e10

    orbital_period = 40.
    host_mass = 1.e27
    host_radius = 1.e6
    host_gravity = 10.
    host_density = 5000.
    host_spin_period = 10.
    host_moi = 1.e5
    host_obliquity = 0.1

    target_radius = 1.e6
    target_mass = 1.e24
    target_gravity = 10.
    target_density = 5000.
    target_moi = 1.e5
    target_spin_period = 10.
    target_obliquity = 0.1
    eccentricity = 0.1
    rheology = 'andrade'
    andrade_inputs = (0.2, 1.)
    dissipation_results = \
        quick_dual_body_tidal_dissipation(
            radii=(host_radius, target_radius),
            masses=(host_mass, target_mass),
            gravities=(host_gravity, target_gravity),
            densities=(host_density, target_density),
            mois=(host_moi, target_moi),
            spin_periods=(host_spin_period, target_spin_period),
            obliquities=(host_obliquity, target_obliquity),
            tidal_scales=(1., 1.),
            viscosities=(viscosity, viscosity), shear_moduli=(shear, shear),
            rheologies=rheology, complex_compliance_inputs=andrade_inputs,
            eccentricity=eccentricity, orbital_period=orbital_period,
            max_tidal_order_l=3, eccentricity_truncation_lvl=6,
            use_obliquity=True, fixed_qs=(10000., 120.)
            )

    assert type(dissipation_results) == dict
    for obj_name in ['host', 'secondary']:
        world_results = dissipation_results[obj_name]

        assert type(world_results['tidal_heating']) == np.ndarray
        assert type(world_results['dUdM']) == np.ndarray
        assert type(world_results['dUdw']) == np.ndarray
        assert type(world_results['dUdO']) == np.ndarray
        assert type(world_results['love_number_by_orderl']) in [Dict, dict]
        assert type(world_results['negative_imk_by_orderl']) in [Dict, dict]
        assert type(world_results['effective_q_by_orderl']) in [Dict, dict]
        assert len(world_results['effective_q_by_orderl']) == 2
        assert type(world_results['love_number_by_orderl'][2]) == np.ndarray
        assert type(world_results['negative_imk_by_orderl'][2]) == np.ndarray
        assert type(world_results['effective_q_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['eccentricity_derivative']) == np.ndarray
    assert type(dissipation_results['semi_major_axis_derivative']) == np.ndarray


def test_quick_tidal_dissipation_andrade_visco_array_higher_l_from_dict():
    """ This will test the quick tidal dissipation (parameters pulled from a dict) calculator -
        using viscosity array with the Andrade rheology. Max l = 3, eccentricity trunc = 4 """

    secondary_dict = {
        'radius'         : 1.e6,
        'mass'           : 1.e24,
        'gravity_surface': 10.,
        'density_bulk'   : 5000.,
        'moi'            : 1.e10
        }

    host_dict = {
        'radius'         : 1.e7,
        'mass'           : 1.e27,
        'gravity_surface': 20.,
        'density_bulk'   : 5000.,
        'moi'            : 1.e10
        }

    viscosity = np.logspace(10., 40., 10)
    shear = 1.e10

    orbital_period = 40.
    host_spin_period = 10.
    host_obliquity = 0.1

    target_spin_period = 10.
    target_obliquity = 0.1
    eccentricity = 0.1
    rheology = 'andrade'
    andrade_inputs = (0.2, 1.)
    dissipation_results = \
        dual_dissipation_from_dict_or_world_instance(
            host_dict, secondary_dict,
            spin_periods=(host_spin_period, target_spin_period),
            obliquities=(host_obliquity, target_obliquity),
            tidal_scales=(1., 1.),
            viscosities=(viscosity, viscosity), shear_moduli=(shear, shear),
            rheologies=rheology, complex_compliance_inputs=andrade_inputs,
            eccentricity=eccentricity, orbital_period=orbital_period,
            max_tidal_order_l=3, eccentricity_truncation_lvl=6,
            use_obliquity=True, fixed_qs=(10000., 120.)
            )

    assert type(dissipation_results) == dict
    for obj_name in ['host', 'secondary']:
        world_results = dissipation_results[obj_name]

        assert type(world_results['tidal_heating']) == np.ndarray
        assert type(world_results['dUdM']) == np.ndarray
        assert type(world_results['dUdw']) == np.ndarray
        assert type(world_results['dUdO']) == np.ndarray
        assert type(world_results['love_number_by_orderl']) in [Dict, dict]
        assert type(world_results['negative_imk_by_orderl']) in [Dict, dict]
        assert type(world_results['effective_q_by_orderl']) in [Dict, dict]
        assert len(world_results['effective_q_by_orderl']) == 2
        assert type(world_results['love_number_by_orderl'][2]) == np.ndarray
        assert type(world_results['negative_imk_by_orderl'][2]) == np.ndarray
        assert type(world_results['effective_q_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['eccentricity_derivative']) == np.ndarray
    assert type(dissipation_results['semi_major_axis_derivative']) == np.ndarray


def test_quick_tidal_dissipation_andrade_visco_array_higher_l_from_world_instance():
    """ This will test the quick tidal dissipation (parameters pulled from a world instance) calculator -
        using viscosity array with the Andrade rheology. Max l = 3, eccentricity trunc = 4 """

    from TidalPy.structures import build_world

    jupiter = build_world('jupiter')
    io = build_world('io_simple')

    viscosity = np.logspace(10., 40., 10)
    shear = 1.e10

    orbital_period = 40.
    host_spin_period = 10.
    host_obliquity = 0.1

    target_spin_period = 10.
    target_obliquity = 0.1
    eccentricity = 0.1
    rheology = 'andrade'
    andrade_inputs = (0.2, 1.)
    dissipation_results = \
        dual_dissipation_from_dict_or_world_instance(
            jupiter, io,
            spin_periods=(host_spin_period, target_spin_period),
            obliquities=(host_obliquity, target_obliquity),
            tidal_scales=(1., 1.),
            viscosities=(viscosity, viscosity), shear_moduli=(shear, shear),
            rheologies=rheology, complex_compliance_inputs=andrade_inputs,
            eccentricity=eccentricity, orbital_period=orbital_period,
            max_tidal_order_l=3, eccentricity_truncation_lvl=6,
            use_obliquity=True, fixed_qs=(10000., 120.)
            )

    assert type(dissipation_results) == dict
    for obj_name in ['host', 'secondary']:
        world_results = dissipation_results[obj_name]

        assert type(world_results['tidal_heating']) == np.ndarray
        assert type(world_results['dUdM']) == np.ndarray
        assert type(world_results['dUdw']) == np.ndarray
        assert type(world_results['dUdO']) == np.ndarray
        assert type(world_results['love_number_by_orderl']) in [Dict, dict]
        assert type(world_results['negative_imk_by_orderl']) in [Dict, dict]
        assert type(world_results['effective_q_by_orderl']) in [Dict, dict]
        assert len(world_results['effective_q_by_orderl']) == 2
        assert type(world_results['love_number_by_orderl'][2]) == np.ndarray
        assert type(world_results['negative_imk_by_orderl'][2]) == np.ndarray
        assert type(world_results['effective_q_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['eccentricity_derivative']) == np.ndarray
    assert type(dissipation_results['semi_major_axis_derivative']) == np.ndarray
