import numpy as np
from numba.typed.typeddict import Dict

import TidalPy


from TidalPy.toolbox.quick_tides import quick_tidal_dissipation, single_dissipation_from_dict_or_world_instance


def test_quick_tidal_dissipation_fixedq_float():
    """ This will test the quick tidal dissipation calculator - using all floats for fixed-q rheology """

    orbital_period = 45.
    host_mass = 1.e27
    target_radius = 1.e6
    target_mass = 1.e24
    target_gravity = 10.
    target_density = 5000.
    target_moi = 1.e5
    eccentricity = 0.1
    obliquity = 0.1
    rheology = 'cpl'
    dissipation_results = \
        quick_tidal_dissipation(
            host_mass, target_radius, target_mass, target_gravity, target_density, target_moi,
            rheology=rheology, eccentricity=eccentricity, obliquity=obliquity,
            orbital_period=orbital_period, spin_period=None,
            max_tidal_order_l=2, eccentricity_truncation_lvl=2,
            use_obliquity=True, tidal_scale=1., fixed_q=120.
            )

    assert type(dissipation_results) == dict
    assert type(dissipation_results['tidal_heating']) in [float, np.float64]
    assert type(dissipation_results['dUdM']) in [float, np.float64]
    assert type(dissipation_results['dUdw']) in [float, np.float64]
    assert type(dissipation_results['dUdO']) in [float, np.float64]
    assert type(dissipation_results['love_number_by_orderl']) in [Dict, dict]
    assert type(dissipation_results['negative_imk_by_orderl']) in [Dict, dict]
    assert type(dissipation_results['effective_q_by_orderl']) in [Dict, dict]
    assert len(dissipation_results['effective_q_by_orderl']) == 1
    assert type(dissipation_results['love_number_by_orderl'][2]) in [complex, np.complex128]
    assert type(dissipation_results['negative_imk_by_orderl'][2]) in [float, np.float64]
    assert type(dissipation_results['effective_q_by_orderl'][2]) in [float, np.float64]


def test_quick_tidal_dissipation_fixedq_orbital_array():
    """ This will test the quick tidal dissipation calculator - using orbital array with the fixed-q rheology """

    orbital_period = np.linspace(10., 40., 10)
    host_mass = 1.e27
    target_radius = 1.e6
    target_mass = 1.e24
    target_gravity = 10.
    target_density = 5000.
    target_moi = 1.e5
    eccentricity = 0.1
    obliquity = 0.1
    rheology = 'cpl'
    dissipation_results = \
        quick_tidal_dissipation(
            host_mass, target_radius, target_mass, target_gravity, target_density, target_moi,
            rheology=rheology, eccentricity=eccentricity, obliquity=obliquity,
            orbital_period=orbital_period, spin_period=None,
            max_tidal_order_l=2, eccentricity_truncation_lvl=2,
            use_obliquity=True, tidal_scale=1., fixed_q=120.
            )

    assert type(dissipation_results) == dict
    assert type(dissipation_results['tidal_heating']) == np.ndarray
    assert type(dissipation_results['dUdM']) == np.ndarray
    assert type(dissipation_results['dUdw']) == np.ndarray
    assert type(dissipation_results['dUdO']) == np.ndarray
    assert type(dissipation_results['love_number_by_orderl']) in [Dict, dict]
    assert type(dissipation_results['negative_imk_by_orderl']) in [Dict, dict]
    assert type(dissipation_results['effective_q_by_orderl']) in [Dict, dict]
    assert len(dissipation_results['effective_q_by_orderl']) == 1
    assert type(dissipation_results['love_number_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['negative_imk_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['effective_q_by_orderl'][2]) == np.ndarray


def test_quick_tidal_dissipation_fixedq_orbital_array_NSR():
    """ This will test the quick tidal dissipation calculator - using orbital array with the fixed-q rheology - NSR """

    orbital_period = np.linspace(10., 40., 10)
    spin_period = np.linspace(5., 10., 10)
    host_mass = 1.e27
    target_radius = 1.e6
    target_mass = 1.e24
    target_gravity = 10.
    target_density = 5000.
    target_moi = 1.e5
    eccentricity = 0.1
    obliquity = 0.1
    rheology = 'cpl'
    dissipation_results = \
        quick_tidal_dissipation(
            host_mass, target_radius, target_mass, target_gravity, target_density, target_moi,
            rheology=rheology, eccentricity=eccentricity, obliquity=obliquity,
            orbital_period=orbital_period, spin_period=spin_period,
            max_tidal_order_l=2, eccentricity_truncation_lvl=2,
            use_obliquity=True, tidal_scale=1., fixed_q=120.
            )

    assert type(dissipation_results) == dict
    assert type(dissipation_results['tidal_heating']) == np.ndarray
    assert type(dissipation_results['dUdM']) == np.ndarray
    assert type(dissipation_results['dUdw']) == np.ndarray
    assert type(dissipation_results['dUdO']) == np.ndarray
    assert type(dissipation_results['love_number_by_orderl']) in [Dict, dict]
    assert type(dissipation_results['negative_imk_by_orderl']) in [Dict, dict]
    assert type(dissipation_results['effective_q_by_orderl']) in [Dict, dict]
    assert len(dissipation_results['effective_q_by_orderl']) == 1
    assert type(dissipation_results['love_number_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['negative_imk_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['effective_q_by_orderl'][2]) == np.ndarray


def test_quick_tidal_dissipation_maxwell_orbital_array_NSR():
    """ This will test the quick tidal dissipation calculator - using orbital array with the Maxwell rheology - NSR """

    orbital_period = np.linspace(10., 40., 10)
    spin_period = np.linspace(5., 10., 10)
    host_mass = 1.e27
    target_radius = 1.e6
    target_mass = 1.e24
    target_gravity = 10.
    target_density = 5000.
    target_moi = 1.e5
    eccentricity = 0.1
    obliquity = 0.1
    viscosity = 1.e22
    shear = 1.e10
    rheology = 'maxwell'
    dissipation_results = \
        quick_tidal_dissipation(
            host_mass, target_radius, target_mass, target_gravity, target_density, target_moi,
            rheology=rheology, eccentricity=eccentricity, obliquity=obliquity,
            orbital_period=orbital_period, spin_period=spin_period,
            viscosity=viscosity, shear_modulus=shear,
            max_tidal_order_l=2, eccentricity_truncation_lvl=2,
            use_obliquity=True, tidal_scale=1., fixed_q=120.
            )

    assert type(dissipation_results) == dict
    assert type(dissipation_results['tidal_heating']) == np.ndarray
    assert type(dissipation_results['dUdM']) == np.ndarray
    assert type(dissipation_results['dUdw']) == np.ndarray
    assert type(dissipation_results['dUdO']) == np.ndarray
    assert type(dissipation_results['love_number_by_orderl']) in [Dict, dict]
    assert type(dissipation_results['negative_imk_by_orderl']) in [Dict, dict]
    assert type(dissipation_results['effective_q_by_orderl']) in [Dict, dict]
    assert len(dissipation_results['effective_q_by_orderl']) == 1
    assert type(dissipation_results['love_number_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['negative_imk_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['effective_q_by_orderl'][2]) == np.ndarray


def test_quick_tidal_dissipation_andrade_orbital_array_NSR():
    """ This will test the quick tidal dissipation calculator - using orbital array with the Andrade rheology - NSR """

    orbital_period = np.linspace(10., 40., 10)
    spin_period = np.linspace(5., 10., 10)
    host_mass = 1.e27
    target_radius = 1.e6
    target_mass = 1.e24
    target_gravity = 10.
    target_density = 5000.
    target_moi = 1.e5
    eccentricity = 0.1
    obliquity = 0.1
    viscosity = 1.e22
    shear = 1.e10
    rheology = 'andrade'
    andrade_inputs = (0.2, 1.)
    dissipation_results = \
        quick_tidal_dissipation(
            host_mass, target_radius, target_mass, target_gravity, target_density, target_moi,
            rheology=rheology, eccentricity=eccentricity, obliquity=obliquity,
            orbital_period=orbital_period, spin_period=spin_period,
            viscosity=viscosity, shear_modulus=shear,
            max_tidal_order_l=2, eccentricity_truncation_lvl=2,
            complex_compliance_inputs=andrade_inputs,
            use_obliquity=True, tidal_scale=1., fixed_q=120.
            )

    assert type(dissipation_results) == dict
    assert type(dissipation_results['tidal_heating']) == np.ndarray
    assert type(dissipation_results['dUdM']) == np.ndarray
    assert type(dissipation_results['dUdw']) == np.ndarray
    assert type(dissipation_results['dUdO']) == np.ndarray
    assert type(dissipation_results['love_number_by_orderl']) in [Dict, dict]
    assert type(dissipation_results['negative_imk_by_orderl']) in [Dict, dict]
    assert type(dissipation_results['effective_q_by_orderl']) in [Dict, dict]
    assert len(dissipation_results['effective_q_by_orderl']) == 1
    assert type(dissipation_results['love_number_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['negative_imk_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['effective_q_by_orderl'][2]) == np.ndarray


def test_quick_tidal_dissipation_andrade_visco_array_NSR():
    """ This will test the quick tidal dissipation calculator - using viscosity array with the Andrade rheology - NSR """

    orbital_period = 40.
    spin_period = 10.
    host_mass = 1.e27
    target_radius = 1.e6
    target_mass = 1.e24
    target_gravity = 10.
    target_density = 5000.
    target_moi = 1.e5
    eccentricity = 0.1
    obliquity = 0.1
    viscosity = np.logspace(10., 40., 10)
    shear = 1.e10
    rheology = 'andrade'
    andrade_inputs = (0.2, 1.)
    dissipation_results = \
        quick_tidal_dissipation(
            host_mass, target_radius, target_mass, target_gravity, target_density, target_moi,
            rheology=rheology, eccentricity=eccentricity, obliquity=obliquity,
            orbital_period=orbital_period, spin_period=spin_period,
            viscosity=viscosity, shear_modulus=shear,
            max_tidal_order_l=2, eccentricity_truncation_lvl=2,
            complex_compliance_inputs=andrade_inputs,
            use_obliquity=True, tidal_scale=1., fixed_q=120.
            )

    assert type(dissipation_results) == dict
    assert type(dissipation_results['tidal_heating']) == np.ndarray
    assert type(dissipation_results['dUdM']) == np.ndarray
    assert type(dissipation_results['dUdw']) == np.ndarray
    assert type(dissipation_results['dUdO']) == np.ndarray
    assert type(dissipation_results['love_number_by_orderl']) in [Dict, dict]
    assert type(dissipation_results['negative_imk_by_orderl']) in [Dict, dict]
    assert type(dissipation_results['effective_q_by_orderl']) in [Dict, dict]
    assert len(dissipation_results['effective_q_by_orderl']) == 1
    assert type(dissipation_results['love_number_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['negative_imk_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['effective_q_by_orderl'][2]) == np.ndarray


def test_quick_tidal_dissipation_andrade_visco_array_higher_l():
    """ This will test the quick tidal dissipation calculator - using viscosity array with the Andrade rheology
        Max l = 3, eccentricity trunc = 4 """

    orbital_period = 40.
    spin_period = 10.
    host_mass = 1.e27
    target_radius = 1.e6
    target_mass = 1.e24
    target_gravity = 10.
    target_density = 5000.
    target_moi = 1.e5
    eccentricity = 0.1
    obliquity = 0.1
    viscosity = np.logspace(10., 40., 10)
    shear = 1.e10
    rheology = 'andrade'
    andrade_inputs = (0.2, 1.)
    dissipation_results = \
        quick_tidal_dissipation(
            host_mass, target_radius, target_mass, target_gravity, target_density, target_moi,
            rheology=rheology, eccentricity=eccentricity, obliquity=obliquity,
            orbital_period=orbital_period, spin_period=spin_period,
            viscosity=viscosity, shear_modulus=shear,
            max_tidal_order_l=3, eccentricity_truncation_lvl=4,
            complex_compliance_inputs=andrade_inputs,
            use_obliquity=True, tidal_scale=1., fixed_q=120.
            )

    assert type(dissipation_results) == dict
    assert type(dissipation_results['tidal_heating']) == np.ndarray
    assert type(dissipation_results['dUdM']) == np.ndarray
    assert type(dissipation_results['dUdw']) == np.ndarray
    assert type(dissipation_results['dUdO']) == np.ndarray
    assert type(dissipation_results['love_number_by_orderl']) in [Dict, dict]
    assert type(dissipation_results['negative_imk_by_orderl']) in [Dict, dict]
    assert type(dissipation_results['effective_q_by_orderl']) in [Dict, dict]
    assert len(dissipation_results['effective_q_by_orderl']) == 2
    assert type(dissipation_results['love_number_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['negative_imk_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['effective_q_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['love_number_by_orderl'][3]) == np.ndarray
    assert type(dissipation_results['negative_imk_by_orderl'][3]) == np.ndarray
    assert type(dissipation_results['effective_q_by_orderl'][3]) == np.ndarray


def test_quick_tidal_dissipation_andrade_visco_array_higher_l_from_dict():
    """ This will test the quick tidal dissipation (parameters pulled from a dict) calculator -
        using viscosity array with the Andrade rheology. Max l = 3, eccentricity trunc = 4 """

    host_dict = {
        'mass': 1.e27
        }

    secondary_dict = {
        'radius'         : 1.e6,
        'mass'           : 1.e24,
        'gravity_surface': 10.,
        'density_bulk'   : 5000.,
        'moi'            : 1.e10
        }

    orbital_period = 40.
    spin_period = 10.
    eccentricity = 0.1
    obliquity = 0.1
    viscosity = np.logspace(10., 40., 10)
    shear = 1.e10
    rheology = 'andrade'
    andrade_inputs = (0.2, 1.)
    kwargs = dict(
        rheology=rheology, eccentricity=eccentricity, obliquity=obliquity,
        orbital_period=orbital_period, spin_period=spin_period,
        viscosity=viscosity, shear_modulus=shear,
        max_tidal_order_l=3, eccentricity_truncation_lvl=4,
        complex_compliance_inputs=andrade_inputs,
        use_obliquity=True, tidal_scale=1., fixed_q=120.
    )
    dissipation_results = \
        single_dissipation_from_dict_or_world_instance(host_dict, secondary_dict, **kwargs)

    assert type(dissipation_results) == dict
    assert type(dissipation_results['tidal_heating']) == np.ndarray
    assert type(dissipation_results['dUdM']) == np.ndarray
    assert type(dissipation_results['dUdw']) == np.ndarray
    assert type(dissipation_results['dUdO']) == np.ndarray
    assert type(dissipation_results['love_number_by_orderl']) in [Dict, dict]
    assert type(dissipation_results['negative_imk_by_orderl']) in [Dict, dict]
    assert type(dissipation_results['effective_q_by_orderl']) in [Dict, dict]
    assert len(dissipation_results['effective_q_by_orderl']) == 2
    assert type(dissipation_results['love_number_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['negative_imk_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['effective_q_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['love_number_by_orderl'][3]) == np.ndarray
    assert type(dissipation_results['negative_imk_by_orderl'][3]) == np.ndarray
    assert type(dissipation_results['effective_q_by_orderl'][3]) == np.ndarray


def test_quick_tidal_dissipation_andrade_visco_array_higher_l_from_world_instance():
    """ This will test the quick tidal dissipation (parameters pulled from a world instance) calculator -
        using viscosity array with the Andrade rheology. Max l = 3, eccentricity trunc = 4 """

    from TidalPy.structures import build_world
    jupiter = build_world('jupiter')
    io = build_world('io_simple')

    orbital_period = 40.
    spin_period = 10.
    eccentricity = 0.1
    obliquity = 0.1
    viscosity = np.logspace(10., 40., 10)
    shear = 1.e10
    rheology = 'andrade'
    andrade_inputs = (0.2, 1.)
    kwargs = dict(
        rheology=rheology, eccentricity=eccentricity, obliquity=obliquity,
        orbital_period=orbital_period, spin_period=spin_period,
        viscosity=viscosity, shear_modulus=shear,
        max_tidal_order_l=3, eccentricity_truncation_lvl=4,
        complex_compliance_inputs=andrade_inputs,
        use_obliquity=True, tidal_scale=1., fixed_q=120.
    )

    dissipation_results = \
        single_dissipation_from_dict_or_world_instance(jupiter, io, **kwargs)

    assert type(dissipation_results) == dict
    assert type(dissipation_results['tidal_heating']) == np.ndarray
    assert type(dissipation_results['dUdM']) == np.ndarray
    assert type(dissipation_results['dUdw']) == np.ndarray
    assert type(dissipation_results['dUdO']) == np.ndarray
    assert type(dissipation_results['love_number_by_orderl']) in [Dict, dict]
    assert type(dissipation_results['negative_imk_by_orderl']) in [Dict, dict]
    assert type(dissipation_results['effective_q_by_orderl']) in [Dict, dict]
    assert len(dissipation_results['effective_q_by_orderl']) == 2
    assert type(dissipation_results['love_number_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['negative_imk_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['effective_q_by_orderl'][2]) == np.ndarray
    assert type(dissipation_results['love_number_by_orderl'][3]) == np.ndarray
    assert type(dissipation_results['negative_imk_by_orderl'][3]) == np.ndarray
    assert type(dissipation_results['effective_q_by_orderl'][3]) == np.ndarray
