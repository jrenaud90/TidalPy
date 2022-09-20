""" Tests for calculating the radial solution across the interior of a planet
"""

import TidalPy
import numpy as np
from TidalPy.constants import G
from TidalPy.tides.multilayer.numerical_int.solver import tidal_y_solver

TidalPy.config['stream_level'] = 'ERROR'
TidalPy.use_disk = False
TidalPy.reinit()

# Model planet - 2layers
density_array = 5000. * np.ones(10)
radius_array = np.linspace(0., 1.e6, 11)
R = radius_array[-1]
volume_array = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
mass_array = volume_array * density_array
planet_mass = sum(mass_array)
bulk_density = planet_mass / ((4. / 3.) * np.pi * R**3)
mass_below = np.asarray([np.sum(mass_array[:i + 1]) for i in range(10)])
gravity_array = G * mass_below / (radius_array[1:]**2)
shear_array = 5.e10 * np.ones(10, dtype=np.complex128)
bulk_array = 10.e10 * np.ones(10, dtype=np.float64)
radius_array_to_use = radius_array[1:]
frequency = 2. * np.pi / (86400. * 1.)


def test_calculate_homogen():
    """ Test the solution calculation for homogeneous planet """

    # Test static
    tidal_y = tidal_y_solver(
        model_name='homogeneous_solid',
        radius=radius_array_to_use, shear_modulus=shear_array, bulk_modulus=bulk_array,
        density=density_array, gravity=gravity_array, frequency=frequency,
        is_solid_by_layer=[True], is_static_by_layer=[True],
        indices_by_layer=[np.ones(radius_array_to_use.shape, dtype=np.bool)],
        order_l=2,
        surface_boundary_condition=None, solve_load_numbers=False,
        use_kamata=False,
        integrator='scipy', integration_method=None,
        integration_rtol=1.0e-6, integration_atol=1.0e-8,
        verbose=False, nondimensionalize=False, planet_bulk_density=bulk_density
        )

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex128, complex]

    # Test dynamic
    tidal_y = tidal_y_solver(
        model_name='homogeneous_solid',
        radius=radius_array_to_use, shear_modulus=shear_array, bulk_modulus=bulk_array,
        density=density_array, gravity=gravity_array, frequency=frequency,
        is_solid_by_layer=[True], is_static_by_layer=[False],
        indices_by_layer=[np.ones(radius_array_to_use.shape, dtype=np.bool)],
        order_l=2,
        surface_boundary_condition=None, solve_load_numbers=False,
        use_kamata=False,
        integrator='scipy', integration_method=None,
        integration_rtol=1.0e-6, integration_atol=1.0e-8,
        verbose=False, nondimensionalize=False, planet_bulk_density=bulk_density
        )

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex128, complex]


ls_layer_0_indx = radius_array_to_use <= R / 2.
ls_layer_1_indx = radius_array_to_use > R / 2.
shear_array_ls = np.copy(shear_array)
shear_array_ls[ls_layer_0_indx] = 0. + 0.j


def test_calculate_ls():
    """ Test the solution calculation for liquid-solid planet structure """

    # Test all static
    tidal_y = tidal_y_solver(
        model_name='liquid_solid',
        radius=radius_array_to_use, shear_modulus=shear_array_ls, bulk_modulus=bulk_array,
        density=density_array, gravity=gravity_array, frequency=frequency,
        is_solid_by_layer=[False, True], is_static_by_layer=[True, True],
        indices_by_layer=[ls_layer_0_indx, ls_layer_1_indx],
        order_l=2,
        surface_boundary_condition=None, solve_load_numbers=False,
        use_kamata=False,
        integrator='scipy', integration_method=None,
        integration_rtol=1.0e-6, integration_atol=1.0e-8,
        verbose=False, nondimensionalize=False, planet_bulk_density=bulk_density
        )

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex128, complex]

    # Test all dynamic
    tidal_y = tidal_y_solver(
        model_name='liquid_solid',
        radius=radius_array_to_use, shear_modulus=shear_array_ls, bulk_modulus=bulk_array,
        density=density_array, gravity=gravity_array, frequency=frequency,
        is_solid_by_layer=[False, True], is_static_by_layer=[False, False],
        indices_by_layer=[ls_layer_0_indx, ls_layer_1_indx],
        order_l=2,
        surface_boundary_condition=None, solve_load_numbers=False,
        use_kamata=False,
        integrator='scipy', integration_method=None,
        integration_rtol=1.0e-6, integration_atol=1.0e-8,
        verbose=False, nondimensionalize=False, planet_bulk_density=bulk_density
        )

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex128, complex]

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex128, complex]

    # Test mix liq=Static
    tidal_y = tidal_y_solver(
        model_name='liquid_solid',
        radius=radius_array_to_use, shear_modulus=shear_array_ls, bulk_modulus=bulk_array,
        density=density_array, gravity=gravity_array, frequency=frequency,
        is_solid_by_layer=[False, True], is_static_by_layer=[True, False],
        indices_by_layer=[ls_layer_0_indx, ls_layer_1_indx],
        order_l=2,
        surface_boundary_condition=None, solve_load_numbers=False,
        use_kamata=False,
        integrator='scipy', integration_method=None,
        integration_rtol=1.0e-6, integration_atol=1.0e-8,
        verbose=False, nondimensionalize=False, planet_bulk_density=bulk_density
        )

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex128, complex]

    # Test mix liq=Dynamic
    tidal_y = tidal_y_solver(
        model_name='liquid_solid',
        radius=radius_array_to_use, shear_modulus=shear_array_ls, bulk_modulus=bulk_array,
        density=density_array, gravity=gravity_array, frequency=frequency,
        is_solid_by_layer=[False, True], is_static_by_layer=[False, True],
        indices_by_layer=[ls_layer_0_indx, ls_layer_1_indx],
        order_l=2,
        surface_boundary_condition=None, solve_load_numbers=False,
        use_kamata=False,
        integrator='scipy', integration_method=None,
        integration_rtol=1.0e-6, integration_atol=1.0e-8,
        verbose=False, nondimensionalize=False, planet_bulk_density=bulk_density
        )

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex128, complex]


def test_calculate_ls_load_numbers():
    """ Test the solution calculation for liquid-solid planet structure """

    # Test all static
    tidal_y, load_y = tidal_y_solver(
        model_name='liquid_solid',
        radius=radius_array_to_use, shear_modulus=shear_array_ls, bulk_modulus=bulk_array,
        density=density_array, gravity=gravity_array, frequency=frequency,
        is_solid_by_layer=[False, True], is_static_by_layer=[True, True],
        indices_by_layer=[ls_layer_0_indx, ls_layer_1_indx],
        order_l=2,
        surface_boundary_condition=None, solve_load_numbers=True,
        use_kamata=False,
        integrator='scipy', integration_method=None,
        integration_rtol=1.0e-6, integration_atol=1.0e-8,
        verbose=False, nondimensionalize=False, planet_bulk_density=bulk_density
        )

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex128, complex]
    assert load_y.shape == (6, 10)
    assert load_y.dtype in [np.complex128, complex]

    # Test all dynamic
    tidal_y, load_y = tidal_y_solver(
        model_name='liquid_solid',
        radius=radius_array_to_use, shear_modulus=shear_array_ls, bulk_modulus=bulk_array,
        density=density_array, gravity=gravity_array, frequency=frequency,
        is_solid_by_layer=[False, True], is_static_by_layer=[False, False],
        indices_by_layer=[ls_layer_0_indx, ls_layer_1_indx],
        order_l=2,
        surface_boundary_condition=None, solve_load_numbers=True,
        use_kamata=False,
        integrator='scipy', integration_method=None,
        integration_rtol=1.0e-6, integration_atol=1.0e-8,
        verbose=False, nondimensionalize=False, planet_bulk_density=bulk_density
        )

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex128, complex]
    assert load_y.shape == (6, 10)
    assert load_y.dtype in [np.complex128, complex]

    # Test mix liq=Static
    tidal_y, load_y = tidal_y_solver(
        model_name='liquid_solid',
        radius=radius_array_to_use, shear_modulus=shear_array_ls, bulk_modulus=bulk_array,
        density=density_array, gravity=gravity_array, frequency=frequency,
        is_solid_by_layer=[False, True], is_static_by_layer=[True, False],
        indices_by_layer=[ls_layer_0_indx, ls_layer_1_indx],
        order_l=2,
        surface_boundary_condition=None, solve_load_numbers=True,
        use_kamata=False,
        integrator='scipy', integration_method=None,
        integration_rtol=1.0e-6, integration_atol=1.0e-8,
        verbose=False, nondimensionalize=False, planet_bulk_density=bulk_density
        )

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex128, complex]
    assert load_y.shape == (6, 10)
    assert load_y.dtype in [np.complex128, complex]

    # Test mix liq=Dynamic
    tidal_y, load_y = tidal_y_solver(
        model_name='liquid_solid',
        radius=radius_array_to_use, shear_modulus=shear_array_ls, bulk_modulus=bulk_array,
        density=density_array, gravity=gravity_array, frequency=frequency,
        is_solid_by_layer=[False, True], is_static_by_layer=[False, True],
        indices_by_layer=[ls_layer_0_indx, ls_layer_1_indx],
        order_l=2,
        surface_boundary_condition=None, solve_load_numbers=True,
        use_kamata=False,
        integrator='scipy', integration_method=None,
        integration_rtol=1.0e-6, integration_atol=1.0e-8,
        verbose=False, nondimensionalize=False, planet_bulk_density=bulk_density
        )

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex128, complex]
    assert load_y.shape == (6, 10)
    assert load_y.dtype in [np.complex128, complex]

sls_layer_0_indx = radius_array_to_use <= R / 3.
sls_layer_1_indx = np.logical_and(radius_array_to_use > R / 3., radius_array_to_use <= 2. * R / 3.)
sls_layer_2_indx = radius_array_to_use > 2. * R / 3.
shear_array_sls = np.copy(shear_array)
shear_array_sls[sls_layer_1_indx] = 0. + 0.j


def test_calculate_sls():
    """ Test the solution calculation for solid-liquid-solid planet structure """

    # Test all static
    tidal_y = tidal_y_solver(
        model_name='solid_liquid_solid',
        radius=radius_array_to_use, shear_modulus=shear_array_sls, bulk_modulus=bulk_array,
        density=density_array, gravity=gravity_array, frequency=frequency,
        is_solid_by_layer=[True, False, True], is_static_by_layer=[True, True, True],
        indices_by_layer=[sls_layer_0_indx, sls_layer_1_indx, sls_layer_2_indx],
        order_l=2,
        surface_boundary_condition=None, solve_load_numbers=False,
        use_kamata=False,
        integrator='scipy', integration_method=None,
        integration_rtol=1.0e-6, integration_atol=1.0e-8,
        verbose=False, nondimensionalize=False, planet_bulk_density=bulk_density
        )

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex128, complex]

    # Test all dynamic
    tidal_y = tidal_y_solver(
        model_name='solid_liquid_solid',
        radius=radius_array_to_use, shear_modulus=shear_array_sls, bulk_modulus=bulk_array,
        density=density_array, gravity=gravity_array, frequency=frequency,
        is_solid_by_layer=[True, False, True], is_static_by_layer=[False, False, False],
        indices_by_layer=[sls_layer_0_indx, sls_layer_1_indx, sls_layer_2_indx],
        order_l=2,
        surface_boundary_condition=None, solve_load_numbers=False,
        use_kamata=False,
        integrator='scipy', integration_method=None,
        integration_rtol=1.0e-6, integration_atol=1.0e-8,
        verbose=False, nondimensionalize=False, planet_bulk_density=bulk_density
        )

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex128, complex]

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex128, complex]

    # Test mix liq=Static
    tidal_y = tidal_y_solver(
        model_name='solid_liquid_solid',
        radius=radius_array_to_use, shear_modulus=shear_array_sls, bulk_modulus=bulk_array,
        density=density_array, gravity=gravity_array, frequency=frequency,
        is_solid_by_layer=[True, False, True], is_static_by_layer=[False, True, False],
        indices_by_layer=[sls_layer_0_indx, sls_layer_1_indx, sls_layer_2_indx],
        order_l=2,
        surface_boundary_condition=None, solve_load_numbers=False,
        use_kamata=False,
        integrator='scipy', integration_method=None,
        integration_rtol=1.0e-6, integration_atol=1.0e-8,
        verbose=False, nondimensionalize=False, planet_bulk_density=bulk_density
        )

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex128, complex]

    # Test mix liq=Dynamic
    tidal_y = tidal_y_solver(
        model_name='solid_liquid_solid',
        radius=radius_array_to_use, shear_modulus=shear_array_sls, bulk_modulus=bulk_array,
        density=density_array, gravity=gravity_array, frequency=frequency,
        is_solid_by_layer=[True, False, True], is_static_by_layer=[True, False, True],
        indices_by_layer=[sls_layer_0_indx, sls_layer_1_indx, sls_layer_2_indx],
        order_l=2,
        surface_boundary_condition=None, solve_load_numbers=False,
        use_kamata=False,
        integrator='scipy', integration_method=None,
        integration_rtol=1.0e-6, integration_atol=1.0e-8,
        verbose=False, nondimensionalize=False, planet_bulk_density=bulk_density
        )

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex128, complex]


ssls_layer_0_indx = radius_array_to_use <= R / 4.
ssls_layer_1_indx = np.logical_and(radius_array_to_use > R / 4., radius_array_to_use <= 2. * R / 4.)
ssls_layer_2_indx = np.logical_and(radius_array_to_use > 2. * R / 4., radius_array_to_use <= 3. * R / 4.)
ssls_layer_3_indx = radius_array_to_use > 3. * R / 4.
shear_array_ssls = np.copy(shear_array)
shear_array_ssls[ssls_layer_2_indx] = 0. + 0.j


def test_calculate_ssls():
    """ Test the solution calculation for solid-solid-liquid-solid planet structure """

    # Test all static
    tidal_y = tidal_y_solver(
        model_name='solid_solid_liquid_solid',
        radius=radius_array_to_use, shear_modulus=shear_array_ssls, bulk_modulus=bulk_array,
        density=density_array, gravity=gravity_array, frequency=frequency,
        is_solid_by_layer=[True, True, False, True], is_static_by_layer=[True, True, True, True],
        indices_by_layer=[ssls_layer_0_indx, ssls_layer_1_indx, ssls_layer_2_indx, ssls_layer_3_indx],
        order_l=2,
        surface_boundary_condition=None, solve_load_numbers=False,
        use_kamata=False,
        integrator='scipy', integration_method=None,
        integration_rtol=1.0e-6, integration_atol=1.0e-8,
        verbose=False, nondimensionalize=False, planet_bulk_density=bulk_density
        )

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex128, complex]

    # Test all dynamic
    tidal_y = tidal_y_solver(
        model_name='solid_solid_liquid_solid',
        radius=radius_array_to_use, shear_modulus=shear_array_ssls, bulk_modulus=bulk_array,
        density=density_array, gravity=gravity_array, frequency=frequency,
        is_solid_by_layer=[True, True, False, True], is_static_by_layer=[False, False, False, False],
        indices_by_layer=[ssls_layer_0_indx, ssls_layer_1_indx, ssls_layer_2_indx, ssls_layer_3_indx],
        order_l=2,
        surface_boundary_condition=None, solve_load_numbers=False,
        use_kamata=False,
        integrator='scipy', integration_method=None,
        integration_rtol=1.0e-6, integration_atol=1.0e-8,
        verbose=False, nondimensionalize=False, planet_bulk_density=bulk_density
        )

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex128, complex]

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex128, complex]

    # Test mix liq=Static
    tidal_y = tidal_y_solver(
        model_name='solid_solid_liquid_solid',
        radius=radius_array_to_use, shear_modulus=shear_array_ssls, bulk_modulus=bulk_array,
        density=density_array, gravity=gravity_array, frequency=frequency,
        is_solid_by_layer=[True, True, False, True], is_static_by_layer=[False, False, True, False],
        indices_by_layer=[ssls_layer_0_indx, ssls_layer_1_indx, ssls_layer_2_indx, ssls_layer_3_indx],
        order_l=2,
        surface_boundary_condition=None, solve_load_numbers=False,
        use_kamata=False,
        integrator='scipy', integration_method=None,
        integration_rtol=1.0e-6, integration_atol=1.0e-8,
        verbose=False, nondimensionalize=False, planet_bulk_density=bulk_density
        )

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex128, complex]

    # Test mix liq=Dynamic
    tidal_y = tidal_y_solver(
        model_name='solid_solid_liquid_solid',
        radius=radius_array_to_use, shear_modulus=shear_array_ssls, bulk_modulus=bulk_array,
        density=density_array, gravity=gravity_array, frequency=frequency,
        is_solid_by_layer=[True, True, False, True], is_static_by_layer=[True, False, True, True],
        indices_by_layer=[ssls_layer_0_indx, ssls_layer_1_indx, ssls_layer_2_indx, ssls_layer_3_indx],
        order_l=2,
        surface_boundary_condition=None, solve_load_numbers=False,
        use_kamata=False,
        integrator='scipy', integration_method=None,
        integration_rtol=1.0e-6, integration_atol=1.0e-8,
        verbose=False, nondimensionalize=False, planet_bulk_density=bulk_density
        )

    assert tidal_y.shape == (6, 10)
    assert tidal_y.dtype in [np.complex128, complex]

# def test_calculate_ssls_incomp():
#     """ Test the solution calculation for solid-solid-liquid-solid planet structure using the incompressible assumption """
#
#     # Test all static
#     tidal_y = tidal_y_solver(
#         model_name='solid_solid_liquid_solid',
#         radius=radius_array_to_use, shear_modulus=shear_array_ssls, bulk_modulus=bulk_array,
#         density=density_array, gravity=gravity_array, frequency=frequency,
#         is_solid_by_layer=[True, True, False, True], is_static_by_layer=[True, True, True, True],
#         indices_by_layer=[ssls_layer_0_indx, ssls_layer_1_indx, ssls_layer_2_indx, ssls_layer_3_indx],
#         order_l=2,
#         surface_boundary_condition=None, solve_load_numbers=False,
#         use_kamata=False,
#         integrator='scipy', integration_method=None,
#         integration_rtol=1.0e-6, integration_atol=1.0e-8,
#         verbose=False, nondimensionalize=False, planet_bulk_density=bulk_density, incompressible=True
#         )
#
#     assert tidal_y.shape == (6, 10)
#     assert tidal_y.dtype in [np.complex128, complex]
#
#     # Test all dynamic
#     tidal_y = tidal_y_solver(
#         model_name='solid_solid_liquid_solid',
#         radius=radius_array_to_use, shear_modulus=shear_array_ssls, bulk_modulus=bulk_array,
#         density=density_array, gravity=gravity_array, frequency=frequency,
#         is_solid_by_layer=[True, True, False, True], is_static_by_layer=[False, False, False, False],
#         indices_by_layer=[ssls_layer_0_indx, ssls_layer_1_indx, ssls_layer_2_indx, ssls_layer_3_indx],
#         order_l=2,
#         surface_boundary_condition=None, solve_load_numbers=False,
#         use_kamata=False,
#         integrator='scipy', integration_method=None,
#         integration_rtol=1.0e-6, integration_atol=1.0e-8,
#         verbose=False, nondimensionalize=False, planet_bulk_density=bulk_density, incompressible=True
#         )
#
#     assert tidal_y.shape == (6, 10)
#     assert tidal_y.dtype in [np.complex128, complex]
#
#     assert tidal_y.shape == (6, 10)
#     assert tidal_y.dtype in [np.complex128, complex]
#
#     # Test mix liq=Static
#     tidal_y = tidal_y_solver(
#         model_name='solid_solid_liquid_solid',
#         radius=radius_array_to_use, shear_modulus=shear_array_ssls, bulk_modulus=bulk_array,
#         density=density_array, gravity=gravity_array, frequency=frequency,
#         is_solid_by_layer=[True, True, False, True], is_static_by_layer=[False, False, True, False],
#         indices_by_layer=[ssls_layer_0_indx, ssls_layer_1_indx, ssls_layer_2_indx, ssls_layer_3_indx],
#         order_l=2,
#         surface_boundary_condition=None, solve_load_numbers=False,
#         use_kamata=False,
#         integrator='scipy', integration_method=None,
#         integration_rtol=1.0e-6, integration_atol=1.0e-8,
#         verbose=False, nondimensionalize=False, planet_bulk_density=bulk_density, incompressible=True
#         )
#
#     assert tidal_y.shape == (6, 10)
#     assert tidal_y.dtype in [np.complex128, complex]
#
#     # Test mix liq=Dynamic
#     tidal_y = tidal_y_solver(
#         model_name='solid_solid_liquid_solid',
#         radius=radius_array_to_use, shear_modulus=shear_array_ssls, bulk_modulus=bulk_array,
#         density=density_array, gravity=gravity_array, frequency=frequency,
#         is_solid_by_layer=[True, True, False, True], is_static_by_layer=[True, False, True, True],
#         indices_by_layer=[ssls_layer_0_indx, ssls_layer_1_indx, ssls_layer_2_indx, ssls_layer_3_indx],
#         order_l=2,
#         surface_boundary_condition=None, solve_load_numbers=False,
#         use_kamata=False,
#         integrator='scipy', integration_method=None,
#         integration_rtol=1.0e-6, integration_atol=1.0e-8,
#         verbose=False, nondimensionalize=False, planet_bulk_density=bulk_density, incompressible=True
#         )
#
#     assert tidal_y.shape == (6, 10)
#     assert tidal_y.dtype in [np.complex128, complex]
