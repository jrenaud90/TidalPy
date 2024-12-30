""" Tests for extracting useful information out of a multilayer tidal propagation
"""

import pytest
import numpy as np

import TidalPy


from TidalPy.constants import G
from TidalPy.tides.multilayer.heating import calc_radial_volumetric_tidal_heating, calc_radial_volumetric_tidal_heating_from_rs_solution
from TidalPy.tides.multilayer.sensitivity import calc_sensitivity_to_shear, calc_sensitivity_to_bulk
from TidalPy.RadialSolver import radial_solver
from TidalPy.utilities.conversions import orbital_motion2semi_a

# Model planet - 2layers
N = 10
planet_bulk_density = 5000.
density_array = planet_bulk_density * np.ones(N)
radius_array = np.linspace(0., 1.e6, N)
upper_radius_bylayer_array = np.asarray((radius_array[-1],), dtype=np.float64)
volume_array = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
mass_array = volume_array * density_array[1:]
planet_mass = sum(mass_array)
mass_below = np.asarray([np.sum(mass_array[:i + 1]) for i in range(N-1)])
gravity_array = G * mass_below / (radius_array[1:]**2)
shear_array = 5.e10 * np.ones(N, dtype=np.complex128)
bulk_array = 5.e12 * np.ones(N, dtype=np.complex128)
host_mass = 10. * planet_mass
orbital_freq = (2. * np.pi / (86400. * 6.))
semi_major_axis = orbital_motion2semi_a(orbital_freq, host_mass, planet_mass)
eccentricity = 0.01


@pytest.mark.parametrize('degree_l', (2, 3))
@pytest.mark.parametrize('use_prop_matrix', (True, False))
@pytest.mark.parametrize('do_shear', (True, False))
def test_radial_sensitivity_functions(degree_l, use_prop_matrix, do_shear):
    
    if use_prop_matrix:
        radial_solution = radial_solver(
            radius_array,
            density_array,
            bulk_array,
            shear_array,
            orbital_freq,
            planet_bulk_density,
            ('solid',),
            (True,),
            (True,),
            upper_radius_bylayer_array,
            degree_l = degree_l,
            use_prop_matrix = use_prop_matrix,
            perform_checks = True
            )
    else:
        radial_solution = radial_solver(
            radius_array,
            density_array,
            bulk_array,
            shear_array,
            orbital_freq,
            planet_bulk_density,
            ('solid',),
            (False,),
            (False,),
            upper_radius_bylayer_array,
            degree_l = degree_l,
            use_prop_matrix = use_prop_matrix,
            perform_checks = True
            )
    assert radial_solution.success

    # Find sensitivities
    if do_shear:
        sensitivity_array = calc_sensitivity_to_shear(
            radial_solution.result, 
            radius_array,
            shear_array,
            bulk_array,
            degree_l = degree_l,
            perform_checks = True
            )
    else:
        sensitivity_array = calc_sensitivity_to_bulk(
            radial_solution.result, 
            radius_array,
            shear_array,
            bulk_array,
            degree_l = degree_l,
            perform_checks = True
            )


    assert sensitivity_array.dtype == np.float64
    assert sensitivity_array.shape == (radius_array.size,)


@pytest.mark.parametrize('degree_l', (2, 3))
@pytest.mark.parametrize('use_prop_matrix', (True, False))
@pytest.mark.parametrize('shear_only', (True, False))
def test_calc_radial_volumetric_tidal_heating(degree_l, use_prop_matrix, shear_only):
    
    if use_prop_matrix:
        radial_solution = radial_solver(
            radius_array,
            density_array,
            bulk_array,
            shear_array,
            orbital_freq,
            planet_bulk_density,
            ('solid',),
            (True,),
            (True,),
            upper_radius_bylayer_array,
            degree_l = degree_l,
            use_prop_matrix = use_prop_matrix,
            perform_checks = True
            )
    else:
        radial_solution = radial_solver(
            radius_array,
            density_array,
            bulk_array,
            shear_array,
            orbital_freq,
            planet_bulk_density,
            ('solid',),
            (False,),
            (False,),
            upper_radius_bylayer_array,
            degree_l = degree_l,
            use_prop_matrix = use_prop_matrix,
            perform_checks = True
            )
    assert radial_solution.success

    # Find sensitivities
    shear_sensitivity_array = calc_sensitivity_to_shear(
        radial_solution.result, 
        radius_array,
        shear_array,
        bulk_array,
        degree_l = degree_l,
        perform_checks = True
        )
    if shear_only:
        bulk_sensitivity_array = np.zeros(shear_sensitivity_array.size, dtype=np.float64, order='C')
    else:
        bulk_sensitivity_array = calc_sensitivity_to_bulk(
            radial_solution.result, 
            radius_array,
            shear_array,
            bulk_array,
            degree_l = degree_l,
            perform_checks = True
            )
        
    # Calculate heating
    eccentricity=0.1
    orbital_frequency=0.001
    semi_major_axis=1.0e5
    host_mass=1.0e23
    volumetric_heating = calc_radial_volumetric_tidal_heating(
        eccentricity,
        orbital_frequency,
        semi_major_axis,
        host_mass,
        radius_array,
        shear_sensitivity_array,
        shear_array,
        bulk_sensitivity_array,
        bulk_array,
        degree_l = degree_l,
        G_to_use = G,
        perform_checks = True
        )


    assert volumetric_heating.dtype == np.float64
    assert volumetric_heating.shape == (radius_array.size,)
    assert np.all(~np.isnan(volumetric_heating))

@pytest.mark.parametrize('degree_l', (2, 3))
@pytest.mark.parametrize('use_prop_matrix', (True, False))
@pytest.mark.parametrize('shear_only', (True, False))
def test_calc_radial_volumetric_tidal_heating(degree_l, use_prop_matrix, shear_only):
    
    if use_prop_matrix:
        radial_solution = radial_solver(
            radius_array,
            density_array,
            bulk_array,
            shear_array,
            orbital_freq,
            planet_bulk_density,
            ('solid',),
            (True,),
            (True,),
            upper_radius_bylayer_array,
            degree_l = degree_l,
            use_prop_matrix = use_prop_matrix,
            perform_checks = True
            )
    else:
        radial_solution = radial_solver(
            radius_array,
            density_array,
            bulk_array,
            shear_array,
            orbital_freq,
            planet_bulk_density,
            ('solid',),
            (False,),
            (False,),
            upper_radius_bylayer_array,
            degree_l = degree_l,
            use_prop_matrix = use_prop_matrix,
            perform_checks = True
            )
    assert radial_solution.success
        
    # Calculate heating
    eccentricity=0.1
    orbital_frequency=0.001
    semi_major_axis=1.0e5
    host_mass=1.0e23
    volumetric_heating = calc_radial_volumetric_tidal_heating_from_rs_solution(
        eccentricity,
        orbital_frequency,
        semi_major_axis,
        host_mass,
        radial_solution,
        perform_checks = True
        )

    assert volumetric_heating.dtype == np.float64
    assert volumetric_heating.shape == (radius_array.size,)
    assert np.all(~np.isnan(volumetric_heating))
