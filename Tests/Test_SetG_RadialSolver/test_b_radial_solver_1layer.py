import pytest

import numpy as np

import TidalPy
TidalPy.test_mode()


from TidalPy.RadialSolver import radial_solver

frequency = 1.0 / (86400.0 * 6.0)
N = 100
radius_array = np.linspace(0.1, 6000.e3, N)
density_array = 3500. * np.ones_like(radius_array)
bulk_modulus_array = 1.0e11 * np.ones_like(radius_array)
viscosity_array = 1.0e20 * np.ones_like(radius_array)
shear_array = 5.0e10 * np.ones_like(radius_array)

from TidalPy.rheology.models import MaxwellRheology
complex_shear_modulus_array = np.empty(N, dtype=np.complex128)
max_rho = MaxwellRheology()
max_rho.vectorize_modulus_viscosity(frequency, shear_array, viscosity_array, complex_shear_modulus_array)


from TidalPy.utilities.spherical_helper import calculate_mass_gravity_arrays
volume_array, mass_array, gravity_array = calculate_mass_gravity_arrays(radius_array, density_array)

planet_bulk_density = np.sum(mass_array) / np.sum(volume_array)
upper_radius_by_layer = (radius_array[-1],)

@pytest.mark.parametrize('is_solid', (True, False))
@pytest.mark.parametrize('is_static', (True, False))
@pytest.mark.parametrize('is_incompressible', (True, False))
def test_radial_solver_1layer(is_solid, is_static, is_incompressible):

    is_solid_by_layer = (is_solid,)
    is_static_by_layer = (is_static,)
    is_incompressible_by_layer = (is_incompressible,)
    
    out = radial_solver(
        radius_array, density_array, gravity_array, bulk_modulus_array, complex_shear_modulus_array,
        frequency, planet_bulk_density, 
        is_solid_by_layer, is_static_by_layer, is_incompressible_by_layer, upper_radius_by_layer,
        degree_l=2, solve_for=None, use_kamata=False,
        integration_method=2, integration_rtol=1.0e-4, integration_atol=1.0e-12,
        scale_rtols_by_layer_type=True,
        max_num_steps=500_000, expected_size=250, max_step=0,
        limit_solution_to_radius=True, verbose=False, nondimensionalize=False)

