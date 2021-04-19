import numpy as np
import matplotlib.pyplot as plt


from TidalPy.constants import G
from TidalPy.tools.conversions import orbital_motion2semi_a
from TidalPy.rheology.complex_compliance.compliance_models import maxwell_array, sundberg_array
# Load TidalPy's multilayer functions
from TidalPy.tides.multilayer.numerical_int.initial_solution import solid_solutions
from TidalPy.utilities.graphics.multilayer.yplot import yplot
from TidalPy.tides.multilayer.numerical_int.radial_functions_solid import radial_func_derivatives_general, radial_func_derivatives_general_liq

full_domain = True

# Planet Structure
R_planet = 1600.e3
rho_mantle = 3300.
# Two different core density models, lower value assumes Fe / Fe-S core composition.
rho_core = 8000.
rho_average = 3500.

use_incompressibility = False # To compare to RN08

viscosity_mantle = 1.e20
shear_mantle = 6.68250e10
if use_incompressibility:
    bulk_mantle = 1.e50
else:
    bulk_mantle = 1.2210e11

# From RN08 "Rigidity and viscosity have been reduced by a factor of 1e6 and 1e9 respectively from the overlying ice shell"
viscosity_core = 1.e13 / 1.e9
shear_core = 1.e-1

if use_incompressibility:
    bulk_core = 1.e50
else:
    bulk_core = 2.88e11

# Core radius is determined from the core density and planet radius
R_core = ((rho_average - rho_mantle) / (rho_core - rho_mantle))**(1 / 3) * R_planet

# Orbital properties
planet_mass = 1.08e20
host_mass = 5.683e26
eccentricity = 0.0045
orbital_freq_HZ = 5.308e-5  # This is reported in Hz in RN08. We will need to convert them to rad s-1
orbital_freq = 2. * np.pi * orbital_freq_HZ
semi_major_axis = orbital_motion2semi_a(orbital_freq, host_mass, planet_mass)

# Setup homogeneous domain
N = 100
radius_array = np.linspace(0., R_planet, N+1)[1:]
volume_array = np.zeros_like(radius_array)
volume_array[0] = (4. / 3.) * np.pi * radius_array[0]**3
volume_array[1:] = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)

density_array = rho_mantle * np.ones_like(radius_array)
density_array[radius_array <= R_core] = rho_core

# Calculate masses, mass below, gravities
mass_array = volume_array * density_array

masses_below = np.zeros_like(radius_array)
masses_below[0] = mass_array[0]
masses_below[1:] = np.asarray([sum(mass_array[0:i+1]) for i in range(1, N)])

gravity_array = G * masses_below / radius_array**2

# Calculate rheological properties
viscosity_array = viscosity_mantle * np.ones_like(radius_array)
viscosity_array[radius_array <= R_core] = viscosity_core

shear_array = shear_mantle * np.ones_like(radius_array)
shear_array[radius_array <= R_core] = shear_core
bulk_array = bulk_mantle * np.ones_like(radius_array)
bulk_array[radius_array <= R_core] = bulk_core

# Use the Maxwell rheology
complex_compliance_array = maxwell_array(orbital_freq, shear_array**(-1), viscosity_array)
complex_shear_array = complex_compliance_array**(-1)

interp_funcs = (lambda x: np.interp(x, radius_array, complex_shear_array),
           lambda x: np.interp(x, radius_array, bulk_array),
           lambda x: np.interp(x, radius_array, density_array),
           lambda x: np.interp(x, radius_array, gravity_array))

def deriv_sol(x, y):
    shear_modulus, bulk_modulus, density, gravity = interp_funcs

    result = radial_func_derivatives_general(x, y, shear_modulus(x), bulk_modulus(x), density(x), gravity(x), orbital_freq)

    return result

def deriv_liq(x, y):
    shear_modulus, bulk_modulus, density, gravity = interp_funcs

    result = radial_func_derivatives_general_liq(x, y, shear_modulus(x), bulk_modulus(x), density(x), gravity(x), orbital_freq)

    return result


initial_solutions = solid_solutions(radius_array[radius_array > R_core],
                                    bulk_array[radius_array > R_core],
                                    complex_shear_array[radius_array > R_core],
                                    density_array[radius_array > R_core],
                                    gravity_array[radius_array > R_core], orbital_freq, order_l=2)


# IVP Approach
from scipy.integrate import solve_ivp
span = (radius_array[radius_array > R_core][0], R_planet)

solution_ys = list()
solution_rs = list()
for sn in range(3):
    print(f'Solving Solution {sn}')
    initial_value = initial_solutions[sn][:, 0]

    # Apply bottom BC (is this right??!)
    initial_value[0] = 0.
    initial_value[2] = 0.
    initial_value[4] = 0.

    solution = solve_ivp(deriv_sol, span, initial_value, t_eval=radius_array[radius_array > R_core])
    if solution.status == 0:
        print(f'Solution {sn} solved.')
        solution_rs.append(solution.t)
        solution_ys.append(solution.y)
    else:
        print(solution.message)

# Build solution matrix at surface
sol_surf_mtx = np.asarray([
    [solution_ys[0][1, -1], solution_ys[1][1, -1], solution_ys[2][1, -1]],
    [solution_ys[0][3, -1], solution_ys[1][3, -1], solution_ys[2][3, -1]],
    [solution_ys[0][5, -1], solution_ys[1][5, -1], solution_ys[2][5, -1]]
])
sol_surf_mtx_inv = np.linalg.inv(sol_surf_mtx)
surf_bc = np.zeros(3, dtype=np.complex128)
surf_bc[2] = 5. / R_planet
Q_vector = sol_surf_mtx_inv @ surf_bc

tidal_y_sol = Q_vector[0] * solution_ys[0] + Q_vector[1] * solution_ys[1] + Q_vector[2] * solution_ys[2]
tidal_y_ivp = np.copy(tidal_y_sol)
tidal_y_ivp_r = np.copy(solution_rs[-1])

solution_ys.append(tidal_y_ivp)
solution_rs.append(tidal_y_ivp_r)

if len(solution_ys) != 0:
    yplot(solution_ys, solution_rs)

# BVP Approach
from scipy.integrate import solve_bvp
surf_condition = (5. / R_planet)
def bc_residual(ya, yb):
    return (ya[0], yb[1], ya[2], yb[3], ya[4], yb[5] - surf_condition)

solution_ys = list()
solution_rs = list()
for sn in range(3):

    initial_guess = initial_solutions[sn]
    solution = solve_bvp(deriv_sol, bc_residual, radius_array[radius_array > R_core], initial_guess, tol=1e-3, verbose=2)

    print(f'Solving Solution {sn}')
    #     solution = solve_ivp(deriv_sol, span, initial_solutions[sn][:, 0])
    if solution.status == 0:
        print(f'Solution {sn} solved.')
        solution_rs.append(solution.x)
        solution_ys.append(solution.y)
    else:
        print(solution.message)

solution_ys.append(tidal_y_ivp)
solution_rs.append(tidal_y_ivp_r)

if len(solution_ys) != 0:
    yplot(solution_ys, solution_rs)