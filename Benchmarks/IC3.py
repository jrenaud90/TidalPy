import numpy as np
import matplotlib.pyplot as plt


from TidalPy.constants import G
from TidalPy.tools.conversions import orbital_motion2semi_a
from TidalPy.rheology.complex_compliance.compliance_models import maxwell_array, sundberg_array
# Load TidalPy's multilayer functions
from TidalPy.tides.multilayer.numerical_int.initial_solution_dynamic import solid_guess_takeuchi, solid_guess_kamata
from TidalPy.utilities.graphics.multilayer.yplot import yplot
from TidalPy.tides.multilayer.numerical_int.radial_derivatives_dynamic import radial_derivatives_solid_general, radial_derivatives_liquid_general
from TidalPy.tides.multilayer.numerical_int.radial_derivatives_static import radial_derivatives_liquid_general as radial_derivatives_liquid_general_static
from scipy.integrate import solve_ivp
from TidalPy.utilities.numpy_helper import find_nearest

# Switches and Integration Variables
use_ts72_initial = False
use_incompressibility = False # To compare to RN08
static_core = True
integration_tol = 1.e-8
integration_method = 'RK45'
plot_tobie = True
plot_roberts = False

# Planet Structure
R_planet = 1600.e3
bulk_density = 3500.

viscosity_mantle = 1.e20
shear_mantle = 6.68250e10
if use_incompressibility:
    bulk_mantle = 1.e50
else:
    bulk_mantle = 1.2210e11

# From RN08 "Rigidity and viscosity have been reduced by a factor of 1e6 and 1e9 respectively from the overlying ice shell"
viscosity_OC = 0.
shear_OC = 0.
bulk_OC = 2.88e11

viscosity_IC = 1.e22
shear_IC = 6.68250e11
bulk_IC = 3.2210e11

# Orbital properties
planet_mass = 1.08e20
host_mass = 5.683e26
eccentricity = 0.0045
orbital_freq_HZ = 2. * np.pi * 5.308e-5  # This is reported in Hz in RN08. We will need to convert them to rad s-1
orbital_freq_TB_match = 2. * np.pi / (86400. * .9)
freq_europa = 2.04793e-05
freq_titan = 4.55938e-06

orbital_freq = 2. * np.pi / (86400. * .03)
semi_major_axis = orbital_motion2semi_a(orbital_freq, host_mass, planet_mass)

# Setup homogeneous domain
N = 300
radius_array = np.linspace(0., R_planet, N+1)[1:]
volume_array = np.zeros_like(radius_array)
volume_array[0] = (4. / 3.) * np.pi * radius_array[0]**3
volume_array[1:] = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)

mantle_density = 3100.

# Core
core_liq_density = 7000.
core_sol_density = 13000.
liq_r_frac = 0.4
sol_r_frac = 0.1

# Determine the radius structure based on the density
R_core = liq_r_frac * R_planet
R_core_sol = sol_r_frac * R_planet
span_IC = (radius_array[0], R_core_sol)
span_OC = (radius_array[radius_array > R_core_sol][0], R_core)
span_mantle = (radius_array[radius_array > R_core][0], radius_array[-1])

# Radius Index
core_sol_index = radius_array <= R_core_sol
core_liq_index = np.logical_and(radius_array <= R_core, radius_array > R_core_sol)
mantle_index = radius_array > R_core
R_IC_index = find_nearest(radius_array, R_core_sol)
R_OC_index = find_nearest(radius_array, R_core)

# Set the density for each layer
density_array = mantle_density * np.ones_like(radius_array)
density_array[core_liq_index] = core_liq_density
density_array[core_sol_index] = core_sol_density

# Radius Index

# Calculate masses, mass below, gravities
mass_array = volume_array * density_array

masses_below = np.zeros_like(radius_array)
masses_below[0] = mass_array[0]
masses_below[1:] = np.asarray([sum(mass_array[0:i+1]) for i in range(1, N)])

gravity_array = G * masses_below / radius_array**2

# Calculate rheological properties
viscosity_array = viscosity_mantle * np.ones_like(radius_array)
viscosity_array[core_sol_index] = viscosity_IC
viscosity_array[core_liq_index] = viscosity_OC

shear_array = shear_mantle * np.ones_like(radius_array)
shear_array[core_sol_index] = shear_IC
shear_array[core_liq_index] = shear_OC

bulk_array = bulk_mantle * np.ones_like(radius_array)
bulk_array[core_sol_index] = bulk_IC
bulk_array[core_liq_index] = bulk_OC

# Use the Maxwell rheology
complex_compliance_array = maxwell_array(orbital_freq, shear_array**(-1), viscosity_array)
complex_shear_array = complex_compliance_array**(-1)

# Solve Rocky Inner Core
interp_funcs = (
        lambda x: np.interp(x, radius_array, complex_shear_array),
        lambda x: np.interp(x, radius_array, bulk_array),
        lambda x: np.interp(x, radius_array, density_array),
        lambda x: np.interp(x, radius_array, gravity_array)
)

def deriv_sol(x, y):
    shear_modulus, bulk_modulus, density, gravity = interp_funcs

    result = radial_derivatives_solid_general(x, y, shear_modulus(x), bulk_modulus(x), density(x), gravity(x), orbital_freq)

    return result

def deriv_liq(x, y):
    shear_modulus, bulk_modulus, density, gravity = interp_funcs

    if static_core:
        result = radial_derivatives_liquid_general_static(x, y, density(x), gravity(x))
    else:
        result = radial_derivatives_liquid_general(x, y, bulk_modulus(x), density(x), gravity(x), orbital_freq)

    return result

# Get initial solutions for the IC
if use_ts72_initial:
    initial_solutions_solid_IC = solid_guess_takeuchi(
                    radius_array[core_sol_index],
                    complex_shear_array[core_sol_index],
                    bulk_array[core_sol_index],
                    density_array[core_sol_index], orbital_freq, order_l=2)
else:
    initial_solutions_solid_IC = solid_guess_kamata(
            radius_array[core_sol_index],
            complex_shear_array[core_sol_index],
            bulk_array[core_sol_index],
            density_array[core_sol_index], orbital_freq, order_l=2)

# Solve for rocky IC
solution_ys_IC = list()
solution_rs_IC = list()
for sn in range(3):
    print(f'Solving IC Solution {sn}')
    initial_value = initial_solutions_solid_IC[sn][:, 0]

    # Boundary conditions must be set manually
    # Roberts & Nimmo (2008): Liquid innermost zone.
    # initial_value[1] = 1.
    # initial_value[3] = 1.
    # initial_value[5] = 1.

    # Henning and Hurford (2014) Solid innermost
    # initial_value[0] = 1.
    # initial_value[1] = 1.
    # initial_value[2] = 1.

    # BC for Tobie 2005 benchmark (From IcyDwarf)
    # initial_value[0] = 0.05
    # initial_value[1] = 0.01
    # initial_value[5] = 1.

    # Use the initial guess function.
    pass

    # Make sure the initial value remains complex.
    initial_value = np.asarray(initial_value, dtype=np.complex128)
    solution = solve_ivp(deriv_sol, span_IC, initial_value,
                         method=integration_method, rtol=integration_tol,
                         t_eval=radius_array[core_sol_index])
    if solution.status == 0:
        print(f'Solution {sn} solved.')
        solution_rs_IC.append(solution.t)
        solution_ys_IC.append(solution.y)
    else:
        print(solution.message)

# Solve Outer Core
solution_ys_OC = list()
solution_rs_OC = list()
if static_core:
    OC_sn_range = range(1)
else:
    OC_sn_range = range(2)

for sn in OC_sn_range:
    print(f'Solving OC Solution {sn+1} of {len(OC_sn_range)}')

    if static_core:
        initial_value = np.zeros(2, dtype=np.complex128)
    else:
        initial_value = np.zeros(4, dtype=np.complex128)

    # Boundary conditions are based on the inner core results.
    if sn == 0:
        if static_core:
            liquid_density = density_array[R_IC_index + 1]
            interface_gravity = gravity_array[R_IC_index]

            # initial_value[0] = solution_ys_IC[0][4, -1]
            # initial_value[1] = solution_ys_IC[0][5, -1] - \
            #                    (4. * np.pi * G * liquid_density / interface_gravity) * initial_value[0]

            y4_frac_1 = solution_ys_IC[0][3, -1] / solution_ys_IC[2][3, -1]
            y4_frac_2 = solution_ys_IC[1][3, -1] / solution_ys_IC[2][3, -1]

            # gamma_j = (y_2j - f_j y_23) - rho( g(y_1j - f_j y_13) - (y_5j - f_j y_53))
            gamma_1 = (solution_ys_IC[0][1, -1] - y4_frac_1 * solution_ys_IC[2][1, -1]) - \
                      liquid_density * (interface_gravity * (solution_ys_IC[0][0, -1] - y4_frac_1 * solution_ys_IC[2][0, -1]) -
                                (solution_ys_IC[0][4, -1] - y4_frac_1 * solution_ys_IC[2][4, -1]))
            gamma_2 = (solution_ys_IC[1][1, -1] - y4_frac_2 * solution_ys_IC[2][1, -1]) - \
                      liquid_density * (interface_gravity * (solution_ys_IC[1][0, -1] - y4_frac_2 * solution_ys_IC[2][0, -1]) -
                                (solution_ys_IC[1][4, -1] - y4_frac_2 * solution_ys_IC[2][4, -1]))

            initial_value[0] = solution_ys_IC[0][4, -1] - (gamma_1 / gamma_2) * solution_ys_IC[1][4, -1] - \
                               (y4_frac_1 - (gamma_1 / gamma_2) * y4_frac_2) * solution_ys_IC[2][4, -1]

            # initial_value[1] = solution_ys_IC[0][5, -1] - \
            #                    (4. * np.pi * G * liquid_density / interface_gravity) * initial_value[0]

            y_7_IC_0 = solution_ys_IC[0][5, -1] + (4. * np.pi * G / interface_gravity) * solution_ys_IC[0][1, -1]
            y_7_IC_1 = solution_ys_IC[1][5, -1] + (4. * np.pi * G / interface_gravity) * solution_ys_IC[1][1, -1]
            y_7_IC_2 = solution_ys_IC[2][5, -1] + (4. * np.pi * G / interface_gravity) * solution_ys_IC[2][1, -1]

            initial_value[1] = y_7_IC_0 - (gamma_1 / gamma_2) * y_7_IC_1 - \
                               (y4_frac_1 - (gamma_1 / gamma_2) * y4_frac_2) * y_7_IC_2
        else:
            sol1_frac = solution_ys_IC[0][3, -1] / solution_ys_IC[2][3, -1]
            initial_value[0] = solution_ys_IC[0][0, -1] - sol1_frac * solution_ys_IC[2][0, -1]
            initial_value[1] = solution_ys_IC[0][1, -1] - sol1_frac * solution_ys_IC[2][1, -1]
            initial_value[2] = solution_ys_IC[0][4, -1] - sol1_frac * solution_ys_IC[2][4, -1]
            initial_value[3] = solution_ys_IC[0][5, -1] - sol1_frac * solution_ys_IC[2][5, -1]

    elif sn == 1:
        sol2_frac = solution_ys_IC[1][3, -1] / solution_ys_IC[2][3, -1]
        initial_value[0] = solution_ys_IC[1][0, -1] - sol2_frac * solution_ys_IC[2][0, -1]
        initial_value[1] = solution_ys_IC[1][1, -1] - sol2_frac * solution_ys_IC[2][1, -1]
        initial_value[2] = solution_ys_IC[1][4, -1] - sol2_frac * solution_ys_IC[2][4, -1]
        initial_value[3] = solution_ys_IC[1][5, -1] - sol2_frac * solution_ys_IC[2][5, -1]

    # Make sure the initial value remains complex.
    initial_value = np.asarray(initial_value, dtype=np.complex128)
    solution = solve_ivp(deriv_liq, span_OC, initial_value,
                         method=integration_method, rtol=integration_tol,
                         t_eval=radius_array[core_liq_index])
    if solution.status == 0:
        print(f'Solution {sn} solved.')
        solution_rs_OC.append(solution.t)
        solution_ys_OC.append(solution.y)
    else:
        print(solution.message)

# Solve for Mantle
solution_ys_mantle = list()
solution_rs_mantle = list()
for sn in range(3):
    print(f'Solving Mantle Solution {sn}')
    initial_value = initial_solutions_solid_IC[sn][:, 0]

    # Boundary conditions are based on the outer core results.
    if static_core:
        # From Saito74
        liquid_density = density_array[R_OC_index]
        interface_gravity = gravity_array[R_OC_index]
        if sn == 0:
            initial_value[0] = 0.
            initial_value[1] = -liquid_density * solution_ys_OC[0][0, -1]
            initial_value[4] = solution_ys_OC[0][0, -1]
            initial_value[5] = solution_ys_OC[0][1, -1] + \
                               (4. * np.pi * G * liquid_density / interface_gravity) * \
                               solution_ys_OC[0][0, -1]

            initial_value[2] = 0.
            initial_value[3] = 0.

        elif sn == 1:
            initial_value[0] = 1.
            initial_value[1] = liquid_density * interface_gravity * initial_value[0]
            initial_value[5] = -4. * np.pi * G * liquid_density * initial_value[0]

            initial_value[2] = 0.
            initial_value[3] = 0.
            initial_value[4] = 0.

        elif sn == 2:
            initial_value[0] = 0.
            initial_value[1] = 0.
            initial_value[2] = 1.
            initial_value[3] = 0.
            initial_value[4] = 0.
            initial_value[5] = 0.
    else:
        if sn == 0:
            initial_value[0] = solution_ys_OC[0][0, -1]
            initial_value[1] = solution_ys_OC[0][1, -1]
            initial_value[4] = solution_ys_OC[0][2, -1]
            initial_value[5] = solution_ys_OC[0][3, -1]

            initial_value[2] = 0.
            initial_value[3] = 0.

        elif sn == 1:
            initial_value[0] = solution_ys_OC[1][0, -1]
            initial_value[1] = solution_ys_OC[1][1, -1]
            initial_value[4] = solution_ys_OC[1][2, -1]
            initial_value[5] = solution_ys_OC[1][3, -1]

            initial_value[2] = 0.
            initial_value[3] = 0.

        elif sn == 2:
            initial_value[0] = 0.
            initial_value[1] = 0.
            initial_value[2] = 1.
            initial_value[3] = 0.
            initial_value[4] = 0.
            initial_value[5] = 0.

    # Make sure the initial value remains complex.
    initial_value = np.asarray(initial_value, dtype=np.complex128)
    solution = solve_ivp(deriv_sol, span_mantle, initial_value,
                         method=integration_method, rtol=integration_tol,
                         t_eval=radius_array[mantle_index])
    if solution.status == 0:
        print(f'Solution {sn} solved.')
        solution_rs_mantle.append(solution.t)
        solution_ys_mantle.append(solution.y)
    else:
        print(solution.message)

# Solve unknowns from top down.
# Build solution matrix at surface
sol_surf_mtx = np.asarray([
    [solution_ys_mantle[0][1, -1], solution_ys_mantle[1][1, -1], solution_ys_mantle[2][1, -1]],
    [solution_ys_mantle[0][3, -1], solution_ys_mantle[1][3, -1], solution_ys_mantle[2][3, -1]],
    [solution_ys_mantle[0][5, -1], solution_ys_mantle[1][5, -1], solution_ys_mantle[2][5, -1]]
])
sol_surf_mtx_inv = np.linalg.inv(sol_surf_mtx)
surf_bc = np.zeros(3, dtype=np.complex128)
surf_bc[2] = 5. / R_planet
Q_mantle_vector = sol_surf_mtx_inv @ surf_bc

if static_core:
    # Solve for the outer core Qs
    Q_OC_vector = np.zeros(1, dtype=np.complex128)
    Q_OC_vector[0] = Q_mantle_vector[0]

    # Solve for inner core Qs
    Q_IC_vector = np.zeros(3, dtype=np.complex128)

    y4_frac_1 = solution_ys_IC[0][3, -1] / solution_ys_IC[2][3, -1]
    y4_frac_2 = solution_ys_IC[1][3, -1] / solution_ys_IC[2][3, -1]

    g_IC = gravity_array[R_IC_index]      # Interface gravity
    rho_IC = density_array[R_IC_index+1]  # Liquid density
    print(g_IC, rho_IC)

    # gamma_j = (y_2j - f_j y_23) - rho( g(y_1j - f_j y_13) - (y_5j - f_j y_53))
    gamma_1 = (solution_ys_IC[0][1, -1] - y4_frac_1 * solution_ys_IC[2][1, -1]) - \
              rho_IC * (g_IC * (solution_ys_IC[0][0, -1] - y4_frac_1 * solution_ys_IC[2][0, -1]) -
                        (solution_ys_IC[0][4, -1] - y4_frac_1 * solution_ys_IC[2][4, -1]))
    gamma_2 = (solution_ys_IC[1][1, -1] - y4_frac_2 * solution_ys_IC[2][1, -1]) - \
              rho_IC * (g_IC * (solution_ys_IC[1][0, -1] - y4_frac_2 * solution_ys_IC[2][0, -1]) -
                        (solution_ys_IC[1][4, -1] - y4_frac_2 * solution_ys_IC[2][4, -1]))

    Q_IC_vector[0] = Q_OC_vector[0]
    Q_IC_vector[1] = (-gamma_1 / gamma_2) * Q_IC_vector[0]
    Q_IC_vector[2] = -y4_frac_1 * Q_IC_vector[0] - y4_frac_2 * Q_IC_vector[1]

    # Solve for total planet y's
    tidal_y_OC = Q_OC_vector[0] * solution_ys_OC[0]

    core_ys = (
        np.nan * np.ones_like(tidal_y_OC[0, :]),
        np.nan * np.ones_like(tidal_y_OC[0, :]),
        np.nan * np.ones_like(tidal_y_OC[0, :]),
        np.nan * np.ones_like(tidal_y_OC[0, :]),
        tidal_y_OC[0, :],
        np.nan * np.ones_like(tidal_y_OC[0, :])
    )

else:
    # Solve for the outer core Qs
    Q_OC_vector = np.zeros(2, dtype=np.complex128)
    Q_OC_vector[0] = Q_mantle_vector[0]
    Q_OC_vector[1] = Q_mantle_vector[1]

    # Solve for inner core Qs
    Q_IC_vector = np.zeros(3, dtype=np.complex128)
    Q_IC_vector[0] = Q_OC_vector[0]
    Q_IC_vector[1] = Q_OC_vector[1]
    y4_frac_1 = solution_ys_IC[0][3, -1] / solution_ys_IC[2][3, -1]
    y4_frac_2 = solution_ys_IC[1][3, -1] / solution_ys_IC[2][3, -1]
    Q_IC_vector[2] = -y4_frac_1 * Q_IC_vector[0] - y4_frac_2 * Q_IC_vector[1]

    # Solve for total planet y's
    tidal_y_IC = Q_IC_vector[0] * solution_ys_IC[0] + Q_IC_vector[1] * solution_ys_IC[1] + \
                 Q_IC_vector[2] * solution_ys_IC[2]

    tidal_y_OC = Q_OC_vector[0] * solution_ys_OC[0] + Q_OC_vector[1] * solution_ys_OC[1]

    tidal_y_mantle = Q_mantle_vector[0] * solution_ys_mantle[0] + Q_mantle_vector[1] * solution_ys_mantle[1] + \
                     Q_mantle_vector[2] * solution_ys_mantle[2]

    # Outer core is missing two y's, fix that now.
    y3_OC = \
        (1. / (orbital_freq**2 * density_array[core_liq_index] * radius_array[core_liq_index])) * \
        (density_array[core_liq_index] * gravity_array[core_liq_index] * tidal_y_OC[0, :] -
         tidal_y_OC[1, :] - density_array[core_liq_index] * tidal_y_OC[2, :])

    core_ys = (
        tidal_y_OC[0, :],
        tidal_y_OC[1, :],
        y3_OC,
        np.nan * np.ones_like(tidal_y_OC[0, :]),
        tidal_y_OC[2, :],
        tidal_y_OC[3, :]
    )
tidal_y_OC_full = np.vstack(core_ys)

# Solve for total planet y's
tidal_y_IC = Q_IC_vector[0] * solution_ys_IC[0] + Q_IC_vector[1] * solution_ys_IC[1] + \
             Q_IC_vector[2] * solution_ys_IC[2]

tidal_y_mantle = Q_mantle_vector[0] * solution_ys_mantle[0] + Q_mantle_vector[1] * solution_ys_mantle[1] + \
                 Q_mantle_vector[2] * solution_ys_mantle[2]

# Combine solutions for all layers
tidal_y = np.concatenate((tidal_y_IC, core_ys, tidal_y_mantle), axis=1)

tidal_y[4, :] = -np.imag(tidal_y[4, :] - 1)

yplot(tidal_y, radius_array, use_tobie_limits=True)