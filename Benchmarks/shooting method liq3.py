import numpy as np
import matplotlib.pyplot as plt


from TidalPy.constants import G
from TidalPy.tools.conversions import orbital_motion2semi_a
from TidalPy.utilities.numpy_helper import find_nearest
from TidalPy.rheology.complex_compliance.compliance_models import maxwell_array, sundberg_array
# Load TidalPy's multilayer functions
from TidalPy.tides.multilayer.numerical_int.initial_solution_dynamic import \
    solid_guess_takeuchi, liquid_guess_takeuchi, solid_guess_kamata, liquid_guess_kamata
from TidalPy.tides.multilayer.numerical_int.initial_solution_static import liquid_guess_saito
from TidalPy.utilities.graphics.multilayer.yplot import yplot
from TidalPy.tides.multilayer.numerical_int.radial_derivatives_static import radial_derivatives_liquid_general as radial_derivatives_liquid_general_static
from TidalPy.tides.multilayer.numerical_int.radial_derivatives_dynamic import radial_derivatives_solid_general, radial_derivatives_liquid_general

# Switches and Integration Variables
use_ts72_initial = False
static_core = True
use_incompressibility = False # To compare to RN08
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

if use_ts72_initial:
    solid_solutions = solid_guess_takeuchi
    liquid_solutions = liquid_guess_takeuchi
    liquid_solutions_static = liquid_guess_saito
else:
    solid_solutions = solid_guess_kamata
    liquid_solutions = liquid_guess_kamata
    liquid_solutions_static = liquid_guess_saito

# From RN08 "Rigidity and viscosity have been reduced by a factor of 1e6 and 1e9 respectively from the overlying ice shell"
viscosity_core = 0.
shear_core = 0.

if use_incompressibility:
    bulk_core = 1.e50
else:
    bulk_core = 2.88e11

# Orbital properties
planet_mass = 1.08e20
host_mass = 5.683e26
eccentricity = 0.0045
orbital_freq_HZ = 2. * np.pi * 5.308e-5  # This is reported in Hz in RN08. We will need to convert them to rad s-1
orbital_freq_TB_match = 2. * np.pi / (86400. * 0.9)
freq_europa = 2.04793e-05
freq_titan = 4.55938e-06

orbital_freq = orbital_freq_TB_match
semi_major_axis = orbital_motion2semi_a(orbital_freq, host_mass, planet_mass)

# Setup homogeneous domain
N = 100
radius_array = np.linspace(0., R_planet, N+1)[1:]
volume_array = np.zeros_like(radius_array)
volume_array[0] = (4. / 3.) * np.pi * radius_array[0]**3
volume_array[1:] = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)

# Setup 3 models based on the three layer structures used in Tobie 2005
models = {'TB05-Homogen': (None, 3500.),
          'TB05-Core1': (5150, 3300.),
          'TB05-Core2': (8000, 3300.)}
model_radii = list()
model_names = list()
model_results = list()

for model_name, (core_density, mantle_density) in models.items():

    print(f'\n\nWorking on model {model_name}.\n------------------')

    # Determine the radius structure based on the density
    if core_density is None:
        R_core = None
    else:
        R_core = ((bulk_density - mantle_density) / (core_density - mantle_density))**(1 / 3) * R_planet

    # Set the density for each layer
    density_array = mantle_density * np.ones_like(radius_array)
    if R_core is not None:
        density_array[radius_array <= R_core] = core_density

    # Calculate masses, mass below, gravities
    mass_array = volume_array * density_array

    masses_below = np.zeros_like(radius_array)
    masses_below[0] = mass_array[0]
    masses_below[1:] = np.asarray([sum(mass_array[0:i+1]) for i in range(1, N)])

    gravity_array = G * masses_below / radius_array**2

    # Calculate rheological properties
    viscosity_array = viscosity_mantle * np.ones_like(radius_array)
    if R_core is not None:
        viscosity_array[radius_array <= R_core] = viscosity_core

    shear_array = shear_mantle * np.ones_like(radius_array)
    bulk_array = bulk_mantle * np.ones_like(radius_array)
    if R_core is not None:
        shear_array[radius_array <= R_core] = shear_core
        bulk_array[radius_array <= R_core] = bulk_core

    # Use the Maxwell rheology
    complex_compliance_array = maxwell_array(orbital_freq, shear_array**(-1), viscosity_array)
    complex_shear_array = complex_compliance_array**(-1)

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

        mu = shear_modulus(x)
        K = bulk_modulus(x)
        rho = density(x)
        g = gravity(x)

        if static_core:
            result = radial_derivatives_liquid_general_static(x, y, rho, g)
        else:
            result = radial_derivatives_liquid_general(x, y, K, rho, g, orbital_freq)

        return result

    if R_core is None:
        mantle_radii = radius_array >= 0.
        core_radii = None
        R_core_index = None
    else:
        mantle_radii = radius_array > R_core
        core_radii = radius_array <= R_core
        R_core_index = find_nearest(radius_array, R_core) - 1

    # Get initial solutions
    initial_solutions_solid = \
        solid_solutions(radius_array[mantle_radii],
                        complex_shear_array[mantle_radii],
                        bulk_array[mantle_radii],
                        density_array[mantle_radii], orbital_freq, order_l=2)

    initial_solutions_liquid = None
    if R_core is not None:
        if static_core:
            initial_solutions_liquid = \
                liquid_solutions_static(
                        radius_array[core_radii],
                        bulk_array[core_radii],
                        density_array[core_radii], order_l=2)
            # make a list to match the dynamic solution type
            initial_solutions_liquid = (initial_solutions_liquid,)
        else:
            initial_solutions_liquid = \
                liquid_solutions(radius_array[core_radii],
                                 bulk_array[core_radii],
                                 density_array[core_radii], orbital_freq, order_l=2)

    # IVP Approach
    from scipy.integrate import solve_ivp

    span_core = None
    if R_core is not None:
        span_core = (radius_array[0], R_core)
    span_mantle = (radius_array[mantle_radii][0], R_planet)

    solution_ys_core = None
    solution_rs_core = None
    if R_core is not None:
        # Solve core first
        solution_ys_core = list()
        solution_rs_core = list()
        if static_core:
            solution_range = range(1)
        else:
            solution_range = range(2)

        for sn in solution_range:
            print(f'Solving Core Solution {sn+1} of {len(solution_range)}')
            initial_value = initial_solutions_liquid[sn][:, 0]

            # Apply bottom boundary conditions
            # Roberts & Nimmo (2008): Liquid innermost zone.
            # initial_value[1] = 1.
            # initial_value[3] = 1.

            # Henning and Hurford (2014) Solid innermost
            # initial_value[0] = 1.
            # initial_value[1] = 1.

            # # BC for Tobie 2005 benchmark (From IcyDwarf)
            # initial_value[0] = 0.05
            # initial_value[1] = 0.01
            # initial_value[3] = -1.

            # # Bottom BC Based on Tobie 2005
            # initial_value[0] = 0.0
            # initial_value[1] = 0.0
            # initial_value[2] = 0.0

            # Use the initial guess function.
            pass

            # Make sure the initial value remains complex.
            initial_value = np.asarray(initial_value, dtype=np.complex128)
            solution = solve_ivp(deriv_liq, span_core, initial_value,
                                 method=integration_method, rtol=integration_tol,
                                 t_eval=radius_array[radius_array <= R_core])
            if solution.status == 0:
                print(f'Solution {sn} solved.')
                solution_rs_core.append(solution.t)
                solution_ys_core.append(solution.y)
            else:
                print(solution.message)

    # Solve mantle
    solution_ys_mantle = list()
    solution_rs_mantle = list()
    for sn in range(3):
        print(f'Solving Mantle Solution {sn+1}')
        initial_value = initial_solutions_solid[sn][:, 0]

        # Apply bottom BC (is this right??!)
        if R_core is None:
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

        else:
            # Boundary conditions are based on the core results.
            if static_core:
                # From Saito74
                if sn == 0:
                    initial_value[0] = 0.
                    initial_value[1] = -density_array[R_core_index] * solution_ys_core[0][0, -1]
                    initial_value[4] = solution_ys_core[0][0, -1]
                    initial_value[5] = solution_ys_core[0][1, -1] + \
                        (4. * np.pi * G * density_array[R_core_index] / gravity_array[R_core_index]) * \
                        solution_ys_core[0][0, -1]

                    initial_value[2] = 0.
                    initial_value[3] = 0.

                elif sn == 1:
                    initial_value[0] = 1.
                    initial_value[1] = density_array[R_core_index] * gravity_array[R_core_index] * initial_value[0]
                    initial_value[5] = -4. * np.pi * G * density_array[R_core_index] * initial_value[0]

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
                    initial_value[0] = solution_ys_core[0][0, -1]
                    initial_value[1] = solution_ys_core[0][1, -1]
                    initial_value[4] = solution_ys_core[0][2, -1]
                    initial_value[5] = solution_ys_core[0][3, -1]

                    initial_value[2] = 0.
                    initial_value[3] = 0.

                elif sn == 1:
                    initial_value[0] = solution_ys_core[1][0, -1]
                    initial_value[1] = solution_ys_core[1][1, -1]
                    initial_value[4] = solution_ys_core[1][2, -1]
                    initial_value[5] = solution_ys_core[1][3, -1]

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
                             t_eval=radius_array[mantle_radii])
        if solution.status == 0:
            print(f'Solution {sn} solved.')
            solution_rs_mantle.append(solution.t)
            solution_ys_mantle.append(solution.y)
        else:
            print(solution.message)

    # Build solution matrix at surface
    sol_surf_mtx = np.asarray([
        [solution_ys_mantle[0][1, -1], solution_ys_mantle[1][1, -1], solution_ys_mantle[2][1, -1]],
        [solution_ys_mantle[0][3, -1], solution_ys_mantle[1][3, -1], solution_ys_mantle[2][3, -1]],
        [solution_ys_mantle[0][5, -1], solution_ys_mantle[1][5, -1], solution_ys_mantle[2][5, -1]]
    ])
    sol_surf_mtx_inv = np.linalg.inv(sol_surf_mtx)
    surf_bc = np.zeros(3, dtype=np.complex128)
    surf_bc[2] = 5. / R_planet
    Q_vector = sol_surf_mtx_inv @ surf_bc

    # Combine Solutions
    tidal_y_sol_core = None
    y3_core = None
    core_ys = None
    tidal_y_sol_core_full = None
    if R_core is not None:
        if static_core:
            tidal_y_sol_core = Q_vector[0] * solution_ys_core[0]
            core_ys = (
                np.nan * np.ones_like(tidal_y_sol_core[0, :]),
                np.nan * np.ones_like(tidal_y_sol_core[0, :]),
                np.nan * np.ones_like(tidal_y_sol_core[0, :]),
                np.nan * np.ones_like(tidal_y_sol_core[0, :]),
                tidal_y_sol_core[0, :],
                np.nan * np.ones_like(tidal_y_sol_core[0, :])
            )

        else:
            tidal_y_sol_core = Q_vector[0] * solution_ys_core[0] + Q_vector[1] * solution_ys_core[1]
            y3_core = \
                (1. / (orbital_freq**2 * density_array[core_radii] * radius_array[core_radii])) * \
                    (density_array[core_radii] * gravity_array[core_radii] * tidal_y_sol_core[0, :] -
                     tidal_y_sol_core[1, :] - density_array[core_radii] * tidal_y_sol_core[2, :])
            core_ys = (
                tidal_y_sol_core[0, :],
                tidal_y_sol_core[1, :],
                y3_core,
                np.zeros_like(y3_core),
                tidal_y_sol_core[2, :],
                tidal_y_sol_core[3, :]
            )
        tidal_y_sol_core_full = np.vstack(core_ys)

    tidal_y_sol_mantle = Q_vector[0] * solution_ys_mantle[0] + Q_vector[1] * solution_ys_mantle[1] + Q_vector[2] * solution_ys_mantle[2]

    if R_core is None:
        tidal_y = tidal_y_sol_mantle
    else:
        tidal_y = np.concatenate((tidal_y_sol_core_full, tidal_y_sol_mantle), axis=1)

    model_radii.append(radius_array)
    model_names.append(model_name)
    model_results.append(tidal_y)

yplot(model_results, model_radii, labels=model_names,
      use_tobie_limits=True, plot_tobie=plot_tobie, plot_roberts=plot_roberts)

# BVP Approach
# from scipy.integrate import solve_bvp
# surf_condition = (5. / R_planet)
# def bc_residual(ya, yb):
#     return (ya[0], yb[1], ya[2], yb[3], ya[4], yb[5] - surf_condition)
#
# solution_ys = list()
# solution_rs = list()
# for sn in range(3):
#
#     initial_guess = initial_solutions[sn]
#     solution = solve_bvp(deriv_sol, bc_residual, radius_array[radius_array > R_core], initial_guess, tol=1e-3, verbose=2)
#
#     print(f'Solving Solution {sn}')
#     #     solution = solve_ivp(deriv_sol, span, initial_solutions[sn][:, 0])
#     if solution.status == 0:
#         print(f'Solution {sn} solved.')
#         solution_rs.append(solution.x)
#         solution_ys.append(solution.y)
#     else
#         print(solution.message)
#
# solution_ys.append(tidal_y_ivp)
# solution_rs.append(tidal_y_ivp_r)
#
# if len(solution_ys) != 0:
#     yplot(solution_ys, solution_rs)