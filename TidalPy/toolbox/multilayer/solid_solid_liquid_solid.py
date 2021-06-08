from typing import Tuple

import numpy as np

from .odes import dynamic_solid_ode, static_solid_ode, dynamic_liquid_ode, static_liquid_ode
from ...constants import G
from ...exceptions import IntegrationFailed, AttributeNotSetError
from ...tides.multilayer.nondimensional import non_dimensionalize_physicals, re_dimensionalize_radial_func
from ...tides.multilayer.numerical_int import find_initial_guess
from ...tides.multilayer.numerical_int.interfaces import find_interface_func
from ...utilities.integration.integrate import rk_integrator
from ...utilities.performance import njit

TidalYSolType = Tuple[Tuple[np.ndarray, np.ndarray, np.ndarray],
                      Tuple[np.ndarray, ...],
                      Tuple[np.ndarray, np.ndarray, np.ndarray]]


@njit(cacheable=True)
def convergence_ssls_static_liq(tidal_y_solutions_by_layer: TidalYSolType, surface_solution: np.ndarray,
                                gravity_array_layer2: np.ndarray, density_array_layer2: np.ndarray) -> np.ndarray:
    """ Determine the radial solution convergence for a planet with solid-solid-liquid-solid structure.
    A static liquid layer is assumed.

    Used in the numerical shooting method.

    Parameters
    ----------
    tidal_y_solutions_by_layer : TidalYSolType
        Radial functions solved via integration separated by layer.
    surface_solution : np.ndarray
        Surface boundary condition used to find the constants at the top-most layer.
    gravity_array_layer2 : np.ndarray
        Acceleration due to gravity within the liquid layer [m s-2]
    density_array_layer2 : np.ndarray
        Density within the liquid layer [kg m-3]

    Returns
    -------
    tidal_y : np.ndarray
        Collapsed radial functions for the entire planet. Scaled by the correct constants.
        This will be a 6 x N ndarray for the six radial functions.

    """

    # Pull out data
    tidal_y_layer0 = tidal_y_solutions_by_layer[0]
    tidal_y_layer1 = tidal_y_solutions_by_layer[1]
    tidal_y_layer2 = tidal_y_solutions_by_layer[2]
    tidal_y_layer3 = tidal_y_solutions_by_layer[3]
    liquid_density_interface_2 = density_array_layer2[0]
    gravity_interface_2 = gravity_array_layer2[0]

    # Build solution matrix at surface
    sol_surf_mtx = np.asarray([
        [tidal_y_layer3[0][1, -1], tidal_y_layer3[1][1, -1], tidal_y_layer3[2][1, -1]],
        [tidal_y_layer3[0][3, -1], tidal_y_layer3[1][3, -1], tidal_y_layer3[2][3, -1]],
        [tidal_y_layer3[0][5, -1], tidal_y_layer3[1][5, -1], tidal_y_layer3[2][5, -1]]
    ])
    sol_surf_mtx_inv = np.linalg.inv(sol_surf_mtx)
    C_layer3_vector = sol_surf_mtx_inv @ surface_solution

    # Solve for the liquid layer's constants
    C_layer2_vector = np.zeros(1, dtype=np.complex128)
    C_layer2_vector[0] = C_layer3_vector[0]

    # Solve for layer 1's constants
    C_layer1_vector = np.zeros(3, dtype=np.complex128)

    y4_frac_1 = tidal_y_layer1[0][3, -1] / tidal_y_layer1[2][3, -1]
    y4_frac_2 = tidal_y_layer1[1][3, -1] / tidal_y_layer1[2][3, -1]

    gamma_1 = (tidal_y_layer1[0][1, -1] - y4_frac_1 * tidal_y_layer1[2][1, -1]) - \
              liquid_density_interface_2 * (gravity_interface_2 * (
            tidal_y_layer1[0][0, -1] - y4_frac_1 * tidal_y_layer1[2][0, -1]) -
                                            (tidal_y_layer1[0][4, -1] - y4_frac_1 * tidal_y_layer1[2][4, -1]))
    gamma_2 = (tidal_y_layer1[1][1, -1] - y4_frac_2 * tidal_y_layer1[2][1, -1]) - \
              liquid_density_interface_2 * (gravity_interface_2 * (
            tidal_y_layer1[1][0, -1] - y4_frac_2 * tidal_y_layer1[2][0, -1]) -
                                            (tidal_y_layer1[1][4, -1] - y4_frac_2 * tidal_y_layer1[2][4, -1]))

    C_layer1_vector[0] = C_layer2_vector[0]
    C_layer1_vector[1] = (-gamma_1 / gamma_2) * C_layer1_vector[0]
    C_layer1_vector[2] = -y4_frac_1 * C_layer1_vector[0] - y4_frac_2 * C_layer1_vector[1]

    # Solve for the innermost layers constants
    C_layer0_vector = np.zeros(3, dtype=np.complex128)
    C_layer0_vector[0] = C_layer1_vector[0]
    C_layer0_vector[1] = C_layer1_vector[1]
    C_layer0_vector[2] = C_layer1_vector[1]

    # Solve for the liquid layer's y's
    tidal_y_layer2 = C_layer2_vector[0] * tidal_y_layer2[0]

    layer2_ys = (
        np.nan * np.ones_like(tidal_y_layer2[0, :]),
        np.nan * np.ones_like(tidal_y_layer2[0, :]),
        np.nan * np.ones_like(tidal_y_layer2[0, :]),
        np.nan * np.ones_like(tidal_y_layer2[0, :]),
        tidal_y_layer2[0, :],
        np.nan * np.ones_like(tidal_y_layer2[0, :])
    )

    tidal_y_layer2_full = np.vstack(layer2_ys)

    # Solve for total planet y's
    tidal_y_layer0 = C_layer0_vector[0] * tidal_y_layer0[0] + C_layer0_vector[1] * tidal_y_layer0[1] + \
                     C_layer0_vector[2] * tidal_y_layer0[2]

    tidal_y_layer1 = C_layer1_vector[0] * tidal_y_layer1[0] + C_layer1_vector[1] * tidal_y_layer1[1] + \
                     C_layer1_vector[2] * tidal_y_layer1[2]

    tidal_y_layer3 = C_layer3_vector[0] * tidal_y_layer3[0] + C_layer3_vector[1] * tidal_y_layer3[1] + \
                     C_layer3_vector[2] * tidal_y_layer3[2]

    # Combine solutions for all layers
    tidal_y = np.concatenate((tidal_y_layer0, tidal_y_layer1, tidal_y_layer2_full, tidal_y_layer3), axis=1)

    return tidal_y


@njit(cacheable=True)
def convergence_ssls_dynamic_liq(tidal_y_solutions_by_layer: TidalYSolType, surface_solution: np.ndarray,
                                 gravity_array_layer2: np.ndarray, density_array_layer2: np.ndarray,
                                 radius_array_layer2: np.ndarray, orbital_freq: float) -> np.ndarray:
    """ Determine the radial solution convergence for a planet with solid-solid-liquid-solid structure.
    A dynamic liquid layer is assumed.

    Used in the numerical shooting method.

    Parameters
    ----------
    tidal_y_solutions_by_layer : TidalYSolType
        Radial functions solved via integration separated by layer.
    surface_solution : np.ndarray
        Surface boundary condition used to find the constants at the top-most layer.
    gravity_array_layer2 : np.ndarray
        Acceleration due to gravity within the liquid layer [m s-2]
    density_array_layer2 : np.ndarray
        Density within the liquid layer [kg m-3]
    radius_array_layer2 : np.ndarray
        Radius array for the liquid layer [m]
    orbital_freq : float
        Forcing frequency [rad s-1]

    Returns
    -------
    tidal_y : np.ndarray
        Collapsed radial functions for the entire planet. Scaled by the correct constants.
        This will be a 6 x N ndarray for the six radial functions.

    """

    # Pull out data
    tidal_y_layer0 = tidal_y_solutions_by_layer[0]
    tidal_y_layer1 = tidal_y_solutions_by_layer[1]
    tidal_y_layer2 = tidal_y_solutions_by_layer[2]
    tidal_y_layer3 = tidal_y_solutions_by_layer[3]

    # Build solution matrix at surface
    sol_surf_mtx = np.asarray([
        [tidal_y_layer3[0][1, -1], tidal_y_layer3[1][1, -1], tidal_y_layer3[2][1, -1]],
        [tidal_y_layer3[0][3, -1], tidal_y_layer3[1][3, -1], tidal_y_layer3[2][3, -1]],
        [tidal_y_layer3[0][5, -1], tidal_y_layer3[1][5, -1], tidal_y_layer3[2][5, -1]]
    ])
    sol_surf_mtx_inv = np.linalg.inv(sol_surf_mtx)
    C_layer3_vector = sol_surf_mtx_inv @ surface_solution

    # Solve for the liquid layer constants
    C_layer2_vector = np.zeros(2, dtype=np.complex128)
    C_layer2_vector[0] = C_layer3_vector[0]
    C_layer2_vector[1] = C_layer3_vector[1]

    # Solve for the layer 1 constants
    C_layer1_vector = np.zeros(3, dtype=np.complex128)
    C_layer1_vector[0] = C_layer2_vector[0]
    C_layer1_vector[1] = C_layer2_vector[1]
    y4_frac_1_layer_1 = tidal_y_layer1[0][3, -1] / tidal_y_layer1[2][3, -1]
    y4_frac_2_layer_1 = tidal_y_layer1[1][3, -1] / tidal_y_layer1[2][3, -1]
    C_layer1_vector[2] = -y4_frac_1_layer_1 * C_layer1_vector[0] - y4_frac_2_layer_1 * C_layer1_vector[1]

    # Solve for the innermost layers constants
    C_layer0_vector = np.zeros(3, dtype=np.complex128)
    C_layer0_vector[0] = C_layer1_vector[0]
    C_layer0_vector[1] = C_layer1_vector[1]
    C_layer0_vector[2] = C_layer1_vector[1]

    # Solve for the liquid layer's y's
    tidal_y_layer2 = C_layer2_vector[0] * tidal_y_layer2[0] + C_layer2_vector[1] * tidal_y_layer2[1]

    # Liquid layer is missing two y's, fix that now.
    y3_layer2 = \
        (1. / (orbital_freq**2 * density_array_layer2 * radius_array_layer2)) * \
        (density_array_layer2 * gravity_array_layer2 * tidal_y_layer2[0, :] -
         tidal_y_layer2[1, :] - density_array_layer2 * tidal_y_layer2[2, :])

    layer2_ys = (
        tidal_y_layer2[0, :],
        tidal_y_layer2[1, :],
        y3_layer2,
        np.nan * np.ones_like(tidal_y_layer2[0, :]),
        tidal_y_layer2[2, :],
        tidal_y_layer2[3, :]
    )

    tidal_y_layer2_full = np.vstack(layer2_ys)

    # Solve for total planet y's
    tidal_y_layer0 = C_layer0_vector[0] * tidal_y_layer0[0] + C_layer0_vector[1] * tidal_y_layer0[1] + \
                     C_layer0_vector[2] * tidal_y_layer0[2]

    tidal_y_layer1 = C_layer1_vector[0] * tidal_y_layer1[0] + C_layer1_vector[1] * tidal_y_layer1[1] + \
                     C_layer1_vector[2] * tidal_y_layer1[2]

    tidal_y_layer3 = C_layer3_vector[0] * tidal_y_layer3[0] + C_layer3_vector[1] * tidal_y_layer3[1] + \
                     C_layer3_vector[2] * tidal_y_layer3[2]

    # Combine solutions for all layers
    tidal_y = np.concatenate((tidal_y_layer0, tidal_y_layer1, tidal_y_layer2_full, tidal_y_layer3), axis=1)

    return tidal_y


def calculate_ssls(radius: np.ndarray, shear_modulus: np.ndarray, bulk_modulus: np.ndarray,
                   density: np.ndarray, gravity: np.ndarray, frequency: float,
                   interface_1_radius: float, interface_2_radius: float, interface_3_radius: float,
                   layer_0_static: bool = False, layer_1_static: bool = False,
                   layer_2_static: bool = True, layer_3_static: bool = False,
                   surface_boundary_condition: np.ndarray = None,
                   order_l: int = 2, use_kamata: bool = True,
                   use_julia: bool = False,
                   use_numba_integrator: bool = False,
                   verbose: bool = False,
                   int_rtol: float = 1.0e-6, int_atol: float = 1.0e-6,
                   scipy_int_method: str = 'RK45', julia_int_method: str = 'Tsit5',
                   non_dimensionalize: bool = False,
                   planet_bulk_density: float = None) -> Tuple[np.ndarray, np.ndarray]:
    """ Calculate the radial solution for a planet that has a three layer structure: Solid-Liquid-Solid.

    Parameters
    ----------
    radius : np.ndarray
        Full planet radius array [m]
    shear_modulus : np.ndarray
        Full planet shear modulus (can be complex) at each `radius` [Pa]
    bulk_modulus : np.ndarray
        Full planet bulk modulus (can be complex) at each `radius` [Pa]
    density : np.ndarray
        Full planet density at each `radius` [kg m-3]
    gravity : np.ndarray
        Full planet gravity at each `radius` [m s-2]
    frequency : float
        Forcing frequency [rad s-1]
    interface_1_radius : float
        Radius of the first (Solid-Solid) radius [m]
    interface_2_radius : float
        Radius of the first (Solid-Liquid) radius [m]
    interface_3_radius : float
        Radius of the first (Liquid-Solid) radius [m]
    layer_0_static : bool = False
        If True, layer 0 will be treated under the static tidal assumption (w=0).
    layer_1_static : bool = False
        If True, layer 1 will be treated under the static tidal assumption (w=0).
    layer_2_static : bool = True
        If True, layer 2 will be treated under the static tidal assumption (w=0).
    layer_3_static : bool = False
        If True, layer 3 will be treated under the static tidal assumption (w=0).
    surface_boundary_condition: np.ndarray = None
        The surface boundary condition, for tidal solutions y2, y4, y6, = (0, 0, (2l+1)/R)
            Tidal solution will be the default if `None` is provided.
    order_l : int = 2
        Tidal harmonic order.
    use_kamata : bool = True
        If True, the Kamata+2015 initial conditions will be used at the base of layer 0.
        Otherwise, the Takeuchi & Saito 1972 initial conditions will be used.
    use_julia : bool = False
        If True, the Julia `diffeqpy` integration tools will be used.
        Otherwise, `scipy.integrate.solve_ivp` or TidalPy's numba-safe integrator will be used.
    use_numba_integrator : bool = False
        If True, TidalPy's numba-safe RK-based integrator will be used.
        Otherwise, `scipy.integrate.solve_ivp` or Julia `diffeqpy` integrator will be used.
    verbose : bool = False
        If True, the function will print some information to console during calculation (may cause a slow down).
    int_rtol : float = 1.0e-6
        Integration relative error.
    int_atol : float = 1.0e-4
        Integration absolute error.
    scipy_int_method : str = 'RK45'
        Integration method for the Scipy integration scheme.
        See options here (note some do not work for complex numbers):
            https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
    julia_int_method : str = 'Tsit5'
        Integration method for the Julia integration scheme.
        See options here (note some do not work for complex numbers):
            `TidalPy.utilities.julia_helper.integration_methods.py`
    non_dimensionalize : bool = False
        If True, integration will use dimensionless variables. These will be converted back before output is given to
        the user.
    planet_bulk_density : float = None
        Must be provided if non_dimensionalize is True. Bulk density of the planet.

    Returns
    -------
    tidal_y : np.ndarray
        The radial solution throughout the entire planet.

    """
    # Figure out the index of the boundaries and geometry of the problem
    layer_0_indx = radius <= interface_1_radius
    layer_1_indx = np.logical_and(radius > interface_1_radius, radius <= interface_2_radius)
    layer_2_indx = np.logical_and(radius > interface_2_radius, radius <= interface_3_radius)
    layer_3_indx = radius > interface_3_radius

    # Non-dimensionalize inputs
    planet_radius = radius[-1]
    if non_dimensionalize:
        if planet_bulk_density is None:
            raise AttributeNotSetError('Planet bulk modulus must be provided if non-dimensionalize is True.')

        radius, gravity, density, shear_modulus, bulk_modulus, frequency, G_to_use = \
            non_dimensionalize_physicals(radius, gravity, density, shear_modulus, bulk_modulus, frequency,
                                         mean_radius=planet_radius, bulk_density=planet_bulk_density)
    else:
        G_to_use = G

    radial_span_0 = (radius[layer_0_indx][0], radius[layer_0_indx][-1])
    radial_span_1 = (radius[layer_1_indx][0], radius[layer_1_indx][-1])
    radial_span_2 = (radius[layer_2_indx][0], radius[layer_2_indx][-1])
    radial_span_3 = (radius[layer_3_indx][0], radius[layer_3_indx][-1])

    radial_solve_0 = radius[layer_0_indx]
    radial_solve_1 = radius[layer_1_indx]
    radial_solve_2 = radius[layer_2_indx]
    radial_solve_3 = radius[layer_3_indx]

    # Find solution at the top of the planet -- this is dependent on the forcing type.
    #     Tides (default here) follow the (y2, y4, y6) = (0, 0, (2l+1)/R) rule
    if surface_boundary_condition is None:
        surface_boundary_condition = np.zeros(3, dtype=np.complex128)
        surface_boundary_condition[2] = (2. * order_l + 1.) / radius[-1]

    # Initial (base) guess will be for a solid layer
    initial_value_func = find_initial_guess(is_kamata=use_kamata, is_solid=True, is_dynamic=(not layer_0_static))
    if layer_0_static:
        initial_value_tuple = initial_value_func(radius[0], shear_modulus[0], bulk_modulus[0], density[0],
                                                 order_l=order_l, G_to_use=G_to_use)
    else:
        initial_value_tuple = initial_value_func(radius[0], shear_modulus[0], bulk_modulus[0], density[0],
                                                 frequency, order_l=order_l, G_to_use=G_to_use)

    # Find the differential equation
    if layer_0_static:
        radial_derivative_layer_0 = static_solid_ode
        derivative_inputs_layer_0 = (radius[layer_0_indx], shear_modulus[layer_0_indx], bulk_modulus[layer_0_indx],
                                     density[layer_0_indx], gravity[layer_0_indx],
                                     order_l, G_to_use)
    else:
        radial_derivative_layer_0 = dynamic_solid_ode
        derivative_inputs_layer_0 = (radius[layer_0_indx], shear_modulus[layer_0_indx], bulk_modulus[layer_0_indx],
                                     density[layer_0_indx], gravity[layer_0_indx], frequency,
                                     order_l, G_to_use)

    if layer_1_static:
        radial_derivative_layer_1 = static_solid_ode
        derivative_inputs_layer_1 = (radius[layer_1_indx], shear_modulus[layer_1_indx], bulk_modulus[layer_1_indx],
                                     density[layer_1_indx], gravity[layer_1_indx],
                                     order_l, G_to_use)
    else:
        radial_derivative_layer_1 = dynamic_solid_ode
        derivative_inputs_layer_1 = (radius[layer_1_indx], shear_modulus[layer_1_indx], bulk_modulus[layer_1_indx],
                                     density[layer_1_indx], gravity[layer_1_indx], frequency,
                                     order_l, G_to_use)

    if layer_2_static:
        radial_derivative_layer_2 = static_liquid_ode
        derivative_inputs_layer_2 = (radius[layer_2_indx],
                                     density[layer_2_indx], gravity[layer_2_indx],
                                     order_l, G_to_use)
    else:
        radial_derivative_layer_2 = dynamic_liquid_ode
        derivative_inputs_layer_2 = (radius[layer_2_indx], bulk_modulus[layer_2_indx],
                                     density[layer_2_indx], gravity[layer_2_indx], frequency,
                                     order_l, G_to_use)
    if layer_3_static:
        radial_derivative_layer_3 = static_solid_ode
        derivative_inputs_layer_3 = (radius[layer_3_indx], shear_modulus[layer_3_indx], bulk_modulus[layer_3_indx],
                                     density[layer_3_indx], gravity[layer_3_indx],
                                     order_l, G_to_use)
    else:
        radial_derivative_layer_3 = dynamic_solid_ode
        derivative_inputs_layer_3 = (radius[layer_3_indx], shear_modulus[layer_3_indx], bulk_modulus[layer_3_indx],
                                     density[layer_3_indx], gravity[layer_3_indx], frequency,
                                     order_l, G_to_use)

    # Find interfaces
    interface_1_func = find_interface_func(lower_layer_is_solid=True, lower_layer_is_static=layer_0_static,
                                           upper_layer_is_solid=True, upper_layer_is_static=layer_1_static)

    interface_2_func = find_interface_func(lower_layer_is_solid=True, lower_layer_is_static=layer_1_static,
                                           upper_layer_is_solid=False, upper_layer_is_static=layer_2_static,
                                           liquid_density=density[layer_2_indx][-1],
                                           interface_gravity=gravity[layer_2_indx][-1])

    interface_3_func = find_interface_func(lower_layer_is_solid=False, lower_layer_is_static=layer_2_static,
                                           upper_layer_is_solid=True, upper_layer_is_static=layer_3_static,
                                           liquid_density=density[layer_2_indx][-1],
                                           interface_gravity=gravity[layer_2_indx][-1])

    solutions_by_layer = [list(), list(), list(), list()]
    for layer_i in range(4):
        # Solve the inner-most layer first using the initial conditions we found above
        if layer_i == 0:
            radial_span = radial_span_0
            radial_solve = radial_solve_0
            derivatives = radial_derivative_layer_0
            initial_values_to_use = initial_value_tuple
            derivative_inputs = derivative_inputs_layer_0
        elif layer_i == 1:
            radial_span = radial_span_1
            radial_solve = radial_solve_1
            derivatives = radial_derivative_layer_1
            derivative_inputs = derivative_inputs_layer_1
            # Initial values are based on the previous layer's results
            initial_values_to_use = interface_1_func(solutions_by_layer[0])
        elif layer_i == 2:
            radial_span = radial_span_2
            radial_solve = radial_solve_2
            derivatives = radial_derivative_layer_2
            derivative_inputs = derivative_inputs_layer_2
            # Initial values are based on the previous layer's results
            initial_values_to_use = interface_2_func(solutions_by_layer[1])

            if layer_2_static:
                # There is only one solution provided by the interface func. Make it a list so the loop later on
                #   works.
                initial_values_to_use = (initial_values_to_use,)
        else:
            radial_span = radial_span_3
            radial_solve = radial_solve_3
            derivatives = radial_derivative_layer_3
            derivative_inputs = derivative_inputs_layer_3
            # Initial values are based on the previous layer's results
            if layer_2_static:
                # The interface function will only expect one solution, but they are stored in a tuple no matter what.
                #   pull that out first.
                liq_sols = solutions_by_layer[2][0]
            else:
                liq_sols = solutions_by_layer[2]

            initial_values_to_use = interface_3_func(liq_sols)

        if use_julia:
            def diffeq_julia(u, p, r):
                output = derivatives(r, u, *p)
                return list(output)

            # Import Julia's Diffeqpy and reinit the problem
            from ...utilities.julia_helper.integration_methods import get_julia_solver
            ode, solver = get_julia_solver(julia_int_method)

            if verbose:
                print(f'Solving Layer {layer_i + 1} (with SciPy, using {scipy_int_method})...')

            for solution_num, initial_values in enumerate(initial_values_to_use):
                problem = ode.ODEProblem(diffeq_julia, initial_values, radial_span, derivative_inputs)
                solution = ode.solve(problem, solver(), abstol=int_atol, reltol=int_rtol)

                # Julia does not have the same t_eval. There is the "saveat" keyword but can cause issues.
                #    So perform an interpolation for the desired radii
                u_T = np.transpose(solution.u)
                u = np.zeros((u_T.shape[0], radial_solve.size), dtype=np.complex128)
                for i in range(u_T.shape[0]):
                    u[i, :] = np.interp(radial_solve, solution.t, u_T[i, :])
                solutions_by_layer[layer_i].append(u)

            if verbose:
                print('\nIntegration Done!')

        elif use_numba_integrator:

            if scipy_int_method == 'RK23':
                rk_method = 0
            elif scipy_int_method == 'RK45':
                rk_method = 1
            else:
                raise NotImplementedError

            if verbose:
                print(f"Solving Layer {layer_i + 1} (with TidalPy's Numba integrator, using {scipy_int_method})...")

            for solution_num, initial_values in enumerate(initial_values_to_use):

                # TODO: The copy below makes this work for numba on multiprocessing due to contiguous arrays... wtf!
                ts, ys, status, message, success = \
                    rk_integrator(derivatives, radial_span, np.copy(initial_values),
                                  args=derivative_inputs,
                                  rk_method=rk_method,
                                  t_eval_N=radial_solve.size, t_eval_log=False, use_teval=True,
                                  rtol=int_rtol, atol=int_atol, verbose=False)

                if status != 0:
                    raise IntegrationFailed(f'Integration Solution Failed for {layer_i} at solution #{solution_num}.'
                                            f'\n\t{message}')

                solutions_by_layer[layer_i].append(ys)

            if verbose:
                print('\nIntegration Done!')

        else:
            from scipy.integrate import solve_ivp

            if verbose:
                print(f'Solving Layer {layer_i + 1} (with SciPy, using {scipy_int_method})...')

            for solution_num, initial_values in enumerate(initial_values_to_use):
                solution = solve_ivp(
                        derivatives, radial_span, initial_values, t_eval=radial_solve, args=derivative_inputs,
                        method=scipy_int_method, vectorized=False, rtol=int_rtol, atol=int_atol)

                if solution.status != 0:
                    raise IntegrationFailed(f'Integration Solution Failed for {layer_i} at solution #{solution_num}.'
                                            f'\n\t{solution.message}')

                solutions_by_layer[layer_i].append(solution.y)

            if verbose:
                print('\nIntegration Done!')

        # Done with layer
        solutions_by_layer[layer_i] = tuple(solutions_by_layer[layer_i])

    solutions_by_layer = tuple(solutions_by_layer)

    # Find the convergence based on the surface boundary conditions
    if verbose:
        print('Solving convergence...')

    if layer_2_static:
        tidal_y = \
            convergence_ssls_static_liq(solutions_by_layer, surface_boundary_condition,
                                        gravity[layer_2_indx], density[layer_2_indx])
    else:
        tidal_y = \
            convergence_ssls_dynamic_liq(solutions_by_layer, surface_boundary_condition,
                                         gravity[layer_2_indx], density[layer_2_indx],
                                         radius[layer_2_indx], frequency)

    if verbose:
        print('Done!')

    if non_dimensionalize:
        if verbose:
            print('Redimensionalizing Radial Functions.')
        tidal_y = re_dimensionalize_radial_func(tidal_y, planet_radius, planet_bulk_density)

    # Now that tidal_y has been found, we can find the radial derivatives which are used in some calculations.
    tidal_y_derivative = np.zeros_like(tidal_y)
    tidal_y_derivative[:, layer_0_indx] = np.stack(
            radial_derivative_layer_0(radius[layer_0_indx], tidal_y[:, layer_0_indx], *derivative_inputs_layer_0)
    )
    tidal_y_derivative[:, layer_1_indx] = np.stack(
            radial_derivative_layer_1(radius[layer_1_indx], tidal_y[:, layer_1_indx], *derivative_inputs_layer_1)
    )
    tidal_y_derivative[:, layer_3_indx] = np.stack(
            radial_derivative_layer_3(radius[layer_3_indx], tidal_y[:, layer_3_indx], *derivative_inputs_layer_3)
    )
    # Layer 2 is liquid so that complicates the calculation slightly
    if layer_2_static:
        # TODO: This is not correct, we could pull out dy_5/dr but it is not used in subsequent calculations,
        #  so let's just leave all derivatives as zero (which is done in the initialization above).
        pass
    else:
        tidal_y_liq = np.zeros((4, len(radius[layer_2_indx])), dtype=np.complex128)
        tidal_y_liq[0, :] = tidal_y[0, layer_2_indx]
        tidal_y_liq[1, :] = tidal_y[1, layer_2_indx]
        tidal_y_liq[2, :] = tidal_y[4, layer_2_indx]
        tidal_y_liq[3, :] = tidal_y[5, layer_2_indx]
        liq_derivatives = np.stack(
                radial_derivative_layer_2(radius[layer_2_indx], tidal_y_liq, *derivative_inputs_layer_2)
        )
        tidal_y_derivative[0, layer_2_indx] = liq_derivatives[0, :]
        tidal_y_derivative[1, layer_2_indx] = liq_derivatives[1, :]
        tidal_y_derivative[4, layer_2_indx] = liq_derivatives[2, :]
        tidal_y_derivative[5, layer_2_indx] = liq_derivatives[3, :]

    return tidal_y, tidal_y_derivative
