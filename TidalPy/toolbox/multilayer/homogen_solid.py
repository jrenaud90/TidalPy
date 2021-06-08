from typing import Tuple

import numpy as np

from .odes import static_solid_ode, dynamic_solid_ode
from ...constants import G
from ...exceptions import IntegrationFailed, AttributeNotSetError
from ...tides.multilayer.nondimensional import non_dimensionalize_physicals, re_dimensionalize_radial_func
from ...tides.multilayer.numerical_int import find_initial_guess
from ...utilities.integration.integrate import rk_integrator
from ...utilities.performance import njit

MAX_DATA_SIZE = 2000


@njit(cacheable=True)
def convergence_solid(solid_solutions: Tuple[np.ndarray, np.ndarray, np.ndarray], surface_solution: np.ndarray):
    """ Determine the radial solution convergence for a homogeneous planet.

        Used in the numerical shooting method.

        Parameters
        ----------
        solid_solutions : Tuple[np.ndarray, np.ndarray, np.ndarray]
            Radial functions solved via integration for an entire homogeneous planet.
        surface_solution : np.ndarray
            Surface boundary condition used to find the constants at the top-most layer.

        Returns
        -------
        tidal_y : np.ndarray
            Collapsed radial functions for the entire planet. Scaled by the correct constants.
            This will be a 6 x N ndarray for the six radial functions.

        """

    # Build solution matrix at surface
    sol_surf_mtx = np.asarray([
        [solid_solutions[0][1, -1], solid_solutions[1][1, -1], solid_solutions[2][1, -1]],
        [solid_solutions[0][3, -1], solid_solutions[1][3, -1], solid_solutions[2][3, -1]],
        [solid_solutions[0][5, -1], solid_solutions[1][5, -1], solid_solutions[2][5, -1]]
    ])
    sol_surf_mtx_inv = np.linalg.inv(sol_surf_mtx)
    C_vector = sol_surf_mtx_inv @ surface_solution

    # Solve for total planet y's
    tidal_y = C_vector[0] * solid_solutions[0] + C_vector[1] * solid_solutions[1] + \
              C_vector[2] * solid_solutions[2]

    return tidal_y


def calculate_homogen_solid(radius: np.ndarray, shear_modulus: np.ndarray, bulk_modulus: np.ndarray,
                            density: np.ndarray, gravity: np.ndarray, frequency: float,
                            order_l: int = 2, use_static: bool = False,
                            surface_boundary_condition: np.ndarray = None,
                            use_kamata: bool = True, use_julia: bool = False,
                            use_numba_integrator: bool = False,
                            verbose: bool = False,
                            int_rtol: float = 1.0e-6, int_atol: float = 1.0e-4, scipy_int_method: str = 'RK45',
                            julia_int_method: str = 'Tsit5',
                            non_dimensionalize: bool = False,
                            planet_bulk_density: float = None) -> Tuple[np.ndarray, np.ndarray]:
    """ Calculate the radial solution for a homogeneous, solid planet.

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
    use_static : bool = False
        If True, the planet will be treated under the static tidal assumption (w=0).
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

    # Determine Geometry
    radial_span = (radius[0], radius[-1])

    # Find solution at the top of the planet -- this is dependent on the forcing type.
    #     Tides (default here) follow the (y2, y4, y6) = (0, 0, (2l+1)/R) rule
    if surface_boundary_condition is None:
        surface_boundary_condition = np.zeros(3, dtype=np.complex128)
        surface_boundary_condition[2] = (2. * order_l + 1.) / radius[-1]

    # Initial (base) guess will be for a solid layer
    is_dynamic = not use_static
    initial_value_func = find_initial_guess(is_kamata=use_kamata, is_solid=True, is_dynamic=is_dynamic)
    if use_static:
        initial_value_tuple = initial_value_func(radius[0], shear_modulus[0], bulk_modulus[0], density[0],
                                                 order_l=order_l, G_to_use=G_to_use)
    else:
        initial_value_tuple = initial_value_func(radius[0], shear_modulus[0], bulk_modulus[0], density[0],
                                                 frequency, order_l=order_l, G_to_use=G_to_use)

    # Find the differential equation
    if use_static:
        radial_derivative = static_solid_ode
        derivative_inputs = (radius, shear_modulus, bulk_modulus, density, gravity, order_l, G_to_use)
    else:
        radial_derivative = dynamic_solid_ode
        derivative_inputs = (radius, shear_modulus, bulk_modulus, density, gravity, frequency, order_l, G_to_use)

    solutions = list()
    for solution_num, initial_values in enumerate(initial_value_tuple):
        if use_julia:
            def diffeq_julia(u, p, r):
                output = radial_derivative(r, u, *p)

                return list(output)

            # Import Julia's Diffeqpy and reinit the problem
            from ...utilities.julia_helper.integration_methods import get_julia_solver
            ode, solver = get_julia_solver(julia_int_method)

            problem = ode.ODEProblem(diffeq_julia, initial_values, radial_span, derivative_inputs)
            if verbose:
                print(f'Solving (with Julia, using {julia_int_method})...')

            solution = ode.solve(problem, solver(), abstol=int_atol, reltol=int_rtol)

            # Julia does not have the same t_eval. There is the "saveat" keyword but can cause issues.
            #    So perform an interpolation for the desired radii
            u_T = np.transpose(solution.u)
            y = np.zeros((u_T.shape[0], radius.size), dtype=np.complex128)
            for i in range(u_T.shape[0]):
                y[i, :] = np.interp(radius, solution.t, u_T[i, :])

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
                print(f"Solving solution {solution_num} (with TidalPy's Numba integrator, using {scipy_int_method})...")

            ts, ys, status, message, success = \
                rk_integrator(radial_derivative, radial_span, initial_values,
                              args=derivative_inputs,
                              rk_method=rk_method,
                              t_eval_N=radius.size, t_eval_log=False, use_teval=True,
                              rtol=int_rtol, atol=int_atol, verbose=False)

            if status != 0:
                raise IntegrationFailed(f'Integration Solution Failed for solution {solution_num}.'
                                        f'\n\t{message}')
            y = ys

            if verbose:
                print('\nIntegration Done!')

        else:
            from scipy.integrate import solve_ivp

            if verbose:
                print(f'Solving (with SciPy, using {scipy_int_method})...')

            solution = solve_ivp(radial_derivative, radial_span, initial_values, t_eval=radius, args=derivative_inputs,
                                 method=scipy_int_method, vectorized=False, rtol=int_rtol, atol=int_atol)

            if solution.status != 0:
                raise IntegrationFailed(
                        f'Integration Solution Failed for homogeneous model at solution #{solution_num}.'
                        f'\n\t{solution.message}')

            y = solution.y

            if verbose:
                print('\nIntegration Done!')

        solutions.append(y)

    # Find the convergence based on the surface boundary conditions
    if verbose:
        print('Solving convergence...')

    tidal_y = convergence_solid(solutions, surface_boundary_condition)

    if verbose:
        print('Done!')

    if non_dimensionalize:
        if verbose:
            print('Redimensionalizing Radial Functions.')
        tidal_y = re_dimensionalize_radial_func(tidal_y, planet_radius, planet_bulk_density)

    # Now that tidal_y has been found, we can find the radial derivatives which are used in some calculations.
    tidal_y_derivative = np.stack(radial_derivative(radius, tidal_y, *derivative_inputs))

    return tidal_y, tidal_y_derivative
