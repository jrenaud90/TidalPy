from typing import List

import numpy as np

from .collapse import find_collapse_func
from .derivatives import known_multilayer_odes
from .initial_conditions import find_initial_guess
from .interfaces import find_interface_func
from ..nondimensional import non_dimensionalize_physicals, re_dimensionalize_radial_func
from ....constants import G
from ....exceptions import AttributeNotSetError, IntegrationFailed
from ....utilities.integration.rk_integrator import rk_integrate


def tidal_y_solver(
    model_name: str,
    radius: np.ndarray, shear_modulus: np.ndarray, bulk_modulus: np.ndarray,
    density: np.ndarray, gravity: np.ndarray, frequency: float,
    is_solid_by_layer: List[bool], is_static_by_layer: List[bool], indices_by_layer: List[np.ndarray],
    order_l: int = 2,
    surface_boundary_condition: np.ndarray = None, solve_load_numbers: bool = False,
    use_kamata: bool = False,
    use_julia: bool = False, use_numba_integrator: bool = False,
    int_rtol: float = 1.0e-8, int_atol: float = 1.0e-12,
    scipy_int_method: str = 'RK45', julia_int_method: str = 'Tsit5',
    verbose: bool = False, nondimensionalize: bool = True, planet_bulk_density: float = None
    ) -> np.ndarray:
    """ Calculate the radial solution for a homogeneous, solid planet.

    Parameters
    ----------
    model_name : str
        Interior model name (based on layer structure and desired collapse functions)
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
    is_solid_by_layer : List[bool]
        Flags for if each layer is solid (True) or liquid (False)
    is_static_by_layer : List[bool]
        Flags for if each layer is static (True) or dynamic (False)
    indices_by_layer : List[np.ndarray]
        Numpy array of booleans for the index of each layer (based on the radius array)
    order_l : int = 2
        Tidal harmonic order.
    surface_boundary_condition : np.ndarray = None
        The surface boundary condition, for tidal solutions y2, y4, y6, = (0, 0, (2l+1)/R)
            Tidal solution or load solution will be the default if `None` is provided.
    solve_load_numbers : bool = False
        If True, then the load solution will be used instead of tidal if surface_boundary_condition = None.
    use_kamata : bool = False
        If True, the Kamata+2015 initial conditions will be used at the base of layer 0.
        Otherwise, the Takeuchi & Saito 1972 initial conditions will be used.
    use_julia : bool = False
        If True, the Julia `diffeqpy` integration tools will be used.
        Otherwise, `scipy.integrate.solve_ivp` or TidalPy's numba-safe integrator will be used.
    use_numba_integrator : bool = False
        If True, TidalPy's numba-safe RK-based integrator will be used.
        Otherwise, `scipy.integrate.solve_ivp` or Julia `diffeqpy` integrator will be used.
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
    verbose: bool = False
        If True, the function will print some information to console during calculation (may cause a slow down).
    nondimensionalize : bool = False
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
    if nondimensionalize:
        if planet_bulk_density is None:
            raise AttributeNotSetError('Planet bulk density must be provided if non-dimensionalize is True.')

        radius, gravity, density, shear_modulus, bulk_modulus, frequency, G_to_use = \
            non_dimensionalize_physicals(
                radius, gravity, density, shear_modulus, bulk_modulus, frequency,
                mean_radius=planet_radius, bulk_density=planet_bulk_density
                )
    else:
        G_to_use = G

    # Find solution at the top of the planet -- this is dependent on the forcing type.
    #     Tides (default here) follow the (y2, y4, y6) = (0, 0, (2l+1)/R) rule
    if surface_boundary_condition is None:
        if solve_load_numbers:
            # TODO
            raise NotImplementedError
        else:
            surface_boundary_condition = np.zeros(3, dtype=np.complex128)
            surface_boundary_condition[2] = (2. * order_l + 1.) / radius[-1]

    # Determine layer structure
    num_layers = len(is_solid_by_layer)
    num_interfaces = num_layers - 1

    # Find the differential equation for each layer
    radial_odes = []
    ode_inputs = []
    bottom_interfaces = []
    bottom_interface_inputs = []
    num_liquid_layers = 0
    liquid_indices = None
    is_liquid_static = None
    for layer_i in range(num_layers):

        # Get the indices for this layer
        layer_indices = indices_by_layer[layer_i]

        # The known radial function ODEs are stored by is the layer solid, and is the layer static (vs dynamic)
        layer_is_solid = is_solid_by_layer[layer_i]
        layer_is_static = is_static_by_layer[layer_i]
        if frequency == 0.:
            layer_is_static = True
        layer_ode = known_multilayer_odes[(layer_is_solid, layer_is_static)]

        if not layer_is_solid:
            if num_liquid_layers > 0:
                # See TidalPy.tides.multilayer.numerical_int.collapse for more detail.
                raise NotImplementedError(
                    'The current solution collapse functions can not handle multiple liquid layers.'
                    )
            num_liquid_layers += 1
            liquid_indices = layer_indices
            is_liquid_static = layer_is_static

        # Dynamic functions require the frequency and solid layers require the shear modulus.
        if layer_is_solid:
            if layer_is_static:
                ode_input = (radius[layer_indices], shear_modulus[layer_indices], bulk_modulus[layer_indices],
                             density[layer_indices], gravity[layer_indices], order_l, G_to_use)
            else:
                ode_input = (radius[layer_indices], shear_modulus[layer_indices], bulk_modulus[layer_indices],
                             density[layer_indices], gravity[layer_indices], frequency, order_l, G_to_use)
        else:
            if layer_is_static:
                ode_input = (radius[layer_indices], density[layer_indices], gravity[layer_indices], order_l, G_to_use)
            else:
                ode_input = (radius[layer_indices], bulk_modulus[layer_indices], density[layer_indices],
                             gravity[layer_indices], frequency, order_l, G_to_use)

        # Store model info for this layer
        radial_odes.append(layer_ode)
        ode_inputs.append(ode_input)

        # Find the initial solution at the center of the planet
        if layer_i == 0:
            is_dynamic = not layer_is_static
            initial_value_tuple = \
                find_initial_guess(
                    use_kamata, layer_is_solid, is_dynamic,
                    radius[0], shear_modulus[0], bulk_modulus[0], density[0], frequency,
                    order_l=order_l, G_to_use=G_to_use
                    )

        # Determine the top interface for this layer
        if num_interfaces == 0:
            # There are no internal interfaces
            bottom_interfaces.append(None)
            bottom_interface_inputs.append(tuple())
        elif layer_i == 0:
            # The first layer has no bottom interface.
            bottom_interfaces.append(None)
            bottom_interface_inputs.append(tuple())
        else:
            # Real interface, need to figure out what type it is.
            layer_below_is_solid = is_solid_by_layer[layer_i - 1]
            layer_below_is_static = is_static_by_layer[layer_i - 1]

            # For interfaces between liquid layers we need the liquid density at the interface.
            # TODO It is a bit unclear what happens at the interface between two liquid layers of different densities
            #   This is not a common problem in TidalPy applications so we will preferentially choose the bottom
            #   layer's density in those circumstances.
            liquid_density = None
            if not layer_below_is_solid:
                # Choose the density at the top of the layer below (at the interface)
                liquid_density = density[indices_by_layer[layer_i - 1]][-1]
            elif not layer_is_solid:
                # Choose the density at the bottom of this layer (at the interface)
                liquid_density = density[layer_indices][0]
            interface_gravity = gravity[layer_indices][0]

            interface_func, interface_input = find_interface_func(
                lower_layer_is_solid=layer_below_is_solid, lower_layer_is_static=layer_below_is_static,
                upper_layer_is_solid=layer_is_solid, upper_layer_is_static=layer_is_static,
                liquid_density=liquid_density, interface_gravity=interface_gravity, G_to_use=G_to_use
                )

            bottom_interfaces.append(interface_func)
            bottom_interface_inputs.append(interface_input)

    # Determine function to collapse the multiple solutions into a single solution for the planet.
    collapse_function, collapse_input = \
        find_collapse_func(
            model_name, surface_boundary_condition,
            is_liquid_static=is_liquid_static, radius_array=radius, gravity_array=gravity,
            density_array=density, liquid_layer_indices=liquid_indices, frequency=frequency
            )

    # Ready to solve the viscoelastic-gravitational problem for each layer, obtaining multiple solutions per layer
    #    which will later be collapsed via a linear combination (subjected to boundary conditions) into one solution
    #    for the entire planet.

    # Initialize empty lists for each layer. This sub lists will contain the multiple solutions.
    # OPT: Until the solutions are collapsed they are totally independent so there is room to multiprocess here.
    #    however we can not multiprocess the layers as the initial conditions for a layer depend on the layer below.
    solutions_by_layer = [list() for layer_i in range(num_layers)]
    for layer_i in range(num_layers):

        # Pull out layer specific information
        layer_indices = indices_by_layer[layer_i]
        layer_radii = radius[layer_indices]
        radial_span = (layer_radii[0], layer_radii[-1])
        diffeq = radial_odes[layer_i]
        diffeq_input = ode_inputs[layer_i]
        bottom_interface = bottom_interfaces[layer_i]
        bottom_interface_input = bottom_interface_inputs[layer_i]

        # Determine initial conditions at the base of the layer.
        if layer_i == 0:
            # The conditions at the base of the layer were found earlier in this function using the initial conditions.
            initial_values_to_use = initial_value_tuple
        else:
            # Initial values are based on the previous layer's results and the interface function.
            initial_values_to_use = bottom_interface(solutions_by_layer[layer_i - 1], *bottom_interface_input)

        # Start integration routine
        if use_julia:
            def diffeq_julia(u, p, r):
                # Julia integrator flips the order of the variables for the differential equation.
                output = diffeq(r, u, *p)
                return list(output)

            # Import Julia's Diffeqpy and reinit the problem
            from ....utilities.julia_helper.integration_methods import get_julia_solver
            ode, solver = get_julia_solver(julia_int_method)

            if verbose:
                print(f'Solving Layer {layer_i + 1} (with SciPy, using {scipy_int_method})...')

            for solution_num, initial_values in enumerate(initial_values_to_use):
                problem = ode.ODEProblem(diffeq_julia, initial_values, radial_span, diffeq_input)
                solution = ode.solve(problem, solver(), abstol=int_atol, reltol=int_rtol)

                # Julia does not have the same t_eval. There is the "saveat" keyword but can cause issues.
                #    So perform an interpolation for the desired radii
                u_T = np.transpose(solution.u)
                u = np.zeros((u_T.shape[0], layer_radii.size), dtype=np.complex128)
                for i in range(u_T.shape[0]):
                    u[i, :] = np.interp(layer_radii, solution.t, u_T[i, :])
                solutions_by_layer[layer_i].append(u)

            if verbose:
                print('\nIntegration Done!')

        elif use_numba_integrator:

            if scipy_int_method.lower() in ['rk23']:
                rk_method = 0
            elif scipy_int_method.lower() in ['rk45']:
                rk_method = 1
            else:
                raise NotImplementedError

            if verbose:
                print(f"Solving Layer {layer_i + 1} (with TidalPy's Numba integrator, using {scipy_int_method})...")

            for solution_num, initial_values in enumerate(initial_values_to_use):

                ts, ys, success, message, = \
                    rk_integrate(
                        diffeq, radial_span, initial_values,
                        args=diffeq_input,
                        rk_method=rk_method,
                        t_eval=layer_radii,
                        rtol=int_rtol, atol=int_atol
                        )

                if not success:
                    raise IntegrationFailed(
                        f'Integration Solution Failed for {layer_i} at solution #{solution_num}.'
                        f'\n\t{message}'
                        )

                solutions_by_layer[layer_i].append(ys)

            if verbose:
                print('\nIntegration Done!')

        else:
            from scipy.integrate import solve_ivp

            if verbose:
                print(f'Solving Layer {layer_i + 1} (with SciPy, using {scipy_int_method})...')

            for solution_num, initial_values in enumerate(initial_values_to_use):
                solution = solve_ivp(
                    diffeq, radial_span, initial_values, t_eval=layer_radii, args=diffeq_input,
                    method=scipy_int_method, vectorized=False, rtol=int_rtol, atol=int_atol
                    )

                if solution.status != 0:
                    raise IntegrationFailed(
                        f'Integration Solution Failed for {layer_i} at solution #{solution_num}.'
                        f'\n\t{solution.message}'
                        )

                solutions_by_layer[layer_i].append(solution.y)

            if verbose:
                print('\nIntegration Done!')

        # Done with layer, turn the inner list into a tuple.
        solutions_by_layer[layer_i] = tuple(solutions_by_layer[layer_i])

    # Done with all layers, turn the outer list into a tuple.
    solutions_by_layer = tuple(solutions_by_layer)

    # Collapse the multiple solutions per layer into a single solution for the planet based on the boundary condition.
    if verbose:
        print('Collapsing solutions...')

    tidal_y = collapse_function(solutions_by_layer, *collapse_input)

    if verbose:
        print('Done!')

    if nondimensionalize:
        if verbose:
            print('Re-dimensionalizing Radial Functions.')
        tidal_y = re_dimensionalize_radial_func(tidal_y, planet_radius, planet_bulk_density)

    return tidal_y
