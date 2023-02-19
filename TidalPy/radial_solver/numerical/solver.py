from typing import List, Union, Tuple

import numpy as np

from TidalPy.constants import G
from TidalPy.exceptions import AttributeNotSetError, IntegrationFailed
from TidalPy.utilities.integration import get_integrator, _nb2cy, cyrk_solver
from TidalPy.utilities.performance import nbList

from .collapse import collapse_solutions
from .derivatives import known_multilayer_odes
from .initial import find_initial_guess
from .interfaces import find_interface_func
from ..nondimensional import non_dimensionalize_physicals, re_dimensionalize_radial_func


def radial_solver(
        radius: np.ndarray, shear_modulus: np.ndarray, bulk_modulus: np.ndarray,
        density: np.ndarray, gravity: np.ndarray, frequency: float, planet_bulk_density: float,
        is_solid_by_layer: Union[List[bool], Tuple[bool, ...]],
        is_static_by_layer: Union[List[bool], Tuple[bool, ...]],
        indices_by_layer: Union[List[np.ndarray], Tuple[np.ndarray, ...]],
        order_l: int = 2,
        surface_boundary_condition: np.ndarray = None, solve_load_numbers: bool = False,
        use_kamata: bool = False,
        integrator: str = 'scipy', integration_method: str = None,
        integration_rtol: float = 1.0e-6, integration_atol: float = 1.0e-8,
        verbose: bool = False, nondimensionalize: bool = True, incompressible: bool = False
        ) -> np.ndarray:
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
    planet_bulk_density : float
        Bulk density of the planet [kg m-3]
    is_solid_by_layer : Union[List[bool], Tuple[bool, ...]]
        Flags for if each layer is solid (True) or liquid (False)
    is_static_by_layer : Union[List[bool], Tuple[bool, ...]]
        Flags for if each layer is static (True) or dynamic (False)
    indices_by_layer : Union[List[np.ndarray], Tuple[np.ndarray, ...]]
        Numpy array of booleans for the index of each layer (based on the radius array)
    order_l : int = 2
        Tidal harmonic order.
    surface_boundary_condition : np.ndarray = None
        The surface boundary condition, for tidal solutions y2, y4, y6, = (0, 0, (2l+1)/R)
            Tidal solution or load solution will be the default if `None` is provided.
    solve_load_numbers : bool = False
        If True, then the load Love numbers will be calculated alongside the tidal. This changes the output signature
         of this function. If _only_ the load Love numbers are required then it will be more efficient to set this to
         False and change the `surface_boundary_condition` to the appropriate format.
    use_kamata : bool = False
        If True, the Kamata+2015 initial conditions will be used at the base of layer 0.
        Otherwise, the Takeuchi & Saito 1972 initial conditions will be used.
    integrator : str = 'scipy'
        Integrator used for solving the system of ODE's. See `TidalPy.utilities.integration` for more information.
        Depending on which packages are installed, the available options are:
            `scipy`      : SciPy's solve_ivp method.
            `cython`     : CyRK's cython-based cyrk_ode method.
            `numba`      : CyRK's numba-based nbrk_ode method.
            `julia`      : Diffeqpy's Julia-based DifferentialEquations method.
            'numbalsoda' : NumbaLSODA package.
    integration_method : str = None
        Integration method used in conjunction with the chosen integrator. If None, then the default for each integrator
        will be used (usually RK45)
    integration_rtol : float = 1.0e-6
        Integration relative error.
    integration_atol : float = 1.0e-8
        Integration absolute error.
    verbose: bool = False
        If True, the function will print some information to console during calculation (may cause a slow-down).
    nondimensionalize : bool = False
        If True, integration will use dimensionless variables. These will be converted back before output is given to
        the user.
    incompressible : bool = False
        If `True`, the incompressible assumption will be used.

    Returns
    -------
    tidal_y : np.ndarray
        The radial solution throughout the entire planet.
    (optional) load_y : np.ndarray
        The radial load solution throughout the entire planet.

    """

    # Find integrator function
    integrator, integrator_method = get_integrator(integrator, integration_method)

    # Nondimensionalize inputs
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

    # Find boundary condition at the top of the planet -- this is dependent on the forcing type.
    #     Tides (default here) follow the (y2, y4, y6) = (0, 0, (2l+1)/R) rule
    if surface_boundary_condition is None:
        surface_boundary_condition = np.zeros(3, dtype=np.complex128)
        if nondimensionalize:
            surface_boundary_condition[2] = (2. * order_l + 1.) / 1.
        else:
            surface_boundary_condition[2] = (2. * order_l + 1.) / radius[-1]

    # Determine layer structure
    num_layers = len(is_solid_by_layer)
    num_interfaces = num_layers - 1

    # Find the differential equation and other parameters for each layer and the interfaces between layers
    radial_odes = []
    ode_inputs = []
    bottom_interfaces = []
    bottom_interface_inputs = []
    gravity_at_interfaces = []
    liquid_density_at_interfaces = []
    for layer_i in range(num_layers):

        # Get the indices for this layer
        layer_indices = indices_by_layer[layer_i]

        # The known radial function ODEs are stored by is the layer solid, and is the layer static (vs dynamic)
        layer_is_solid = is_solid_by_layer[layer_i]
        layer_is_static = is_static_by_layer[layer_i]
        if frequency == 0.:
            # There will be divide by zero errors in dynamic layers.
            # TODO: Perhaps an error should be thrown instead?
            layer_is_static = True
        layer_ode = known_multilayer_odes[(layer_is_solid, layer_is_static)]

        # Dynamic functions require the frequency and solid layers require the shear modulus.
        if layer_is_solid:
            if layer_is_static:
                ode_input = (radius[layer_indices], shear_modulus[layer_indices], bulk_modulus[layer_indices],
                             density[layer_indices], gravity[layer_indices], order_l, G_to_use, incompressible)
            else:
                ode_input = (radius[layer_indices], shear_modulus[layer_indices], bulk_modulus[layer_indices],
                             density[layer_indices], gravity[layer_indices], frequency, order_l, G_to_use,
                             incompressible)
        else:
            if layer_is_static:
                ode_input = (radius[layer_indices], density[layer_indices], gravity[layer_indices], order_l, G_to_use,
                             incompressible)
            else:
                ode_input = (radius[layer_indices], bulk_modulus[layer_indices], density[layer_indices],
                             gravity[layer_indices], frequency, order_l, G_to_use, incompressible)

        # If cython solver is used, convert function type now
        if integrator is cyrk_solver:
            layer_ode = _nb2cy(layer_ode, use_njit=True, cache_njit=True)

        # Store model info for this layer
        radial_odes.append(layer_ode)
        ode_inputs.append(ode_input)

        # Find the initial solution at the center of the planet
        if layer_i == 0:
            initial_value_tuple = \
                find_initial_guess(
                        layer_is_solid, layer_is_static, incompressible, use_kamata,
                        radius[0], shear_modulus[0], bulk_modulus[0], density[0], frequency,
                        order_l=order_l, G_to_use=G_to_use
                        )

        # Determine the bottom interface for this layer
        if num_interfaces == 0:
            # There are no internal interfaces
            bottom_interfaces.append(None)
            bottom_interface_inputs.append(tuple())
            gravity_at_interfaces.append(0.)
            liquid_density_at_interfaces.append(0.)
        elif layer_i == 0:
            # The first layer has no bottom interface.
            bottom_interfaces.append(None)
            bottom_interface_inputs.append(tuple())
            gravity_at_interfaces.append(0.)
            liquid_density_at_interfaces.append(0.)
        else:
            # Real interface, need to figure out what type it is.
            layer_below_is_solid = is_solid_by_layer[layer_i - 1]
            layer_below_is_static = is_static_by_layer[layer_i - 1]

            # For interfaces between liquid layers we need the liquid density at the interface.
            # TODO It is a bit unclear what happens at the interface between two liquid layers of different densities
            #   This is not a common problem in TidalPy applications so we will preferentially choose the bottom
            #   layer's density in those circumstances.
            if layer_below_is_solid:
                if layer_is_solid:
                    # Both layers are solid, no liquid density is needed.
                    static_liquid_density = None
                else:
                    # This layer is liquid. Grab the density at the base of this layer.
                    static_liquid_density = density[layer_indices][0]
            else:
                if layer_is_solid:
                    # Layer below is liquid. Grab the density at the top of the below layer.
                    static_liquid_density = density[indices_by_layer[layer_i - 1]][-1]
                else:
                    # If both layers are liquid, choose the density of the static layer.
                    if layer_is_static:
                        if layer_below_is_static:
                            # TODO: Both layers are static. Not sure what to do here so just choose this layer's value.
                            static_liquid_density = density[layer_indices][0]
                        else:
                            # Layer below is dynamic. Grab this layer's density,
                            static_liquid_density = density[layer_indices][0]
                    else:
                        if layer_below_is_static:
                            # Layer below is static, grab its density.
                            static_liquid_density = density[indices_by_layer[layer_i - 1]][-1]
                        else:
                            # Both layers are dynamic. Liquid density will not be needed.
                            static_liquid_density = None
            # For the gravity, take the average between the bottom of this layer and the top of the layer below.
            gravity_bot_this_layer = gravity[layer_indices][0]
            gravity_top_lower_layer = gravity[indices_by_layer[layer_i - 1]][-1]
            interface_gravity = 0.5 * (gravity_bot_this_layer + gravity_top_lower_layer)

            # Record interface parameters
            gravity_at_interfaces.append(interface_gravity)
            liquid_density_to_store = static_liquid_density
            if liquid_density_to_store is None:
                liquid_density_to_store = 0.
            liquid_density_at_interfaces.append(liquid_density_to_store)

            # Find interface function and input
            interface_func, interface_input = find_interface_func(
                    lower_layer_is_solid=layer_below_is_solid, lower_layer_is_static=layer_below_is_static,
                    upper_layer_is_solid=layer_is_solid, upper_layer_is_static=layer_is_static,
                    static_liquid_density=static_liquid_density, interface_gravity=interface_gravity, G_to_use=G_to_use
                    )

            bottom_interfaces.append(interface_func)
            bottom_interface_inputs.append(interface_input)

    # Ready to solve the viscoelastic-gravitational problem for each layer, obtaining multiple solutions per layer
    #    which will later be collapsed via a linear combination (subjected to boundary conditions) into one solution
    #    for the entire planet.

    # Initialize empty lists for each layer. This sub lists will contain the multiple solutions.
    # OPT: Until the solutions are collapsed they are totally independent so there is room to multiprocess here.
    #    however we can not multiprocess the layers as the initial conditions for a layer depend on the layer below.
    solutions_by_layer = [list() for _ in range(num_layers)]
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
            lower_layer_top_solutions = nbList()
            for lower_layer_solution in solutions_by_layer[layer_i - 1]:
                lower_layer_top_solutions.append(lower_layer_solution[:, -1])
            initial_values_to_use = bottom_interface(lower_layer_top_solutions, *bottom_interface_input)

        if verbose:
            print(f'Solving Layer {layer_i + 1} (using {integrator} with {integrator_method})...')

        # Integrate over each solution
        for solution_num, initial_values in enumerate(initial_values_to_use):
            # Convert initial values to floats
            num_initial_values = initial_values.size
            initial_values_float = np.empty(2 * num_initial_values, dtype=np.float64)
            for y_i in range(num_initial_values):
                initial_values_float[2 * y_i] = np.real(initial_values[y_i])
                initial_values_float[2 * y_i + 1] = np.imag(initial_values[y_i])

            # Start integration routine
            radius_domain, y_results_float, success, message = \
                integrator(
                    diffeq, radial_span, initial_values_float, args=diffeq_input,
                    rtol=integration_rtol, atol=integration_atol, method=integrator_method, t_eval=layer_radii)

            if not success:
                raise IntegrationFailed(
                    f'Integration Solution Failed for layer {layer_i} at solution #{solution_num}.'
                    f'\n\t{message}')

            # Convert floats back to complex
            len_r = radius_domain.size
            y_results = np.empty((num_initial_values, len_r), dtype=np.complex128)
            for y_i in range(num_initial_values):
                y_results[y_i, :] = y_results_float[2 * y_i, :] + 1.0j * y_results_float[2 * y_i + 1, :]

            # Store result for layer
            solutions_by_layer[layer_i].append(y_results)

        # solutions_by_layer[layer_i] = tuple(solutions_by_layer[layer_i])

        if verbose:
            print('\nIntegration Done!')

    # Done with all layers, turn the outer list into a tuple.
    solutions_by_layer = tuple(solutions_by_layer)

    # Collapse the multiple solutions per layer into a single solution for the planet based on the boundary condition.
    if verbose:
        print('Collapsing solutions...')

    tidal_y = collapse_solutions(
            solutions_by_layer,
            is_solid_by_layer, is_static_by_layer, indices_by_layer,
            surface_boundary_condition,
            radius, density, gravity,
            gravity_at_interfaces, liquid_density_at_interfaces,
            gravity[-1], frequency, G_to_use=G_to_use
            )

    if nondimensionalize:
        if verbose:
            print('Re-dimensionalizing Radial Functions.')
        tidal_y = re_dimensionalize_radial_func(tidal_y, planet_radius, planet_bulk_density)

    if solve_load_numbers:
        # In additional to the main calculation, also calculate the radial solution for surface loading.
        if verbose:
            print('Collapsing load solutions...')

        # Find boundary condition at the top of the planet for surface loading
        surface_loading_boundary_condition = np.zeros(3, dtype=np.complex128)
        if nondimensionalize:
            surface_loading_boundary_condition[0] = -(2. * order_l + 1.) / 3.
            surface_loading_boundary_condition[2] = (2. * order_l + 1.) / 1.
        else:
            surface_loading_boundary_condition[0] = -(2. * order_l + 1.) * planet_bulk_density / 3.
            surface_loading_boundary_condition[2] = (2. * order_l + 1.) / radius[-1]

        load_y = collapse_solutions(
                solutions_by_layer,
                is_solid_by_layer, is_static_by_layer, indices_by_layer,
                surface_loading_boundary_condition,
                radius, density, gravity,
                gravity_at_interfaces, liquid_density_at_interfaces,
                gravity[-1], frequency, G_to_use=G_to_use
                )
        if nondimensionalize:
            if verbose:
                print('Re-dimensionalizing Radial (Load) Functions.')
            load_y = re_dimensionalize_radial_func(load_y, planet_radius, planet_bulk_density)

        output_ = (tidal_y, load_y)
    else:
        output_ = tidal_y

    if verbose:
        print('Done!')

    return output_
