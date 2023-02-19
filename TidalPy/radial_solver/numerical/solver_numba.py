from typing import Tuple, List, Union

import numpy as np

from TidalPy.utilities.integration import _nbrk_ode, cyrk_installed
from TidalPy.utilities.performance import nbList, njit

from .collapse import collapse_solutions
from .derivatives import dynamic_liquid_ode, dynamic_solid_ode, static_liquid_ode, static_solid_ode
from .initial import find_initial_guess
from .interfaces import (
    interface_LDy_LDy, interface_LDy_SDy, interface_LDy_SSt, interface_LSt_LSt, interface_LSt_SDy,
    interface_LSt_SSt, interface_SDy_LDy, interface_SDy_LSt, interface_SDy_SDy, interface_SDy_SSt,
    interface_SSt_LDy, interface_LSt_LDy, interface_LDy_LSt,
    interface_SSt_LSt, interface_SSt_SDy, interface_SSt_SSt)
from ..nondimensional import non_dimensionalize_physicals, re_dimensionalize_radial_func


@njit(cacheable=True)
def _get_initial_values(
        layer_i, layer_is_solid, layer_is_static, layer_below_is_solid, layer_below_is_static,
        layer_below_y_solutions,
        radius_at_bottom, shear_mod_at_bottom, bulk_mod_at_bottom, density_at_bottom, gravity_at_bottom,
        frequency, interface_gravity, liquid_density,
        order_l, G_to_use, use_kamata, incompressible
        ):
    initial_values_to_use = nbList([np.empty((6, 0), dtype=np.complex128)])
    if layer_i == 0:
        # Find the initial solution at the center of the planet
        initial_values_to_use = \
            find_initial_guess(
                    layer_is_solid, layer_is_static, incompressible, use_kamata,
                    radius_at_bottom, shear_mod_at_bottom, bulk_mod_at_bottom, density_at_bottom, frequency,
                    order_l=order_l, G_to_use=G_to_use
                    )
    else:
        # Find solutions from the top of the layer below.
        layer_below_y_solutions_top = nbList([layer_below_y[:, -1] for layer_below_y in layer_below_y_solutions])

        # For interfaces between liquid layers we need the liquid density at the interface.
        if layer_below_is_solid:
            if layer_is_solid:
                if layer_below_is_static:
                    if layer_is_static:
                        # Solid-Solid, Static-Static
                        initial_values_to_use = interface_SSt_SSt(layer_below_y_solutions_top)
                    else:
                        # Solid-Solid, Static-Dynamic
                        initial_values_to_use = interface_SSt_SDy(layer_below_y_solutions_top)
                else:
                    if layer_is_static:
                        # Solid-Solid, Dynamic-Static
                        initial_values_to_use = interface_SDy_SSt(layer_below_y_solutions_top)
                    else:
                        # Solid-Solid, Dynamic-Dynamic
                        initial_values_to_use = interface_SDy_SDy(layer_below_y_solutions_top)
            else:
                if layer_below_is_static:
                    if layer_is_static:
                        # Solid-Liquid, Static-Static
                        initial_values_to_use = interface_SSt_LSt(
                                layer_below_y_solutions_top,
                                interface_gravity, liquid_density, G_to_use
                                )
                    else:
                        # Solid-Liquid, Static-Dynamic
                        initial_values_to_use = interface_SSt_LDy(layer_below_y_solutions_top)
                else:
                    if layer_is_static:
                        # Solid-Liquid, Dynamic-Static
                        initial_values_to_use = interface_SDy_LSt(
                                layer_below_y_solutions_top,
                                interface_gravity, liquid_density, G_to_use
                                )
                    else:
                        # Solid-Liquid, Dynamic-Dynamic
                        initial_values_to_use = interface_SDy_LDy(layer_below_y_solutions_top)
        else:
            if layer_is_solid:
                if layer_below_is_static:
                    if layer_is_static:
                        # Liquid-Solid, Static-Static
                        initial_values_to_use = interface_LSt_SSt(
                                layer_below_y_solutions_top,
                                interface_gravity, liquid_density, G_to_use
                                )
                    else:
                        # Liquid-Solid, Static-Dynamic
                        initial_values_to_use = interface_LSt_SDy(
                                layer_below_y_solutions_top,
                                interface_gravity, liquid_density, G_to_use
                                )
                else:
                    if layer_is_static:
                        # Liquid-Solid, Dynamic-Static
                        initial_values_to_use = interface_LDy_SSt(layer_below_y_solutions_top)
                    else:
                        # Liquid-Solid, Dynamic-Dynamic
                        initial_values_to_use = interface_LDy_SDy(layer_below_y_solutions_top)
            else:
                if layer_below_is_static:
                    if layer_is_static:
                        # Liquid-Liquid, Static-Static
                        initial_values_to_use = interface_LSt_LSt(layer_below_y_solutions_top)
                    else:
                        # Liquid-Liquid, Static-Dynamic
                        initial_values_to_use = interface_LSt_LDy(
                                layer_below_y_solutions_top,
                                interface_gravity, liquid_density, G_to_use
                                )
                else:
                    if layer_is_static:
                        # Liquid-Liquid, Dynamic-Static
                        initial_values_to_use = interface_LDy_LSt(
                                layer_below_y_solutions_top,
                                interface_gravity, liquid_density, G_to_use
                                )
                    else:
                        # Liquid-Liquid, Dynamic-Dynamic
                        initial_values_to_use = interface_LDy_LDy(layer_below_y_solutions_top)

    return initial_values_to_use


@njit(cacheable=True)
def _single_layer_integrate(
        layer_is_solid, layer_is_static, radial_span, initial_values_to_use, frequency,
        radius_array, shear_modulus_array, bulk_modulus_array, density_array,
        gravity_array, order_l, G_to_use, incompressible, integration_rtol, integration_atol, rk_method
        ):

    len_r = radius_array.size
    solution_num = 0
    layer_solutions = nbList()
    for initial_values in initial_values_to_use:
        # Convert initial values from n complex numbers to 2n floats
        num_initial_values = initial_values.size
        initial_values_float = np.empty(2 * num_initial_values, dtype=np.float64)
        for y_i in range(num_initial_values):
            initial_values_float[2 * y_i] = np.real(initial_values[y_i])
            initial_values_float[2 * y_i + 1] = np.imag(initial_values[y_i])

        if layer_is_solid:
            if layer_is_static:
                ts, ys, success, message, = \
                    _nbrk_ode(
                            static_solid_ode, radial_span, initial_values_float,
                            args=(radius_array, shear_modulus_array, bulk_modulus_array,
                                  density_array, gravity_array, order_l, G_to_use, incompressible),
                            rk_method=rk_method,
                            t_eval=radius_array,
                            rtol=integration_rtol, atol=integration_atol
                            )
            else:
                ts, ys, success, message, = \
                    _nbrk_ode(
                            dynamic_solid_ode, radial_span, initial_values_float,
                            args=(radius_array, shear_modulus_array, bulk_modulus_array,
                                  density_array, gravity_array, frequency, order_l, G_to_use, incompressible),
                            rk_method=rk_method,
                            t_eval=radius_array,
                            rtol=integration_rtol, atol=integration_atol
                            )
        else:
            if layer_is_static:
                ts, ys, success, message, = \
                    _nbrk_ode(
                            static_liquid_ode, radial_span, initial_values_float,
                            args=(radius_array, density_array,
                                  gravity_array, order_l, G_to_use, incompressible),
                            rk_method=rk_method,
                            t_eval=radius_array,
                            rtol=integration_rtol, atol=integration_atol
                            )
            else:
                ts, ys, success, message, = \
                    _nbrk_ode(
                            dynamic_liquid_ode, radial_span, initial_values_float,
                            args=(radius_array, bulk_modulus_array, density_array,
                                  gravity_array, frequency, order_l, G_to_use, incompressible),
                            rk_method=rk_method,
                            t_eval=radius_array,
                            rtol=integration_rtol, atol=integration_atol
                            )

        if not success:
            print(message)
            raise Exception('Integration Solution Failed.')

        # Convert floats back to complex
        ys_complex = np.empty((num_initial_values, len_r), dtype=np.complex128)
        for y_i in range(num_initial_values):
            for r_i in range(len_r):
                ys_complex[y_i, r_i] = ys[2 * y_i, r_i] + 1.0j * ys[2 * y_i + 1, r_i]

        # Store solution
        layer_solutions.append(ys_complex)
        solution_num += 1

    return layer_solutions


@njit(cacheable=True)
def radial_solver(
        radius: np.ndarray, shear_modulus: np.ndarray, bulk_modulus: np.ndarray,
        density: np.ndarray, gravity: np.ndarray, frequency: float, planet_bulk_density: float,
        is_solid_by_layer: Union[List[bool], Tuple[bool, ...]],
        is_static_by_layer: Union[List[bool], Tuple[bool, ...]],
        indices_by_layer: Union[List[np.ndarray], Tuple[np.ndarray, ...]],
        order_l: int = 2,
        surface_boundary_condition: np.ndarray = None, solve_load_numbers: bool = False,
        use_kamata: bool = False,
        integration_rtol: float = 1.0e-8, integration_atol: float = 1.0e-12, integration_method: int = 1,
        verbose: bool = False, nondimensionalize: bool = True,
        incompressible: bool = False
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
        If True, then the load solution will be used instead of tidal if surface_boundary_condition = None.
        OPT: Unlike the non-numba version, this function is not able to calculate the tidal and load numbers
         simultaneously. This is only due to the return signature would be different like the non-numba function.
    use_kamata : bool = False
        If True, the Kamata+2015 initial conditions will be used at the base of layer 0.
        Otherwise, the Takeuchi & Saito 1972 initial conditions will be used.
    integration_rtol : float = 1.0e-6
        Integration relative error.
    integration_atol : float = 1.0e-4
        Integration absolute error.
    integration_method : int = 1
    The type of RK method used for integration
        0 = RK23
        1 = RK45
        2 = DOP
    verbose: bool = False
        If True, the function will print some information to console during calculation (may cause a slow down).
    nondimensionalize : bool = False
        If True, integration will use dimensionless variables. These will be converted back before output is given to
        the user.
    incompressible : bool = False
        If `True`, the incompressible assumption will be used.

    Returns
    -------
    tidal_y : np.ndarray
        The radial solution throughout the entire planet.

    """

    if not cyrk_installed:
        raise Exception('CyRK package is required for the numba-based solver.')

    # Non-dimensionalize inputs
    planet_radius = radius[-1]
    if nondimensionalize:
        if planet_bulk_density is None:
            raise Exception('Planet bulk density must be provided if non-dimensionalize is True.')

        radius, gravity, density, shear_modulus, bulk_modulus, frequency, G_to_use = \
            non_dimensionalize_physicals(
                    radius, gravity, density, shear_modulus, bulk_modulus, frequency,
                    mean_radius=planet_radius, bulk_density=planet_bulk_density
                    )
    else:
        # Gravitational Constant
        G_to_use = 6.6743e-11

    # Find solution at the top of the planet -- this is dependent on the forcing type.
    #     Tides (default here) follow the (y2, y4, y6) = (0, 0, (2l+1)/R) rule
    if surface_boundary_condition is None:
        if solve_load_numbers:
            if planet_bulk_density is None:
                raise Exception('Planet bulk density must be provided if calculating load Love numbers.')
            # Based on Eq. 6 of Beuthe (2015; Icarus)
            surface_boundary_condition = np.zeros(3, dtype=np.complex128)
            if nondimensionalize:
                surface_boundary_condition[0] = -(2. * order_l + 1.) * 1. / 3.
                surface_boundary_condition[2] = (2. * order_l + 1.) / 1.
            else:
                surface_boundary_condition[0] = -(2. * order_l + 1.) * planet_bulk_density / 3.
                surface_boundary_condition[2] = (2. * order_l + 1.) / radius[-1]
        else:
            surface_boundary_condition = np.zeros(3, dtype=np.complex128)
            if nondimensionalize:
                surface_boundary_condition[2] = (2. * order_l + 1.) / 1.
            else:
                surface_boundary_condition[2] = (2. * order_l + 1.) / radius[-1]

    # Determine layer structure
    num_layers = len(is_solid_by_layer)
    num_interfaces = num_layers - 1
    gravity_at_interfaces = nbList()
    liquid_density_at_interfaces = nbList()

    # Find the differential equation for each layer
    num_liquid_layers = 0

    # Dummy variables for compile
    layer_below_ys = nbList([np.empty((6, 1), dtype=np.complex128)])

    # Ready to solve the viscoelastic-gravitational problem for each layer, obtaining multiple solutions per layer
    #    which will later be collapsed via a linear combination (subjected to boundary conditions) into one solution
    #    for the entire planet.

    # Initialize empty lists for each layer. This sub lists will contain the multiple solutions.
    # OPT: Until the solutions are collapsed they are totally independent so there is room to multiprocess here.
    #    however we can not multiprocess the layers as the initial conditions for a layer depend on the layer below.
    solutions_by_layer = nbList()
    for layer_i in range(num_layers):

        # Pull out layer specific information
        layer_indices = indices_by_layer[layer_i]
        layer_radii = radius[layer_indices]
        radial_span = (layer_radii[0], layer_radii[-1])

        layer_density = density[layer_indices]
        layer_gravity = gravity[layer_indices]
        layer_shear = shear_modulus[layer_indices]
        layer_bulk = bulk_modulus[layer_indices]

        # The known radial function ODEs are stored by is the layer solid, and is the layer static (vs dynamic)
        layer_is_solid = is_solid_by_layer[layer_i]
        layer_is_static = is_static_by_layer[layer_i]
        if frequency == 0.:
            layer_is_static = True

        # Determine initial conditions at the base of the layer.
        if layer_i == 0:
            layer_below_is_solid = False
            layer_below_is_static = False
            layer_below_top_density = 0.
        else:
            layer_below_is_solid = is_solid_by_layer[layer_i - 1]
            layer_below_is_static = is_static_by_layer[layer_i - 1]
            layer_below_top_density = density[layer_indices - 1][-1]

        # Determine the bottom interface for this layer
        if num_interfaces == 0:
            # There are no internal interfaces
            interface_gravity = 0.
            static_liquid_density = 0.
            gravity_at_interfaces.append(interface_gravity)
            liquid_density_at_interfaces.append(static_liquid_density)
        elif layer_i == 0:
            # The first layer has no bottom interface.
            interface_gravity = 0.
            static_liquid_density = 0.
            gravity_at_interfaces.append(interface_gravity)
            liquid_density_at_interfaces.append(static_liquid_density)
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
                    static_liquid_density = 0.
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
                            static_liquid_density = 0.
            # For the gravity, take the average between the bottom of this layer and the top of the layer below.
            gravity_bot_this_layer = gravity[layer_indices][0]
            gravity_top_lower_layer = gravity[indices_by_layer[layer_i - 1]][-1]
            interface_gravity = 0.5 * (gravity_bot_this_layer + gravity_top_lower_layer)

            # Record interface parameters
            gravity_at_interfaces.append(interface_gravity)
            liquid_density_at_interfaces.append(static_liquid_density)

        initial_values_to_use = \
            _get_initial_values(
                    layer_i, layer_is_solid, layer_is_static, layer_below_is_solid, layer_below_is_static,
                    layer_below_ys,
                    layer_radii[0], layer_shear[0], layer_bulk[0], layer_density[0], layer_gravity[0],
                    frequency, interface_gravity, static_liquid_density,
                    order_l, G_to_use, use_kamata, incompressible
                    )

        # Start integration routine
        if verbose:
            msg = 'Solving Layer ' + str(layer_i + 1) + " (with TidalPy's Numba integrator)..."
            print(msg)

        layer_below_ys = \
            _single_layer_integrate(
                    layer_is_solid, layer_is_static, radial_span, initial_values_to_use, frequency,
                    layer_radii, layer_shear,
                    layer_bulk, layer_density,
                    layer_gravity, order_l, G_to_use,
                    incompressible, integration_rtol, integration_atol, integration_method
                    )

        # Add solutions to outer list
        solutions_by_layer.append(layer_below_ys)

        if verbose:
            print('\nIntegration Done!')

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

    if verbose:
        print('Done!')

    if nondimensionalize:
        if verbose:
            print('Re-dimensionalizing Radial Functions.')
        tidal_y = re_dimensionalize_radial_func(tidal_y, planet_radius, planet_bulk_density)

    return tidal_y
