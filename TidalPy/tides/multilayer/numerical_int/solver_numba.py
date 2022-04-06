from typing import List as ListType

import numpy as np

from .collapse import (collapse_homogen_solid, collapse_ls_dynamic_liq, collapse_ls_static_liq,
                       collapse_sls_dynamic_liq,
                       collapse_sls_static_liq, collapse_ssls_dynamic_liq, collapse_ssls_static_liq)
from .derivatives import dynamic_liquid_ode, dynamic_solid_ode, static_liquid_ode, static_solid_ode
from .initial_conditions import find_initial_guess
from .interfaces import (interface_LDy_LDy, interface_LDy_SDy, interface_LDy_SSt, interface_LSt_LSt, interface_LSt_SDy,
                         interface_LSt_SSt, interface_SDy_LDy, interface_SDy_LSt, interface_SDy_SDy, interface_SDy_SSt,
                         interface_SSt_LDy,
                         interface_SSt_LSt, interface_SSt_SDy, interface_SSt_SSt)
from ..nondimensional import non_dimensionalize_physicals, re_dimensionalize_radial_func
from ....constants import G
from ....utilities.integration.rk_integrator import rk_integrate
from ....utilities.performance import njit, nbList


@njit(cacheable=True)
def _get_initial_values(
    layer_i, layer_is_solid, layer_is_static, layer_below_is_solid, layer_below_is_static,
    layer_below_top_density, last_below_ys,
    radius_array, shear_mod_array, bulk_mod_array, density_array, gravity_array, frequency,
    order_l, G_to_use, use_kamata
    ):
    initial_values_to_use = nbList([np.empty(6, dtype=np.complex128)])
    if layer_i == 0:
        # Find the initial solution at the center of the planet
        is_dynamic = not layer_is_static
        initial_values_to_use = \
            find_initial_guess(
                use_kamata, layer_is_solid, is_dynamic,
                radius_array[0], shear_mod_array[0], bulk_mod_array[0], density_array[0], frequency,
                order_l=order_l, G_to_use=G_to_use
                )
    else:
        # For interfaces between liquid layers we need the liquid density at the interface.
        # TODO It is a bit unclear what happens at the interface between two liquid layers of different densities
        #   This is not a common problem in TidalPy applications so we will preferentially choose the bottom
        #   layer's density in those circumstances
        liquid_density = None
        if not layer_below_is_solid:
            # Choose the density at the top of the layer below (at the interface)
            liquid_density = layer_below_top_density
        elif not layer_is_solid:
            # Choose the density at the bottom of this layer (at the interface)
            liquid_density = density_array[0]
        interface_gravity = gravity_array[0]
        if layer_below_is_solid:
            if layer_is_solid:
                if layer_below_is_static:
                    if layer_is_static:
                        # Solid-Solid, Static-Static
                        initial_values_to_use = interface_SSt_SSt(last_below_ys)
                    else:
                        # Solid-Solid, Static-Dynamic
                        initial_values_to_use = interface_SSt_SDy(last_below_ys)
                else:
                    if layer_is_static:
                        # Solid-Solid, Dynamic-Static
                        initial_values_to_use = interface_SDy_SSt(last_below_ys)
                    else:
                        # Solid-Solid, Dynamic-Dynamic
                        initial_values_to_use = interface_SDy_SDy(last_below_ys)
            else:
                if layer_below_is_static:
                    if layer_is_static:
                        # Solid-Liquid, Static-Static
                        initial_values_to_use = interface_SSt_LSt(
                            last_below_ys,
                            interface_gravity, liquid_density, G_to_use
                            )
                    else:
                        # Solid-Liquid, Static-Dynamic
                        initial_values_to_use = interface_SSt_LDy(last_below_ys)
                else:
                    if layer_is_static:
                        # Solid-Liquid, Dynamic-Static
                        initial_values_to_use = interface_SDy_LSt(
                            last_below_ys,
                            interface_gravity, liquid_density, G_to_use
                            )
                    else:
                        # Solid-Liquid, Dynamic-Dynamic
                        initial_values_to_use = interface_SDy_LDy(last_below_ys)
        else:
            if layer_is_solid:
                if layer_below_is_static:
                    if layer_is_static:
                        # Liquid-Solid, Static-Static
                        initial_values_to_use = interface_LSt_SSt(
                            last_below_ys,
                            interface_gravity, liquid_density, G_to_use
                            )
                    else:
                        # Liquid-Solid, Static-Dynamic
                        initial_values_to_use = interface_LSt_SDy(
                            last_below_ys,
                            interface_gravity, liquid_density, G_to_use
                            )
                else:
                    if layer_is_static:
                        # Liquid-Solid, Dynamic-Static
                        initial_values_to_use = interface_LDy_SSt(last_below_ys)
                    else:
                        # Liquid-Solid, Dynamic-Dynamic
                        initial_values_to_use = interface_LDy_SDy(last_below_ys)
            else:
                if layer_below_is_static:

                    if layer_is_static:
                        # Liquid-Liquid, Static-Static
                        initial_values_to_use = interface_LSt_LSt(last_below_ys)
                    else:
                        # Liquid-Liquid, Static-Dynamic
                        initial_values_to_use = interface_LSt_LSt(last_below_ys)
                        raise Exception('Not Implemented')
                else:
                    if layer_is_static:
                        # Liquid-Liquid, Dynamic-Static
                        initial_values_to_use = interface_LSt_LSt(last_below_ys)
                        raise Exception('Not Implemented')
                    else:
                        # Liquid-Liquid, Dynamic-Dynamic
                        initial_values_to_use = interface_LDy_LDy(last_below_ys)

    return initial_values_to_use


@njit(cacheable=True)
def _single_layer_integrate(
    layer_is_solid, layer_is_static, radial_span, initial_values_to_use, frequency,
    radius_array, shear_modulus_array, bulk_modulus_array, density_array,
    gravity_array, order_l, G_to_use, rk_method, int_rtol, int_atol
    ):
    solution_num = 0
    layer_solutions = nbList()
    for initial_values in initial_values_to_use:
        initial_values_copy = initial_values.copy()

        if layer_is_solid:
            if layer_is_static:
                ts, ys, success, message, = \
                    rk_integrate(
                        static_solid_ode, radial_span, initial_values_copy,
                        args=(radius_array, shear_modulus_array, bulk_modulus_array,
                              density_array, gravity_array, order_l, G_to_use),
                        rk_method=rk_method,
                        t_eval=radius_array,
                        rtol=int_rtol, atol=int_atol
                        )
            else:
                ts, ys, success, message, = \
                    rk_integrate(
                        dynamic_solid_ode, radial_span, initial_values_copy,
                        args=(radius_array, shear_modulus_array, bulk_modulus_array,
                              density_array, gravity_array, frequency, order_l, G_to_use),
                        rk_method=rk_method,
                        t_eval=radius_array,
                        rtol=int_rtol, atol=int_atol
                        )
        else:
            if layer_is_static:
                ts, ys, success, message, = \
                    rk_integrate(
                        static_liquid_ode, radial_span, initial_values_copy,
                        args=(radius_array, density_array,
                              gravity_array, order_l, G_to_use),
                        rk_method=rk_method,
                        t_eval=radius_array,
                        rtol=int_rtol, atol=int_atol
                        )
            else:
                ts, ys, success, message, = \
                    rk_integrate(
                        dynamic_liquid_ode, radial_span, initial_values_copy,
                        args=(radius_array, bulk_modulus_array, density_array,
                              gravity_array, frequency, order_l, G_to_use),
                        rk_method=rk_method,
                        t_eval=radius_array,
                        rtol=int_rtol, atol=int_atol
                        )

        if not success:
            raise Exception('Integration Solution Failed.')

        layer_solutions.append(ys)
        solution_num += 1

    return layer_solutions


@njit(cacheable=True)
def tidal_y_solver(
    model_name: str,
    radius: np.ndarray, shear_modulus: np.ndarray, bulk_modulus: np.ndarray,
    density: np.ndarray, gravity: np.ndarray, frequency: float,
    is_solid_by_layer: ListType[bool], is_static_by_layer: ListType[bool], indices_by_layer: ListType[np.ndarray],
    order_l: int = 2,
    surface_boundary_condition: np.ndarray = None, solve_load_numbers: bool = False,
    use_kamata: bool = False,
    int_rtol: float = 1.0e-8, int_atol: float = 1.0e-12, rk_method: int = 1,
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
    int_rtol : float = 1.0e-6
        Integration relative error.
    int_atol : float = 1.0e-4
        Integration absolute error.
    rk_method : int = 1
    The type of RK method used for integration
        0 = RK23
        1 = RK45
        # TODO NotImplemented 2 = DOP
    verbose: bool = False
        If True, the function will print some information to console during calculation (may cause a slow down).
    nondimensionalize : bool = False
        If True, integration will use dimensionless variables. These will be converted back before output is given to
        the user.
    planet_bulk_density : float = None
        Must be provided if nondimensionalize is True. Bulk density of the planet.

    Returns
    -------
    tidal_y : np.ndarray
        The radial solution throughout the entire planet.

    """

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
        G_to_use = G

    # Find solution at the top of the planet -- this is dependent on the forcing type.
    #     Tides (default here) follow the (y2, y4, y6) = (0, 0, (2l+1)/R) rule
    if surface_boundary_condition is None:
        if solve_load_numbers:
            # TODO
            raise Exception('Not Implemented.')
        else:
            surface_boundary_condition = np.zeros(3, dtype=np.complex128)
            surface_boundary_condition[2] = (2. * order_l + 1.) / radius[-1]

    # Determine layer structure
    num_layers = len(is_solid_by_layer)
    num_interfaces = num_layers - 1

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

        if not layer_is_solid:
            if num_liquid_layers > 0:
                # See TidalPy.tides.multilayer.numerical_int.collapse for more detail.
                raise Exception('The current solution collapse functions can not handle multiple liquid layers.')
            num_liquid_layers += 1
            liquid_indices = layer_indices
            is_liquid_static = layer_is_static

        # Determine initial conditions at the base of the layer.
        if layer_i == 0:
            layer_below_is_solid = False
            layer_below_is_static = False
            layer_below_top_density = 0.
        else:
            layer_below_is_solid = is_solid_by_layer[layer_i - 1]
            layer_below_is_static = is_static_by_layer[layer_i - 1]
            layer_below_top_density = density[layer_indices - 1][-1]

        initial_values_to_use = \
            _get_initial_values(
                layer_i, layer_is_solid, layer_is_static, layer_below_is_solid, layer_below_is_static,
                layer_below_top_density, layer_below_ys,
                layer_radii, layer_shear, layer_bulk, layer_density, layer_gravity, frequency,
                order_l, G_to_use, use_kamata
                )

        # Start integration routine
        if verbose:
            print(f"Solving Layer {layer_i + 1} (with TidalPy's Numba integrator)...")

        layer_below_ys = \
            _single_layer_integrate(
                layer_is_solid, layer_is_static, radial_span, initial_values_to_use, frequency,
                layer_radii, layer_shear, layer_bulk, layer_density,
                layer_gravity, order_l, G_to_use, rk_method, int_rtol, int_atol
                )

        # Add solutions to outer list
        solutions_by_layer.append(layer_below_ys)

        if verbose:
            print('\nIntegration Done!')

    # Collapse the multiple solutions per layer into a single solution for the planet based on the boundary condition.
    if verbose:
        print('Collapsing solutions...')

    # Determine function to collapse the multiple solutions into a single solution for the planet.
    if model_name == 'homogeneous_solid':
        # No liquid model
        tidal_y = collapse_homogen_solid(solutions_by_layer, surface_boundary_condition)
    else:
        liquid_gravity = gravity[liquid_indices]
        liquid_density = density[liquid_indices]
        liquid_radii = radius[liquid_indices]
        if model_name == 'liquid_solid':
            if is_liquid_static:
                # The static liquid case for a liquid-solid model requires no other input.
                tidal_y = collapse_ls_static_liq(solutions_by_layer, surface_boundary_condition)
            else:
                tidal_y = collapse_ls_dynamic_liq(
                    solutions_by_layer, surface_boundary_condition,
                    liquid_gravity, liquid_density, liquid_radii, frequency
                    )
        elif model_name == 'solid_liquid_solid':
            if is_liquid_static:
                tidal_y = collapse_sls_static_liq(
                    solutions_by_layer, surface_boundary_condition,
                    liquid_gravity, liquid_density
                    )
            else:
                tidal_y = collapse_sls_dynamic_liq(
                    solutions_by_layer, surface_boundary_condition,
                    liquid_gravity, liquid_density, liquid_radii, frequency
                    )
        elif model_name == 'solid_solid_liquid_solid':
            if is_liquid_static:
                tidal_y = collapse_ssls_static_liq(
                    solutions_by_layer, surface_boundary_condition,
                    liquid_gravity, liquid_density
                    )
            else:
                tidal_y = collapse_ssls_dynamic_liq(
                    solutions_by_layer, surface_boundary_condition,
                    liquid_gravity, liquid_density, liquid_radii, frequency
                    )
        else:
            raise Exception('Unknown Model.')

    if verbose:
        print('Done!')

    if nondimensionalize:
        if verbose:
            print('Re-dimensionalizing Radial Functions.')
        tidal_y = re_dimensionalize_radial_func(tidal_y, planet_radius, planet_bulk_density)

    return tidal_y
