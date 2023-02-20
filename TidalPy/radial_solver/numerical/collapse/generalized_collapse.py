""" Contains generalized collapse function to find the final solution for the radial functions y within a planet.

References
----------
TS72  : Takeuchi & Saito 1972 (Seismology; DOI: 10.1016/B978-0-12-460811-5.50010-6)
S74   : Saito 1974 (JPE; DOI: 10.4294/jpe1952.22.123)
KMN15 : Kamata et al 2015 (JGR; DOI: 10.1002/2015JE004821)
KTC21 : Kervazo et al 2021 (A&A; DOI: 10.1051/0004-6361/202039433)

"""

from typing import Tuple, Union, List

import numpy as np

from TidalPy.constants import G
from TidalPy.utilities.performance import njit, nbList

from .surface_condition import solid_surface, static_liquid_surface, dynamic_liquid_surface


@njit(cacheable=True)
def collapse_solutions(
        y_solutions_by_layer: Tuple[Tuple[np.ndarray, np.ndarray, np.ndarray]],
        is_solid_by_layer: Union[List[bool], Tuple[bool, ...]],
        is_static_by_layer: Union[List[bool], Tuple[bool, ...]],
        indices_by_layer: Union[List[np.ndarray], Tuple[np.ndarray, ...]],
        surface_boundary_condition: np.ndarray,
        radius_array: np.ndarray, density_array: np.ndarray, gravity_array: np.ndarray,
        gravity_at_interfaces: Union[List[float], Tuple[float, ...]],
        liquid_density_at_interfaces: Union[List[float], Tuple[float, ...]],
        gravity_at_surface: float,
        frequency: float,
        G_to_use: float = G
        ):
    """ Function to solve for unknown constants of integration in order to find the final viscoelastic-gravitational
     solution throughout a planet.

    Used in the numerical shooting method.

    References
    ----------
    TS72  : Takeuchi & Saito 1972 (Seismology; DOI: 10.1016/B978-0-12-460811-5.50010-6)
    S74   : Saito 1974 (JPE; DOI: 10.4294/jpe1952.22.123)
    KMN15 : Kamata et al 2015 (JGR; DOI: 10.1002/2015JE004821)
    KTC21 : Kervazo et al 2021 (A&A; DOI: 10.1051/0004-6361/202039433)

    Parameters
    ----------
    y_solutions_by_layer : Tuple[np.ndarray, np.ndarray, np.ndarray]
        Radial functions solved via integration for an entire homogeneous planet.
    is_solid_by_layer : Union[List[bool], Tuple[bool, ...]]
        Flags for if each layer is solid (True) or liquid (False)
    is_static_by_layer : Union[List[bool], Tuple[bool, ...]]
        Flags for if each layer is static (True) or dynamic (False)
    indices_by_layer : Union[List[np.ndarray], Tuple[np.ndarray, ...]]
    surface_boundary_condition : np.ndarray
        Surface boundary condition used to find the constants at the top-most layer.
    radius_array : np.ndarray
        Numpy array of radii throughout the planet [m]
    density_array : np.ndarray
        Numpy array of density throughout the planet [kg m-3]
    gravity_array : np.ndarray
        Numpy array of acceleration due to gravity throughout the planet [m s-2]
    gravity_at_interfaces : Union[List[float], Tuple[float, ...]]
        Acceleration due to gravity at each layer's interface [m s-2]
    liquid_density_at_interfaces : Union[List[float], Tuple[float, ...]]
        Density of the static liquid layer at each interface (0 if interface of two solids) [kg m-3]
    gravity_at_surface : float
        The acceleration due to gravity at the planet's surface [m s-2]
    frequency : float
        Forcing frequency [rad s-1]
    G_to_use : float = G
        Gravitational constant. Can be provided by user if other functions have had their values nondimensionalized.

    Returns
    -------
    total_y : np.ndarray
        Collapsed radial functions for the entire planet. Scaled by the correct constants.
        This will be a 6 x N ndarray for the six radial functions, where N is the number of radial shells.

    """
    # Determine number of layers
    num_layers = len(y_solutions_by_layer)

    # Make empty variables that will be set in the loop
    layer_above_y = y_solutions_by_layer[0][0]  # Pull out the first layer's first solution.
    radius_len = indices_by_layer[0].size
    layer_above_is_solid = True
    layer_above_is_static = True
    layer_above_constants = np.zeros(3, dtype=np.complex128)

    # Create storage for final results. Put a fake result into it to start out with so that numba can compile.
    total_y = np.empty((6, radius_len), dtype=np.complex128)

    # Work from the surface to the core
    for layer_i_bot in range(num_layers):
        # Convert so we are working from surface down
        layer_i = num_layers - (layer_i_bot + 1)

        # Pull out some layer information
        is_solid = is_solid_by_layer[layer_i]
        is_static = is_static_by_layer[layer_i]
        y_solutions = y_solutions_by_layer[layer_i]
        layer_index = indices_by_layer[layer_i]
        layer_density = density_array[layer_index]
        layer_gravity = gravity_array[layer_index]
        layer_radius = radius_array[layer_index]
        y_surface_solutions = [y_sol[:, -1] for y_sol in y_solutions]
        num_solutions = len(y_surface_solutions)

        # Find y's in layer based on the multiple solutions and global constants found via boundary conditions.
        if layer_i_bot == 0:
            # Working on surface layer
            # Determine the uppermost layer's solution based on boundary conditions
            if is_solid:
                # Surface is solid
                constant_vector = solid_surface(y_surface_solutions, surface_boundary_condition)
            elif is_static:
                # Surface is a static liquid layer
                constant_vector = static_liquid_surface(
                        y_surface_solutions, surface_boundary_condition,
                        gravity_at_surface, G_to_use)
            else:
                # Surface is a dynamic liquid layer
                constant_vector = dynamic_liquid_surface(y_surface_solutions, surface_boundary_condition)
        else:
            # Working on interior layers. Will need to find the constants of integration based on the layer above.

            # Interfaces are defined at the bottom of the layer in question. However, this function is calculating
            #  the transition at the top of each layer as it works its way down.
            #  So, for interface values, we actually need the ones of the layer above us.
            gravity_at_interface = gravity_at_interfaces[layer_i + 1]
            liquid_density_at_interface = liquid_density_at_interfaces[layer_i + 1]

            if is_solid:
                if layer_above_is_solid:
                    # Both layers are solid. Constants are the same.
                    constant_vector = np.copy(layer_above_constants)
                else:
                    constant_vector = np.empty(3, dtype=np.complex128)

                    # Create some helper functions that will be needed
                    y4_frac_1 = -y_surface_solutions[0][3] / y_surface_solutions[2][3]
                    y4_frac_2 = -y_surface_solutions[1][3] / y_surface_solutions[2][3]

                    if layer_above_is_static:
                        # Need to find 3 solid constants from 1 liquid constant
                        # S74, Page 131
                        constant_vector[0] = layer_above_constants[0]
                        # Derived by JPR based on Eq 21 (2nd line) of S74
                        gamma_1 = \
                            (y_surface_solutions[0][1] + y4_frac_1 * y_surface_solutions[2][1]) - \
                            (liquid_density_at_interface *
                             (gravity_at_interface * (
                                     y_surface_solutions[0][0] + y4_frac_1 * y_surface_solutions[2][0]) -
                              (y_surface_solutions[0][4] + y4_frac_1 * y_surface_solutions[2][4])
                              )
                             )
                        gamma_2 = \
                            (y_surface_solutions[1][1] + y4_frac_2 * y_surface_solutions[2][1]) - \
                            (liquid_density_at_interface *
                             (gravity_at_interface * (
                                     y_surface_solutions[1][0] + y4_frac_2 * y_surface_solutions[2][0]) -
                              (y_surface_solutions[1][4] + y4_frac_2 * y_surface_solutions[2][4])
                              )
                             )

                        constant_vector[1] = (-gamma_1 / gamma_2) * constant_vector[0]
                        # TS72, Eq. 142 (utilizes y_4 = 0)
                        constant_vector[2] = y4_frac_1 * constant_vector[0] + y4_frac_2 * constant_vector[1]

                    else:
                        # Need to find 3 solid constants from 2 liquid constants
                        # TS72, Eq. 144
                        constant_vector[0] = layer_above_constants[0]
                        constant_vector[1] = layer_above_constants[1]
                        # TS72, Eq. 142 (utilizes y_4 = 0)
                        constant_vector[2] = y4_frac_1 * constant_vector[0] + y4_frac_2 * constant_vector[1]
            else:
                if is_static:
                    if not layer_above_is_solid:
                        # Liquid layer above
                        if layer_above_is_static:
                            # Both layers are static liquids. Constants are the same.
                            constant_vector = np.copy(layer_above_constants)
                        else:
                            # Dynamic liquid above
                            # JPR decided to follow a similar approach as Eq. 20 in S74:
                            #   Treat the lower static liquid as normal.
                            #   The upper dynamic liquid layer is treated like the solid layer in Eq. 20 except
                            #    that y_3 is undefined as is "set 3" solution mentioned in that text.
                            constant_vector = np.empty(1, dtype=np.complex128)
                            constant_vector[0] = layer_above_constants[0]
                    else:
                        # Solid layer above
                        # Based on S74. The constant in this layer is just equal to the constant in solution 1 of the
                        #  layer above.
                        constant_vector = np.empty(1, dtype=np.complex128)
                        constant_vector[0] = layer_above_constants[0]
                else:
                    if not layer_above_is_solid:
                        # Liquid layer above
                        if layer_above_is_static:
                            constant_vector = np.empty(2, dtype=np.complex128)
                            # Need to find 2 liquid (dynamic) constants from 1 liquid (static) constant
                            # S74, Page 131
                            constant_vector[0] = layer_above_constants[0]
                            # Derived by JPR based on Eq 21 (2nd line) of S74
                            # Pull out ys
                            # # Solution 1
                            lower_s1y1 = y_surface_solutions[0][0]
                            lower_s1y2 = y_surface_solutions[0][1]
                            lower_s1y5 = y_surface_solutions[0][2]
                            lower_s1y6 = y_surface_solutions[0][3]
                            # # Solution 2
                            lower_s2y1 = y_surface_solutions[1][0]
                            lower_s2y2 = y_surface_solutions[1][1]
                            lower_s2y5 = y_surface_solutions[1][2]
                            lower_s2y6 = y_surface_solutions[1][3]
                            # lambda_j = (y_2j - rho * ( g * y_1j - y_5j))
                            lambda_1 = lower_s1y2 - liquid_density_at_interface * \
                                       (gravity_at_interface * lower_s1y1 - lower_s1y5)
                            lambda_2 = lower_s2y2 - liquid_density_at_interface * \
                                       (gravity_at_interface * lower_s2y1 - lower_s2y5)
                            constant_vector[1] = (-lambda_1 / lambda_2) * constant_vector[0]
                        else:
                            # Both layers are dynamic liquids. Constants are the same.
                            constant_vector = np.copy(layer_above_constants)
                    else:
                        # Solid layer above
                        # TS72 Eqs. 148-149
                        constant_vector = np.empty(2, dtype=np.complex128)
                        constant_vector[0] = layer_above_constants[0]
                        constant_vector[1] = layer_above_constants[1]

        # Solve for y in layer
        y_at_layer = constant_vector[0] * y_solutions[0]
        if num_solutions > 1:
            for sol_i in range(num_solutions):
                if sol_i == 0:
                    # Already started with this value so skip it now.
                    continue
                y_at_layer += constant_vector[sol_i] * y_solutions[sol_i]

        # Convert y solution into correct shape & solve for any dependent values.
        if is_solid:
            # All radial functions for a solid layer are already present and in the correct order.
            y_at_layer_expanded = y_at_layer
        else:
            # Get layer's radial shape to build nan arrays.
            shape = y_at_layer[0, :].shape
            if is_static:
                # A static liquid layer's y_1-y_4 & y_6 are undefined. We could pass y_7 in place of y_6 but that may
                #  cause confusion elsewhere.
                layer_y_tuple = (
                    np.full(shape, np.nan, dtype=np.complex128),
                    np.full(shape, np.nan, dtype=np.complex128),
                    np.full(shape, np.nan, dtype=np.complex128),
                    np.full(shape, np.nan, dtype=np.complex128),
                    y_at_layer[0, :],
                    np.full(shape, np.nan, dtype=np.complex128)
                    )
                y_at_layer_expanded = np.vstack(layer_y_tuple)
            else:
                # Dynamic liquid layer's y3 is not uniquely determined through integration, but is finite and can be
                #  determined using the other radial functions.
                # A dynamic liquid layer will be missing two y's, fix that now.
                y3_dynamic_liq = \
                    (1. / (frequency**2 * layer_density * layer_radius)) * \
                    (layer_density * layer_gravity * y_at_layer[0, :] -
                     y_at_layer[1, :] - layer_density * y_at_layer[2, :])

                # A dynamic liquid layer's y_4 is undefined. Other y's need to shift around to make the final y-matrix
                layer_y_tuple = (
                    y_at_layer[0, :],
                    y_at_layer[1, :],
                    y3_dynamic_liq,
                    np.full(shape, np.nan, dtype=np.complex128),
                    y_at_layer[2, :],
                    y_at_layer[3, :]
                    )
                y_at_layer_expanded = np.vstack(layer_y_tuple)

        # Store result and reference for next layer's calculation
        total_y[:, layer_index] = y_at_layer_expanded
        layer_above_y = y_at_layer
        layer_above_constants = constant_vector
        layer_above_is_solid = is_solid
        layer_above_is_static = is_static

    return total_y
