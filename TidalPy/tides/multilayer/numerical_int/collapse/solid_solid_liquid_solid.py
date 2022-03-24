from typing import Tuple

import numpy as np

from .....utilities.performance import njit

TidalYSolType = Tuple[Tuple[np.ndarray, np.ndarray, np.ndarray],
                      Tuple[np.ndarray, ...],
                      Tuple[np.ndarray, np.ndarray, np.ndarray]]


@njit(cacheable=True)
def collapse_ssls_static_liq(
    tidal_y_solutions_by_layer: TidalYSolType, surface_solution: np.ndarray,
    liquid_gravity_array: np.ndarray, liquid_density_array: np.ndarray
    ) -> np.ndarray:
    """ Determine the radial solution convergence for a planet with solid-solid-liquid-solid structure.
    A static liquid layer is assumed.

    Used in the numerical shooting method.

    Parameters
    ----------
    tidal_y_solutions_by_layer : TidalYSolType
        Radial functions solved via integration separated by layer.
    surface_solution : np.ndarray
        Surface boundary condition used to find the constants at the top-most layer.
    liquid_gravity_array : np.ndarray
        Acceleration due to gravity within the liquid layer [m s-2]
    liquid_density_array : np.ndarray
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
    liquid_density_interface_2 = liquid_density_array[0]
    gravity_interface_2 = liquid_gravity_array[0]

    # Build solution matrix at surface
    sol_surf_mtx = np.asarray(
        [
            [tidal_y_layer3[0][1, -1], tidal_y_layer3[1][1, -1], tidal_y_layer3[2][1, -1]],
            [tidal_y_layer3[0][3, -1], tidal_y_layer3[1][3, -1], tidal_y_layer3[2][3, -1]],
            [tidal_y_layer3[0][5, -1], tidal_y_layer3[1][5, -1], tidal_y_layer3[2][5, -1]]
            ]
        )
    sol_surf_mtx_inv = np.linalg.inv(sol_surf_mtx)
    C_layer3_vector = sol_surf_mtx_inv @ surface_solution

    # Solve for the liquid layer's constants
    C_layer2_vector = np.empty(1, dtype=np.complex128)
    C_layer2_vector[0] = C_layer3_vector[0]

    # Solve for layer 1's constants
    C_layer1_vector = np.empty(3, dtype=np.complex128)

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
    C_layer0_vector = np.empty(3, dtype=np.complex128)
    C_layer0_vector[0] = C_layer1_vector[0]
    C_layer0_vector[1] = C_layer1_vector[1]
    C_layer0_vector[2] = C_layer1_vector[2]

    # Solve for the liquid layer's y's
    tidal_y_layer2 = C_layer2_vector[0] * tidal_y_layer2[0]

    shape = tidal_y_layer2[0, :].shape
    layer2_ys = (
        np.full(shape, np.nan, dtype=np.complex128),
        np.full(shape, np.nan, dtype=np.complex128),
        np.full(shape, np.nan, dtype=np.complex128),
        np.full(shape, np.nan, dtype=np.complex128),
        tidal_y_layer2[0, :],
        np.full(shape, np.nan, dtype=np.complex128),
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
def collapse_ssls_dynamic_liq(
    tidal_y_solutions_by_layer: TidalYSolType, surface_solution: np.ndarray,
    liquid_gravity_array: np.ndarray, liquid_density_array: np.ndarray,
    liquid_radii: np.ndarray, frequency: float
    ) -> np.ndarray:
    """ Determine the radial solution convergence for a planet with solid-solid-liquid-solid structure.
    A dynamic liquid layer is assumed.

    Used in the numerical shooting method.

    Parameters
    ----------
    tidal_y_solutions_by_layer : TidalYSolType
        Radial functions solved via integration separated by layer.
    surface_solution : np.ndarray
        Surface boundary condition used to find the constants at the top-most layer.
    liquid_gravity_array : np.ndarray
        Acceleration due to gravity within the liquid layer [m s-2]
    liquid_density_array : np.ndarray
        Density within the liquid layer [kg m-3]
    liquid_radii : np.ndarray
        Radius array for the liquid layer [m]
    frequency : float
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
    sol_surf_mtx = np.asarray(
        [
            [tidal_y_layer3[0][1, -1], tidal_y_layer3[1][1, -1], tidal_y_layer3[2][1, -1]],
            [tidal_y_layer3[0][3, -1], tidal_y_layer3[1][3, -1], tidal_y_layer3[2][3, -1]],
            [tidal_y_layer3[0][5, -1], tidal_y_layer3[1][5, -1], tidal_y_layer3[2][5, -1]]
            ]
        )
    sol_surf_mtx_inv = np.linalg.inv(sol_surf_mtx)
    C_layer3_vector = sol_surf_mtx_inv @ surface_solution

    # Solve for the liquid layer constants
    C_layer2_vector = np.empty(2, dtype=np.complex128)
    C_layer2_vector[0] = C_layer3_vector[0]
    C_layer2_vector[1] = C_layer3_vector[1]

    # Solve for the layer 1 constants
    C_layer1_vector = np.empty(3, dtype=np.complex128)
    C_layer1_vector[0] = C_layer2_vector[0]
    C_layer1_vector[1] = C_layer2_vector[1]
    y4_frac_1_layer_1 = tidal_y_layer1[0][3, -1] / tidal_y_layer1[2][3, -1]
    y4_frac_2_layer_1 = tidal_y_layer1[1][3, -1] / tidal_y_layer1[2][3, -1]
    C_layer1_vector[2] = -y4_frac_1_layer_1 * C_layer1_vector[0] - y4_frac_2_layer_1 * C_layer1_vector[1]

    # Solve for the innermost layers constants
    C_layer0_vector = np.empty(3, dtype=np.complex128)
    C_layer0_vector[0] = C_layer1_vector[0]
    C_layer0_vector[1] = C_layer1_vector[1]
    C_layer0_vector[2] = C_layer1_vector[1]

    # Solve for the liquid layer's y's
    tidal_y_layer2 = C_layer2_vector[0] * tidal_y_layer2[0] + C_layer2_vector[1] * tidal_y_layer2[1]

    # Liquid layer is missing two y's, fix that now.
    y3_layer2 = \
        (1. / (frequency**2 * liquid_density_array * liquid_radii)) * \
        (liquid_radii * liquid_gravity_array * tidal_y_layer2[0, :] -
         tidal_y_layer2[1, :] - liquid_density_array * tidal_y_layer2[2, :])

    shape = tidal_y_layer2[0, :].shape
    layer2_ys = (
        tidal_y_layer2[0, :],
        tidal_y_layer2[1, :],
        y3_layer2,
        np.full(shape, np.nan, dtype=np.complex128),
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
