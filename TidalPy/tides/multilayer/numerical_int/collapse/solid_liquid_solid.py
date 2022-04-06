from typing import Tuple

import numpy as np

from .....utilities.performance import njit

TidalYSolType = Tuple[Tuple[np.ndarray, np.ndarray, np.ndarray],
                      Tuple[np.ndarray, ...],
                      Tuple[np.ndarray, np.ndarray, np.ndarray]]


@njit(cacheable=True)
def collapse_sls_static_liq(
    tidal_y_solutions_by_layer: TidalYSolType, surface_solution: np.ndarray,
    liquid_gravity_array: np.ndarray, liquid_density_array: np.ndarray
    ) -> np.ndarray:
    """ Determine the radial solution convergence for a planet with solid-liquid-solid structure.
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
    liquid_density_interface_1 = liquid_density_array[0]
    gravity_interface_1 = liquid_gravity_array[0]

    # Build solution matrix at surface
    sol_surf_mtx = np.asarray(
        [
            [tidal_y_layer2[0][1, -1], tidal_y_layer2[1][1, -1], tidal_y_layer2[2][1, -1]],
            [tidal_y_layer2[0][3, -1], tidal_y_layer2[1][3, -1], tidal_y_layer2[2][3, -1]],
            [tidal_y_layer2[0][5, -1], tidal_y_layer2[1][5, -1], tidal_y_layer2[2][5, -1]]
            ]
        )
    sol_surf_mtx_inv = np.linalg.inv(sol_surf_mtx)
    C_layer2_vector = sol_surf_mtx_inv @ surface_solution

    # Solve for the outer core Qs
    C_layer1_vector = np.empty(1, dtype=np.complex128)
    C_layer1_vector[0] = C_layer2_vector[0]

    # Solve for inner core Qs
    C_layer0_vector = np.empty(3, dtype=np.complex128)

    y4_frac_1 = tidal_y_layer0[0][3, -1] / tidal_y_layer0[2][3, -1]
    y4_frac_2 = tidal_y_layer0[1][3, -1] / tidal_y_layer0[2][3, -1]

    # gamma_j = (y_2j - f_j y_23) - rho( g(y_1j - f_j y_13) - (y_5j - f_j y_53))
    gamma_1 = (tidal_y_layer0[0][1, -1] - y4_frac_1 * tidal_y_layer0[2][1, -1]) - \
              liquid_density_interface_1 * (gravity_interface_1 * (
            tidal_y_layer0[0][0, -1] - y4_frac_1 * tidal_y_layer0[2][0, -1]) -
                                            (tidal_y_layer0[0][4, -1] - y4_frac_1 * tidal_y_layer0[2][4, -1]))
    gamma_2 = (tidal_y_layer0[1][1, -1] - y4_frac_2 * tidal_y_layer0[2][1, -1]) - \
              liquid_density_interface_1 * (gravity_interface_1 * (
            tidal_y_layer0[1][0, -1] - y4_frac_2 * tidal_y_layer0[2][0, -1]) -
                                            (tidal_y_layer0[1][4, -1] - y4_frac_2 * tidal_y_layer0[2][4, -1]))

    # Find inner core constants
    C_layer0_vector[0] = C_layer1_vector[0]
    C_layer0_vector[1] = (-gamma_1 / gamma_2) * C_layer0_vector[0]
    C_layer0_vector[2] = -y4_frac_1 * C_layer0_vector[0] - y4_frac_2 * C_layer0_vector[1]

    # Solve for the liquid layer's y's
    tidal_y_layer1 = C_layer1_vector[0] * tidal_y_layer1[0]

    shape = tidal_y_layer1[0, :].shape
    layer1_ys = (
        np.full(shape, np.nan, dtype=np.complex128),
        np.full(shape, np.nan, dtype=np.complex128),
        np.full(shape, np.nan, dtype=np.complex128),
        np.full(shape, np.nan, dtype=np.complex128),
        tidal_y_layer1[0, :],
        np.full(shape, np.nan, dtype=np.complex128),
        )

    tidal_y_layer1_full = np.vstack(layer1_ys)

    # Solve for total planet y's
    tidal_y_layer0 = C_layer0_vector[0] * tidal_y_layer0[0] + C_layer0_vector[1] * tidal_y_layer0[1] + \
                     C_layer0_vector[2] * tidal_y_layer0[2]

    tidal_y_layer2 = C_layer2_vector[0] * tidal_y_layer2[0] + C_layer2_vector[1] * tidal_y_layer2[1] + \
                     C_layer2_vector[2] * tidal_y_layer2[2]

    # Combine solutions for all layers
    tidal_y = np.concatenate((tidal_y_layer0, tidal_y_layer1_full, tidal_y_layer2), axis=1)

    return tidal_y


@njit(cacheable=True)
def collapse_sls_dynamic_liq(
    tidal_y_solutions_by_layer: TidalYSolType, surface_solution: np.ndarray,
    liquid_gravity_array: np.ndarray, liquid_density_array: np.ndarray,
    liquid_radius_array: np.ndarray, frequency: float
    ) -> np.ndarray:
    """ Determine the radial solution convergence for a planet with solid-liquid-solid structure.
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
    liquid_radius_array : np.ndarray
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

    # Build solution matrix at surface
    sol_surf_mtx = np.asarray(
        [
            [tidal_y_layer2[0][1, -1], tidal_y_layer2[1][1, -1], tidal_y_layer2[2][1, -1]],
            [tidal_y_layer2[0][3, -1], tidal_y_layer2[1][3, -1], tidal_y_layer2[2][3, -1]],
            [tidal_y_layer2[0][5, -1], tidal_y_layer2[1][5, -1], tidal_y_layer2[2][5, -1]]
            ]
        )
    sol_surf_mtx_inv = np.linalg.inv(sol_surf_mtx)
    C_layer2_vector = sol_surf_mtx_inv @ surface_solution

    # Solve for the outer core Qs
    C_layer1_vector = np.empty(2, dtype=np.complex128)
    C_layer1_vector[0] = C_layer2_vector[0]
    C_layer1_vector[1] = C_layer2_vector[1]

    # Solve for inner core Qs
    C_layer0_vector = np.empty(3, dtype=np.complex128)
    C_layer0_vector[0] = C_layer1_vector[0]
    C_layer0_vector[1] = C_layer1_vector[1]
    y4_frac_1 = tidal_y_layer0[0][3, -1] / tidal_y_layer0[2][3, -1]
    y4_frac_2 = tidal_y_layer0[1][3, -1] / tidal_y_layer0[2][3, -1]
    C_layer0_vector[2] = -y4_frac_1 * C_layer0_vector[0] - y4_frac_2 * C_layer0_vector[1]

    # Solve for the liquid layer's y's
    tidal_y_layer1 = C_layer1_vector[0] * tidal_y_layer1[0] + C_layer1_vector[1] * tidal_y_layer1[1]

    # Outer core is missing two y's, fix that now.
    y3_layer1 = \
        (1. / (frequency**2 * liquid_density_array * liquid_radius_array)) * \
        (liquid_density_array * liquid_gravity_array * tidal_y_layer1[0, :] -
         tidal_y_layer1[1, :] - liquid_density_array * tidal_y_layer1[2, :])

    shape = tidal_y_layer1[0, :].shape
    layer1_ys = (
        tidal_y_layer1[0, :],
        tidal_y_layer1[1, :],
        y3_layer1,
        np.full(shape, np.nan, dtype=np.complex128),
        tidal_y_layer1[2, :],
        tidal_y_layer1[3, :]
        )

    tidal_y_layer1_full = np.vstack(layer1_ys)

    # Solve for total planet y's
    tidal_y_layer0 = C_layer0_vector[0] * tidal_y_layer0[0] + C_layer0_vector[1] * tidal_y_layer0[1] + \
                     C_layer0_vector[2] * tidal_y_layer0[2]

    tidal_y_layer2 = C_layer2_vector[0] * tidal_y_layer2[0] + C_layer2_vector[1] * tidal_y_layer2[1] + \
                     C_layer2_vector[2] * tidal_y_layer2[2]

    # Combine solutions for all layers
    tidal_y = np.concatenate((tidal_y_layer0, tidal_y_layer1_full, tidal_y_layer2), axis=1)

    return tidal_y
