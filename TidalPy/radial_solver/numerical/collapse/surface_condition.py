""" Helper functions to solve for the radial solution surface boundary conditions at a world's surface.

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
from TidalPy.utilities.performance import njit


@njit(cacheable=True)
def solid_surface(
        y_solutions_at_surface: Union[List[np.ndarray], Tuple[np.ndarray, np.ndarray, np.ndarray]],
        surface_boundary_condition: np.ndarray
        ) -> np.ndarray:
    """ Calculate unknown radial function constants for a solid surface layer.

    Methods
    -------
    At the surface,
    y_2 = S_1
    y_4 = S_4
    y_6 = S_6

    y_2 = \sum_{i=0}^{3} y_{2,i} C_{i}
    y_4 = \sum_{i=0}^{3} y_{4,i} C_{i}
    y_6 = \sum_{i=0}^{3} y_{6,i} C_{i}

    To find the unknown C's, must solve the linear system of equations
    \matrix{y} \vector{C} = \vector{S}

    Where S are surface solutions that vary based on if the problem is tidal, loading, etc.

    References
    ----------
    Eq. B.37 in KTC21
    Eq. 16 in KMN15

    Parameters
    ----------
    y_solutions_at_surface : Union[List[np.ndarray], Tuple[np.ndarray, np.ndarray, np.ndarray]]
        Radial functions solved via integration for an entire homogeneous planet. Limited to the uppermost layer's
        solution at the very top.
    surface_boundary_condition : np.ndarray
        Surface boundary condition used to find the constants in the uppermost layer.

    Returns
    -------
    C_vector : np.ndarray
        Calculated constants for the uppermost layer's solution.

    """

    # Build y-solution matrix to be applied to the surface.
    sol_surf_mtx_solid = np.asarray(
            (
                (y_solutions_at_surface[0][1], y_solutions_at_surface[1][1], y_solutions_at_surface[2][1]),
                (y_solutions_at_surface[0][3], y_solutions_at_surface[1][3], y_solutions_at_surface[2][3]),
                (y_solutions_at_surface[0][5], y_solutions_at_surface[1][5], y_solutions_at_surface[2][5])
                )
            )

    # Solve for the upper layer's (solid) constants by inverting the matrix and applying to boundary conditions.
    sol_surf_mtx_inv = np.linalg.inv(sol_surf_mtx_solid)
    C_vector = sol_surf_mtx_inv @ surface_boundary_condition

    return C_vector


@njit(cacheable=True)
def dynamic_liquid_surface(
        y_solutions_at_surface: Union[List[np.ndarray], Tuple[np.ndarray, np.ndarray]],
        surface_boundary_condition: np.ndarray
        ) -> np.ndarray:
    """ Calculate unknown radial function constants for a dynamic liquid surface layer.

    Methods
    -------
    Unlike the solid layer, a liquid layer's y_4 is undefined. That leads to one less boundary condition and one
     less solution (2 total).

    At the surface,
    y_2 = S_1
    y_6 = S_6

    y_2 = \sum_{i=0}^{2} y_{2,i} C_{i}
    y_6 = \sum_{i=0}^{2} y_{6,i} C_{i}

    To find the unknown C's, must solve the linear system of equations
    \matrix{y} \vector{C} = \vector{S}

    Where S are surface solutions that vary based on if the problem is tidal, loading, etc.

    References
    ----------
    Eq. B.38 in KTC21
    Eq. 17 in KMN15

    Parameters
    ----------
    y_solutions_at_surface : Union[List[np.ndarray], Tuple[np.ndarray, np.ndarray]]
        Radial functions solved via integration for an entire homogeneous planet. Limited to the uppermost layer's
        solution at the very top.
    surface_boundary_condition : np.ndarray
        Surface boundary condition used to find the constants in the uppermost layer.

    Returns
    -------
    C_vector : np.ndarray
        Calculated constants for the uppermost layer's solution.

    """

    # The surface boundary condition will still have 3 members. Drop the one related to y_4
    surface_boundary_condition_to_use = np.empty(2, dtype=np.complex128)
    surface_boundary_condition_to_use[0] = surface_boundary_condition[0]  # y_2
    surface_boundary_condition_to_use[1] = surface_boundary_condition[2]  # y_6

    # Build y-solution matrix to be applied to the surface.
    # Note: for a dynamic liquid, y_2 and y_6 are held at indices 1 and 3 respectively
    sol_surf_mtx_dynamic_liquid = np.asarray(
            (
                (y_solutions_at_surface[0][1], y_solutions_at_surface[1][1]),
                (y_solutions_at_surface[0][3], y_solutions_at_surface[1][3])
                )
            )

    # Solve for the upper layer's (solid) constants by inverting the matrix and applying to boundary conditions.
    sol_surf_mtx_inv = np.linalg.inv(sol_surf_mtx_dynamic_liquid)
    C_vector = sol_surf_mtx_inv @ surface_boundary_condition_to_use

    return C_vector


@njit(cacheable=True)
def static_liquid_surface(
        y_solutions_at_surface: Union[List[np.ndarray], Tuple[np.ndarray]],
        surface_boundary_condition: np.ndarray,
        gravity_at_surface: float,
        G_to_use: float = G,
        ) -> np.ndarray:
    """ Calculate unknown radial function constants for a static liquid surface layer.

    Methods
    -------
    Unlike the dynamic liquid layer, a static liquid layer's y_2 is undefined. That leads to one less boundary condition
     and one less solution (1 total).

    At the surface,
    y_7 = S_7

    y_7 = y_{7,0} C_{0}

    To find the unknown C's, must solve the linear system of equations
    \matrix{y} \vector{C} = \vector{S}

    Where S are surface solutions that vary based on if the problem is tidal, loading, etc.

    References
    ----------
    Eq. 17, 10 in S74

    Parameters
    ----------
    y_solutions_at_surface : Union[List[np.ndarray], Tuple[np.ndarray]]
        Radial functions solved via integration for an entire homogeneous planet. Limited to the uppermost layer's
        solution at the very top.
    surface_boundary_condition : np.ndarray
        Surface boundary condition used to find the constants in the uppermost layer.
    gravity_at_surface : float
        The acceleration due to gravity at the planet's surface [m s-2]
    G_to_use : float = G
        Gravitational constant. Can be provided by user if other functions have had their values nondimensionalized.

    Returns
    -------
    C_vector : np.ndarray
        Calculated constants for the uppermost layer's solution.

    """

    # The surface boundary condition will still have 3 members. Drop the one related to y_4 and solve for the new one
    #   related to y_7. Found by the surface conditions for y_2 and y_6 and Eq. 17 from S74.
    surface_boundary_condition_to_use = np.empty(1, dtype=np.complex128)
    # y_7 = y_6 + (4 pi G / g) y_2
    surface_boundary_condition_to_use[0] = \
        surface_boundary_condition[2] + \
        surface_boundary_condition[0] * (4. * np.pi * G_to_use / gravity_at_surface)

    # Build y-solution matrix to be applied to the surface.
    # Note: for a static liquid layer, y_7 held in index 1 (index 0 is y_5).
    sol_surf_mtx_static_liquid = np.asarray(
            (
                (y_solutions_at_surface[0][1],),
                )
            )

    # Solve for the upper layer's (solid) constants by inverting the matrix and applying to boundary conditions.
    sol_surf_mtx_inv = np.linalg.inv(sol_surf_mtx_static_liquid)
    C_vector = sol_surf_mtx_inv @ surface_boundary_condition_to_use

    return C_vector
