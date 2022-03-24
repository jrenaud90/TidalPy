from typing import Tuple

import numpy as np

from .....utilities.performance import njit


@njit(cacheable=True)
def collapse_homogen_solid(solid_solutions: Tuple[Tuple[np.ndarray, np.ndarray, np.ndarray]], surface_solution: np.ndarray):
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
    # The homogenous planet only has 1 layer but the solution is still stored in a list as if it had more than one.
    #   so we need to pull out this "layer".
    solid_solutions = solid_solutions[0]

    # Build solution matrix at surface
    sol_surf_mtx = np.asarray(
        [
            [solid_solutions[0][1, -1], solid_solutions[1][1, -1], solid_solutions[2][1, -1]],
            [solid_solutions[0][3, -1], solid_solutions[1][3, -1], solid_solutions[2][3, -1]],
            [solid_solutions[0][5, -1], solid_solutions[1][5, -1], solid_solutions[2][5, -1]]
            ]
        )
    sol_surf_mtx_inv = np.linalg.inv(sol_surf_mtx)
    C_vector = sol_surf_mtx_inv @ surface_solution

    # Solve for total planet y's
    tidal_y = C_vector[0] * solid_solutions[0] + C_vector[1] * solid_solutions[1] + \
              C_vector[2] * solid_solutions[2]

    return tidal_y
