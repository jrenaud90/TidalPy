""" Propagation of tidal solution using the fundamental matrix

References
----------
SVC16 : Sabadini, Vermeerson, & Cambiotti (2016, DOI: 10.1007/978-94-017-7552-6)
HH14  : Henning & Hurford (2014, DOI: 10.1088/0004-637X/789/1/30)
ID    : IcyDwarf Code by Marc Neveu (https://github.com/MarcNeveu/IcyDwarf/blob/master/IcyDwarf/Thermal.h)
"""

import numpy as np

from ...utilities.performance import njit


@njit(cacheable=True)
def propagate(fundamental_matrix: np.ndarray, fundamental_matrix_inverse: np.ndarray,
              central_boundary_condition: np.ndarray, world_radius: float,
              order_l: int = 2):
    """ This function will propagate the incompressible tidal equations, via the fundamental matrix, through a world's
        shells.

    References
    ----------
    See appendix of HH14 and lines ~2355 of thermal.h in ID

    Parameters
    ----------
    fundamental_matrix : np.ndarray
        Fundamental (incompressible) Matrix (6 x 6 x N); See fundamental.py
    fundamental_matrix_inverse : np.ndarray
        Inverse of the fundamental (incompressible) Matrix (6 x 6 x N); See fundamental.py
    central_boundary_condition : np.ndarray
        Boundary conditions (In (6 x 3) Matrix) of the tidal problem at the inner surface.
    world_radius : float
        Radius of the world [m]
    order_l : int = 2
        Tidal harmonic order

    Returns
    -------
    tidal_y : np.ndarray
        Matrix [6 x N] of tidal solutions. See decompression.py on how useful information is extracted.

    """
    num_shells = fundamental_matrix.shape[2]

    # Find propagation matrix
    propagation_mtx = np.zeros((6, 6, num_shells), dtype=np.complex128)
    for i in range(num_shells):
        if i == 0:
            # Skip central shell - boundary conditions provided.
            continue

        propagation_mtx[:, :, i] = fundamental_matrix[:, :, i] @ fundamental_matrix_inverse[:, :, i-1]

    # Create aggregate matrix based on central boundary conditions
    #    Aggregate matrix has 6 rows but only 3 columns.
    aggregate_matrix = np.zeros((6, 3, num_shells), dtype=np.complex128)
    #    It is initialized with the central boundary condition matrix.
    aggregate_matrix[:, :, 0] = central_boundary_condition

    # Step through the planet's shells
    for i in range(num_shells):
        if i == 0:
            # Skip central shell - boundary conditions provided.
            continue

        aggregate_matrix[:, :, i] = propagation_mtx[:, :, i] @ aggregate_matrix[:, :, i-1]

    # Surface condition matrix is a 3x3 matrix of the top-most shell of the aggregate matrix's rows [3, 4, 6]
    surface_matrix = np.vstack( (aggregate_matrix[2, :, -1], aggregate_matrix[3, :, -1], aggregate_matrix[5, :, -1]) )

    # The surface boundary conditions are (always?) static and only based on the radius and order-l
    surface_bc = np.zeros((3,), dtype=np.complex128)
    surface_bc[2] = -1. * (2. * order_l + 1) * world_radius

    # Invert the surface matrix and solve using the surface boundary condition
    surface_matrix_inv = np.linalg.inv(surface_matrix)
    surface_solution = surface_matrix_inv @ surface_bc

    # Using the aggregate matrix, solve for the tidal "y"s
    tidal_y = np.zeros((6, num_shells), dtype=np.complex128)
    for i in range(num_shells):
        tidal_y[:, i] = aggregate_matrix[:, :, i] @ surface_solution

    return tidal_y
