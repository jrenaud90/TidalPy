""" Propagation of tidal solution using the fundamental matrix

References
----------
SVC16 : Sabadini, Vermeerson, & Cambiotti (2016, DOI: 10.1007/978-94-017-7552-6)
HH14  : Henning & Hurford (2014, DOI: 10.1088/0004-637X/789/1/30)
ID    : IcyDwarf Code by Marc Neveu (https://github.com/MarcNeveu/IcyDwarf/blob/master/IcyDwarf/Thermal.h)
B13   : Beuthe (2013, DOI: 10.1016/j.icarus.2012.11.020)
"""

import numpy as np

from ....utilities.performance import njit


# @njit(cacheable=True)
def propagate(fundamental_matrix: np.ndarray, fundamental_matrix_inverse: np.ndarray, derivative_matrix: np.ndarray,
              inner_boundary_condition: np.ndarray, world_radius: float,
              order_l: int = 2):
    """ This function will propagate the incompressible tidal equations, via the fundamental matrix, through a world or
    layers sub-shells.

    Assumptions:
    - Each shell is homogeneous (However, it can have different properties from the shells above/below it).

    References
    ----------
    See appendix of HH14 and lines ~2355 of thermal.h in ID

    Parameters
    ----------
    fundamental_matrix : np.ndarray
        Fundamental (incompressible) Matrix (6 x 6 x N); See fundamental.py
    fundamental_matrix_inverse : np.ndarray
        Inverse of the fundamental (incompressible) Matrix (6 x 6 x N); See fundamental.py
    derivative_matrix : np.ndarray
        Derivative matrix, A, defined by the function dy/dr = A dot y
    inner_boundary_condition : np.ndarray
        Boundary condition Matrix (6 x 3) of the tidal problem at the inner surface.
    world_radius : float
        Radius of the world [m]
    order_l : int = 2
        Tidal harmonic order

    Returns
    -------
    tidal_y : np.ndarray
        Matrix [6 x N] of tidal solutions. See decompression.py on how useful information is extracted.
    tidal_y_derivative : np.ndarray
        Matrix [6 x N] of the derivative of the tidal solutions with respect to radius.

    """
    num_shells = fundamental_matrix.shape[2]

    # Find propagation matrix
    #    Prop matrix has 6 rows but only 3 columns.
    propagation_mtx = np.zeros((6, 3, num_shells), dtype=np.complex128)
    #    It is initialized with the inner boundary condition matrix.
    propagation_mtx[:, :, 0] = inner_boundary_condition

    # Step through the planet's shells
    for i in range(num_shells):
        if i == 0:
            # Skip central shell - boundary conditions provided.
            continue

        # Dot product of between the fundamental matrix and (the dot product between the inverse fundamental matrix one
        #    layer below and the propagation matrix one layer below).
        propagation_mtx[:, :, i] = \
            fundamental_matrix[:, :, i] @ (fundamental_matrix_inverse[:, :, i-1] @ propagation_mtx[:, :, i-1])

    # Surface condition matrix is a 3x3 matrix of the top-most shell of the aggregate matrix's rows [3, 4, 6]
    surface_matrix = np.vstack( (propagation_mtx[2, :, -1], propagation_mtx[3, :, -1], propagation_mtx[5, :, -1]) )

    # The surface boundary conditions are (always?) static and only based on the radius and order-l
    surface_bc = np.zeros((3,), dtype=np.complex128)
    surface_bc[2] = -1. * (2. * order_l + 1.) / world_radius

    # Invert the surface matrix and solve using the surface boundary condition
    surface_matrix_inv = np.linalg.inv(surface_matrix)
    surface_solution = surface_matrix_inv @ surface_bc

    # Using the aggregate matrix, solve for the tidal "y"s
    tidal_y_sv = np.zeros((6, num_shells), dtype=np.complex128)
    tidal_y_derivative_sv = np.zeros((6, num_shells), dtype=np.complex128)
    for i in range(num_shells):
        tidal_y_sv[:, i] = propagation_mtx[:, :, i] @ surface_solution

        # Calculate the derivatives of the tidal solution with radius
        tidal_y_derivative_sv[:, i] = derivative_matrix[:, :, i] @ tidal_y_sv[:, i]

    # As discussed in B13 (discussed near their equation 7), SVC16 (and the earlier 2004 book) use a different
    #    convention for tidal_y than is used by Takeuchi and Saito (1972). Since a good chunk of the field follows the
    #    latter, we will do the same. Below are the conversions from SVC16 to TS72
    tidal_y = np.zeros_like(tidal_y_sv)
    tidal_y[0, :] = tidal_y_sv[0, :]
    tidal_y[1, :] = tidal_y_sv[2, :]  # Flip y3 and y2
    tidal_y[2, :] = tidal_y_sv[1, :]  # Flip y3 and y2
    tidal_y[3, :] = tidal_y_sv[3, :]
    tidal_y[4, :] = tidal_y_sv[4, :] * -1.
    tidal_y[5, :] = tidal_y_sv[5, :] * -1.

    # Likewise, take convert the derivatives to match the TS72 format
    tidal_y_derivative = np.zeros_like(tidal_y_derivative_sv)
    tidal_y_derivative[0, :] = tidal_y_derivative_sv[0, :]
    tidal_y_derivative[1, :] = tidal_y_derivative_sv[2, :]
    tidal_y_derivative[2, :] = tidal_y_derivative_sv[1, :]
    tidal_y_derivative[3, :] = tidal_y_derivative_sv[3, :]
    tidal_y_derivative[4, :] = tidal_y_derivative_sv[4, :] * -1.
    tidal_y_derivative[5, :] = tidal_y_derivative_sv[5, :] * -1.

    return tidal_y, tidal_y_derivative
