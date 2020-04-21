from scipy.constants import G, pi
import numpy as np


#@njit
def build_fundamental_matrix(radius_array: np.ndarray, shear_array: np.ndarray,
                             density_array: np.ndarray, gravity_array: np.ndarray, order_l: int = 2):
    """ Construct the fundamental matrix and its inverse

    See Sabadini, Vermeerson, & Cambiotti (2016) (Hereafter cited as SVC)
    
    Parameters
    ----------
    radius_array : np.ndarray
    shear_array : np.ndarray
    density_array : np.ndarray
    gravity_array : np.ndarray
    order_l : int
        Tidal Harmonic Degree Index

    Returns
    -------

    """

    # This function manually defines the inverse matrix which is about 15--30% faster than a version that uses
    #    np.linalg.inv() to calculate the inverse matrix.

    num_shells = radius_array.shape[0]

    # First index: Rows of the propagation matrix
    # Second index: Columns of the propagation matrix
    # Third index: layer index (set outside of function)
    Y = np.zeros((6, 6, num_shells), dtype=np.complex)

    # Optimizations
    rgp = radius_array * gravity_array * density_array

    # Build Fundamental Matrix (zeros do not need to be specifically stated as they were put in at initialization)
    #     Eq. 2.42 in SVC
    ## Column 1
    Y[0, 0, :] = order_l * radius_array**(order_l + 1) / (2. * (2. * order_l + 3.))
    Y[1, 0, :] = (order_l + 3) * radius_array**(order_l + 1) / (2 * (2 * order_l + 3) * (order_l + 1))
    Y[2, 0, :] = (order_l * rgp + 2 * (order_l**2 - order_l - 3) * shear_array) * radius_array**order_l \
                 / (2 * (2 * order_l + 3))
    Y[3, 0, :] = order_l * (order_l + 2) * shear_array * radius_array**order_l\
                 / ((2 * order_l + 3) * (order_l + 1))
    # Y[4, 0, :] = np.zeros(shape)
    Y[5, 0, :] = 2 * pi * G * density_array * order_l * radius_array**(order_l + 1) / (2 * order_l + 3)

    ## Column 2
    Y[0, 1, :] = radius_array**(order_l - 1)
    Y[1, 1, :] = radius_array**(order_l - 1) / order_l
    Y[2, 1, :] = (rgp + 2 * (order_l - 1) * shear_array) * radius_array**(order_l - 2)
    Y[3, 1, :] = 2 * (order_l - 1) * shear_array * radius_array**(order_l - 2) / order_l
    # Y[4, 1, :] = np.zeros(shape)
    Y[5, 1, :] = 4 * pi * G * density_array * radius_array**(order_l - 1)

    ## Column 3
    # Y[0, 2, :] = np.zeros(shape)
    # Y[1, 2, :] = np.zeros(shape)
    Y[2, 2, :] = -density_array * radius_array**order_l  # TODO: Should be negative? SVC changes sign conventions (I flipped back)
    # Y[3, 2, :] = np.zeros(shape)
    Y[4, 2, :] = -radius_array**order_l  # TODO: Should be negative? SVC changes sign conventions (I flipped back)
    Y[5, 2, :] = -(2 * order_l + 1) * radius_array**(order_l - 1)  # TODO: Should be negative? SVC changes sign conventions (I flipped back)

    ## Column 4
    Y[0, 3, :] = (order_l + 1) * radius_array**(-order_l) / (2 * (2 * order_l - 1))
    Y[1, 3, :] = (2 - order_l) * radius_array**(-order_l) / (2 * order_l * (2 * order_l - 1))
    Y[2, 3, :] = ((order_l + 1) * rgp - 2 * (order_l**2 + 3*order_l - 1) * shear_array)\
                 / (2 * (2 * order_l - 1) * radius_array**(order_l + 1))
    Y[3, 3, :] = (order_l**2 - 1) * shear_array / (order_l * (2 * order_l - 1) * radius_array**(order_l + 1))
    # Y[4, 3, :] = np.zeros(shape)
    Y[5, 3, :] = 2 * pi * G * density_array * (order_l + 1) / ((2 * order_l - 1) * radius_array**order_l)

    ## Column 5
    Y[0, 4, :] = radius_array**(-order_l - 2)
    Y[1, 4, :] = -radius_array**(-order_l - 2) / (order_l + 1)
    Y[2, 4, :] = (rgp - 2 * (order_l + 2) * shear_array) / radius_array**(order_l + 3)
    Y[3, 4, :] = 2*(order_l + 2)*shear_array/((order_l + 1)*radius_array**(order_l + 3))
    # Y[4, 4, :] = np.zeros(radius_array.shape)
    Y[5, 4, :] = 4 * pi * G * density_array / radius_array**(order_l + 2)

    ## Column 6
    # Y[0, 5, :] = np.zeros(shape)
    # Y[1, 5, :] = np.zeros(shape)
    Y[2, 5, :] = -density_array / radius_array**(order_l + 1) # TODO: Should be negative? SVC changes sign conventions (I flipped back)
    # Y[3, 5, :] = np.zeros(shape)
    Y[4, 5, :] = -1 / radius_array**(order_l + 1) # TODO: Should be negative? SVC changes sign conventions (I flipped back)
    # Y[5, 5, :] = np.zeros(shape)

    # Calculate Inverse Analytically, See Eqs. 2.45--2.47 in SVC
    #     Inverse matrix is found as Y^-1(r) = D(r)Y^bar(r)

    # D = a diagonal matrix with is non-zero components equal to...
    D_diag = np.zeros((6, 6, num_shells), dtype=np.complex)
    D_diag[0, 0, :] = (order_l + 1) / radius_array**(order_l + 1)
    D_diag[1, 1, :] = order_l * (order_l + 1) / (2 * (2 * order_l - 1) * radius_array**(order_l - 1))
    D_diag[2, 2, :] = +1 / radius_array**(order_l - 1) # TODO: Should be negative? SVC changes sign conventions (I flipped back)
    D_diag[3, 3, :] = order_l * radius_array**order_l
    D_diag[4, 4, :] = radius_array**(order_l + 2) * order_l * (order_l + 1) / (2 * (2 * order_l + 3))
    D_diag[5, 5, :] = -radius_array**(order_l + 1) # TODO: Should be negative? SVC changes sign conventions (I flipped back)

    D = (1 / (2 * order_l + 1)) * D_diag

    # Y^bar is a matrix that has a similar format to the fundamental matrix
    Yinv_bar = np.zeros(Y.shape, dtype=np.complex)

    # Optimizations
    rgp_s = rgp / shear_array
    r_s = radius_array / shear_array
    pr_s = density_array * r_s

    ## Column 1
    Yinv_bar[0, 0, :] = rgp_s - 2 * (order_l + 2)
    Yinv_bar[1, 0, :] = -rgp_s + 2 * (order_l**2 + 3*order_l - 1) / (order_l + 1)
    Yinv_bar[2, 0, :] = 4 * pi * G * density_array
    Yinv_bar[3, 0, :] = rgp_s + 2 * (order_l-1)
    Yinv_bar[4, 0, :] = -rgp_s - 2 * (order_l**2 - order_l - 3) / order_l
    Yinv_bar[5, 0, :] = 4 * pi * G * density_array * radius_array

    ## Column 2
    Yinv_bar[0, 1, :] = 2 * order_l * (order_l + 2)
    Yinv_bar[1, 1, :] = -2 * (order_l**2 - 1)
    # Yinv_bar[2, 1, :] = np.zeros(shape)
    Yinv_bar[3, 1, :] = 2 * (order_l**2 - 1)
    Yinv_bar[4, 1, :] = -2 * order_l * (order_l + 2)
    # Yinv_bar[5, 1, :] = np.zeros(shape)

    ## Column 3
    Yinv_bar[0, 2, :] = -r_s
    Yinv_bar[1, 2, :] = r_s
    # Yinv_bar[2, 2, :] = np.zeros(shape)
    Yinv_bar[3, 2, :] = -r_s
    Yinv_bar[4, 2, :] = r_s
    # Yinv_bar[5, 2, :] = np.zeros(shape)

    ## Column 4
    Yinv_bar[0, 3, :] = order_l * r_s
    Yinv_bar[1, 3, :] = (2 - order_l) * r_s
    # Yinv_bar[2, 3, :] = np.zeros(shape)
    Yinv_bar[3, 3, :] = -(order_l + 1) * r_s
    Yinv_bar[4, 3, :] = (order_l + 3) * r_s
    # Yinv_bar[5, 3, :] = np.zeros(shape)

    ## Column 5
    Yinv_bar[0, 4, :] = pr_s
    Yinv_bar[1, 4, :] = -pr_s
    # Yinv_bar[2, 4, :] = np.zeros(shape)
    Yinv_bar[3, 4, :] = pr_s
    Yinv_bar[4, 4, :] = -pr_s
    Yinv_bar[5, 4, :] = 2 * order_l + 1

    ## Column 6
    # Yinv_bar[0, 5, :] = np.zeros(shape)
    # Yinv_bar[1, 5, :] = np.zeros(shape)
    Yinv_bar[2, 5, :] = -np.ones(num_shells)
    # Yinv_bar[3, 5, :] = np.zeros(shape)
    # Yinv_bar[4, 5, :] = np.zeros(shape)
    Yinv_bar[5, 5, :] = -radius_array

    # Step through each shell index
    Yinv = np.zeros((6, 6, num_shells), dtype=np.complex)
    Y_reduced_shifted = np.zeros((6, 6, num_shells), dtype=np.complex)
    for shell_i in range(num_shells):
        # Calculate the fundamental matrix inverse
        Yinv[:, :, shell_i] = np.matmul(D[:, :, shell_i], Yinv_bar[:, :, shell_i])

        # If shell is not core (innermost) shell then we can calculate the shifted matrix.
        if shell_i > 0:
            Y_reduced_shifted[:, :, shell_i] = np.matmul(Y[:, :, shell_i], Yinv[:, :, shell_i - 1])

    return Y, Yinv, Y_reduced_shifted

#@njit
def find_propagator_matrix(Y_reduced_shifted: np.ndarray, seed_core: np.ndarray) -> np.ndarray:
    """ Build the propagator matrix out of the fundamental matrix and the seed matrix at the core

    See appendix of Henning & Hurford 2014

    Parameters
    ----------
    Y_reduced_shifted : np.ndarray
        Combination of the propagator matrix and its shifted inverse = Y_i @ Y_(i-1)^(-1)
    seed_core : np.ndarray
        Seed matrix at the core (aggregate matrix is built from bottom up)

    Returns
    -------
    propagator_matrix : np.ndarray
        Propagator matrix at each shell

    """

    # We are trying to solve for B where B_i = Y_i Y_(i-1)^(-1) B_(i-1)

    # The shape of Y should be 6 x 6 x number_of_layers
    num_shells = Y_reduced_shifted.shape[2]  # This should be the same for B_core
    seed_shape = seed_core.shape

    B = np.zeros((seed_shape[0], seed_shape[1], num_shells), dtype=np.complex)

    # TODO: is there a more efficient way to do this with out the nested for-loops?
    #    I think it may be better to invert the matrix and solve it that way - but initial tests lead to lots of Singular Matrices when inverting
    for row_i in range(seed_shape[0]):
        for col_i in range(seed_shape[1]):
            for shell_i in range(num_shells):

                if shell_i == 0:
                    B[row_i, col_i, 0] = seed_core[row_i, col_i]
                else:
                    # Recursive otherwise
                    B[row_i, col_i, shell_i] = Y_reduced_shifted[row_i, col_i, shell_i] * B[row_i, col_i, shell_i-1]

    return B