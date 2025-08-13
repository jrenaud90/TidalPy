# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

import numpy as np
cimport numpy as cnp

from TidalPy.exceptions import ArgumentException, UnknownModelError
from TidalPy.constants cimport d_NAN_DBL


cdef int cf_get_surface_bc(
        double* boundary_conditions_ptr,
        int* bc_model_ptr,
        size_t num_bcs,
        double radius_to_use,
        double bulk_density_to_use,
        double degree_l_dbl,
        ) noexcept nogil:
    """
    cf_get_surface_bc

    Populates a boundary condition array at a planet's surface based on the specified boundary condition models.

    Parameters
    ----------
    boundary_conditions_ptr : double*
        Pointer to an array where the boundary condition values will be written. The array must have space for 15 elements 
        (5 max allowed models * 3 surface bc's per model). Unused elements are set to `NaN`.
    bc_model_ptr : int*
        Pointer to an array indicating the types of boundary conditions to apply for each model. Each value must correspond 
        to a valid boundary condition type:
            - 0: Free surface
            - 1: Tidal potential
            - 2: Loading potential
    num_bcs : size_t
        Number of boundary condition models specified in `bc_model_ptr`. Must be between 1 and 5 (inclusive).
    radius_to_use : double
        Radius of the planet used for the boundary condition calculations.
    bulk_density_to_use : double
        Bulk density of the planet used for the boundary condition calculations.
    degree_l_dbl : double
        Degree of the spherical harmonic for the boundary condition calculations.

    Returns
    -------
    int
        Return code indicating the status of the function:
            - 0: Success
            - -1: `num_bcs` exceeds 5
            - -2: `num_bcs` is less than or equal to 0
            - -3: Invalid or not implemented boundary condition type found in `bc_model_ptr`

    Notes
    -----
    - The function supports three types of boundary conditions:
        - 0: **Free Surface**: All boundary condition values are set to 0.
        - 1: **Tidal Potential**: Derived from the degree of spherical harmonic and the planet's radius.
        - 2: **Loading Potential**: Uses bulk density, degree of spherical harmonic, and radius.
    - If `num_bcs` is greater than 5 or less than or equal to 0, the function returns an error code.
    - The maximum number of boundary condition models (`num_bcs`) is 5. Unused slots in the output array are initialized to `NaN`.
    - See Eq. 6 in Beuthe (2015) and Eq. 9 in Saito (1974) for loading potential calculations.

    References
    ----------
    - Beuthe (2015)
    - Saito (1974)

    Raises
    ------
    None
        The function does not explicitly raise exceptions. Errors are indicated via return codes.

    Examples
    --------
    # Example usage in a Cython environment:
    >>> cdef double boundary_conditions[15]
    >>> cdef int bc_models[3] = [0, 1, 2]
    >>> cdef int result
    >>> result = cf_get_surface_bc(
            &boundary_conditions[0], &bc_models[0], 
            num_bcs=3, radius_to_use=6.37e6, 
            bulk_density_to_use=5515, degree_l_dbl=2.0
        )
    >>> if result == 0:
    ...     print("Boundary conditions calculated successfully")
    """

    # `num_bcs` should equal the length of `bc_model_ptr`
    if num_bcs > 5:
        return -1
    elif num_bcs <= 0:
        return -2
    cdef size_t i, j

    # Inititalize all boundary conditions to NaN
    # 15 = 5 (max_num_solutions) * 3 (number of surface conditions)
    for i in range(15):
        boundary_conditions_ptr[i] = d_NAN_DBL
    
    for j in range(num_bcs):
        if bc_model_ptr[j] == 0:
            # Free Surface
            boundary_conditions_ptr[j * 3 + 0] = 0.
            boundary_conditions_ptr[j * 3 + 1] = 0.
            boundary_conditions_ptr[j * 3 + 2] = 0.
        elif bc_model_ptr[j] == 1:
            # Tidal Potential
            boundary_conditions_ptr[j * 3 + 0] = 0.
            boundary_conditions_ptr[j * 3 + 1] = 0.
            boundary_conditions_ptr[j * 3 + 2] = (2. * degree_l_dbl + 1.) / radius_to_use
        elif bc_model_ptr[j] == 2:
            # Loading Potential
            # See Eq. 6 in Beuthe (2015) and Eq. 9 of Saito (1974)
            boundary_conditions_ptr[j * 3 + 0] = (-1. / 3.) * (2. * degree_l_dbl + 1.) * bulk_density_to_use
            boundary_conditions_ptr[j * 3 + 1] = 0.
            boundary_conditions_ptr[j * 3 + 2] = (2. * degree_l_dbl + 1.) / radius_to_use
        else:
            return -3
    return 0


def get_surface_bc(
        int[::1] bc_model_view,
        double radius_to_use,
        double bulk_density_to_use,
        double degree_l_dbl,
        ):

    cdef int* bc_model_ptr = &bc_model_view[0]
    cdef size_t num_bcs = bc_model_view.size

    # Build output array
    cdef cnp.ndarray[cnp.float64_t, ndim=1] boundary_conditions_arr = np.empty(15, dtype=np.float64)
    cdef double[::1] boundary_conditions_view = boundary_conditions_arr
    cdef double* boundary_conditions_ptr = &boundary_conditions_view[0]

    # Call cythonized function
    cdef int error_code
    error_code = cf_get_surface_bc(
        boundary_conditions_ptr,
        bc_model_ptr,
        num_bcs,
        radius_to_use,
        bulk_density_to_use,
        degree_l_dbl,
        )
    
    if error_code < 0:
        if error_code == -1:
            raise ArgumentException(f"Unsupported number of boundaries conditions encountered."
                "Provided: {num_bcs} when maximum supported is 5.")
        elif error_code == -2:
            raise ArgumentException(f"Unsupported number of boundaries conditions encountered."
                "Provided: {num_bcs} when minimum supported is 1.")
        elif error_code == -3:
            raise UnknownModelError(f"Unknown boundary condition model. Supported models are:\n"
                "\t0: Free Surface.\n"
                "\t1: Tidal Potential.\n"
                "\t2: Loading Potential.\n")
        else:
            raise RuntimeError("Unknown error encountered.")

    return boundary_conditions_arr
