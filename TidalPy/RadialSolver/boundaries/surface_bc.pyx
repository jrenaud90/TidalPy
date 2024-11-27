# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

import numpy as np
cimport numpy as np

from libc.stdio cimport printf
from libc.stdlib cimport exit

from TidalPy.constants cimport d_NAN_DBL


cdef void cf_get_surface_bc(
    double* boundary_conditions_ptr,
    int* bc_model_ptr,
    size_t num_bcs,
    double radius_to_use,
    double bulk_density_to_use,
    double degree_l_dbl,
    ) noexcept nogil:
    """Find the surface boundary condition. """

    # `num_bcs` should equal the length of `bc_model_ptr`
    if num_bcs > 5:
        printf(
            "Unsupported number of boundaries conditions encountered."
            " Provided: %d when maximum supported is 5.", num_bcs)
        exit(-1)
    elif num_bcs <= 0:
        printf(
            "Unsupported number of boundaries conditions encountered."
            " Provided: %d when minimum supported is 1.", num_bcs)
        exit(-1)
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
            boundary_conditions_ptr[j * 3 + 0] = (-1. / 3.) * (2. * degree_l_dbl + 1.) * bulk_density_to_use
            boundary_conditions_ptr[j * 3 + 1] = 0.
            boundary_conditions_ptr[j * 3 + 2] = (2. * degree_l_dbl + 1.) / radius_to_use
        else:
            printf(
                "Unknown boundary condition model: %d. Supported models are:\n"
                " 0: Free Surface.\n"
                " 1: Tidal Potential.\n"
                " 2: Loading Potential.\n", bc_model_ptr[j])
            exit(-1)


def get_surface_bc(
    int[::1] bc_model_view,
    double radius_to_use,
    double bulk_density_to_use,
    double degree_l_dbl,
    ):

    cdef int* bc_model_ptr = &bc_model_view[0]
    cdef size_t num_bcs = bc_model_view.size

    # Build output array
    cdef np.ndarray[np.float64_t, ndim=1] boundary_conditions_arr = np.empty(15, dtype=np.float64)
    cdef double[::1] boundary_conditions_view = boundary_conditions_arr
    cdef double* boundary_conditions_ptr = &boundary_conditions_view[0]
    cf_get_surface_bc(
        boundary_conditions_ptr,
        bc_model_ptr,
        num_bcs,
        radius_to_use,
        bulk_density_to_use,
        degree_l_dbl,
        )
    
    return boundary_conditions_arr
