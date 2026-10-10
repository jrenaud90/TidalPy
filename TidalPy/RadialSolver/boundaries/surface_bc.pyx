# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

import numpy as np
cimport numpy as cnp
cnp.import_array()


def get_surface_bc(
        int[::1] bc_model_view,
        double radius_to_use,
        double bulk_density_to_use,
        double degree_l_dbl):
    """
    Get surface boundary conditions for the radial solver.

    Parameters
    ----------
    bc_model_view : int[::1]
        Array of BC model types: 0=free surface, 1=tidal potential, 2=loading potential.
    radius_to_use : float
        Planet radius [m].
    bulk_density_to_use : float
        Planet bulk density [kg m-3].
    degree_l_dbl : float
        Tidal harmonic degree.

    Returns
    -------
    boundary_conditions : ndarray
        Array of shape (15,) with boundary condition values.
    """
    cdef int* bc_model_ptr = &bc_model_view[0]
    cdef size_t num_bcs = bc_model_view.size

    # Build output array
    cdef cnp.ndarray[cnp.float64_t, ndim=1] boundary_conditions_arr = np.empty(15, dtype=np.float64)
    cdef double[::1] boundary_conditions_view = boundary_conditions_arr
    cdef double* boundary_conditions_ptr = &boundary_conditions_view[0]

    cdef int error_code
    error_code = c_get_surface_bc(
        boundary_conditions_ptr,
        bc_model_ptr,
        num_bcs,
        radius_to_use,
        bulk_density_to_use,
        degree_l_dbl,
        )

    if error_code < 0:
        if error_code == -1:
            raise ValueError(f"Unsupported number of boundary conditions. Provided: {num_bcs}, max supported: 5.")
        elif error_code == -2:
            raise ValueError(f"Unsupported number of boundary conditions. Provided: {num_bcs}, min supported: 1.")
        elif error_code == -3:
            raise ValueError("Unknown boundary condition model. Supported: 0=Free, 1=Tidal, 2=Loading.")
        else:
            raise RuntimeError("Unknown error encountered.")

    return boundary_conditions_arr
