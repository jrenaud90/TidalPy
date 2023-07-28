# distutils: language = c++
import cython
import numpy as np
cimport numpy as np
np.import_array()

from TidalPy.utilities.performance.array.interp cimport interp, interp_complex
from TidalPy.radial_solver.numerical.derivatives.radial_derivatives_dynamic_x cimport radial_derivatives_solid_general_x

cdef double G = 6.67430e-11

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef void dynamic_solid_ode(double radius, double[:] y_vector, double[:] dy_vector,
                             double[:] radius_array, double complex[:] shear_modulus_array,
                             double[:] bulk_modulus_array, double[:] density_array, double[:] gravity_array,
                             double frequency, int order_l = 2, double G_to_use = G) nogil:
    """ A njit-safe radial derivative function for static, solid layers.

    Parameters
    ----------
    radius : float
        Requested radius to at which the radial derivatives are calculated [m]
    y_vector : np.ndarray
        The radial functions at (or near) the provided radius. Used to estimate the derivatives.
    radius_array : np.ndarray
        Array of radius_array for the interior of a planet or layer where this radial derivative function is valid.
        The bottom most radius_array of a planet should not be equal to zero [m]
    shear_modulus_array : np.ndarray
        Shear modulus at each `radius_array` [Pa] (can be complex for shear dissipation)
    bulk_modulus_array : np.ndarray
        Bulk modulus at each `radius_array` [Pa] (can be complex for bulk dissipation)
    density_array : np.ndarray
        Density at each `radius_array` [kg m-3]
    gravity_array : np.ndarray
        Acceleration due to gravity calculated at each `radius_array` [m s-2]
    frequency : float
        Forcing frequency [rad s-1]
    order_l : int = 2
        Tidal harmonic order
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.
    incompressible : bool = False
        If `True`, the incompressible assumption will be used.

    Returns
    -------
    solid_dyn_derivatives : np.ndarray
        The radial derivatives estimated using the provided radius for a solid layer under the dynamic assumption.
    """

    # These physical parameters are generally given in a coarser resolution than what is required during integration.
    #    Therefore, we interpolate their values at the integrator's requested `radius`.
    cdef double complex shear_modulus = interp_complex(radius, radius_array, shear_modulus_array)
    cdef double density = interp(radius, radius_array, density_array)
    cdef double gravity = interp(radius, radius_array, gravity_array)
    cdef double bulk_modulus = interp(radius, radius_array, bulk_modulus_array)

    radial_derivatives_solid_general_x(
        radius, y_vector, dy_vector, shear_modulus, bulk_modulus, density, gravity, frequency,
        order_l, G_to_use
        )
