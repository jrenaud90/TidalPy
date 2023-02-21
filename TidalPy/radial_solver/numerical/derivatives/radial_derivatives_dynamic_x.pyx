# distutils: language = c++
""" Functions to calculate the Rayleigh wave propagation differential equations

These functions are the most general and allow for:
    - Compressibility
    - Dynamic (non-static) Tides
    - Liquid layer propagation

References
----------
KMN15 : Kamata+ (2015; JGR-P; DOI: 10.1002/2015JE004821)
B15   : Beuthe+ (2015; Icarus; DOI: 10.1016/j.icarus.2015.06.008)
TMS05 : Tobie+ (2005; Icarus; DOI: 10.1016/j.icarus.2005.04.006)
S74   : Saito (1974; J. Phy. Earth; DOI: 10.4294/jpe1952.22.123)
TS72  : Takeuchi, H., and M. Saito (1972), Seismic surface waves, Methods Comput. Phys., 11, 217–295.
"""
import cython
import numpy as np
cimport numpy as np
np.import_array()
from libc.math cimport isnan, M_PI

# Get machine precision.
cdef double EPS
EPS = np.finfo(dtype=np.float64).eps

cdef double G = 6.67430e-11


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef void radial_derivatives_solid_general_x(
    double radius,
    double[:] radial_functions,
    double[:] dy,
    double complex shear_modulus, double bulk_modulus, double density,
    double gravity, double frequency,
    int order_l=2, double G_to_use=G
    ) nogil:
    """ Calculates the radial derivative of the radial functions in the most general form - for solid layers.

    Allows for compressibility and dynamic tides.
    Tidal harmonic l is allowed to be an integer >= 2.

    OPT: This and the other radial derivative functions are called many times during integration.
        As of TidalPy V0.3.4 this function call for floats is around 2.3 micro seconds.

    References
    ----------
    KMN15; B15; TS72

    Parameters
    ----------
    radius : FloatArray
        Radius where the radial functions are calculated. [m]
    radial_functions : np.ndarray
        Tuple of radial functions for a solid layer broken up into real and imaginary portions.
        (y1_real, y1_imag, y2_real, y2_imag, y3_real, y3_imag, y4_real, y4_imag, y5_real, y5_imag, y6_real, y6_imag)
    shear_modulus : NumArray
        Shear modulus (can be complex for dissipation) at `radius` [Pa]
    bulk_modulus : NumArray
        Bulk modulus (can be complex for dissipation) at `radius` [Pa]
    density : FloatArray
        Density at `radius` [kg m-3]
    gravity : FloatArray
        Acceleration due to gravity at `radius` [m s-2]
    frequency : FloatArray
        Forcing frequency (for spin-synchronous tides this is the orbital motion) [rad s-1]
    order_l : int = 2
        Tidal harmonic order.
    G_to_use : float = G
        Newton's gravitational constant. This can be provided as in its MKS units or dimensionless to match the other
        inputs.

    Returns
    -------
    radial_derivatives : np.ndarray
        The derivatives of the radial functions for a solid layer (dynamic assumption)

    """

    cdef double y1_real = radial_functions[0]
    cdef double y1_imag = radial_functions[1]
    cdef double y2_real = radial_functions[2]
    cdef double y2_imag = radial_functions[3]
    cdef double y3_real = radial_functions[4]
    cdef double y3_imag = radial_functions[5]
    cdef double y4_real = radial_functions[6]
    cdef double y4_imag = radial_functions[7]
    cdef double y5_real = radial_functions[8]
    cdef double y5_imag = radial_functions[9]
    cdef double y6_real = radial_functions[10]
    cdef double y6_imag = radial_functions[11]

    # Convert floats to complex
    cdef double complex y1 = y1_real + 1.0j * y1_imag
    cdef double complex y2 = y2_real + 1.0j * y2_imag
    cdef double complex y3 = y3_real + 1.0j * y3_imag
    cdef double complex y4 = y4_real + 1.0j * y4_imag
    cdef double complex y5 = y5_real + 1.0j * y5_imag
    cdef double complex y6 = y6_real + 1.0j * y6_imag

    # Convert compressibility parameters (the first lame parameter can be complex)
    cdef double complex lame = (<double complex>bulk_modulus - (2. / 3.) * shear_modulus)

    # Optimizations
    cdef double degree_l_flt = <double>order_l
    cdef double lp1              = degree_l_flt + 1.
    cdef double lm1              = degree_l_flt - 1.
    cdef double llp1             = degree_l_flt * lp1
    cdef double complex lame_2mu         = lame + 2. * shear_modulus
    cdef double complex lame_2mu_inverse = 1. / lame_2mu
    cdef double r_inverse        = 1. / radius
    cdef double complex two_shear_r_inv  = 2. * shear_modulus * r_inverse
    cdef double density_gravity  = density * gravity
    cdef double dynamic_term     = -frequency * frequency * density * radius
    cdef double grav_term        = 4. * M_PI * G_to_use * density
    cdef double complex y1_y3_term       = 2. * y1 - llp1 * y3

    # See Eq. 82 in TS72 or Eqs. 4--9 in KMN15 or Eqs. 13--18 in B15
    #   Note: There appears to be a missing factor of mu^2 in some of the terms in KMN15.
    # dy2 and dy4 contain all three of: dynamic, viscoelastic, and gravitational terms.
    cdef double complex dy1
    cdef double complex dy2
    cdef double complex dy3
    cdef double complex dy4
    cdef double complex dy5
    cdef double complex dy6

    dy1 = lame_2mu_inverse * (
            y1_y3_term * -lame * r_inverse +
            y2
    )

    dy2 = r_inverse * (
            y1 * (dynamic_term - 2. * density_gravity) +
            y2 * -2. +
            y4 * llp1 +
            y5 * density * lp1 +
            y6 * -density * radius +
            dy1 * 2. * lame +
            y1_y3_term * (2. * (lame + shear_modulus) * r_inverse - density_gravity)
    )

    dy3 = \
        y1 * -r_inverse + \
        y3 * r_inverse + \
        y4 * (1. / shear_modulus)

    dy4 = r_inverse * (
            y1 * (density_gravity + two_shear_r_inv) +
            y3 * (dynamic_term - two_shear_r_inv) +
            y4 * -3. +
            y5 * -density +
            dy1 * -lame +
            y1_y3_term * -lame_2mu * r_inverse
    )

    dy5 = \
        y1 * grav_term + \
        y5 * -lp1 * r_inverse + \
        y6

    dy6 = r_inverse * (
            y1 * grav_term * lm1 +
            y6 * lm1 +
            y1_y3_term * grav_term
    )

    # Convert back to floats
    dy[0] = dy1.real
    dy[1] = dy1.imag
    dy[2] = dy2.real
    dy[3] = dy2.imag
    dy[4] = dy3.real
    dy[5] = dy3.imag
    dy[6] = dy4.real
    dy[7] = dy4.imag
    dy[8] = dy5.real
    dy[9] = dy5.imag
    dy[10] = dy6.real
    dy[11] = dy6.imag


# @njit(cacheable=True)
# def radial_derivatives_liquid_general(
#     radius: FloatArray,
#     radial_functions: np.ndarray,
#     bulk_modulus: NumArray, density: FloatArray,
#     gravity: FloatArray, frequency: FloatArray,
#     order_l: int = 2, G_to_use: float = G
#     ) -> np.ndarray:
#     """ Calculates the radial derivative of the radial functions in the most general form - for liquid layers.
#
#     Allows for compressibility and dynamic tides.
#     Tidal harmonic l is allowed to be an integer >= 2.
#
#     References
#     ----------
#     KMN15; B15; TS72
#
#     Parameters
#     ----------
#     radius : FloatArray
#         Radius where the radial functions are calculated. [m]
#     radial_functions : np.ndarray
#         Tuple of radial functions for a solid layer broken up into real and imaginary portions.
#         (y1_real, y1_imag, y2_real, y2_imag, y5_real, y5_imag, y6_real, y6_imag)
#     bulk_modulus : NumArray
#         Bulk modulus (can be complex for dissipation) at `radius` [Pa]
#     density : FloatArray
#         Density at `radius` [kg m-3]
#     gravity : FloatArray
#         Acceleration due to gravity at `radius` [m s-2]
#     frequency : FloatArray
#         Forcing frequency (for spin-synchronous tides this is the orbital motion) [rad s-1]
#     order_l : int = 2
#         Tidal harmonic order.
#     G_to_use : float = G
#         Newton's gravitational constant. This can be provided as in its MKS units or dimensionless to match the other
#         inputs.
#
#     Returns
#     -------
#     radial_derivatives : np.ndarray
#         The derivatives of the radial functions for a liquid layer (dynamic assumption)
#
#     """
#
#     # For the dynamic version, y4 = 0 always in a liquid layer and y3 is defined by y1, y2, and y5 analytically
#     y1_real = radial_functions[0]
#     y1_imag = radial_functions[1]
#     y2_real = radial_functions[2]
#     y2_imag = radial_functions[3]
#     y5_real = radial_functions[4]
#     y5_imag = radial_functions[5]
#     y6_real = radial_functions[6]
#     y6_imag = radial_functions[7]
#
#     # Convert floats to complex
#     y1 = y1_real + 1.0j * y1_imag
#     y2 = y2_real + 1.0j * y2_imag
#     y5 = y5_real + 1.0j * y5_imag
#     y6 = y6_real + 1.0j * y6_imag
#
#     # Convert compressibility parameters (the first lame parameter can be complex)
#     # For the liquid layer it is assumed that the shear modulus is zero so the lame parameter simply
#     #    equals the bulk modulus
#     lame = bulk_modulus
#
#     # Optimizations
#     lp1               = order_l + 1.
#     lm1               = order_l - 1.
#     llp1              = order_l * lp1
#     lame_inverse      = 1. / lame
#     r_inverse         = 1. / radius
#     density_gravity   = density * gravity
#     f2                = frequency * frequency
#     dynamic_term_no_r = -f2 * density
#     dynamic_term      = dynamic_term_no_r * radius
#     grav_term         = 4. * pi * G_to_use * density
#
#     # y3 derivative is undetermined for a liquid layer, but we can calculate its value which is still used in the
#     #   other derivatives.
#     y3         = (1. / dynamic_term) * (y2 - density_gravity * y1 + density * y5)
#     y1_y3_term = 2. * y1 - llp1 * y3
#
#     # Eqs. 11--14 in KMN15 equations look like they don't match TS72 because they applied the rheology already.
#     #    and substituted y3.
#     # We will use TS72 eq. 87 to allow for a generic rheology and bulk dissipation.
#     # # dy2 contain all three of: dynamic, viscoelastic, and gravitational terms.
#     dy1 = \
#         y2 * lame_inverse - \
#         y1_y3_term * r_inverse
#
#     # TODO: In the solid version there is a [2. * (lame + shear_modulus) * r_inverse] coefficient for y1_y3_term
#     #   In TS72 the first term is gone. Shouldn't Lame + mu = Lame = Bulk for liquid layer?
#     dy2 = \
#         y1 * (dynamic_term_no_r - 2. * density_gravity * r_inverse) + \
#         y5 * density * lp1 * r_inverse - \
#         y6 * density - \
#         y1_y3_term * density_gravity * r_inverse
#
#     dy5 = \
#         y1 * grav_term - \
#         y5 * lp1 * r_inverse + \
#         y6
#
#     dy6 = r_inverse * (
#             lm1 * (y1 * grav_term + y6) +
#             y1_y3_term * grav_term
#     )
#
#     # dy1 = y1 * (-2. * r_inverse + llp1 * gravity / (frequency**2 * radius**2)) + \
#     #     y2 * (lame_inverse - llp1 / (frequency**2 * density * radius**2)) - \
#     #     y5 * llp1 / (frequency**2 * radius**2)
#     #
#     # dy2 = y1 * (-frequency**2 * density - 4. * density * gravity / radius + llp1 * density * gravity**2 / (frequency**2 * radius**2)) - \
#     #     y2 * llp1 * gravity / (frequency**2 * radius**2) + \
#     #     y5 * (lp1 * density / radius) * (1. - order_l * gravity / (frequency**2 * radius)) - \
#     #     y6 * density
#     #
#     # dy5 = y1 * (4. * pi * G_to_use * density) - \
#     #     y5 * (lp1 / radius) + \
#     #     y6
#     #
#     # dy6 = y1 * (4. * pi * lp1 * G_to_use * density / radius) * (1. - order_l * gravity / (frequency**2 * radius)) + \
#     #     y2 * (4. * pi * llp1 * G_to_use / (frequency**2 * radius**2)) + \
#     #     y5 * (4. * pi * llp1 * density * G_to_use / (frequency**2 * radius**2)) + \
#     #     y6 * (lm1 / radius)
#
#     # Build output
#     dy = np.empty(8, dtype=np.float64)
#
#     # Convert back to floats
#     dy[0] = np.real(dy1)
#     dy[1] = np.imag(dy1)
#     dy[2] = np.real(dy2)
#     dy[3] = np.imag(dy2)
#     dy[4] = np.real(dy5)
#     dy[5] = np.imag(dy5)
#     dy[6] = np.real(dy6)
#     dy[7] = np.imag(dy6)
#
#     return dy
