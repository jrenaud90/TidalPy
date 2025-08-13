# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from TidalPy.utilities.math.complex cimport cf_build_dblcmplx, cf_cinv
from TidalPy.constants cimport d_EPS_DBL
from libc.math cimport fabs
from libcpp cimport bool as cpp_bool
from libcpp.cmath cimport ilogb, ldexp
from libcpp.limits cimport numeric_limits

from libc.float cimport DBL_MIN_EXP

cdef double d_EPS_DBL_10 = 10.0 * d_EPS_DBL

cdef void cf_solid_dynamic_compressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused) noexcept nogil:
    """

    References
    ----------
    KMN15; B15; TS72
    """

    # Recast the additional arguments
    cdef RadialSolverDiffeqArgStruct* rs_args_ptr = <RadialSolverDiffeqArgStruct*>args_ptr

    # Update Equation of State at this radius value

    cdef double[9] eos_array
    cdef double* eos_array_ptr = &eos_array[0]
    
    # Call equation of state solution for this layer.
    rs_args_ptr.eos_solution_ptr.call(rs_args_ptr.layer_index, radius, eos_array_ptr)

    # Pull out results.
    # The EOS stores 9 doubles:
    #   0: Gravity
    #   1: Pressure
    #   2: Mass
    #   3: Moment of Inertia
    #   4: Density
    #   5: Shear Mod (real)
    #   6: Shear Mod (imag)
    #   7: Bulk Mod (real)
    #   8: Bulk Mod (imag)
    cdef double gravity                  = eos_array_ptr[0]
    cdef double density                  = eos_array_ptr[4]
    cdef double complex shear_modulus    = cf_build_dblcmplx(eos_array_ptr[5], eos_array_ptr[6])
    cdef double complex bulk_modulus     = cf_build_dblcmplx(eos_array_ptr[7], eos_array_ptr[8])

    # Pull out y values
    cdef double complex y1 = cf_build_dblcmplx(y_ptr[0], y_ptr[1])
    cdef double complex y2 = cf_build_dblcmplx(y_ptr[2], y_ptr[3])
    cdef double complex y3 = cf_build_dblcmplx(y_ptr[4], y_ptr[5])
    cdef double complex y4 = cf_build_dblcmplx(y_ptr[6], y_ptr[7])
    cdef double complex y5 = cf_build_dblcmplx(y_ptr[8], y_ptr[9])
    cdef double complex y6 = cf_build_dblcmplx(y_ptr[10], y_ptr[11])

    # Convert compressibility parameters (the first lame parameter can be complex)
    cdef double complex lame = (bulk_modulus - (2. / 3.) * shear_modulus)

    # Optimizations
    cdef double r_inverse                = 1. / radius
    cdef double density_gravity          = density * gravity
    cdef double dynamic_term             = -rs_args_ptr.frequency * rs_args_ptr.frequency * density * radius
    cdef double grav_term                = rs_args_ptr.grav_coeff * density
    cdef double complex lame_2mu         = lame + 2. * shear_modulus
    cdef double complex lame_2mu_inverse = 1. / lame_2mu
    cdef double complex two_shear_r_inv  = 2. * shear_modulus * r_inverse
    cdef double complex y1_y3_term       = 2. * y1 - rs_args_ptr.llp1 * y3

    # See Eq. 82 in TS72 or Eqs. 4--9 in KMN15 or Eqs. 13--18 in B15
    #   Note: There appears to be a missing factor of mu^2 in some of the terms in KMN15.
    # dy2 and dy4 contain all three of: dynamic, viscoelastic, and gravitational terms.
    cdef double complex dy1 = \
        lame_2mu_inverse * (
            y1_y3_term * -lame * r_inverse +
            y2
        )

    cdef double complex dy2 = \
        r_inverse * (
            y1 * (dynamic_term - 2. * density_gravity) +
            y2 * -2. +
            y4 * rs_args_ptr.llp1 +
            y5 * density * rs_args_ptr.lp1 +
            y6 * -density * radius +
            dy1 * 2. * lame +
            y1_y3_term * (2. * (lame + shear_modulus) * r_inverse - density_gravity)
        )

    cdef double complex dy3 = \
        y1 * -r_inverse + \
        y3 * r_inverse + \
        y4 * (1. / shear_modulus)

    cdef double complex dy4 = \
        r_inverse * (
            y1 * (density_gravity + two_shear_r_inv) +
            y3 * (dynamic_term - two_shear_r_inv) +
            y4 * -3. +
            y5 * -density +
            dy1 * -lame +
            y1_y3_term * -lame_2mu * r_inverse
        )

    cdef double complex dy5 = \
        y1 * grav_term + \
        y5 * -rs_args_ptr.lp1 * r_inverse + \
        y6

    cdef double complex dy6 = \
        r_inverse * (
            y1 * grav_term * rs_args_ptr.lm1 +
            y6 * rs_args_ptr.lm1 +
            y1_y3_term * grav_term
        )

    # Convert back to floats
    dy_ptr[0]  = dy1.real
    dy_ptr[1]  = dy1.imag
    dy_ptr[2]  = dy2.real
    dy_ptr[3]  = dy2.imag
    dy_ptr[4]  = dy3.real
    dy_ptr[5]  = dy3.imag
    dy_ptr[6]  = dy4.real
    dy_ptr[7]  = dy4.imag
    dy_ptr[8]  = dy5.real
    dy_ptr[9]  = dy5.imag
    dy_ptr[10] = dy6.real
    dy_ptr[11] = dy6.imag


cdef void cf_solid_dynamic_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused) noexcept nogil:
    """

    References
    ----------
    KMN15; B15; TS72
    """

    # Recast the additional arguments
    cdef RadialSolverDiffeqArgStruct* rs_args_ptr = <RadialSolverDiffeqArgStruct*>args_ptr

    # Update Equation of State at this radius value
    cdef double[9] eos_array 
    cdef double* eos_array_ptr = &eos_array[0]
    # Call equation of state solution for this layer.
    rs_args_ptr.eos_solution_ptr.call(rs_args_ptr.layer_index, radius, eos_array_ptr)
    # Pull out results.
    # The EOS stores 9 doubles:
    #   0: Gravity
    #   1: Pressure
    #   2: Mass
    #   3: Moment of Inertia
    #   4: Density
    #   5: Shear Mod (real)
    #   6: Shear Mod (imag)
    #   7: Bulk Mod (real)
    #   8: Bulk Mod (imag)
    cdef double gravity                  = eos_array_ptr[0]
    cdef double density                  = eos_array_ptr[4]
    cdef double complex shear_modulus    = cf_build_dblcmplx(eos_array_ptr[5], eos_array_ptr[6])
    cdef double complex bulk_modulus     = cf_build_dblcmplx(eos_array_ptr[7], eos_array_ptr[8])

    # Pull out y values
    cdef double complex y1 = cf_build_dblcmplx(y_ptr[0], y_ptr[1])
    cdef double complex y2 = cf_build_dblcmplx(y_ptr[2], y_ptr[3])
    cdef double complex y3 = cf_build_dblcmplx(y_ptr[4], y_ptr[5])
    cdef double complex y4 = cf_build_dblcmplx(y_ptr[6], y_ptr[7])
    cdef double complex y5 = cf_build_dblcmplx(y_ptr[8], y_ptr[9])
    cdef double complex y6 = cf_build_dblcmplx(y_ptr[10], y_ptr[11])

    # Optimizations
    cdef double r_inverse               = 1. / radius
    cdef double density_gravity         = density * gravity
    cdef double dynamic_term            = -rs_args_ptr.frequency * rs_args_ptr.frequency * density * radius
    cdef double grav_term               = rs_args_ptr.grav_coeff * density
    cdef double complex two_shear_r_inv = 2. * shear_modulus * r_inverse
    cdef double complex y1_y3_term      = 2. * y1 - rs_args_ptr.llp1 * y3

    # See Eq. 82 in TS72 or Eqs. 4--9 in KMN15 or Eqs. 13--18 in B15
    #   Note: There appears to be a missing factor of mu^2 in some of the terms in KMN15.
    # dy2 and dy4 contain all three of: dynamic, viscoelastic, and gravitational terms.
    cdef double complex dy1 = \
        y1_y3_term * -1. * r_inverse

    cdef double complex dy2 = \
        r_inverse * (
            y1 * (dynamic_term + 12. * shear_modulus * r_inverse - 4. * density_gravity) +
            y3 * rs_args_ptr.llp1 * (density_gravity - 6. * shear_modulus * r_inverse) +
            y4 * rs_args_ptr.llp1 +
            y5 * density * rs_args_ptr.lp1 +
            y6 * -density * radius
        )

    cdef double complex dy3 = \
        y1 * -r_inverse + \
        y3 * r_inverse + \
        y4 * (1. / shear_modulus)

    cdef double complex dy4 = \
        r_inverse * (
            y1 * (density_gravity - 3. * two_shear_r_inv) +
            y2 * -1. +
            y3 * (dynamic_term + two_shear_r_inv * (2. * rs_args_ptr.llp1 - 1.)) +
            y4 * -3. +
            y5 * -density
        )

    cdef double complex dy5 = \
        y1 * grav_term + \
        y5 * -rs_args_ptr.lp1 * r_inverse + \
        y6

    cdef double complex dy6 = \
        r_inverse * (
            y1 * grav_term * rs_args_ptr.lm1 +
            y6 * rs_args_ptr.lm1 +
            y1_y3_term * grav_term
        )

    # Convert back to floats
    dy_ptr[0]  = dy1.real
    dy_ptr[1]  = dy1.imag
    dy_ptr[2]  = dy2.real
    dy_ptr[3]  = dy2.imag
    dy_ptr[4]  = dy3.real
    dy_ptr[5]  = dy3.imag
    dy_ptr[6]  = dy4.real
    dy_ptr[7]  = dy4.imag
    dy_ptr[8]  = dy5.real
    dy_ptr[9]  = dy5.imag
    dy_ptr[10] = dy6.real
    dy_ptr[11] = dy6.imag


cdef void cf_solid_static_compressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused) noexcept nogil:

    # Recast the additional arguments
    cdef RadialSolverDiffeqArgStruct* rs_args_ptr = <RadialSolverDiffeqArgStruct*>args_ptr

    # Update Equation of State at this radius value
    cdef double[9] eos_array 
    cdef double* eos_array_ptr = &eos_array[0]
    # Call equation of state solution for this layer.
    rs_args_ptr.eos_solution_ptr.call(rs_args_ptr.layer_index, radius, eos_array_ptr)
    # Pull out results.
    # The EOS stores 9 doubles:
    #   0: Gravity
    #   1: Pressure
    #   2: Mass
    #   3: Moment of Inertia
    #   4: Density
    #   5: Shear Mod (real)
    #   6: Shear Mod (imag)
    #   7: Bulk Mod (real)
    #   8: Bulk Mod (imag)
    cdef double gravity                  = eos_array_ptr[0]
    cdef double density                  = eos_array_ptr[4]
    cdef double complex shear_modulus    = cf_build_dblcmplx(eos_array_ptr[5], eos_array_ptr[6])
    cdef double complex bulk_modulus     = cf_build_dblcmplx(eos_array_ptr[7], eos_array_ptr[8])

    # Pull out y values
    cdef double complex y1 = cf_build_dblcmplx(y_ptr[0], y_ptr[1])
    cdef double complex y2 = cf_build_dblcmplx(y_ptr[2], y_ptr[3])
    cdef double complex y3 = cf_build_dblcmplx(y_ptr[4], y_ptr[5])
    cdef double complex y4 = cf_build_dblcmplx(y_ptr[6], y_ptr[7])
    cdef double complex y5 = cf_build_dblcmplx(y_ptr[8], y_ptr[9])
    cdef double complex y6 = cf_build_dblcmplx(y_ptr[10], y_ptr[11])

    # Convert compressibility parameters (the first lame parameter can be complex)
    cdef double complex lame = (bulk_modulus - (2. / 3.) * shear_modulus)

    # Optimizations
    cdef double r_inverse                = 1. / radius
    cdef double density_gravity          = density * gravity
    cdef double grav_term                = rs_args_ptr.grav_coeff * density
    cdef double complex lame_2mu         = lame + 2. * shear_modulus
    cdef double complex lame_2mu_inverse = 1. / lame_2mu
    cdef double complex two_shear_r_inv  = 2. * shear_modulus * r_inverse
    cdef double complex y1_y3_term       = 2. * y1 - rs_args_ptr.llp1 * y3

    # See Eq. 82 in TS72 or Eqs. 4--9 in KMN15 or Eqs. 13--18 in B15
    #   Note: There appears to be a missing factor of mu^2 in some of the terms in KMN15.
    # The static case just sets all frequency_to_use dependence in these equations to zero.
    # dy2 and dy4 contain: viscoelastic, and gravitational terms.
    cdef double complex dy1 = \
        lame_2mu_inverse * (
            y1_y3_term * -lame * r_inverse +
            y2
        )

    cdef double complex dy2 = \
        r_inverse * (
            y1 * -2. * density_gravity +
            y2 * -2. +
            y4 * rs_args_ptr.llp1 +
            y5 * density * rs_args_ptr.lp1 +
            y6 * -density * radius +
            dy1 * 2. * lame +
            y1_y3_term * (2. * (lame + shear_modulus) * r_inverse - density_gravity)
        )

    cdef double complex dy3 = \
        y1 * -r_inverse + \
        y3 * r_inverse + \
        y4 * (1. / shear_modulus)

    cdef double complex dy4 = \
        r_inverse * (
            y1 * (density_gravity + two_shear_r_inv) +
            y3 * -two_shear_r_inv +
            y4 * -3. +
            y5 * -density +
            dy1 * -lame +
            y1_y3_term * -lame_2mu * r_inverse
        )

    cdef double complex dy5 = \
        y1 * grav_term + \
        y5 * -rs_args_ptr.lp1 * r_inverse + \
        y6

    cdef double complex dy6 = \
        r_inverse * (
            y1 * grav_term * rs_args_ptr.lm1 +
            y6 * rs_args_ptr.lm1 +
            y1_y3_term * grav_term
        )

    # Convert back to floats
    dy_ptr[0]  = dy1.real
    dy_ptr[1]  = dy1.imag
    dy_ptr[2]  = dy2.real
    dy_ptr[3]  = dy2.imag
    dy_ptr[4]  = dy3.real
    dy_ptr[5]  = dy3.imag
    dy_ptr[6]  = dy4.real
    dy_ptr[7]  = dy4.imag
    dy_ptr[8]  = dy5.real
    dy_ptr[9]  = dy5.imag
    dy_ptr[10] = dy6.real
    dy_ptr[11] = dy6.imag


cdef void cf_solid_static_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused) noexcept nogil:

    # Recast the additional arguments
    cdef RadialSolverDiffeqArgStruct* rs_args_ptr = <RadialSolverDiffeqArgStruct*>args_ptr

    # Update Equation of State at this radius value
    cdef double[9] eos_array 
    cdef double* eos_array_ptr = &eos_array[0]
    # Call equation of state solution for this layer.
    rs_args_ptr.eos_solution_ptr.call(rs_args_ptr.layer_index, radius, eos_array_ptr)
    # Pull out results.
    # The EOS stores 9 doubles:
    #   0: Gravity
    #   1: Pressure
    #   2: Mass
    #   3: Moment of Inertia
    #   4: Density
    #   5: Shear Mod (real)
    #   6: Shear Mod (imag)
    #   7: Bulk Mod (real)
    #   8: Bulk Mod (imag)
    cdef double gravity                  = eos_array_ptr[0]
    cdef double density                  = eos_array_ptr[4]
    cdef double complex shear_modulus    = cf_build_dblcmplx(eos_array_ptr[5], eos_array_ptr[6])
    cdef double complex bulk_modulus     = cf_build_dblcmplx(eos_array_ptr[7], eos_array_ptr[8])

    # Pull out y values
    cdef double complex y1 = cf_build_dblcmplx(y_ptr[0], y_ptr[1])
    cdef double complex y2 = cf_build_dblcmplx(y_ptr[2], y_ptr[3])
    cdef double complex y3 = cf_build_dblcmplx(y_ptr[4], y_ptr[5])
    cdef double complex y4 = cf_build_dblcmplx(y_ptr[6], y_ptr[7])
    cdef double complex y5 = cf_build_dblcmplx(y_ptr[8], y_ptr[9])
    cdef double complex y6 = cf_build_dblcmplx(y_ptr[10], y_ptr[11])

    # Optimizations
    cdef double r_inverse               = 1. / radius
    cdef double density_gravity         = density * gravity
    cdef double grav_term               = rs_args_ptr.grav_coeff * density
    cdef double complex two_shear_r_inv = 2. * shear_modulus * r_inverse
    cdef double complex y1_y3_term      = 2. * y1 - rs_args_ptr.llp1 * y3

    # Solve for radial derivatives
    cdef double complex dy1 = \
        -1. * y1_y3_term * r_inverse

    cdef double complex dy2 = \
        r_inverse * (
            y1 * (12. * shear_modulus * r_inverse - 4. * density_gravity) +
            y3 * rs_args_ptr.llp1 * (density_gravity - 6. * shear_modulus * r_inverse) +
            y4 * rs_args_ptr.llp1 +
            y5 * density * rs_args_ptr.lp1 +
            y6 * -density * radius
        )

    cdef double complex dy3 = \
        y1 * -r_inverse + \
        y3 * r_inverse + \
        y4 * (1. / shear_modulus)

    cdef double complex dy4 = \
        r_inverse * (
            y1 * (density_gravity - 3. * two_shear_r_inv) +
            y2 * -1. +
            y3 * (two_shear_r_inv * (2. * rs_args_ptr.llp1 - 1.)) +
            y4 * -3. +
            y5 * -density
        )

    cdef double complex dy5 = \
        y1 * grav_term + \
        y5 * -rs_args_ptr.lp1 * r_inverse + \
        y6

    cdef double complex dy6 = \
        r_inverse * (
            y1 * grav_term * rs_args_ptr.lm1 +
            y6 * rs_args_ptr.lm1 +
            y1_y3_term * grav_term
        )

    # Convert back to floats
    dy_ptr[0]  = dy1.real
    dy_ptr[1]  = dy1.imag
    dy_ptr[2]  = dy2.real
    dy_ptr[3]  = dy2.imag
    dy_ptr[4]  = dy3.real
    dy_ptr[5]  = dy3.imag
    dy_ptr[6]  = dy4.real
    dy_ptr[7]  = dy4.imag
    dy_ptr[8]  = dy5.real
    dy_ptr[9]  = dy5.imag
    dy_ptr[10] = dy6.real
    dy_ptr[11] = dy6.imag

cdef cpp_bool equal_zero_within_ulps(double x, size_t n) noexcept nogil:
    #// Since `epsilon()` is the gap size (ULP, unit in the last place)
    #// of floating-point numbers in interval [1, 2), we can scale it to
    #// the gap size in interval [2^e, 2^{e+1}), where `e` is the exponent
    #// of `x` and `y`.
 
    #// If `x` and `y` have different gap sizes (which means they have
    #// different exponents), we take the smaller one. Taking the bigger
    #// one is also reasonable, I guess.
    cdef double m = min(fabs(x), 0.0)
 
    # // Subnormal numbers have fixed exponent, which is `min_exponent - 1`.
    cdef int exp
    if m < numeric_limits[double].min():
        exp = <int>DBL_MIN_EXP - 1
    else:
        exp = ilogb(m)
 
    #// We consider `x` and `y` equal if the difference between them is
    #// within `n` ULPs.
    return fabs(x) <= n * ldexp(numeric_limits[double].epsilon(), exp)

cdef void cf_liquid_dynamic_compressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused) noexcept nogil:

    # Recast the additional arguments
    cdef RadialSolverDiffeqArgStruct* rs_args_ptr = <RadialSolverDiffeqArgStruct*>args_ptr

    # Update Equation of State at this radius value
    cdef double[9] eos_array 
    cdef double* eos_array_ptr = &eos_array[0]
    # Call equation of state solution for this layer.
    rs_args_ptr.eos_solution_ptr.call(rs_args_ptr.layer_index, radius, eos_array_ptr)
    # Pull out results.
    # The EOS stores 9 doubles:
    #   0: Gravity
    #   1: Pressure
    #   2: Mass
    #   3: Moment of Inertia
    #   4: Density
    #   5: Shear Mod (real)
    #   6: Shear Mod (imag)
    #   7: Bulk Mod (real)
    #   8: Bulk Mod (imag)
    cdef double gravity                  = eos_array_ptr[0]
    cdef double density                  = eos_array_ptr[4]
    cdef double complex shear_modulus    = cf_build_dblcmplx(eos_array_ptr[5], eos_array_ptr[6])
    cdef double complex bulk_modulus     = cf_build_dblcmplx(eos_array_ptr[7], eos_array_ptr[8])

    # Pull out y values
    # For the dynamic version, y4 = 0 always in a liquid layer and y3 is defined by y1, y2, and y5 analytically
    cdef double complex y1 = cf_build_dblcmplx(y_ptr[0], y_ptr[1])
    cdef double complex y2 = cf_build_dblcmplx(y_ptr[2], y_ptr[3])
    cdef double complex y5 = cf_build_dblcmplx(y_ptr[4], y_ptr[5])
    cdef double complex y6 = cf_build_dblcmplx(y_ptr[6], y_ptr[7])
    # Optimizations
    cdef double r_inverse         = 1. / radius
    cdef double density_gravity   = density * gravity
    cdef double f2                = rs_args_ptr.frequency * rs_args_ptr.frequency
    cdef double dynamic_term_no_r = -f2 * density
    cdef double dynamic_term      = dynamic_term_no_r * radius
    cdef double grav_term         = rs_args_ptr.grav_coeff * density

    # Check if dynamic term is close to zero. It will always be negative so compare to negative eps
    # if dynamic_term < EPS_100:
    #     # TODO: is faking this okay?
    #     dynamic_term = EPS_100

    # For the liquid layer it is assumed that the shear modulus is zero so the lame parameter simply
    #    equals the bulk modulus. Until bulk dissipation is considered, it will always be real-valued
    cdef double complex lame_inverse = cf_cinv(bulk_modulus)

    # y3 derivative is undetermined for a liquid layer, but we can calculate its value which is still used in the
    #   other derivatives.
    
    # TODO: There are issues with this summation and the floating point errors. If you sum term 1 + 2 + 3 you get one
    # Answer, if you sum 1 + 3 + 2 you get another. In either case the sum is close to zero in a lot of cases
    # So the differences here can cause sign changes or large swings in outcomes. 
    # cdef double complex y3 = (1. / dynamic_term) * (y2 - density_gravity * y1 + density * y5)

    # After some testing, found that dividing out the y2 and then checking if the sum is < eps may help.
    # cdef double complex y3 = (y2 / dynamic_term) * (1.0 - density_gravity * y1/y2 + density * y5/y2)
    # cdef double complex y1y2   = y1 / y2
    # cdef double complex y5y2   = y5 / y2
    cdef double complex coeff_r  = (rs_args_ptr.llp1 / (f2 * radius))
    # cdef double complex y3_sum = (y2 + density * y5)
    # y3_sum -= density_gravity * y1

    # if cf_cabs(y3_sum) < d_EPS_DBL:
    #     y3 = cf_build_dblcmplx(0.0, 0.0)
    #     y1_y3_term = 2. * y1
    # else:
    # y3 = coeff * y3_sum
    cdef double complex y1_y3_term = \
        (y1 * (2.0 - gravity * coeff_r)) + \
        (y2 * coeff_r / density) + \
        (y5 * coeff_r)

    # Eqs. 11--14 in KMN15 equations look like they don't match TS72 because they applied the rheology already.
    #    and substituted y3.
    # We will use TS72 eq. 87 to allow for a generic rheology and bulk dissipation.
    # dy2 contain all three of: dynamic, viscoelastic, and gravitational terms.
    cdef double complex dy1 = \
        y2 * lame_inverse - \
        y1_y3_term * r_inverse

    # TODO: In the solid version there is a [2. * (lame + shear_modulus) * r_inverse] coefficient for y1_y3_term
    #   In TS72 the first term is gone. Shouldn't Lame + mu = Lame = Bulk for liquid layer?
    cdef double complex dy2 = \
        y1 * (dynamic_term_no_r - 2. * density_gravity * r_inverse) + \
        y5 * density * rs_args_ptr.lp1 * r_inverse - \
        y6 * density - \
        y1_y3_term * density_gravity * r_inverse

    cdef double complex dy5 = \
        y1 * grav_term - \
        y5 * rs_args_ptr.lp1 * r_inverse + \
        y6

    cdef double complex dy6 = \
        r_inverse * (
            rs_args_ptr.lm1 * (y1 * grav_term + y6) +
            y1_y3_term * grav_term
        )

    # Convert back to floats
    dy_ptr[0] = dy1.real
    dy_ptr[1] = dy1.imag
    dy_ptr[2] = dy2.real
    dy_ptr[3] = dy2.imag
    dy_ptr[4] = dy5.real
    dy_ptr[5] = dy5.imag
    dy_ptr[6] = dy6.real
    dy_ptr[7] = dy6.imag


cdef void cf_liquid_dynamic_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused) noexcept nogil:

    # Recast the additional arguments
    cdef RadialSolverDiffeqArgStruct* rs_args_ptr = <RadialSolverDiffeqArgStruct*>args_ptr

    # Update Equation of State at this radius value
    cdef double[9] eos_array 
    cdef double* eos_array_ptr = &eos_array[0]
    # Call equation of state solution for this layer.
    rs_args_ptr.eos_solution_ptr.call(rs_args_ptr.layer_index, radius, eos_array_ptr)
    # Pull out results.
    # The EOS stores 9 doubles:
    #   0: Gravity
    #   1: Pressure
    #   2: Mass
    #   3: Moment of Inertia
    #   4: Density
    #   5: Shear Mod (real)
    #   6: Shear Mod (imag)
    #   7: Bulk Mod (real)
    #   8: Bulk Mod (imag)
    cdef double gravity                  = eos_array_ptr[0]
    cdef double density                  = eos_array_ptr[4]
    cdef double complex shear_modulus    = cf_build_dblcmplx(eos_array_ptr[5], eos_array_ptr[6])
    cdef double complex bulk_modulus     = cf_build_dblcmplx(eos_array_ptr[7], eos_array_ptr[8])

    # Pull out y values
    # For the dynamic version, y4 = 0 always in a liquid layer and y3 is defined by y1, y2, and y5 analytically
    cdef double complex y1 = cf_build_dblcmplx(y_ptr[0], y_ptr[1])
    cdef double complex y2 = cf_build_dblcmplx(y_ptr[2], y_ptr[3])
    cdef double complex y5 = cf_build_dblcmplx(y_ptr[4], y_ptr[5])
    cdef double complex y6 = cf_build_dblcmplx(y_ptr[6], y_ptr[7])

    # Optimizations
    cdef double r_inverse       = 1. / radius
    cdef double density_gravity = density * gravity
    cdef double dynamic_term    = -rs_args_ptr.frequency * rs_args_ptr.frequency * density * radius
    cdef double grav_term       = rs_args_ptr.grav_coeff * density

    # Check if dynamic term is close to zero. It will always be negative so compare to negative eps
    # if dynamic_term < EPS_100:
    #     # TODO: is faking this okay?
    #     dynamic_term = EPS_100

    # y3 derivative is undetermined for a liquid layer, but we can calculate its value which is still used in the
    #   other derivatives.
    cdef double complex y3 = (1. / dynamic_term) * (y2 + density * y5 - density_gravity * y1)
    cdef double complex y1_y3_term = 2. * y1 - rs_args_ptr.llp1 * y3

    cdef double complex dy1 = \
        y1_y3_term * -r_inverse

    cdef double complex dy2 = \
        r_inverse * (
            y1 * (dynamic_term - 2. * density_gravity) +
            y5 * density * rs_args_ptr.lp1 +
            y6 * -density * radius +
            # TODO: In the solid version there is a [2. * (lame + shear_modulus) * r_inverse] coefficient for y1_y3_term
            #   In TS72 the first term is gone. Shouldn't Lame + mu = Lame = Bulk for liquid layer?
            y1_y3_term * -density_gravity
        )

    cdef double complex dy5 = \
        y1 * grav_term + \
        y5 * -rs_args_ptr.lp1 * r_inverse + \
        y6

    cdef double complex dy6 = \
        r_inverse * (
            y1 * grav_term * rs_args_ptr.lm1 +
            y6 * rs_args_ptr.lm1 +
            y1_y3_term * grav_term
        )

    # Convert back to floats
    dy_ptr[0] = dy1.real
    dy_ptr[1] = dy1.imag
    dy_ptr[2] = dy2.real
    dy_ptr[3] = dy2.imag
    dy_ptr[4] = dy5.real
    dy_ptr[5] = dy5.imag
    dy_ptr[6] = dy6.real
    dy_ptr[7] = dy6.imag


cdef void cf_liquid_static_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused) noexcept nogil:

    # Recast the additional arguments
    cdef RadialSolverDiffeqArgStruct* rs_args_ptr = <RadialSolverDiffeqArgStruct*>args_ptr

    # Update Equation of State at this radius value
    cdef double[9] eos_array 
    cdef double* eos_array_ptr = &eos_array[0]
    # Call equation of state solution for this layer.
    rs_args_ptr.eos_solution_ptr.call(rs_args_ptr.layer_index, radius, eos_array_ptr)
    # Pull out results.
    # The EOS stores 9 doubles:
    #   0: Gravity
    #   1: Pressure
    #   2: Mass
    #   3: Moment of Inertia
    #   4: Density
    #   5: Shear Mod (real)
    #   6: Shear Mod (imag)
    #   7: Bulk Mod (real)
    #   8: Bulk Mod (imag)
    cdef double gravity                  = eos_array_ptr[0]
    cdef double density                  = eos_array_ptr[4]
    cdef double complex shear_modulus    = cf_build_dblcmplx(eos_array_ptr[5], eos_array_ptr[6])
    cdef double complex bulk_modulus     = cf_build_dblcmplx(eos_array_ptr[7], eos_array_ptr[8])

    # For the static liquid version, only y5 and y7 are defined.
    cdef double complex y5 = cf_build_dblcmplx(y_ptr[0], y_ptr[1])
    cdef double complex y7 = cf_build_dblcmplx(y_ptr[2], y_ptr[3])

    # Optimizations
    cdef double r_inverse = 1. / radius
    cdef double grav_term = rs_args_ptr.grav_coeff * density / gravity

    # See Eq. 18 in S75
    cdef double complex dy5 = \
        y5 * (grav_term - rs_args_ptr.lp1 * r_inverse) + \
        y7

    cdef double complex dy7 = \
        y5 * 2. * rs_args_ptr.lm1 * r_inverse * grav_term + \
        y7 * (rs_args_ptr.lm1 * r_inverse - grav_term)

    # Convert back to floats; stor in output array
    dy_ptr[0] = dy5.real
    dy_ptr[1] = dy5.imag
    dy_ptr[2] = dy7.real
    dy_ptr[3] = dy7.imag


cdef DiffeqFuncType cf_find_layer_diffeq(
        int layer_type,
        int layer_is_static,
        int layer_is_incomp) noexcept nogil:
    
    if layer_type == 0:
        # Solid
        if layer_is_static == 1:
            # Solid-Static
            if layer_is_incomp == 1:
                # Solid-Static-Incomp
                return cf_solid_static_incompressible
            else:
                # Solid-Static-Comp
                return cf_solid_static_compressible
        else:
            # Solid-Dynamic
            if layer_is_incomp:
                # Solid-Dynamic-Incomp
                return cf_solid_dynamic_incompressible
                
            else:
                # Solid-Dynamic-Comp
                return cf_solid_dynamic_compressible
    else:
        # Liquid
        if layer_is_static == 1:
            # Liquid-Static
            if layer_is_incomp == 1:
                # Liquid-Static-Incomp
                return cf_liquid_static_incompressible
            else:
                # Liquid-Static-Comp
                # TODO: Compressible is the same as the incompressible version; Check if this is true.
                return cf_liquid_static_incompressible
        else:
            # Liquid-Dynamic
            if layer_is_incomp:
                # Liquid-Dynamic-Incomp
                return cf_liquid_dynamic_incompressible
                
            else:
                # Liquid-Dynamic-Comp
                return cf_liquid_dynamic_compressible
