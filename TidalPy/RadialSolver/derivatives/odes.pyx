# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from TidalPy.utilities.math.complex cimport cf_build_dblcmplx


cdef void cf_solid_dynamic_compressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        const void* void_args_ptr,
        PreEvalFunc unused) noexcept nogil:
    """

    References
    ----------
    KMN15; B15; TS72
    """

    # Recast the additional arguments
    cdef RadialSolverDiffeqArgStruct* args_ptr = <RadialSolverDiffeqArgStruct*>void_args_ptr

    # Update Equation of State at this radius value

    cdef double[7] eos_array
    cdef double* eos_array_ptr = &eos_array[0]
    # Call equation of state solution for this layer.
    args_ptr.eos_solution_sptr.get().call(radius, eos_array_ptr)

    # Pull out results.
    # The EOS stores 7 doubles:
    #   0: Gravity
    #   1: Pressure
    #   2: Density
    #   3: Shear Mod (real)
    #   4: Shear Mod (imag)
    #   5: Bulk Mod (real)
    #   6: Bulk Mod (imag)
    cdef double gravity                  = eos_array_ptr[0]
    cdef double density                  = eos_array_ptr[2]
    cdef double complex shear_modulus    = cf_build_dblcmplx(eos_array_ptr[3], eos_array_ptr[4])
    cdef double complex bulk_modulus     = cf_build_dblcmplx(eos_array_ptr[5], eos_array_ptr[6])

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
    cdef double dynamic_term             = -args_ptr.frequency * args_ptr.frequency * density * radius
    cdef double grav_term                = args_ptr.grav_coeff * density
    cdef double complex lame_2mu         = lame + 2. * shear_modulus
    cdef double complex lame_2mu_inverse = 1. / lame_2mu
    cdef double complex two_shear_r_inv  = 2. * shear_modulus * r_inverse
    cdef double complex y1_y3_term       = 2. * y1 - args_ptr.llp1 * y3

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
            y4 * args_ptr.llp1 +
            y5 * density * args_ptr.lp1 +
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
        y5 * -args_ptr.lp1 * r_inverse + \
        y6

    cdef double complex dy6 = \
        r_inverse * (
            y1 * grav_term * args_ptr.lm1 +
            y6 * args_ptr.lm1 +
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
        const void* void_args_ptr,
        PreEvalFunc unused) noexcept nogil:
    """

    References
    ----------
    KMN15; B15; TS72
    """

    # Recast the additional arguments
    cdef RadialSolverDiffeqArgStruct* args_ptr = <RadialSolverDiffeqArgStruct*>void_args_ptr

    # Update Equation of State at this radius value
    cdef double[7] eos_array 
    cdef double* eos_array_ptr = &eos_array[0]
    # Call equation of state solution for this layer.
    args_ptr.eos_solution_sptr.get().call(radius, eos_array_ptr)
    # Pull out results.
    # The EOS stores 7 doubles:
    #   0: Gravity
    #   1: Pressure
    #   2: Density
    #   3: Shear Mod (real)
    #   4: Shear Mod (imag)
    #   5: Bulk Mod (real)
    #   6: Bulk Mod (imag)
    cdef double gravity                  = eos_array_ptr[0]
    cdef double density                  = eos_array_ptr[2]
    cdef double complex shear_modulus    = cf_build_dblcmplx(eos_array_ptr[3], eos_array_ptr[4])
    cdef double complex bulk_modulus     = cf_build_dblcmplx(eos_array_ptr[5], eos_array_ptr[6])

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
    cdef double dynamic_term            = -args_ptr.frequency * args_ptr.frequency * density * radius
    cdef double grav_term               = args_ptr.grav_coeff * density
    cdef double complex two_shear_r_inv = 2. * shear_modulus * r_inverse
    cdef double complex y1_y3_term      = 2. * y1 - args_ptr.llp1 * y3

    # See Eq. 82 in TS72 or Eqs. 4--9 in KMN15 or Eqs. 13--18 in B15
    #   Note: There appears to be a missing factor of mu^2 in some of the terms in KMN15.
    # dy2 and dy4 contain all three of: dynamic, viscoelastic, and gravitational terms.
    cdef double complex dy1 = \
        y1_y3_term * -1. * r_inverse

    cdef double complex dy2 = \
        r_inverse * (
            y1 * (dynamic_term + 12. * shear_modulus * r_inverse - 4. * density_gravity) +
            y3 * args_ptr.llp1 * (density_gravity - 6. * shear_modulus * r_inverse) +
            y4 * args_ptr.llp1 +
            y5 * density * args_ptr.lp1 +
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
            y3 * (dynamic_term + two_shear_r_inv * (2. * args_ptr.llp1 - 1.)) +
            y4 * -3. +
            y5 * -density
        )

    cdef double complex dy5 = \
        y1 * grav_term + \
        y5 * -args_ptr.lp1 * r_inverse + \
        y6

    cdef double complex dy6 = \
        r_inverse * (
            y1 * grav_term * args_ptr.lm1 +
            y6 * args_ptr.lm1 +
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
        const void* void_args_ptr,
        PreEvalFunc unused) noexcept nogil:

    # Recast the additional arguments
    cdef RadialSolverDiffeqArgStruct* args_ptr = <RadialSolverDiffeqArgStruct*>void_args_ptr

    # Update Equation of State at this radius value
    cdef double[7] eos_array 
    cdef double* eos_array_ptr = &eos_array[0]
    # Call equation of state solution for this layer.
    args_ptr.eos_solution_sptr.get().call(radius, eos_array_ptr)
    # Pull out results.
    # The EOS stores 7 doubles:
    #   0: Gravity
    #   1: Pressure
    #   2: Density
    #   3: Shear Mod (real)
    #   4: Shear Mod (imag)
    #   5: Bulk Mod (real)
    #   6: Bulk Mod (imag)
    cdef double gravity                  = eos_array_ptr[0]
    cdef double density                  = eos_array_ptr[2]
    cdef double complex shear_modulus    = cf_build_dblcmplx(eos_array_ptr[3], eos_array_ptr[4])
    cdef double complex bulk_modulus     = cf_build_dblcmplx(eos_array_ptr[5], eos_array_ptr[6])

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
    cdef double grav_term                = args_ptr.grav_coeff * density
    cdef double complex lame_2mu         = lame + 2. * shear_modulus
    cdef double complex lame_2mu_inverse = 1. / lame_2mu
    cdef double complex two_shear_r_inv  = 2. * shear_modulus * r_inverse
    cdef double complex y1_y3_term       = 2. * y1 - args_ptr.llp1 * y3

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
            y4 * args_ptr.llp1 +
            y5 * density * args_ptr.lp1 +
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
        y5 * -args_ptr.lp1 * r_inverse + \
        y6

    cdef double complex dy6 = \
        r_inverse * (
            y1 * grav_term * args_ptr.lm1 +
            y6 * args_ptr.lm1 +
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
        const void* void_args_ptr,
        PreEvalFunc unused) noexcept nogil:

    # Recast the additional arguments
    cdef RadialSolverDiffeqArgStruct* args_ptr = <RadialSolverDiffeqArgStruct*>void_args_ptr

    # Update Equation of State at this radius value
    cdef double[7] eos_array 
    cdef double* eos_array_ptr = &eos_array[0]
    # Call equation of state solution for this layer.
    args_ptr.eos_solution_sptr.get().call(radius, eos_array_ptr)
    # Pull out results.
    # The EOS stores 7 doubles:
    #   0: Gravity
    #   1: Pressure
    #   2: Density
    #   3: Shear Mod (real)
    #   4: Shear Mod (imag)
    #   5: Bulk Mod (real)
    #   6: Bulk Mod (imag)
    cdef double gravity                  = eos_array_ptr[0]
    cdef double density                  = eos_array_ptr[2]
    cdef double complex shear_modulus    = cf_build_dblcmplx(eos_array_ptr[3], eos_array_ptr[4])
    cdef double complex bulk_modulus     = cf_build_dblcmplx(eos_array_ptr[5], eos_array_ptr[6])

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
    cdef double grav_term               = args_ptr.grav_coeff * density
    cdef double complex two_shear_r_inv = 2. * shear_modulus * r_inverse
    cdef double complex y1_y3_term      = 2. * y1 - args_ptr.llp1 * y3

    # Solve for radial derivatives
    cdef double complex dy1 = \
        -1. * y1_y3_term * r_inverse

    cdef double complex dy2 = \
        r_inverse * (
            y1 * (12. * shear_modulus * r_inverse - 4. * density_gravity) +
            y3 * args_ptr.llp1 * (density_gravity - 6. * shear_modulus * r_inverse) +
            y4 * args_ptr.llp1 +
            y5 * density * args_ptr.lp1 +
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
            y3 * (two_shear_r_inv * (2. * args_ptr.llp1 - 1.)) +
            y4 * -3. +
            y5 * -density
        )

    cdef double complex dy5 = \
        y1 * grav_term + \
        y5 * -args_ptr.lp1 * r_inverse + \
        y6

    cdef double complex dy6 = \
        r_inverse * (
            y1 * grav_term * args_ptr.lm1 +
            y6 * args_ptr.lm1 +
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


cdef void cf_liquid_dynamic_compressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        const void* void_args_ptr,
        PreEvalFunc unused) noexcept nogil:

    # Recast the additional arguments
    cdef RadialSolverDiffeqArgStruct* args_ptr = <RadialSolverDiffeqArgStruct*>void_args_ptr

    # Update Equation of State at this radius value
    cdef double[7] eos_array 
    cdef double* eos_array_ptr = &eos_array[0]
    # Call equation of state solution for this layer.
    args_ptr.eos_solution_sptr.get().call(radius, eos_array_ptr)
    # Pull out results.
    # The EOS stores 7 doubles:
    #   0: Gravity
    #   1: Pressure
    #   2: Density
    #   3: Shear Mod (real)
    #   4: Shear Mod (imag)
    #   5: Bulk Mod (real)
    #   6: Bulk Mod (imag)
    cdef double gravity                  = eos_array_ptr[0]
    cdef double density                  = eos_array_ptr[2]
    cdef double complex shear_modulus    = cf_build_dblcmplx(eos_array_ptr[3], eos_array_ptr[4])
    cdef double complex bulk_modulus     = cf_build_dblcmplx(eos_array_ptr[5], eos_array_ptr[6])

    # Pull out y values
    # For the dynamic version, y4 = 0 always in a liquid layer and y3 is defined by y1, y2, and y5 analytically
    cdef double complex y1 = cf_build_dblcmplx(y_ptr[0], y_ptr[1])
    cdef double complex y2 = cf_build_dblcmplx(y_ptr[2], y_ptr[3])
    cdef double complex y5 = cf_build_dblcmplx(y_ptr[4], y_ptr[5])
    cdef double complex y6 = cf_build_dblcmplx(y_ptr[6], y_ptr[7])

    # Optimizations
    cdef double r_inverse         = 1. / radius
    cdef double density_gravity   = density * gravity
    cdef double f2                = args_ptr.frequency * args_ptr.frequency
    cdef double dynamic_term_no_r = -f2 * density
    cdef double dynamic_term      = dynamic_term_no_r * radius
    cdef double grav_term         = args_ptr.grav_coeff * density

    # Check if dynamic term is close to zero. It will always be negative so compare to negative eps
    # if dynamic_term < EPS_100:
    #     # TODO: is faking this okay?
    #     dynamic_term = EPS_100

    # For the liquid layer it is assumed that the shear modulus is zero so the lame parameter simply
    #    equals the bulk modulus. Until bulk dissipation is considered, it will always be real-valued
    cdef double complex lame_inverse = 1. / bulk_modulus

    # y3 derivative is undetermined for a liquid layer, but we can calculate its value which is still used in the
    #   other derivatives.
    y3 = (1. / dynamic_term) * (y2 - density_gravity * y1 + density * y5)
    y1_y3_term = 2. * y1 - args_ptr.llp1 * y3

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
        y5 * density * args_ptr.lp1 * r_inverse - \
        y6 * density - \
        y1_y3_term * density_gravity * r_inverse

    cdef double complex dy5 = \
        y1 * grav_term - \
        y5 * args_ptr.lp1 * r_inverse + \
        y6

    cdef double complex dy6 = \
        r_inverse * (
            args_ptr.lm1 * (y1 * grav_term + y6) +
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
        const void* void_args_ptr,
        PreEvalFunc unused) noexcept nogil:

    # Recast the additional arguments
    cdef RadialSolverDiffeqArgStruct* args_ptr = <RadialSolverDiffeqArgStruct*>void_args_ptr

    # Update Equation of State at this radius value
    cdef double[7] eos_array 
    cdef double* eos_array_ptr = &eos_array[0]
    # Call equation of state solution for this layer.
    args_ptr.eos_solution_sptr.get().call(radius, eos_array_ptr)
    # Pull out results.
    # The EOS stores 7 doubles:
    #   0: Gravity
    #   1: Pressure
    #   2: Density
    #   3: Shear Mod (real)
    #   4: Shear Mod (imag)
    #   5: Bulk Mod (real)
    #   6: Bulk Mod (imag)
    cdef double gravity                  = eos_array_ptr[0]
    cdef double density                  = eos_array_ptr[2]
    cdef double complex shear_modulus    = cf_build_dblcmplx(eos_array_ptr[3], eos_array_ptr[4])
    cdef double complex bulk_modulus     = cf_build_dblcmplx(eos_array_ptr[5], eos_array_ptr[6])

    # Pull out y values
    # For the dynamic version, y4 = 0 always in a liquid layer and y3 is defined by y1, y2, and y5 analytically
    cdef double complex y1 = cf_build_dblcmplx(y_ptr[0], y_ptr[1])
    cdef double complex y2 = cf_build_dblcmplx(y_ptr[2], y_ptr[3])
    cdef double complex y5 = cf_build_dblcmplx(y_ptr[4], y_ptr[5])
    cdef double complex y6 = cf_build_dblcmplx(y_ptr[6], y_ptr[7])

    # Optimizations
    cdef double r_inverse       = 1. / radius
    cdef double density_gravity = density * gravity
    cdef double dynamic_term    = -args_ptr.frequency * args_ptr.frequency * density * radius
    cdef double grav_term       = args_ptr.grav_coeff * density

    # Check if dynamic term is close to zero. It will always be negative so compare to negative eps
    # if dynamic_term < EPS_100:
    #     # TODO: is faking this okay?
    #     dynamic_term = EPS_100

    # y3 derivative is undetermined for a liquid layer, but we can calculate its value which is still used in the
    #   other derivatives.
    cdef double complex y3 = (1. / dynamic_term) * (y2 + density * y5 - density_gravity * y1)
    cdef double complex y1_y3_term = 2. * y1 - args_ptr.llp1 * y3

    cdef double complex dy1 = \
        y1_y3_term * -r_inverse

    cdef double complex dy2 = \
        r_inverse * (
            y1 * (dynamic_term - 2. * density_gravity) +
            y5 * density * args_ptr.lp1 +
            y6 * -density * radius +
            # TODO: In the solid version there is a [2. * (lame + shear_modulus) * r_inverse] coefficient for y1_y3_term
            #   In TS72 the first term is gone. Shouldn't Lame + mu = Lame = Bulk for liquid layer?
            y1_y3_term * -density_gravity
        )

    cdef double complex dy5 = \
        y1 * grav_term + \
        y5 * -args_ptr.lp1 * r_inverse + \
        y6

    cdef double complex dy6 = \
        r_inverse * (
            y1 * grav_term * args_ptr.lm1 +
            y6 * args_ptr.lm1 +
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
        const void* void_args_ptr,
        PreEvalFunc unused) noexcept nogil:

    # Recast the additional arguments
    cdef RadialSolverDiffeqArgStruct* args_ptr = <RadialSolverDiffeqArgStruct*>void_args_ptr

    # Update Equation of State at this radius value
    cdef double[7] eos_array 
    cdef double* eos_array_ptr = &eos_array[0]
    # Call equation of state solution for this layer.
    args_ptr.eos_solution_sptr.get().call(radius, eos_array_ptr)
    # Pull out results.
    # The EOS stores 7 doubles:
    #   0: Gravity
    #   1: Pressure
    #   2: Density
    #   3: Shear Mod (real)
    #   4: Shear Mod (imag)
    #   5: Bulk Mod (real)
    #   6: Bulk Mod (imag)
    cdef double gravity                  = eos_array_ptr[0]
    cdef double density                  = eos_array_ptr[2]
    cdef double complex shear_modulus    = cf_build_dblcmplx(eos_array_ptr[3], eos_array_ptr[4])
    cdef double complex bulk_modulus     = cf_build_dblcmplx(eos_array_ptr[5], eos_array_ptr[6])

    # For the static liquid version, only y5 and y7 are defined.
    cdef double complex y5 = cf_build_dblcmplx(y_ptr[0], y_ptr[1])
    cdef double complex y7 = cf_build_dblcmplx(y_ptr[2], y_ptr[3])

    # Optimizations
    cdef double r_inverse = 1. / radius
    cdef double grav_term = args_ptr.grav_coeff * density / gravity

    # See Eq. 18 in S75
    cdef double complex dy5 = \
        y5 * (grav_term - args_ptr.lp1 * r_inverse) + \
        y7

    cdef double complex dy7 = \
        y5 * 2. * args_ptr.lm1 * r_inverse * grav_term + \
        y7 * (args_ptr.lm1 * r_inverse - grav_term)

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
