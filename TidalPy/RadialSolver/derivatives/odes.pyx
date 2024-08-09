# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from CyRK.cy.cysolverNew cimport PreEvalFunc

from TidalPy.utilities.constants_x cimport PI_DBL
from TidalPy.utilities.math.complex cimport cf_build_dblcmplx
from TidalPy.RadialSolver.derivatives.common cimport ODEInput
from TidalPy.RadialSolver.eos.common cimport EOSOutput


cdef void solid_dynamic_compressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        const void* input_args,
        PreEvalFunc eos_function) noexcept nogil:

    # Pull out additional input arguments by casting the input args to the correct structure
    cdef ODEInput* additional_input_ptr = <ODEInput*>input_args
    cdef double grav_coeff = 4. * PI_DBL * additional_input_ptr.G_to_use
    cdef double lp1  = additional_input_ptr.degree_l + 1.0
    cdef double llp1 = additional_input_ptr.degree_l * lp1
    cdef double lm1  = additional_input_ptr.degree_l - 1.0

    # This function depends on both shear and bulk modulus so tell the EOS it needs to update both.
    additional_input_ptr.update_bulk  = True
    additional_input_ptr.update_shear = True

    # Tell the EOS function where to find the gravity and pressure data in the y values
    additional_input_ptr.gravity_index  = 12
    additional_input_ptr.pressure_index = 13

    # Pull out y information
    # For the static liquid version, only y5 and y7 are defined.
    cdef double complex y1 = cf_build_dblcmplx(y_ptr[0], y_ptr[1])
    cdef double complex y2 = cf_build_dblcmplx(y_ptr[2], y_ptr[3])
    cdef double complex y3 = cf_build_dblcmplx(y_ptr[4], y_ptr[5])
    cdef double complex y4 = cf_build_dblcmplx(y_ptr[6], y_ptr[7])
    cdef double complex y5 = cf_build_dblcmplx(y_ptr[8], y_ptr[9])
    cdef double complex y6 = cf_build_dblcmplx(y_ptr[10], y_ptr[11])

    # Pull out other diffeq parameters
    cdef double gravity = y_ptr[12]

    # Update viscoelastic parameters using the user-provided equation of state
    cdef EOSOutput eos_output
    cdef EOSOutput* eos_output_ptr = &eos_output 
    eos_function(eos_output_ptr, radius, y_ptr, input_args)

    # Convert compressibility parameters (the first lame parameter can be complex)
    cdef double complex lame
    lame = (eos_output.bulk_modulus - (2. / 3.) * eos_output.shear_modulus)

    # Optimizations
    cdef double r_inverse, density_gravity, dynamic_term, grav_term
    cdef double complex lame_2mu, lame_2mu_inverse, two_shear_r_inv, y1_y3_term
    lame_2mu         = lame + 2. * eos_output.shear_modulus
    lame_2mu_inverse = 1. / lame_2mu
    r_inverse        = 1. / radius
    two_shear_r_inv  = 2. * eos_output.shear_modulus * r_inverse
    density_gravity  = eos_output.density * gravity
    dynamic_term     = -additional_input_ptr.frequency * additional_input_ptr.frequency * eos_output.density * radius
    grav_term        = grav_coeff * eos_output.density
    y1_y3_term       = 2. * y1 - llp1 * y3

    # See Eq. 82 in TS72 or Eqs. 4--9 in KMN15 or Eqs. 13--18 in B15
    #   Note: There appears to be a missing factor of mu^2 in some of the terms in KMN15.
    # dy2 and dy4 contain all three of: dynamic, viscoelastic, and gravitational terms.
    cdef double complex dy1, dy2, dy3, dy4, dy5, dy6

    dy1 = lame_2mu_inverse * (
            y1_y3_term * -lame * r_inverse +
            y2
    )

    dy2 = r_inverse * (
            y1 * (dynamic_term - 2. * density_gravity) +
            y2 * -2. +
            y4 * llp1 +
            y5 * eos_output.density * lp1 +
            y6 * -eos_output.density * radius +
            dy1 * 2. * lame +
            y1_y3_term * (2. * (lame + eos_output.shear_modulus) * r_inverse - density_gravity)
    )

    dy3 = \
        y1 * -r_inverse + \
        y3 * r_inverse + \
        y4 * (1. / eos_output.shear_modulus)

    dy4 = r_inverse * (
            y1 * (density_gravity + two_shear_r_inv) +
            y3 * (dynamic_term - two_shear_r_inv) +
            y4 * -3. +
            y5 * -eos_output.density +
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

    # Solve for gravity
    dy_ptr[12] = grav_coeff * eos_output.density - 2.0 * gravity * r_inverse

    # Solve for pressure
    dy_ptr[13] = -eos_output.density * gravity


cdef void solid_dynamic_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        const void* input_args,
        PreEvalFunc eos_function) noexcept nogil:

    # Pull out additional input arguments by casting the input args to the correct structure
    cdef ODEInput* additional_input_ptr = <ODEInput*>input_args
    cdef double grav_coeff = 4. * PI_DBL * additional_input_ptr.G_to_use
    cdef double lp1  = additional_input_ptr.degree_l + 1.0
    cdef double llp1 = additional_input_ptr.degree_l * lp1
    cdef double lm1  = additional_input_ptr.degree_l - 1.0

    # This function depends on shear modulus so tell the EOS it needs to update it and not bulk.
    additional_input_ptr.update_bulk  = False
    additional_input_ptr.update_shear = True

    # Tell the EOS function where to find the gravity and pressure data in the y values
    additional_input_ptr.gravity_index  = 12
    additional_input_ptr.pressure_index = 13

    # Pull out y information
    # For the static liquid version, only y5 and y7 are defined.
    cdef double complex y1 = cf_build_dblcmplx(y_ptr[0], y_ptr[1])
    cdef double complex y2 = cf_build_dblcmplx(y_ptr[2], y_ptr[3])
    cdef double complex y3 = cf_build_dblcmplx(y_ptr[4], y_ptr[5])
    cdef double complex y4 = cf_build_dblcmplx(y_ptr[6], y_ptr[7])
    cdef double complex y5 = cf_build_dblcmplx(y_ptr[8], y_ptr[9])
    cdef double complex y6 = cf_build_dblcmplx(y_ptr[10], y_ptr[11])

    # Pull out other diffeq parameters
    cdef double gravity = y_ptr[12]

    # Update viscoelastic parameters using the user-provided equation of state
    cdef EOSOutput eos_output
    cdef EOSOutput* eos_output_ptr = &eos_output 
    eos_function(eos_output_ptr, radius, y_ptr, input_args)

    # Optimizations
    cdef double r_inverse, density_gravity, dynamic_term, grav_term
    cdef double complex two_shear_r_inv, y1_y3_term
    r_inverse       = 1. / radius
    two_shear_r_inv = 2. * eos_output.shear_modulus * r_inverse
    density_gravity = eos_output.density * gravity
    dynamic_term    = -additional_input_ptr.frequency * additional_input_ptr.frequency * eos_output.density * radius
    grav_term       = grav_coeff * eos_output.density
    y1_y3_term      = 2. * y1 - llp1 * y3

    # See Eq. 82 in TS72 or Eqs. 4--9 in KMN15 or Eqs. 13--18 in B15
    #   Note: There appears to be a missing factor of mu^2 in some of the terms in KMN15.
    # dy2 and dy4 contain all three of: dynamic, viscoelastic, and gravitational terms.
    cdef double complex dy1, dy2, dy3, dy4, dy5, dy6

    dy1 = y1_y3_term * -1. * r_inverse

    dy2 = r_inverse * (
            y1 * (dynamic_term + 12. * eos_output.shear_modulus * r_inverse - 4. * density_gravity) +
            y3 * llp1 * (density_gravity - 6. * eos_output.shear_modulus * r_inverse) +
            y4 * llp1 +
            y5 * eos_output.density * lp1 +
            y6 * -eos_output.density * radius
    )

    dy3 = \
        y1 * -r_inverse + \
        y3 * r_inverse + \
        y4 * (1. / eos_output.shear_modulus)

    dy4 = r_inverse * (
            y1 * (density_gravity - 3. * two_shear_r_inv) +
            y2 * -1. +
            y3 * (dynamic_term + two_shear_r_inv * (2. * llp1 - 1.)) +
            y4 * -3. +
            y5 * -eos_output.density
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

    # Solve for gravity
    dy_ptr[12] = grav_coeff * eos_output.density - 2.0 * gravity * r_inverse

    # Solve for pressure
    dy_ptr[13] = -eos_output.density * gravity


cdef void solid_static_compressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        const void* input_args,
        PreEvalFunc eos_function) noexcept nogil:

    # Pull out additional input arguments by casting the input args to the correct structure
    cdef ODEInput* additional_input_ptr = <ODEInput*>input_args
    cdef double grav_coeff = 4. * PI_DBL * additional_input_ptr.G_to_use
    cdef double lp1  = additional_input_ptr.degree_l + 1.0
    cdef double llp1 = additional_input_ptr.degree_l * lp1
    cdef double lm1  = additional_input_ptr.degree_l - 1.0

    # This function depends on shear and bulk modulus so tell the EOS it needs to update both.
    additional_input_ptr.update_bulk  = True
    additional_input_ptr.update_shear = True

    # Tell the EOS function where to find the gravity and pressure data in the y values
    additional_input_ptr.gravity_index  = 12
    additional_input_ptr.pressure_index = 13

    # Pull out y information
    # For the static liquid version, only y5 and y7 are defined.
    cdef double complex y1 = cf_build_dblcmplx(y_ptr[0], y_ptr[1])
    cdef double complex y2 = cf_build_dblcmplx(y_ptr[2], y_ptr[3])
    cdef double complex y3 = cf_build_dblcmplx(y_ptr[4], y_ptr[5])
    cdef double complex y4 = cf_build_dblcmplx(y_ptr[6], y_ptr[7])
    cdef double complex y5 = cf_build_dblcmplx(y_ptr[8], y_ptr[9])
    cdef double complex y6 = cf_build_dblcmplx(y_ptr[10], y_ptr[11])

    # Pull out other diffeq parameters
    cdef double gravity = y_ptr[12]

    # Update viscoelastic parameters using the user-provided equation of state
    cdef EOSOutput eos_output
    cdef EOSOutput* eos_output_ptr = &eos_output 
    eos_function(eos_output_ptr, radius, y_ptr, input_args)

    # Convert compressibility parameters (the first lame parameter can be complex)
    cdef double complex lame
    lame = (eos_output.bulk_modulus - (2. / 3.) * eos_output.shear_modulus)

    # Optimizations
    cdef double r_inverse, density_gravity, grav_term
    cdef double complex lame_2mu, lame_2mu_inverse, two_shear_r_inv, y1_y3_term

    lame_2mu         = lame + 2. * eos_output.shear_modulus
    lame_2mu_inverse = 1. / lame_2mu
    r_inverse        = 1. / radius
    two_shear_r_inv  = 2. * eos_output.shear_modulus * r_inverse
    density_gravity  = eos_output.density * gravity
    grav_term        = grav_coeff * eos_output.density
    y1_y3_term       = 2. * y1 - llp1 * y3

    # See Eq. 82 in TS72 or Eqs. 4--9 in KMN15 or Eqs. 13--18 in B15
    #   Note: There appears to be a missing factor of mu^2 in some of the terms in KMN15.
    # The static case just sets all frequency_to_use dependence in these equations to zero.
    # dy2 and dy4 contain: viscoelastic, and gravitational terms.
    cdef double complex dy1, dy2, dy3, dy4, dy5, dy6

    dy1 = lame_2mu_inverse * (
            y1_y3_term * -lame * r_inverse +
            y2
    )

    dy2 = r_inverse * (
            y1 * -2. * density_gravity +
            y2 * -2. +
            y4 * llp1 +
            y5 * eos_output.density * lp1 +
            y6 * -eos_output.density * radius +
            dy1 * 2. * lame +
            y1_y3_term * (2. * (lame + eos_output.shear_modulus) * r_inverse - density_gravity)
    )

    dy3 = \
        y1 * -r_inverse + \
        y3 * r_inverse + \
        y4 * (1. / eos_output.shear_modulus)

    dy4 = r_inverse * (
            y1 * (density_gravity + two_shear_r_inv) +
            y3 * -two_shear_r_inv +
            y4 * -3. +
            y5 * -eos_output.density +
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

    # Solve for gravity
    dy_ptr[12] = grav_coeff * eos_output.density - 2.0 * gravity * r_inverse

    # Solve for pressure
    dy_ptr[13] = -eos_output.density * gravity


cdef void solid_static_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        const void* input_args,
        PreEvalFunc eos_function) noexcept nogil:

    # Pull out additional input arguments by casting the input args to the correct structure
    cdef ODEInput* additional_input_ptr = <ODEInput*>input_args
    cdef double grav_coeff = 4. * PI_DBL * additional_input_ptr.G_to_use
    cdef double lp1  = additional_input_ptr.degree_l + 1.0
    cdef double llp1 = additional_input_ptr.degree_l * lp1
    cdef double lm1  = additional_input_ptr.degree_l - 1.0

    # This function does not depend on bulk modulus so tell the EOS there is no reason to update it.
    additional_input_ptr.update_bulk  = False
    additional_input_ptr.update_shear = True

    # Tell the EOS function where to find the gravity and pressure data in the y values
    additional_input_ptr.gravity_index  = 12
    additional_input_ptr.pressure_index = 13

    # Pull out y information
    # For the static liquid version, only y5 and y7 are defined.
    cdef double complex y1 = cf_build_dblcmplx(y_ptr[0], y_ptr[1])
    cdef double complex y2 = cf_build_dblcmplx(y_ptr[2], y_ptr[3])
    cdef double complex y3 = cf_build_dblcmplx(y_ptr[4], y_ptr[5])
    cdef double complex y4 = cf_build_dblcmplx(y_ptr[6], y_ptr[7])
    cdef double complex y5 = cf_build_dblcmplx(y_ptr[8], y_ptr[9])
    cdef double complex y6 = cf_build_dblcmplx(y_ptr[10], y_ptr[11])

    # Pull out other diffeq parameters
    cdef double gravity = y_ptr[12]

    # Update viscoelastic parameters using the user-provided equation of state
    cdef EOSOutput eos_output
    cdef EOSOutput* eos_output_ptr = &eos_output 
    eos_function(eos_output_ptr, radius, y_ptr, input_args)

    # Optimizations
    cdef double r_inverse, density_gravity, grav_term
    cdef double complex two_shear_r_inv, y1_y3_term

    r_inverse       = 1. / radius
    two_shear_r_inv = 2. * eos_output.shear_modulus * r_inverse
    density_gravity = eos_output.density * gravity
    grav_term       = grav_coeff * eos_output.density
    y1_y3_term      = 2. * y1 - llp1 * y3

    cdef double complex dy1, dy2, dy3, dy4, dy5, dy6

    dy1 = y1_y3_term * -1. * r_inverse

    dy2 = r_inverse * (
            y1 * (12. * eos_output.shear_modulus * r_inverse - 4. * density_gravity) +
            y3 * llp1 * (density_gravity - 6. * eos_output.shear_modulus * r_inverse) +
            y4 * llp1 +
            y5 * eos_output.density * lp1 +
            y6 * -eos_output.density * radius
    )

    dy3 = \
        y1 * -r_inverse + \
        y3 * r_inverse + \
        y4 * (1. / eos_output.shear_modulus)

    dy4 = r_inverse * (
            y1 * (density_gravity - 3. * two_shear_r_inv) +
            y2 * -1. +
            y3 * (two_shear_r_inv * (2. * llp1 - 1.)) +
            y4 * -3. +
            y5 * -eos_output.density
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

    # Solve for gravity
    dy_ptr[12] = grav_coeff * eos_output.density - 2.0 * gravity * r_inverse

    # Solve for pressure
    dy_ptr[13] = -eos_output.density * gravity



cdef void liquid_dynamic_compressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        const void* input_args,
        PreEvalFunc eos_function) noexcept nogil:

    # Pull out additional input arguments by casting the input args to the correct structure
    cdef ODEInput* additional_input_ptr = <ODEInput*>input_args
    cdef double grav_coeff = 4. * PI_DBL * additional_input_ptr.G_to_use
    cdef double lp1  = additional_input_ptr.degree_l + 1.0
    cdef double llp1 = additional_input_ptr.degree_l * lp1
    cdef double lm1  = additional_input_ptr.degree_l - 1.0

    # This function does not depend on shear modulus so tell the EOS there is no reason to update it.
    additional_input_ptr.update_bulk  = True
    additional_input_ptr.update_shear = False

    # Tell the EOS function where to find the gravity and pressure data in the y values
    additional_input_ptr.gravity_index  = 8
    additional_input_ptr.pressure_index = 9

    # Pull out y information
    # For the static liquid version, only y5 and y7 are defined.
    cdef double complex y1 = cf_build_dblcmplx(y_ptr[0], y_ptr[1])
    cdef double complex y2 = cf_build_dblcmplx(y_ptr[2], y_ptr[3])
    cdef double complex y5 = cf_build_dblcmplx(y_ptr[4], y_ptr[5])
    cdef double complex y6 = cf_build_dblcmplx(y_ptr[6], y_ptr[7])

    # Pull out other diffeq parameters
    cdef double gravity = y_ptr[8]

    # Update viscoelastic parameters using the user-provided equation of state
    cdef EOSOutput eos_output
    cdef EOSOutput* eos_output_ptr = &eos_output 
    eos_function(eos_output_ptr, radius, y_ptr, input_args)

    # Optimizations
    cdef double r_inverse, density_gravity, f2, dynamic_term, dynamic_term_no_r, grav_term

    r_inverse         = 1. / radius
    density_gravity   = eos_output.density * gravity
    f2                = additional_input_ptr.frequency * additional_input_ptr.frequency
    dynamic_term_no_r = -f2 * eos_output.density
    dynamic_term      = dynamic_term_no_r * radius
    grav_term         = grav_coeff * eos_output.density

    # Check if dynamic term is close to zero. It will always be negative so compare to negative eps
    # if dynamic_term < EPS_100:
    #     # TODO: is faking this okay?
    #     dynamic_term = EPS_100

    # Until bulk dissipation is considered, lame_inverse will always be real-valued for a liquid layer.
    cdef double lame_inverse

    # For the liquid layer it is assumed that the shear modulus is zero so the lame parameter simply
    #    equals the bulk modulus. Until bulk dissipation is considered, it will always be real-valued
    lame_inverse = 1. / eos_output.bulk_modulus

    # y3 derivative is undetermined for a liquid layer, but we can calculate its value which is still used in the
    #   other derivatives.
    y3 = (1. / dynamic_term) * (y2 - density_gravity * y1 + eos_output.density * y5)
    y1_y3_term = 2. * y1 - llp1 * y3

    # Eqs. 11--14 in KMN15 equations look like they don't match TS72 because they applied the rheology already.
    #    and substituted y3.
    # We will use TS72 eq. 87 to allow for a generic rheology and bulk dissipation.
    # dy2 contain all three of: dynamic, viscoelastic, and gravitational terms.
    cdef double complex dy1, dy2, dy5, dy6

    dy1 = \
        y2 * lame_inverse - \
        y1_y3_term * r_inverse

    # TODO: In the solid version there is a [2. * (lame + shear_modulus) * r_inverse] coefficient for y1_y3_term
    #   In TS72 the first term is gone. Shouldn't Lame + mu = Lame = Bulk for liquid layer?
    dy2 = \
        y1 * (dynamic_term_no_r - 2. * density_gravity * r_inverse) + \
        y5 * eos_output.density * lp1 * r_inverse - \
        y6 * eos_output.density - \
        y1_y3_term * density_gravity * r_inverse

    dy5 = \
        y1 * grav_term - \
        y5 * lp1 * r_inverse + \
        y6

    dy6 = r_inverse * (
            lm1 * (y1 * grav_term + y6) +
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

    # Solve for gravity
    dy_ptr[8] = grav_coeff * eos_output.density - 2.0 * gravity * r_inverse

    # Solve for pressure
    dy_ptr[9] = -eos_output.density * gravity


cdef void liquid_dynamic_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        const void* input_args,
        PreEvalFunc eos_function) noexcept nogil:

    # Pull out additional input arguments by casting the input args to the correct structure
    cdef ODEInput* additional_input_ptr = <ODEInput*>input_args
    cdef double grav_coeff = 4. * PI_DBL * additional_input_ptr.G_to_use
    cdef double lp1  = additional_input_ptr.degree_l + 1.0
    cdef double llp1 = additional_input_ptr.degree_l * lp1
    cdef double lm1  = additional_input_ptr.degree_l - 1.0

    # This function does not depend on bulk or shear modulus so tell the EOS there is no reason to update it.
    additional_input_ptr.update_bulk  = False
    additional_input_ptr.update_shear = False

    # Tell the EOS function where to find the gravity and pressure data in the y values
    additional_input_ptr.gravity_index  = 8
    additional_input_ptr.pressure_index = 9

    # Pull out y information
    # For the static liquid version, only y5 and y7 are defined.
    cdef double complex y1 = cf_build_dblcmplx(y_ptr[0], y_ptr[1])
    cdef double complex y2 = cf_build_dblcmplx(y_ptr[2], y_ptr[3])
    cdef double complex y5 = cf_build_dblcmplx(y_ptr[4], y_ptr[5])
    cdef double complex y6 = cf_build_dblcmplx(y_ptr[6], y_ptr[7])

    # Pull out other diffeq parameters
    cdef double gravity = y_ptr[8]

    # Update viscoelastic parameters using the user-provided equation of state
    cdef EOSOutput eos_output
    cdef EOSOutput* eos_output_ptr = &eos_output 
    eos_function(eos_output_ptr, radius, y_ptr, input_args)

    # Optimizations
    cdef double r_inverse, density_gravity, dynamic_term, grav_term

    r_inverse       = 1. / radius
    density_gravity = eos_output.density * gravity
    dynamic_term    = -additional_input_ptr.frequency * additional_input_ptr.frequency * eos_output.density * radius
    grav_term       = grav_coeff * eos_output.density

    # Check if dynamic term is close to zero. It will always be negative so compare to negative eps
    # if dynamic_term < EPS_100:
    #     # TODO: is faking this okay?
    #     dynamic_term = EPS_100

    # y3 derivative is undetermined for a liquid layer, but we can calculate its value which is still used in the
    #   other derivatives.
    cdef double complex y3, y1_y3_term
    y3         = (1. / dynamic_term) * (y2 + eos_output.density * y5 - density_gravity * y1)
    y1_y3_term = 2. * y1 - llp1 * y3

    cdef double complex dy1, dy2, dy5, dy6

    dy1 = y1_y3_term * -r_inverse

    dy2 = r_inverse * (
            y1 * (dynamic_term - 2. * density_gravity) +
            y5 * eos_output.density * lp1 +
            y6 * -eos_output.density * radius +
            # TODO: In the solid version there is a [2. * (lame + shear_modulus) * r_inverse] coefficient for y1_y3_term
            #   In TS72 the first term is gone. Shouldn't Lame + mu = Lame = Bulk for liquid layer?
            y1_y3_term * -density_gravity
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
    dy_ptr[0] = dy1.real
    dy_ptr[1] = dy1.imag
    dy_ptr[2] = dy2.real
    dy_ptr[3] = dy2.imag
    dy_ptr[4] = dy5.real
    dy_ptr[5] = dy5.imag
    dy_ptr[6] = dy6.real
    dy_ptr[7] = dy6.imag

    # Solve for gravity
    dy_ptr[8] = grav_coeff * eos_output.density - 2.0 * gravity * r_inverse

    # Solve for pressure
    dy_ptr[9] = -eos_output.density * gravity


cdef void liquid_static_compressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        const void* input_args,
        PreEvalFunc eos_function) noexcept nogil:

    # Using the static assumption results in the same function for liquid compressible and incompressible.
    # TODO: Check that this is actually valid.
    liquid_static_incompressible(dy_ptr, radius, y_ptr, input_args, eos_function)

cdef void liquid_static_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        const void* input_args,
        PreEvalFunc eos_function) noexcept nogil:

    # Pull out additional input arguments by casting the input args to the correct structure
    cdef ODEInput* additional_input_ptr = <ODEInput*>input_args
    cdef double grav_coeff = 4. * PI_DBL * additional_input_ptr.G_to_use
    cdef double lp1 = additional_input_ptr.degree_l + 1.0
    cdef double lm1 = additional_input_ptr.degree_l - 1.0

    # This function does not depend on bulk or shear modulus so tell the EOS there is no reason to update it.
    additional_input_ptr.update_bulk  = False
    additional_input_ptr.update_shear = False

    # Tell the EOS function where to find the gravity and pressure data in the y values
    additional_input_ptr.gravity_index  = 4
    additional_input_ptr.pressure_index = 5

    # Pull out y information
    # For the static liquid version, only y5 and y7 are defined.
    cdef double complex y5 = cf_build_dblcmplx(y_ptr[0], y_ptr[1])
    cdef double complex y7 = cf_build_dblcmplx(y_ptr[2], y_ptr[3])

    # Pull out other diffeq parameters
    cdef double gravity = y_ptr[4]

    # Update viscoelastic parameters using the user-provided equation of state
    cdef EOSOutput eos_output
    cdef EOSOutput* eos_output_ptr = &eos_output 
    eos_function(eos_output_ptr, radius, y_ptr, input_args)

    # Optimizations
    cdef double r_inverse, grav_term

    r_inverse = 1. / radius
    grav_term = grav_coeff * eos_output.density / gravity

    # See Eq. 18 in S75
    cdef double complex dy5, dy7

    dy5 = \
        y5 * (grav_term - lp1 * r_inverse) + \
        y7

    dy7 = \
        y5 * 2. * lm1 * r_inverse * grav_term + \
        y7 * (lm1 * r_inverse - grav_term)

    # Convert back to floats
    dy_ptr[0] = dy5.real
    dy_ptr[1] = dy5.imag
    dy_ptr[2] = dy7.real
    dy_ptr[3] = dy7.imag

    # Solve for gravity
    dy_ptr[4] = grav_coeff * eos_output.density - 2.0 * gravity * r_inverse

    # Solve for pressure
    dy_ptr[5] = -eos_output.density * gravity
