# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

import numpy as np

from TidalPy.radial_solver.numerical.derivatives.ode_base_x cimport RadialSolverBase

cdef double EPS = np.finfo(np.float64).eps
cdef double EPS_10 = EPS * 10.
cdef double EPS_100 = EPS * 100.


cdef class SolidDynamicCompressible(RadialSolverBase):
    cdef void diffeq(self) noexcept nogil:
        """

        References
        ----------
        KMN15; B15; TS72
        """

        # Note that `t_now` is the current "radius" not time.
        cdef double radius
        radius = self.t_now

        # Update interpolation
        self.update_interp(update_bulk=True, update_shear=True)

        # Pull out y values
        cdef double y1_real, y2_real, y3_real, y4_real, y5_real, y6_real
        cdef double y1_imag, y2_imag, y3_imag, y4_imag, y5_imag, y6_imag
        y1_real = self.y_ptr[0]
        y1_imag = self.y_ptr[1]
        y2_real = self.y_ptr[2]
        y2_imag = self.y_ptr[3]
        y3_real = self.y_ptr[4]
        y3_imag = self.y_ptr[5]
        y4_real = self.y_ptr[6]
        y4_imag = self.y_ptr[7]
        y5_real = self.y_ptr[8]
        y5_imag = self.y_ptr[9]
        y6_real = self.y_ptr[10]
        y6_imag = self.y_ptr[11]

        # Convert floats to complex
        cdef double complex y1, y2, y3, y4, y5, y6
        y1 = y1_real + 1.0j * y1_imag
        y2 = y2_real + 1.0j * y2_imag
        y3 = y3_real + 1.0j * y3_imag
        y4 = y4_real + 1.0j * y4_imag
        y5 = y5_real + 1.0j * y5_imag
        y6 = y6_real + 1.0j * y6_imag

        # Convert compressibility parameters (the first lame parameter can be complex)
        cdef double complex lame
        lame = (<double complex> self.bulk_modulus - (2. / 3.) * self.shear_modulus)

        # Optimizations
        cdef double r_inverse, density_gravity, dynamic_term, grav_term
        cdef double complex lame_2mu, lame_2mu_inverse, two_shear_r_inv, y1_y3_term
        lame_2mu         = lame + 2. * self.shear_modulus
        lame_2mu_inverse = 1. / lame_2mu
        r_inverse        = 1. / radius
        two_shear_r_inv  = 2. * self.shear_modulus * r_inverse
        density_gravity  = self.density * self.gravity
        dynamic_term     = -self.frequency * self.frequency * self.density * radius
        grav_term        = self.grav_coeff * self.density
        y1_y3_term       = 2. * y1 - self.llp1 * y3

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
                y4 * self.llp1 +
                y5 * self.density * self.lp1 +
                y6 * -self.density * radius +
                dy1 * 2. * lame +
                y1_y3_term * (2. * (lame + self.shear_modulus) * r_inverse - density_gravity)
        )

        dy3 = \
            y1 * -r_inverse + \
            y3 * r_inverse + \
            y4 * (1. / self.shear_modulus)

        dy4 = r_inverse * (
                y1 * (density_gravity + two_shear_r_inv) +
                y3 * (dynamic_term - two_shear_r_inv) +
                y4 * -3. +
                y5 * -self.density +
                dy1 * -lame +
                y1_y3_term * -lame_2mu * r_inverse
        )

        dy5 = \
            y1 * grav_term + \
            y5 * -self.lp1 * r_inverse + \
            y6

        dy6 = r_inverse * (
                y1 * grav_term * self.lm1 +
                y6 * self.lm1 +
                y1_y3_term * grav_term
        )

        # Convert back to floats
        self.dy_ptr[0]  = dy1.real
        self.dy_ptr[1]  = dy1.imag
        self.dy_ptr[2]  = dy2.real
        self.dy_ptr[3]  = dy2.imag
        self.dy_ptr[4]  = dy3.real
        self.dy_ptr[5]  = dy3.imag
        self.dy_ptr[6]  = dy4.real
        self.dy_ptr[7]  = dy4.imag
        self.dy_ptr[8]  = dy5.real
        self.dy_ptr[9]  = dy5.imag
        self.dy_ptr[10] = dy6.real
        self.dy_ptr[11] = dy6.imag


cdef class SolidDynamicIncompressible(RadialSolverBase):

    cdef void diffeq(self) noexcept nogil:
        """

        References
        ----------
        KMN15; B15; TS72
        """

        # Note that `t_now` is the current "radius" not time.
        cdef double radius
        radius = self.t_now

        # Update interpolation
        self.update_interp(update_bulk=True, update_shear=True)

        # Pull out y values
        cdef double y1_real, y2_real, y3_real, y4_real, y5_real, y6_real
        cdef double y1_imag, y2_imag, y3_imag, y4_imag, y5_imag, y6_imag
        y1_real = self.y_ptr[0]
        y1_imag = self.y_ptr[1]
        y2_real = self.y_ptr[2]
        y2_imag = self.y_ptr[3]
        y3_real = self.y_ptr[4]
        y3_imag = self.y_ptr[5]
        y4_real = self.y_ptr[6]
        y4_imag = self.y_ptr[7]
        y5_real = self.y_ptr[8]
        y5_imag = self.y_ptr[9]
        y6_real = self.y_ptr[10]
        y6_imag = self.y_ptr[11]

        # Convert floats to complex
        cdef double complex y1, y2, y3, y4, y5, y6
        y1 = y1_real + 1.0j * y1_imag
        y2 = y2_real + 1.0j * y2_imag
        y3 = y3_real + 1.0j * y3_imag
        y4 = y4_real + 1.0j * y4_imag
        y5 = y5_real + 1.0j * y5_imag
        y6 = y6_real + 1.0j * y6_imag

        # Optimizations
        cdef double r_inverse, density_gravity, dynamic_term, grav_term
        cdef double complex two_shear_r_inv, y1_y3_term
        r_inverse       = 1. / radius
        two_shear_r_inv = 2. * self.shear_modulus * r_inverse
        density_gravity = self.density * self.gravity
        dynamic_term    = -self.frequency * self.frequency * self.density * radius
        grav_term       = self.grav_coeff * self.density
        y1_y3_term      = 2. * y1 - self.llp1 * y3

        # See Eq. 82 in TS72 or Eqs. 4--9 in KMN15 or Eqs. 13--18 in B15
        #   Note: There appears to be a missing factor of mu^2 in some of the terms in KMN15.
        # dy2 and dy4 contain all three of: dynamic, viscoelastic, and gravitational terms.
        cdef double complex dy1, dy2, dy3, dy4, dy5, dy6

        dy1 = y1_y3_term * -1. * r_inverse

        dy2 = r_inverse * (
                y1 * (dynamic_term + 12. * self.shear_modulus * r_inverse - 4. * density_gravity) +
                y3 * self.llp1 * (density_gravity - 6. * self.shear_modulus * r_inverse) +
                y4 * self.llp1 +
                y5 * self.density * self.lp1 +
                y6 * -self.density * radius
        )

        dy3 = \
            y1 * -r_inverse + \
            y3 * r_inverse + \
            y4 * (1. / self.shear_modulus)

        dy4 = r_inverse * (
                y1 * (density_gravity - 3. * two_shear_r_inv) +
                y2 * -1. +
                y3 * (dynamic_term + two_shear_r_inv * (2. * self.llp1 - 1.)) +
                y4 * -3. +
                y5 * -self.density
        )

        dy5 = \
            y1 * grav_term + \
            y5 * -self.lp1 * r_inverse + \
            y6

        dy6 = r_inverse * (
                y1 * grav_term * self.lm1 +
                y6 * self.lm1 +
                y1_y3_term * grav_term
        )

        # Convert back to floats
        self.dy_ptr[0]  = dy1.real
        self.dy_ptr[1]  = dy1.imag
        self.dy_ptr[2]  = dy2.real
        self.dy_ptr[3]  = dy2.imag
        self.dy_ptr[4]  = dy3.real
        self.dy_ptr[5]  = dy3.imag
        self.dy_ptr[6]  = dy4.real
        self.dy_ptr[7]  = dy4.imag
        self.dy_ptr[8]  = dy5.real
        self.dy_ptr[9]  = dy5.imag
        self.dy_ptr[10] = dy6.real
        self.dy_ptr[11] = dy6.imag


cdef class SolidStaticCompressible(RadialSolverBase):
    cdef void diffeq(self) noexcept nogil:
        # Note that `t_now` is the current "radius" not time.
        cdef double radius
        radius = self.t_now

        # Update interpolation
        self.update_interp(update_bulk=True, update_shear=True)

        # Pull out y values
        cdef double y1_real, y2_real, y3_real, y4_real, y5_real, y6_real
        cdef double y1_imag, y2_imag, y3_imag, y4_imag, y5_imag, y6_imag

        y1_real = self.y_ptr[0]
        y1_imag = self.y_ptr[1]
        y2_real = self.y_ptr[2]
        y2_imag = self.y_ptr[3]
        y3_real = self.y_ptr[4]
        y3_imag = self.y_ptr[5]
        y4_real = self.y_ptr[6]
        y4_imag = self.y_ptr[7]
        y5_real = self.y_ptr[8]
        y5_imag = self.y_ptr[9]
        y6_real = self.y_ptr[10]
        y6_imag = self.y_ptr[11]

        # Convert floats to complex
        cdef double complex y1, y2, y3, y4, y5, y6

        y1 = y1_real + 1.0j * y1_imag
        y2 = y2_real + 1.0j * y2_imag
        y3 = y3_real + 1.0j * y3_imag
        y4 = y4_real + 1.0j * y4_imag
        y5 = y5_real + 1.0j * y5_imag
        y6 = y6_real + 1.0j * y6_imag

        # Convert compressibility parameters (the first lame parameter can be complex)
        cdef double complex lame
        lame = (<double complex>self.bulk_modulus - (2. / 3.) * self.shear_modulus)

        # Optimizations
        cdef double r_inverse, density_gravity, grav_term
        cdef double complex lame_2mu, lame_2mu_inverse, two_shear_r_inv, y1_y3_term

        lame_2mu         = lame + 2. * self.shear_modulus
        lame_2mu_inverse = 1. / lame_2mu
        r_inverse        = 1. / radius
        two_shear_r_inv  = 2. * self.shear_modulus * r_inverse
        density_gravity  = self.density * self.gravity
        grav_term        = self.grav_coeff * self.density
        y1_y3_term       = 2. * y1 - self.llp1 * y3

        # See Eq. 82 in TS72 or Eqs. 4--9 in KMN15 or Eqs. 13--18 in B15
        #   Note: There appears to be a missing factor of mu^2 in some of the terms in KMN15.
        # The static case just sets all frequency dependence in these equations to zero.
        # dy2 and dy4 contain: viscoelastic, and gravitational terms.
        cdef double complex dy1, dy2, dy3, dy4, dy5, dy6

        dy1 = lame_2mu_inverse * (
                y1_y3_term * -lame * r_inverse +
                y2
        )

        dy2 = r_inverse * (
                y1 * -2. * density_gravity +
                y2 * -2. +
                y4 * self.llp1 +
                y5 * self.density * self.lp1 +
                y6 * -self.density * radius +
                dy1 * 2. * lame +
                y1_y3_term * (2. * (lame + self.shear_modulus) * r_inverse - density_gravity)
        )

        dy3 = \
            y1 * -r_inverse + \
            y3 * r_inverse + \
            y4 * (1. / self.shear_modulus)

        dy4 = r_inverse * (
                y1 * (density_gravity + two_shear_r_inv) +
                y3 * -two_shear_r_inv +
                y4 * -3. +
                y5 * -self.density +
                dy1 * -lame +
                y1_y3_term * -lame_2mu * r_inverse
        )

        dy5 = \
            y1 * grav_term + \
            y5 * -self.lp1 * r_inverse + \
            y6

        dy6 = r_inverse * (
                y1 * grav_term * self.lm1 +
                y6 * self.lm1 +
                y1_y3_term * grav_term
        )

        # Convert back to floats
        self.dy_ptr[0]  = dy1.real
        self.dy_ptr[1]  = dy1.imag
        self.dy_ptr[2]  = dy2.real
        self.dy_ptr[3]  = dy2.imag
        self.dy_ptr[4]  = dy3.real
        self.dy_ptr[5]  = dy3.imag
        self.dy_ptr[6]  = dy4.real
        self.dy_ptr[7]  = dy4.imag
        self.dy_ptr[8]  = dy5.real
        self.dy_ptr[9]  = dy5.imag
        self.dy_ptr[10] = dy6.real
        self.dy_ptr[11] = dy6.imag


cdef class SolidStaticIncompressible(RadialSolverBase):
    cdef void diffeq(self) noexcept nogil:
        # Note that `t_now` is the current "radius" not time.
        cdef double radius
        radius = self.t_now

        # Update interpolation
        self.update_interp(update_bulk=False, update_shear=True)

        # Pull out y values
        cdef double y1_real, y2_real, y3_real, y4_real, y5_real, y6_real
        cdef double y1_imag, y2_imag, y3_imag, y4_imag, y5_imag, y6_imag

        y1_real = self.y_ptr[0]
        y1_imag = self.y_ptr[1]
        y2_real = self.y_ptr[2]
        y2_imag = self.y_ptr[3]
        y3_real = self.y_ptr[4]
        y3_imag = self.y_ptr[5]
        y4_real = self.y_ptr[6]
        y4_imag = self.y_ptr[7]
        y5_real = self.y_ptr[8]
        y5_imag = self.y_ptr[9]
        y6_real = self.y_ptr[10]
        y6_imag = self.y_ptr[11]

        # Convert floats to complex
        cdef double complex y1, y2, y3, y4, y5, y6

        y1 = y1_real + 1.0j * y1_imag
        y2 = y2_real + 1.0j * y2_imag
        y3 = y3_real + 1.0j * y3_imag
        y4 = y4_real + 1.0j * y4_imag
        y5 = y5_real + 1.0j * y5_imag
        y6 = y6_real + 1.0j * y6_imag

        # Optimizations
        cdef double r_inverse, density_gravity, grav_term
        cdef double complex two_shear_r_inv, y1_y3_term

        r_inverse = 1. / radius
        two_shear_r_inv = 2. * self.shear_modulus * r_inverse
        density_gravity = self.density * self.gravity
        grav_term = self.grav_coeff * self.density
        y1_y3_term = 2. * y1 - self.llp1 * y3

        cdef double complex dy1, dy2, dy3, dy4, dy5, dy6

        dy1 = y1_y3_term * -1. * r_inverse

        dy2 = r_inverse * (
                y1 * (12. * self.shear_modulus * r_inverse - 4. * density_gravity) +
                y3 * self.llp1 * (density_gravity - 6. * self.shear_modulus * r_inverse) +
                y4 * self.llp1 +
                y5 * self.density * self.lp1 +
                y6 * -self.density * radius
        )

        dy3 = \
            y1 * -r_inverse + \
            y3 * r_inverse + \
            y4 * (1. / self.shear_modulus)

        dy4 = r_inverse * (
                y1 * (density_gravity - 3. * two_shear_r_inv) +
                y2 * -1. +
                y3 * (two_shear_r_inv * (2. * self.llp1 - 1.)) +
                y4 * -3. +
                y5 * -self.density
        )

        dy5 = \
            y1 * grav_term + \
            y5 * -self.lp1 * r_inverse + \
            y6

        dy6 = r_inverse * (
                y1 * grav_term * self.lm1 +
                y6 * self.lm1 +
                y1_y3_term * grav_term
        )

        # Convert back to floats
        self.dy_ptr[0]  = dy1.real
        self.dy_ptr[1]  = dy1.imag
        self.dy_ptr[2]  = dy2.real
        self.dy_ptr[3]  = dy2.imag
        self.dy_ptr[4]  = dy3.real
        self.dy_ptr[5]  = dy3.imag
        self.dy_ptr[6]  = dy4.real
        self.dy_ptr[7]  = dy4.imag
        self.dy_ptr[8]  = dy5.real
        self.dy_ptr[9]  = dy5.imag
        self.dy_ptr[10] = dy6.real
        self.dy_ptr[11] = dy6.imag


cdef class LiquidDynamicCompressible(RadialSolverBase):
    cdef void diffeq(self) noexcept nogil:
        # Note that `t_now` is the current "radius" not time.
        cdef double radius
        radius = self.t_now

        # Update interpolation
        self.update_interp(update_bulk=True, update_shear=False)

        # For the dynamic version, y4 = 0 always in a liquid layer and y3 is defined by y1, y2, and y5 analytically
        # Pull out y values
        cdef double y1_real, y2_real, y5_real, y6_real
        cdef double y1_imag, y2_imag, y5_imag, y6_imag

        y1_real = self.y_ptr[0]
        y1_imag = self.y_ptr[1]
        y2_real = self.y_ptr[2]
        y2_imag = self.y_ptr[3]
        y5_real = self.y_ptr[4]
        y5_imag = self.y_ptr[5]
        y6_real = self.y_ptr[6]
        y6_imag = self.y_ptr[7]

        # Convert floats to complex
        cdef double complex y1, y2, y5, y6

        y1 = y1_real + 1.0j * y1_imag
        y2 = y2_real + 1.0j * y2_imag
        y5 = y5_real + 1.0j * y5_imag
        y6 = y6_real + 1.0j * y6_imag

        # Optimizations
        cdef double r_inverse, density_gravity, f2, dynamic_term, dynamic_term_no_r, grav_term

        r_inverse         = 1. / radius
        density_gravity   = self.density * self.gravity
        f2                = self.frequency * self.frequency
        dynamic_term_no_r = -f2 * self.density
        dynamic_term      = dynamic_term_no_r * radius
        grav_term         = self.grav_coeff * self.density

        # Check if dynamic term is close to zero. It will always be negative so compare to negative eps
        # if dynamic_term < EPS_100:
        #     # TODO: is faking this okay?
        #     dynamic_term = EPS_100

        # Until bulk dissipation is considered, lame_inverse will always be real-valued for a liquid layer.
        cdef double lame_inverse

        # For the liquid layer it is assumed that the shear modulus is zero so the lame parameter simply
        #    equals the bulk modulus. Until bulk dissipation is considered, it will always be real-valued
        lame_inverse = 1. / self.bulk_modulus

        # y3 derivative is undetermined for a liquid layer, but we can calculate its value which is still used in the
        #   other derivatives.
        y3 = (1. / dynamic_term) * (y2 - density_gravity * y1 + self.density * y5)
        y1_y3_term = 2. * y1 - self.llp1 * y3

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
            y5 * self.density * self.lp1 * r_inverse - \
            y6 * self.density - \
            y1_y3_term * density_gravity * r_inverse

        dy5 = \
            y1 * grav_term - \
            y5 * self.lp1 * r_inverse + \
            y6

        dy6 = r_inverse * (
                self.lm1 * (y1 * grav_term + y6) +
                y1_y3_term * grav_term
        )


        # Convert back to floats
        self.dy_ptr[0] = dy1.real
        self.dy_ptr[1] = dy1.imag
        self.dy_ptr[2] = dy2.real
        self.dy_ptr[3] = dy2.imag
        self.dy_ptr[4] = dy5.real
        self.dy_ptr[5] = dy5.imag
        self.dy_ptr[6] = dy6.real
        self.dy_ptr[7] = dy6.imag


cdef class LiquidDynamicIncompressible(RadialSolverBase):
    cdef void diffeq(self) noexcept nogil:
        # Note that `t_now` is the current "radius" not time.
        cdef double radius
        radius = self.t_now

        # Update interpolation
        self.update_interp(update_bulk=False, update_shear=False)

        # Pull out y values
        cdef double y1_real, y2_real, y5_real, y6_real
        cdef double y1_imag, y2_imag, y5_imag, y6_imag

        y1_real = self.y_ptr[0]
        y1_imag = self.y_ptr[1]
        y2_real = self.y_ptr[2]
        y2_imag = self.y_ptr[3]
        y5_real = self.y_ptr[4]
        y5_imag = self.y_ptr[5]
        y6_real = self.y_ptr[6]
        y6_imag = self.y_ptr[7]

        # Convert floats to complex
        cdef double complex y1, y2, y5, y6

        y1 = y1_real + 1.0j * y1_imag
        y2 = y2_real + 1.0j * y2_imag
        y5 = y5_real + 1.0j * y5_imag
        y6 = y6_real + 1.0j * y6_imag

        # Optimizations
        cdef double r_inverse, density_gravity, dynamic_term, grav_term

        r_inverse = 1. / radius
        density_gravity = self.density * self.gravity
        dynamic_term = -self.frequency * self.frequency * self.density * radius
        grav_term = self.grav_coeff * self.density

        # Check if dynamic term is close to zero. It will always be negative so compare to negative eps
        # if dynamic_term < EPS_100:
        #     # TODO: is faking this okay?
        #     dynamic_term = EPS_100

        # y3 derivative is undetermined for a liquid layer, but we can calculate its value which is still used in the
        #   other derivatives.
        cdef double complex y3, y1_y3_term
        y3 = (1. / dynamic_term) * (y2 + self.density * y5 - density_gravity * y1)
        y1_y3_term = 2. * y1 - self.llp1 * y3

        cdef double complex dy1, dy2, dy5, dy6

        dy1 = y1_y3_term * -r_inverse

        dy2 = r_inverse * (
                y1 * (dynamic_term - 2. * density_gravity) +
                y5 * self.density * self.lp1 +
                y6 * -self.density * radius +
                # TODO: In the solid version there is a [2. * (lame + shear_modulus) * r_inverse] coefficient for y1_y3_term
                #   In TS72 the first term is gone. Shouldn't Lame + mu = Lame = Bulk for liquid layer?
                y1_y3_term * -density_gravity
        )

        dy5 = \
            y1 * grav_term + \
            y5 * -self.lp1 * r_inverse + \
            y6

        dy6 = r_inverse * (
                y1 * grav_term * self.lm1 +
                y6 * self.lm1 +
                y1_y3_term * grav_term
        )

        # Convert back to floats
        self.dy_ptr[0] = dy1.real
        self.dy_ptr[1] = dy1.imag
        self.dy_ptr[2] = dy2.real
        self.dy_ptr[3] = dy2.imag
        self.dy_ptr[4] = dy5.real
        self.dy_ptr[5] = dy5.imag
        self.dy_ptr[6] = dy6.real
        self.dy_ptr[7] = dy6.imag


cdef class LiquidStaticCompressible(RadialSolverBase):
    cdef void diffeq(self) noexcept nogil:
        # Note that `t_now` is the current "radius" not time.
        cdef double radius
        radius = self.t_now

        # Update interpolation
        self.update_interp(update_bulk=False, update_shear=False)

        # Pull out y values
        cdef double y5_real, y5_imag
        cdef double y7_real, y7_imag

        # For the static liquid version, only y5 and y7 are defined.
        y5_real = self.y_ptr[0]
        y5_imag = self.y_ptr[1]
        y7_real = self.y_ptr[2]
        y7_imag = self.y_ptr[3]

        # Convert floats to complex
        cdef double complex y5, y7

        y5 = y5_real + 1.0j * y5_imag
        y7 = y7_real + 1.0j * y7_imag

        # Optimizations
        cdef double r_inverse, grav_term

        r_inverse = 1. / radius
        grav_term = self.grav_coeff * self.density / self.gravity

        # See Eq. 18 in S75
        cdef double complex dy5, dy7

        dy5 = \
            y5 * (grav_term - self.lp1 * r_inverse) + \
            y7

        dy7 = \
            y5 * 2. * self.lm1 * r_inverse * grav_term + \
            y7 * (self.lm1 * r_inverse - grav_term)

        # Convert back to floats
        self.dy_ptr[0] = dy5.real
        self.dy_ptr[1] = dy5.imag
        self.dy_ptr[2] = dy7.real
        self.dy_ptr[3] = dy7.imag


cdef class LiquidStaticIncompressible(RadialSolverBase):
    cdef void diffeq(self) noexcept nogil:
        # Note that `t_now` is the current "radius" not time.
        cdef double radius
        radius = self.t_now

        # Update interpolation
        self.update_interp(update_bulk=False, update_shear=False)

        # Pull out y values
        cdef double y5_real, y5_imag
        cdef double y7_real, y7_imag

        # For the static liquid version, only y5 and y7 are defined.
        y5_real = self.y_ptr[0]
        y5_imag = self.y_ptr[1]
        y7_real = self.y_ptr[2]
        y7_imag = self.y_ptr[3]

        # Convert floats to complex
        cdef double complex y5, y7

        y5 = y5_real + 1.0j * y5_imag
        y7 = y7_real + 1.0j * y7_imag

        # Optimizations
        cdef double r_inverse, grav_term

        r_inverse = 1. / radius
        grav_term = self.grav_coeff * self.density / self.gravity

        # See Eq. 18 in S75
        cdef double complex dy5, dy7

        dy5 = \
            y5 * (grav_term - self.lp1 * r_inverse) + \
            y7

        dy7 = \
            y5 * 2. * self.lm1 * r_inverse * grav_term + \
            y7 * (self.lm1 * r_inverse - grav_term)

        # Convert back to floats
        self.dy_ptr[0] = dy5.real
        self.dy_ptr[1] = dy5.imag
        self.dy_ptr[2] = dy7.real
        self.dy_ptr[3] = dy7.imag


cdef RadialSolverBase build_solver(
        bool_cpp_t is_solid,
        bool_cpp_t is_static,
        bool_cpp_t is_incomp,

        # RadialSolverBase Inputs
        Py_ssize_t num_slices,
        double* radius_array_ptr,
        double* density_array_ptr,
        double* gravity_array_ptr,
        double* bulk_modulus_array_ptr,
        double complex* shear_modulus_array_ptr,
        double frequency,
        unsigned int degree_l,
        double G_to_use,

        # Regular CySolver Inputs
        (double, double) t_span,
        const double[::1] y0,
        double* rtols,
        double* atols,
        unsigned char rk_method,
        double max_step,
        Py_ssize_t max_num_steps,
        Py_ssize_t expected_size,

        # Additional optional arguments for RadialSolver class
        bool_cpp_t limit_solution_to_radius
        ):

    cdef RadialSolverBase solver

    if is_solid:
        if is_static:
            if is_incomp:
                solver = SolidStaticIncompressible(
                    frequency,
                    degree_l,
                    G_to_use,
                    t_span,
                    y0,
                    rk_method=rk_method,
                    max_step=max_step,
                    max_num_steps=max_num_steps,
                    expected_size=expected_size,
                    call_first_reset=False,
                    auto_solve=False
                    )
            else:
                solver = SolidStaticCompressible(
                    frequency,
                    degree_l,
                    G_to_use,
                    t_span,
                    y0,
                    rk_method=rk_method,
                    max_step=max_step,
                    max_num_steps=max_num_steps,
                    expected_size=expected_size,
                    call_first_reset=False,
                    auto_solve=False
                    )
        else:
            if is_incomp:
                solver = SolidDynamicIncompressible(
                    frequency,
                    degree_l,
                    G_to_use,
                    t_span,
                    y0,
                    rk_method=rk_method,
                    max_step=max_step,
                    max_num_steps=max_num_steps,
                    expected_size=expected_size,
                    call_first_reset=False,
                    auto_solve=False
                    )
            else:
                solver = SolidDynamicCompressible(
                    frequency,
                    degree_l,
                    G_to_use,
                    t_span,
                    y0,
                    rk_method=rk_method,
                    max_step=max_step,
                    max_num_steps=max_num_steps,
                    expected_size=expected_size,
                    call_first_reset=False,
                    auto_solve=False
                    )
    else:
        if is_static:
            if is_incomp:
                solver = LiquidStaticIncompressible(
                    frequency,
                    degree_l,
                    G_to_use,
                    t_span,
                    y0,
                    rk_method=rk_method,
                    max_step=max_step,
                    max_num_steps=max_num_steps,
                    expected_size=expected_size,
                    call_first_reset=False,
                    auto_solve=False
                    )
            else:
                solver = LiquidStaticCompressible(
                    frequency,
                    degree_l,
                    G_to_use,
                    t_span,
                    y0,
                    rk_method=rk_method,
                    max_step=max_step,
                    max_num_steps=max_num_steps,
                    expected_size=expected_size,
                    call_first_reset=False,
                    auto_solve=False
                    )
        else:
            if is_incomp:
                solver = LiquidDynamicIncompressible(
                    frequency,
                    degree_l,
                    G_to_use,
                    t_span,
                    y0,
                    rk_method=rk_method,
                    max_step=max_step,
                    max_num_steps=max_num_steps,
                    expected_size=expected_size,
                    call_first_reset=False,
                    auto_solve=False
                    )
            else:
                solver = LiquidDynamicCompressible(
                    frequency,
                    degree_l,
                    G_to_use,
                    t_span,
                    y0,
                    rk_method=rk_method,
                    max_step=max_step,
                    max_num_steps=max_num_steps,
                    expected_size=expected_size,
                    call_first_reset=False,
                    auto_solve=False
                    )

    # Install non-python objects
    solver.install_pointers(
        num_slices,
        radius_array_ptr,
        density_array_ptr,
        gravity_array_ptr,
        bulk_modulus_array_ptr,
        shear_modulus_array_ptr,
        atols,
        rtols,
        limit_solution_to_radius=limit_solution_to_radius,
        call_first_reset=False,
        auto_solve=False
        )

    return solver