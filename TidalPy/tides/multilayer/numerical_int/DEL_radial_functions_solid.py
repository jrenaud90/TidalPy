from ....constants import pi, G
import numpy as np


def radial_func_derivatives_general(radius, radial_functions, shear_modulus, bulk_modulus, density, gravity, frequency,
                                    order_l: int = 2):
    """"""

    y1, y2, y3, y4, y5, y6 = radial_functions

    # Convert compressibility parameters
    poisson = (3. * bulk_modulus - 2. * shear_modulus) / (6. * bulk_modulus + 2. * shear_modulus)
    poisson_frac = (1. + poisson) / (1. - poisson)
    compressibility = (1. - 2. * poisson) / (1. - poisson)
    laplace_eigenvalue = (order_l - 1.) * (order_l + 2.)

    # Optimizations
    dynamic_term = frequency * frequency
    r_inverse = 1. / radius
    r2_inverse = r_inverse * r_inverse

    # See Eqs. 13 - 18 in B15
    # # dy2 and dy4 contain all three of: dynamic, viscoelastic, and gravitational terms.
    dy1 = compressibility / (2. * shear_modulus) * y2 - r_inverse * (1. - compressibility) * \
          (2. * y1 - order_l * (order_l + 1.) * y3)

    dy2 = -dynamic_term * density * y1 - (2. * r_inverse) * compressibility * y2 + \
          r2_inverse * (2. * shear_modulus * poisson_frac - density * gravity * radius) * \
              (2. * y1 - order_l * (order_l + 1.) * y3) + \
          order_l * (order_l + 1.) * r_inverse * y4 - \
          density * (y6 - (order_l + 1.) * r_inverse * y5 + 2. * r_inverse * gravity * y1)

    dy3 = (1. / shear_modulus) * y4 + r_inverse * (y3 - y1)

    dy4 = -dynamic_term * density * y3 - r_inverse * (1. - compressibility) * y2 - \
          (2. * r2_inverse) * shear_modulus * \
              (poisson_frac * y1 - ((laplace_eigenvalue + 1. + poisson) / (1. - poisson)) * y3) - \
          3. * r_inverse * y4 - density * r_inverse * (y5 - gravity * y1)

    dy5 = y6 + 4. * pi * G * density * y1 - (order_l + 1.) * r_inverse * y5

    dy6 = (order_l - 1.) * r_inverse * y6 + \
          4. * pi * G * density * r_inverse * (order_l + 1.) * (y1 - order_l * y3)

    return (dy1, dy2, dy3, dy4, dy5, dy6)


def radial_func_derivatives_general_liq(radius, radial_functions, shear_modulus, bulk_modulus, density, gravity, frequency,
                                    order_l: int = 2):
    """"""

    y1, y2, y5, y6 = radial_functions

    # Convert compressibility parameters
    lame_1 = bulk_modulus - (2. / 3.) * shear_modulus

    # Optimizations
    dynamic_term = frequency * frequency
    r_inverse = 1. / radius
    r2_inverse = r_inverse * r_inverse
    llp1 = order_l * (order_l + 1.)

    # See Eqs. 13 - 18 in B15
    # # dy2 and dy4 contain all three of: dynamic, viscoelastic, and gravitational terms.
    dy1 = (-2. * r_inverse + r2_inverse * llp1 * gravity / dynamic_term) * y1 + \
          ((1. / lame_1) - r2_inverse * llp1 / (dynamic_term * density)) * y2 - \
          (llp1 * r2_inverse / dynamic_term) * y5

    dy2 = (-dynamic_term * density - 4. * density * gravity * r_inverse + llp1 * density * gravity**2 * r2_inverse / dynamic_term) * y1 - \
          (llp1 * gravity * r2_inverse / dynamic_term) * y2 + \
          (order_l + 1) * density * r_inverse * (1 - order_l * gravity * r_inverse / dynamic_term) * y5 - \
          density * y6

    dy5 = 4. * pi * G * density * y1 - (order_l + 1) * r_inverse * y5 + y6

    dy6 = 4. * pi * (order_l + 1.) * G * density * r_inverse * (1 - order_l * gravity * r_inverse / dynamic_term) * y1 + \
          (4. * pi * llp1 * G * r2_inverse / dynamic_term) * y2 + \
          (4. * pi * llp1 * density * G * r2_inverse / dynamic_term) * y5 + \
          (order_l - 1.) * r_inverse * y6

    return (dy1, dy2, dy5, dy6)

# def radial_func_derivatives_general_liq_y7(radius, radial_functions, shear_modulus, bulk_modulus, density, gravity, frequency,
#                                     order_l: int = 2):
#     """"""
#
#     y1, y2, y5, y7 = radial_functions
#
#     # Convert compressibility parameters
#     poisson = (3. * bulk_modulus - 2. * shear_modulus) / (6. * bulk_modulus + 2. * shear_modulus)
#     poisson_frac = (1. + poisson) / (1. - poisson)
#     compressibility = (1. - 2. * poisson) / (1. - poisson)
#     lame_1 = bulk_modulus - (2. / 3.) * shear_modulus
#     laplace_eigenvalue = (order_l - 1.) * (order_l + 2.)
#
#     llp1 = order_l * (order_l + 1.)
#
#     # Optimizations
#     dynamic_term = frequency * frequency
#     r_inverse = 1. / radius
#     r2_inverse = r_inverse * r_inverse
#
#     # See Eqs. 13 - 18 in B15
#     # # dy2 and dy4 contain all three of: dynamic, viscoelastic, and gravitational terms.
#     dy1 = (-2. * r_inverse + r2_inverse * llp1 * gravity / dynamic_term) * y1 + \
#           ((1. / lame_1) - r2_inverse * llp1 / (dynamic_term * density)) * y2 - \
#           (llp1 * r2_inverse / dynamic_term) * y5
#
#     dy2 = (-dynamic_term * density - 4. * density * gravity * r_inverse + llp1 * density * gravity**2 * r2_inverse / dynamic_term) * y1 - \
#           (llp1 * gravity * r2_inverse / dynamic_term) * y2 + \
#           (order_l + 1) * density * r_inverse * (1 - order_l * gravity * r_inverse / dynamic_term) * y5 - \
#           density * y6
#
#     dy5 = 4. * pi * G * density * y1 - (order_l + 1) * r_inverse * y5 + y6
#
#     dy6 = 4. * pi * (order_l + 1.) * G * density * r_inverse * (1 - order_l * gravity * r_inverse / dynamic_term) * y1 + \
#           (4. * pi * llp1 * G * r2_inverse / dynamic_term) * y2 + \
#           (4. * pi * llp1 * density * G * r2_inverse / dynamic_term) * y5 + \
#           (order_l - 1.) * r_inverse * y6
#
#     return (dy1, dy2, dy5, dy6)