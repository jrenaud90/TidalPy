import numpy as np

from ....constants import pi, G
from ....utilities.performance import njit

@njit(cacheable=True)
def z_calc(x_squared, order_l, init_l=20):

    high_z = x_squared / (2 * init_l + 3)
    for l_fake in range(order_l-1, init_l-1):
        l = init_l - l_fake
        # OPT: The above range is nicely written as range(order_l, init_l)[::-1]; but njit does not support this atm.
        high_z = x_squared / ((2 * l + 1) - high_z)

    return high_z

def solid_solutions(radius, bulk_modulus, shear_modulus, density, gravity, frequency, order_l: int = 2):

    # Convert compressibility parameters
    poisson = (3. * bulk_modulus - 2. * shear_modulus) / (6. * bulk_modulus + 2. * shear_modulus)
    poisson_frac = (1. + poisson) / (1. - poisson)
    compressibility = (1. - 2. * poisson) / (1. - poisson)
    lame_1 = bulk_modulus - (2. / 3.) * shear_modulus
    laplace_eigenvalue = (order_l - 1.) * (order_l + 2.)

    # Constants
    w2 = frequency * frequency
    alpha_2 = lame_1 + 2. * shear_modulus / density
    beta_2 = shear_modulus / density
    gamma = 4. * pi * G * density / 3.

    k2_quad_pos = (w2 / beta_2) + (w2 + 4. * gamma) / alpha_2
    k2_quad_neg = (w2 / beta_2) - (w2 + 4. * gamma) / alpha_2
    k2_quad = k2_quad_neg**2 + ((4. * order_l * (order_l + 1) * gamma**2) / (alpha_2 * beta_2))

    k2_pos = (1. / 2.) * (k2_quad_pos + np.sqrt(k2_quad))
    k2_neg = (1. / 2.) * (k2_quad_pos - np.sqrt(k2_quad))

    f_k2_pos = (beta_2 * k2_pos - w2) / gamma
    f_k2_neg = (beta_2 * k2_neg - w2) / gamma

    h_k2_pos = f_k2_pos - (order_l + 1)
    h_k2_neg = f_k2_neg - (order_l + 1)

    z_k2_pos = z_calc(k2_pos * radius**2, order_l=order_l)
    z_k2_neg = z_calc(k2_neg * radius**2, order_l=order_l)

    # Appendix B, Kamata+2015
    y1_s1 = -f_k2_pos * z_k2_pos / radius
    y1_s2 = -f_k2_neg * z_k2_neg / radius
    y2_s1 = -density * f_k2_pos * alpha_2 * k2_pos + (2. * shear_modulus / (radius**2)) * \
                                                     (2 * f_k2_pos + order_l * (order_l + 1)) * z_k2_pos
    y2_s2 = -density * f_k2_neg * alpha_2 * k2_neg + (2. * shear_modulus / (radius**2)) * \
                                                     (2 * f_k2_neg + order_l * (order_l + 1)) * z_k2_neg
    y3_s1 = (1. / radius) * z_k2_pos
    y3_s2 = (1. / radius) * z_k2_neg
    y4_s1 = shear_modulus * k2_pos - (2. * shear_modulus / (radius**2)) * (f_k2_pos + 1.) * z_k2_pos
    y4_s2 = shear_modulus * k2_neg - (2. * shear_modulus / (radius**2)) * (f_k2_neg + 1.) * z_k2_neg
    y5_s1 = 3. * gamma * f_k2_pos - h_k2_pos * (order_l * gamma - w2)
    y5_s2 = 3. * gamma * f_k2_neg - h_k2_neg * (order_l * gamma - w2)
    y6_s1 = (2. * order_l + 1) * y5_s1 / radius
    y6_s2 = (2. * order_l + 1) * y5_s2 / radius

    y1_s3 = order_l / radius
    y2_s3 = 2. * shear_modulus * order_l * (order_l - 1) / (radius**2)
    y3_s3 = (1. / radius)
    y4_s3 = 2. * shear_modulus * (order_l - 1.) / (radius**2)
    y5_s3 = order_l * gamma - w2
    y6_s3 = ((2 * order_l + 1) * y5_s3 / radius) - (3. * order_l * gamma / radius)

    tidaly_s1 = np.vstack((y1_s1, y2_s1, y3_s1, y4_s1, y5_s1, y6_s1))
    tidaly_s2 = np.vstack((y1_s2, y2_s2, y3_s2, y4_s2, y5_s2, y6_s2))
    tidaly_s3 = np.vstack((y1_s3, y2_s3, y3_s3, y4_s3, y5_s3, y6_s3))

    return (tidaly_s1, tidaly_s2, tidaly_s3)

def liquid_solutions(radius, bulk_modulus, shear_modulus, density, gravity, frequency, order_l: int = 2):

    # Convert compressibility parameters
    poisson = (3. * bulk_modulus - 2. * shear_modulus) / (6. * bulk_modulus + 2. * shear_modulus)
    poisson_frac = (1. + poisson) / (1. - poisson)
    compressibility = (1. - 2. * poisson) / (1. - poisson)
    lame_1 = bulk_modulus - (2. / 3.) * shear_modulus
    laplace_eigenvalue = (order_l - 1.) * (order_l + 2.)

    # Constants
    gamma = (4. * pi * G * density / 3.)
    w2 = frequency**2
    alpha_2 = lame_1 / density
    k2 = (1 / alpha_2) * (w2 + 4. * gamma - (order_l * (order_l + 1) * gamma**2 / w2))
    f = -w2 / gamma
    h = f - (order_l + 1)

    # Appendix B, Kamata+2015
    y1_s1 = -(f / radius) * z_calc(k2 * radius**2, order_l=order_l)
    y1_s2 = order_l / radius

    y2_s1 = -density * (f * (w2 + 4 * gamma) + order_l * (order_l + 1) * gamma)
    y2_s2 = 0. * radius

    y5_s1 = 3. * gamma * f - h * (order_l * gamma - w2)
    y5_s2 = order_l * gamma - w2

    y6_s1 = (2 * order_l + 1) * y5_s1 / radius
    y6_s2 = ((2 * order_l + 1) * y5_s2 / radius) - ((3. * order_l * gamma) / radius)

    tidaly_s1 = np.vstack((y1_s1, y2_s1, y5_s1, y6_s1))
    tidaly_s2 = np.vstack((y1_s2, y2_s2, y5_s2, y6_s2))

    return (tidaly_s1, tidaly_s2)