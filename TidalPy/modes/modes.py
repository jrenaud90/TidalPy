from numba import njit
import numpy as np

def spin_sync_modes(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray):
    """

    :param orbital_freq:
    :param spin_freq:
    :param eccentricity:
    :param inclination:
    :return:
    """

    modes = tuple(orbital_freq)
    heating_coeffs = tuple(7. * eccentricity**2 + inclination**2)
    torque_coeffs = tuple(12. * eccentricity**2)

    return modes, heating_coeffs, torque_coeffs


def nsr_modes(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray):
    """

    :param orbital_freq:
    :param spin_freq:
    :param eccentricity:
    :param inclination:
    :return:
    """
    e2 = eccentricity**2
    i2 = inclination**2

    modes = (
            orbital_freq,
            -spin_freq,
            2. * orbital_freq - spin_freq,
            orbital_freq - 2. * spin_freq,
            3. * orbital_freq - 2. * spin_freq,
            2. * orbital_freq - 2. * spin_freq
    )
    heating_coeffs = (
         3. / 4. * e2 * modes[0],
         1. / 2. * i2 * modes[1],
         1. / 2. * i2 * modes[2],
         1. / 8. * e2 * modes[3],
        49. / 8. * e2 * modes[4],
        (0.5 - 0.5 * i2 - (5./2) * e2) * modes[5]


    )



            tuple(7. * eccentricity**2 + inclination**2)
    torque_coeffs = tuple(12. * eccentricity**2)

    return modes, heating_coeffs, torque_coeffs