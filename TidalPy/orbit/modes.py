import numpy as np

from ..performance import njit


@njit
def spin_sync_modes(orbital_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray):
    """

    :param orbital_freq:
    :param eccentricity:
    :param inclination:
    :return:
    """

    modes = (orbital_freq,)
    heating_coeffs = (7. * eccentricity**2 + inclination**2,)
    ztorque_coeffs = (12. * eccentricity**2,)

    return modes, heating_coeffs, ztorque_coeffs

@njit
def nsr_modes(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray):
    """

    These should all have the same shape!!

    :param orbital_freq:
    :param spin_freq:
    :param eccentricity:
    :param inclination:
    :return:
    """

    # TODO: add an explicit check to ensure all of these have the same shape? does njit like type checks like this?

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
         0.75 * e2 * modes[0],
         -0.5 * i2 * modes[1],
         0.5 * i2 * modes[2],
         0.125 * e2 * modes[3],
         6.125 * e2 * modes[4],
         (0.5 - 0.5 * i2 - 2.5 * e2) * modes[5]
    )
    ztorque_coeffs = (
        np.zeros_like(orbital_freq),
        -0.5 * i2,
        0.5 * i2,
        0.25 * e2,
        12.25 * e2,
        (1.0 - i2 - 5.0 * e2)
    )

    return modes, heating_coeffs, ztorque_coeffs