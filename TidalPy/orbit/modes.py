import numpy as np

from ..performance import njit


# Minimum tidal mode in [rad s-1]
MODE_ZERO_TOL = 1.e-10


@njit
def spin_sync_modes(orbital_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray):
    """

    :param orbital_freq:
    :param eccentricity:
    :param inclination:
    :return:
    """

    modes = (orbital_freq,)
    freqs = (np.abs(orbital_freq),)
    heating_coeffs = ((7. * eccentricity**2 + inclination**2) * freqs[0],)
    ztorque_coeffs = (12. * eccentricity**2 * np.ones_like(freqs[0]),)

    return modes, freqs, heating_coeffs, ztorque_coeffs


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

    unchecked_modes = [
        orbital_freq,
        spin_freq,
        2. * orbital_freq - spin_freq,
        orbital_freq - 2. * spin_freq,
        3. * orbital_freq - 2. * spin_freq,
        2. * orbital_freq - 2. * spin_freq
    ]

    modes = list()
    for mode in unchecked_modes:
        mode[np.abs(mode) < MODE_ZERO_TOL] = 0.
        modes.append(mode)

    freqs = [np.abs(mode) for mode in modes]
    signs = [np.sign(mode) for mode in modes]

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    heating_coeffs = (
        # The sign will cancel with the Im(k2) function (which is an odd function) so no need to add them here.
        (3. / 4.) * e2 * freqs[0],
        (1. / 2.) * i2 * freqs[1],
        (1. / 2.) * i2 * freqs[2],
        (1. / 8.) * e2 * freqs[3],
        (49. / 8.) * e2 * freqs[4],
        (0.5 - 0.5 * i2 - 2.5 * e2) * freqs[5]
    )

    ztorque_coeffs = (
        # Don't need signs on these first two as the signs will always be positive
        np.zeros_like(orbital_freq),
        -0.5 * i2,
        0.5 * i2 * signs[2],
        0.25 * e2 * signs[3],
        12.25 * e2 * signs[4],
        (1.0 - i2 - 5.0 * e2) * signs[5]
    )

    return modes, freqs, heating_coeffs, ztorque_coeffs
