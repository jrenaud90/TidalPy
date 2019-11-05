import numpy as np

from TidalPy.performance import njit


# Minimum tidal mode in [rad s-1]
MODE_ZERO_TOL = 1.e-12

@njit
def spin_sync_modes(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray,
                      inclination: np.ndarray):
    """ Syncrhonous Rotation Tidal Modes (for heating and torque) for l=2, e^2, I^2

    These should all have the same shape!!

    :param orbital_freq:
    :param spin_freq:
    :param eccentricity:
    :param inclination:
    :return:
    """
    e2 = eccentricity ** 2
    i2 = inclination ** 2

    modes_coeffs = (
        (
            orbital_freq,
            7. * e2 + i2,
            12. * e2
        ),
    )

    heating_coeffs = list()
    ztorque_coeffs = list()
    freqs = list()
    modes = list()

    m_i = 0
    for mode, coeff, torque_coeff in modes_coeffs:
        sgn = np.sign(mode)
        freq = np.abs(mode)

        heating_coeff = freq * coeff
        ztorque_coeff = torque_coeff * sgn

        freq[freq < MODE_ZERO_TOL] = 0.
        mode[freq < MODE_ZERO_TOL] = 0.
        heating_coeff[freq < MODE_ZERO_TOL] = 0.
        ztorque_coeff[freq < MODE_ZERO_TOL] = 0.

        freqs.append(freq)
        modes.append(mode)
        heating_coeffs.append(heating_coeff)
        ztorque_coeffs.append(ztorque_coeff)

        m_i += 1

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return modes, freqs, heating_coeffs, ztorque_coeffs


@njit
def spin_sync_modes_4(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray):
    """ Syncrhonous Rotation Tidal Modes (for heating and torque) for l=2, e^4, I^4

    These should all have the same shape!!

    :param orbital_freq:
    :param spin_freq:
    :param eccentricity:
    :param inclination:
    :return:
    """
    e2 = eccentricity**2
    e4 = eccentricity**4
    i2 = inclination**2
    i4 = inclination**4
    i2e2 = e2 * i2

    modes_coeffs = (
        (
            orbital_freq,
            7. * e2 - (101./4.) * e4 + i2 - (13./12.) * i4 - (19./2.) * i2e2,
            12. * e2 - (215./4.) * e4 + (1./4.) * i4 - 16. * i2e2
        ),
        (
            2. * orbital_freq,
            (605./16.) * e4 + (5./16.) * i4 + (29./4.) * i2e2,
            (289./4.) * e4 - (1./4.) * i4 + 5. * i2e2
        ),
    )

    heating_coeffs = list()
    ztorque_coeffs = list()
    freqs = list()
    modes = list()

    m_i = 0
    for mode, coeff, torque_coeff in modes_coeffs:
        sgn = np.sign(mode)
        freq = np.abs(mode)

        heating_coeff = freq * coeff
        ztorque_coeff = torque_coeff * sgn

        freq[freq < MODE_ZERO_TOL] = 0.
        mode[freq < MODE_ZERO_TOL] = 0.
        heating_coeff[freq < MODE_ZERO_TOL] = 0.
        ztorque_coeff[freq < MODE_ZERO_TOL] = 0.

        freqs.append(freq)
        modes.append(mode)
        heating_coeffs.append(heating_coeff)
        ztorque_coeffs.append(ztorque_coeff)

        m_i += 1

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return modes, freqs, heating_coeffs, ztorque_coeffs


@njit
def spin_sync_modes_6(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray):
    """ Syncrhonous Rotation Tidal Modes (for heating and torque) for l=2, e^6, I^6

    These should all have the same shape!!

    :param orbital_freq:
    :param spin_freq:
    :param eccentricity:
    :param inclination:
    :return:
    """
    e2 = eccentricity**2
    e4 = eccentricity**4
    e6 = eccentricity**6
    i2 = inclination**2
    i4 = inclination**4
    i6 = inclination**6
    i2e2 = e2 * i2
    i2e4 = i2 * e4
    i4e2 = i4 * e2

    modes_coeffs = (
        (
            orbital_freq,
            7. * e2 - (101./4.) * e4 + (551./12.) * e6 + i2 - (13./12.) * i4 + (739./1440.) * i6 - (19./2.) * i2e2
            + (1003./32.) * i2e4 + (1097./192.) * i4e2,
            12. * e2 - (215./4.) * e4 + (8239./96.) * e6 + (1./4.) * i4 - (19./96.) * i6 - 16. * i2e2
            + (1831./32.) * i2e4 + (433./48.) * i4e2
        ),
        (
            2. * orbital_freq,
            (605./16.) * e4 - (3847./24.) * e6 + (5./16.) * i4 - (5./24.) * i6 + (29./4.) * i2e2 - (1049./16.) * i2e4
            - (43./6.) * i4e2,
            (289./4.) * e4 - (1955./6.) * e6 - (1./4.) * i4 + (1./6.) * i6 + 5. * i2e2 - (1627./16.) * i2e4
            - (209./48.) * i4e2
        ),
        (
            3. * orbital_freq,
            (2855./18.) * e6 + (1./32.) * i6 + (1237./32.) * i2e4 + (165./64.) * i4e2,
            (9917./32.) * e6 - (1./32.) * i6 + (1075./32.) * i2e4 - (9./16.) * i4e2
        )
    )

    heating_coeffs = list()
    ztorque_coeffs = list()
    freqs = list()
    modes = list()

    m_i = 0
    for mode, coeff, torque_coeff in modes_coeffs:
        sgn = np.sign(mode)
        freq = np.abs(mode)

        heating_coeff = freq * coeff
        ztorque_coeff = torque_coeff * sgn

        freq[freq < MODE_ZERO_TOL] = 0.
        mode[freq < MODE_ZERO_TOL] = 0.
        heating_coeff[freq < MODE_ZERO_TOL] = 0.
        ztorque_coeff[freq < MODE_ZERO_TOL] = 0.

        freqs.append(freq)
        modes.append(mode)
        heating_coeffs.append(heating_coeff)
        ztorque_coeffs.append(ztorque_coeff)

        m_i += 1

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return modes, freqs, heating_coeffs, ztorque_coeffs

@njit
def nsr_modes(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray):
    """ Non-Syncrhonous Rotation Tidal Modes (for heating and torque) for l=2, e^2, I^2

    These should all have the same shape!!

    :param orbital_freq:
    :param spin_freq:
    :param eccentricity:
    :param inclination:
    :return:
    """
    e2 = eccentricity ** 2
    i2 = inclination ** 2

    modes_coeffs = (
        (
            orbital_freq,
            (3. / 4.) * e2,
            0.
        ),
        (
            spin_freq,
            (1. / 2.) * i2,
            -1.
        ),
        (
            2. * orbital_freq - spin_freq,
            (1. / 2.) * i2,
            1.
        ),
        (
            orbital_freq - 2. * spin_freq,
            (1. / 8.) * e2,
            2.
        ),
        (
            2. * orbital_freq - 2. * spin_freq,
            (1. / 2.) - (5. / 2.) * e2 - (1. / 2.) * i2,
            2.
        ),
        (
            3. * orbital_freq - 2. * spin_freq,
            (49. / 8.) * e2,
            2.
        ),
    )

    heating_coeffs = list()
    ztorque_coeffs = list()
    freqs = list()
    modes = list()

    m_i = 0
    for mode, coeff, torque_multi in modes_coeffs:
        sgn = np.sign(mode)
        freq = np.abs(mode)

        heating_coeff = freq * coeff
        ztorque_coeff = coeff * torque_multi * sgn

        freq[freq < MODE_ZERO_TOL] = 0.
        mode[freq < MODE_ZERO_TOL] = 0.
        heating_coeff[freq < MODE_ZERO_TOL] = 0.
        ztorque_coeff[freq < MODE_ZERO_TOL] = 0.

        freqs.append(freq)
        modes.append(mode)
        heating_coeffs.append(heating_coeff)
        ztorque_coeffs.append(ztorque_coeff)

        m_i += 1

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return modes, freqs, heating_coeffs, ztorque_coeffs


@njit
def nsr_modes_4(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray):
    """ Non-Syncrhonous Rotation Tidal Modes (for heating and torque) for l=2, e^4, I^4

    These should all have the same shape!!

    :param orbital_freq:
    :param spin_freq:
    :param eccentricity:
    :param inclination:
    :return:
    """
    e2 = eccentricity ** 2
    e4 = eccentricity ** 4
    i2 = inclination ** 2
    i4 = inclination ** 4
    i2e2 = e2 * i2

    modes_coeffs = (
        (
            orbital_freq,
            (3. / 4.) * e2 + (27. / 16.) * e4 - (9. / 4.) * i2e2,
            0.
        ),
        (
            2. * orbital_freq,
            (27. / 16.) * e4 + (3. / 16.) * i4,
            0.
        ),
        (
            spin_freq,
            (1. / 2.) * i2 - (2. / 3.) * i4 + (3. / 2.) * i2e2,
            -1.
        ),
        (
            2. * spin_freq,
            (1. / 8.) * i4,
            -2.
        ),
        (
            orbital_freq - spin_freq,
            (5. / 4.) * i2e2,
            1.
        ),
        (
            2. * orbital_freq - spin_freq,
            (1. / 2.) * i2 - (5. / 12.) * i4 - (5. / 2.) * i2e2,
            1.
        ),
        (
            3. * orbital_freq - spin_freq,
            (49. / 8.) * i2e2,
            1.
        ),
        (
            orbital_freq - 2. * spin_freq,
            (1. / 8.) * e2 - (1. / 32.) * e4 - (1. / 8.) * i2e2,
            2.
        ),
        (
            2. * orbital_freq - 2. * spin_freq,
            (1. / 2.) - (5. / 2.) * e2 + (63. / 16.) * e4 - (1. / 2.) * i2 + (11. / 48.) * i4 + (5. / 2.) * i2e2,
            2.
        ),
        (
            3. * orbital_freq - 2. * spin_freq,
            (49. / 8.) * e2 - (861. / 32.) * e4 - (49. / 8.) * i2e2,
            2.
        ),
        (
            4. * orbital_freq - 2. * spin_freq,
            (289. / 8.) * e4,
            2.
        ),
        (
            -1. * orbital_freq - spin_freq,
            (9. / 8.) * i2e2,
            1.
        )
    )

    heating_coeffs = list()
    ztorque_coeffs = list()
    freqs = list()
    modes = list()

    m_i = 0
    for mode, coeff, torque_multi in modes_coeffs:
        sgn = np.sign(mode)
        freq = np.abs(mode)

        heating_coeff = freq * coeff
        ztorque_coeff = coeff * torque_multi * sgn

        freq[freq < MODE_ZERO_TOL] = 0.
        mode[freq < MODE_ZERO_TOL] = 0.
        heating_coeff[freq < MODE_ZERO_TOL] = 0.
        ztorque_coeff[freq < MODE_ZERO_TOL] = 0.

        freqs.append(freq)
        modes.append(mode)
        heating_coeffs.append(heating_coeff)
        ztorque_coeffs.append(ztorque_coeff)

        m_i += 1

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return modes, freqs, heating_coeffs, ztorque_coeffs


@njit
def nsr_modes_6(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray):
    """ Non-Syncrhonous Rotation Tidal Modes (for heating and torque) for l=2, e^6, I^6

    These should all have the same shape!!

    :param orbital_freq:
    :param spin_freq:
    :param eccentricity:
    :param inclination:
    :return:
    """
    e2 = eccentricity**2
    e4 = eccentricity**4
    e6 = eccentricity**6
    i2 = inclination**2
    i4 = inclination**4
    i6 = inclination**6
    i2e2 = e2 * i2
    i2e4 = i2 * e4
    i4e2 = i4 * e2

    modes_coeffs = (
        (
            orbital_freq,
            (3./4.) * e2 + (27./16.) * e4 + (765./256.) * e6 - (9./4.) * i2e2 - (81./16.) * i2e4 + (159./64.) * i4e2,
            0.
        ),
        (
            2. * orbital_freq,
            (27./16.) * e4 + (21./8.) * e6 - (81./16.) * i2e4 + (3./16.) * i4 - (15./16.) * i4e2 - (1./8.) * i6,
            0.
        ),
        (
            3. * orbital_freq,
            (2809./768.) * e6 + (147./64.) * i4e2,
            0.
        ),
        (
            spin_freq,
            (1./2.) * i2 - (2./3.) * i4 + (16./45.) * i6 + (3./2.) * i2e2 + 3. * i2e4 - 2. * i4e2,
            -1.
        ),
        (
            2. * spin_freq,
            (1./8.) * i4 - (1./12.) * i6 + (3./8.) * i4e2,
            -2.
         ),
        (
            orbital_freq - spin_freq,
            (5./4.) * i2e2 + (5./2.) * i2e4 - (77./48.) * i4e2,
            1.
        ),
        (
            2. * orbital_freq - spin_freq,
            (1./2.) * i2 - (5./12.) * i4 + (227./1440.) * i6 - (5./2.) * i2e2 + (207./32.) * i2e4 + (25./12.) * i4e2,
            1.
        ),
        (
            3. * orbital_freq - spin_freq,
            (49./8.) * i2e2 - (861./32.) * i2e4 - (245./48.) * i4e2,
            1.
        ),
        (
            4.* orbital_freq - spin_freq,
            (289./8.) * i2e4,
            1.
        ),
        (
            orbital_freq - 2. * spin_freq,
            (1./8.) * e2 - (1./32.) * e4 + (13./1536.) * e6 - (1./8.) * i2e2 + (1./32.) * i2e4 + (65./192.) * i4e2,
            2.
        ),
        (
            2. * orbital_freq - 2. * spin_freq,
            (1./2.) - (5./2.) * e2 + (63./16.) * e4 - (79./36.) * e6 - (1./2.) * i2 + (11./48.) * i4 - (23./360.) * i6
            + (5./2.) * i2e2 - (63./16.) * i2e4 - (55./48.) * i4e2,
            2.
        ),
        (
            3. * orbital_freq - 2. * spin_freq,
            (49./8.) * e2 - (861./32.) * e4 + (21975./512.) * e6 - (49./8.) * i2e2 + (861./32.) * i2e4
            + (539./192.) * i4e2,
            2.
        ),
        (
            4. * orbital_freq - 2. * spin_freq,
            (289./8.) * e4 - (1955./12.) * e6 - (289./8.) * i2e4,
            2.
        ),
        (
            5. * orbital_freq - 2. * spin_freq,
            (714025./4608.) * e6,
            2.
        ),
        (
            -1. * orbital_freq - spin_freq,
            (9./8.) * i2e2 + (81./32.) * i2e4 - (3./2.) * i4e2,
            1.
        ),
        (
            -1. * orbital_freq - 2. * spin_freq,
            (1./4608.) * e6 + (9./32.) * i4e2,
            2.
        ),
        (
            -2. * orbital_freq - 1. * spin_freq,
            (1./32.) * i6 + (81./32.) * i2e4,
            1.
        )
    )

    heating_coeffs = list()
    ztorque_coeffs = list()
    freqs = list()
    modes = list()

    m_i = 0
    for mode, coeff, torque_multi in modes_coeffs:
        sgn = np.sign(mode)
        freq = np.abs(mode)

        heating_coeff = freq * coeff
        ztorque_coeff = coeff * torque_multi * sgn

        freq[freq < MODE_ZERO_TOL] = 0.
        mode[freq < MODE_ZERO_TOL] = 0.
        heating_coeff[freq < MODE_ZERO_TOL] = 0.
        ztorque_coeff[freq < MODE_ZERO_TOL] = 0.

        freqs.append(freq)
        modes.append(mode)
        heating_coeffs.append(heating_coeff)
        ztorque_coeffs.append(ztorque_coeff)

        m_i += 1

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return modes, freqs, heating_coeffs, ztorque_coeffs