import numpy as np

from TidalPy.performance import njit
from . import MODE_ZERO_TOL

# TODO: A lot of the mode list construction is static and could take place outside of the functions at the expense of
#   readibility.

@njit
def spin_sync_modes(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray,
                    inclination: np.ndarray):
    """ Synchronous Rotation Tidal Modes (for heating and torque) for l=2, e^2, I^2

    These should all have the same shape!!
    """
    e2 = eccentricity ** 2
    i2 = inclination ** 2

    modes_coeffs = (
        (
            'n',
            orbital_freq,
            7. * e2 + i2,
            12. * e2
        ),
    )

    heating_coeffs = list()
    ztorque_coeffs = list()
    freqs = list()
    modes = list()
    mode_names = list()

    for mode_name, mode, coeff, torque_coeff in modes_coeffs:
        sgn = np.sign(mode)
        freq = np.abs(mode)

        heating_coeff = freq * coeff
        ztorque_coeff = torque_coeff * sgn

        freq_too_low = freq < MODE_ZERO_TOL
        freq[freq_too_low] = 0.
        mode[freq_too_low] = 0.
        heating_coeff[freq_too_low] = 0.
        ztorque_coeff[freq_too_low] = 0.

        freqs.append(freq)
        modes.append(mode)
        mode_names.append(mode_name)
        heating_coeffs.append(heating_coeff)
        ztorque_coeffs.append(ztorque_coeff)

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return mode_names, modes, freqs, heating_coeffs, ztorque_coeffs


@njit
def spin_sync_modes_4(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray):
    """ Synchronous Rotation Tidal Modes (for heating and torque) for l=2, e^4, I^4

    These should all have the same shape!!
    """
    e2 = eccentricity**2
    e4 = eccentricity**4
    i2 = inclination**2
    i4 = inclination**4
    i2e2 = e2 * i2

    modes_coeffs = (
        (
            'n',
            orbital_freq,
            7. * e2 - (101./4.) * e4 + i2 - (13./12.) * i4 - (19./2.) * i2e2,
            12. * e2 - (215./4.) * e4 + (1./4.) * i4 - 16. * i2e2
        ),
        (
            '2n',
            2. * orbital_freq,
            (605./16.) * e4 + (5./16.) * i4 + (29./4.) * i2e2,
            (289./4.) * e4 - (1./4.) * i4 + 5. * i2e2
        ),
    )

    heating_coeffs = list()
    ztorque_coeffs = list()
    freqs = list()
    modes = list()
    mode_names = list()

    for mode_name, mode, coeff, torque_coeff in modes_coeffs:
        sgn = np.sign(mode)
        freq = np.abs(mode)

        heating_coeff = freq * coeff
        ztorque_coeff = torque_coeff * sgn

        freq_too_low = freq < MODE_ZERO_TOL
        freq[freq_too_low] = 0.
        mode[freq_too_low] = 0.
        heating_coeff[freq_too_low] = 0.
        ztorque_coeff[freq_too_low] = 0.

        freqs.append(freq)
        modes.append(mode)
        mode_names.append(mode_name)
        heating_coeffs.append(heating_coeff)
        ztorque_coeffs.append(ztorque_coeff)

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return mode_names, modes, freqs, heating_coeffs, ztorque_coeffs


@njit
def spin_sync_modes_6(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray):
    """ Synchronous Rotation Tidal Modes (for heating and torque) for l=2, e^6, I^6

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
            'n',
            orbital_freq,
            7. * e2 - (101./4.) * e4 + (551./12.) * e6 + i2 - (13./12.) * i4 + (739./1440.) * i6 - (19./2.) * i2e2
            + (1003./32.) * i2e4 + (1097./192.) * i4e2,
            12. * e2 - (215./4.) * e4 + (8239./96.) * e6 + (1./4.) * i4 - (19./96.) * i6 - 16. * i2e2
            + (1831./32.) * i2e4 + (433./48.) * i4e2
        ),
        (
            '2n',
            2. * orbital_freq,
            (605./16.) * e4 - (3847./24.) * e6 + (5./16.) * i4 - (5./24.) * i6 + (29./4.) * i2e2 - (1049./16.) * i2e4
            - (43./6.) * i4e2,
            (289./4.) * e4 - (1955./6.) * e6 - (1./4.) * i4 + (1./6.) * i6 + 5. * i2e2 - (1627./16.) * i2e4
            - (209./48.) * i4e2
        ),
        (
            '3n',
            3. * orbital_freq,
            (2855./18.) * e6 + (1./32.) * i6 + (1237./32.) * i2e4 + (165./64.) * i4e2,
            (9917./32.) * e6 - (1./32.) * i6 + (1075./32.) * i2e4 - (9./16.) * i4e2
        )
    )

    heating_coeffs = list()
    ztorque_coeffs = list()
    freqs = list()
    modes = list()
    mode_names = list()

    for mode_name, mode, coeff, torque_coeff in modes_coeffs:
        sgn = np.sign(mode)
        freq = np.abs(mode)

        heating_coeff = freq * coeff
        ztorque_coeff = torque_coeff * sgn

        freq_too_low = freq < MODE_ZERO_TOL
        freq[freq_too_low] = 0.
        mode[freq_too_low] = 0.
        heating_coeff[freq_too_low] = 0.
        ztorque_coeff[freq_too_low] = 0.

        freqs.append(freq)
        modes.append(mode)
        mode_names.append(mode_name)
        heating_coeffs.append(heating_coeff)
        ztorque_coeffs.append(ztorque_coeff)

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return mode_names, modes, freqs, heating_coeffs, ztorque_coeffs

@njit
def spin_sync_modes_8(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray):
    """ Synchronous Rotation Tidal Modes (for heating and torque) for l=2, e^8, I^8

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
    e8 = eccentricity**8
    i2 = inclination**2
    i4 = inclination**4
    i6 = inclination**6
    i8 = inclination**8
    i2e2 = i2 * e2
    i2e4 = i2 * e4
    i2e6 = i2 * e6
    i4e2 = i4 * e2
    i6e2 = i6 * e2
    i4e4 = i4 * e4

    modes_coeffs = (
        (
            'n',
            orbital_freq,
            7. * e2 - (101./4.) * e4 + (551./12.) * e6 - (65773./2304.) * e8 +
            i2 - (13./12.) * i4 + (739./1440.) * i6 - (2773./20160.) * i8 -
            (19./2.) * i2e2 + (1003./32.) * i2e4 - (52013./1152.) * i2e6 +
            (1097./192.) * i4e2 - (12977./768.) * i4e4 - (707./360.) * i6e2,
            12. * e2 - (215. / 4.) * e4 + (8239. / 96.) * e6 - (305561. / 4608.) * e8 +
            (1. / 4.) * i4 - (19. / 96.) * i6 + (21. / 320.) * i8 -
            16 * i2e2 + (1831. / 32.) * i2e4 - (25655. / 288.) * i2e6 +
            (433. / 48.) * i4e2 - (5483. / 192.) * i4e4 - (4339. / 1440.) * i6e2
        ),
        (
            '2n',
            2. * orbital_freq,
            (605./16.) * e4 - (3847./24.) * e6 + (339187./1152.) * e8 +
            (5./16.) * i4 - (5./24.) * i6 + (1./16.) * i8 +
            (29./4.) * i2e2 - (1049./16.) * i2e4 + (14971./72.) * i2e6 -
            (43./6.) * i4e2 + (16633./384.) * i4e4 + (1121./360.) * i6e2,
            (289. / 4.) * e4 - (1955. / 6.) * e6 + (83551. / 144.) * e8 -
            (1. / 4.) * i4 + (1. / 6.) * i6 - (1. / 20.) * i8 +
            5. * i2e2 - (1627. / 16.) * i2e4 + (53245. / 144.) * i2e6 -
            (209. / 48.) * i4e2 + (11023. / 192.) * i4e4 + (935. / 576.) * i6e2
        ),
        (
            '3n',
            3. * orbital_freq,
            (2855./18.) * e6 - (1709915./2304.) * e8 + (1./32.) * i6 - (1./64.) * i8 +
            (1237./32.) * i2e4 - (374291./1152.) * i2e6 + (165./64.) * i4e2 - (32975./768.) * i4e4 - (15./8.) * i6e2,

            (9917. / 32.) * e6 - (2290303. / 1536.) * e8 - (1. / 32.) * i6 + (1. / 64.) * i8 +
            (1075. / 32.) * i2e4 - (45769. / 96.) * i2e6 - (9. / 16.) * i4e2 - (5375. / 192.) * i4e4 + (17./32.) * i6e2
        ),
        (
            '4n',
            4. * orbital_freq,
            (2592379./4608.) * e8 + (1./512.) * i8 + (369653./2304.) * i2e6 + (1815./128.) * i4e4 + (49./128.) * i6e2,
            (2556797. / 2304.) * e8 - (1. / 256.) * i8 + (86093. / 576.) * i2e6 - (81. / 64.) * i4e4 - (49./128.) * i6e2
        )
    )

    heating_coeffs = list()
    ztorque_coeffs = list()
    freqs = list()
    modes = list()
    mode_names = list()

    for mode_name, mode, coeff, torque_coeff in modes_coeffs:
        sgn = np.sign(mode)
        freq = np.abs(mode)

        heating_coeff = freq * coeff
        ztorque_coeff = torque_coeff * sgn

        freq_too_low = freq < MODE_ZERO_TOL
        freq[freq_too_low] = 0.
        mode[freq_too_low] = 0.
        heating_coeff[freq_too_low] = 0.
        ztorque_coeff[freq_too_low] = 0.

        freqs.append(freq)
        modes.append(mode)
        mode_names.append(mode_name)
        heating_coeffs.append(heating_coeff)
        ztorque_coeffs.append(ztorque_coeff)

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return mode_names, modes, freqs, heating_coeffs, ztorque_coeffs

@njit
def nsr_modes(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray):
    """ Non-Synchronous Rotation Tidal Modes (for heating and torque) for l=2, e^2, I^2

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
            'n',
            orbital_freq,
            (3. / 4.) * e2,
            0.
        ),
        (
            'o',
            spin_freq,
            (1. / 2.) * i2,
            -1.
        ),
        (
            '2n-o',
            2. * orbital_freq - spin_freq,
            (1. / 2.) * i2,
            1.
        ),
        (
            'n-2o',
            orbital_freq - 2. * spin_freq,
            (1. / 8.) * e2,
            2.
        ),
        (
            '2n-2o',
            2. * orbital_freq - 2. * spin_freq,
            (1. / 2.) - (5. / 2.) * e2 - (1. / 2.) * i2,
            2.
        ),
        (
            '3n-2o',
            3. * orbital_freq - 2. * spin_freq,
            (49. / 8.) * e2,
            2.
        ),
    )

    heating_coeffs = list()
    ztorque_coeffs = list()
    freqs = list()
    modes = list()
    mode_names = list()

    for mode_name, mode, coeff, torque_multi in modes_coeffs:
        sgn = np.sign(mode)
        freq = np.abs(mode)

        heating_coeff = freq * coeff
        ztorque_coeff = coeff * torque_multi * sgn

        freq_too_low = freq < MODE_ZERO_TOL
        freq[freq_too_low] = 0.
        mode[freq_too_low] = 0.
        heating_coeff[freq_too_low] = 0.
        ztorque_coeff[freq_too_low] = 0.

        freqs.append(freq)
        modes.append(mode)
        mode_names.append(mode_name)
        heating_coeffs.append(heating_coeff)
        ztorque_coeffs.append(ztorque_coeff)

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return mode_names, modes, freqs, heating_coeffs, ztorque_coeffs


@njit
def nsr_modes_4(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray):
    """ Non-Synchronous Rotation Tidal Modes (for heating and torque) for l=2, e^4, I^4

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
            'n',
            orbital_freq,
            (3. / 4.) * e2 + (27. / 16.) * e4 - (9. / 4.) * i2e2,
            0.
        ),
        (
            '2n',
            2. * orbital_freq,
            (27. / 16.) * e4 + (3. / 16.) * i4,
            0.
        ),
        (
            'o',
            spin_freq,
            (1. / 2.) * i2 - (2. / 3.) * i4 + (3. / 2.) * i2e2,
            -1.
        ),
        (
            '2o',
            2. * spin_freq,
            (1. / 8.) * i4,
            -2.
        ),
        (
            'n-o',
            orbital_freq - spin_freq,
            (5. / 4.) * i2e2,
            1.
        ),
        (
            '2n-o',
            2. * orbital_freq - spin_freq,
            (1. / 2.) * i2 - (5. / 12.) * i4 - (5. / 2.) * i2e2,
            1.
        ),
        (
            '3n-o',
            3. * orbital_freq - spin_freq,
            (49. / 8.) * i2e2,
            1.
        ),
        (
            'n-2o',
            orbital_freq - 2. * spin_freq,
            (1. / 8.) * e2 - (1. / 32.) * e4 - (1. / 8.) * i2e2,
            2.
        ),
        (
            '2n-2o',
            2. * orbital_freq - 2. * spin_freq,
            (1. / 2.) - (5. / 2.) * e2 + (63. / 16.) * e4 - (1. / 2.) * i2 + (11. / 48.) * i4 + (5. / 2.) * i2e2,
            2.
        ),
        (
            '3n-2o',
            3. * orbital_freq - 2. * spin_freq,
            (49. / 8.) * e2 - (861. / 32.) * e4 - (49. / 8.) * i2e2,
            2.
        ),
        (
            '4n-2o',
            4. * orbital_freq - 2. * spin_freq,
            (289. / 8.) * e4,
            2.
        ),
        (
            '-n-o',
            -1. * orbital_freq - spin_freq,
            (9. / 8.) * i2e2,
            1.
        )
    )

    heating_coeffs = list()
    ztorque_coeffs = list()
    freqs = list()
    modes = list()
    mode_names = list()

    for mode_name, mode, coeff, torque_multi in modes_coeffs:
        sgn = np.sign(mode)
        freq = np.abs(mode)

        heating_coeff = freq * coeff
        ztorque_coeff = coeff * torque_multi * sgn

        freq_too_low = freq < MODE_ZERO_TOL
        freq[freq_too_low] = 0.
        mode[freq_too_low] = 0.
        heating_coeff[freq_too_low] = 0.
        ztorque_coeff[freq_too_low] = 0.

        freqs.append(freq)
        modes.append(mode)
        mode_names.append(mode_name)
        heating_coeffs.append(heating_coeff)
        ztorque_coeffs.append(ztorque_coeff)

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return mode_names, modes, freqs, heating_coeffs, ztorque_coeffs


@njit
def nsr_modes_6(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray):
    """ Non-Synchronous Rotation Tidal Modes (for heating and torque) for l=2, e^6, I^6

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
            'n',
            orbital_freq,
            (3./4.) * e2 + (27./16.) * e4 + (765./256.) * e6 - (9./4.) * i2e2 - (81./16.) * i2e4 + (159./64.) * i4e2,
            0.
        ),
        (
            '2n',
            2. * orbital_freq,
            (27./16.) * e4 + (21./8.) * e6 - (81./16.) * i2e4 + (3./16.) * i4 - (15./16.) * i4e2 - (1./8.) * i6,
            0.
        ),
        (
            '3n',
            3. * orbital_freq,
            (2809./768.) * e6 + (147./64.) * i4e2,
            0.
        ),
        (
            'o',
            spin_freq,
            (1./2.) * i2 - (2./3.) * i4 + (16./45.) * i6 + (3./2.) * i2e2 + 3. * i2e4 - 2. * i4e2,
            -1.
        ),
        (
            '2o',
            2. * spin_freq,
            (1./8.) * i4 - (1./12.) * i6 + (3./8.) * i4e2,
            -2.
         ),
        (
            'n-o',
            orbital_freq - spin_freq,
            (5./4.) * i2e2 + (5./2.) * i2e4 - (77./48.) * i4e2,
            1.
        ),
        (
            '2n-o',
            2. * orbital_freq - spin_freq,
            (1./2.) * i2 - (5./12.) * i4 + (227./1440.) * i6 - (5./2.) * i2e2 + (207./32.) * i2e4 + (25./12.) * i4e2,
            1.
        ),
        (
            '3n-o',
            3. * orbital_freq - spin_freq,
            (49./8.) * i2e2 - (861./32.) * i2e4 - (245./48.) * i4e2,
            1.
        ),
        (
            '4n-o',
            4.* orbital_freq - spin_freq,
            (289./8.) * i2e4,
            1.
        ),
        (
            'n-2o',
            orbital_freq - 2. * spin_freq,
            (1./8.) * e2 - (1./32.) * e4 + (13./1536.) * e6 - (1./8.) * i2e2 + (1./32.) * i2e4 + (65./192.) * i4e2,
            2.
        ),
        (
            '2n-2o',
            2. * orbital_freq - 2. * spin_freq,
            (1./2.) - (5./2.) * e2 + (63./16.) * e4 - (79./36.) * e6 - (1./2.) * i2 + (11./48.) * i4 - (23./360.) * i6
            + (5./2.) * i2e2 - (63./16.) * i2e4 - (55./48.) * i4e2,
            2.
        ),
        (
            '3n-2o',
            3. * orbital_freq - 2. * spin_freq,
            (49./8.) * e2 - (861./32.) * e4 + (21975./512.) * e6 - (49./8.) * i2e2 + (861./32.) * i2e4
            + (539./192.) * i4e2,
            2.
        ),
        (
            '4n-2o',
            4. * orbital_freq - 2. * spin_freq,
            (289./8.) * e4 - (1955./12.) * e6 - (289./8.) * i2e4,
            2.
        ),
        (
            '5n-2o',
            5. * orbital_freq - 2. * spin_freq,
            (714025./4608.) * e6,
            2.
        ),
        (
            '-n-o',
            -1. * orbital_freq - spin_freq,
            (9./8.) * i2e2 + (81./32.) * i2e4 - (3./2.) * i4e2,
            1.
        ),
        (
            '-n-2o',
            -1. * orbital_freq - 2. * spin_freq,
            (1./4608.) * e6 + (9./32.) * i4e2,
            2.
        ),
        (
            '-2n-o',
            -2. * orbital_freq - 1. * spin_freq,
            (1./32.) * i6 + (81./32.) * i2e4,
            1.
        )
    )

    heating_coeffs = list()
    ztorque_coeffs = list()
    freqs = list()
    modes = list()
    mode_names = list()

    for mode_name, mode, coeff, torque_multi in modes_coeffs:
        sgn = np.sign(mode)
        freq = np.abs(mode)

        heating_coeff = freq * coeff
        ztorque_coeff = coeff * torque_multi * sgn

        freq_too_low = freq < MODE_ZERO_TOL
        freq[freq_too_low] = 0.
        mode[freq_too_low] = 0.
        heating_coeff[freq_too_low] = 0.
        ztorque_coeff[freq_too_low] = 0.

        freqs.append(freq)
        modes.append(mode)
        mode_names.append(mode_name)
        heating_coeffs.append(heating_coeff)
        ztorque_coeffs.append(ztorque_coeff)

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return mode_names, modes, freqs, heating_coeffs, ztorque_coeffs

@njit
def nsr_modes_8(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray):
    """ Non-Synchronous Rotation Tidal Modes (for heating and torque) for l=2, e^8, I^8

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
    e8 = eccentricity**8
    i2 = inclination**2
    i4 = inclination**4
    i6 = inclination**6
    i8 = inclination**8
    i2e2 = e2 * i2
    i2e4 = i2 * e4
    i4e2 = i4 * e2
    i6e2 = i6 * e2
    i4e4 = i4 * e4
    i2e6 = i2 * e6

    modes_coeffs = (
        (
            'n',
            orbital_freq,
            (3./4.) * e2 + (27./16.) * e4 + (765./256.) * e6 + (28211./6144.) * e8 - (9./4.) * i2e2 - (81./16.) * i2e4 -
            (2295./256.) * i2e6 + (159./64.) * i4e2 + (1401./256.) * i4e4 - (201./160.) * i6e2,
            0.
        ),
        (
            '2n',
            2. * orbital_freq,
            (27./16.) * e4 + (21./8.) * e6 + (1661./384.) * e8 - (81./16.) * i2e4 - (63./8.) * i2e6 + (3./16.) * i4 -
            (15./16.) * i4e2 + (891./128.) * i4e4 - (1./8.) * i6 + (5./8.) * i6e2 + (3./80.) * i8,
            0.
        ),
        (
            '3n',
            3. * orbital_freq,
            (2809./768.) * e6 + (6943./2048.) * e8 - (2809./256.) * i2e6 + (147./64.) * i4e2 - (2583./256.) * i4e4 -
            (49./32.) * i6e2,
            0.
        ),
        (
            '4n',
            4. * orbital_freq,
            (5929./768.) * e8 + (867. / 64.) * i4e4,
            0.
        ),
        (
            'o',
            spin_freq,
            (1./2.) * i2 - (2./3.) * i4 + (16./45.) * i6 - (32./325.) * i8 + (3./2.) * i2e2 + 3. * i2e4 + 5. * i2e6 -
            2. * i4e2 - 4. * i4e4 + (16./15.) * i6e2,
            -1.
        ),
        (
            '2o',
            2. * spin_freq,
            (1./8.) * i4 - (1./12.) * i6 + (1./40.) * i8 + (3./8.) * i4e2 + (3./4.) * i4e4 - (1./4.) * i6e2,
            -2.
         ),
        (
            'n-o',
            orbital_freq - spin_freq,
            (5./4.) * i2e2 + (5./2.) * i2e4 + (3449./768.) * i2e6 - (77./48.) * i4e2 - (643./192.) * i4e4 +
            (967./1152.) * i6e2,
            1.
        ),
        (
            '2n-o',
            2. * orbital_freq - spin_freq,
            (1./2.) * i2 - (5./12.) * i4 + (227./1440.) * i6 - (145./4032.) * i8 - (5./2.) * i2e2 + (207./32.) * i2e4 +
            (251./144.) * i2e6 + (25./12.) * i4e2 - (213./32.) * i4e4 - (227./288.) * i6e2,
            1.
        ),
        (
            '3n-o',
            3. * orbital_freq - spin_freq,
            (49./8.) * i2e2 - (861./32.) * i2e4 + (1549./32.) * i2e6 - (245./48.) * i4e2 + (1435./64.) * i4e4 +
            (11123./5760.) * i6e2,
            1.
        ),
        (
            '4n-o',
            4.* orbital_freq - spin_freq,
            (289./8.) * i2e4 - (1955./12.) * i2e6 - (1445./48.) * i4e4,
            1.
        ),
        (
            '5n-o',
            5.* orbital_freq - spin_freq,
            (714025./4608.) * i2e6,
            1.
        ),
        (
            'n-2o',
            orbital_freq - 2. * spin_freq,
            (1./8.) * e2 - (1./32.) * e4 + (13./1536.) * e6 + (305./36864.) * e8 - (1./8.) * i2e2 + (1./32.) * i2e4 -
            (13./1536.) * i2e6 + (65./192.) * i4e2 + (475./768.) * i4e4 - (293./1440.) * i6e2,
            2.
        ),
        (
            '2n-2o',
            2. * orbital_freq - 2. * spin_freq,
            (1./2.) - (5./2.) * e2 + (63./16.) * e4 - (79./36.) * e6 + (3697./4608.) * e8 - (1./2.) * i2 +
            (11./48.) * i4 - (23./360.) * i6 + (1957./161280.) * i8 + (5./2.) * i2e2 - (63./16.) * i2e4 +
            (79./36.) * i2e6 - (55./48.) * i4e2 + (39./16.) * i4e4 + (23./72.) * i6e2,
            2.
        ),
        (
            '3n-2o',
            3. * orbital_freq - 2. * spin_freq,
            (49./8.) * e2 - (861./32.) * e4 + (21975./512.) * e6 - (135771./4096.) * e8 - (49./8.) * i2e2 +
            (861./32.) * i2e4 - (21975./512.) * i2e6 + (539./192.) * i4e2 - (3157./256.) * i4e4 - (1127./1440.) * i6e2,
            2.
        ),
        (
            '4n-2o',
            4. * orbital_freq - 2. * spin_freq,
            (289./8.) * e4 - (1955./12.) * e6 + (83551./288.) * e8 - (289./8.) * i2e4 + (1955./12.) * i2e6 +
            (3179./192.) * i4e4,
            2.
        ),
        (
            '5n-2o',
            5. * orbital_freq - 2. * spin_freq,
            (714025./4608.) * e6 - (27483625./36864.) * e8 - (714025./4608.) * i2e6,
            2.
        ),
        (
            '6n-2o',
            6. * orbital_freq - 2. * spin_freq,
            (284089./512.) * e8,
            2.
        ),
        (
            '-n-o',
            -1. * orbital_freq - spin_freq,
            (9./8.) * i2e2 + (81./32.) * i2e4 + (1291./288.) * i2e6 - (3./2.) * i4e2 - (27./8.) * i4e4 +
            (517./640.) * i6e2,
            1.
        ),
        (
            '-2n-o',
            -2. * orbital_freq - 1. * spin_freq,
            (1./32.) * i6 - (1./64.) * i8 + (81./32.) * i2e4 + (63./16.) * i2e6 - (27./8.) * i4e4 - (5./32.) * i6e2,
            1.
        ),
        (
            '-3n-o',
            -3. * orbital_freq - 1. * spin_freq,
            (2809./512.) * i2e6 + (49./128.) * i6e2,
            1.
        ),
        (
            '-n-2o',
            -1. * orbital_freq - 2. * spin_freq,
            (1./4608.) * e6 + (11./36864.) * e8 - (1./4608.) * i2e6 + (9./32.) * i4e2 + (81./128.) * i4e4 -
            (3./16.) * i6e2,
            2.
        ),
        (
            '-2n-2o',
            -2. * orbital_freq - 2. * spin_freq,
            (1./1152.) * e8 + (81./128.) * i4e4 + (1./512) * i8,
            2.
        )
    )

    heating_coeffs = list()
    ztorque_coeffs = list()
    freqs = list()
    modes = list()
    mode_names = list()

    for mode_name, mode, coeff, torque_multi in modes_coeffs:
        sgn = np.sign(mode)
        freq = np.abs(mode)

        heating_coeff = freq*coeff
        ztorque_coeff = coeff*torque_multi*sgn

        freq_too_low = freq < MODE_ZERO_TOL
        freq[freq_too_low] = 0.
        mode[freq_too_low] = 0.
        heating_coeff[freq_too_low] = 0.
        ztorque_coeff[freq_too_low] = 0.

        # if freq < MODE_ZERO_TOL:
        #     freq = 0.
        #     heating_coeff = 0.
        #     ztorque_coeff = 0.

        freqs.append(freq)
        modes.append(mode)
        mode_names.append(mode_name)
        heating_coeffs.append(heating_coeff)
        ztorque_coeffs.append(ztorque_coeff)

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return mode_names, modes, freqs, heating_coeffs, ztorque_coeffs


@njit
def nsr_modes_10(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray):
    """ Non-Synchronous Rotation Tidal Modes (for heating and torque) for l=2, e^10, I^10

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
    e8 = eccentricity**8
    e10 = eccentricity**10
    i2 = inclination**2
    i4 = inclination**4
    i6 = inclination**6
    i8 = inclination**8
    i10 = inclination**10
    i2e2 = i2 * e2
    i2e4 = i2 * e4
    i2e6 = i2 * e6
    i2e8 = i2 * e8
    i4e2 = i4 * e2
    i4e4 = i4 * e4
    i4e6 = i4 * e6
    i6e2 = i6 * e2
    i6e4 = i6 * e4
    i8e2 = i8 * e2

    modes_coeffs = (
        (
            'n',
            orbital_freq,
            793*i8e2/2240 - 1759*i6e4/640 - 201*i6e2/160 + 29845*i4e6/3072 + 1401*i4e4/256 + 159*i4e2/64
            - 28211*i2e8/2048 - 2295*i2e6/256 - 81*i2e4/16 - 9*i2e2/4 + 1064887*e10/163840 + 28211*e8/6144
            + 765*e6/256 + 27*e4/16 + 3*e2/4,
            0.
        ),
        (
            '2n',
            2. * orbital_freq,
            - 17*i10/2520 - 3*i8e2/16 + 3*i8/80 - 1197*i6e4/320 + 5*i6e2/8 - i6/8 + 1483*i4e6/192 + 891*i4e4/128
            - 15*i4e2/16 + 3*i4/16 - 1661*i2e8/128 - 63*i2e6/8 - 81*i2e4/16 + 3919*e10/640 + 1661*e8/384 + 21*e6/8
            + 27*e4/16,
            0.
        ),
        (
            '3n',
            3. * orbital_freq,
            + 147*i8e2/320 + 861*i6e4/128 - 49*i6e2/32 + 343843*i4e6/12288 - 2583*i4e4/256 + 147*i4e2/64
            - 20829*i2e8/2048 - 2809*i2e6/256 + 2006627*e10/327680 + 6943*e8/2048 + 2809*e6/768,
            0.
        ),
        (
            '4n',
            4. * orbital_freq,
            - 289*i6e4/32 - 1955*i4e6/32 + 867*i4e4/64 - 5929*i2e8/256 + 3311*e10/1280 + 5929*e8/768,
            0.
        ),
        (
            '5n',
            5. * orbital_freq,
            + 714025*i4e6/12288 + 1047843*e10/65536,
            0.
        ),
        (
            'o',
            spin_freq,
            + 256*i10/14175 - 32*i8e2/105 - 32*i8/315 + 32*i6e4/15 + 16*i6e2/15 + 16*i6/45 - 20*i4e6/3 - 4*i4e4
            - 2*i4e2 - 2*i4/3 + 15*i2e8/2 + 5*i2e6 + 3*i2e4 + 3*i2e2/2 + i2/2,
            -1.
        ),
        (
            '2o',
            2. * spin_freq,
            - 17*i10/3780 + 3*i8e2/40 + i8/40 - i6e4/2 - i6e2/4 - i6/12 + 5*i4e6/4 + 3*i4e4/4 + 3*i4e2/8 + i4/8,
            -2.
         ),
        (
            'n-o',
            orbital_freq - spin_freq,
            - 19157*i8e2/80640 + 8249*i6e4/4608 + 967*i6e2/1152 - 55145*i4e6/9216 - 643*i4e4/192 - 77*i4e2/48
            + 63551*i2e8/9216 + 3449*i2e6/768 + 5*i2e4/2 + 5*i2e2/4,
            1.
        ),
        (
            '2n-o',
            2. * orbital_freq - spin_freq,
            + 40277*i10/7257600 + 725*i8e2/4032 - 145*i8/4032 + 3893*i6e4/1280 - 227*i6e2/288 + 227*i6/1440
            - 739*i4e6/216 - 213*i4e4/32 + 25*i4e2/12 - 5*i4/12 + 33595*i2e8/4608 + 251*i2e6/144 + 207*i2e4/32
            - 5*i2e2/2 + i2/2,
            1.
        ),
        (
            '3n-o',
            3. * orbital_freq - spin_freq,
            - 1015*i8e2/2304 - 65149*i6e4/7680 + 11123*i6e2/5760 - 132347*i4e6/3072 + 1435*i4e4/64 - 245*i4e2/48
            - 57471*i2e8/2048 + 1549*i2e6/32 - 861*i2e4/32 + 49*i2e2/8,
            1.
        ),
        (
            '4n-o',
            4.* orbital_freq - spin_freq,
            + 65603*i6e4/5760 + 9775*i4e6/72 - 1445*i4e4/48 + 1390177*i2e8/4608 - 1955*i2e6/12 + 289*i2e4/8,
            1.
        ),
        (
            '5n-o',
            5.* orbital_freq - spin_freq,
            - 3570125*i4e6/27648 - 27483625*i2e8/36864 + 714025*i2e6/4608,
            1.
        ),
        (
            '6n-o',
            6.* orbital_freq - spin_freq,
            284089*i2e8/512,
            1.
        ),
        (
            'n-2o',
            orbital_freq - 2. * spin_freq,
            + 7649*i8e2/129024 - 2407*i6e4/5760 - 293*i6e2/1440 + 41453*i4e6/36864 + 475*i4e4/768 + 65*i4e2/192
            - 305*i2e8/36864 - 13*i2e6/1536 + i2e4/32 - i2e2/8 + 1733*e10/327680 + 305*e8/36864 + 13*e6/1536
            - e4/32 + e2/8,
            2.
        ),
        (
            '2n-2o',
            2. * orbital_freq - 2. * spin_freq,
            - 12107*i10/7257600 - 1957*i8e2/32256 + 1957*i8/161280 - 37*i6e4/40 + 23*i6e2/72 - 23*i6/360
            - 37*i4e6/1728 + 39*i4e4/16 - 55*i4e2/48 + 11*i4/48 - 3697*i2e8/4608 + 79*i2e6/36 - 63*i2e4/16 + 5*i2e2/2
            - i2/2 - 11041*e10/38400 + 3697*e8/4608 - 79*e6/36 + 63*e4/16 - 5*e2/2 + 1/2,
            2.
        ),
        (
            '3n-2o',
            3. * orbital_freq - 2. * spin_freq,
            + 13699*i8e2/92160 + 6601*i6e4/1920 - 1127*i6e2/1440 + 86193*i4e6/4096 - 3157*i4e4/256 + 539*i4e2/192
            + 135771*i2e8/4096 - 21975*i2e6/512 + 861*i2e4/32 - 49*i2e2/8 + 5568309*e10/327680 - 135771*e8/4096
            + 21975*e6/512 - 861*e4/32 + 49*e2/8,
            2.
        ),
        (
            '4n-2o',
            4. * orbital_freq - 2. * spin_freq,
            - 6647*i6e4/1440 - 21505*i4e6/288 + 3179*i4e4/192 - 83551*i2e8/288 + 1955*i2e6/12 - 289*i2e4/8
            - 134209*e10/480 + 83551*e8/288 - 1955*e6/12 + 289*e4/8,
            2.
        ),
        (
            '5n-2o',
            5. * orbital_freq - 2. * spin_freq,
            + 7854275*i4e6/110592 + 27483625*i2e8/36864 - 714025*i2e6/4608 + 587225375*e10/393216 - 27483625*e8/36864
            + 714025*e6/4608,
            2.
        ),
        (
            '6n-2o',
            6. * orbital_freq - 2. * spin_freq,
            - 284089*i2e8/512 - 7369791*e10/2560 + 284089*e8/512,
            2.
        ),
        (
            '7n-2o',
            7. * orbital_freq - 2. * spin_freq,
            52142352409*e10/29491200,
            2.
        ),
        (
            '-n-o',
            -1. * orbital_freq - spin_freq,
            - 2083*i8e2/8960 + 4603*i6e4/2560 + 517*i6e2/640 - 165245*i4e6/27648 - 27*i4e4/8 - 3*i4e2/2
            + 126955*i2e8/18432 + 1291*i2e6/288 + 81*i2e4/32 + 9*i2e2/8,
            1.
        ),
        (
            '-2n-o',
            -2. * orbital_freq - 1. * spin_freq,
            + 9*i10/2560 + 5*i8e2/64 - i8/64 + 2619*i6e4/1280 - 5*i6e2/32 + i6/32 - 21*i4e6/4 - 27*i4e4/8
            + 14951*i2e8/2304 + 63*i2e6/16 + 81*i2e4/32,
            1.
        ),
        (
            '-3n-o',
            -3. * orbital_freq - 1. * spin_freq,
            - 49*i8e2/256 - 861*i6e4/512 + 49*i6e2/128 - 2809*i4e6/384 + 20829*i2e8/4096 + 2809*i2e6/512,
            1.
        ),
        (
            '-4n-o',
            -4. * orbital_freq - 1. * spin_freq,
            + 289*i6e4/128 + 5929*i2e8/512,
            1.
        ),
        (
            '-n-2o',
            -1. * orbital_freq - 2. * spin_freq,
            + 581*i8e2/10240 - 27*i6e4/64 - 3*i6e2/16 + 123941*i4e6/110592 + 81*i4e4/128 + 9*i4e2/32 - 11*i2e8/36864
            - i2e6/4608 + 619*e10/1966080 + 11*e8/36864 + e6/4608,
            2.
        ),
        (
            '-2n-2o',
            -2. * orbital_freq - 2. * spin_freq,
            - i10/1536 - 5*i8e2/512 + i8/512 - 27*i6e4/64 + 63*i4e6/64 + 81*i4e4/128 - i2e8/1152 + 7*e10/5760 + e8/1152,
            2.
        ),
        (
            '-3n-2o',
            -3. * orbital_freq - 2. * spin_freq,
            + 49*i8e2/2048 + 2809*i4e6/2048 + 6561*e10/3276800,
            2.
        )
    )

    heating_coeffs = list()
    ztorque_coeffs = list()
    freqs = list()
    modes = list()
    mode_names = list()

    for mode_name, mode, coeff, torque_multi in modes_coeffs:
        sgn = np.sign(mode)
        freq = np.abs(mode)

        heating_coeff = freq * coeff
        ztorque_coeff = coeff * torque_multi * sgn

        freq_too_low = freq < MODE_ZERO_TOL
        freq[freq_too_low] = 0.
        mode[freq_too_low] = 0.
        heating_coeff[freq_too_low] = 0.
        ztorque_coeff[freq_too_low] = 0.

        freqs.append(freq)
        modes.append(mode)
        mode_names.append(mode_name)
        heating_coeffs.append(heating_coeff)
        ztorque_coeffs.append(ztorque_coeff)

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return mode_names, modes, freqs, heating_coeffs, ztorque_coeffs

@njit
def nsr_modes_12(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray):
    """ Non-Synchronous Rotation Tidal Modes (for heating and torque) for l=2, e^12, I^12

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
    e8 = eccentricity**8
    e10 = eccentricity**10
    e12 = eccentricity**12
    i2 = inclination**2
    i4 = inclination**4
    i6 = inclination**6
    i8 = inclination**8
    i10 = inclination**10
    i12 = inclination**12
    i2e2 = i2 * e2
    i2e4 = i2 * e4
    i2e6 = i2 * e6
    i2e8 = i2 * e8
    i2e10 = i2 * e10
    i4e2 = i4 * e2
    i4e4 = i4 * e4
    i4e6 = i4 * e6
    i4e8 = i4 * e8
    i6e2 = i6 * e2
    i6e4 = i6 * e4
    i6e6 = i6 * e6
    i8e2 = i8 * e2
    i8e4 = i8 * e4
    i10e2 = i10 * e2

    modes_coeffs = (
        (
            'n',
            orbital_freq,
            -3161*i10e2/50400 + 6927*i8e4/8960 + 793*i8e2/2240 - 22501*i6e6/4608 - 1759*i6e4/640 - 201*i6e2/160
            + 61137*i4e8/4096 + 29845*i4e6/3072 + 1401*i4e4/256 + 159*i4e2/64 - 3194661*i2e10/163840 - 28211*i2e8/2048
            - 2295*i2e6/256 - 81*i2e4/16 - 9*i2e2/4 + 85553819*e12/9830400 + 1064887*e10/163840 + 28211*e8/6144
            + 765*e6/256 + 27*e4/16 + 3*e2/4,
            0.
        ),
        (
            '2n',
            2. * orbital_freq,
            31*i12/37800 + 17*i10e2/504 - 17*i10/2520 + 4797*i8e4/4480 - 3*i8e2/16 + 3*i8/80 - 5399*i6e6/1440
            - 1197*i6e4/320 + 5*i6e2/8 - i6/8 + 58543*i4e8/4096 + 1483*i4e6/192 + 891*i4e4/128 - 15*i4e2/16 + 3*i4/16
            - 11757*i2e10/640 - 1661*i2e8/128 - 63*i2e6/8 - 81*i2e4/16 + 505601*e12/61440 + 3919*e10/640 + 1661*e8/384
            + 21*e6/8 + 27*e4/16,
            0.
        ),
        (
            '3n',
            3. * orbital_freq,
            -119*i10e2/1440 - 2583*i8e4/1280 + 147*i8e2/320 - 1539439*i6e6/92160 + 861*i6e4/128 - 49*i6e2/32
            - 46277*i4e8/32768 + 343843*i4e6/12288 - 2583*i4e4/256 + 147*i4e2/64 - 6019881*i2e10/327680
            - 20829*i2e8/2048 - 2809*i2e6/256 + 30355211*e12/3932160 + 2006627*e10/327680 + 6943*e8/2048 + 2809*e6/768,
            0.
        ),
        (
            '4n',
            4. * orbital_freq,
            867*i8e4/320 + 1955*i6e6/48 - 289*i6e4/32 + 411281*i4e8/3072 - 1955*i4e6/32 + 867*i4e4/64
            - 9933*i2e10/1280 - 5929*i2e8/256 + 708871*e12/76800 + 3311*e10/1280 + 5929*e8/768,
            0.
        ),
        (
            '5n',
            5. * orbital_freq,
            -714025*i6e6/18432 - 27483625*i4e8/98304 + 714025*i4e6/12288 - 3143529*i2e10/65536 - 982439*e12/262144
            + 1047843*e10/65536,
            0.
        ),
        (
            '6n',
            6. * orbital_freq,
            852267*i4e8/4096 + 10029889*e12/307200,
            0.
        ),
        (
            'o',
            spin_freq,
            -1024*i12/467775 + 256*i10e2/4725 + 256*i10/14175 - 64*i8e4/105 - 32*i8e2/105 - 32*i8/315 + 32*i6e6/9
            + 32*i6e4/15 + 16*i6e2/15 + 16*i6/45 - 10*i4e8 - 20*i4e6/3 - 4*i4e4 - 2*i4e2 - 2*i4/3 + 21*i2e10/2
            + 15*i2e8/2 + 5*i2e6 + 3*i2e4 + 3*i2e2/2 + i2/2,
            -1.
        ),
        (
            '2o',
            2. * spin_freq,
            31*i12/56700 - 17*i10e2/1260 - 17*i10/3780 + 3*i8e4/20 + 3*i8e2/40 + i8/40 - 5*i6e6/6 - i6e4/2 - i6e2/4
            - i6/12 + 15*i4e8/8 + 5*i4e6/4 + 3*i4e4/4 + 3*i4e2/8 + i4/8,
            -2.
         ),
        (
            'n-o',
            orbital_freq - spin_freq,
            6971*i10e2/165888 - 165163*i8e4/322560 - 19157*i8e2/80640 + 1764043*i6e6/552960 + 8249*i6e4/4608
            + 967*i6e2/1152 - 2032717*i4e8/221184 - 55145*i4e6/9216 - 643*i4e4/192 - 77*i4e2/48 + 1598197*i2e10/163840
            + 63551*i2e8/9216 + 3449*i2e6/768 + 5*i2e4/2 + 5*i2e2/4,
            1.
        ),
        (
            '2n-o',
            2. * orbital_freq - spin_freq,
            -59123*i12/95800320 - 40277*i10e2/1451520 + 40277*i10/7257600 - 14291*i8e4/17920 + 725*i8e2/4032
            - 145*i8/4032 + 54643*i6e6/25920 + 3893*i6e4/1280 - 227*i6e2/288 + 227*i6/1440 - 257669*i4e8/27648
            - 739*i4e6/216 - 213*i4e4/32 + 25*i4e2/12 - 5*i4/12 + 341669*i2e10/38400 + 33595*i2e8/4608 + 251*i2e6/144
            + 207*i2e4/32 - 5*i2e2/2 + i2/2,
            1.
        ),
        (
            '3n-o',
            3. * orbital_freq - spin_freq,
            281939*i10e2/4147200 + 5945*i8e4/3072 - 1015*i8e2/2304 + 6426533*i6e6/368640 - 65149*i6e4/7680
            + 11123*i6e2/5760 + 170741*i4e8/8192 - 132347*i4e6/3072 + 1435*i4e4/64 - 245*i4e2/48
            + 17156499*i2e10/655360 - 57471*i2e8/2048 + 1549*i2e6/32 - 861*i2e4/32 + 49*i2e2/8,
            1.
        ),
        (
            '4n-o',
            4.* orbital_freq - spin_freq,
            -41905*i8e4/16128 - 88757*i6e6/1728 + 65603*i6e4/5760 - 888871*i4e8/3456 + 9775*i4e6/72 - 1445*i4e4/48
            - 423509*i2e10/1536 + 1390177*i2e8/4608 - 1955*i2e6/12 + 289*i2e4/8,
            1.
        ),
        (
            '5n-o',
            5.* orbital_freq - spin_freq,
            32416735*i6e6/663552 + 137418125*i4e8/221184 - 3570125*i4e6/27648 + 298327981*i2e10/196608
            - 27483625*i2e8/36864 + 714025*i2e6/4608,
            1.
        ),
        (
            '6n-o',
            6.* orbital_freq - spin_freq,
            -1420445*i4e8/3072 - 7369791*i2e10/2560 + 284089*i2e8/512,
            1.
        ),
        (
            '7n-o',
            7.* orbital_freq - spin_freq,
            52142352409*i2e10/29491200,
            1.
        ),
        (
            'n-2o',
            orbital_freq - 2. * spin_freq,
            -305867*i10e2/29030400 + 64927*i8e4/516096 + 7649*i8e2/129024 - 206849*i6e6/276480 - 2407*i6e4/5760
            - 293*i6e2/1440 + 1526749*i4e8/884736 + 41453*i4e6/36864 + 475*i4e4/768 + 65*i4e2/192 - 1733*i2e10/327680
            - 305*i2e8/36864 - 13*i2e6/1536 + i2e4/32 - i2e2/8 + 277229*e12/58982400 + 1733*e10/327680 + 305*e8/36864
            + 13*e6/1536 - e4/32 + e2/8,
            2.
        ),
        (
            '2n-2o',
            2. * orbital_freq - 2. * spin_freq,
            330367*i12/1916006400 + 12107*i10e2/1451520 - 12107*i10/7257600 + 4549*i8e4/20480 - 1957*i8e2/32256
            + 1957*i8/161280 - 4871*i6e6/12960 - 37*i6e4/40 + 23*i6e2/72 - 23*i6/360 + 220055*i4e8/110592
            - 37*i4e6/1728 + 39*i4e4/16 - 55*i4e2/48 + 11*i4/48 + 11041*i2e10/38400 - 3697*i2e8/4608 + 79*i2e6/36
            - 63*i2e4/16 + 5*i2e2/2 - i2/2 + 34471*e12/552960 - 11041*e10/38400 + 3697*e8/4608 - 79*e6/36 + 63*e4/16
            - 5*e2/2 + 1/2,
            2.
        ),
        (
            '3n-2o',
            3. * orbital_freq - 2. * spin_freq,
            -84749*i10e2/4147200 - 80237*i8e4/122880 + 13699*i8e2/92160 - 39313*i6e6/6144 + 6601*i6e4/1920
            - 1127*i6e2/1440 - 456169*i4e8/32768 + 86193*i4e6/4096 - 3157*i4e4/256 + 539*i4e2/192
            - 5568309*i2e10/327680 + 135771*i2e8/4096 - 21975*i2e6/512 + 861*i2e4/32 - 49*i2e2/8
            - 47982879*e12/6553600 + 5568309*e10/327680 - 135771*e8/4096 + 21975*e6/512 - 861*e4/32 + 49*e2/8,
            2.
        ),
        (
            '4n-2o',
            4. * orbital_freq - 2. * spin_freq,
            565573*i8e4/645120 + 8993*i6e6/432 - 6647*i6e4/1440 + 7512571*i4e8/55296 - 21505*i4e6/288 + 3179*i4e4/192
            + 134209*i2e10/480 - 83551*i2e8/288 + 1955*i2e6/12 - 289*i2e4/8 + 8421731*e12/46080 - 134209*e10/480
            + 83551*e8/288 - 1955*e6/12 + 289*e4/8,
            2.
        ),
        (
            '5n-2o',
            5. * orbital_freq - 2. * spin_freq,
            -3284515*i6e6/165888 - 302319875*i4e8/884736 + 7854275*i4e6/110592 - 587225375*i2e10/393216
            + 27483625*i2e8/36864 - 714025*i2e6/4608 - 72670996375*e12/42467328 + 587225375*e10/393216
            - 27483625*e8/36864 + 714025*e6/4608,
            2.
        ),
        (
            '6n-2o',
            6. * orbital_freq - 2. * spin_freq,
            3124979*i4e8/12288 + 7369791*i2e10/2560 - 284089*i2e8/512 + 659870313*e12/102400 - 7369791*e10/2560
            + 284089*e8/512,
            2.
        ),
        (
            '7n-2o',
            7. * orbital_freq - 2. * spin_freq,
            -52142352409*i2e10/29491200 - 140254152605*e12/14155776 + 52142352409*e10/29491200,
            2.
        ),
        (
            '8n-2o',
            8. * orbital_freq - 2. * spin_freq,
            5383010161*e12/1036800,
            2.
        ),
        (
            '-n-o',
            -1. * orbital_freq - spin_freq,
            133907*i10e2/3225600 - 18397*i8e4/35840 - 2083*i8e2/8960 + 5288671*i6e6/1658880 + 4603*i6e4/2560
            + 517*i6e2/640 - 2031247*i4e8/221184 - 165245*i4e6/27648 - 27*i4e4/8 - 3*i4e2/2 + 3833717*i2e10/393216
            + 126955*i2e8/18432 + 1291*i2e6/288 + 81*i2e4/32 + 9*i2e2/8,
            1.
        ),
        (
            '-2n-o',
            -2. * orbital_freq - 1. * spin_freq,
            -463*i12/967680 - 9*i10e2/512 + 9*i10/2560 - 11421*i8e4/17920 + 5*i8e2/64 - i8/64 + 1921*i6e6/720
            + 2619*i6e4/1280 - 5*i6e2/32 + i6/32 - 59801*i4e8/6912 - 21*i4e6/4 - 27*i4e4/8 + 105827*i2e10/11520
            + 14951*i2e8/2304 + 63*i2e6/16 + 81*i2e4/32,
            1.
        ),
        (
            '-3n-o',
            -3. * orbital_freq - 1. * spin_freq,
            441*i10e2/10240 + 861*i8e4/1024 - 49*i8e2/256 + 2427083*i6e6/368640 - 861*i6e4/512 + 49*i6e2/128
            - 6943*i4e8/1024 - 2809*i4e6/384 + 15052983*i2e10/1638400 + 20829*i2e8/4096 + 2809*i2e6/512,
            1.
        ),
        (
            '-4n-o',
            -4. * orbital_freq - 1. * spin_freq,
            -289*i8e4/256 - 1955*i6e6/192 + 289*i6e4/128 - 5929*i4e8/384 + 9933*i2e10/2560 + 5929*i2e8/512,
            1.
        ),
        (
            '-5n-o',
            -5. * orbital_freq - 1. * spin_freq,
            714025*i6e6/73728 + 3143529*i2e10/131072,
            1.
        ),
        (
            '-n-2o',
            -1. * orbital_freq - 2. * spin_freq,
            -737*i10e2/71680 + 5179*i8e4/40960 + 581*i8e2/10240 - 619673*i6e6/829440 - 27*i6e4/64 - 3*i6e2/16
            + 1523515*i4e8/884736 + 123941*i4e6/110592 + 81*i4e4/128 + 9*i4e2/32 - 619*i2e10/1966080 - 11*i2e8/36864
            - i2e6/4608 + 62617*e12/212336640 + 619*e10/1966080 + 11*e8/36864 + e6/4608,
            2.
        ),
        (
            '-2n-2o',
            -2. * orbital_freq - 2. * spin_freq,
            19*i12/184320 + 5*i10e2/1536 - i10/1536 + 2907*i8e4/20480 - 5*i8e2/512 + i8/512 - 21*i6e6/32
            - 27*i6e4/64 + 22429*i4e8/13824 + 63*i4e6/64 + 81*i4e4/128 - 7*i2e10/5760 - i2e8/1152 + 949*e12/691200
            + 7*e10/5760 + e8/1152,
            2.
        ),
        (
            '-3n-2o',
            -3. * orbital_freq - 2. * spin_freq,
            -49*i10e2/6144 - 861*i8e4/8192 + 49*i8e2/2048 - 2809*i6e6/3072 + 20829*i4e8/16384 + 2809*i4e6/2048
            - 6561*i2e10/3276800 + 6561*e12/2621440 + 6561*e10/3276800,
            2.
        ),
        (
            '-4n-2o',
            -4. * orbital_freq - 2. * spin_freq,
            289*i8e4/2048 + 5929*i4e8/2048 + 8*e12/2025,
            2.
        )
    )

    heating_coeffs = list()
    ztorque_coeffs = list()
    freqs = list()
    modes = list()
    mode_names = list()

    for mode_name, mode, coeff, torque_multi in modes_coeffs:
        sgn = np.sign(mode)
        freq = np.abs(mode)

        heating_coeff = freq * coeff
        ztorque_coeff = coeff * torque_multi * sgn

        freq_too_low = freq < MODE_ZERO_TOL
        freq[freq_too_low] = 0.
        mode[freq_too_low] = 0.
        heating_coeff[freq_too_low] = 0.
        ztorque_coeff[freq_too_low] = 0.

        freqs.append(freq)
        modes.append(mode)
        mode_names.append(mode_name)
        heating_coeffs.append(heating_coeff)
        ztorque_coeffs.append(ztorque_coeff)

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return mode_names, modes, freqs, heating_coeffs, ztorque_coeffs


@njit
def nsr_modes_16(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray):
    """ Non-Synchronous Rotation Tidal Modes (for heating and torque) for l=2, e^16, I^16

    These should all have the same shape!!

    :param orbital_freq:
    :param spin_freq:
    :param eccentricity:
    :param inclination:
    :return:
    """
    # Conversion from input_heating_nsr_16.dat
    # Eccentricity max order = 16, Inclination max order = 16.
    # Total Number of Unique Modes = 42

    i2 = inclination**2
    i4 = inclination**4
    i6 = inclination**6
    i8 = inclination**8
    i10 = inclination**10
    i12 = inclination**12
    i14 = inclination**14
    i16 = inclination**16
    e2 = eccentricity**2
    e4 = eccentricity**4
    e6 = eccentricity**6
    e8 = eccentricity**8
    e10 = eccentricity**10
    e12 = eccentricity**12
    e14 = eccentricity**14
    e16 = eccentricity**16
    i2e2 = i2*e2
    i2e4 = i2*e4
    i2e6 = i2*e6
    i2e8 = i2*e8
    i2e10 = i2*e10
    i2e12 = i2*e12
    i2e14 = i2*e14
    i4e2 = i4*e2
    i4e4 = i4*e4
    i4e6 = i4*e6
    i4e8 = i4*e8
    i4e10 = i4*e10
    i4e12 = i4*e12
    i6e2 = i6*e2
    i6e4 = i6*e4
    i6e6 = i6*e6
    i6e8 = i6*e8
    i6e10 = i6*e10
    i8e2 = i8*e2
    i8e4 = i8*e4
    i8e6 = i8*e6
    i8e8 = i8*e8
    i10e2 = i10*e2
    i10e4 = i10*e4
    i10e6 = i10*e6
    i12e2 = i12*e2
    i12e4 = i12*e4
    i14e2 = i14*e2

    modes_coeffs = (
        (
            '9n-o',
            9.*orbital_freq - spin_freq,
            (147483366698529./10276044800.)*i2e14,
            1.
        ),
        (
            '10n-2o',
            10.*orbital_freq - 2.*spin_freq,
            (402063787225./10616832.)*e16,
            2.
        ),
        (
            '8n',
            8.*orbital_freq,
            (5383010161./2764800.)*i4e12 + (31801945561./240844800.)*e16,
            0.
        ),
        (
            '-7n-o',
            -7.*orbital_freq - spin_freq,
            (52142352409./471859200.)*i6e10 + (186702632281./1887436800.)*i2e14,
            1.
        ),
        (
            '9n-2o',
            9.*orbital_freq - 2.*spin_freq,
            (-147483366698529./10276044800.)*i2e14 - (3064748517300717./32883343360.)*e16 + (
                        147483366698529./10276044800.)*e14,
            2.
        ),
        (
            '8n-o',
            8.*orbital_freq - spin_freq,
            (-5383010161./1244160.)*i4e12 - (56914314263./1814400.)*i2e14 + (5383010161./1036800.)*i2e12,
            1.
        ),
        (
            '-6n-2o',
            -6.*orbital_freq - 2.*spin_freq,
            (284089./131072.)*i8e8 + (10029889./819200.)*i4e12 + (531441./40140800.)*e16,
            2.
        ),
        (
            '8n-2o',
            8.*orbital_freq - 2.*spin_freq,
            (59213111771./24883200.)*i4e12 + (56914314263./1814400.)*i2e14 - (5383010161./1036800.)*i2e12 + (
                        578658802849./6773760.)*e16 - (56914314263./1814400.)*e14 + (5383010161./1036800.)*e12,
            2.
        ),
        (
            '7n',
            7.*orbital_freq,
            (-52142352409./117964800.)*i6e10 - (140254152605./37748736.)*i4e12 + (52142352409./78643200.)*i4e10 - (
                        186702632281./943718400.)*i2e14 - (1328179895713./15099494400.)*e16 + (
                        186702632281./2831155200.)*e14,
            0.
        ),
        (
            '7n-o',
            7.*orbital_freq - spin_freq,
            (11836313996843./21233664000.)*i6e10 + (701270763025./84934656.)*i4e12 - (52142352409./35389440.)*i4e10 + (
                        104746088252857./4246732800.)*i2e14 - (140254152605./14155776.)*i2e12 + (
                        52142352409./29491200.)*i2e10,
            1.
        ),
        (
            '-6n-o',
            -6.*orbital_freq - spin_freq,
            (-284089./16384.)*i8e8 - (7369791./40960.)*i6e10 + (284089./8192.)*i6e8 - (10029889./153600.)*i4e12 - (
                        6926229./179200.)*i2e14 + (10029889./204800.)*i2e12,
            1.
        ),
        (
            '-5n-2o',
            -5.*orbital_freq - 2.*spin_freq,
            (-714025./3538944.)*i10e6 - (27483625./9437184.)*i8e8 + (714025./1179648.)*i8e6 - (
                        1047843./262144.)*i6e10 - (2947317./2097152.)*i4e12 + (3143529./524288.)*i4e10 - (
                        244140625./33294385152.)*i2e14 + (2685546875./532710162432.)*e16 + (
                        244140625./33294385152.)*e14,
            2.
        ),
        (
            '7n-2o',
            7.*orbital_freq - 2.*spin_freq,
            (-1199274105407./5308416000.)*i6e10 - (1542795678655./339738624.)*i4e12 + (
                        573565876499./707788800.)*i4e10 - (417304029320899./16986931200.)*i2e14 + (
                        140254152605./14155776.)*i2e12 - (52142352409./29491200.)*i2e10 - (
                        9997999669389767./271790899200.)*e16 + (417304029320899./16986931200.)*e14 - (
                        140254152605./14155776.)*e12 + (52142352409./29491200.)*e10,
            2.
        ),
        (
            '6n-o',
            6.*orbital_freq - spin_freq,
            (-41192905./1032192.)*i8e8 - (557647519./614400.)*i6e10 + (64488203./368640.)*i6e8 - (
                        3339471121./614400.)*i4e12 + (2456597./1024.)*i4e10 - (1420445./3072.)*i4e8 - (
                        613914843./71680.)*i2e14 + (265954103./40960.)*i2e12 - (7369791./2560.)*i2e10 + (
                        284089./512.)*i2e8,
            1.
        ),
        (
            '6n',
            6.*orbital_freq,
            (852267./20480.)*i8e8 + (7369791./10240.)*i6e10 - (284089./2048.)*i6e8 + (6199609931./2457600.)*i4e12 - (
                        22109373./20480.)*i4e10 + (852267./4096.)*i4e8 + (6926229./89600.)*i2e14 - (
                        10029889./102400.)*i2e12 + (65845785./1605632.)*e16 - (2308743./89600.)*e14 + (
                        10029889./307200.)*e12,
            0.
        ),
        (
            '-5n-o',
            -5.*orbital_freq - spin_freq,
            (142805./131072.)*i10e6 + (27483625./1179648.)*i8e8 - (714025./147456.)*i8e6 + (
                        3472622491./31457280.)*i6e10 - (27483625./589824.)*i6e8 + (714025./73728.)*i6e6 + (
                        982439./131072.)*i4e12 - (1047843./32768.)*i4e10 + (213067123621./8323596288.)*i2e14 - (
                        2947317./524288.)*i2e12 + (3143529./131072.)*i2e10,
            1.
        ),
        (
            '6n-2o',
            6.*orbital_freq - 2.*spin_freq,
            (555962173./41287680.)*i8e8 + (56501731./153600.)*i6e10 - (6534047./92160.)*i6e8 + (
                        242955437./81920.)*i4e12 - (27022567./20480.)*i4e10 + (3124979./12288.)*i4e8 + (
                        3055721757./358400.)*i2e14 - (659870313./102400.)*i2e12 + (7369791./2560.)*i2e10 - (
                        284089./512.)*i2e8 + (180544031973./22937600.)*e16 - (3055721757./358400.)*e14 + (
                        659870313./102400.)*e12 - (7369791./2560.)*e10 + (284089./512.)*e8,
            2.
        ),
        (
            '5n-o',
            5.*orbital_freq - spin_freq,
            (1150351397./668860416.)*i10e6 + (3985125625./74317824.)*i8e8 - (103533625./9289728.)*i8e6 + (
                        138128620669./283115520.)*i6e10 - (1247756575./5308416.)*i6e8 + (32416735./663552.)*i6e6 + (
                        365264843291./254803968.)*i4e12 - (3011571571./2359296.)*i4e10 + (137418125./221184.)*i4e8 - (
                        3570125./27648.)*i4e6 + (1637783339303./1189085184.)*i2e14 - (18227432263./10616832.)*i2e12 + (
                        298327981./196608.)*i2e10 - (27483625./36864.)*i2e8 + (714025./4608.)*i2e6,
            1.
        ),
        (
            '-4n-2o',
            -4.*orbital_freq - 2.*spin_freq,
            (5491./737280.)*i12e4 + (1955./9216.)*i10e6 - (289./6144.)*i10e4 + (631199./368640.)*i8e8 - (
                        1955./3072.)*i8e6 + (289./2048.)*i8e4 - (3311./5120.)*i6e10 - (5929./3072.)*i6e8 + (
                        34469153./9953280.)*i4e12 + (9933./10240.)*i4e10 + (5929./2048.)*i4e8 - (8./2025.)*i2e14 - (
                        8./2025.)*i2e12 + (23./4725.)*e16 + (8./2025.)*e14 + (8./2025.)*e12,
            2.
        ),
        (
            '5n',
            5.*orbital_freq,
            (-2427685./1161216.)*i10e6 - (5496725./98304.)*i8e8 + (142805./12288.)*i8e6 - (
                        3141504103./7864320.)*i6e10 + (27483625./147456.)*i6e8 - (714025./18432.)*i6e6 - (
                        74050340731./113246208.)*i4e12 + (641713211./1048576.)*i4e10 - (27483625./98304.)*i4e8 + (
                        714025./12288.)*i4e6 - (13524196093./264241152.)*i2e14 + (2947317./262144.)*i2e12 - (
                        3143529./65536.)*i2e10 + (28491858875./4227858432.)*e16 + (13524196093./792723456.)*e14 - (
                        982439./262144.)*e12 + (1047843./65536.)*e10,
            0.
        ),
        (
            '-4n-o',
            -4.*orbital_freq - spin_freq,
            (-133807./3870720.)*i12e4 - (1173./1024.)*i10e6 + (2601./10240.)*i10e4 - (526171./46080.)*i8e8 + (
                        1955./384.)*i8e6 - (289./256.)*i8e4 - (565093./38400.)*i6e10 + (607483./23040.)*i6e8 - (
                        1955./192.)*i6e6 + (289./128.)*i6e4 - (57428791./3110400.)*i4e12 - (3311./640.)*i4e10 - (
                        5929./384.)*i4e8 + (5450899./414720.)*i2e14 + (11486987./829440.)*i2e12 + (
                        9933./2560.)*i2e10 + (5929./512.)*i2e8,
            1.
        ),
        (
            '5n-2o',
            5.*orbital_freq - 2.*spin_freq,
            (-345788027./668860416.)*i10e6 - (10757090825./594542592.)*i8e8 + (279469385./74317824.)*i8e6 - (
                        2757820247./14155776.)*i6e10 + (126424675./1327104.)*i6e8 - (3284515./165888.)*i6e6 - (
                        800813356187./1019215872.)*i4e12 + (6516062647./9437184.)*i4e10 - (302319875./884736.)*i4e8 + (
                        7854275./110592.)*i4e6 - (6429415592375./4756340736.)*i2e14 + (72670996375./42467328.)*i2e12 - (
                        587225375./393216.)*i2e10 + (27483625./36864.)*i2e8 - (714025./4608.)*i2e6 - (
                        21274828753525./25367150592.)*e16 + (6429415592375./4756340736.)*e14 - (
                        72670996375./42467328.)*e12 + (587225375./393216.)*e10 - (27483625./36864.)*e8 + (
                        714025./4608.)*e6,
            2.
        ),
        (
            '4n-o',
            4.*orbital_freq - spin_freq,
            (-17086547./383201280.)*i12e4 - (15748307./8709120.)*i10e6 + (11640053./29030400.)*i10e4 - (
                        67404683./2903040.)*i8e8 + (283475./24192.)*i8e6 - (41905./16128.)*i8e4 - (
                        1180475./13824.)*i6e10 + (20673629./207360.)*i6e8 - (88757./1728.)*i6e6 + (
                        65603./5760.)*i6e4 - (236062631./1382400.)*i4e12 + (1312291./5760.)*i4e10 - (
                        888871./3456.)*i4e8 + (9775./72.)*i4e6 - (1445./48.)*i4e4 - (1196078171./14515200.)*i2e14 + (
                        90597149./460800.)*i2e12 - (423509./1536.)*i2e10 + (1390177./4608.)*i2e8 - (1955./12.)*i2e6 + (
                        289./8.)*i2e4,
            1.
        ),
        (
            '4n',
            4.*orbital_freq,
            (8959./151200.)*i12e4 + (6647./3024.)*i10e6 - (4913./10080.)*i10e4 + (1166083./46080.)*i8e8 - (
                        391./32.)*i8e6 + (867./320.)*i8e4 + (840647./12800.)*i6e10 - (217949./2560.)*i6e8 + (
                        1955./48.)*i6e6 - (289./32.)*i6e4 + (544861901./5529600.)*i4e12 - (493793./5120.)*i4e10 + (
                        411281./3072.)*i4e8 - (1955./32.)*i4e6 + (867./64.)*i4e4 - (3027367./115200.)*i2e14 - (
                        708871./25600.)*i2e12 - (9933./1280.)*i2e10 - (5929./256.)*i2e8 + (124914751./10321920.)*e16 + (
                        3027367./345600.)*e14 + (708871./76800.)*e12 + (3311./1280.)*e10 + (5929./768.)*e8,
            0.
        ),
        (
            '-3n-2o',
            -3.*orbital_freq - 2.*spin_freq,
            (-7./55296.)*i14e2 - (5453./983040.)*i12e4 + (931./737280.)*i12e2 - (13034509./123863040.)*i10e6 + (
                        287./8192.)*i10e4 - (49./6144.)*i10e2 + (654201./5242880.)*i8e8 + (289651./655360.)*i8e6 - (
                        861./8192.)*i8e4 + (49./2048.)*i8e2 - (100348117./65536000.)*i6e10 - (6943./8192.)*i6e8 - (
                        2809./3072.)*i6e6 + (60734479./20971520.)*i4e12 + (60222867./26214400.)*i4e10 + (
                        20829./16384.)*i4e8 + (2809./2048.)*i4e6 - (4284333./1468006400.)*i2e14 - (
                        6561./2621440.)*i2e12 - (6561./3276800.)*i2e10 + (9428157./3355443200.)*e16 + (
                        4284333./1468006400.)*e14 + (6561./2621440.)*e12 + (6561./3276800.)*e10,
            2.
        ),
        (
            '-3n-o',
            -3.*orbital_freq - spin_freq,
            (1211./2211840.)*i14e2 + (18983./737280.)*i12e4 - (3241./552960.)*i12e2 + (928873373./1857945600.)*i10e6 - (
                        7749./40960.)*i10e4 + (441./10240.)*i10e2 + (36691./13762560.)*i8e8 - (
                        12674957./5160960.)*i8e6 + (861./1024.)*i8e4 - (49./256.)*i8e2 + (
                        5972707919./786432000.)*i6e10 + (1518251./983040.)*i6e8 + (2427083./368640.)*i6e6 - (
                        861./512.)*i6e4 + (49./128.)*i6e2 - (242874493./15728640.)*i4e12 - (
                        16055203./1310720.)*i4e10 - (6943./1024.)*i4e8 - (2809./384.)*i4e6 + (
                        5540569119./367001600.)*i2e14 + (7590443./655360.)*i2e12 + (15052983./1638400.)*i2e10 + (
                        20829./4096.)*i2e8 + (2809./512.)*i2e6,
            1.
        ),
        (
            '-2o',
            -2.*spin_freq,
            (5461./1702701000.)*i16 - (1./6930.)*i14e2 - (1./20790.)*i14 + (31./9450.)*i12e4 + (31./18900.)*i12e2 + (
                        31./56700.)*i12 - (17./378.)*i10e6 - (17./630.)*i10e4 - (17./1260.)*i10e2 - (17./3780.)*i10 + (
                        3./8.)*i8e8 + (1./4.)*i8e6 + (3./20.)*i8e4 + (3./40.)*i8e2 + (1./40.)*i8 - (7./4.)*i6e10 - (
                        5./4.)*i6e8 - (5./6.)*i6e6 - (1./2.)*i6e4 - (1./4.)*i6e2 - (1./12.)*i6 + (7./2.)*i4e12 + (
                        21./8.)*i4e10 + (15./8.)*i4e8 + (5./4.)*i4e6 + (3./4.)*i4e4 + (3./8.)*i4e2 + (1./8.)*i4,
            2.
        ),
        (
            '4n-2o',
            4.*orbital_freq - 2.*spin_freq,
            (95476063./7664025600.)*i12e4 + (4733837./8709120.)*i10e6 - (3498923./29030400.)*i10e4 + (
                        176956279./23224320.)*i8e8 - (765187./193536.)*i8e6 + (565573./645120.)*i8e4 + (
                        24247471./691200.)*i6e10 - (16173799./414720.)*i6e8 + (8993./432.)*i6e6 - (6647./1440.)*i6e4 + (
                        241167361./2764800.)*i4e12 - (2344199./18432.)*i4e10 + (7512571./55296.)*i4e8 - (
                        21505./288.)*i4e6 + (3179./192.)*i4e4 + (346700573./3628800.)*i2e14 - (
                        8421731./46080.)*i2e12 + (134209./480.)*i2e10 - (83551./288.)*i2e8 + (1955./12.)*i2e6 - (
                        289./8.)*i2e4 + (173370469./4147200.)*e16 - (346700573./3628800.)*e14 + (
                        8421731./46080.)*e12 - (134209./480.)*e10 + (83551./288.)*e8 - (1955./12.)*e6 + (289./8.)*e4,
            2.
        ),
        (
            '-n-o',
            -orbital_freq - spin_freq,
            (34420297./77491814400.)*i14e2 - (9411719./851558400.)*i12e4 - (1074041./212889600.)*i12e2 + (
                        1354163761./8360755200.)*i10e6 + (130757./1433600.)*i10e4 + (133907./3225600.)*i10e2 - (
                        260044601./185794560.)*i8e8 - (4231445./4644864.)*i8e6 - (18397./35840.)*i8e4 - (
                        2083./8960.)*i8e2 + (1962921403./283115520.)*i6e10 + (13001251./2654208.)*i6e8 + (
                        5288671./1658880.)*i6e6 + (4603./2560.)*i6e4 + (517./640.)*i6e2 - (
                        110879314849./6370099200.)*i4e12 - (153346823./11796480.)*i4e10 - (2031247./221184.)*i4e8 - (
                        165245./27648.)*i4e6 - (27./8.)*i4e4 - (3./2.)*i4e2 + (499264587977./29727129600.)*i2e14 + (
                        13860031763./1061683200.)*i2e12 + (3833717./393216.)*i2e10 + (126955./18432.)*i2e8 + (
                        1291./288.)*i2e6 + (81./32.)*i2e4 + (9./8.)*i2e2,
            1.
        ),
        (
            'n-o',
            orbital_freq - spin_freq,
            (5654153./12680478720.)*i14e2 - (84639041./7664025600.)*i12e4 - (9732799./1916006400.)*i12e2 + (
                        451481413./2786918400.)*i10e6 + (2115311./23224320.)*i10e4 + (6971./165888.)*i10e2 - (
                        260104871./185794560.)*i8e8 - (1411001./1548288.)*i8e6 - (165163./322560.)*i8e4 - (
                        19157./80640.)*i8e2 + (3272128931./471859200.)*i6e10 + (65033009./13271040.)*i6e8 + (
                        1764043./552960.)*i6e6 + (8249./4608.)*i6e4 + (967./1152.)*i6e2 - (
                        6161261113./353894400.)*i4e12 - (25565953./1966080.)*i4e10 - (2032717./221184.)*i4e8 - (
                        55145./9216.)*i4e6 - (643./192.)*i4e4 - (77./48.)*i4e2 + (13316818727./792723456.)*i2e14 + (
                        962827./73728.)*i2e12 + (1598197./163840.)*i2e10 + (63551./9216.)*i2e8 + (3449./768.)*i2e6 + (
                        5./2.)*i2e4 + (5./4.)*i2e2,
            1.
        ),
        (
            '3n-o',
            3.*orbital_freq - spin_freq,
            (8988527./14233190400.)*i14e2 + (2424043./72990720.)*i12e4 - (413861./54743040.)*i12e2 + (
                        1253268323./1857945600.)*i10e6 - (1651357./5529600.)*i10e4 + (281939./4147200.)*i10e2 + (
                        18592061./13762560.)*i8e8 - (21684707./5160960.)*i8e6 + (5945./3072.)*i8e4 - (
                        1015./2304.)*i8e2 + (1870083469./157286400.)*i6e10 - (6718523./983040.)*i6e8 + (
                        6426533./368640.)*i6e6 - (65149./7680.)*i6e4 + (11123./5760.)*i6e2 - (
                        14687593./1572864.)*i4e12 - (17307023./655360.)*i4e10 + (170741./8192.)*i4e8 - (
                        132347./3072.)*i4e6 + (1435./64.)*i4e4 - (245./48.)*i4e2 + (924913403./52428800.)*i2e14 + (
                        55810297./13107200.)*i2e12 + (17156499./655360.)*i2e10 - (57471./2048.)*i2e8 + (
                        1549./32.)*i2e6 - (861./32.)*i2e4 + (49./8.)*i2e2,
            1.
        ),
        (
            '-2n-o',
            -2.*orbital_freq - spin_freq,
            (-437./141926400.)*i16 - (173./774144.)*i14e2 + (173./3870720.)*i14 - (234181./15769600.)*i12e4 + (
                        463./193536.)*i12e2 - (463./967680.)*i12 + (7337./57600.)*i10e6 + (85381./716800.)*i10e4 - (
                        9./512.)*i10e2 + (9./2560.)*i10 - (61884479./46448640.)*i8e8 - (1057./1440.)*i8e6 - (
                        11421./17920.)*i8e4 + (5./64.)*i8e2 - (1./64.)*i8 + (108469001./16588800.)*i6e10 + (
                        15401609./3317760.)*i6e8 + (1921./720.)*i6e6 + (2619./1280.)*i6e4 - (5./32.)*i6e2 + (
                        1./32.)*i6 - (853261./51840.)*i4e12 - (423287./34560.)*i4e10 - (59801./6912.)*i4e8 - (
                        21./4.)*i4e6 - (27./8.)*i4e4 + (230112373./14515200.)*i2e14 + (68263727./5529600.)*i2e12 + (
                        105827./11520.)*i2e10 + (14951./2304.)*i2e8 + (63./16.)*i2e6 + (81./32.)*i2e4,
            1.
        ),
        (
            '3n',
            3.*orbital_freq,
            (-7./7920.)*i14e2 - (1271./28800.)*i12e4 + (217./21600.)*i12e2 - (25451359./29030400.)*i10e6 + (
                        697./1920.)*i10e4 - (119./1440.)*i10e2 - (3193577./3440640.)*i8e8 + (6321823./1290240.)*i8e6 - (
                        2583./1280.)*i8e4 + (147./320.)*i8e2 - (112073569./7864320.)*i6e10 + (675737./245760.)*i6e8 - (
                        1539439./92160.)*i6e6 + (861./128.)*i6e4 - (49./32.)*i6e2 + (7028958283./314572800.)*i4e12 + (
                        688791973./26214400.)*i4e10 - (46277./32768.)*i4e8 + (343843./12288.)*i4e6 - (
                        2583./256.)*i4e4 + (147./64.)*i4e2 - (3165427449./104857600.)*i2e14 - (
                        30355211./1310720.)*i2e12 - (6019881./327680.)*i2e10 - (20829./2048.)*i2e8 - (
                        2809./256.)*i2e6 + (145747579353./11744051200.)*e16 + (1055142483./104857600.)*e14 + (
                        30355211./3932160.)*e12 + (2006627./327680.)*e10 + (6943./2048.)*e8 + (2809./768.)*e6,
            0.
        ),
        (
            '-n-2o',
            -orbital_freq - 2.*spin_freq,
            (-337./3041280.)*i14e2 + (285031./103219200.)*i12e4 + (3601./2867200.)*i12e2 - (
                        337187791./8360755200.)*i10e6 - (19549./860160.)*i10e4 - (737./71680.)*i10e2 + (
                        102383837./297271296.)*i8e8 + (41647601./185794560.)*i8e6 + (5179./40960.)*i8e4 + (
                        581./10240.)*i8e2 - (575053217./353894400.)*i6e10 - (7617223./6635520.)*i6e8 - (
                        619673./829440.)*i6e6 - (27./64.)*i6e4 - (3./16.)*i6e2 + (83161756003./25480396800.)*i4e12 + (
                        23002921./9437184.)*i4e10 + (1523515./884736.)*i4e8 + (123941./110592.)*i4e6 + (
                        81./128.)*i4e4 + (9./32.)*i4e2 - (31398887./118908518400.)*i2e14 - (62617./212336640.)*i2e12 - (
                        619./1966080.)*i2e10 - (11./36864.)*i2e8 - (1./4608.)*i2e6 + (147400583./634178764800.)*e16 + (
                        31398887./118908518400.)*e14 + (62617./212336640.)*e12 + (619./1966080.)*e10 + (
                        11./36864.)*e8 + (1./4608.)*e6,
            2.
        ),
        (
            '-2n-2o',
            -2.*orbital_freq - 2.*spin_freq,
            (457./619315200.)*i16 + (5./96768.)*i14e2 - (1./96768.)*i14 + (20527./5734400.)*i12e4 - (
                        19./36864.)*i12e2 + (19./184320.)*i12 - (1129./34560.)*i10e6 - (3999./143360.)*i10e4 + (
                        5./1536.)*i10e2 - (1./1536.)*i10 + (121207039./371589120.)*i8e8 + (2173./11520.)*i8e6 + (
                        2907./20480.)*i8e4 - (5./512.)*i8e2 + (1./512.)*i8 - (396839./259200.)*i6e10 - (
                        112129./103680.)*i6e8 - (21./32.)*i6e6 - (27./64.)*i6e4 + (204810161./66355200.)*i4e12 + (
                        79379./34560.)*i4e10 + (22429./13824.)*i4e8 + (63./64.)*i4e6 + (81./128.)*i4e4 - (
                        2417./1814400.)*i2e14 - (949./691200.)*i2e12 - (7./5760.)*i2e10 - (1./1152.)*i2e8 + (
                        22601./18579456.)*e16 + (2417./1814400.)*e14 + (949./691200.)*e12 + (7./5760.)*e10 + (
                        1./1152.)*e8,
            2.
        ),
        (
            '-o',
            -spin_freq,
            (-8192./638512875.)*i16 + (8192./14189175.)*i14e2 + (8192./42567525.)*i14 - (2048./155925.)*i12e4 - (
                        1024./155925.)*i12e2 - (1024./467775.)*i12 + (512./2835.)*i10e6 + (512./4725.)*i10e4 + (
                        256./4725.)*i10e2 + (256./14175.)*i10 - (32./21.)*i8e8 - (64./63.)*i8e6 - (64./105.)*i8e4 - (
                        32./105.)*i8e2 - (32./315.)*i8 + (112./15.)*i6e10 + (16./3.)*i6e8 + (32./9.)*i6e6 + (
                        32./15.)*i6e4 + (16./15.)*i6e2 + (16./45.)*i6 - (56./3.)*i4e12 - 14*i4e10 - 10*i4e8 - (
                        20./3.)*i4e6 - 4*i4e4 - 2*i4e2 - (2./3.)*i4 + 18*i2e14 + 14*i2e12 + (21./2.)*i2e10 + (
                        15./2.)*i2e8 + 5*i2e6 + 3*i2e4 + (3./2.)*i2e2 + (1./2.)*i2,
            1.
        ),
        (
            '2n-o',
            2.*orbital_freq - spin_freq,
            (-3490169./1046139494400.)*i16 - (8988527./34871316480.)*i14e2 + (8988527./174356582400.)*i14 - (
                        6787897./425779200.)*i12e4 + (59123./19160064.)*i12e2 - (59123./95800320.)*i12 + (
                        15397573./130636800.)*i10e6 + (871763./6451200.)*i10e4 - (40277./1451520.)*i10e2 + (
                        40277./7257600.)*i10 - (63912689./46448640.)*i8e8 - (233029./362880.)*i8e6 - (
                        14291./17920.)*i8e4 + (725./4032.)*i8e2 - (145./4032.)*i8 + (178083313./27648000.)*i6e10 + (
                        3229435./663552.)*i6e8 + (54643./25920.)*i6e6 + (3893./1280.)*i6e4 - (227./288.)*i6e2 + (
                        227./1440.)*i6 - (54777263./3317760.)*i4e12 - (110659./9216.)*i4e10 - (257669./27648.)*i4e8 - (
                        739./216.)*i4e6 - (213./32.)*i4e4 + (25./12.)*i4e2 - (5./12.)*i4 + (
                        1072027567./67737600.)*i2e14 + (13720169./1105920.)*i2e12 + (341669./38400.)*i2e10 + (
                        33595./4608.)*i2e8 + (251./144.)*i2e6 + (207./32.)*i2e4 - (5./2.)*i2e2 + (1./2.)*i2,
            1.
        ),
        (
            'n',
            orbital_freq,
            (-50521./75675600.)*i14e2 + (110287./6652800.)*i12e4 + (4211./554400.)*i12e2 - (353141./1451520.)*i10e6 - (
                        27599./201600.)*i10e4 - (3161./50400.)*i10e2 + (2723191./1290240.)*i8e8 + (
                        29543./21504.)*i8e6 + (6927./8960.)*i8e4 + (793./2240.)*i8e2 - (417490789./39321600.)*i6e10 - (
                        691367./92160.)*i6e8 - (22501./4608.)*i6e6 - (1759./640.)*i6e4 - (201./160.)*i6e2 + (
                        80083677791./2831155200.)*i4e12 + (22151853./1048576.)*i4e10 + (61137./4096.)*i4e8 + (
                        29845./3072.)*i4e6 + (1401./256.)*i4e4 + (159./64.)*i4e2 - (24654653741./734003200.)*i2e14 - (
                        85553819./3276800.)*i2e12 - (3194661./163840.)*i2e10 - (28211./2048.)*i2e8 - (
                        2295./256.)*i2e6 - (81./16.)*i2e4 - (9./4.)*i2e2 + (689312857627./49325015040.)*e16 + (
                        24654653741./2202009600.)*e14 + (85553819./9830400.)*e12 + (1064887./163840.)*e10 + (
                        28211./6144.)*e8 + (765./256.)*e6 + (27./16.)*e4 + (3./4.)*e2,
            0.
        ),
        (
            '3n-2o',
            3.*orbital_freq - 2.*spin_freq,
            (-27269./161740800.)*i14e2 - (13545047./1459814400.)*i12e4 + (2312569./1094860800.)*i12e2 - (
                        7949713./41287680.)*i10e6 + (496387./5529600.)*i10e4 - (84749./4147200.)*i10e2 - (
                        60573773./110100480.)*i8e8 + (18110321./13762560.)*i8e6 - (80237./122880.)*i8e4 + (
                        13699./92160.)*i8e2 - (12131629./3276800.)*i6e10 + (832621./245760.)*i6e8 - (
                        39313./6144.)*i6e6 + (6601./1920.)*i6e4 - (1127./1440.)*i6e2 - (1510073./3276800.)*i4e12 + (
                        13218507./1310720.)*i4e10 - (456169./32768.)*i4e8 + (86193./4096.)*i4e6 - (3157./256.)*i4e4 + (
                        539./192.)*i4e2 - (534226163./209715200.)*i2e14 + (47982879./6553600.)*i2e12 - (
                        5568309./327680.)*i2e10 + (135771./4096.)*i2e8 - (21975./512.)*i2e6 + (861./32.)*i2e4 - (
                        49./8.)*i2e2 - (3973253733./4697620480.)*e16 + (534226163./209715200.)*e14 - (
                        47982879./6553600.)*e12 + (5568309./327680.)*e10 - (135771./4096.)*e8 + (21975./512.)*e6 - (
                        861./32.)*e4 + (49./8.)*e2,
            2.
        ),
        (
            'n-2o',
            orbital_freq - 2.*spin_freq,
            (-9734839./87178291200.)*i14e2 + (16904269./6131220480.)*i12e4 + (1951667./1532805120.)*i12e2 - (
                        112442683./2786918400.)*i10e6 - (2631733./116121600.)*i10e4 - (305867./29030400.)*i10e2 + (
                        512160559./1486356480.)*i8e8 + (13892933./61931520.)*i8e6 + (64927./516096.)*i8e4 + (
                        7649./129024.)*i8e2 - (95879689./58982400.)*i6e10 - (1524797./1327104.)*i6e8 - (
                        206849./276480.)*i6e6 - (2407./5760.)*i6e4 - (293./1440.)*i6e2 + (
                        924591149./283115520.)*i4e12 + (19187029./7864320.)*i4e10 + (1526749./884736.)*i4e8 + (
                        41453./36864.)*i4e6 + (475./768.)*i4e4 + (65./192.)*i4e2 - (165285343./39636172800.)*i2e14 - (
                        277229./58982400.)*i2e12 - (1733./327680.)*i2e10 - (305./36864.)*i2e8 - (13./1536.)*i2e6 + (
                        1./32.)*i2e4 - (1./8.)*i2e2 + (49450862117./13317754060800.)*e16 + (
                        165285343./39636172800.)*e14 + (277229./58982400.)*e12 + (1733./327680.)*e10 + (
                        305./36864.)*e8 + (13./1536.)*e6 - (1./32.)*e4 + (1./8.)*e2,
            2.
        ),
        (
            '2n',
            2.*orbital_freq,
            (5461./1135134000.)*i16 + (1./2772.)*i14e2 - (1./13860.)*i14 + (1219./52800.)*i12e4 - (31./7560.)*i12e2 + (
                        31./37800.)*i12 - (83719./453600.)*i10e6 - (2133./11200.)*i10e4 + (17./504.)*i10e2 - (
                        17./2520.)*i10 + (2625169./1290240.)*i8e8 + (1003./960.)*i8e6 + (4797./4480.)*i8e4 - (
                        3./16.)*i8e2 + (3./80.)*i8 - (4598441./460800.)*i6e10 - (665537./92160.)*i6e8 - (
                        5399./1440.)*i6e6 - (1197./320.)*i6e4 + (5./8.)*i6e2 - (1./8.)*i6 + (
                        197256941./7372800.)*i4e12 + (6103337./307200.)*i4e10 + (58543./4096.)*i4e8 + (
                        1483./192.)*i4e6 + (891./128.)*i4e4 - (15./16.)*i4e2 + (3./16.)*i4 - (
                        25565893./806400.)*i2e14 - (505601./20480.)*i2e12 - (11757./640.)*i2e10 - (1661./128.)*i2e8 - (
                        63./8.)*i2e6 - (81./16.)*i2e4 + (678544541./51609600.)*e16 + (25565893./2419200.)*e14 + (
                        505601./61440.)*e12 + (3919./640.)*e10 + (1661./384.)*e8 + (21./8.)*e6 + (27./16.)*e4,
            0.
        ),
        (
            '2n-2o',
            2.*orbital_freq - 2.*spin_freq,
            (72518377./83691159552000.)*i16 + (27269./396264960.)*i14e2 - (27269./1981324800.)*i14 + (
                        7026553./1703116800.)*i12e4 - (330367./383201280.)*i12e2 + (330367./1916006400.)*i12 - (
                        3670267./130636800.)*i10e6 - (231629./6451200.)*i10e4 + (12107./1451520.)*i10e2 - (
                        12107./7257600.)*i10 + (25557005./74317824.)*i8e8 + (416933./2903040.)*i8e6 + (
                        4549./20480.)*i8e4 - (1957./32256.)*i8e2 + (1957./161280.)*i8 - (10327357./6912000.)*i6e10 - (
                        981971./829440.)*i6e8 - (4871./12960.)*i6e6 - (37./40.)*i6e4 + (23./72.)*i6e2 - (
                        23./360.)*i6 + (20666431./6635520.)*i4e12 + (1994809./921600.)*i4e10 + (
                        220055./110592.)*i4e8 - (37./1728.)*i4e6 + (39./16.)*i4e4 - (55./48.)*i4e2 + (11./48.)*i4 + (
                        1739939./67737600.)*i2e14 - (34471./552960.)*i2e12 + (11041./38400.)*i2e10 - (
                        3697./4608.)*i2e8 + (79./36.)*i2e6 - (63./16.)*i2e4 + (5./2.)*i2e2 - (1./2.)*i2 - (
                        561889./240844800.)*e16 - (1739939./67737600.)*e14 + (34471./552960.)*e12 - (
                        11041./38400.)*e10 + (3697./4608.)*e8 - (79./36.)*e6 + (63./16.)*e4 - (5./2.)*e2 + (1./2.),
            2.
        )
    )

    heating_coeffs = list()
    ztorque_coeffs = list()
    freqs = list()
    modes = list()
    mode_names = list()

    for mode_name, mode, coeff, torque_multi in modes_coeffs:
        sgn = np.sign(mode)
        freq = np.abs(mode)

        heating_coeff = freq * coeff
        ztorque_coeff = coeff * torque_multi * sgn

        # freq_too_low = freq < MODE_ZERO_TOL
        # freq[freq_too_low] = 0.
        # mode[freq_too_low] = 0.
        # heating_coeff[freq_too_low] = 0.
        # ztorque_coeff[freq_too_low] = 0.

        if freq < MODE_ZERO_TOL:
            freq = 0.
            heating_coeff = 0.
            ztorque_coeff = 0.

        freqs.append(freq)
        modes.append(mode)
        mode_names.append(mode_name)
        heating_coeffs.append(heating_coeff)
        ztorque_coeffs.append(ztorque_coeff)

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return mode_names, modes, freqs, heating_coeffs, ztorque_coeffs


@njit
def nsr_modes_20(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray):
    """ Non-Synchronous Rotation Tidal Modes (for heating and torque) for l=2, e^20, I^20

    These should all have the same shape!!

    :param orbital_freq:
    :param spin_freq:
    :param eccentricity:
    :param inclination:
    :return:
    """

    # Conversion from input_heating_nsr_20.dat
    # Eccentricity max order = 20, Inclination max order = 20.
    # Total Number of Unique Modes = 52

    i2 = inclination**2
    i4 = inclination**4
    i6 = inclination**6
    i8 = inclination**8
    i10 = inclination**10
    i12 = inclination**12
    i14 = inclination**14
    i16 = inclination**16
    i18 = inclination**18
    i20 = inclination**20
    e2 = eccentricity**2
    e4 = eccentricity**4
    e6 = eccentricity**6
    e8 = eccentricity**8
    e10 = eccentricity**10
    e12 = eccentricity**12
    e14 = eccentricity**14
    e16 = eccentricity**16
    e18 = eccentricity**18
    e20 = eccentricity**20
    i2e2 = i2*e2
    i2e4 = i2*e4
    i2e6 = i2*e6
    i2e8 = i2*e8
    i2e10 = i2*e10
    i2e12 = i2*e12
    i2e14 = i2*e14
    i2e16 = i2*e16
    i2e18 = i2*e18
    i4e2 = i4*e2
    i4e4 = i4*e4
    i4e6 = i4*e6
    i4e8 = i4*e8
    i4e10 = i4*e10
    i4e12 = i4*e12
    i4e14 = i4*e14
    i4e16 = i4*e16
    i6e2 = i6*e2
    i6e4 = i6*e4
    i6e6 = i6*e6
    i6e8 = i6*e8
    i6e10 = i6*e10
    i6e12 = i6*e12
    i6e14 = i6*e14
    i8e2 = i8*e2
    i8e4 = i8*e4
    i8e6 = i8*e6
    i8e8 = i8*e8
    i8e10 = i8*e10
    i8e12 = i8*e12
    i10e2 = i10*e2
    i10e4 = i10*e4
    i10e6 = i10*e6
    i10e8 = i10*e8
    i10e10 = i10*e10
    i12e2 = i12*e2
    i12e4 = i12*e4
    i12e6 = i12*e6
    i12e8 = i12*e8
    i14e2 = i14*e2
    i14e4 = i14*e4
    i14e6 = i14*e6
    i16e2 = i16*e2
    i16e4 = i16*e4
    i18e2 = i18*e2

    modes_coeffs = (
        (
            '11n-o',
            11.*orbital_freq - spin_freq,
            (6648821549377771726369./69039237051187200.)*i2e18,
            1.
        ),
        (
            '12n-2o',
            12.*orbital_freq - 2.*spin_freq,
            (3816001995797209./16056320000.)*e20,
            2.
        ),
        (
            '10n',
            10.*orbital_freq,
            (402063787225./28311552.)*i4e16 + (58301303109841./112368549888.)*e20,
            0.
        ),
        (
            '-9n-o',
            -9.*orbital_freq - spin_freq,
            (147483366698529./164416716800.)*i6e14 + (4143587919679849./10522669875200.)*i2e18,
            1.
        ),
        (
            '11n-2o',
            11.*orbital_freq - 2.*spin_freq,
            (-6648821549377771726369./69039237051187200.)*i2e18 - (
                        987278781529197450645517./1380784741023744000.)*e20 + (
                        6648821549377771726369./69039237051187200.)*e18,
            2.
        ),
        (
            '10n-o',
            10.*orbital_freq - spin_freq,
            (-2010318936125./63700992.)*i4e16 - (176187983600875./668860416.)*i2e18 + (402063787225./10616832.)*i2e16,
            1.
        ),
        (
            '-8n-2o',
            -8.*orbital_freq - 2.*spin_freq,
            (5383010161./265420800.)*i8e12 + (31801945561./642252800.)*i4e16 + (8388608./200930625.)*e20,
            2.
        ),
        (
            '10n-2o',
            10.*orbital_freq - 2.*spin_freq,
            (4422701659475./254803968.)*i4e16 + (176187983600875./668860416.)*i2e18 - (
                        402063787225./10616832.)*i2e16 + (285521789807747375./337105649664.)*e20 - (
                        176187983600875./668860416.)*e18 + (402063787225./10616832.)*e16,
            2.
        ),
        (
            '9n',
            9.*orbital_freq,
            (-147483366698529./41104179200.)*i6e14 - (9194245551902151./263066746880.)*i4e16 + (
                        442450100095587./82208358400.)*i4e14 - (4143587919679849./5261334937600.)*i2e18 - (
                        9482058568573459./15032385536000.)*e20 + (4143587919679849./15784004812800.)*e18,
            0.
        ),
        (
            '9n-o',
            9.*orbital_freq - spin_freq,
            (3719858248951787./822083584000.)*i6e14 + (1021582839100239./13153337344.)*i4e16 - (
                        49161122232843./4110417920.)*i4e14 + (1457889300753667049./5261334937600.)*i2e18 - (
                        3064748517300717./32883343360.)*i2e16 + (147483366698529./10276044800.)*i2e14,
            1.
        ),
        (
            '-8n-o',
            -8.*orbital_freq - spin_freq,
            (-5383010161./33177600.)*i8e12 - (56914314263./29030400.)*i6e14 + (5383010161./16588800.)*i6e12 - (
                        31801945561./120422400.)*i4e16 - (1606626956771./4335206400.)*i2e18 + (
                        31801945561./160563200.)*i2e16,
            1.
        ),
        (
            '-7n-2o',
            -7.*orbital_freq - 2.*spin_freq,
            (-52142352409./22649241600.)*i10e10 - (140254152605./3623878656.)*i8e12 + (
                        52142352409./7549747200.)*i8e10 - (186702632281./11324620800.)*i6e14 - (
                        1328179895713./40265318400.)*i4e16 + (186702632281./7549747200.)*i4e14 - (
                        33232930569601./1408964021452800.)*i2e18 - (33232930569601./28179280429056000.)*e20 + (
                        33232930569601./1408964021452800.)*e18,
            2.
        ),
        (
            '9n-2o',
            9.*orbital_freq - 2.*spin_freq,
            (-376901937118463./205520896000.)*i6e14 - (11237411230102629./263066746880.)*i4e16 + (
                        540772344561273./82208358400.)*i4e14 - (415947859083950607./1503238553600.)*i2e18 + (
                        3064748517300717./32883343360.)*i2e16 - (147483366698529./10276044800.)*i2e14 - (
                        107205587887490347401./210453397504000.)*e20 + (415947859083950607./1503238553600.)*e18 - (
                        3064748517300717./32883343360.)*e16 + (147483366698529./10276044800.)*e14,
            2.
        ),
        (
            '8n-o',
            8.*orbital_freq - spin_freq,
            (-156107294669./418037760.)*i8e12 - (12919549337701./1306368000.)*i6e14 + (
                        1221943306547./746496000.)*i6e12 - (4741268850403./66355200.)*i4e16 + (
                        56914314263./2177280.)*i4e14 - (5383010161./1244160.)*i4e12 - (
                        50221734739989587./351151718400.)*i2e18 + (371200286353507./4335206400.)*i2e16 - (
                        56914314263./1814400.)*i2e14 + (5383010161./1036800.)*i2e12,
            1.
        ),
        (
            '8n',
            8.*orbital_freq,
            (5383010161./13824000.)*i8e12 + (56914314263./7257600.)*i6e14 - (5383010161./4147200.)*i6e12 + (
                        93825684332719./2890137600.)*i4e16 - (56914314263./4838400.)*i4e14 + (
                        5383010161./2764800.)*i4e12 + (1606626956771./2167603200.)*i2e18 - (
                        31801945561./80281600.)*i2e16 + (583590180249631./1755758592000.)*e20 - (
                        1606626956771./6502809600.)*e18 + (31801945561./240844800.)*e16,
            0.
        ),
        (
            '-7n-o',
            -7.*orbital_freq - spin_freq,
            (52142352409./4194304000.)*i10e10 + (140254152605./452984832.)*i8e12 - (52142352409./943718400.)*i8e10 + (
                        2182111894332367./1358954496000.)*i6e14 - (140254152605./226492416.)*i6e12 + (
                        52142352409./471859200.)*i6e10 + (1328179895713./7549747200.)*i4e16 - (
                        186702632281./1415577600.)*i4e14 + (121515502749179129./704482010726400.)*i2e18 - (
                        1328179895713./10066329600.)*i2e16 + (186702632281./1887436800.)*i2e14,
            1.
        ),
        (
            '8n-2o',
            8.*orbital_freq - 2.*spin_freq,
            (10534550885077./83607552000.)*i8e12 + (1309029228049./326592000.)*i6e14 - (
                        123809233703./186624000.)*i6e12 + (2039454943618921./52022476800.)*i4e16 - (
                        626057456893./43545600.)*i4e14 + (59213111771./24883200.)*i4e12 + (
                        391340609035087./2743372800.)*i2e18 - (578658802849./6773760.)*i2e16 + (
                        56914314263./1814400.)*i2e14 - (5383010161./1036800.)*i2e12 + (
                        74152056390168773./438939648000.)*e20 - (391340609035087./2743372800.)*e18 + (
                        578658802849./6773760.)*e16 - (56914314263./1814400.)*e14 + (5383010161./1036800.)*e12,
            2.
        ),
        (
            '7n-o',
            7.*orbital_freq - spin_freq,
            (300019646853899./15288238080000.)*i10e10 + (2905264589675./4076863488.)*i8e12 - (
                        216018317123./1698693120.)*i8e10 + (95588340385394921./12230590464000.)*i6e14 - (
                        6367538528267./2038431744.)*i6e12 + (11836313996843./21233664000.)*i6e10 + (
                        50276885204422843./1630745395200.)*i4e16 - (2099962736128727./101921587200.)*i4e14 + (
                        701270763025./84934656.)*i4e12 - (52142352409./35389440.)*i4e10 + (
                        1219852130257107097./31310311587840.)*i2e18 - (5016930263287009./135895449600.)*i2e16 + (
                        104746088252857./4246732800.)*i2e14 - (140254152605./14155776.)*i2e12 + (
                        52142352409./29491200.)*i2e10,
            1.
        ),
        (
            '-6n-2o',
            -6.*orbital_freq - 2.*spin_freq,
            (5397691./47185920.)*i12e8 + (2456597./655360.)*i10e10 - (284089./393216.)*i10e8 + (
                        3620308013./131072000.)*i8e12 - (7369791./655360.)*i8e10 + (284089./131072.)*i8e8 + (
                        2308743./358400.)*i6e14 - (10029889./1228800.)*i6e12 + (1235095623./80281600.)*i4e16 - (
                        6926229./716800.)*i4e14 + (10029889./819200.)*i4e12 - (177147./40140800.)*i2e18 - (
                        531441./40140800.)*i2e16 + (18128043./1605632000.)*e20 + (177147./40140800.)*e18 + (
                        531441./40140800.)*e16,
            2.
        ),
        (
            '7n',
            7.*orbital_freq,
            (-126631427279./5308416000.)*i10e10 - (28050830521./37748736.)*i8e12 + (52142352409./393216000.)*i8e10 - (
                        707704620843857./113246208000.)*i6e14 + (140254152605./56623104.)*i6e12 - (
                        52142352409./117964800.)*i6e10 - (2041039146624199./144955146240.)*i4e16 + (
                        427012566199511./45298483200.)*i4e14 - (140254152605./37748736.)*i4e12 + (
                        52142352409./78643200.)*i4e10 - (999990833612299./2899102924800.)*i2e18 + (
                        1328179895713./5033164800.)*i2e16 - (186702632281./943718400.)*i2e14 - (
                        105829510207034821./1565515579392000.)*e20 + (999990833612299./8697308774400.)*e18 - (
                        1328179895713./15099494400.)*e16 + (186702632281./2831155200.)*e14,
            0.
        ),
        (
            '-6n-o',
            -6.*orbital_freq - spin_freq,
            (-131533207./247726080.)*i12e8 - (66328119./3276800.)*i10e10 + (2556801./655360.)*i10e8 - (
                        218129754931./1032192000.)*i8e12 + (7369791./81920.)*i8e10 - (284089./16384.)*i8e8 - (
                        16066659729./28672000.)*i6e14 + (32261815669./73728000.)*i6e12 - (7369791./40960.)*i6e10 + (
                        284089./8192.)*i6e8 - (188156121./2293760.)*i4e16 + (2308743./44800.)*i4e14 - (
                        10029889./153600.)*i4e12 - (188201579./16056320.)*i2e18 + (4939496757./80281600.)*i2e16 - (
                        6926229./179200.)*i2e14 + (10029889./204800.)*i2e12,
            1.
        ),
        (
            '7n-2o',
            7.*orbital_freq - 2.*spin_freq,
            (-90183922945109./15288238080000.)*i10e10 - (7842210761371./32614907904.)*i8e12 + (
                        14577511952059./339738624000.)*i8e10 - (9648402385096547./3057647616000.)*i6e14 + (
                        645169101983./509607936.)*i6e12 - (1199274105407./5308416000.)*i6e10 - (
                        110193161506392943./6522981580800.)*i4e16 + (4600426264673063./407686348800.)*i4e14 - (
                        1542795678655./339738624.)*i4e12 + (573565876499./707788800.)*i4e10 - (
                        1518065224694500853./39137889484800.)*i2e18 + (9997999669389767./271790899200.)*i2e16 - (
                        417304029320899./16986931200.)*i2e14 + (140254152605./14155776.)*i2e12 - (
                        52142352409./29491200.)*i2e10 - (373553969127047054087./11741366845440000.)*e20 + (
                        1518065224694500853./39137889484800.)*e18 - (9997999669389767./271790899200.)*e16 + (
                        417304029320899./16986931200.)*e14 - (140254152605./14155776.)*e12 + (
                        52142352409./29491200.)*e10,
            2.
        ),
        (
            '6n-o',
            6.*orbital_freq - spin_freq,
            (-16796193947./24524881920.)*i12e8 - (98944357369./3096576000.)*i10e10 + (
                        11442252653./1857945600.)*i10e8 - (488676583261./1032192000.)*i8e12 + (
                        71241313./344064.)*i8e10 - (41192905./1032192.)*i8e8 - (46716086489./17203200.)*i6e14 + (
                        30471642527./14745600.)*i6e12 - (557647519./614400.)*i6e10 + (64488203./368640.)*i6e8 - (
                        426537070737./64225280.)*i4e16 + (5129809483./716800.)*i4e14 - (3339471121./614400.)*i4e12 + (
                        2456597./1024.)*i4e10 - (1420445./3072.)*i4e8 - (911835584033./160563200.)*i2e18 + (
                        1273685091561./160563200.)*i2e16 - (613914843./71680.)*i2e14 + (265954103./40960.)*i2e12 - (
                        7369791./2560.)*i2e10 + (284089./512.)*i2e8,
            1.
        ),
        (
            '6n',
            6.*orbital_freq,
            (8806759./9676800.)*i12e8 + (41762149./1075200.)*i10e10 - (4829513./645120.)*i10e8 + (
                        128587026311./258048000.)*i8e12 - (22109373./102400.)*i8e10 + (852267./20480.)*i8e8 + (
                        15580284537./7168000.)*i6e14 - (30677093207./18432000.)*i6e12 + (7369791./10240.)*i6e10 - (
                        284089./2048.)*i6e8 + (158505203589./51380224.)*i4e16 - (9407274543./2867200.)*i4e14 + (
                        6199609931./2457600.)*i4e12 - (22109373./20480.)*i4e10 + (852267./4096.)*i4e8 + (
                        941362189./40140800.)*i2e18 - (197537355./1605632.)*i2e16 + (6926229./89600.)*i2e14 - (
                        10029889./102400.)*i2e12 + (168060374177./6422528000.)*e20 - (941362189./120422400.)*e18 + (
                        65845785./1605632.)*e16 - (2308743./89600.)*e14 + (10029889./307200.)*e12,
            0.
        ),
        (
            '-5n-2o',
            -5.*orbital_freq - 2.*spin_freq,
            (-714025./222953472.)*i14e6 - (104437775./679477248.)*i12e8 + (2713295./84934656.)*i12e6 - (
                        22832994493./10569646080.)*i10e10 + (27483625./28311552.)*i10e8 - (714025./3538944.)*i10e6 - (
                        378633873203./54358179840.)*i8e12 + (3539684443./503316480.)*i8e10 - (
                        27483625./9437184.)*i8e8 + (714025./1179648.)*i8e6 - (5113269170029./1198597865472.)*i6e14 + (
                        982439./1048576.)*i6e12 - (1047843./262144.)*i6e10 + (32339308979875./12785043898368.)*i4e16 + (
                        5114831670029./799065243648.)*i4e14 - (2947317./2097152.)*i4e12 + (3143529./524288.)*i4e10 - (
                        323974609375./43834436222976.)*i2e18 - (2685546875./532710162432.)*i2e16 - (
                        244140625./33294385152.)*i2e14 + (21589111328125./3682092642729984.)*e20 + (
                        323974609375./43834436222976.)*e18 + (2685546875./532710162432.)*e16 + (
                        244140625./33294385152.)*e14,
            2.
        ),
        (
            '-5n-o',
            -5.*orbital_freq - spin_freq,
            (24705265./1783627776.)*i14e6 + (2544983675./3567255552.)*i12e8 - (66118715./445906944.)*i12e6 + (
                        66745437623./5872025600.)*i10e10 - (5496725./1048576.)*i10e8 + (142805./131072.)*i10e6 + (
                        2597809820069./47563407360.)*i8e12 - (22698870589./440401920.)*i8e10 + (
                        27483625./1179648.)*i8e8 - (714025./147456.)*i8e6 + (1201932519871./11705057280.)*i6e14 - (
                        376936218611./3397386240.)*i6e12 + (3472622491./31457280.)*i6e10 - (27483625./589824.)*i6e8 + (
                        714025./73728.)*i6e6 - (43093118353375./3196260974592.)*i4e16 - (
                        6817415533997./199766310912.)*i4e14 + (982439./131072.)*i4e12 - (1047843./32768.)*i4e10 + (
                        7139291827499245./306841053560832.)*i2e18 + (2693823437125./266355081216.)*i2e16 + (
                        213067123621./8323596288.)*i2e14 - (2947317./524288.)*i2e12 + (3143529./131072.)*i2e10,
            1.
        ),
        (
            '6n-2o',
            6.*orbital_freq - 2.*spin_freq,
            (93853630663./490497638400.)*i12e8 + (29742019879./3096576000.)*i10e10 - (3439465523./1857945600.)*i10e8 + (
                        87439097251./550502400.)*i8e12 - (4807560329./68812800.)*i8e10 + (555962173./41287680.)*i8e8 + (
                        23565724717./21504000.)*i6e14 - (851525863./1024000.)*i6e12 + (56501731./153600.)*i6e10 - (
                        6534047./92160.)*i6e8 + (4653717222807./1284505600.)*i4e16 - (449280721./114688.)*i4e14 + (
                        242955437./81920.)*i4e12 - (27022567./20480.)*i4e10 + (3124979./12288.)*i4e8 + (
                        25998653133./4587520.)*i2e18 - (180544031973./22937600.)*i2e16 + (3055721757./358400.)*i2e14 - (
                        659870313./102400.)*i2e12 + (7369791./2560.)*i2e10 - (284089./512.)*i2e8 + (
                        21810011108497./6422528000.)*e20 - (25998653133./4587520.)*e18 + (
                        180544031973./22937600.)*e16 - (3055721757./358400.)*e14 + (659870313./102400.)*e12 - (
                        7369791./2560.)*e10 + (284089./512.)*e8,
            2.
        ),
        (
            '5n-o',
            5.*orbital_freq - spin_freq,
            (19747793819./1236054048768.)*i14e6 + (324982872175./353158299648.)*i12e8 - (
                        8443060015./44144787456.)*i12e6 + (24887762328139./1426902220800.)*i10e10 - (
                        44278318565./5350883328.)*i10e8 + (1150351397./668860416.)*i10e6 + (
                        53175396894371./428070666240.)*i8e12 - (63578891293./566231040.)*i8e10 + (
                        3985125625./74317824.)*i8e8 - (103533625./9289728.)*i8e6 + (
                        332893400993369./749123665920.)*i6e14 - (16618547307749./30576476160.)*i6e12 + (
                        138128620669./283115520.)*i6e10 - (1247756575./5308416.)*i6e8 + (32416735./663552.)*i6e6 + (
                        104322729928625./152202903552.)*i4e16 - (33120820080571./28538044416.)*i4e14 + (
                        365264843291./254803968.)*i4e12 - (3011571571./2359296.)*i4e10 + (137418125./221184.)*i4e8 - (
                        3570125./27648.)*i4e6 + (52764829492826545./115065395085312.)*i2e18 - (
                        10509201011825./12683575296.)*i2e16 + (1637783339303./1189085184.)*i2e14 - (
                        18227432263./10616832.)*i2e12 + (298327981./196608.)*i2e10 - (27483625./36864.)*i2e8 + (
                        714025./4608.)*i2e6,
            1.
        ),
        (
            '5n',
            5.*orbital_freq,
            (-142805./6386688.)*i14e6 - (34079695./27869184.)*i12e8 + (885391./3483648.)*i12e6 - (
                        53137321943./2477260800.)*i10e10 + (93444325./9289728.)*i10e8 - (2427685./1161216.)*i10e6 - (
                        515522960797./3963617280.)*i8e12 + (4380222557./36700160.)*i8e10 - (5496725./98304.)*i8e8 + (
                        142805./12288.)*i8e6 - (121792981737923./332943851520.)*i6e14 + (
                        368554049063./849346560.)*i6e12 - (3141504103./7864320.)*i6e10 + (27483625./147456.)*i6e8 - (
                        714025./18432.)*i6e6 - (207827804192825./710280216576.)*i4e16 + (
                        24964480332551./44392513536.)*i4e14 - (74050340731./113246208.)*i4e12 + (
                        641713211./1048576.)*i4e10 - (27483625./98304.)*i4e8 + (714025./12288.)*i4e6 - (
                        2447539096445./52613349376.)*i2e18 - (28491858875./1409286144.)*i2e16 - (
                        13524196093./264241152.)*i2e14 + (2947317./262144.)*i2e12 - (3143529./65536.)*i2e10 + (
                        321698310085813./21917218111488.)*e20 + (2447539096445./157840048128.)*e18 + (
                        28491858875./4227858432.)*e16 + (13524196093./792723456.)*e14 - (982439./262144.)*e12 + (
                        1047843./65536.)*e10,
            0.
        ),
        (
            '-4n-o',
            -4.*orbital_freq - spin_freq,
            (-126293./567705600.)*i16e4 - (67643./4644864.)*i14e6 + (49997./15482880.)*i14e4 - (
                        228744469./696729600.)*i12e8 + (181033./1161216.)*i12e6 - (133807./3870720.)*i12e4 - (
                        50479829./27648000.)*i10e10 + (40776779./16588800.)*i10e8 - (1173./1024.)*i10e6 + (
                        2601./10240.)*i10e4 - (178180572349./20901888000.)*i8e12 + (610501./76800.)*i8e10 - (
                        526171./46080.)*i8e8 + (1955./384.)*i8e6 - (289./256.)*i8e4 + (352576639./104509440.)*i6e14 + (
                        6351002783./298598400.)*i6e12 - (565093./38400.)*i6e10 + (607483./23040.)*i6e8 - (
                        1955./192.)*i6e6 + (289./128.)*i6e4 - (224884235./9289728.)*i4e16 - (
                        27251423./1555200.)*i4e14 - (57428791./3110400.)*i4e12 - (3311./640.)*i4e10 - (
                        5929./384.)*i4e8 + (527512887821./25082265600.)*i2e18 + (5622671123./309657600.)*i2e16 + (
                        5450899./414720.)*i2e14 + (11486987./829440.)*i2e12 + (9933./2560.)*i2e10 + (5929./512.)*i2e8,
            1.
        ),
        (
            '-4n-2o',
            -4.*orbital_freq - 2.*spin_freq,
            (132073./2477260800.)*i16e4 + (1955./580608.)*i14e6 - (289./387072.)*i14e4 + (3205931./44236800.)*i12e8 - (
                        7429./221184.)*i12e6 + (5491./737280.)*i12e4 + (67413./204800.)*i10e10 - (
                        177649./368640.)*i10e8 + (1955./9216.)*i10e6 - (289./6144.)*i10e4 + (
                        47029973497./33443020800.)*i8e12 - (551849./614400.)*i8e10 + (631199./368640.)*i8e8 - (
                        1955./3072.)*i8e6 + (289./2048.)*i8e4 - (408788753./186624000.)*i6e14 - (
                        861466681./373248000.)*i6e12 - (3311./5120.)*i6e10 - (5929./3072.)*i6e8 + (
                        2410254527./530841600.)*i4e16 + (16356793./4976640.)*i4e14 + (34469153./9953280.)*i4e12 + (
                        9933./10240.)*i4e10 + (5929./2048.)*i4e8 - (68./15309.)*i2e18 - (23./4725.)*i2e16 - (
                        8./2025.)*i2e14 - (8./2025.)*i2e12 + (1733531./428652000.)*e20 + (68./15309.)*e18 + (
                        23./4725.)*e16 + (8./2025.)*e14 + (8./2025.)*e12,
            2.
        ),
        (
            '5n-2o',
            5.*orbital_freq - 2.*spin_freq,
            (-59909993./14046068736.)*i14e6 - (363187309615./1412633198592.)*i12e8 + (
                        9435611887./176579149824.)*i12e6 - (1483470394961./285380444160.)*i10e10 + (
                        13309769915./5350883328.)*i10e8 - (345788027./668860416.)*i10e6 - (
                        143179710059539./3424565329920.)*i8e12 + (1187224185659./31708938240.)*i8e10 - (
                        10757090825./594542592.)*i8e8 + (279469385./74317824.)*i8e6 - (
                        30305618313947./171228266496.)*i6e14 + (335718979387./1528823808.)*i6e12 - (
                        2757820247./14155776.)*i6e10 + (126424675./1327104.)*i6e8 - (3284515./165888.)*i6e6 - (
                        232484555909525./608811614208.)*i4e16 + (71453878105147./114152177664.)*i4e14 - (
                        800813356187./1019215872.)*i4e12 + (6516062647./9437184.)*i4e10 - (302319875./884736.)*i4e8 + (
                        7854275./110592.)*i4e6 - (2044426346565875./4696546738176.)*i2e18 + (
                        21274828753525./25367150592.)*i2e16 - (6429415592375./4756340736.)*i2e14 + (
                        72670996375./42467328.)*i2e12 - (587225375./393216.)*i2e10 + (27483625./36864.)*i2e8 - (
                        714025./4608.)*i2e6 - (180939012603859375./920523160682496.)*e20 + (
                        2044426346565875./4696546738176.)*e18 - (21274828753525./25367150592.)*e16 + (
                        6429415592375./4756340736.)*e14 - (72670996375./42467328.)*e12 + (587225375./393216.)*e10 - (
                        27483625./36864.)*e8 + (714025./4608.)*e6,
            2.
        ),
        (
            '4n-o',
            4.*orbital_freq - spin_freq,
            (-1008658841./4184557977600.)*i16e4 - (3514514057./209227898880.)*i14e6 + (
                        2597684303./697426329600.)*i14e4 - (28195995361./68976230400.)*i12e8 + (
                        23117093./114960384.)*i12e6 - (17086547./383201280.)*i12e4 - (1032284497./348364800.)*i10e10 + (
                        3802316939./1045094400.)*i10e8 - (15748307./8709120.)*i10e6 + (11640053./29030400.)*i10e4 - (
                        4118001859./258048000.)*i8e12 + (93487253./4838400.)*i8e10 - (67404683./2903040.)*i8e8 + (
                        283475./24192.)*i8e6 - (41905./16128.)*i8e4 - (54287697463./2612736000.)*i6e14 + (
                        1243549381./18432000.)*i6e12 - (1180475./13824.)*i6e10 + (20673629./207360.)*i6e8 - (
                        88757./1728.)*i6e6 + (65603./5760.)*i6e4 - (8227071409./139345920.)*i4e16 + (
                        1352054623./21772800.)*i4e14 - (236062631./1382400.)*i4e12 + (1312291./5760.)*i4e10 - (
                        888871./3456.)*i4e8 + (9775./72.)*i4e6 - (1445./48.)*i4e4 + (
                        286613322889./58525286400.)*i2e18 + (55698476441./928972800.)*i2e16 - (
                        1196078171./14515200.)*i2e14 + (90597149./460800.)*i2e12 - (423509./1536.)*i2e10 + (
                        1390177./4608.)*i2e8 - (1955./12.)*i2e6 + (289./8.)*i2e4,
            1.
        ),
        (
            '4n',
            4.*orbital_freq,
            (1578229./4540536000.)*i16e4 + (391./16632.)*i14e6 - (289./55440.)*i14e4 + (4005557./7257600.)*i12e8 - (
                        12121./45360.)*i12e6 + (8959./151200.)*i12e4 + (43084901./12096000.)*i10e10 - (
                        32966741./7257600.)*i10e8 + (6647./3024.)*i10e6 - (4913./10080.)*i10e4 + (
                        3473775227./193536000.)*i8e12 - (1519219./76800.)*i8e10 + (1166083./46080.)*i8e8 - (
                        391./32.)*i8e6 + (867./320.)*i8e4 + (43440269./4536000.)*i6e14 - (
                        2520154657./41472000.)*i6e12 + (840647./12800.)*i6e10 - (217949./2560.)*i6e8 + (
                        1955./48.)*i6e6 - (289./32.)*i6e4 + (973378663./17694720.)*i4e16 - (222487./30240.)*i4e14 + (
                        544861901./5529600.)*i4e12 - (493793./5120.)*i4e10 + (411281./3072.)*i4e8 - (1955./32.)*i4e6 + (
                        867./64.)*i4e4 - (2170376447./51609600.)*i2e18 - (124914751./3440640.)*i2e16 - (
                        3027367./115200.)*i2e14 - (708871./25600.)*i2e12 - (9933./1280.)*i2e10 - (5929./256.)*i2e8 + (
                        4872964014199./292626432000.)*e20 + (2170376447./154828800.)*e18 + (
                        124914751./10321920.)*e16 + (3027367./345600.)*e14 + (708871./76800.)*e12 + (
                        3311./1280.)*e10 + (5929./768.)*e8,
            0.
        ),
        (
            '-3n-o',
            -3.*orbital_freq - spin_freq,
            (2743907./1366386278400.)*i18e2 + (17917./108134400.)*i16e4 - (3059./81100800.)*i16e2 + (
                        265500849863./44635285094400.)*i14e6 - (7093./2949120.)*i14e4 + (1211./2211840.)*i14e2 + (
                        441708911./46714060800.)*i12e8 - (7981790359./122624409600.)*i12e6 + (18983./737280.)*i12e4 - (
                        3241./552960.)*i12e2 + (1788735836969./3963617280000.)*i10e10 - (
                        244699459./4954521600.)*i10e8 + (928873373./1857945600.)*i10e6 - (7749./40960.)*i10e4 + (
                        441./10240.)*i10e2 - (8018055157./3774873600.)*i8e12 - (5279234101./2202009600.)*i8e10 + (
                        36691./13762560.)*i8e8 - (12674957./5160960.)*i8e6 + (861./1024.)*i8e4 - (49./256.)*i8e2 + (
                        159918691141./14680064000.)*i6e14 + (14679665557./1887436800.)*i6e12 + (
                        5972707919./786432000.)*i6e10 + (1518251./983040.)*i6e8 + (2427083./368640.)*i6e6 - (
                        861./512.)*i6e4 + (49./128.)*i6e2 - (1166090629989./46976204800.)*i4e16 - (
                        59095119603./2936012800.)*i4e14 - (242874493./15728640.)*i4e12 - (16055203./1310720.)*i4e10 - (
                        6943./1024.)*i4e8 - (2809./384.)*i4e6 + (29595673932843./1315333734400.)*i2e18 + (
                        218654367579./11744051200.)*i2e16 + (5540569119./367001600.)*i2e14 + (
                        7590443./655360.)*i2e12 + (15052983./1638400.)*i2e10 + (20829./4096.)*i2e8 + (2809./512.)*i2e6,
            1.
        ),
        (
            '-3n-2o',
            -3.*orbital_freq - 2.*spin_freq,
            (-3437./7007109120.)*i18e2 - (18737./471859200.)*i16e4 + (3199./353894400.)*i16e2 - (
                        1927729./1362493440.)*i14e6 + (41./73728.)*i14e4 - (7./55296.)*i14e2 - (
                        50464553./39636172800.)*i12e8 + (220689271./14863564800.)*i12e6 - (5453./983040.)*i12e4 + (
                        931./737280.)*i12e2 - (138402957053./1321205760000.)*i10e10 - (852013./330301440.)*i10e8 - (
                        13034509./123863040.)*i10e6 + (287./8192.)*i10e4 - (49./6144.)*i10e2 + (
                        517152035./939524096.)*i8e12 + (15434868243./29360128000.)*i8e10 + (654201./5242880.)*i8e8 + (
                        289651./655360.)*i8e6 - (861./8192.)*i8e4 + (49./2048.)*i8e2 - (
                        73870922661./29360128000.)*i6e14 - (303602411./157286400.)*i6e12 - (
                        100348117./65536000.)*i6e10 - (6943./8192.)*i6e8 - (2809./3072.)*i6e6 + (
                        874727465481./187904819200.)*i4e16 + (44331693507./11744051200.)*i4e14 + (
                        60734479./20971520.)*i4e12 + (60222867./26214400.)*i4e10 + (20829./16384.)*i4e8 + (
                        2809./2048.)*i4e6 - (270241029./105226698752.)*i2e18 - (9428157./3355443200.)*i2e16 - (
                        4284333./1468006400.)*i2e14 - (6561./2621440.)*i2e12 - (6561./3276800.)*i2e10 + (
                        591397458471./263066746880000.)*e20 + (270241029./105226698752.)*e18 + (
                        9428157./3355443200.)*e16 + (4284333./1468006400.)*e14 + (6561./2621440.)*e12 + (
                        6561./3276800.)*e10,
            2.
        ),
        (
            '-2o',
            -2.*spin_freq,
            (73./10337827500.)*i20 - (257./510810300.)*i18e2 - (257./1532430900.)*i18 + (5461./283783500.)*i16e4 + (
                        5461./567567000.)*i16e2 + (5461./1702701000.)*i16 - (1./2079.)*i14e6 - (1./3465.)*i14e4 - (
                        1./6930.)*i14e2 - (1./20790.)*i14 + (31./3780.)*i12e8 + (31./5670.)*i12e6 + (
                        31./9450.)*i12e4 + (31./18900.)*i12e2 + (31./56700.)*i12 - (17./180.)*i10e10 - (
                        17./252.)*i10e8 - (17./378.)*i10e6 - (17./630.)*i10e4 - (17./1260.)*i10e2 - (17./3780.)*i10 + (
                        7./10.)*i8e12 + (21./40.)*i8e10 + (3./8.)*i8e8 + (1./4.)*i8e6 + (3./20.)*i8e4 + (
                        3./40.)*i8e2 + (1./40.)*i8 - 3*i6e14 - (7./3.)*i6e12 - (7./4.)*i6e10 - (5./4.)*i6e8 - (
                        5./6.)*i6e6 - (1./2.)*i6e4 - (1./4.)*i6e2 - (1./12.)*i6 + (45./8.)*i4e16 + (9./2.)*i4e14 + (
                        7./2.)*i4e12 + (21./8.)*i4e10 + (15./8.)*i4e8 + (5./4.)*i4e6 + (3./4.)*i4e4 + (3./8.)*i4e2 + (
                        1./8.)*i4,
            2.
        ),
        (
            '4n-2o',
            4.*orbital_freq - 2.*spin_freq,
            (20957810953./334764638208000.)*i16e4 + (10662179./2377589760.)*i14e6 - (7880741./7925299200.)*i14e4 + (
                        31096144609./275904921600.)*i12e8 - (129173497./2299207680.)*i12e6 + (
                        95476063./7664025600.)*i12e4 + (1564078403./1741824000.)*i10e10 - (
                        1120408397./1045094400.)*i10e8 + (4733837./8709120.)*i10e6 - (3498923./29030400.)*i10e4 + (
                        95268679979./18579456000.)*i8e12 - (51027533./7741440.)*i8e10 + (176956279./23224320.)*i8e8 - (
                        765187./193536.)*i8e6 + (565573./645120.)*i8e4 + (13087364543./1306368000.)*i6e14 - (
                        21283933./829440.)*i6e12 + (24247471./691200.)*i6e10 - (16173799./414720.)*i6e8 + (
                        8993./432.)*i6e6 - (6647./1440.)*i6e4 + (264182891963./11147673600.)*i4e16 - (
                        7055240243./174182400.)*i4e14 + (241167361./2764800.)*i4e12 - (2344199./18432.)*i4e10 + (
                        7512571./55296.)*i4e8 - (21505./288.)*i4e6 + (3179./192.)*i4e4 + (
                        2949969133./182891520.)*i2e18 - (173370469./4147200.)*i2e16 + (346700573./3628800.)*i2e14 - (
                        8421731./46080.)*i2e12 + (134209./480.)*i2e10 - (83551./288.)*i2e8 + (1955./12.)*i2e6 - (
                        289./8.)*i2e4 + (38506861007131./7023034368000.)*e20 - (2949969133./182891520.)*e18 + (
                        173370469./4147200.)*e16 - (346700573./3628800.)*e14 + (8421731./46080.)*e12 - (
                        134209./480.)*e10 + (83551./288.)*e8 - (1955./12.)*e6 + (289./8.)*e4,
            2.
        ),
        (
            '-n-o',
            -orbital_freq - spin_freq,
            (8823166687./5690998849536000.)*i18e2 - (22303343./344408064000.)*i16e4 - (6263489./211341312000.)*i16e2 + (
                        346689850051./200858782924800.)*i14e6 + (301124023./309967257600.)*i14e4 + (
                        34420297./77491814400.)*i14e2 - (133152775507./4414478745600.)*i12e8 - (
                        2166757111./110361968640.)*i12e6 - (9411719./851558400.)*i12e4 - (1074041./212889600.)*i12e2 + (
                        71790682171./203843174400.)*i10e10 + (665745497./2675441664.)*i10e8 + (
                        1354163761./8360755200.)*i10e6 + (130757./1433600.)*i10e4 + (133907./3225600.)*i10e2 - (
                        28386502727599./10701766656000.)*i8e12 - (39259718513./19818086400.)*i8e10 - (
                        260044601./185794560.)*i8e8 - (4231445./4644864.)*i8e6 - (18397./35840.)*i8e4 - (
                        2083./8960.)*i8e2 + (511253620507703./42807066624000.)*i6e14 + (
                        7096471588921./764411904000.)*i6e12 + (1962921403./283115520.)*i6e10 + (
                        13001251./2654208.)*i6e8 + (5288671./1658880.)*i6e6 + (4603./2560.)*i6e4 + (517./640.)*i6e2 - (
                        148892609051513./5327101624320.)*i4e16 - (15976372618603./713451110400.)*i4e14 - (
                        110879314849./6370099200.)*i4e12 - (153346823./11796480.)*i4e10 - (2031247./221184.)*i4e8 - (
                        165245./27648.)*i4e6 - (27./8.)*i4e4 - (3./2.)*i4e2 + (
                        294008809328516923./11506539508531200.)*i2e18 + (46529133791863./2219625676800.)*i2e16 + (
                        499264587977./29727129600.)*i2e14 + (13860031763./1061683200.)*i2e12 + (
                        3833717./393216.)*i2e10 + (126955./18432.)*i2e8 + (1291./288.)*i2e6 + (81./32.)*i2e4 + (
                        9./8.)*i2e2,
            1.
        ),
        (
            'n-o',
            orbital_freq - spin_freq,
            (15901071061./10243797929164800.)*i18e2 - (5418367139./83691159552000.)*i16e4 - (
                        88775803./2988969984000.)*i16e2 + (115570856383./66952927641600.)*i14e6 + (
                        541784093./557941063680.)*i14e4 + (5654153./12680478720.)*i14e2 - (
                        133162540717./4414478745600.)*i12e8 - (144467303./7357464576.)*i12e6 - (
                        84639041./7664025600.)*i12e4 - (9732799./1916006400.)*i12e2 + (
                        119658011123./339738624000.)*i10e10 + (16645807439./66886041600.)*i10e8 + (
                        451481413./2786918400.)*i10e6 + (2115311./23224320.)*i10e4 + (6971./165888.)*i10e2 - (
                        350474324849./132120576000.)*i8e12 - (13087909301./6606028800.)*i8e10 - (
                        260104871./185794560.)*i8e8 - (1411001./1548288.)*i8e6 - (165163./322560.)*i8e4 - (
                        19157./80640.)*i8e2 + (1363455696287./114152177664.)*i6e14 + (
                        17524253851./1887436800.)*i6e12 + (3272128931./471859200.)*i6e10 + (
                        65033009./13271040.)*i6e8 + (1764043./552960.)*i6e6 + (8249./4608.)*i6e4 + (967./1152.)*i6e2 - (
                        63817740372059./2283043553280.)*i4e16 - (5326231634771./237817036800.)*i4e14 - (
                        6161261113./353894400.)*i4e12 - (25565953./1966080.)*i4e10 - (2032717./221184.)*i4e8 - (
                        55145./9216.)*i4e6 - (643./192.)*i4e4 - (77./48.)*i4e2 + (
                        49007498831132419./1917756584755200.)*i2e18 + (69805289550263./3329438515200.)*i2e16 + (
                        13316818727./792723456.)*i2e14 + (962827./73728.)*i2e12 + (1598197./163840.)*i2e10 + (
                        63551./9216.)*i2e8 + (3449./768.)*i2e6 + (5./2.)*i2e4 + (5./4.)*i2e2,
            1.
        ),
        (
            '3n-o',
            3.*orbital_freq - spin_freq,
            (2195943977./1045285502976000.)*i18e2 + (143096929./797058662400.)*i16e4 - (3490169./85399142400.)*i16e2 + (
                        291777280313./44635285094400.)*i14e6 - (368529607./132843110400.)*i14e4 + (
                        8988527./14233190400.)*i14e2 + (871197841./46714060800.)*i12e8 - (
                        9441589609./122624409600.)*i12e6 + (2424043./72990720.)*i12e4 - (413861./54743040.)*i12e2 + (
                        412540295467./792723456000.)*i10e10 - (912783293./4954521600.)*i10e8 + (
                        1253268323./1857945600.)*i10e6 - (1651357./5529600.)*i10e4 + (281939./4147200.)*i10e2 - (
                        48254503933./26424115200.)*i8e12 - (6801059227./2202009600.)*i8e10 + (
                        18592061./13762560.)*i8e8 - (21684707./5160960.)*i8e6 + (5945./3072.)*i8e4 - (
                        1015./2304.)*i8e2 + (1524246268151./132120576000.)*i6e14 + (55926589319./9437184000.)*i6e12 + (
                        1870083469./157286400.)*i6e10 - (6718523./983040.)*i6e8 + (6426533./368640.)*i6e6 - (
                        65149./7680.)*i6e4 + (11123./5760.)*i6e2 - (1132870187049./46976204800.)*i4e16 - (
                        27994550407./1258291200.)*i4e14 - (14687593./1572864.)*i4e12 - (17307023./655360.)*i4e10 + (
                        170741./8192.)*i4e8 - (132347./3072.)*i4e6 + (1435./64.)*i4e4 - (245./48.)*i4e2 + (
                        17066354521071./751619276800.)*i2e18 + (208688234697./11744051200.)*i2e16 + (
                        924913403./52428800.)*i2e14 + (55810297./13107200.)*i2e12 + (17156499./655360.)*i2e10 - (
                        57471./2048.)*i2e8 + (1549./32.)*i2e6 - (861./32.)*i2e4 + (49./8.)*i2e2,
            1.
        ),
        (
            '-2n-o',
            -2.*orbital_freq - spin_freq,
            (-166711./23911759872000.)*i20 - (2743907./3347646382080.)*i18e2 + (2743907./16738231910400.)*i18 - (
                        46081037./516612096000.)*i16e4 + (437./28385280.)*i16e2 - (437./141926400.)*i16 + (
                        16524163./12454041600.)*i14e6 + (7612757./5740134400.)*i14e4 - (173./774144.)*i14e2 + (
                        173./3870720.)*i14 - (31824774733./1103619686400.)*i12e8 - (3644933./239500800.)*i12e6 - (
                        234181./15769600.)*i12e4 + (463./193536.)*i12e2 - (463./967680.)*i12 + (
                        27801816239./83607552000.)*i10e10 + (567377177./2388787200.)*i10e8 + (7337./57600.)*i10e6 + (
                        85381./716800.)*i10e4 - (9./512.)*i10e2 + (9./2560.)*i10 - (1996744519./796262400.)*i8e12 - (
                        434201279./232243200.)*i8e10 - (61884479./46448640.)*i8e8 - (1057./1440.)*i8e6 - (
                        11421./17920.)*i8e4 + (5./64.)*i8e2 - (1./64.)*i8 + (1649647706809./146313216000.)*i6e14 + (
                        17472864227./1990656000.)*i6e12 + (108469001./16588800.)*i6e10 + (15401609./3317760.)*i6e8 + (
                        1921./720.)*i6e6 + (2619./1280.)*i6e4 - (5./32.)*i6e2 + (1./32.)*i6 - (
                        73285635553./2786918400.)*i4e16 - (115052561./5443200.)*i4e14 - (853261./51840.)*i4e12 - (
                        423287./34560.)*i4e10 - (59801./6912.)*i4e8 - (21./4.)*i4e6 - (27./8.)*i4e4 + (
                        1400928769637./58525286400.)*i2e18 + (18321832657./928972800.)*i2e16 + (
                        230112373./14515200.)*i2e14 + (68263727./5529600.)*i2e12 + (105827./11520.)*i2e10 + (
                        14951./2304.)*i2e8 + (63./16.)*i2e6 + (81./32.)*i2e4,
            1.
        ),
        (
            '3n',
            3.*orbital_freq,
            (-257./83397600.)*i18e2 - (223901./864864000.)*i16e4 + (5461./92664000.)*i16e2 - (
                        408042079./43589145600.)*i14e6 + (41./10560.)*i14e4 - (7./7920.)*i14e2 - (
                        53550377./2554675200.)*i12e8 + (101969503./958003200.)*i12e6 - (1271./28800.)*i12e4 + (
                        217./21600.)*i12e2 - (9012556853./12386304000.)*i10e10 + (1894991./11059200.)*i10e8 - (
                        25451359./29030400.)*i10e6 + (697./1920.)*i10e4 - (119./1440.)*i10e2 + (
                        99039786343./33030144000.)*i8e12 + (1607718319./393216000.)*i8e10 - (3193577./3440640.)*i8e8 + (
                        6321823./1290240.)*i8e6 - (2583./1280.)*i8e4 + (147./320.)*i8e2 - (
                        250635583131./14680064000.)*i6e14 - (5086224779./471859200.)*i6e12 - (
                        112073569./7864320.)*i6e10 + (675737./245760.)*i6e8 - (1539439./92160.)*i6e6 + (
                        861./128.)*i6e4 - (49./32.)*i6e2 + (537105236547./13421772800.)*i4e16 + (
                        197651733117./5872025600.)*i4e14 + (7028958283./314572800.)*i4e12 + (
                        688791973./26214400.)*i4e10 - (46277./32768.)*i4e8 + (343843./12288.)*i4e6 - (
                        2583./256.)*i4e4 + (147./64.)*i4e2 - (8454941691423./187904819200.)*i2e18 - (
                        437242738059./11744051200.)*i2e16 - (3165427449./104857600.)*i2e14 - (
                        30355211./1310720.)*i2e12 - (6019881./327680.)*i2e10 - (20829./2048.)*i2e8 - (
                        2809./256.)*i2e6 + (466859091912363./26306674688000.)*e20 + (
                        2818313897141./187904819200.)*e18 + (145747579353./11744051200.)*e16 + (
                        1055142483./104857600.)*e14 + (30355211./3932160.)*e12 + (2006627./327680.)*e10 + (
                        6943./2048.)*e8 + (2809./768.)*e6,
            0.
        ),
        (
            '-n-2o',
            -orbital_freq - 2.*spin_freq,
            (-8644781./22317642547200.)*i18e2 + (802969931./49594761216000.)*i16e4 + (
                        91760309./12398690304000.)*i16e2 - (10831066667./25107347865600.)*i14e6 - (
                        20681./85155840.)*i14e4 - (337./3041280.)*i14e2 + (2418511621./321052999680.)*i12e8 + (
                        10822609811./2207239372800.)*i12e6 + (285031./103219200.)*i12e4 + (3601./2867200.)*i12e2 - (
                        625699035023./7134511104000.)*i10e10 - (592059833./9555148800.)*i10e8 - (
                        337187791./8360755200.)*i10e6 - (19549./860160.)*i10e4 - (737./71680.)*i10e2 + (
                        55884570305471./85614133248000.)*i8e12 + (15457945133./31708938240.)*i8e10 + (
                        102383837./297271296.)*i8e8 + (41647601./185794560.)*i8e6 + (5179./40960.)*i8e4 + (
                        581./10240.)*i8e2 - (59911530765031./21403533312000.)*i6e14 - (
                        83159752259./38220595200.)*i6e12 - (575053217./353894400.)*i6e10 - (7617223./6635520.)*i6e8 - (
                        619673./829440.)*i6e6 - (27./64.)*i6e4 - (3./16.)*i6e2 + (
                        558354764522761./106542032486400.)*i4e16 + (11982507105883./2853804441600.)*i4e14 + (
                        83161756003./25480396800.)*i4e12 + (23002921./9437184.)*i4e10 + (1523515./884736.)*i4e8 + (
                        123941./110592.)*i4e6 + (81./128.)*i4e4 + (9./32.)*i4e2 - (
                        167431204877./821895679180800.)*i2e18 - (147400583./634178764800.)*i2e16 - (
                        31398887./118908518400.)*i2e14 - (62617./212336640.)*i2e12 - (619./1966080.)*i2e10 - (
                        11./36864.)*i2e8 - (1./4608.)*i2e6 + (20572630185001./115065395085312000.)*e20 + (
                        167431204877./821895679180800.)*e18 + (147400583./634178764800.)*e16 + (
                        31398887./118908518400.)*e14 + (62617./212336640.)*e12 + (619./1966080.)*e10 + (
                        11./36864.)*e8 + (1./4608.)*e6,
            2.
        ),
        (
            '-2n-2o',
            -2.*orbital_freq - 2.*spin_freq,
            (164573./95647039488000.)*i20 + (491./2452488192.)*i18e2 - (491./12262440960.)*i18 + (
                        20249169./918421504000.)*i16e4 - (457./123863040.)*i16e2 + (457./619315200.)*i16 - (
                        8027./23950080.)*i14e6 - (1537./4730880.)*i14e4 + (5./96768.)*i14e2 - (1./96768.)*i14 + (
                        31728878989./4414478745600.)*i12e8 + (16051./4147200.)*i12e6 + (20527./5734400.)*i12e4 - (
                        19./36864.)*i12e2 + (19./184320.)*i12 - (6919375601./83607552000.)*i10e10 - (
                        985676273./16721510400.)*i10e8 - (1129./34560.)*i10e6 - (3999./143360.)*i10e4 + (
                        5./1536.)*i10e2 - (1./1536.)*i10 + (137600820757./222953472000.)*i8e12 + (
                        122014537./265420800.)*i8e10 + (121207039./371589120.)*i8e8 + (2173./11520.)*i8e6 + (
                        2907./20480.)*i8e4 - (5./512.)*i8e2 + (1./512.)*i8 - (3451617919./1306368000.)*i6e14 - (
                        1023929333./497664000.)*i6e12 - (396839./259200.)*i6e10 - (112129./103680.)*i6e8 - (
                        21./32.)*i6e6 - (27./64.)*i6e4 + (6871040387./1393459200.)*i4e16 + (
                        690385459./174182400.)*i4e14 + (204810161./66355200.)*i4e12 + (79379./34560.)*i4e10 + (
                        22429./13824.)*i4e8 + (63./64.)*i4e6 + (81./128.)*i4e4 - (128441./119439360.)*i2e18 - (
                        22601./18579456.)*i2e16 - (2417./1814400.)*i2e14 - (949./691200.)*i2e12 - (7./5760.)*i2e10 - (
                        1./1152.)*i2e8 + (3291434567./3511517184000.)*e20 + (128441./119439360.)*e18 + (
                        22601./18579456.)*e16 + (2417./1814400.)*e14 + (949./691200.)*e12 + (7./5760.)*e10 + (
                        1./1152.)*e8,
            2.
        ),
        (
            '-o',
            -spin_freq,
            (-262144./9280784638125.)*i20 + (65536./32564156625.)*i18e2 + (65536./97692469875.)*i18 - (
                        16384./212837625.)*i16e4 - (8192./212837625.)*i16e2 - (8192./638512875.)*i16 + (
                        16384./8513505.)*i14e6 + (16384./14189175.)*i14e4 + (8192./14189175.)*i14e2 + (
                        8192./42567525.)*i14 - (1024./31185.)*i12e8 - (2048./93555.)*i12e6 - (2048./155925.)*i12e4 - (
                        1024./155925.)*i12e2 - (1024./467775.)*i12 + (256./675.)*i10e10 + (256./945.)*i10e8 + (
                        512./2835.)*i10e6 + (512./4725.)*i10e4 + (256./4725.)*i10e2 + (256./14175.)*i10 - (
                        128./45.)*i8e12 - (32./15.)*i8e10 - (32./21.)*i8e8 - (64./63.)*i8e6 - (64./105.)*i8e4 - (
                        32./105.)*i8e2 - (32./315.)*i8 + (64./5.)*i6e14 + (448./45.)*i6e12 + (112./15.)*i6e10 + (
                        16./3.)*i6e8 + (32./9.)*i6e6 + (32./15.)*i6e4 + (16./15.)*i6e2 + (
                        16./45.)*i6 - 30*i4e16 - 24*i4e14 - (56./3.)*i4e12 - 14*i4e10 - 10*i4e8 - (
                        20./3.)*i4e6 - 4*i4e4 - 2*i4e2 - (2./3.)*i4 + (55./2.)*i2e18 + (
                        45./2.)*i2e16 + 18*i2e14 + 14*i2e12 + (21./2.)*i2e10 + (15./2.)*i2e8 + 5*i2e6 + 3*i2e4 + (
                        3./2.)*i2e2 + (1./2.)*i2,
            1.
        ),
        (
            '2n-o',
            2.*orbital_freq - spin_freq,
            (-3479571749./486580401635328000.)*i20 - (2195943977./2560949482291200.)*i18e2 + (
                        2195943977./12804747411456000.)*i18 - (424145803./4649508864000.)*i16e4 + (
                        3490169./209227898880.)*i16e2 - (3490169./1046139494400.)*i16 + (
                        4046247103./3138418483200.)*i14e6 + (213914633./154983628800.)*i14e4 - (
                        8988527./34871316480.)*i14e2 + (8988527./174356582400.)*i14 - (
                        32444130643./1103619686400.)*i12e8 - (125282063./8622028800.)*i12e6 - (
                        6787897./425779200.)*i12e4 + (59123./19160064.)*i12e2 - (59123./95800320.)*i12 + (
                        45786897463./139345920000.)*i10e10 + (813559357./3344302080.)*i10e8 + (
                        15397573./130636800.)*i10e6 + (871763./6451200.)*i10e4 - (40277./1451520.)*i10e2 + (
                        40277./7257600.)*i10 - (14004087071./5573836800.)*i8e12 - (142872011./77414400.)*i8e10 - (
                        63912689./46448640.)*i8e8 - (233029./362880.)*i8e6 - (14291./17920.)*i8e4 + (
                        725./4032.)*i8e2 - (145./4032.)*i8 + (549378057479./48771072000.)*i6e14 + (
                        3502573193./398131200.)*i6e12 + (178083313./27648000.)*i6e10 + (3229435./663552.)*i6e8 + (
                        54643./25920.)*i6e6 + (3893./1280.)*i6e4 - (227./288.)*i6e2 + (227./1440.)*i6 - (
                        37995684851./1445068800.)*i4e16 - (8581440353./406425600.)*i4e14 - (
                        54777263./3317760.)*i4e12 - (110659./9216.)*i4e10 - (257669./27648.)*i4e8 - (739./216.)*i4e6 - (
                        213./32.)*i4e4 + (25./12.)*i4e2 - (5./12.)*i4 + (6302286505933./263363788800.)*i2e18 + (
                        2374624949./120422400.)*i2e16 + (1072027567./67737600.)*i2e14 + (13720169./1105920.)*i2e12 + (
                        341669./38400.)*i2e10 + (33595./4608.)*i2e8 + (251./144.)*i2e6 + (207./32.)*i2e4 - (
                        5./2.)*i2e2 + (1./2.)*i2,
            1.
        ),
        (
            'n',
            orbital_freq,
            (-89809./38594556000.)*i18e2 + (1764047./18162144000.)*i16e4 + (202073./4540536000.)*i16e2 - (
                        1881127./726485760.)*i14e6 - (147013./100900800.)*i14e4 - (50521./75675600.)*i14e2 + (
                        43359671./958003200.)*i12e8 + (1411189./47900160.)*i12e6 + (110287./6652800.)*i12e4 + (
                        4211./554400.)*i12e2 - (2184040423./4128768000.)*i10e10 - (3616829./9676800.)*i10e8 - (
                        353141./1451520.)*i10e6 - (27599./201600.)*i10e4 - (3161./50400.)*i10e2 + (
                        396322412057./99090432000.)*i8e12 + (328883377./110100480.)*i8e10 + (2723191./1290240.)*i8e8 + (
                        29543./21504.)*i8e6 + (6927./8960.)*i8e4 + (793./2240.)*i8e2 - (
                        5436680684213./297271296000.)*i6e14 - (301860389467./21233664000.)*i6e12 - (
                        417490789./39321600.)*i6e10 - (691367./92160.)*i6e8 - (22501./4608.)*i6e6 - (
                        1759./640.)*i6e4 - (201./160.)*i6e2 + (80652231656077./1775700541440.)*i4e16 + (
                        1442363150713./39636172800.)*i4e14 + (80083677791./2831155200.)*i4e12 + (
                        22151853./1048576.)*i4e10 + (61137./4096.)*i4e8 + (29845./3072.)*i4e6 + (1401./256.)*i4e4 + (
                        159./64.)*i4e2 - (725941889609009./14205604331520.)*i2e18 - (
                        689312857627./16441671680.)*i2e16 - (24654653741./734003200.)*i2e14 - (
                        85553819./3276800.)*i2e12 - (3194661./163840.)*i2e10 - (28211./2048.)*i2e8 - (
                        2295./256.)*i2e6 - (81./16.)*i2e4 - (9./4.)*i2e2 + (
                        2343932511383522191./115065395085312000.)*e20 + (725941889609009./42616812994560.)*e18 + (
                        689312857627./49325015040.)*e16 + (24654653741./2202009600.)*e14 + (85553819./9830400.)*e12 + (
                        1064887./163840.)*e10 + (28211./6144.)*e8 + (765./256.)*e6 + (27./16.)*e4 + (3./4.)*e2,
            0.
        ),
        (
            '3n-2o',
            3.*orbital_freq - 2.*spin_freq,
            (-561142037./1045285502976000.)*i18e2 - (2973253457./63764692992000.)*i16e4 + (
                        72518377./6831931392000.)*i16e2 - (635755327./371960709120.)*i14e6 + (
                        1118029./1509580800.)*i14e4 - (27269./161740800.)*i14e2 - (7678265783./1307993702400.)*i12e8 + (
                        3400796131./163499212800.)*i12e6 - (13545047./1459814400.)*i12e4 + (
                        2312569./1094860800.)*i12e2 - (12271338893./88080384000.)*i10e10 + (
                        45900997./707788800.)*i10e8 - (7949713./41287680.)*i10e6 + (496387./5529600.)*i10e4 - (
                        84749./4147200.)*i10e2 + (6732799983./16777216000.)*i8e12 + (5118548329./5872025600.)*i8e10 - (
                        60573773./110100480.)*i8e8 + (18110321./13762560.)*i8e6 - (80237./122880.)*i8e4 + (
                        13699./92160.)*i8e2 - (107250025219./37748736000.)*i6e14 - (24438221./24576000.)*i6e12 - (
                        12131629./3276800.)*i6e10 + (832621./245760.)*i6e8 - (39313./6144.)*i6e6 + (
                        6601./1920.)*i6e4 - (1127./1440.)*i6e2 + (114520355859./26843545600.)*i4e16 + (
                        24869052487./5033164800.)*i4e14 - (1510073./3276800.)*i4e12 + (13218507./1310720.)*i4e10 - (
                        456169./32768.)*i4e8 + (86193./4096.)*i4e6 - (3157./256.)*i4e4 + (539./192.)*i4e2 - (
                        6258845529./30064771072.)*i2e18 + (3973253733./4697620480.)*i2e16 - (
                        534226163./209715200.)*i2e14 + (47982879./6553600.)*i2e12 - (5568309./327680.)*i2e10 + (
                        135771./4096.)*i2e8 - (21975./512.)*i2e6 + (861./32.)*i2e4 - (49./8.)*i2e2 - (
                        8388292638491./105226698752000.)*e20 + (6258845529./30064771072.)*e18 - (
                        3973253733./4697620480.)*e16 + (534226163./209715200.)*e14 - (47982879./6553600.)*e12 + (
                        5568309./327680.)*e10 - (135771./4096.)*e8 + (21975./512.)*e6 - (861./32.)*e4 + (49./8.)*e2,
            2.
        ),
        (
            'n-2o',
            orbital_freq - 2.*spin_freq,
            (-19888199957./51218989645824000.)*i18e2 + (4333885291./267811710566400.)*i16e4 + (
                        7656277./1030045040640.)*i16e2 - (277756067./643778150400.)*i14e6 - (
                        84613961./348713164800.)*i14e4 - (9734839./87178291200.)*i14e2 + (
                        133037669869./17657914982400.)*i12e8 + (3608378063./735746457600.)*i12e6 + (
                        16904269./6131220480.)*i12e4 + (1951667./1532805120.)*i12e2 - (
                        208590407867./2378170368000.)*i10e10 - (118442971./1911029760.)*i10e8 - (
                        112442683./2786918400.)*i10e6 - (2631733./116121600.)*i10e4 - (305867./29030400.)*i10e2 + (
                        248410000001./380507258880.)*i8e12 + (25768594261./52848230400.)*i8e10 + (
                        512160559./1486356480.)*i8e8 + (13892933./61931520.)*i8e6 + (64927./516096.)*i8e4 + (
                        7649./129024.)*i8e2 - (19974071093099./7134511104000.)*i6e14 - (
                        23105907397./10616832000.)*i6e12 - (95879689./58982400.)*i6e10 - (1524797./1327104.)*i6e8 - (
                        206849./276480.)*i6e6 - (2407./5760.)*i6e4 - (293./1440.)*i6e2 + (
                        1675574203516897./319626097459200.)*i4e16 + (114167772709./27179089920.)*i4e14 + (
                        924591149./283115520.)*i4e12 + (19187029./7864320.)*i4e10 + (1526749./884736.)*i4e8 + (
                        41453./36864.)*i4e6 + (475./768.)*i4e4 + (65./192.)*i4e2 - (
                        12842565048623./3835513169510400.)*i2e18 - (49450862117./13317754060800.)*i2e16 - (
                        165285343./39636172800.)*i2e14 - (277229./58982400.)*i2e12 - (1733./327680.)*i2e10 - (
                        305./36864.)*i2e8 - (13./1536.)*i2e6 + (1./32.)*i2e4 - (1./8.)*i2e2 + (
                        2104490540764777./690392370511872000.)*e20 + (12842565048623./3835513169510400.)*e18 + (
                        49450862117./13317754060800.)*e16 + (165285343./39636172800.)*e14 + (277229./58982400.)*e12 + (
                        1733./327680.)*e10 + (305./36864.)*e8 + (13./1536.)*e6 - (1./32.)*e4 + (1./8.)*e2,
            2.
        ),
        (
            '2n',
            2.*orbital_freq,
            (73./6891885000.)*i20 + (257./204324120.)*i18e2 - (257./1021620600.)*i18 + (45511./336336000.)*i16e4 - (
                        5461./227026800.)*i16e2 + (5461./1135134000.)*i16 - (63659./32432400.)*i14e6 - (
                        3103./1528800.)*i14e4 + (1./2772.)*i14e2 - (1./13860.)*i14 + (3801619./87091200.)*i12e8 + (
                        334343./14968800.)*i12e6 + (1219./52800.)*i12e4 - (31./7560.)*i12e2 + (31./37800.)*i12 - (
                        72153913./145152000.)*i10e10 - (1162633./3225600.)*i10e8 - (83719./453600.)*i10e6 - (
                        2133./11200.)*i10e4 + (17./504.)*i10e2 - (17./2520.)*i10 + (976317787./258048000.)*i8e12 + (
                        10060853./3584000.)*i8e10 + (2625169./1290240.)*i8e8 + (1003./960.)*i8e6 + (
                        4797./4480.)*i8e4 - (3./16.)*i8e2 + (3./80.)*i8 - (11692772467./677376000.)*i6e14 - (
                        29743849./2211840.)*i6e12 - (4598441./460800.)*i6e10 - (665537./92160.)*i6e8 - (
                        5399./1440.)*i6e6 - (1197./320.)*i6e4 + (5./8.)*i6e2 - (1./8.)*i6 + (
                        741056268847./17340825600.)*i4e16 + (9306367493./270950400.)*i4e14 + (
                        197256941./7372800.)*i4e12 + (6103337./307200.)*i4e10 + (58543./4096.)*i4e8 + (
                        1483./192.)*i4e6 + (891./128.)*i4e4 - (15./16.)*i4e2 + (3./16.)*i4 - (
                        51883919761./1083801600.)*i2e18 - (678544541./17203200.)*i2e16 - (25565893./806400.)*i2e14 - (
                        505601./20480.)*i2e12 - (11757./640.)*i2e10 - (1661./128.)*i2e8 - (63./8.)*i2e6 - (
                        81./16.)*i2e4 + (6351400670087./334430208000.)*e20 + (51883919761./3251404800.)*e18 + (
                        678544541./51609600.)*e16 + (25565893./2419200.)*e14 + (505601./61440.)*e12 + (
                        3919./640.)*e10 + (1661./384.)*e8 + (21./8.)*e6 + (27./16.)*e4,
            0.
        ),
        (
            '2n-2o',
            2.*orbital_freq - 2.*spin_freq,
            (17616175987./9731608032706560000.)*i20 + (561142037./2560949482291200.)*i18e2 - (
                        561142037./12804747411456000.)*i18 + (131962651./5722472448000.)*i16e4 - (
                        72518377./16738231910400.)*i16e2 + (72518377./83691159552000.)*i16 - (
                        124902599./392302310400.)*i14e6 - (973879./2767564800.)*i14e4 + (27269./396264960.)*i14e2 - (
                        27269./1981324800.)*i14 + (6508373303./882895749120.)*i12e8 + (
                        122391503./34488115200.)*i12e6 + (7026553./1703116800.)*i12e4 - (330367./383201280.)*i12e2 + (
                        330367./1916006400.)*i12 - (11379001513./139345920000.)*i10e10 - (
                        1020649199./16721510400.)*i10e8 - (3670267./130636800.)*i10e6 - (231629./6451200.)*i10e4 + (
                        12107./1451520.)*i10e2 - (12107./7257600.)*i10 + (27588572527./44590694400.)*i8e12 + (
                        1400534183./3096576000.)*i8e10 + (25557005./74317824.)*i8e8 + (416933./2903040.)*i8e6 + (
                        4549./20480.)*i8e4 - (1957./32256.)*i8e2 + (1957./161280.)*i8 - (
                        32173006583./12192768000.)*i6e14 - (102780619./49766400.)*i6e12 - (10327357./6912000.)*i6e10 - (
                        981971./829440.)*i6e8 - (4871./12960.)*i6e6 - (37./40.)*i6e4 + (23./72.)*i6e2 - (
                        23./360.)*i6 + (28492689943./5780275200.)*i4e16 + (6423465707./1625702400.)*i4e14 + (
                        20666431./6635520.)*i4e12 + (1994809./921600.)*i4e10 + (220055./110592.)*i4e8 - (
                        37./1728.)*i4e6 + (39./16.)*i4e4 - (55./48.)*i4e2 + (11./48.)*i4 + (
                        459927151./75246796800.)*i2e18 + (561889./240844800.)*i2e16 + (1739939./67737600.)*i2e14 - (
                        34471./552960.)*i2e12 + (11041./38400.)*i2e10 - (3697./4608.)*i2e8 + (79./36.)*i2e6 - (
                        63./16.)*i2e4 + (5./2.)*i2e2 - (1./2.)*i2 - (475843828001./105345515520000.)*e20 - (
                        459927151./75246796800.)*e18 - (561889./240844800.)*e16 - (1739939./67737600.)*e14 + (
                        34471./552960.)*e12 - (11041./38400.)*e10 + (3697./4608.)*e8 - (79./36.)*e6 + (63./16.)*e4 - (
                        5./2.)*e2 + (1./2.),
            2.
        )
    )

    heating_coeffs = list()
    ztorque_coeffs = list()
    freqs = list()
    modes = list()
    mode_names = list()

    for mode_name, mode, coeff, torque_multi in modes_coeffs:
        sgn = np.sign(mode)
        freq = np.abs(mode)

        heating_coeff = freq * coeff
        ztorque_coeff = coeff * torque_multi * sgn

        freq_too_low = freq < MODE_ZERO_TOL
        freq[freq_too_low] = 0.
        mode[freq_too_low] = 0.
        heating_coeff[freq_too_low] = 0.
        ztorque_coeff[freq_too_low] = 0.

        # if freq < MODE_ZERO_TOL:
        #     freq = 0.
        #     heating_coeff = 0.
        #     ztorque_coeff = 0.

        freqs.append(freq)
        modes.append(mode)
        mode_names.append(mode_name)
        heating_coeffs.append(heating_coeff)
        ztorque_coeffs.append(ztorque_coeff)

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return mode_names, modes, freqs, heating_coeffs, ztorque_coeffs