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