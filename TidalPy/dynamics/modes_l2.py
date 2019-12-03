from typing import Tuple, List

import numpy as np

from TidalPy.performance import njit
from . import MODE_ZERO_TOL

ModeOutput = Tuple[List[str, ...], List[np.ndarray, ...], List[np.ndarray, ...], List[np.ndarray, ...],
                   List[np.ndarray, ...]]


# TODO: A lot of the mode list construction is static and could take place outside of the functions at the expense of
#   readibility.

@njit
def spin_sync_modes(orbital_frequency: np.ndarray, spin_frequency: np.ndarray, eccentricity: np.ndarray,
                    inclination: np.ndarray) -> ModeOutput:
    """ Synchronous Rotation Tidal Modes (for heating and torque) for l=2, e^2, I^2

    These should all have the same shape!!

    Parameters
    ----------
    orbital_frequency : np.ndarray
        Orbital frequency (mean motion) [rads s-1]
    spin_frequency : np.ndarray
        Planet's spin frequency (inverse of period) [rads s-1]
    eccentricity : np.ndarray
        Orbital eccentricity
    inclination : np.ndarray
        Orbital inclination (relative to the orbital plane of the satellite) [rads]
        Note: this is not relative to the host planet's equator

    Returns
    -------
    mode_names : List[str, ...]
        List of easy to read names for each mode.
            "n" == orbital motion frequency
            "o" == spin frequency
    modes : List[np.ndarray, ...]
        List of calculated modes [rads s-1]
    freqs : List[np.ndarray, ...]
        List of calculated frequencies (absolute value of modes) [rads s-1]
    heating_subterms: List[np.ndarray, ...]
        List of tidal heating coefficients [N m]
    ztorque_subterms: List[np.ndarray, ...]
        List of tidal torque coefficients [N m]
    """

    e2 = eccentricity ** 2
    i2 = inclination ** 2

    modes_coeffs = (
        (
            'n',
            orbital_frequency,
            7. * e2 + i2,
            12. * e2
        ),
    )

    heating_subterms = list()
    ztorque_subterms = list()
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
        heating_subterms.append(heating_coeff)
        ztorque_subterms.append(ztorque_coeff)

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return mode_names, modes, freqs, heating_subterms, ztorque_subterms


@njit
def spin_sync_modes_4(orbital_frequency: np.ndarray, spin_frequency: np.ndarray, eccentricity: np.ndarray,
                      inclination: np.ndarray) -> ModeOutput:
    """ Synchronous Rotation Tidal Modes (for heating and torque) for l=2, e^4, I^4

    These should all have the same shape!!

    Parameters
    ----------
    orbital_frequency : np.ndarray
        Orbital frequency (mean motion) [rads s-1]
    spin_frequency : np.ndarray
        Planet's spin frequency (inverse of period) [rads s-1]
    eccentricity : np.ndarray
        Orbital eccentricity
    inclination : np.ndarray
        Orbital inclination (relative to the orbital plane of the satellite) [rads]
        Note: this is not relative to the host planet's equator

    Returns
    -------
    mode_names : List[str, ...]
        List of easy to read names for each mode.
            "n" == orbital motion frequency
            "o" == spin frequency
    modes : List[np.ndarray, ...]
        List of calculated modes [rads s-1]
    freqs : List[np.ndarray, ...]
        List of calculated frequencies (absolute value of modes) [rads s-1]
    heating_subterms: List[np.ndarray, ...]
        List of tidal heating coefficients [N m]
    ztorque_subterms: List[np.ndarray, ...]
        List of tidal torque coefficients [N m]
    """

    e2 = eccentricity**2
    e4 = eccentricity**4
    i2 = inclination**2
    i4 = inclination**4
    i2e2 = e2 * i2

    modes_coeffs = (
        (
            'n',
            orbital_frequency,
            7. * e2 - (101./4.) * e4 + i2 - (13./12.) * i4 - (19./2.) * i2e2,
            12. * e2 - (215./4.) * e4 + (1./4.) * i4 - 16. * i2e2
        ),
        (
            '2n',
            2. * orbital_frequency,
            (605./16.) * e4 + (5./16.) * i4 + (29./4.) * i2e2,
            (289./4.) * e4 - (1./4.) * i4 + 5. * i2e2
        ),
    )

    heating_subterms = list()
    ztorque_subterms = list()
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
        heating_subterms.append(heating_coeff)
        ztorque_subterms.append(ztorque_coeff)

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return mode_names, modes, freqs, heating_subterms, ztorque_subterms


@njit
def spin_sync_modes_6(orbital_frequency: np.ndarray, spin_frequency: np.ndarray, eccentricity: np.ndarray,
                      inclination: np.ndarray) -> ModeOutput:
    """ Synchronous Rotation Tidal Modes (for heating and torque) for l=2, e^6, I^6

    These should all have the same shape!!

    Parameters
    ----------
    orbital_frequency : np.ndarray
        Orbital frequency (mean motion) [rads s-1]
    spin_frequency : np.ndarray
        Planet's spin frequency (inverse of period) [rads s-1]
    eccentricity : np.ndarray
        Orbital eccentricity
    inclination : np.ndarray
        Orbital inclination (relative to the orbital plane of the satellite) [rads]
        Note: this is not relative to the host planet's equator

    Returns
    -------
    mode_names : List[str, ...]
        List of easy to read names for each mode.
            "n" == orbital motion frequency
            "o" == spin frequency
    modes : List[np.ndarray, ...]
        List of calculated modes [rads s-1]
    freqs : List[np.ndarray, ...]
        List of calculated frequencies (absolute value of modes) [rads s-1]
    heating_subterms: List[np.ndarray, ...]
        List of tidal heating coefficients [N m]
    ztorque_subterms: List[np.ndarray, ...]
        List of tidal torque coefficients [N m]
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
            orbital_frequency,
            7. * e2 - (101./4.) * e4 + (551./12.) * e6 + i2 - (13./12.) * i4 + (739./1440.) * i6 - (19./2.) * i2e2
            + (1003./32.) * i2e4 + (1097./192.) * i4e2,
            12. * e2 - (215./4.) * e4 + (8239./96.) * e6 + (1./4.) * i4 - (19./96.) * i6 - 16. * i2e2
            + (1831./32.) * i2e4 + (433./48.) * i4e2
        ),
        (
            '2n',
            2. * orbital_frequency,
            (605./16.) * e4 - (3847./24.) * e6 + (5./16.) * i4 - (5./24.) * i6 + (29./4.) * i2e2 - (1049./16.) * i2e4
            - (43./6.) * i4e2,
            (289./4.) * e4 - (1955./6.) * e6 - (1./4.) * i4 + (1./6.) * i6 + 5. * i2e2 - (1627./16.) * i2e4
            - (209./48.) * i4e2
        ),
        (
            '3n',
            3. * orbital_frequency,
            (2855./18.) * e6 + (1./32.) * i6 + (1237./32.) * i2e4 + (165./64.) * i4e2,
            (9917./32.) * e6 - (1./32.) * i6 + (1075./32.) * i2e4 - (9./16.) * i4e2
        )
    )

    heating_subterms = list()
    ztorque_subterms = list()
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
        heating_subterms.append(heating_coeff)
        ztorque_subterms.append(ztorque_coeff)

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return mode_names, modes, freqs, heating_subterms, ztorque_subterms

@njit
def nsr_modes(orbital_frequency: np.ndarray, spin_frequency: np.ndarray, eccentricity: np.ndarray,
              inclination: np.ndarray) -> ModeOutput:
    """ Non-Synchronous Rotation Tidal Modes (for heating and torque) for l=2, e^2, I^2

    These should all have the same shape!!

    Parameters
    ----------
    orbital_frequency : np.ndarray
        Orbital frequency (mean motion) [rads s-1]
    spin_frequency : np.ndarray
        Planet's spin frequency (inverse of period) [rads s-1]
    eccentricity : np.ndarray
        Orbital eccentricity
    inclination : np.ndarray
        Orbital inclination (relative to the orbital plane of the satellite) [rads]
        Note: this is not relative to the host planet's equator

    Returns
    -------
    mode_names : List[str, ...]
        List of easy to read names for each mode.
            "n" == orbital motion frequency
            "o" == spin frequency
    modes : List[np.ndarray, ...]
        List of calculated modes [rads s-1]
    freqs : List[np.ndarray, ...]
        List of calculated frequencies (absolute value of modes) [rads s-1]
    heating_subterms: List[np.ndarray, ...]
        List of tidal heating coefficients [N m]
    ztorque_subterms: List[np.ndarray, ...]
        List of tidal torque coefficients [N m]
    """

    e2 = eccentricity ** 2
    i2 = inclination ** 2

    modes_coeffs = (
        (
            'n',
            orbital_frequency,
            (3. / 4.) * e2,
            0.
        ),
        (
            'o',
            spin_frequency,
            (1. / 2.) * i2,
            -1.
        ),
        (
            '2n-o',
            2. * orbital_frequency - spin_frequency,
            (1. / 2.) * i2,
            1.
        ),
        (
            'n-2o',
            orbital_frequency - 2. * spin_frequency,
            (1. / 8.) * e2,
            2.
        ),
        (
            '2n-2o',
            2. * orbital_frequency - 2. * spin_frequency,
            (1. / 2.) - (5. / 2.) * e2 - (1. / 2.) * i2,
            2.
        ),
        (
            '3n-2o',
            3. * orbital_frequency - 2. * spin_frequency,
            (49. / 8.) * e2,
            2.
        ),
    )

    heating_subterms = list()
    ztorque_subterms = list()
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
        heating_subterms.append(heating_coeff)
        ztorque_subterms.append(ztorque_coeff)

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return mode_names, modes, freqs, heating_subterms, ztorque_subterms


@njit
def nsr_modes_4(orbital_frequency: np.ndarray, spin_frequency: np.ndarray, eccentricity: np.ndarray,
                inclination: np.ndarray) -> ModeOutput:
    """ Non-Synchronous Rotation Tidal Modes (for heating and torque) for l=2, e^4, I^4

    These should all have the same shape!!
    
    Parameters
    ----------
    orbital_frequency : np.ndarray
        Orbital frequency (mean motion) [rads s-1]
    spin_frequency : np.ndarray
        Planet's spin frequency (inverse of period) [rads s-1]
    eccentricity : np.ndarray
        Orbital eccentricity
    inclination : np.ndarray
        Orbital inclination (relative to the orbital plane of the satellite) [rads]
        Note: this is not relative to the host planet's equator

    Returns
    -------
    mode_names : List[str, ...]
        List of easy to read names for each mode.
            "n" == orbital motion frequency
            "o" == spin frequency
    modes : List[np.ndarray, ...]
        List of calculated modes [rads s-1]
    freqs : List[np.ndarray, ...]
        List of calculated frequencies (absolute value of modes) [rads s-1]
    heating_subterms: List[np.ndarray, ...]
        List of tidal heating coefficients [N m]
    ztorque_subterms: List[np.ndarray, ...]
        List of tidal torque coefficients [N m]
    """
    
    e2 = eccentricity ** 2
    e4 = eccentricity ** 4
    i2 = inclination ** 2
    i4 = inclination ** 4
    i2e2 = e2 * i2

    modes_coeffs = (
        (
            'n',
            orbital_frequency,
            (3. / 4.) * e2 + (27. / 16.) * e4 - (9. / 4.) * i2e2,
            0.
        ),
        (
            '2n',
            2. * orbital_frequency,
            (27. / 16.) * e4 + (3. / 16.) * i4,
            0.
        ),
        (
            'o',
            spin_frequency,
            (1. / 2.) * i2 - (2. / 3.) * i4 + (3. / 2.) * i2e2,
            -1.
        ),
        (
            '2o',
            2. * spin_frequency,
            (1. / 8.) * i4,
            -2.
        ),
        (
            'n-o',
            orbital_frequency - spin_frequency,
            (5. / 4.) * i2e2,
            1.
        ),
        (
            '2n-o',
            2. * orbital_frequency - spin_frequency,
            (1. / 2.) * i2 - (5. / 12.) * i4 - (5. / 2.) * i2e2,
            1.
        ),
        (
            '3n-o',
            3. * orbital_frequency - spin_frequency,
            (49. / 8.) * i2e2,
            1.
        ),
        (
            'n-2o',
            orbital_frequency - 2. * spin_frequency,
            (1. / 8.) * e2 - (1. / 32.) * e4 - (1. / 8.) * i2e2,
            2.
        ),
        (
            '2n-2o',
            2. * orbital_frequency - 2. * spin_frequency,
            (1. / 2.) - (5. / 2.) * e2 + (63. / 16.) * e4 - (1. / 2.) * i2 + (11. / 48.) * i4 + (5. / 2.) * i2e2,
            2.
        ),
        (
            '3n-2o',
            3. * orbital_frequency - 2. * spin_frequency,
            (49. / 8.) * e2 - (861. / 32.) * e4 - (49. / 8.) * i2e2,
            2.
        ),
        (
            '4n-2o',
            4. * orbital_frequency - 2. * spin_frequency,
            (289. / 8.) * e4,
            2.
        ),
        (
            '-n-o',
            -1. * orbital_frequency - spin_frequency,
            (9. / 8.) * i2e2,
            1.
        )
    )

    heating_subterms = list()
    ztorque_subterms = list()
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
        heating_subterms.append(heating_coeff)
        ztorque_subterms.append(ztorque_coeff)

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return mode_names, modes, freqs, heating_subterms, ztorque_subterms


@njit
def nsr_modes_6(orbital_frequency: np.ndarray, spin_frequency: np.ndarray, eccentricity: np.ndarray,
                inclination: np.ndarray) -> ModeOutput:
    """ Non-Synchronous Rotation Tidal Modes (for heating and torque) for l=2, e^6, I^6

    These should all have the same shape!!

    Parameters
    ----------
    orbital_frequency : np.ndarray
        Orbital frequency (mean motion) [rads s-1]
    spin_frequency : np.ndarray
        Planet's spin frequency (inverse of period) [rads s-1]
    eccentricity : np.ndarray
        Orbital eccentricity
    inclination : np.ndarray
        Orbital inclination (relative to the orbital plane of the satellite) [rads]
        Note: this is not relative to the host planet's equator

    Returns
    -------
    mode_names : List[str, ...]
        List of easy to read names for each mode.
            "n" == orbital motion frequency
            "o" == spin frequency
    modes : List[np.ndarray, ...]
        List of calculated modes [rads s-1]
    freqs : List[np.ndarray, ...]
        List of calculated frequencies (absolute value of modes) [rads s-1]
    heating_subterms: List[np.ndarray, ...]
        List of tidal heating coefficients [N m]
    ztorque_subterms: List[np.ndarray, ...]
        List of tidal torque coefficients [N m]
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
            orbital_frequency,
            (3./4.) * e2 + (27./16.) * e4 + (765./256.) * e6 - (9./4.) * i2e2 - (81./16.) * i2e4 + (159./64.) * i4e2,
            0.
        ),
        (
            '2n',
            2. * orbital_frequency,
            (27./16.) * e4 + (21./8.) * e6 - (81./16.) * i2e4 + (3./16.) * i4 - (15./16.) * i4e2 - (1./8.) * i6,
            0.
        ),
        (
            '3n',
            3. * orbital_frequency,
            (2809./768.) * e6 + (147./64.) * i4e2,
            0.
        ),
        (
            'o',
            spin_frequency,
            (1./2.) * i2 - (2./3.) * i4 + (16./45.) * i6 + (3./2.) * i2e2 + 3. * i2e4 - 2. * i4e2,
            -1.
        ),
        (
            '2o',
            2. * spin_frequency,
            (1./8.) * i4 - (1./12.) * i6 + (3./8.) * i4e2,
            -2.
         ),
        (
            'n-o',
            orbital_frequency - spin_frequency,
            (5./4.) * i2e2 + (5./2.) * i2e4 - (77./48.) * i4e2,
            1.
        ),
        (
            '2n-o',
            2. * orbital_frequency - spin_frequency,
            (1./2.) * i2 - (5./12.) * i4 + (227./1440.) * i6 - (5./2.) * i2e2 + (207./32.) * i2e4 + (25./12.) * i4e2,
            1.
        ),
        (
            '3n-o',
            3. * orbital_frequency - spin_frequency,
            (49./8.) * i2e2 - (861./32.) * i2e4 - (245./48.) * i4e2,
            1.
        ),
        (
            '4n-o',
            4.* orbital_frequency - spin_frequency,
            (289./8.) * i2e4,
            1.
        ),
        (
            'n-2o',
            orbital_frequency - 2. * spin_frequency,
            (1./8.) * e2 - (1./32.) * e4 + (13./1536.) * e6 - (1./8.) * i2e2 + (1./32.) * i2e4 + (65./192.) * i4e2,
            2.
        ),
        (
            '2n-2o',
            2. * orbital_frequency - 2. * spin_frequency,
            (1./2.) - (5./2.) * e2 + (63./16.) * e4 - (79./36.) * e6 - (1./2.) * i2 + (11./48.) * i4 - (23./360.) * i6
            + (5./2.) * i2e2 - (63./16.) * i2e4 - (55./48.) * i4e2,
            2.
        ),
        (
            '3n-2o',
            3. * orbital_frequency - 2. * spin_frequency,
            (49./8.) * e2 - (861./32.) * e4 + (21975./512.) * e6 - (49./8.) * i2e2 + (861./32.) * i2e4
            + (539./192.) * i4e2,
            2.
        ),
        (
            '4n-2o',
            4. * orbital_frequency - 2. * spin_frequency,
            (289./8.) * e4 - (1955./12.) * e6 - (289./8.) * i2e4,
            2.
        ),
        (
            '5n-2o',
            5. * orbital_frequency - 2. * spin_frequency,
            (714025./4608.) * e6,
            2.
        ),
        (
            '-n-o',
            -1. * orbital_frequency - spin_frequency,
            (9./8.) * i2e2 + (81./32.) * i2e4 - (3./2.) * i4e2,
            1.
        ),
        (
            '-n-2o',
            -1. * orbital_frequency - 2. * spin_frequency,
            (1./4608.) * e6 + (9./32.) * i4e2,
            2.
        ),
        (
            '-2n-o',
            -2. * orbital_frequency - 1. * spin_frequency,
            (1./32.) * i6 + (81./32.) * i2e4,
            1.
        )
    )

    heating_subterms = list()
    ztorque_subterms = list()
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
        heating_subterms.append(heating_coeff)
        ztorque_subterms.append(ztorque_coeff)

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return mode_names, modes, freqs, heating_subterms, ztorque_subterms