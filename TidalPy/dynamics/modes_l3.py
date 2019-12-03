import numpy as np

from TidalPy.performance import njit
from . import MODE_ZERO_TOL
from .modes_l2 import ModeOutput

# TODO: A lot of the mode list construction is static and could take place outside of the functions at the expense of
#   readibility.

@njit
def spin_sync_modes(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray,
                    inclination: np.ndarray) -> ModeOutput:
    """ Synchronous Rotation Tidal Modes (for heating and torque) for l=3, e^2, I^2

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
    i2 = inclination**2

    modes_coeffs = (
        (
            'n',
            orbital_freq,
            (40./3.) * e2 + 2. * i2,
            32. * e2
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
def spin_sync_modes_4(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray,
                      inclination: np.ndarray) -> ModeOutput:
    """ Synchronous Rotation Tidal Modes (for heating and torque) for l=3, e^4, I^4

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
            orbital_freq,
            (40./3.) * e2 - (262./3.) * e4 + 2. * i2 - (53./12.) * i4 - 32. * i2e2,
            32. * e2 - 269. * e4 + (5./4.) * i4 - 76. * i2e2
        ),
        (
            '2n',
            2. * orbital_freq,
            (2795./24.) * e4 + (19./16) * i4 + 23. * i2e2,
            (651./2.) * e4 - (7./8.) * i4 + 30. * i2e2
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
def spin_sync_modes_6(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray,
                      inclination: np.ndarray) -> ModeOutput:
    """ Synchronous Rotation Tidal Modes (for heating and torque) for l=3, e^6, I^6

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
    i2e2 = i2 * e2
    i2e4 = i2 * e4
    i4e2 = i4 * e2

    modes_coeffs = (
        (
            'n',
            orbital_freq,
            (40./3.) * e2 - (262./3.) * e4 + (46241./144.) * e6 + 2. * i2 - (53./12.) * i4 + (12631./2880.) * i6
            - 32. * i2e2 + (11833./64.) * i2e4 + (1829./48.) * i4e2,
            32. * e2 - 269. * e4 + (7397./8.) * e6 + (5./4.) * i4 - (85./48.) * i6 - 76. * i2e2 + (15809./32.) * i2e4
            + (238./3.) * i4e2
        ),
        (
            '2n',
            2. * orbital_freq,
            (2795./24.) * e4 - (7511./9.) * e6 + (19./16) * i4 - (169./96.) * i6 + 23. * i2e2 - (10959./32.) * i2e4
            - (515./12.) * i4e2,
            (651./2.) * e4 - (10113./4.) * e6 - (7./8.) * i4 + (67./48.) * i6 + 30. * i2e2 - (1623./2.) * i2e4
            - (97./2.) * i4e2
        ),
        (
            '3n',
            3. * orbital_freq,
            (105695./144.) * e6 + (55./192.) * i6 + (12263./64.) * i2e4 + (241./16.) * i4e2,
            (50783./24.) * e6 - (5./16.) * i6 + (10005./32.) * i2e4 + 4. * i4e2
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
def nsr_modes(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray,
              inclination: np.ndarray) -> ModeOutput:
    """ Non-Synchronous Rotation Tidal Modes (for heating and torque) for l=3, e^2, I^2

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
    i2 = inclination**2

    modes_coeffs = (
        (
            'n',
            orbital_freq,
            (3./4.) * i2,
            0.
        ),
        (
            'o',
            spin_freq,
            (1./4.) * e2,
            -1.
        ),
        (
            'n-o',
            orbital_freq - spin_freq,
            (1./4.) + e2 - (11./8.) * i2,
            1.
        ),
        (
            '2n-o',
            2. * orbital_freq - spin_freq,
            (9./4.) * e2,
            1.
        ),
        (
            'n-2o',
            orbital_freq - 2. * spin_freq,
            (5./8.) * i2,
            2.
        ),
        (
            '3n-2o',
            3. * orbital_freq - 2. * spin_freq,
            (5./8.) * i2,
            2.
        ),
        (
            '2n-3o',
            2. * orbital_freq - 3. * spin_freq,
            (5./12.) * e2,
            3.
        ),
        (
            '3n-3o',
            3. * orbital_freq - 3. * spin_freq,
            (5./12.) - 5. * e2 - (5./8) * i2,
            3.
        ),
        (
            '4n-3o',
            4. * orbital_freq - 3. * spin_freq,
            (125./12.) * e2,
            3.
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
def nsr_modes_4(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray,
                inclination: np.ndarray) -> ModeOutput:
    """ Non-Synchronous Rotation Tidal Modes (for heating and torque) for l=3, e^4, I^4

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
    i2e2 = i2 * e2

    modes_coeffs = (
        (
            'n',
            orbital_freq,
            (3./4.) * i2 - (17./8.) * i4 + 3. * i2e2,
            0.
        ),
        (
            '2n',
            2. * orbital_freq,
            (27./4.) * i2e2,
            0.
        ),
        (
            'o',
            spin_freq,
            (1./4.) * e2 + (5./4.) * e4 - (11./8.) * i2e2,
            -1.
        ),
        (
            '2o',
            2. * spin_freq,
            (5./8.) * i2e2,
            -2.
         ),
        (
            'n-o',
            orbital_freq - spin_freq,
            (1./4.) + e2 + (367./128.) * e4 - (11./8.) * i2 + (535./192.) * i4 - (11./2.) * i2e2,
            1.
        ),
        (
            '2n-o',
            2. * orbital_freq - spin_freq,
            (9./4.) * e2 + (33./8.) * e4 - (99./8.) * i2e2,
            1.
        ),
        (
            '3n-o',
            3. * orbital_freq - spin_freq,
            (2809./256.) * e4 + (25./64.) * i4,
            1.
        ),
        (
            'n-2o',
            orbital_freq - 2. * spin_freq,
            (5./8.) * i2 - (35./24.) * i4 + (5./2.) * i2e2,
            2.
        ),
        (
            '2n-2o',
            2. * orbital_freq - 2. * spin_freq,
            (25./4) * i2e2,
            2.
        ),
        (
            '3n-2o',
            3. * orbital_freq - 2. * spin_freq,
            (5./8.) * i2 - (5./6.) * i4 - (15./2.) * i2e2,
            2.
        ),
        (
            '4n-2o',
            4. * orbital_freq - 2. * spin_freq,
            (125./8.) * i2e2,
            2.
        ),
        (
            'n-3o',
            orbital_freq - 3. * spin_freq,
            (5./768.) * e4 + (15./64.) * i4,
            3.
        ),
        (
            '2n-3o',
            2. * orbital_freq - 3. * spin_freq,
            (5./12.) * e2 - (25./24.) * e4 - (5./8.) * i2e2,
            3.
        ),
        (
            '3n-3o',
            3. * orbital_freq - 3. * spin_freq,
            (5./12.) - 5. * e2 + (2625./128.) * e4 - (5./8) * i2 + (85./192.) * i4 + (15./2.) * i2e2,
            3.
        ),
        (
            '4n-3o',
            4. * orbital_freq - 3. * spin_freq,
            (125./12.) * e2 - (275./3.) * e4 - (125./8.) * i2e2,
            3.
        ),
        (
            '5n-3o',
            5. * orbital_freq - 3. * spin_freq,
            (80645./768.) * e4,
            3.
        ),
        (
            '-n-o',
            -1. * orbital_freq - spin_freq,
            (121./256.) * e4 + (9./16.) * i4,
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
def nsr_modes_6(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray,
                inclination: np.ndarray) -> ModeOutput:
    """ Non-Synchronous Rotation Tidal Modes (for heating and torque) for l=3, e^6, I^6

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
    i2e2 = i2 * e2
    i2e4 = i2 * e4
    i4e2 = i4 * e2

    modes_coeffs = (
        (
            'n',
            orbital_freq,
            (3./4.) * i2 - (17./8.) * i4 + (2357./960.) * i6 + 3. * i2e2 + (2565./256.) * i2e4 - (17./2.) * i4e2,
            0.
        ),
        (
            '2n',
            2. * orbital_freq,
            (27./4.) * i2e2 + (99./8.) * i2e4 - (153./8.) * i4e2,
            0.
        ),
        (
            '3n',
            3. * orbital_freq,
            (25./192.) * i6 + (8427./256.) * i2e4,
            0.
        ),
        (
            'o',
            spin_freq,
            (1./4.) * e2 + (5./4.) * e4 + (15./4.) * e6 - (11./8.) * i2e2 - (55./8.) * i2e4 + (643./192.) * i4e2,
            -1.
        ),
        (
            '2o',
            2. * spin_freq,
            (5./8.) * i2e2 + (25./8.)*i2e4 - (35./24.) * i4e2,
            -2.
         ),
        (
            '3o',
            3. * spin_freq,
            (15./64.) * i4e2,
            -3.
        ),
        (
            'n-o',
            orbital_freq - spin_freq,
            (1./4.) + e2 + (367./128.) * e4 + (7613./1152.) * e6 - (11./8.) * i2 + (535./192.) * i4
            - (15637./5760.) * i6 - (11./2.) * i2e2 - (4037./256.) * i2e4 + (535./48.) * i4e2,
            1.
        ),
        (
            '2n-o',
            2. * orbital_freq - spin_freq,
            (9./4.) * e2 + (33./8.) * e4 + (611./64.) * e6 - (99./8.) * i2e2 - (363./16.) * i2e4 + (815./32.) * i4e2,
            1.
        ),
        (
            '3n-o',
            3. * orbital_freq - spin_freq,
            (2809./256.) * e4 + (2067./256.) * e6 + (25./64.) * i4 - (175./384.) * i6 - (30899./512.) * i2e4
            - (75./16.) * i4e2,
            1.
        ),
        (
            '4n-o',
            4.* orbital_freq - spin_freq,
            (5929./144.) * e6 + (625./64.) * i4e2,
            1.
        ),
        (
            'n-2o',
            orbital_freq - 2. * spin_freq,
            (5./8.) * i2 - (35./24.) * i4 + (811./576.) * i6 + (5./2.) * i2e2 + (3675./512.) * i2e4 - (35./6.) * i4e2,
            2.
        ),
        (
            '2n-2o',
            2. * orbital_freq - 2. * spin_freq,
            (25./4) * i2e2 + (35./4.) * i2e4 - (335./24.) * i4e2,
            2.
        ),
        (
            '3n-2o',
            3. * orbital_freq - 2. * spin_freq,
            (5./8.) * i2 - (5./6.) * i4 + (301./576.) * i6 - (15./2.) * i2e2 + (29795./512.) * i2e4 + 10. * i4e2,
            2.
        ),
        (
            '4n-2o',
            4. * orbital_freq - 2. * spin_freq,
            (125./8.) * i2e2 - (275./2.) * i2e4 - (125./6.) * i4e2,
            2.
        ),
        (
            '5n-2o',
            5. * orbital_freq - 2. * spin_freq,
            (80645./512.) * i2e4,
            2.
        ),
        (
            'n-3o',
            orbital_freq - 3. * spin_freq,
            (5./768.) * e4 + (5./2304.) * e6 + (15./64.) * i4 - (35./128.) * i6 - (5./512.) * i2e4 + (15./16.) * i4e2,
            3.
        ),
        (
            '2n-3o',
            2. * orbital_freq - 3. * spin_freq,
            (5./12.) * e2 - (25./24.) * e4 + (445./576.) * e6 - (5./8.) * i2e2 + (25./16.) * i2e4 + (245./96.) *i4e2,
            3.
        ),
        (
            '3n-3o',
            3. * orbital_freq - 3. * spin_freq,
            (5./12.) - 5. * e2 + (2625./128.) * e4 - (4445./128.) * e6 - (5./8) * i2 + (85./192.) * i4
            - (227./1152.) * i6 + (15./2.) * i2e2 - (7875./256.) * i2e4 - (85./16.) * i4e2,
            3.
        ),
        (
            '4n-3o',
            4. * orbital_freq - 3. * spin_freq,
            (125./12.) * e2 - (275./3.) * e4 + (44215./144.) * e6 - (125./8.) * i2e2 + (275./2.) * i2e4
            + (2125./192.) * i4e2,
            3.
        ),
        (
            '5n-3o',
            5. * orbital_freq - 3. * spin_freq,
            (80645./768.) * e4 - (1946275./2304.) * e6 - (80645./512.) * i2e4,
            3.
        ),
        (
            '6n-3o',
            6. * orbital_freq - 3. * spin_freq,
            (132845./192.) * e6,
            3.
        ),
        (
            '-n-o',
            -1. * orbital_freq - spin_freq,
            (121./256.) * e4 + (539./256.) * e6 + (9./16.) * i4 - (33./32.) * i6 - (1331./512.) * i2e4 + (9./4.) * i4e2,
            1.
        ),
        (
            '-n-2o',
            -1. * orbital_freq - 2. * spin_freq,
            (605./512.) * i2e4 + (5./32.) * i6,
            2.
        ),
        (
            '-2n-o',
            -2. * orbital_freq - 1. * spin_freq,
            (529./576.) * e6 + (81./16.) * i4e2,
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