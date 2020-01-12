from typing import Tuple

import numpy as np

from ..types import NoneType
from ..performance import njit

OutputType = Tuple[Tuple[str, str, str, np.ndarray, NoneType, np.ndarray, np.ndarray, np.ndarray, np.ndarray], ...]

@njit
def sync_modes_t2(orbital_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray) \
        -> OutputType:
    """ Spin-sync Tidal Modes for truncation level 2, tidal order 2

    Parameters
    ----------
    orbital_freq : np.ndarray
        Planet's orbital mean motion [rad s-1]
    eccentricity : np.ndarray
        Planet's eccentricity
    inclination : np.ndarray
        Planet's inclination [rads]

    Returns
    -------
    mode_data_tuple : OutputType
        List of mode datas. Each mode data contains information to calculate the heating and torques.
    """
    # Tidal Potential Builder for l-order = 2, NSR = False.
    # Max Eccentricity Order = 2
    # Max Inclination Order = 2
    # Max q = 2.
    # Number of unique modes = 2.
    # Number of unique frequencies = 1.

    i2 = inclination**2
    e2 = eccentricity**2

    mode_data_output = (
        (
            '-n',
            'n',
            'None',
            -orbital_freq,
            None,
            (1./2.)*i2 + (1./2.)*e2,
            (-1./4.)*e2,
            (1./4.)*e2,
            (1./2.)*i2 + (1./4.)*e2
        ),
        (
            'n',
            'n',
            'None',
            orbital_freq,
            None,
            (1./2.)*i2 + (13./2.)*e2,
            i2 + (75./4.)*e2,
            i2 + (49./4.)*e2,
            (1./2.)*i2 + (49./4.)*e2
        )
    )

    return mode_data_output

@njit
def sync_modes_t4(orbital_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray) \
        -> OutputType:
    """ Spin-sync Tidal Modes for truncation level 4, tidal order 2

    Parameters
    ----------
    orbital_freq : np.ndarray
        Planet's orbital mean motion [rad s-1]
    eccentricity : np.ndarray
        Planet's eccentricity
    inclination : np.ndarray
        Planet's inclination [rads]

    Returns
    -------
    mode_data_tuple : OutputType
        List of mode datas. Each mode data contains information to calculate the heating and torques.
    """
    # Tidal Potential Builder for l-order = 2, NSR = False.
    # Max Eccentricity Order = 4
    # Max Inclination Order = 4
    # Max q = 3.
    # Number of unique modes = 4.
    # Number of unique frequencies = 2.

    i2 = inclination**2
    i4 = inclination**4
    e2 = eccentricity**2
    e4 = eccentricity**4
    i2e2 = i2*e2

    mode_data_output = (
        (
            '2*n',
            '2*n',
            'None',
            2*orbital_freq,
            None,
            (3./32.)*i4 + (49./8.)*i2e2 + (1183./32.)*e4,
            (3./16.)*i4 + (147./8.)*i2e2 + (2339./16.)*e4,
            (3./16.)*i4 + (49./4.)*i2e2 + (289./4.)*e4,
            (49./8.)*i2e2 + (289./4.)*e4
        ),
        (
            '-2*n',
            '2*n',
            'None',
            -2*orbital_freq,
            None,
            (7./32.)*i4 + (9./8.)*i2e2 + (27./32.)*e4,
            (-3./16.)*i4 + (-9./8.)*i2e2 + (-27./16.)*e4,
            (-3./16.)*i4,
            (1./4.)*i4 + (9./8.)*i2e2
        ),
        (
            '-n',
            'n',
            'None',
            -orbital_freq,
            None,
            (-2./3.)*i4 + (1./4.)*i2e2 + (1./2.)*i2 + (13./16.)*e4 + (1./2.)*e2,
            i2e2 + (-7./8.)*e4 + (-1./4.)*e2,
            (-1./4.)*i2e2 + (-1./16.)*e4 + (1./4.)*e2,
            (-2./3.)*i4 + (5./4.)*i2e2 + (1./2.)*i2 + (-1./16.)*e4 + (1./4.)*e2
        ),
        (
            'n',
            'n',
            'None',
            orbital_freq,
            None,
            (-5./12.)*i4 + (-39./4.)*i2e2 + (1./2.)*i2 + (-417./16.)*e4 + (13./2.)*e2,
            (-5./6.)*i4 + (-49./2.)*i2e2 + i2 + (-639./8.)*e4 + (75./4.)*e2,
            (-5./6.)*i4 + (-69./4.)*i2e2 + i2 + (-861./16.)*e4 + (49./4.)*e2,
            (-5./12.)*i4 + (-59./4.)*i2e2 + (1./2.)*i2 + (-861./16.)*e4 + (49./4.)*e2
        )
    )

    return mode_data_output

@njit
def sync_modes_t6(orbital_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray) \
        -> OutputType:
    """ Spin-sync Tidal Modes for truncation level 6, tidal order 2

    Parameters
    ----------
    orbital_freq : np.ndarray
        Planet's orbital mean motion [rad s-1]
    eccentricity : np.ndarray
        Planet's eccentricity
    inclination : np.ndarray
        Planet's inclination [rads]

    Returns
    -------
    mode_data_tuple : OutputType
        List of mode datas. Each mode data contains information to calculate the heating and torques.
    """
    # Tidal Potential Builder for l-order = 2, NSR = False.
    # Max Eccentricity Order = 6
    # Max Inclination Order = 6
    # Max q = 4.
    # Number of unique modes = 6.
    # Number of unique frequencies = 3.

    i2 = inclination**2
    i4 = inclination**4
    i6 = inclination**6
    e2 = eccentricity**2
    e4 = eccentricity**4
    e6 = eccentricity**6
    i2e2 = i2*e2
    i2e4 = i2*e4
    i4e2 = i4*e2

    mode_data_output = (
        (
            'n',
            'n',
            'None',
            orbital_freq,
            None,
            (227./1440.)*i6 + (785./128.)*i4e2 + (-5./12.)*i4 + (987./32.)*i2e4 + (-39./4.)*i2e2 + (1./2.)*i2 + (
                        5685./128.)*e6 + (-417./16.)*e4 + (13./2.)*e2,
            (227./720.)*i6 + (5311./384.)*i4e2 + (-5./6.)*i4 + (729./8.)*i2e4 + (-49./2.)*i2e2 + i2 + (
                        33345./256.)*e6 + (-639./8.)*e4 + (75./4.)*e2,
            (227./720.)*i6 + (629./64.)*i4e2 + (-5./6.)*i4 + (987./16.)*i2e4 + (-69./4.)*i2e2 + i2 + (
                        21975./256.)*e6 + (-861./16.)*e4 + (49./4.)*e2,
            (227./1440.)*i6 + (739./96.)*i4e2 + (-5./12.)*i4 + (1929./32.)*i2e4 + (-59./4.)*i2e2 + (1./2.)*i2 + (
                        21975./256.)*e6 + (-861./16.)*e4 + (49./4.)*e2
        ),
        (
            '2*n',
            '2*n',
            'None',
            2*orbital_freq,
            None,
            (-1./16.)*i6 + (-535./96.)*i4e2 + (3./32.)*i4 + (-1049./16.)*i2e4 + (49./8.)*i2e2 + (-7757./48.)*e6 + (
                        1183./32.)*e4,
            (-1./8.)*i6 + (-65./4.)*i4e2 + (3./16.)*i4 + (-7369./32.)*i2e4 + (147./8.)*i2e2 + (-15577./24.)*e6 + (
                        2339./16.)*e4,
            (-1./8.)*i6 + (-535./48.)*i4e2 + (3./16.)*i4 + (-2017./16.)*i2e4 + (49./4.)*i2e2 + (-1955./6.)*e6 + (
                        289./4.)*e4,
            (-245./48.)*i4e2 + (-3173./32.)*i2e4 + (49./8.)*i2e2 + (-1955./6.)*e6 + (289./4.)*e4
        ),
        (
            '3*n',
            '3*n',
            'None',
            3*orbital_freq,
            None,
            (147./128.)*i4e2 + (289./8.)*i2e4 + (180613./1152.)*e6,
            (441./128.)*i4e2 + (289./2.)*i2e4 + (1797703./2304.)*e6,
            (147./64.)*i4e2 + (289./4.)*i2e4 + (714025./2304.)*e6,
            (289./8.)*i2e4 + (714025./2304.)*e6
        ),
        (
            '-3*n',
            '3*n',
            'None',
            -3*orbital_freq,
            None,
            (1./32.)*i6 + (183./128.)*i4e2 + (81./32.)*i2e4 + (2107./1152.)*e6,
            (-1./16.)*i6 + (-477./128.)*i4e2 + (-81./16.)*i2e4 + (-12641./2304.)*e6,
            (-1./16.)*i6 + (-147./64.)*i4e2 + (1./2304.)*e6,
            (1./32.)*i6 + (9./16.)*i4e2 + (81./32.)*i2e4 + (1./2304.)*e6
        ),
        (
            '-2*n',
            '2*n',
            'None',
            -2*orbital_freq,
            None,
            (-7./48.)*i6 + (-51./32.)*i4e2 + (7./32.)*i4 + (9./8.)*i2e2 + (21./16.)*e6 + (27./32.)*e4,
            (1./8.)*i6 + (39./16.)*i4e2 + (-3./16.)*i4 + (81./32.)*i2e4 + (-9./8.)*i2e2 + (-21./8.)*e6 + (-27./16.)*e4,
            (1./8.)*i6 + (15./16.)*i4e2 + (-3./16.)*i4,
            (-1./6.)*i6 + (-3./4.)*i4e2 + (1./4.)*i4 + (81./32.)*i2e4 + (9./8.)*i2e2
        ),
        (
            '-n',
            'n',
            'None',
            -orbital_freq,
            None,
            (16./45.)*i6 + (-161./384.)*i4e2 + (-2./3.)*i4 + (1./2.)*i2e4 + (1./4.)*i2e2 + (1./2.)*i2 +
            (577./384.)*e6 + (13./16.)*e4 + (1./2.)*e2,
            (-347./384.)*i4e2 + (41./16.)*i2e4 + i2e2 + (-1141./768.)*e6 + (-7./8.)*e4 + (-1./4.)*e2,
            (13./192.)*i4e2 + (1./16.)*i2e4 + (-1./4.)*i2e2 + (13./768.)*e6 + (-1./16.)*e4 + (1./4.)*e2,
            (16./45.)*i6 + (-127./96.)*i4e2 + (-2./3.)*i4 + (49./16.)*i2e4 + (5./4.)*i2e2 + (1./2.)*i2 +
            (13./768.)*e6 + (-1./16.)*e4 + (1./4.)*e2
        )
    )

    return mode_data_output