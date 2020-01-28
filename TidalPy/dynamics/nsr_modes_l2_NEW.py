from typing import Tuple

import numpy as np

from ..types import NoneType
from ..performance import njit

OutputType = Tuple[Tuple[str, str, str,  np.ndarray, np.ndarray, NoneType, float, float, float], ...]

@njit
def nsr_modes_t2(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray) \
        -> OutputType:
    """ Non-sync Tidal Modes for truncation level 2, tidal order 2

    Parameters
    ----------
    orbital_freq : np.ndarray
        Planet's orbital mean motion [rad s-1]
    spin_freq : np.ndarray
        Planet's spin frequency [rad s-1]
    eccentricity : np.ndarray
        Planet's eccentricity
    inclination : np.ndarray
        Planet's inclination [rads]

    Returns
    -------
    mode_data_tuple : OutputType
        List of mode datas. Each mode data contains information to calculate the heating and torques.
    """
    # Tidal Potential Builder for l-order = 2, NSR = True.
    # Max Eccentricity Order = 2
    # Max Inclination Order = 2
    # Max q = 2.
    # Number of unique modes = 7.
    # Number of unique frequencies = 6.

    i2 = inclination**2
    e2 = eccentricity**2

    mode_data_output = (
        (
            '-n',
            'n',
            '2, 0, 1, -1',
            -orbital_freq,
            (3./8.)*e2,
            None,
            -1.,
            0.,
            0.
        ),
        (
            'n',
            'n',
            '2, 0, 1, 1',
            orbital_freq,
            (3./8.)*e2,
            None,
            1.,
            0.,
            0.
        ),
        (
            '-O + 2*n',
            'O - 2*n',
            '2, 1, 0, 0',
            -spin_freq + 2*orbital_freq,
            (1./2.)*i2,
            None,
            2.,
            2.,
            1.
        ),
        (
            '-O',
            'O',
            '2, 1, 1, 0',
            -spin_freq,
            (1./2.)*i2,
            None,
            0.,
            0.,
            1.
        ),
        (
            '-2*O + n',
            '2*O - n',
            '2, 2, 0, -1',
            -2*spin_freq + orbital_freq,
            (1./8.)*e2,
            None,
            1.,
            2.,
            2.
        ),
        (
            '-2*O + 2*n',
            '2*O - 2*n',
            '2, 2, 0, 0',
            -2*spin_freq + 2*orbital_freq,
            (-1./2.)*i2 + (-5./2.)*e2 + (1./2.),
            None,
            2.,
            2.,
            2.
        ),
        (
            '-2*O + 3*n',
            '2*O - 3*n',
            '2, 2, 0, 1',
            -2*spin_freq + 3*orbital_freq,
            (49./8.)*e2,
            None,
            3.,
            2.,
            2.
        )
    )

    return mode_data_output

@njit
def nsr_modes_t4(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray) \
        -> OutputType:
    """ Non-sync Tidal Modes for truncation level 4, tidal order 2

    Parameters
    ----------
    orbital_freq : np.ndarray
        Planet's orbital mean motion [rad s-1]
    spin_freq : np.ndarray
        Planet's spin frequency [rad s-1]
    eccentricity : np.ndarray
        Planet's eccentricity
    inclination : np.ndarray
        Planet's inclination [rads]

    Returns
    -------
    mode_data_tuple : OutputType
        List of mode datas. Each mode data contains information to calculate the heating and torques.
    """
    # Tidal Potential Builder for l-order = 2, NSR = True.
    # Max Eccentricity Order = 4
    # Max Inclination Order = 4
    # Max q = 3.
    # Number of unique modes = 14.
    # Number of unique frequencies = 12.

    i2 = inclination**2
    i4 = inclination**4
    e2 = eccentricity**2
    e4 = eccentricity**4
    i2e2 = i2*e2

    mode_data_output = (
        (
            '2*n',
            '2*n',
            '2, 0, 0, 0',
            2*orbital_freq,
            (3./32.)*i4,
            None,
            2.,
            2.,
            0.
        ),
        (
            '-2*n',
            '2*n',
            '2, 0, 1, -2',
            -2*orbital_freq,
            (27./32.)*e4,
            None,
            -2.,
            0.,
            0.
        ),
        (
            '-n',
            'n',
            '2, 0, 1, -1',
            -orbital_freq,
            (-9./8.)*i2e2 + (27./32.)*e4 + (3./8.)*e2,
            None,
            -1.,
            0.,
            0.
        ),
        (
            'n',
            'n',
            '2, 0, 1, 1',
            orbital_freq,
            (-9./8.)*i2e2 + (27./32.)*e4 + (3./8.)*e2,
            None,
            1.,
            0.,
            0.
        ),
        (
            '2*n',
            '2*n',
            '2, 0, 1, 2',
            2*orbital_freq,
            (27./32.)*e4,
            None,
            2.,
            0.,
            0.
        ),
        (
            '-2*n',
            '2*n',
            '2, 0, 2, 0',
            -2*orbital_freq,
            (3./32.)*i4,
            None,
            -2.,
            -2.,
            0.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 0, -1',
            -spin_freq + orbital_freq,
            (1./8.)*i2e2,
            None,
            1.,
            2.,
            1.
        ),
        (
            '-O + 2*n',
            'O - 2*n',
            '2, 1, 0, 0',
            -spin_freq + 2*orbital_freq,
            (-5./12.)*i4 + (-5./2.)*i2e2 + (1./2.)*i2,
            None,
            2.,
            2.,
            1.
        ),
        (
            '-O + 3*n',
            'O - 3*n',
            '2, 1, 0, 1',
            -spin_freq + 3*orbital_freq,
            (49./8.)*i2e2,
            None,
            3.,
            2.,
            1.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 1, -1',
            -spin_freq - orbital_freq,
            (9./8.)*i2e2,
            None,
            -1.,
            0.,
            1.
        ),
        (
            '-O',
            'O',
            '2, 1, 1, 0',
            -spin_freq,
            (-2./3.)*i4 + (3./2.)*i2e2 + (1./2.)*i2,
            None,
            0.,
            0.,
            1.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 1, 1',
            -spin_freq + orbital_freq,
            (9./8.)*i2e2,
            None,
            1.,
            0.,
            1.
        ),
        (
            '-2*O + n',
            '2*O - n',
            '2, 2, 0, -1',
            -2*spin_freq + orbital_freq,
            (-1./8.)*i2e2 + (-1./32.)*e4 + (1./8.)*e2,
            None,
            1.,
            2.,
            2.
        ),
        (
            '-2*O + 2*n',
            '2*O - 2*n',
            '2, 2, 0, 0',
            -2*spin_freq + 2*orbital_freq,
            (11./48.)*i4 + (5./2.)*i2e2 + (-1./2.)*i2 + (63./16.)*e4 + (-5./2.)*e2 + (1./2.),
            None,
            2.,
            2.,
            2.
        ),
        (
            '-2*O + 3*n',
            '2*O - 3*n',
            '2, 2, 0, 1',
            -2*spin_freq + 3*orbital_freq,
            (-49./8.)*i2e2 + (-861./32.)*e4 + (49./8.)*e2,
            None,
            3.,
            2.,
            2.
        ),
        (
            '-2*O + 4*n',
            '2*O - 4*n',
            '2, 2, 0, 2',
            -2*spin_freq + 4*orbital_freq,
            (289./8.)*e4,
            None,
            4.,
            2.,
            2.
        ),
        (
            '-2*O',
            '2*O',
            '2, 2, 1, 0',
            -2*spin_freq,
            (1./8.)*i4,
            None,
            0.,
            0.,
            2.
        )
    )

    return mode_data_output

@njit
def nsr_modes_t6(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray) \
        -> OutputType:
    """ Non-sync Tidal Modes for truncation level 6, tidal order 2

    Parameters
    ----------
    orbital_freq : np.ndarray
        Planet's orbital mean motion [rad s-1]
    spin_freq : np.ndarray
        Planet's spin frequency [rad s-1]
    eccentricity : np.ndarray
        Planet's eccentricity
    inclination : np.ndarray
        Planet's inclination [rads]

    Returns
    -------
    mode_data_tuple : OutputType
        List of mode datas. Each mode data contains information to calculate the heating and torques.
    """
    # Tidal Potential Builder for l-order = 2, NSR = True.
    # Max Eccentricity Order = 6
    # Max Inclination Order = 6
    # Max q = 4.
    # Number of unique modes = 20.
    # Number of unique frequencies = 17.

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
            '2, 0, 0, -1',
            orbital_freq,
            (3./128.)*i4e2,
            None,
            1.,
            2.,
            0.
        ),
        (
            '2*n',
            '2*n',
            '2, 0, 0, 0',
            2*orbital_freq,
            (-1./16.)*i6 + (-15./32.)*i4e2 + (3./32.)*i4,
            None,
            2.,
            2.,
            0.
        ),
        (
            '3*n',
            '3*n',
            '2, 0, 0, 1',
            3*orbital_freq,
            (147./128.)*i4e2,
            None,
            3.,
            2.,
            0.
        ),
        (
            '-3*n',
            '3*n',
            '2, 0, 1, -3',
            -3*orbital_freq,
            (2809./1536.)*e6,
            None,
            -3.,
            0.,
            0.
        ),
        (
            '-2*n',
            '2*n',
            '2, 0, 1, -2',
            -2*orbital_freq,
            (-81./32.)*i2e4 + (21./16.)*e6 + (27./32.)*e4,
            None,
            -2.,
            0.,
            0.
        ),
        (
            '-n',
            'n',
            '2, 0, 1, -1',
            -orbital_freq,
            (39./32.)*i4e2 + (-81./32.)*i2e4 + (-9./8.)*i2e2 + (765./512.)*e6 + (27./32.)*e4 + (3./8.)*e2,
            None,
            -1.,
            0.,
            0.
        ),
        (
            'n',
            'n',
            '2, 0, 1, 1',
            orbital_freq,
            (39./32.)*i4e2 + (-81./32.)*i2e4 + (-9./8.)*i2e2 + (765./512.)*e6 + (27./32.)*e4 + (3./8.)*e2,
            None,
            1.,
            0.,
            0.
        ),
        (
            '2*n',
            '2*n',
            '2, 0, 1, 2',
            2*orbital_freq,
            (-81./32.)*i2e4 + (21./16.)*e6 + (27./32.)*e4,
            None,
            2.,
            0.,
            0.
        ),
        (
            '3*n',
            '3*n',
            '2, 0, 1, 3',
            3*orbital_freq,
            (2809./1536.)*e6,
            None,
            3.,
            0.,
            0.
        ),
        (
            '-3*n',
            '3*n',
            '2, 0, 2, -1',
            -3*orbital_freq,
            (147./128.)*i4e2,
            None,
            -3.,
            -2.,
            0.
        ),
        (
            '-2*n',
            '2*n',
            '2, 0, 2, 0',
            -2*orbital_freq,
            (-1./16.)*i6 + (-15./32.)*i4e2 + (3./32.)*i4,
            None,
            -2.,
            -2.,
            0.
        ),
        (
            '-n',
            'n',
            '2, 0, 2, 1',
            -orbital_freq,
            (3./128.)*i4e2,
            None,
            -1.,
            -2.,
            0.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 0, -1',
            -spin_freq + orbital_freq,
            (-5./48.)*i4e2 + (-1./32.)*i2e4 + (1./8.)*i2e2,
            None,
            1.,
            2.,
            1.
        ),
        (
            '-O + 2*n',
            'O - 2*n',
            '2, 1, 0, 0',
            -spin_freq + 2*orbital_freq,
            (227./1440.)*i6 + (25./12.)*i4e2 + (-5./12.)*i4 + (63./16.)*i2e4 + (-5./2.)*i2e2 + (1./2.)*i2,
            None,
            2.,
            2.,
            1.
        ),
        (
            '-O + 3*n',
            'O - 3*n',
            '2, 1, 0, 1',
            -spin_freq + 3*orbital_freq,
            (-245./48.)*i4e2 + (-861./32.)*i2e4 + (49./8.)*i2e2,
            None,
            3.,
            2.,
            1.
        ),
        (
            '-O + 4*n',
            'O - 4*n',
            '2, 1, 0, 2',
            -spin_freq + 4*orbital_freq,
            (289./8.)*i2e4,
            None,
            4.,
            2.,
            1.
        ),
        (
            '-O - 2*n',
            'O + 2*n',
            '2, 1, 1, -2',
            -spin_freq - 2*orbital_freq,
            (81./32.)*i2e4,
            None,
            -2.,
            0.,
            1.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 1, -1',
            -spin_freq - orbital_freq,
            (-3./2.)*i4e2 + (81./32.)*i2e4 + (9./8.)*i2e2,
            None,
            -1.,
            0.,
            1.
        ),
        (
            '-O',
            'O',
            '2, 1, 1, 0',
            -spin_freq,
            (16./45.)*i6 + -2.*i4e2 + (-2./3.)*i4 + 3.*i2e4 + (3./2.)*i2e2 + (1./2.)*i2,
            None,
            0.,
            0.,
            1.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 1, 1',
            -spin_freq + orbital_freq,
            (-3./2.)*i4e2 + (81./32.)*i2e4 + (9./8.)*i2e2,
            None,
            1.,
            0.,
            1.
        ),
        (
            '-O + 2*n',
            'O - 2*n',
            '2, 1, 1, 2',
            -spin_freq + 2*orbital_freq,
            (81./32.)*i2e4,
            None,
            2.,
            0.,
            1.
        ),
        (
            '-O - 2*n',
            'O + 2*n',
            '2, 1, 2, 0',
            -spin_freq - 2*orbital_freq,
            (1./32.)*i6,
            None,
            -2.,
            -2.,
            1.
        ),
        (
            '-2*O - n',
            '2*O + n',
            '2, 2, 0, -3',
            -2*spin_freq - orbital_freq,
            (1./4608.)*e6,
            None,
            -1.,
            2.,
            2.
        ),
        (
            '-2*O + n',
            '2*O - n',
            '2, 2, 0, -1',
            -2*spin_freq + orbital_freq,
            (11./192.)*i4e2 + (1./32.)*i2e4 + (-1./8.)*i2e2 + (13./1536.)*e6 + (-1./32.)*e4 + (1./8.)*e2,
            None,
            1.,
            2.,
            2.
        ),
        (
            '-2*O + 2*n',
            '2*O - 2*n',
            '2, 2, 0, 0',
            -2*spin_freq + 2*orbital_freq,
            (-23./360.)*i6 + (-55./48.)*i4e2 + (11./48.)*i4 + (-63./16.)*i2e4 + (5./2.)*i2e2 + (-1./2.)*i2 + (
                        -79./36.)*e6 + (63./16.)*e4 + (-5./2.)*e2 + (1./2.),
            None,
            2.,
            2.,
            2.
        ),
        (
            '-2*O + 3*n',
            '2*O - 3*n',
            '2, 2, 0, 1',
            -2*spin_freq + 3*orbital_freq,
            (539./192.)*i4e2 + (861./32.)*i2e4 + (-49./8.)*i2e2 + (21975./512.)*e6 + (-861./32.)*e4 + (49./8.)*e2,
            None,
            3.,
            2.,
            2.
        ),
        (
            '-2*O + 4*n',
            '2*O - 4*n',
            '2, 2, 0, 2',
            -2*spin_freq + 4*orbital_freq,
            (-289./8.)*i2e4 + (-1955./12.)*e6 + (289./8.)*e4,
            None,
            4.,
            2.,
            2.
        ),
        (
            '-2*O + 5*n',
            '2*O - 5*n',
            '2, 2, 0, 3',
            -2*spin_freq + 5*orbital_freq,
            (714025./4608.)*e6,
            None,
            5.,
            2.,
            2.
        ),
        (
            '-2*O - n',
            '2*O + n',
            '2, 2, 1, -1',
            -2*spin_freq - orbital_freq,
            (9./32.)*i4e2,
            None,
            -1.,
            0.,
            2.
        ),
        (
            '-2*O',
            '2*O',
            '2, 2, 1, 0',
            -2*spin_freq,
            (-1./12.)*i6 + (3./8.)*i4e2 + (1./8.)*i4,
            None,
            0.,
            0.,
            2.
        ),
        (
            '-2*O + n',
            '2*O - n',
            '2, 2, 1, 1',
            -2*spin_freq + orbital_freq,
            (9./32.)*i4e2,
            None,
            1.,
            0.,
            2.
        )
    )

    return mode_data_output

@njit
def nsr_modes_t8(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray) \
        -> OutputType:
    """ Non-sync Tidal Modes for truncation level 8, tidal order 2

    Parameters
    ----------
    orbital_freq : np.ndarray
        Planet's orbital mean motion [rad s-1]
    spin_freq : np.ndarray
        Planet's spin frequency [rad s-1]
    eccentricity : np.ndarray
        Planet's eccentricity
    inclination : np.ndarray
        Planet's inclination [rads]

    Returns
    -------
    mode_data_tuple : OutputType
        List of mode datas. Each mode data contains information to calculate the heating and torques.
    """
    # Tidal Potential Builder for l-order = 2, NSR = True.
    # Max Eccentricity Order = 8
    # Max Inclination Order = 8
    # Max q = 5.
    # Number of unique modes = 26.
    # Number of unique frequencies = 22.

    i2 = inclination**2
    i4 = inclination**4
    i6 = inclination**6
    i8 = inclination**8
    e2 = eccentricity**2
    e4 = eccentricity**4
    e6 = eccentricity**6
    e8 = eccentricity**8
    i2e2 = i2 * e2
    i2e4 = i2 * e4
    i2e6 = i2 * e6
    i4e2 = i4 * e2
    i4e4 = i4 * e4
    i6e2 = i6 * e2

    mode_data_output = (
        (
            'n',
            'n',
            '2, 0, 0, -1',
            orbital_freq,
            (-1./64.)*i6e2 + (-3./512.)*i4e4 + (3./128.)*i4e2,
            None,
            1.,
            2.,
            0.
        ),
        (
            '2*n',
            '2*n',
            '2, 0, 0, 0',
            2*orbital_freq,
            (3./160.)*i8 + (5./16.)*i6e2 + (-1./16.)*i6 + (189./256.)*i4e4 + (-15./32.)*i4e2 + (3./32.)*i4,
            None,
            2.,
            2.,
            0.
        ),
        (
            '3*n',
            '3*n',
            '2, 0, 0, 1',
            3*orbital_freq,
            (-49./64.)*i6e2 + (-2583./512.)*i4e4 + (147./128.)*i4e2,
            None,
            3.,
            2.,
            0.
        ),
        (
            '4*n',
            '4*n',
            '2, 0, 0, 2',
            4*orbital_freq,
            (867./128.)*i4e4,
            None,
            4.,
            2.,
            0.
        ),
        (
            '-4*n',
            '4*n',
            '2, 0, 1, -4',
            -4*orbital_freq,
            (5929./1536.)*e8,
            None,
            -4.,
            0.,
            0.
        ),
        (
            '-3*n',
            '3*n',
            '2, 0, 1, -3',
            -3*orbital_freq,
            (-2809./512.)*i2e6 + (6943./4096.)*e8 + (2809./1536.)*e6,
            None,
            -3.,
            0.,
            0.
        ),
        (
            '-2*n',
            '2*n',
            '2, 0, 1, -2',
            -2*orbital_freq,
            (351./128.)*i4e4 + (-63./16.)*i2e6 + (-81./32.)*i2e4 + (1661./768.)*e8 + (21./16.)*e6 + (27./32.)*e4,
            None,
            -2.,
            0.,
            0.
        ),
        (
            '-n',
            'n',
            '2, 0, 1, -1',
            -orbital_freq,
            (-49./80.)*i6e2 + (351./128.)*i4e4 + (39./32.)*i4e2 + (-2295./512.)*i2e6 + (-81./32.)*i2e4 + (-9./8.)*i2e2 + (28211./12288.)*e8 + (765./512.)*e6 + (27./32.)*e4 + (3./8.)*e2,
            None,
            -1.,
            0.,
            0.
        ),
        (
            'n',
            'n',
            '2, 0, 1, 1',
            orbital_freq,
            (-49./80.)*i6e2 + (351./128.)*i4e4 + (39./32.)*i4e2 + (-2295./512.)*i2e6 + (-81./32.)*i2e4 + (-9./8.)*i2e2 + (28211./12288.)*e8 + (765./512.)*e6 + (27./32.)*e4 + (3./8.)*e2,
            None,
            1.,
            0.,
            0.
        ),
        (
            '2*n',
            '2*n',
            '2, 0, 1, 2',
            2*orbital_freq,
            (351./128.)*i4e4 + (-63./16.)*i2e6 + (-81./32.)*i2e4 + (1661./768.)*e8 + (21./16.)*e6 + (27./32.)*e4,
            None,
            2.,
            0.,
            0.
        ),
        (
            '3*n',
            '3*n',
            '2, 0, 1, 3',
            3*orbital_freq,
            (-2809./512.)*i2e6 + (6943./4096.)*e8 + (2809./1536.)*e6,
            None,
            3.,
            0.,
            0.
        ),
        (
            '4*n',
            '4*n',
            '2, 0, 1, 4',
            4*orbital_freq,
            (5929./1536.)*e8,
            None,
            4.,
            0.,
            0.
        ),
        (
            '-4*n',
            '4*n',
            '2, 0, 2, -2',
            -4*orbital_freq,
            (867./128.)*i4e4,
            None,
            -4.,
            -2.,
            0.
        ),
        (
            '-3*n',
            '3*n',
            '2, 0, 2, -1',
            -3*orbital_freq,
            (-49./64.)*i6e2 + (-2583./512.)*i4e4 + (147./128.)*i4e2,
            None,
            -3.,
            -2.,
            0.
        ),
        (
            '-2*n',
            '2*n',
            '2, 0, 2, 0',
            -2*orbital_freq,
            (3./160.)*i8 + (5./16.)*i6e2 + (-1./16.)*i6 + (189./256.)*i4e4 + (-15./32.)*i4e2 + (3./32.)*i4,
            None,
            -2.,
            -2.,
            0.
        ),
        (
            '-n',
            'n',
            '2, 0, 2, 1',
            -orbital_freq,
            (-1./64.)*i6e2 + (-3./512.)*i4e4 + (3./128.)*i4e2,
            None,
            -1.,
            -2.,
            0.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 0, -3',
            -spin_freq - orbital_freq,
            (1./4608.)*i2e6,
            None,
            -1.,
            2.,
            1.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 0, -1',
            -spin_freq + orbital_freq,
            (227./5760.)*i6e2 + (5./192.)*i4e4 + (-5./48.)*i4e2 + (13./1536.)*i2e6 + (-1./32.)*i2e4 + (1./8.)*i2e2,
            None,
            1.,
            2.,
            1.
        ),
        (
            '-O + 2*n',
            'O - 2*n',
            '2, 1, 0, 0',
            -spin_freq + 2*orbital_freq,
            (-145./4032.)*i8 + (-227./288.)*i6e2 + (227./1440.)*i6 + (-105./32.)*i4e4 + (25./12.)*i4e2 + (-5./12.)*i4 + (-79./36.)*i2e6 + (63./16.)*i2e4 + (-5./2.)*i2e2 + (1./2.)*i2,
            None,
            2.,
            2.,
            1.
        ),
        (
            '-O + 3*n',
            'O - 3*n',
            '2, 1, 0, 1',
            -spin_freq + 3*orbital_freq,
            (11123./5760.)*i6e2 + (1435./64.)*i4e4 + (-245./48.)*i4e2 + (21975./512.)*i2e6 + (-861./32.)*i2e4 + (49./8.)*i2e2,
            None,
            3.,
            2.,
            1.
        ),
        (
            '-O + 4*n',
            'O - 4*n',
            '2, 1, 0, 2',
            -spin_freq + 4*orbital_freq,
            (-1445./48.)*i4e4 + (-1955./12.)*i2e6 + (289./8.)*i2e4,
            None,
            4.,
            2.,
            1.
        ),
        (
            '-O + 5*n',
            'O - 5*n',
            '2, 1, 0, 3',
            -spin_freq + 5*orbital_freq,
            (714025./4608.)*i2e6,
            None,
            5.,
            2.,
            1.
        ),
        (
            '-O - 3*n',
            'O + 3*n',
            '2, 1, 1, -3',
            -spin_freq - 3*orbital_freq,
            (2809./512.)*i2e6,
            None,
            -3.,
            0.,
            1.
        ),
        (
            '-O - 2*n',
            'O + 2*n',
            '2, 1, 1, -2',
            -spin_freq - 2*orbital_freq,
            (-27./8.)*i4e4 + (63./16.)*i2e6 + (81./32.)*i2e4,
            None,
            -2.,
            0.,
            1.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 1, -1',
            -spin_freq - orbital_freq,
            (4./5.)*i6e2 + (-27./8.)*i4e4 + (-3./2.)*i4e2 + (2295./512.)*i2e6 + (81./32.)*i2e4 + (9./8.)*i2e2,
            None,
            -1.,
            0.,
            1.
        ),
        (
            '-O',
            'O',
            '2, 1, 1, 0',
            -spin_freq,
            (-32./315.)*i8 + (16./15.)*i6e2 + (16./45.)*i6 + -4.*i4e4 + -2.*i4e2 + (-2./3.)*i4 + 5.*i2e6 + 3.*i2e4 + (3./2.)*i2e2 + (1./2.)*i2,
            None,
            0.,
            0.,
            1.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 1, 1',
            -spin_freq + orbital_freq,
            (4./5.)*i6e2 + (-27./8.)*i4e4 + (-3./2.)*i4e2 + (2295./512.)*i2e6 + (81./32.)*i2e4 + (9./8.)*i2e2,
            None,
            1.,
            0.,
            1.
        ),
        (
            '-O + 2*n',
            'O - 2*n',
            '2, 1, 1, 2',
            -spin_freq + 2*orbital_freq,
            (-27./8.)*i4e4 + (63./16.)*i2e6 + (81./32.)*i2e4,
            None,
            2.,
            0.,
            1.
        ),
        (
            '-O + 3*n',
            'O - 3*n',
            '2, 1, 1, 3',
            -spin_freq + 3*orbital_freq,
            (2809./512.)*i2e6,
            None,
            3.,
            0.,
            1.
        ),
        (
            '-O - 3*n',
            'O + 3*n',
            '2, 1, 2, -1',
            -spin_freq - 3*orbital_freq,
            (49./128.)*i6e2,
            None,
            -3.,
            -2.,
            1.
        ),
        (
            '-O - 2*n',
            'O + 2*n',
            '2, 1, 2, 0',
            -spin_freq - 2*orbital_freq,
            (-1./64.)*i8 + (-5./32.)*i6e2 + (1./32.)*i6,
            None,
            -2.,
            -2.,
            1.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 2, 1',
            -spin_freq - orbital_freq,
            (1./128.)*i6e2,
            None,
            -1.,
            -2.,
            1.
        ),
        (
            '-2*O - 2*n',
            '2*O + 2*n',
            '2, 2, 0, -4',
            -2*spin_freq - 2*orbital_freq,
            (1./1152.)*e8,
            None,
            -2.,
            2.,
            2.
        ),
        (
            '-2*O - n',
            '2*O + n',
            '2, 2, 0, -3',
            -2*spin_freq - orbital_freq,
            (-1./4608.)*i2e6 + (11./36864.)*e8 + (1./4608.)*e6,
            None,
            -1.,
            2.,
            2.
        ),
        (
            '-2*O + n',
            '2*O - n',
            '2, 2, 0, -1',
            -2*spin_freq + orbital_freq,
            (-23./1440.)*i6e2 + (-11./768.)*i4e4 + (11./192.)*i4e2 + (-13./1536.)*i2e6 + (1./32.)*i2e4 + (-1./8.)*i2e2 + (305./36864.)*e8 + (13./1536.)*e6 + (-1./32.)*e4 + (1./8.)*e2,
            None,
            1.,
            2.,
            2.
        ),
        (
            '-2*O + 2*n',
            '2*O - 2*n',
            '2, 2, 0, 0',
            -2*spin_freq + 2*orbital_freq,
            (1957./161280.)*i8 + (23./72.)*i6e2 + (-23./360.)*i6 + (231./128.)*i4e4 + (-55./48.)*i4e2 + (11./48.)*i4 + (79./36.)*i2e6 + (-63./16.)*i2e4 + (5./2.)*i2e2 + (-1./2.)*i2 + (3697./4608.)*e8 + (-79./36.)*e6 + (63./16.)*e4 + (-5./2.)*e2 + (1./2.),
            None,
            2.,
            2.,
            2.
        ),
        (
            '-2*O + 3*n',
            '2*O - 3*n',
            '2, 2, 0, 1',
            -2*spin_freq + 3*orbital_freq,
            (-1127./1440.)*i6e2 + (-3157./256.)*i4e4 + (539./192.)*i4e2 + (-21975./512.)*i2e6 + (861./32.)*i2e4 + (-49./8.)*i2e2 + (-135771./4096.)*e8 + (21975./512.)*e6 + (-861./32.)*e4 + (49./8.)*e2,
            None,
            3.,
            2.,
            2.
        ),
        (
            '-2*O + 4*n',
            '2*O - 4*n',
            '2, 2, 0, 2',
            -2*spin_freq + 4*orbital_freq,
            (3179./192.)*i4e4 + (1955./12.)*i2e6 + (-289./8.)*i2e4 + (83551./288.)*e8 + (-1955./12.)*e6 + (289./8.)*e4,
            None,
            4.,
            2.,
            2.
        ),
        (
            '-2*O + 5*n',
            '2*O - 5*n',
            '2, 2, 0, 3',
            -2*spin_freq + 5*orbital_freq,
            (-714025./4608.)*i2e6 + (-27483625./36864.)*e8 + (714025./4608.)*e6,
            None,
            5.,
            2.,
            2.
        ),
        (
            '-2*O + 6*n',
            '2*O - 6*n',
            '2, 2, 0, 4',
            -2*spin_freq + 6*orbital_freq,
            (284089./512.)*e8,
            None,
            6.,
            2.,
            2.
        ),
        (
            '-2*O - 2*n',
            '2*O + 2*n',
            '2, 2, 1, -2',
            -2*spin_freq - 2*orbital_freq,
            (81./128.)*i4e4,
            None,
            -2.,
            0.,
            2.
        ),
        (
            '-2*O - n',
            '2*O + n',
            '2, 2, 1, -1',
            -2*spin_freq - orbital_freq,
            (-3./16.)*i6e2 + (81./128.)*i4e4 + (9./32.)*i4e2,
            None,
            -1.,
            0.,
            2.
        ),
        (
            '-2*O',
            '2*O',
            '2, 2, 1, 0',
            -2*spin_freq,
            (1./40.)*i8 + (-1./4.)*i6e2 + (-1./12.)*i6 + (3./4.)*i4e4 + (3./8.)*i4e2 + (1./8.)*i4,
            None,
            0.,
            0.,
            2.
        ),
        (
            '-2*O + n',
            '2*O - n',
            '2, 2, 1, 1',
            -2*spin_freq + orbital_freq,
            (-3./16.)*i6e2 + (81./128.)*i4e4 + (9./32.)*i4e2,
            None,
            1.,
            0.,
            2.
        ),
        (
            '-2*O + 2*n',
            '2*O - 2*n',
            '2, 2, 1, 2',
            -2*spin_freq + 2*orbital_freq,
            (81./128.)*i4e4,
            None,
            2.,
            0.,
            2.
        ),
        (
            '-2*O - 2*n',
            '2*O + 2*n',
            '2, 2, 2, 0',
            -2*spin_freq - 2*orbital_freq,
            (1./512.)*i8,
            None,
            -2.,
            -2.,
            2.
        )
    )

    return mode_data_output

@njit
def nsr_modes_t10(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray) \
        -> OutputType:
    """ Non-sync Tidal Modes for truncation level 10, tidal order 2

    Parameters
    ----------
    orbital_freq : np.ndarray
        Planet's orbital mean motion [rad s-1]
    spin_freq : np.ndarray
        Planet's spin frequency [rad s-1]
    eccentricity : np.ndarray
        Planet's eccentricity
    inclination : np.ndarray
        Planet's inclination [rads]

    Returns
    -------
    mode_data_tuple : OutputType
        List of mode datas. Each mode data contains information to calculate the heating and torques.
    """
    # Tidal Potential Builder for l-order = 2, NSR = True.
    # Max Eccentricity Order = 10
    # Max Inclination Order = 10
    # Max q = 6.
    # Number of unique modes = 32.
    # Number of unique frequencies = 27.

    i2 = inclination**2
    i4 = inclination**4
    i6 = inclination**6
    i8 = inclination**8
    i10 = inclination**10
    e2 = eccentricity**2
    e4 = eccentricity**4
    e6 = eccentricity**6
    e8 = eccentricity**8
    e10 = eccentricity**10
    i2e2 = i2*e2
    i2e4 = i2*e4
    i2e6 = i2*e6
    i2e8 = i2*e8
    i4e2 = i4*e2
    i4e4 = i4*e4
    i4e6 = i4*e6
    i6e2 = i6*e2
    i6e4 = i6*e4
    i8e2 = i8*e2

    mode_data_output = (
        (
            '-n',
            'n',
            '2, 0, 0, -3',
            -orbital_freq,
            (1./24576.)*i4e6,
            None,
            -1.,
            2.,
            0.
        ),
        (
            'n',
            'n',
            '2, 0, 0, -1',
            orbital_freq,
            (3./640.)*i8e2 + (1./256.)*i6e4 + (-1./64.)*i6e2 + (13./8192.)*i4e6 + (-3./512.)*i4e4 + (3./128.)*i4e2,
            None,
            1.,
            2.,
            0.
        ),
        (
            '2*n',
            '2*n',
            '2, 0, 0, 0',
            2*orbital_freq,
            (-17./5040.)*i10 + (-3./32.)*i8e2 + (3./160.)*i8 + (-63./128.)*i6e4 + (5./16.)*i6e2 + (-1./16.)*i6 + (
                        -79./192.)*i4e6 + (189./256.)*i4e4 + (-15./32.)*i4e2 + (3./32.)*i4,
            None,
            2.,
            2.,
            0.
        ),
        (
            '3*n',
            '3*n',
            '2, 0, 0, 1',
            3*orbital_freq,
            (147./640.)*i8e2 + (861./256.)*i6e4 + (-49./64.)*i6e2 + (65925./8192.)*i4e6 + (-2583./512.)*i4e4 + (
                        147./128.)*i4e2,
            None,
            3.,
            2.,
            0.
        ),
        (
            '4*n',
            '4*n',
            '2, 0, 0, 2',
            4*orbital_freq,
            (-289./64.)*i6e4 + (-1955./64.)*i4e6 + (867./128.)*i4e4,
            None,
            4.,
            2.,
            0.
        ),
        (
            '5*n',
            '5*n',
            '2, 0, 0, 3',
            5*orbital_freq,
            (714025./24576.)*i4e6,
            None,
            5.,
            2.,
            0.
        ),
        (
            '-5*n',
            '5*n',
            '2, 0, 1, -5',
            -5*orbital_freq,
            (1047843./131072.)*e10,
            None,
            -5.,
            0.,
            0.
        ),
        (
            '-4*n',
            '4*n',
            '2, 0, 1, -4',
            -4*orbital_freq,
            (-5929./512.)*i2e8 + (3311./2560.)*e10 + (5929./1536.)*e8,
            None,
            -4.,
            0.,
            0.
        ),
        (
            '-3*n',
            '3*n',
            '2, 0, 1, -3',
            -3*orbital_freq,
            (36517./6144.)*i4e6 + (-20829./4096.)*i2e8 + (-2809./512.)*i2e6 + (2006627./655360.)*e10 + (
                        6943./4096.)*e8 + (2809./1536.)*e6,
            None,
            -3.,
            0.,
            0.
        ),
        (
            '-2*n',
            '2*n',
            '2, 0, 1, -2',
            -2*orbital_freq,
            (-441./320.)*i6e4 + (273./64.)*i4e6 + (351./128.)*i4e4 + (-1661./256.)*i2e8 + (-63./16.)*i2e6 + (
                        -81./32.)*i2e4 + (3919./1280.)*e10 + (1661./768.)*e8 + (21./16.)*e6 + (27./32.)*e4,
            None,
            -2.,
            0.,
            0.
        ),
        (
            '-n',
            'n',
            '2, 0, 1, -1',
            -orbital_freq,
            (193./1120.)*i8e2 + (-441./320.)*i6e4 + (-49./80.)*i6e2 + (9945./2048.)*i4e6 + (351./128.)*i4e4 + (
                        39./32.)*i4e2 + (-28211./4096.)*i2e8 + (-2295./512.)*i2e6 + (-81./32.)*i2e4 + (-9./8.)*i2e2 + (
                        1064887./327680.)*e10 + (28211./12288.)*e8 + (765./512.)*e6 + (27./32.)*e4 + (3./8.)*e2,
            None,
            -1.,
            0.,
            0.
        ),
        (
            'n',
            'n',
            '2, 0, 1, 1',
            orbital_freq,
            (193./1120.)*i8e2 + (-441./320.)*i6e4 + (-49./80.)*i6e2 + (9945./2048.)*i4e6 + (351./128.)*i4e4 + (
                        39./32.)*i4e2 + (-28211./4096.)*i2e8 + (-2295./512.)*i2e6 + (-81./32.)*i2e4 + (-9./8.)*i2e2 + (
                        1064887./327680.)*e10 + (28211./12288.)*e8 + (765./512.)*e6 + (27./32.)*e4 + (3./8.)*e2,
            None,
            1.,
            0.,
            0.
        ),
        (
            '2*n',
            '2*n',
            '2, 0, 1, 2',
            2*orbital_freq,
            (-441./320.)*i6e4 + (273./64.)*i4e6 + (351./128.)*i4e4 + (-1661./256.)*i2e8 + (-63./16.)*i2e6 + (
                        -81./32.)*i2e4 + (3919./1280.)*e10 + (1661./768.)*e8 + (21./16.)*e6 + (27./32.)*e4,
            None,
            2.,
            0.,
            0.
        ),
        (
            '3*n',
            '3*n',
            '2, 0, 1, 3',
            3*orbital_freq,
            (36517./6144.)*i4e6 + (-20829./4096.)*i2e8 + (-2809./512.)*i2e6 + (2006627./655360.)*e10 + (
                        6943./4096.)*e8 + (2809./1536.)*e6,
            None,
            3.,
            0.,
            0.
        ),
        (
            '4*n',
            '4*n',
            '2, 0, 1, 4',
            4*orbital_freq,
            (-5929./512.)*i2e8 + (3311./2560.)*e10 + (5929./1536.)*e8,
            None,
            4.,
            0.,
            0.
        ),
        (
            '5*n',
            '5*n',
            '2, 0, 1, 5',
            5*orbital_freq,
            (1047843./131072.)*e10,
            None,
            5.,
            0.,
            0.
        ),
        (
            '-5*n',
            '5*n',
            '2, 0, 2, -3',
            -5*orbital_freq,
            (714025./24576.)*i4e6,
            None,
            -5.,
            -2.,
            0.
        ),
        (
            '-4*n',
            '4*n',
            '2, 0, 2, -2',
            -4*orbital_freq,
            (-289./64.)*i6e4 + (-1955./64.)*i4e6 + (867./128.)*i4e4,
            None,
            -4.,
            -2.,
            0.
        ),
        (
            '-3*n',
            '3*n',
            '2, 0, 2, -1',
            -3*orbital_freq,
            (147./640.)*i8e2 + (861./256.)*i6e4 + (-49./64.)*i6e2 + (65925./8192.)*i4e6 + (-2583./512.)*i4e4 + (
                        147./128.)*i4e2,
            None,
            -3.,
            -2.,
            0.
        ),
        (
            '-2*n',
            '2*n',
            '2, 0, 2, 0',
            -2*orbital_freq,
            (-17./5040.)*i10 + (-3./32.)*i8e2 + (3./160.)*i8 + (-63./128.)*i6e4 + (5./16.)*i6e2 + (-1./16.)*i6 + (
                        -19./48.)*i4e6 + (189./256.)*i4e4 + (-15./32.)*i4e2 + (3./32.)*i4,
            None,
            -2.,
            -2.,
            0.
        ),
        (
            '-n',
            'n',
            '2, 0, 2, 1',
            -orbital_freq,
            (3./640.)*i8e2 + (1./256.)*i6e4 + (-1./64.)*i6e2 + (13./8192.)*i4e6 + (-3./512.)*i4e4 + (3./128.)*i4e2,
            None,
            -1.,
            -2.,
            0.
        ),
        (
            'n',
            'n',
            '2, 0, 2, 3',
            orbital_freq,
            (1./24576.)*i4e6,
            None,
            1.,
            -2.,
            0.
        ),
        (
            '-O - 2*n',
            'O + 2*n',
            '2, 1, 0, -4',
            -spin_freq - 2*orbital_freq,
            (1./1152.)*i2e8,
            None,
            -2.,
            2.,
            1.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 0, -3',
            -spin_freq - orbital_freq,
            (-5./27648.)*i4e6 + (11./36864.)*i2e8 + (1./4608.)*i2e6,
            None,
            -1.,
            2.,
            1.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 0, -1',
            -spin_freq + orbital_freq,
            (-145./16128.)*i8e2 + (-227./23040.)*i6e4 + (227./5760.)*i6e2 + (-65./9216.)*i4e6 + (5./192.)*i4e4 + (
                        -5./48.)*i4e2 + (305./36864.)*i2e8 + (13./1536.)*i2e6 + (-1./32.)*i2e4 + (1./8.)*i2e2,
            None,
            1.,
            2.,
            1.
        ),
        (
            '-O + 2*n',
            'O - 2*n',
            '2, 1, 0, 0',
            -spin_freq + 2*orbital_freq,
            (40277./7257600.)*i10 + (725./4032.)*i8e2 + (-145./4032.)*i8 + (1589./1280.)*i6e4 + (-227./288.)*i6e2 + (
                        227./1440.)*i6 + (395./216.)*i4e6 + (-105./32.)*i4e4 + (25./12.)*i4e2 + (-5./12.)*i4 + (
                        3697./4608.)*i2e8 + (-79./36.)*i2e6 + (63./16.)*i2e4 + (-5./2.)*i2e2 + (1./2.)*i2,
            None,
            2.,
            2.,
            1.
        ),
        (
            '-O + 3*n',
            'O - 3*n',
            '2, 1, 0, 1',
            -spin_freq + 3*orbital_freq,
            (-1015./2304.)*i8e2 + (-65149./7680.)*i6e4 + (11123./5760.)*i6e2 + (-36625./1024.)*i4e6 + (
                        1435./64.)*i4e4 + (-245./48.)*i4e2 + (-135771./4096.)*i2e8 + (21975./512.)*i2e6 + (
                        -861./32.)*i2e4 + (49./8.)*i2e2,
            None,
            3.,
            2.,
            1.
        ),
        (
            '-O + 4*n',
            'O - 4*n',
            '2, 1, 0, 2',
            -spin_freq + 4*orbital_freq,
            (65603./5760.)*i6e4 + (9775./72.)*i4e6 + (-1445./48.)*i4e4 + (83551./288.)*i2e8 + (-1955./12.)*i2e6 + (
                        289./8.)*i2e4,
            None,
            4.,
            2.,
            1.
        ),
        (
            '-O + 5*n',
            'O - 5*n',
            '2, 1, 0, 3',
            -spin_freq + 5*orbital_freq,
            (-3570125./27648.)*i4e6 + (-27483625./36864.)*i2e8 + (714025./4608.)*i2e6,
            None,
            5.,
            2.,
            1.
        ),
        (
            '-O + 6*n',
            'O - 6*n',
            '2, 1, 0, 4',
            -spin_freq + 6*orbital_freq,
            (284089./512.)*i2e8,
            None,
            6.,
            2.,
            1.
        ),
        (
            '-O - 4*n',
            'O + 4*n',
            '2, 1, 1, -4',
            -spin_freq - 4*orbital_freq,
            (5929./512.)*i2e8,
            None,
            -4.,
            0.,
            1.
        ),
        (
            '-O - 3*n',
            'O + 3*n',
            '2, 1, 1, -3',
            -spin_freq - 3*orbital_freq,
            (-2809./384.)*i4e6 + (20829./4096.)*i2e8 + (2809./512.)*i2e6,
            None,
            -3.,
            0.,
            1.
        ),
        (
            '-O - 2*n',
            'O + 2*n',
            '2, 1, 1, -2',
            -spin_freq - 2*orbital_freq,
            (9./5.)*i6e4 + (-21./4.)*i4e6 + (-27./8.)*i4e4 + (1661./256.)*i2e8 + (63./16.)*i2e6 + (81./32.)*i2e4,
            None,
            -2.,
            0.,
            1.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 1, -1',
            -spin_freq - orbital_freq,
            (-8./35.)*i8e2 + (9./5.)*i6e4 + (4./5.)*i6e2 + (-765./128.)*i4e6 + (-27./8.)*i4e4 + (-3./2.)*i4e2 + (
                        28211./4096.)*i2e8 + (2295./512.)*i2e6 + (81./32.)*i2e4 + (9./8.)*i2e2,
            None,
            -1.,
            0.,
            1.
        ),
        (
            '-O',
            'O',
            '2, 1, 1, 0',
            -spin_freq,
            (256./14175.)*i10 + (-32./105.)*i8e2 + (-32./315.)*i8 + (32./15.)*i6e4 + (16./15.)*i6e2 + (16./45.)*i6 + (
                        -20./3.)*i4e6 + -4.*i4e4 + -2.*i4e2 + (-2./3.)*i4 + (15./2.)*i2e8 + 5.*i2e6 + 3.*i2e4 + (
                        3./2.)*i2e2 + (1./2.)*i2,
            None,
            0.,
            0.,
            1.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 1, 1',
            -spin_freq + orbital_freq,
            (-8./35.)*i8e2 + (9./5.)*i6e4 + (4./5.)*i6e2 + (-765./128.)*i4e6 + (-27./8.)*i4e4 + (-3./2.)*i4e2 + (
                        28211./4096.)*i2e8 + (2295./512.)*i2e6 + (81./32.)*i2e4 + (9./8.)*i2e2,
            None,
            1.,
            0.,
            1.
        ),
        (
            '-O + 2*n',
            'O - 2*n',
            '2, 1, 1, 2',
            -spin_freq + 2*orbital_freq,
            (9./5.)*i6e4 + (-21./4.)*i4e6 + (-27./8.)*i4e4 + (1661./256.)*i2e8 + (63./16.)*i2e6 + (81./32.)*i2e4,
            None,
            2.,
            0.,
            1.
        ),
        (
            '-O + 3*n',
            'O - 3*n',
            '2, 1, 1, 3',
            -spin_freq + 3*orbital_freq,
            (-2809./384.)*i4e6 + (20829./4096.)*i2e8 + (2809./512.)*i2e6,
            None,
            3.,
            0.,
            1.
        ),
        (
            '-O + 4*n',
            'O - 4*n',
            '2, 1, 1, 4',
            -spin_freq + 4*orbital_freq,
            (5929./512.)*i2e8,
            None,
            4.,
            0.,
            1.
        ),
        (
            '-O - 4*n',
            'O + 4*n',
            '2, 1, 2, -2',
            -spin_freq - 4*orbital_freq,
            (289./128.)*i6e4,
            None,
            -4.,
            -2.,
            1.
        ),
        (
            '-O - 3*n',
            'O + 3*n',
            '2, 1, 2, -1',
            -spin_freq - 3*orbital_freq,
            (-49./256.)*i8e2 + (-861./512.)*i6e4 + (49./128.)*i6e2,
            None,
            -3.,
            -2.,
            1.
        ),
        (
            '-O - 2*n',
            'O + 2*n',
            '2, 1, 2, 0',
            -spin_freq - 2*orbital_freq,
            (9./2560.)*i10 + (5./64.)*i8e2 + (-1./64.)*i8 + (63./256.)*i6e4 + (-5./32.)*i6e2 + (1./32.)*i6,
            None,
            -2.,
            -2.,
            1.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 2, 1',
            -spin_freq - orbital_freq,
            (-1./256.)*i8e2 + (-1./512.)*i6e4 + (1./128.)*i6e2,
            None,
            -1.,
            -2.,
            1.
        ),
        (
            '-2*O - 3*n',
            '2*O + 3*n',
            '2, 2, 0, -5',
            -2*spin_freq - 3*orbital_freq,
            (6561./3276800.)*e10,
            None,
            -3.,
            2.,
            2.
        ),
        (
            '-2*O - 2*n',
            '2*O + 2*n',
            '2, 2, 0, -4',
            -2*spin_freq - 2*orbital_freq,
            (-1./1152.)*i2e8 + (7./5760.)*e10 + (1./1152.)*e8,
            None,
            -2.,
            2.,
            2.
        ),
        (
            '-2*O - n',
            '2*O + n',
            '2, 2, 0, -3',
            -2*spin_freq - orbital_freq,
            (11./110592.)*i4e6 + (-11./36864.)*i2e8 + (-1./4608.)*i2e6 + (619./1966080.)*e10 + (11./36864.)*e8 + (
                        1./4608.)*e6,
            None,
            -1.,
            2.,
            2.
        ),
        (
            '-2*O + n',
            '2*O - n',
            '2, 2, 0, -1',
            -2*spin_freq + orbital_freq,
            (1957./645120.)*i8e2 + (23./5760.)*i6e4 + (-23./1440.)*i6e2 + (143./36864.)*i4e6 + (-11./768.)*i4e4 + (
                        11./192.)*i4e2 + (-305./36864.)*i2e8 + (-13./1536.)*i2e6 + (1./32.)*i2e4 + (-1./8.)*i2e2 + (
                        1733./327680.)*e10 + (305./36864.)*e8 + (13./1536.)*e6 + (-1./32.)*e4 + (1./8.)*e2,
            None,
            1.,
            2.,
            2.
        ),
        (
            '-2*O + 2*n',
            '2*O - 2*n',
            '2, 2, 0, 0',
            -2*spin_freq + 2*orbital_freq,
            (-12107./7257600.)*i10 + (-1957./32256.)*i8e2 + (1957./161280.)*i8 + (-161./320.)*i6e4 + (23./72.)*i6e2 + (
                        -23./360.)*i6 + (-869./864.)*i4e6 + (231./128.)*i4e4 + (-55./48.)*i4e2 + (11./48.)*i4 + (
                        -3697./4608.)*i2e8 + (79./36.)*i2e6 + (-63./16.)*i2e4 + (5./2.)*i2e2 + (-1./2.)*i2 + (
                        -11041./38400.)*e10 + (3697./4608.)*e8 + (-79./36.)*e6 + (63./16.)*e4 + (-5./2.)*e2 + (1./2.),
            None,
            2.,
            2.,
            2.
        ),
        (
            '-2*O + 3*n',
            '2*O - 3*n',
            '2, 2, 0, 1',
            -2*spin_freq + 3*orbital_freq,
            (13699./92160.)*i8e2 + (6601./1920.)*i6e4 + (-1127./1440.)*i6e2 + (80575./4096.)*i4e6 + (
                        -3157./256.)*i4e4 + (539./192.)*i4e2 + (135771./4096.)*i2e8 + (-21975./512.)*i2e6 + (
                        861./32.)*i2e4 + (-49./8.)*i2e2 + (5568309./327680.)*e10 + (-135771./4096.)*e8 + (
                        21975./512.)*e6 + (-861./32.)*e4 + (49./8.)*e2,
            None,
            3.,
            2.,
            2.
        ),
        (
            '-2*O + 4*n',
            '2*O - 4*n',
            '2, 2, 0, 2',
            -2*spin_freq + 4*orbital_freq,
            (-6647./1440.)*i6e4 + (-21505./288.)*i4e6 + (3179./192.)*i4e4 + (-83551./288.)*i2e8 + (1955./12.)*i2e6 + (
                        -289./8.)*i2e4 + (-134209./480.)*e10 + (83551./288.)*e8 + (-1955./12.)*e6 + (289./8.)*e4,
            None,
            4.,
            2.,
            2.
        ),
        (
            '-2*O + 5*n',
            '2*O - 5*n',
            '2, 2, 0, 3',
            -2*spin_freq + 5*orbital_freq,
            (7854275./110592.)*i4e6 + (27483625./36864.)*i2e8 + (-714025./4608.)*i2e6 + (587225375./393216.)*e10 + (
                        -27483625./36864.)*e8 + (714025./4608.)*e6,
            None,
            5.,
            2.,
            2.
        ),
        (
            '-2*O + 6*n',
            '2*O - 6*n',
            '2, 2, 0, 4',
            -2*spin_freq + 6*orbital_freq,
            (-284089./512.)*i2e8 + (-7369791./2560.)*e10 + (284089./512.)*e8,
            None,
            6.,
            2.,
            2.
        ),
        (
            '-2*O + 7*n',
            '2*O - 7*n',
            '2, 2, 0, 5',
            -2*spin_freq + 7*orbital_freq,
            (52142352409./29491200.)*e10,
            None,
            7.,
            2.,
            2.
        ),
        (
            '-2*O - 3*n',
            '2*O + 3*n',
            '2, 2, 1, -3',
            -2*spin_freq - 3*orbital_freq,
            (2809./2048.)*i4e6,
            None,
            -3.,
            0.,
            2.
        ),
        (
            '-2*O - 2*n',
            '2*O + 2*n',
            '2, 2, 1, -2',
            -2*spin_freq - 2*orbital_freq,
            (-27./64.)*i6e4 + (63./64.)*i4e6 + (81./128.)*i4e4,
            None,
            -2.,
            0.,
            2.
        ),
        (
            '-2*O - n',
            '2*O + n',
            '2, 2, 1, -1',
            -2*spin_freq - orbital_freq,
            (9./160.)*i8e2 + (-27./64.)*i6e4 + (-3./16.)*i6e2 + (2295./2048.)*i4e6 + (81./128.)*i4e4 + (9./32.)*i4e2,
            None,
            -1.,
            0.,
            2.
        ),
        (
            '-2*O',
            '2*O',
            '2, 2, 1, 0',
            -2*spin_freq,
            (-17./3780.)*i10 + (3./40.)*i8e2 + (1./40.)*i8 + (-1./2.)*i6e4 + (-1./4.)*i6e2 + (-1./12.)*i6 + (
                        5./4.)*i4e6 + (3./4.)*i4e4 + (3./8.)*i4e2 + (1./8.)*i4,
            None,
            0.,
            0.,
            2.
        ),
        (
            '-2*O + n',
            '2*O - n',
            '2, 2, 1, 1',
            -2*spin_freq + orbital_freq,
            (9./160.)*i8e2 + (-27./64.)*i6e4 + (-3./16.)*i6e2 + (2295./2048.)*i4e6 + (81./128.)*i4e4 + (9./32.)*i4e2,
            None,
            1.,
            0.,
            2.
        ),
        (
            '-2*O + 2*n',
            '2*O - 2*n',
            '2, 2, 1, 2',
            -2*spin_freq + 2*orbital_freq,
            (-27./64.)*i6e4 + (63./64.)*i4e6 + (81./128.)*i4e4,
            None,
            2.,
            0.,
            2.
        ),
        (
            '-2*O + 3*n',
            '2*O - 3*n',
            '2, 2, 1, 3',
            -2*spin_freq + 3*orbital_freq,
            (2809./2048.)*i4e6,
            None,
            3.,
            0.,
            2.
        ),
        (
            '-2*O - 3*n',
            '2*O + 3*n',
            '2, 2, 2, -1',
            -2*spin_freq - 3*orbital_freq,
            (49./2048.)*i8e2,
            None,
            -3.,
            -2.,
            2.
        ),
        (
            '-2*O - 2*n',
            '2*O + 2*n',
            '2, 2, 2, 0',
            -2*spin_freq - 2*orbital_freq,
            (-1./1536.)*i10 + (-5./512.)*i8e2 + (1./512.)*i8,
            None,
            -2.,
            -2.,
            2.
        ),
        (
            '-2*O - n',
            '2*O + n',
            '2, 2, 2, 1',
            -2*spin_freq - orbital_freq,
            (1./2048.)*i8e2,
            None,
            -1.,
            -2.,
            2.
        )
    )

    return mode_data_output

@njit
def nsr_modes_t12(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray) \
        -> OutputType:
    """ Non-sync Tidal Modes for truncation level 12, tidal order 2

    Parameters
    ----------
    orbital_freq : np.ndarray
        Planet's orbital mean motion [rad s-1]
    spin_freq : np.ndarray
        Planet's spin frequency [rad s-1]
    eccentricity : np.ndarray
        Planet's eccentricity
    inclination : np.ndarray
        Planet's inclination [rads]

    Returns
    -------
    mode_data_tuple : OutputType
        List of mode datas. Each mode data contains information to calculate the heating and torques.
    """
    # Tidal Potential Builder for l-order = 2, NSR = True.
    # Max Eccentricity Order = 12
    # Max Inclination Order = 12
    # Max q = 7.
    # Number of unique modes = 38.
    # Number of unique frequencies = 32.

    i2 = inclination**2
    i4 = inclination**4
    i6 = inclination**6
    i8 = inclination**8
    i10 = inclination**10
    i12 = inclination**12
    e2 = eccentricity**2
    e4 = eccentricity**4
    e6 = eccentricity**6
    e8 = eccentricity**8
    e10 = eccentricity**10
    e12 = eccentricity**12
    i2e2 = i2*e2
    i2e4 = i2*e4
    i2e6 = i2*e6
    i2e8 = i2*e8
    i2e10 = i2*e10
    i4e2 = i4*e2
    i4e4 = i4*e4
    i4e6 = i4*e6
    i4e8 = i4*e8
    i6e2 = i6*e2
    i6e4 = i6*e4
    i6e6 = i6*e6
    i8e2 = i8*e2
    i8e4 = i8*e4
    i10e2 = i10*e2

    mode_data_output = (
        (
            '-2*n',
            '2*n',
            '2, 0, 0, -4',
            -2*orbital_freq,
            (1./6144.)*i4e8,
            None,
            -2.,
            2.,
            0.
        ),
        (
            '-n',
            'n',
            '2, 0, 0, -3',
            -orbital_freq,
            (-1./36864.)*i6e6 + (11./196608.)*i4e8 + (1./24576.)*i4e6,
            None,
            -1.,
            2.,
            0.
        ),
        (
            'n',
            'n',
            '2, 0, 0, -1',
            orbital_freq,
            (-17./20160.)*i10e2 + (-3./2560.)*i8e4 + (3./640.)*i8e2 + (-13./12288.)*i6e6 + (1./256.)*i6e4 + (
                        -1./64.)*i6e2 + (305./196608.)*i4e8 + (13./8192.)*i4e6 + (-3./512.)*i4e4 + (3./128.)*i4e2,
            None,
            1.,
            2.,
            0.
        ),
        (
            '2*n',
            '2*n',
            '2, 0, 0, 0',
            2*orbital_freq,
            (31./75600.)*i12 + (17./1008.)*i10e2 + (-17./5040.)*i10 + (189./1280.)*i8e4 + (-3./32.)*i8e2 + (
                        3./160.)*i8 + (79./288.)*i6e6 + (-63./128.)*i6e4 + (5./16.)*i6e2 + (-1./16.)*i6 + (
                        3697./24576.)*i4e8 + (-79./192.)*i4e6 + (189./256.)*i4e4 + (-15./32.)*i4e2 + (3./32.)*i4,
            None,
            2.,
            2.,
            0.
        ),
        (
            '3*n',
            '3*n',
            '2, 0, 0, 1',
            3*orbital_freq,
            (-119./2880.)*i10e2 + (-2583./2560.)*i8e4 + (147./640.)*i8e2 + (-21975./4096.)*i6e6 + (861./256.)*i6e4 + (
                        -49./64.)*i6e2 + (-407313./65536.)*i4e8 + (65925./8192.)*i4e6 + (-2583./512.)*i4e4 + (
                        147./128.)*i4e2,
            None,
            3.,
            2.,
            0.
        ),
        (
            '4*n',
            '4*n',
            '2, 0, 0, 2',
            4*orbital_freq,
            (867./640.)*i8e4 + (1955./96.)*i6e6 + (-289./64.)*i6e4 + (83551./1536.)*i4e8 + (-1955./64.)*i4e6 + (
                        867./128.)*i4e4,
            None,
            4.,
            2.,
            0.
        ),
        (
            '5*n',
            '5*n',
            '2, 0, 0, 3',
            5*orbital_freq,
            (-714025./36864.)*i6e6 + (-27483625./196608.)*i4e8 + (714025./24576.)*i4e6,
            None,
            5.,
            2.,
            0.
        ),
        (
            '6*n',
            '6*n',
            '2, 0, 0, 4',
            6*orbital_freq,
            (852267./8192.)*i4e8,
            None,
            6.,
            2.,
            0.
        ),
        (
            '-6*n',
            '6*n',
            '2, 0, 1, -6',
            -6*orbital_freq,
            (10029889./614400.)*e12,
            None,
            -6.,
            0.,
            0.
        ),
        (
            '-5*n',
            '5*n',
            '2, 0, 1, -5',
            -5*orbital_freq,
            (-3143529./131072.)*i2e10 + (-982439./524288.)*e12 + (1047843./131072.)*e10,
            None,
            -5.,
            0.,
            0.
        ),
        (
            '-4*n',
            '4*n',
            '2, 0, 1, -4',
            -4*orbital_freq,
            (77077./6144.)*i4e8 + (-9933./2560.)*i2e10 + (-5929./512.)*i2e8 + (708871./153600.)*e12 + (
                        3311./2560.)*e10 + (5929./1536.)*e8,
            None,
            -4.,
            0.,
            0.
        ),
        (
            '-3*n',
            '3*n',
            '2, 0, 1, -3',
            -3*orbital_freq,
            (-137641./46080.)*i6e6 + (90259./16384.)*i4e8 + (36517./6144.)*i4e6 + (-6019881./655360.)*i2e10 + (
                        -20829./4096.)*i2e8 + (-2809./512.)*i2e6 + (30355211./7864320.)*e12 + (2006627./655360.)*e10 + (
                        6943./4096.)*e8 + (2809./1536.)*e6,
            None,
            -3.,
            0.,
            0.
        ),
        (
            '-2*n',
            '2*n',
            '2, 0, 1, -2',
            -2*orbital_freq,
            (1737./4480.)*i8e4 + (-343./160.)*i6e6 + (-441./320.)*i6e4 + (21593./3072.)*i4e8 + (273./64.)*i4e6 + (
                        351./128.)*i4e4 + (-11757./1280.)*i2e10 + (-1661./256.)*i2e8 + (-63./16.)*i2e6 + (
                        -81./32.)*i2e4 + (505601./122880.)*e12 + (3919./1280.)*e10 + (1661./768.)*e8 + (21./16.)*e6 + (
                        27./32.)*e4,
            None,
            -2.,
            0.,
            0.
        ),
        (
            '-n',
            'n',
            '2, 0, 1, -1',
            -orbital_freq,
            (-769./25200.)*i10e2 + (1737./4480.)*i8e4 + (193./1120.)*i8e2 + (-2499./1024.)*i6e6 + (-441./320.)*i6e4 + (
                        -49./80.)*i6e2 + (366743./49152.)*i4e8 + (9945./2048.)*i4e6 + (351./128.)*i4e4 + (
                        39./32.)*i4e2 + (-3194661./327680.)*i2e10 + (-28211./4096.)*i2e8 + (-2295./512.)*i2e6 + (
                        -81./32.)*i2e4 + (-9./8.)*i2e2 + (85553819./19660800.)*e12 + (1064887./327680.)*e10 + (
                        28211./12288.)*e8 + (765./512.)*e6 + (27./32.)*e4 + (3./8.)*e2,
            None,
            -1.,
            0.,
            0.
        ),
        (
            'n',
            'n',
            '2, 0, 1, 1',
            orbital_freq,
            (-769./25200.)*i10e2 + (1737./4480.)*i8e4 + (193./1120.)*i8e2 + (-2499./1024.)*i6e6 + (-441./320.)*i6e4 + (
                        -49./80.)*i6e2 + (366743./49152.)*i4e8 + (9945./2048.)*i4e6 + (351./128.)*i4e4 + (
                        39./32.)*i4e2 + (-3194661./327680.)*i2e10 + (-28211./4096.)*i2e8 + (-2295./512.)*i2e6 + (
                        -81./32.)*i2e4 + (-9./8.)*i2e2 + (85553819./19660800.)*e12 + (1064887./327680.)*e10 + (
                        28211./12288.)*e8 + (765./512.)*e6 + (27./32.)*e4 + (3./8.)*e2,
            None,
            1.,
            0.,
            0.
        ),
        (
            '2*n',
            '2*n',
            '2, 0, 1, 2',
            2*orbital_freq,
            (1737./4480.)*i8e4 + (-343./160.)*i6e6 + (-441./320.)*i6e4 + (21593./3072.)*i4e8 + (273./64.)*i4e6 + (
                        351./128.)*i4e4 + (-11757./1280.)*i2e10 + (-1661./256.)*i2e8 + (-63./16.)*i2e6 + (
                        -81./32.)*i2e4 + (505601./122880.)*e12 + (3919./1280.)*e10 + (1661./768.)*e8 + (21./16.)*e6 + (
                        27./32.)*e4,
            None,
            2.,
            0.,
            0.
        ),
        (
            '3*n',
            '3*n',
            '2, 0, 1, 3',
            3*orbital_freq,
            (-137641./46080.)*i6e6 + (90259./16384.)*i4e8 + (36517./6144.)*i4e6 + (-6019881./655360.)*i2e10 + (
                        -20829./4096.)*i2e8 + (-2809./512.)*i2e6 + (30355211./7864320.)*e12 + (2006627./655360.)*e10 + (
                        6943./4096.)*e8 + (2809./1536.)*e6,
            None,
            3.,
            0.,
            0.
        ),
        (
            '4*n',
            '4*n',
            '2, 0, 1, 4',
            4*orbital_freq,
            (77077./6144.)*i4e8 + (-9933./2560.)*i2e10 + (-5929./512.)*i2e8 + (708871./153600.)*e12 + (
                        3311./2560.)*e10 + (5929./1536.)*e8,
            None,
            4.,
            0.,
            0.
        ),
        (
            '5*n',
            '5*n',
            '2, 0, 1, 5',
            5*orbital_freq,
            (-3143529./131072.)*i2e10 + (-982439./524288.)*e12 + (1047843./131072.)*e10,
            None,
            5.,
            0.,
            0.
        ),
        (
            '6*n',
            '6*n',
            '2, 0, 1, 6',
            6*orbital_freq,
            (10029889./614400.)*e12,
            None,
            6.,
            0.,
            0.
        ),
        (
            '-6*n',
            '6*n',
            '2, 0, 2, -4',
            -6*orbital_freq,
            (852267./8192.)*i4e8,
            None,
            -6.,
            -2.,
            0.
        ),
        (
            '-5*n',
            '5*n',
            '2, 0, 2, -3',
            -5*orbital_freq,
            (-714025./36864.)*i6e6 + (-27483625./196608.)*i4e8 + (714025./24576.)*i4e6,
            None,
            -5.,
            -2.,
            0.
        ),
        (
            '-4*n',
            '4*n',
            '2, 0, 2, -2',
            -4*orbital_freq,
            (867./640.)*i8e4 + (1955./96.)*i6e6 + (-289./64.)*i6e4 + (83551./1536.)*i4e8 + (-1955./64.)*i4e6 + (
                        867./128.)*i4e4,
            None,
            -4.,
            -2.,
            0.
        ),
        (
            '-3*n',
            '3*n',
            '2, 0, 2, -1',
            -3*orbital_freq,
            (-119./2880.)*i10e2 + (-2583./2560.)*i8e4 + (147./640.)*i8e2 + (-21975./4096.)*i6e6 + (861./256.)*i6e4 + (
                        -49./64.)*i6e2 + (-407313./65536.)*i4e8 + (65925./8192.)*i4e6 + (-2583./512.)*i4e4 + (
                        147./128.)*i4e2,
            None,
            -3.,
            -2.,
            0.
        ),
        (
            '-2*n',
            '2*n',
            '2, 0, 2, 0',
            -2*orbital_freq,
            (31./75600.)*i12 + (17./1008.)*i10e2 + (-17./5040.)*i10 + (189./1280.)*i8e4 + (-3./32.)*i8e2 + (
                        3./160.)*i8 + (19./72.)*i6e6 + (-63./128.)*i6e4 + (5./16.)*i6e2 + (-1./16.)*i6 + (
                        2065./24576.)*i4e8 + (-19./48.)*i4e6 + (189./256.)*i4e4 + (-15./32.)*i4e2 + (3./32.)*i4,
            None,
            -2.,
            -2.,
            0.
        ),
        (
            '-n',
            'n',
            '2, 0, 2, 1',
            -orbital_freq,
            (-17./20160.)*i10e2 + (-3./2560.)*i8e4 + (3./640.)*i8e2 + (-13./12288.)*i6e6 + (1./256.)*i6e4 + (
                        -1./64.)*i6e2 + (305./196608.)*i4e8 + (13./8192.)*i4e6 + (-3./512.)*i4e4 + (3./128.)*i4e2,
            None,
            -1.,
            -2.,
            0.
        ),
        (
            'n',
            'n',
            '2, 0, 2, 3',
            orbital_freq,
            (-1./36864.)*i6e6 + (11./196608.)*i4e8 + (1./24576.)*i4e6,
            None,
            1.,
            -2.,
            0.
        ),
        (
            '2*n',
            '2*n',
            '2, 0, 2, 4',
            2*orbital_freq,
            (1./6144.)*i4e8,
            None,
            2.,
            -2.,
            0.
        ),
        (
            '-O - 3*n',
            'O + 3*n',
            '2, 1, 0, -5',
            -spin_freq - 3*orbital_freq,
            (6561./3276800.)*i2e10,
            None,
            -3.,
            2.,
            1.
        ),
        (
            '-O - 2*n',
            'O + 2*n',
            '2, 1, 0, -4',
            -spin_freq - 2*orbital_freq,
            (-5./6912.)*i4e8 + (7./5760.)*i2e10 + (1./1152.)*i2e8,
            None,
            -2.,
            2.,
            1.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 0, -3',
            -spin_freq - orbital_freq,
            (227./3317760.)*i6e6 + (-55./221184.)*i4e8 + (-5./27648.)*i4e6 + (619./1966080.)*i2e10 + (
                        11./36864.)*i2e8 + (1./4608.)*i2e6,
            None,
            -1.,
            2.,
            1.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 0, -1',
            -spin_freq + orbital_freq,
            (40277./29030400.)*i10e2 + (145./64512.)*i8e4 + (-145./16128.)*i8e2 + (2951./1105920.)*i6e6 + (
                        -227./23040.)*i6e4 + (227./5760.)*i6e2 + (-1525./221184.)*i4e8 + (-65./9216.)*i4e6 + (
                        5./192.)*i4e4 + (-5./48.)*i4e2 + (1733./327680.)*i2e10 + (305./36864.)*i2e8 + (
                        13./1536.)*i2e6 + (-1./32.)*i2e4 + (1./8.)*i2e2,
            None,
            1.,
            2.,
            1.
        ),
        (
            '-O + 2*n',
            'O - 2*n',
            '2, 1, 0, 0',
            -spin_freq + 2*orbital_freq,
            (-59123./95800320.)*i12 + (-40277./1451520.)*i10e2 + (40277./7257600.)*i10 + (-145./512.)*i8e4 + (
                        725./4032.)*i8e2 + (-145./4032.)*i8 + (-17933./25920.)*i6e6 + (1589./1280.)*i6e4 + (
                        -227./288.)*i6e2 + (227./1440.)*i6 + (-18485./27648.)*i4e8 + (395./216.)*i4e6 + (
                        -105./32.)*i4e4 + (25./12.)*i4e2 + (-5./12.)*i4 + (-11041./38400.)*i2e10 + (
                        3697./4608.)*i2e8 + (-79./36.)*i2e6 + (63./16.)*i2e4 + (-5./2.)*i2e2 + (1./2.)*i2,
            None,
            2.,
            2.,
            1.
        ),
        (
            '-O + 3*n',
            'O - 3*n',
            '2, 1, 0, 1',
            -spin_freq + 3*orbital_freq,
            (281939./4147200.)*i10e2 + (5945./3072.)*i8e4 + (-1015./2304.)*i8e2 + (332555./24576.)*i6e6 + (
                        -65149./7680.)*i6e4 + (11123./5760.)*i6e2 + (226285./8192.)*i4e8 + (-36625./1024.)*i4e6 + (
                        1435./64.)*i4e4 + (-245./48.)*i4e2 + (5568309./327680.)*i2e10 + (-135771./4096.)*i2e8 + (
                        21975./512.)*i2e6 + (-861./32.)*i2e4 + (49./8.)*i2e2,
            None,
            3.,
            2.,
            1.
        ),
        (
            '-O + 4*n',
            'O - 4*n',
            '2, 1, 0, 2',
            -spin_freq + 4*orbital_freq,
            (-41905./16128.)*i8e4 + (-88757./1728.)*i6e6 + (65603./5760.)*i6e4 + (-417755./1728.)*i4e8 + (
                        9775./72.)*i4e6 + (-1445./48.)*i4e4 + (-134209./480.)*i2e10 + (83551./288.)*i2e8 + (
                        -1955./12.)*i2e6 + (289./8.)*i2e4,
            None,
            4.,
            2.,
            1.
        ),
        (
            '-O + 5*n',
            'O - 5*n',
            '2, 1, 0, 3',
            -spin_freq + 5*orbital_freq,
            (32416735./663552.)*i6e6 + (137418125./221184.)*i4e8 + (-3570125./27648.)*i4e6 + (
                        587225375./393216.)*i2e10 + (-27483625./36864.)*i2e8 + (714025./4608.)*i2e6,
            None,
            5.,
            2.,
            1.
        ),
        (
            '-O + 6*n',
            'O - 6*n',
            '2, 1, 0, 4',
            -spin_freq + 6*orbital_freq,
            (-1420445./3072.)*i4e8 + (-7369791./2560.)*i2e10 + (284089./512.)*i2e8,
            None,
            6.,
            2.,
            1.
        ),
        (
            '-O + 7*n',
            'O - 7*n',
            '2, 1, 0, 5',
            -spin_freq + 7*orbital_freq,
            (52142352409./29491200.)*i2e10,
            None,
            7.,
            2.,
            1.
        ),
        (
            '-O - 5*n',
            'O + 5*n',
            '2, 1, 1, -5',
            -spin_freq - 5*orbital_freq,
            (3143529./131072.)*i2e10,
            None,
            -5.,
            0.,
            1.
        ),
        (
            '-O - 4*n',
            'O + 4*n',
            '2, 1, 1, -4',
            -spin_freq - 4*orbital_freq,
            (-5929./384.)*i4e8 + (9933./2560.)*i2e10 + (5929./512.)*i2e8,
            None,
            -4.,
            0.,
            1.
        ),
        (
            '-O - 3*n',
            'O + 3*n',
            '2, 1, 1, -3',
            -spin_freq - 3*orbital_freq,
            (2809./720.)*i6e6 + (-6943./1024.)*i4e8 + (-2809./384.)*i4e6 + (6019881./655360.)*i2e10 + (
                        20829./4096.)*i2e8 + (2809./512.)*i2e6,
            None,
            -3.,
            0.,
            1.
        ),
        (
            '-O - 2*n',
            'O + 2*n',
            '2, 1, 1, -2',
            -spin_freq - 2*orbital_freq,
            (-18./35.)*i8e4 + (14./5.)*i6e6 + (9./5.)*i6e4 + (-1661./192.)*i4e8 + (-21./4.)*i4e6 + (-27./8.)*i4e4 + (
                        11757./1280.)*i2e10 + (1661./256.)*i2e8 + (63./16.)*i2e6 + (81./32.)*i2e4,
            None,
            -2.,
            0.,
            1.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 1, -1',
            -spin_freq - orbital_freq,
            (64./1575.)*i10e2 + (-18./35.)*i8e4 + (-8./35.)*i8e2 + (51./16.)*i6e6 + (9./5.)*i6e4 + (4./5.)*i6e2 + (
                        -28211./3072.)*i4e8 + (-765./128.)*i4e6 + (-27./8.)*i4e4 + (-3./2.)*i4e2 + (
                        3194661./327680.)*i2e10 + (28211./4096.)*i2e8 + (2295./512.)*i2e6 + (81./32.)*i2e4 + (
                        9./8.)*i2e2,
            None,
            -1.,
            0.,
            1.
        ),
        (
            '-O',
            'O',
            '2, 1, 1, 0',
            -spin_freq,
            (-1024./467775.)*i12 + (256./4725.)*i10e2 + (256./14175.)*i10 + (-64./105.)*i8e4 + (-32./105.)*i8e2 + (
                        -32./315.)*i8 + (32./9.)*i6e6 + (32./15.)*i6e4 + (16./15.)*i6e2 + (16./45.)*i6 + -10.*i4e8 + (
                        -20./3.)*i4e6 + -4.*i4e4 + -2.*i4e2 + (-2./3.)*i4 + (21./2.)*i2e10 + (
                        15./2.)*i2e8 + 5.*i2e6 + 3.*i2e4 + (3./2.)*i2e2 + (1./2.)*i2,
            None,
            0.,
            0.,
            1.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 1, 1',
            -spin_freq + orbital_freq,
            (64./1575.)*i10e2 + (-18./35.)*i8e4 + (-8./35.)*i8e2 + (51./16.)*i6e6 + (9./5.)*i6e4 + (4./5.)*i6e2 + (
                        -28211./3072.)*i4e8 + (-765./128.)*i4e6 + (-27./8.)*i4e4 + (-3./2.)*i4e2 + (
                        3194661./327680.)*i2e10 + (28211./4096.)*i2e8 + (2295./512.)*i2e6 + (81./32.)*i2e4 + (
                        9./8.)*i2e2,
            None,
            1.,
            0.,
            1.
        ),
        (
            '-O + 2*n',
            'O - 2*n',
            '2, 1, 1, 2',
            -spin_freq + 2*orbital_freq,
            (-18./35.)*i8e4 + (14./5.)*i6e6 + (9./5.)*i6e4 + (-1661./192.)*i4e8 + (-21./4.)*i4e6 + (-27./8.)*i4e4 + (
                        11757./1280.)*i2e10 + (1661./256.)*i2e8 + (63./16.)*i2e6 + (81./32.)*i2e4,
            None,
            2.,
            0.,
            1.
        ),
        (
            '-O + 3*n',
            'O - 3*n',
            '2, 1, 1, 3',
            -spin_freq + 3*orbital_freq,
            (2809./720.)*i6e6 + (-6943./1024.)*i4e8 + (-2809./384.)*i4e6 + (6019881./655360.)*i2e10 + (
                        20829./4096.)*i2e8 + (2809./512.)*i2e6,
            None,
            3.,
            0.,
            1.
        ),
        (
            '-O + 4*n',
            'O - 4*n',
            '2, 1, 1, 4',
            -spin_freq + 4*orbital_freq,
            (-5929./384.)*i4e8 + (9933./2560.)*i2e10 + (5929./512.)*i2e8,
            None,
            4.,
            0.,
            1.
        ),
        (
            '-O + 5*n',
            'O - 5*n',
            '2, 1, 1, 5',
            -spin_freq + 5*orbital_freq,
            (3143529./131072.)*i2e10,
            None,
            5.,
            0.,
            1.
        ),
        (
            '-O - 5*n',
            'O + 5*n',
            '2, 1, 2, -3',
            -spin_freq - 5*orbital_freq,
            (714025./73728.)*i6e6,
            None,
            -5.,
            -2.,
            1.
        ),
        (
            '-O - 4*n',
            'O + 4*n',
            '2, 1, 2, -2',
            -spin_freq - 4*orbital_freq,
            (-289./256.)*i8e4 + (-1955./192.)*i6e6 + (289./128.)*i6e4,
            None,
            -4.,
            -2.,
            1.
        ),
        (
            '-O - 3*n',
            'O + 3*n',
            '2, 1, 2, -1',
            -spin_freq - 3*orbital_freq,
            (441./10240.)*i10e2 + (861./1024.)*i8e4 + (-49./256.)*i8e2 + (21975./8192.)*i6e6 + (-861./512.)*i6e4 + (
                        49./128.)*i6e2,
            None,
            -3.,
            -2.,
            1.
        ),
        (
            '-O - 2*n',
            'O + 2*n',
            '2, 1, 2, 0',
            -spin_freq - 2*orbital_freq,
            (-463./967680.)*i12 + (-9./512.)*i10e2 + (9./2560.)*i10 + (-63./512.)*i8e4 + (5./64.)*i8e2 + (
                        -1./64.)*i8 + (-19./144.)*i6e6 + (63./256.)*i6e4 + (-5./32.)*i6e2 + (1./32.)*i6,
            None,
            -2.,
            -2.,
            1.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 2, 1',
            -spin_freq - orbital_freq,
            (9./10240.)*i10e2 + (1./1024.)*i8e4 + (-1./256.)*i8e2 + (13./24576.)*i6e6 + (-1./512.)*i6e4 + (
                        1./128.)*i6e2,
            None,
            -1.,
            -2.,
            1.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 2, 3',
            -spin_freq + orbital_freq,
            (1./73728.)*i6e6,
            None,
            1.,
            -2.,
            1.
        ),
        (
            '-2*O - 4*n',
            '2*O + 4*n',
            '2, 2, 0, -6',
            -2*spin_freq - 4*orbital_freq,
            (8./2025.)*e12,
            None,
            -4.,
            2.,
            2.
        ),
        (
            '-2*O - 3*n',
            '2*O + 3*n',
            '2, 2, 0, -5',
            -2*spin_freq - 3*orbital_freq,
            (-6561./3276800.)*i2e10 + (6561./2621440.)*e12 + (6561./3276800.)*e10,
            None,
            -3.,
            2.,
            2.
        ),
        (
            '-2*O - 2*n',
            '2*O + 2*n',
            '2, 2, 0, -4',
            -2*spin_freq - 2*orbital_freq,
            (11./27648.)*i4e8 + (-7./5760.)*i2e10 + (-1./1152.)*i2e8 + (949./691200.)*e12 + (7./5760.)*e10 + (
                        1./1152.)*e8,
            None,
            -2.,
            2.,
            2.
        ),
        (
            '-2*O - n',
            '2*O + n',
            '2, 2, 0, -3',
            -2*spin_freq - orbital_freq,
            (-23./829440.)*i6e6 + (121./884736.)*i4e8 + (11./110592.)*i4e6 + (-619./1966080.)*i2e10 + (
                        -11./36864.)*i2e8 + (-1./4608.)*i2e6 + (62617./212336640.)*e12 + (619./1966080.)*e10 + (
                        11./36864.)*e8 + (1./4608.)*e6,
            None,
            -1.,
            2.,
            2.
        ),
        (
            '-2*O + n',
            '2*O - n',
            '2, 2, 0, -1',
            -2*spin_freq + orbital_freq,
            (-12107./29030400.)*i10e2 + (-1957./2580480.)*i8e4 + (1957./645120.)*i8e2 + (-299./276480.)*i6e6 + (
                        23./5760.)*i6e4 + (-23./1440.)*i6e2 + (3355./884736.)*i4e8 + (143./36864.)*i4e6 + (
                        -11./768.)*i4e4 + (11./192.)*i4e2 + (-1733./327680.)*i2e10 + (-305./36864.)*i2e8 + (
                        -13./1536.)*i2e6 + (1./32.)*i2e4 + (-1./8.)*i2e2 + (277229./58982400.)*e12 + (
                        1733./327680.)*e10 + (305./36864.)*e8 + (13./1536.)*e6 + (-1./32.)*e4 + (1./8.)*e2,
            None,
            1.,
            2.,
            2.
        ),
        (
            '-2*O + 2*n',
            '2*O - 2*n',
            '2, 2, 0, 0',
            -2*spin_freq + 2*orbital_freq,
            (330367./1916006400.)*i12 + (12107./1451520.)*i10e2 + (-12107./7257600.)*i10 + (1957./20480.)*i8e4 + (
                        -1957./32256.)*i8e2 + (1957./161280.)*i8 + (1817./6480.)*i6e6 + (-161./320.)*i6e4 + (
                        23./72.)*i6e2 + (-23./360.)*i6 + (40667./110592.)*i4e8 + (-869./864.)*i4e6 + (
                        231./128.)*i4e4 + (-55./48.)*i4e2 + (11./48.)*i4 + (11041./38400.)*i2e10 + (
                        -3697./4608.)*i2e8 + (79./36.)*i2e6 + (-63./16.)*i2e4 + (5./2.)*i2e2 + (-1./2.)*i2 + (
                        34471./552960.)*e12 + (-11041./38400.)*e10 + (3697./4608.)*e8 + (-79./36.)*e6 + (63./16.)*e4 + (
                        -5./2.)*e2 + (1./2.),
            None,
            2.,
            2.,
            2.
        ),
        (
            '-2*O + 3*n',
            '2*O - 3*n',
            '2, 2, 0, 1',
            -2*spin_freq + 3*orbital_freq,
            (-84749./4147200.)*i10e2 + (-80237./122880.)*i8e4 + (13699./92160.)*i8e2 + (-33695./6144.)*i6e6 + (
                        6601./1920.)*i6e4 + (-1127./1440.)*i6e2 + (-497827./32768.)*i4e8 + (80575./4096.)*i4e6 + (
                        -3157./256.)*i4e4 + (539./192.)*i4e2 + (-5568309./327680.)*i2e10 + (135771./4096.)*i2e8 + (
                        -21975./512.)*i2e6 + (861./32.)*i2e4 + (-49./8.)*i2e2 + (-47982879./6553600.)*e12 + (
                        5568309./327680.)*e10 + (-135771./4096.)*e8 + (21975./512.)*e6 + (-861./32.)*e4 + (49./8.)*e2,
            None,
            3.,
            2.,
            2.
        ),
        (
            '-2*O + 4*n',
            '2*O - 4*n',
            '2, 2, 0, 2',
            -2*spin_freq + 4*orbital_freq,
            (565573./645120.)*i8e4 + (8993./432.)*i6e6 + (-6647./1440.)*i6e4 + (919061./6912.)*i4e8 + (
                        -21505./288.)*i4e6 + (3179./192.)*i4e4 + (134209./480.)*i2e10 + (-83551./288.)*i2e8 + (
                        1955./12.)*i2e6 + (-289./8.)*i2e4 + (8421731./46080.)*e12 + (-134209./480.)*e10 + (
                        83551./288.)*e8 + (-1955./12.)*e6 + (289./8.)*e4,
            None,
            4.,
            2.,
            2.
        ),
        (
            '-2*O + 5*n',
            '2*O - 5*n',
            '2, 2, 0, 3',
            -2*spin_freq + 5*orbital_freq,
            (-3284515./165888.)*i6e6 + (-302319875./884736.)*i4e8 + (7854275./110592.)*i4e6 + (
                        -587225375./393216.)*i2e10 + (27483625./36864.)*i2e8 + (-714025./4608.)*i2e6 + (
                        -72670996375./42467328.)*e12 + (587225375./393216.)*e10 + (-27483625./36864.)*e8 + (
                        714025./4608.)*e6,
            None,
            5.,
            2.,
            2.
        ),
        (
            '-2*O + 6*n',
            '2*O - 6*n',
            '2, 2, 0, 4',
            -2*spin_freq + 6*orbital_freq,
            (3124979./12288.)*i4e8 + (7369791./2560.)*i2e10 + (-284089./512.)*i2e8 + (659870313./102400.)*e12 + (
                        -7369791./2560.)*e10 + (284089./512.)*e8,
            None,
            6.,
            2.,
            2.
        ),
        (
            '-2*O + 7*n',
            '2*O - 7*n',
            '2, 2, 0, 5',
            -2*spin_freq + 7*orbital_freq,
            (-52142352409./29491200.)*i2e10 + (-140254152605./14155776.)*e12 + (52142352409./29491200.)*e10,
            None,
            7.,
            2.,
            2.
        ),
        (
            '-2*O + 8*n',
            '2*O - 8*n',
            '2, 2, 0, 6',
            -2*spin_freq + 8*orbital_freq,
            (5383010161./1036800.)*e12,
            None,
            8.,
            2.,
            2.
        ),
        (
            '-2*O - 4*n',
            '2*O + 4*n',
            '2, 2, 1, -4',
            -2*spin_freq - 4*orbital_freq,
            (5929./2048.)*i4e8,
            None,
            -4.,
            0.,
            2.
        ),
        (
            '-2*O - 3*n',
            '2*O + 3*n',
            '2, 2, 1, -3',
            -2*spin_freq - 3*orbital_freq,
            (-2809./3072.)*i6e6 + (20829./16384.)*i4e8 + (2809./2048.)*i4e6,
            None,
            -3.,
            0.,
            2.
        ),
        (
            '-2*O - 2*n',
            '2*O + 2*n',
            '2, 2, 1, -2',
            -2*spin_freq - 2*orbital_freq,
            (81./640.)*i8e4 + (-21./32.)*i6e6 + (-27./64.)*i6e4 + (1661./1024.)*i4e8 + (63./64.)*i4e6 + (81./128.)*i4e4,
            None,
            -2.,
            0.,
            2.
        ),
        (
            '-2*O - n',
            '2*O + n',
            '2, 2, 1, -1',
            -2*spin_freq - orbital_freq,
            (-17./1680.)*i10e2 + (81./640.)*i8e4 + (9./160.)*i8e2 + (-765./1024.)*i6e6 + (-27./64.)*i6e4 + (
                        -3./16.)*i6e2 + (28211./16384.)*i4e8 + (2295./2048.)*i4e6 + (81./128.)*i4e4 + (9./32.)*i4e2,
            None,
            -1.,
            0.,
            2.
        ),
        (
            '-2*O',
            '2*O',
            '2, 2, 1, 0',
            -2*spin_freq,
            (31./56700.)*i12 + (-17./1260.)*i10e2 + (-17./3780.)*i10 + (3./20.)*i8e4 + (3./40.)*i8e2 + (1./40.)*i8 + (
                        -5./6.)*i6e6 + (-1./2.)*i6e4 + (-1./4.)*i6e2 + (-1./12.)*i6 + (15./8.)*i4e8 + (5./4.)*i4e6 + (
                        3./4.)*i4e4 + (3./8.)*i4e2 + (1./8.)*i4,
            None,
            0.,
            0.,
            2.
        ),
        (
            '-2*O + n',
            '2*O - n',
            '2, 2, 1, 1',
            -2*spin_freq + orbital_freq,
            (-17./1680.)*i10e2 + (81./640.)*i8e4 + (9./160.)*i8e2 + (-765./1024.)*i6e6 + (-27./64.)*i6e4 + (
                        -3./16.)*i6e2 + (28211./16384.)*i4e8 + (2295./2048.)*i4e6 + (81./128.)*i4e4 + (9./32.)*i4e2,
            None,
            1.,
            0.,
            2.
        ),
        (
            '-2*O + 2*n',
            '2*O - 2*n',
            '2, 2, 1, 2',
            -2*spin_freq + 2*orbital_freq,
            (81./640.)*i8e4 + (-21./32.)*i6e6 + (-27./64.)*i6e4 + (1661./1024.)*i4e8 + (63./64.)*i4e6 + (81./128.)*i4e4,
            None,
            2.,
            0.,
            2.
        ),
        (
            '-2*O + 3*n',
            '2*O - 3*n',
            '2, 2, 1, 3',
            -2*spin_freq + 3*orbital_freq,
            (-2809./3072.)*i6e6 + (20829./16384.)*i4e8 + (2809./2048.)*i4e6,
            None,
            3.,
            0.,
            2.
        ),
        (
            '-2*O + 4*n',
            '2*O - 4*n',
            '2, 2, 1, 4',
            -2*spin_freq + 4*orbital_freq,
            (5929./2048.)*i4e8,
            None,
            4.,
            0.,
            2.
        ),
        (
            '-2*O - 4*n',
            '2*O + 4*n',
            '2, 2, 2, -2',
            -2*spin_freq - 4*orbital_freq,
            (289./2048.)*i8e4,
            None,
            -4.,
            -2.,
            2.
        ),
        (
            '-2*O - 3*n',
            '2*O + 3*n',
            '2, 2, 2, -1',
            -2*spin_freq - 3*orbital_freq,
            (-49./6144.)*i10e2 + (-861./8192.)*i8e4 + (49./2048.)*i8e2,
            None,
            -3.,
            -2.,
            2.
        ),
        (
            '-2*O - 2*n',
            '2*O + 2*n',
            '2, 2, 2, 0',
            -2*spin_freq - 2*orbital_freq,
            (19./184320.)*i12 + (5./1536.)*i10e2 + (-1./1536.)*i10 + (63./4096.)*i8e4 + (-5./512.)*i8e2 + (1./512.)*i8,
            None,
            -2.,
            -2.,
            2.
        ),
        (
            '-2*O - n',
            '2*O + n',
            '2, 2, 2, 1',
            -2*spin_freq - orbital_freq,
            (-1./6144.)*i10e2 + (-1./8192.)*i8e4 + (1./2048.)*i8e2,
            None,
            -1.,
            -2.,
            2.
        )
    )

    return mode_data_output

@njit
def nsr_modes_t14(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray) \
        -> OutputType:
    """ Non-sync Tidal Modes for truncation level 14, tidal order 2

    Parameters
    ----------
    orbital_freq : np.ndarray
        Planet's orbital mean motion [rad s-1]
    spin_freq : np.ndarray
        Planet's spin frequency [rad s-1]
    eccentricity : np.ndarray
        Planet's eccentricity
    inclination : np.ndarray
        Planet's inclination [rads]

    Returns
    -------
    mode_data_tuple : OutputType
        List of mode datas. Each mode data contains information to calculate the heating and torques.
    """
    # Tidal Potential Builder for l-order = 2, NSR = True.
    # Max Eccentricity Order = 14
    # Max Inclination Order = 14
    # Max q = 8.
    # Number of unique modes = 44.
    # Number of unique frequencies = 37.

    i2 = inclination**2
    i4 = inclination**4
    i6 = inclination**6
    i8 = inclination**8
    i10 = inclination**10
    i12 = inclination**12
    i14 = inclination**14
    e2 = eccentricity**2
    e4 = eccentricity**4
    e6 = eccentricity**6
    e8 = eccentricity**8
    e10 = eccentricity**10
    e12 = eccentricity**12
    e14 = eccentricity**14
    i2e2 = i2*e2
    i2e4 = i2*e4
    i2e6 = i2*e6
    i2e8 = i2*e8
    i2e10 = i2*e10
    i2e12 = i2*e12
    i4e2 = i4*e2
    i4e4 = i4*e4
    i4e6 = i4*e6
    i4e8 = i4*e8
    i4e10 = i4*e10
    i6e2 = i6*e2
    i6e4 = i6*e4
    i6e6 = i6*e6
    i6e8 = i6*e8
    i8e2 = i8*e2
    i8e4 = i8*e4
    i8e6 = i8*e6
    i10e2 = i10*e2
    i10e4 = i10*e4
    i12e2 = i12*e2

    mode_data_output = (
        (
            '-3*n',
            '3*n',
            '2, 0, 0, -5',
            -3*orbital_freq,
            (19683./52428800.)*i4e10,
            None,
            -3.,
            2.,
            0.
        ),
        (
            '-2*n',
            '2*n',
            '2, 0, 0, -4',
            -2*orbital_freq,
            (-1./9216.)*i6e8 + (7./30720.)*i4e10 + (1./6144.)*i4e8,
            None,
            -2.,
            2.,
            0.
        ),
        (
            '-n',
            'n',
            '2, 0, 0, -3',
            -orbital_freq,
            (1./122880.)*i8e6 + (-11./294912.)*i6e8 + (-1./36864.)*i6e6 + (619./10485760.)*i4e10 + (
                        11./196608.)*i4e8 + (1./24576.)*i4e6,
            None,
            -1.,
            2.,
            0.
        ),
        (
            'n',
            'n',
            '2, 0, 0, -1',
            orbital_freq,
            (31./302400.)*i12e2 + (17./80640.)*i10e4 + (-17./20160.)*i10e2 + (13./40960.)*i8e6 + (-3./2560.)*i8e4 + (
                        3./640.)*i8e2 + (-305./294912.)*i6e8 + (-13./12288.)*i6e6 + (1./256.)*i6e4 + (-1./64.)*i6e2 + (
                        5199./5242880.)*i4e10 + (305./196608.)*i4e8 + (13./8192.)*i4e6 + (-3./512.)*i4e4 + (
                        3./128.)*i4e2,
            None,
            1.,
            2.,
            0.
        ),
        (
            '2*n',
            '2*n',
            '2, 0, 0, 0',
            2*orbital_freq,
            (-1./27720.)*i14 + (-31./15120.)*i12e2 + (31./75600.)*i12 + (-17./640.)*i10e4 + (17./1008.)*i10e2 + (
                        -17./5040.)*i10 + (-79./960.)*i8e6 + (189./1280.)*i8e4 + (-3./32.)*i8e2 + (3./160.)*i8 + (
                        -3697./36864.)*i6e8 + (79./288.)*i6e6 + (-63./128.)*i6e4 + (5./16.)*i6e2 + (-1./16.)*i6 + (
                        -11041./204800.)*i4e10 + (3697./24576.)*i4e8 + (-79./192.)*i4e6 + (189./256.)*i4e4 + (
                        -15./32.)*i4e2 + (3./32.)*i4,
            None,
            2.,
            2.,
            0.
        ),
        (
            '3*n',
            '3*n',
            '2, 0, 0, 1',
            3*orbital_freq,
            (217./43200.)*i12e2 + (697./3840.)*i10e4 + (-119./2880.)*i10e2 + (13185./8192.)*i8e6 + (
                        -2583./2560.)*i8e4 + (147./640.)*i8e2 + (135771./32768.)*i6e8 + (-21975./4096.)*i6e6 + (
                        861./256.)*i6e4 + (-49./64.)*i6e2 + (16704927./5242880.)*i4e10 + (-407313./65536.)*i4e8 + (
                        65925./8192.)*i4e6 + (-2583./512.)*i4e4 + (147./128.)*i4e2,
            None,
            3.,
            2.,
            0.
        ),
        (
            '4*n',
            '4*n',
            '2, 0, 0, 2',
            4*orbital_freq,
            (-4913./20160.)*i10e4 + (-391./64.)*i8e6 + (867./640.)*i8e4 + (-83551./2304.)*i6e8 + (1955./96.)*i6e6 + (
                        -289./64.)*i6e4 + (-134209./2560.)*i4e10 + (83551./1536.)*i4e8 + (-1955./64.)*i4e6 + (
                        867./128.)*i4e4,
            None,
            4.,
            2.,
            0.
        ),
        (
            '5*n',
            '5*n',
            '2, 0, 0, 3',
            5*orbital_freq,
            (142805./24576.)*i8e6 + (27483625./294912.)*i6e8 + (-714025./36864.)*i6e6 + (587225375./2097152.)*i4e10 + (
                        -27483625./196608.)*i4e8 + (714025./24576.)*i4e6,
            None,
            5.,
            2.,
            0.
        ),
        (
            '6*n',
            '6*n',
            '2, 0, 0, 4',
            6*orbital_freq,
            (-284089./4096.)*i6e8 + (-22109373./40960.)*i4e10 + (852267./8192.)*i4e8,
            None,
            6.,
            2.,
            0.
        ),
        (
            '7*n',
            '7*n',
            '2, 0, 0, 5',
            7*orbital_freq,
            (52142352409./157286400.)*i4e10,
            None,
            7.,
            2.,
            0.
        ),
        (
            '-7*n',
            '7*n',
            '2, 0, 1, -7',
            -7*orbital_freq,
            (186702632281./5662310400.)*e14,
            None,
            -7.,
            0.,
            0.
        ),
        (
            '-6*n',
            '6*n',
            '2, 0, 1, -6',
            -6*orbital_freq,
            (-10029889./204800.)*i2e12 + (-2308743./179200.)*e14 + (10029889./614400.)*e12,
            None,
            -6.,
            0.,
            0.
        ),
        (
            '-5*n',
            '5*n',
            '2, 0, 1, -5',
            -5*orbital_freq,
            (13621959./524288.)*i4e10 + (2947317./524288.)*i2e12 + (-3143529./131072.)*i2e10 + (
                        13524196093./1585446912.)*e14 + (-982439./524288.)*e12 + (1047843./131072.)*e10,
            None,
            -5.,
            0.,
            0.
        ),
        (
            '-4*n',
            '4*n',
            '2, 0, 1, -4',
            -4*orbital_freq,
            (-290521./46080.)*i6e8 + (43043./10240.)*i4e10 + (77077./6144.)*i4e8 + (-708871./51200.)*i2e12 + (
                        -9933./2560.)*i2e10 + (-5929./512.)*i2e8 + (3027367./691200.)*e14 + (708871./153600.)*e12 + (
                        3311./2560.)*e10 + (5929./1536.)*e8,
            None,
            -4.,
            0.,
            0.
        ),
        (
            '-3*n',
            '3*n',
            '2, 0, 1, -3',
            -3*orbital_freq,
            (542137./645120.)*i8e6 + (-340207./122880.)*i6e8 + (-137641./46080.)*i6e6 + (26086151./2621440.)*i4e10 + (
                        90259./16384.)*i4e8 + (36517./6144.)*i4e6 + (-30355211./2621440.)*i2e12 + (
                        -6019881./655360.)*i2e10 + (-20829./4096.)*i2e8 + (-2809./512.)*i2e6 + (
                        1055142483./209715200.)*e14 + (30355211./7864320.)*e12 + (2006627./655360.)*e10 + (
                        6943./4096.)*e8 + (2809./1536.)*e6,
            None,
            -3.,
            0.,
            0.
        ),
        (
            '-2*n',
            '2*n',
            '2, 0, 1, -2',
            -2*orbital_freq,
            (-769./11200.)*i10e4 + (193./320.)*i8e6 + (1737./4480.)*i8e4 + (-81389./23040.)*i6e8 + (-343./160.)*i6e6 + (
                        -441./320.)*i6e4 + (50947./5120.)*i4e10 + (21593./3072.)*i4e8 + (273./64.)*i4e6 + (
                        351./128.)*i4e4 + (-505601./40960.)*i2e12 + (-11757./1280.)*i2e10 + (-1661./256.)*i2e8 + (
                        -63./16.)*i2e6 + (-81./32.)*i2e4 + (25565893./4838400.)*e14 + (505601./122880.)*e12 + (
                        3919./1280.)*e10 + (1661./768.)*e8 + (21./16.)*e6 + (27./32.)*e4,
            None,
            -2.,
            0.,
            0.
        ),
        (
            '-n',
            'n',
            '2, 0, 1, -1',
            -orbital_freq,
            (439./118800.)*i12e2 + (-769./11200.)*i10e4 + (-769./25200.)*i10e2 + (9843./14336.)*i8e6 + (
                        1737./4480.)*i8e4 + (193./1120.)*i8e2 + (-1382339./368640.)*i6e8 + (-2499./1024.)*i6e6 + (
                        -441./320.)*i6e4 + (-49./80.)*i6e2 + (13843531./1310720.)*i4e10 + (366743./49152.)*i4e8 + (
                        9945./2048.)*i4e6 + (351./128.)*i4e4 + (39./32.)*i4e2 + (-85553819./6553600.)*i2e12 + (
                        -3194661./327680.)*i2e10 + (-28211./4096.)*i2e8 + (-2295./512.)*i2e6 + (-81./32.)*i2e4 + (
                        -9./8.)*i2e2 + (24654653741./4404019200.)*e14 + (85553819./19660800.)*e12 + (
                        1064887./327680.)*e10 + (28211./12288.)*e8 + (765./512.)*e6 + (27./32.)*e4 + (3./8.)*e2,
            None,
            -1.,
            0.,
            0.
        ),
        (
            'n',
            'n',
            '2, 0, 1, 1',
            orbital_freq,
            (439./118800.)*i12e2 + (-769./11200.)*i10e4 + (-769./25200.)*i10e2 + (9843./14336.)*i8e6 + (
                        1737./4480.)*i8e4 + (193./1120.)*i8e2 + (-1382339./368640.)*i6e8 + (-2499./1024.)*i6e6 + (
                        -441./320.)*i6e4 + (-49./80.)*i6e2 + (13843531./1310720.)*i4e10 + (366743./49152.)*i4e8 + (
                        9945./2048.)*i4e6 + (351./128.)*i4e4 + (39./32.)*i4e2 + (-85553819./6553600.)*i2e12 + (
                        -3194661./327680.)*i2e10 + (-28211./4096.)*i2e8 + (-2295./512.)*i2e6 + (-81./32.)*i2e4 + (
                        -9./8.)*i2e2 + (24654653741./4404019200.)*e14 + (85553819./19660800.)*e12 + (
                        1064887./327680.)*e10 + (28211./12288.)*e8 + (765./512.)*e6 + (27./32.)*e4 + (3./8.)*e2,
            None,
            1.,
            0.,
            0.
        ),
        (
            '2*n',
            '2*n',
            '2, 0, 1, 2',
            2*orbital_freq,
            (-769./11200.)*i10e4 + (193./320.)*i8e6 + (1737./4480.)*i8e4 + (-81389./23040.)*i6e8 + (-343./160.)*i6e6 + (
                        -441./320.)*i6e4 + (50947./5120.)*i4e10 + (21593./3072.)*i4e8 + (273./64.)*i4e6 + (
                        351./128.)*i4e4 + (-505601./40960.)*i2e12 + (-11757./1280.)*i2e10 + (-1661./256.)*i2e8 + (
                        -63./16.)*i2e6 + (-81./32.)*i2e4 + (25565893./4838400.)*e14 + (505601./122880.)*e12 + (
                        3919./1280.)*e10 + (1661./768.)*e8 + (21./16.)*e6 + (27./32.)*e4,
            None,
            2.,
            0.,
            0.
        ),
        (
            '3*n',
            '3*n',
            '2, 0, 1, 3',
            3*orbital_freq,
            (542137./645120.)*i8e6 + (-340207./122880.)*i6e8 + (-137641./46080.)*i6e6 + (26086151./2621440.)*i4e10 + (
                        90259./16384.)*i4e8 + (36517./6144.)*i4e6 + (-30355211./2621440.)*i2e12 + (
                        -6019881./655360.)*i2e10 + (-20829./4096.)*i2e8 + (-2809./512.)*i2e6 + (
                        1055142483./209715200.)*e14 + (30355211./7864320.)*e12 + (2006627./655360.)*e10 + (
                        6943./4096.)*e8 + (2809./1536.)*e6,
            None,
            3.,
            0.,
            0.
        ),
        (
            '4*n',
            '4*n',
            '2, 0, 1, 4',
            4*orbital_freq,
            (-290521./46080.)*i6e8 + (43043./10240.)*i4e10 + (77077./6144.)*i4e8 + (-708871./51200.)*i2e12 + (
                        -9933./2560.)*i2e10 + (-5929./512.)*i2e8 + (3027367./691200.)*e14 + (708871./153600.)*e12 + (
                        3311./2560.)*e10 + (5929./1536.)*e8,
            None,
            4.,
            0.,
            0.
        ),
        (
            '5*n',
            '5*n',
            '2, 0, 1, 5',
            5*orbital_freq,
            (13621959./524288.)*i4e10 + (2947317./524288.)*i2e12 + (-3143529./131072.)*i2e10 + (
                        13524196093./1585446912.)*e14 + (-982439./524288.)*e12 + (1047843./131072.)*e10,
            None,
            5.,
            0.,
            0.
        ),
        (
            '6*n',
            '6*n',
            '2, 0, 1, 6',
            6*orbital_freq,
            (-10029889./204800.)*i2e12 + (-2308743./179200.)*e14 + (10029889./614400.)*e12,
            None,
            6.,
            0.,
            0.
        ),
        (
            '7*n',
            '7*n',
            '2, 0, 1, 7',
            7*orbital_freq,
            (186702632281./5662310400.)*e14,
            None,
            7.,
            0.,
            0.
        ),
        (
            '-7*n',
            '7*n',
            '2, 0, 2, -5',
            -7*orbital_freq,
            (52142352409./157286400.)*i4e10,
            None,
            -7.,
            -2.,
            0.
        ),
        (
            '-6*n',
            '6*n',
            '2, 0, 2, -4',
            -6*orbital_freq,
            (-284089./4096.)*i6e8 + (-22109373./40960.)*i4e10 + (852267./8192.)*i4e8,
            None,
            -6.,
            -2.,
            0.
        ),
        (
            '-5*n',
            '5*n',
            '2, 0, 2, -3',
            -5*orbital_freq,
            (142805./24576.)*i8e6 + (27483625./294912.)*i6e8 + (-714025./36864.)*i6e6 + (587225375./2097152.)*i4e10 + (
                        -27483625./196608.)*i4e8 + (714025./24576.)*i4e6,
            None,
            -5.,
            -2.,
            0.
        ),
        (
            '-4*n',
            '4*n',
            '2, 0, 2, -2',
            -4*orbital_freq,
            (-4913./20160.)*i10e4 + (-391./64.)*i8e6 + (867./640.)*i8e4 + (-83551./2304.)*i6e8 + (1955./96.)*i6e6 + (
                        -289./64.)*i6e4 + (-134209./2560.)*i4e10 + (83551./1536.)*i4e8 + (-1955./64.)*i4e6 + (
                        867./128.)*i4e4,
            None,
            -4.,
            -2.,
            0.
        ),
        (
            '-3*n',
            '3*n',
            '2, 0, 2, -1',
            -3*orbital_freq,
            (217./43200.)*i12e2 + (697./3840.)*i10e4 + (-119./2880.)*i10e2 + (13185./8192.)*i8e6 + (
                        -2583./2560.)*i8e4 + (147./640.)*i8e2 + (135771./32768.)*i6e8 + (-21975./4096.)*i6e6 + (
                        861./256.)*i6e4 + (-49./64.)*i6e2 + (16704927./5242880.)*i4e10 + (-407313./65536.)*i4e8 + (
                        65925./8192.)*i4e6 + (-2583./512.)*i4e4 + (147./128.)*i4e2,
            None,
            -3.,
            -2.,
            0.
        ),
        (
            '-2*n',
            '2*n',
            '2, 0, 2, 0',
            -2*orbital_freq,
            (-1./27720.)*i14 + (-31./15120.)*i12e2 + (31./75600.)*i12 + (-17./640.)*i10e4 + (17./1008.)*i10e2 + (
                        -17./5040.)*i10 + (-19./240.)*i8e6 + (189./1280.)*i8e4 + (-3./32.)*i8e2 + (3./160.)*i8 + (
                        -2065./36864.)*i6e8 + (19./72.)*i6e6 + (-63./128.)*i6e4 + (5./16.)*i6e2 + (-1./16.)*i6 + (
                        4079./204800.)*i4e10 + (2065./24576.)*i4e8 + (-19./48.)*i4e6 + (189./256.)*i4e4 + (
                        -15./32.)*i4e2 + (3./32.)*i4,
            None,
            -2.,
            -2.,
            0.
        ),
        (
            '-n',
            'n',
            '2, 0, 2, 1',
            -orbital_freq,
            (31./302400.)*i12e2 + (17./80640.)*i10e4 + (-17./20160.)*i10e2 + (13./40960.)*i8e6 + (-3./2560.)*i8e4 + (
                        3./640.)*i8e2 + (-305./294912.)*i6e8 + (-13./12288.)*i6e6 + (1./256.)*i6e4 + (-1./64.)*i6e2 + (
                        5199./5242880.)*i4e10 + (305./196608.)*i4e8 + (13./8192.)*i4e6 + (-3./512.)*i4e4 + (
                        3./128.)*i4e2,
            None,
            -1.,
            -2.,
            0.
        ),
        (
            'n',
            'n',
            '2, 0, 2, 3',
            orbital_freq,
            (1./122880.)*i8e6 + (-11./294912.)*i6e8 + (-1./36864.)*i6e6 + (619./10485760.)*i4e10 + (
                        11./196608.)*i4e8 + (1./24576.)*i4e6,
            None,
            1.,
            -2.,
            0.
        ),
        (
            '2*n',
            '2*n',
            '2, 0, 2, 4',
            2*orbital_freq,
            (-1./9216.)*i6e8 + (7./30720.)*i4e10 + (1./6144.)*i4e8,
            None,
            2.,
            -2.,
            0.
        ),
        (
            '3*n',
            '3*n',
            '2, 0, 2, 5',
            3*orbital_freq,
            (19683./52428800.)*i4e10,
            None,
            3.,
            -2.,
            0.
        ),
        (
            '-O - 4*n',
            'O + 4*n',
            '2, 1, 0, -6',
            -spin_freq - 4*orbital_freq,
            (8./2025.)*i2e12,
            None,
            -4.,
            2.,
            1.
        ),
        (
            '-O - 3*n',
            'O + 3*n',
            '2, 1, 0, -5',
            -spin_freq - 3*orbital_freq,
            (-2187./1310720.)*i4e10 + (6561./2621440.)*i2e12 + (6561./3276800.)*i2e10,
            None,
            -3.,
            2.,
            1.
        ),
        (
            '-O - 2*n',
            'O + 2*n',
            '2, 1, 0, -4',
            -spin_freq - 2*orbital_freq,
            (227./829440.)*i6e8 + (-7./6912.)*i4e10 + (-5./6912.)*i4e8 + (949./691200.)*i2e12 + (7./5760.)*i2e10 + (
                        1./1152.)*i2e8,
            None,
            -2.,
            2.,
            1.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 0, -3',
            -spin_freq - orbital_freq,
            (-145./9289728.)*i8e6 + (2497./26542080.)*i6e8 + (227./3317760.)*i6e6 + (-619./2359296.)*i4e10 + (
                        -55./221184.)*i4e8 + (-5./27648.)*i4e6 + (62617./212336640.)*i2e12 + (619./1966080.)*i2e10 + (
                        11./36864.)*i2e8 + (1./4608.)*i2e6,
            None,
            -1.,
            2.,
            1.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 0, -1',
            -spin_freq + orbital_freq,
            (-59123./383201280.)*i12e2 + (-40277./116121600.)*i10e4 + (40277./29030400.)*i10e2 + (
                        -1885./3096576.)*i8e6 + (145./64512.)*i8e4 + (-145./16128.)*i8e2 + (13847./5308416.)*i6e8 + (
                        2951./1105920.)*i6e6 + (-227./23040.)*i6e4 + (227./5760.)*i6e2 + (-1733./393216.)*i4e10 + (
                        -1525./221184.)*i4e8 + (-65./9216.)*i4e6 + (5./192.)*i4e4 + (-5./48.)*i4e2 + (
                        277229./58982400.)*i2e12 + (1733./327680.)*i2e10 + (305./36864.)*i2e8 + (13./1536.)*i2e6 + (
                        -1./32.)*i2e4 + (1./8.)*i2e2,
            None,
            1.,
            2.,
            1.
        ),
        (
            '-O + 2*n',
            'O - 2*n',
            '2, 1, 0, 0',
            -spin_freq + 2*orbital_freq,
            (8988527./174356582400.)*i14 + (59123./19160064.)*i12e2 + (-59123./95800320.)*i12 + (
                        40277./921600.)*i10e4 + (-40277./1451520.)*i10e2 + (40277./7257600.)*i10 + (
                        11455./72576.)*i8e6 + (-145./512.)*i8e4 + (725./4032.)*i8e2 + (-145./4032.)*i8 + (
                        839219./3317760.)*i6e8 + (-17933./25920.)*i6e6 + (1589./1280.)*i6e4 + (-227./288.)*i6e2 + (
                        227./1440.)*i6 + (11041./46080.)*i4e10 + (-18485./27648.)*i4e8 + (395./216.)*i4e6 + (
                        -105./32.)*i4e4 + (25./12.)*i4e2 + (-5./12.)*i4 + (34471./552960.)*i2e12 + (
                        -11041./38400.)*i2e10 + (3697./4608.)*i2e8 + (-79./36.)*i2e6 + (63./16.)*i2e4 + (
                        -5./2.)*i2e2 + (1./2.)*i2,
            None,
            2.,
            2.,
            1.
        ),
        (
            '-O + 3*n',
            'O - 3*n',
            '2, 1, 0, 1',
            -spin_freq + 3*orbital_freq,
            (-413861./54743040.)*i12e2 + (-1651357./5529600.)*i10e4 + (281939./4147200.)*i10e2 + (
                        -1062125./344064.)*i8e6 + (5945./3072.)*i8e4 + (-1015./2304.)*i8e2 + (
                        -10273339./983040.)*i6e8 + (332555./24576.)*i6e6 + (-65149./7680.)*i6e4 + (
                        11123./5760.)*i6e2 + (-1856103./131072.)*i4e10 + (226285./8192.)*i4e8 + (-36625./1024.)*i4e6 + (
                        1435./64.)*i4e4 + (-245./48.)*i4e2 + (-47982879./6553600.)*i2e12 + (5568309./327680.)*i2e10 + (
                        -135771./4096.)*i2e8 + (21975./512.)*i2e6 + (-861./32.)*i2e4 + (49./8.)*i2e2,
            None,
            3.,
            2.,
            1.
        ),
        (
            '-O + 4*n',
            'O - 4*n',
            '2, 1, 0, 2',
            -spin_freq + 4*orbital_freq,
            (11640053./29030400.)*i10e4 + (283475./24192.)*i8e6 + (-41905./16128.)*i8e4 + (18966077./207360.)*i6e8 + (
                        -88757./1728.)*i6e6 + (65603./5760.)*i6e4 + (134209./576.)*i4e10 + (-417755./1728.)*i4e8 + (
                        9775./72.)*i4e6 + (-1445./48.)*i4e4 + (8421731./46080.)*i2e12 + (-134209./480.)*i2e10 + (
                        83551./288.)*i2e8 + (-1955./12.)*i2e6 + (289./8.)*i2e4,
            None,
            4.,
            2.,
            1.
        ),
        (
            '-O + 5*n',
            'O - 5*n',
            '2, 1, 0, 3',
            -spin_freq + 5*orbital_freq,
            (-103533625./9289728.)*i8e6 + (-1247756575./5308416.)*i6e8 + (32416735./663552.)*i6e6 + (
                        -2936126875./2359296.)*i4e10 + (137418125./221184.)*i4e8 + (-3570125./27648.)*i4e6 + (
                        -72670996375./42467328.)*i2e12 + (587225375./393216.)*i2e10 + (-27483625./36864.)*i2e8 + (
                        714025./4608.)*i2e6,
            None,
            5.,
            2.,
            1.
        ),
        (
            '-O + 6*n',
            'O - 6*n',
            '2, 1, 0, 4',
            -spin_freq + 6*orbital_freq,
            (64488203./368640.)*i6e8 + (2456597./1024.)*i4e10 + (-1420445./3072.)*i4e8 + (659870313./102400.)*i2e12 + (
                        -7369791./2560.)*i2e10 + (284089./512.)*i2e8,
            None,
            6.,
            2.,
            1.
        ),
        (
            '-O + 7*n',
            'O - 7*n',
            '2, 1, 0, 5',
            -spin_freq + 7*orbital_freq,
            (-52142352409./35389440.)*i4e10 + (-140254152605./14155776.)*i2e12 + (52142352409./29491200.)*i2e10,
            None,
            7.,
            2.,
            1.
        ),
        (
            '-O + 8*n',
            'O - 8*n',
            '2, 1, 0, 6',
            -spin_freq + 8*orbital_freq,
            (5383010161./1036800.)*i2e12,
            None,
            8.,
            2.,
            1.
        ),
        (
            '-O - 6*n',
            'O + 6*n',
            '2, 1, 1, -6',
            -spin_freq - 6*orbital_freq,
            (10029889./204800.)*i2e12,
            None,
            -6.,
            0.,
            1.
        ),
        (
            '-O - 5*n',
            'O + 5*n',
            '2, 1, 1, -5',
            -spin_freq - 5*orbital_freq,
            (-1047843./32768.)*i4e10 + (-2947317./524288.)*i2e12 + (3143529./131072.)*i2e10,
            None,
            -5.,
            0.,
            1.
        ),
        (
            '-O - 4*n',
            'O + 4*n',
            '2, 1, 1, -4',
            -spin_freq - 4*orbital_freq,
            (5929./720.)*i6e8 + (-3311./640.)*i4e10 + (-5929./384.)*i4e8 + (708871./51200.)*i2e12 + (
                        9933./2560.)*i2e10 + (5929./512.)*i2e8,
            None,
            -4.,
            0.,
            1.
        ),
        (
            '-O - 3*n',
            'O + 3*n',
            '2, 1, 1, -3',
            -spin_freq - 3*orbital_freq,
            (-2809./2520.)*i8e6 + (6943./1920.)*i6e8 + (2809./720.)*i6e6 + (-2006627./163840.)*i4e10 + (
                        -6943./1024.)*i4e8 + (-2809./384.)*i4e6 + (30355211./2621440.)*i2e12 + (
                        6019881./655360.)*i2e10 + (20829./4096.)*i2e8 + (2809./512.)*i2e6,
            None,
            -3.,
            0.,
            1.
        ),
        (
            '-O - 2*n',
            'O + 2*n',
            '2, 1, 1, -2',
            -spin_freq - 2*orbital_freq,
            (16./175.)*i10e4 + (-4./5.)*i8e6 + (-18./35.)*i8e4 + (1661./360.)*i6e8 + (14./5.)*i6e6 + (9./5.)*i6e4 + (
                        -3919./320.)*i4e10 + (-1661./192.)*i4e8 + (-21./4.)*i4e6 + (-27./8.)*i4e4 + (
                        505601./40960.)*i2e12 + (11757./1280.)*i2e10 + (1661./256.)*i2e8 + (63./16.)*i2e6 + (
                        81./32.)*i2e4,
            None,
            -2.,
            0.,
            1.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 1, -1',
            -spin_freq - orbital_freq,
            (-256./51975.)*i12e2 + (16./175.)*i10e4 + (64./1575.)*i10e2 + (-51./56.)*i8e6 + (-18./35.)*i8e4 + (
                        -8./35.)*i8e2 + (28211./5760.)*i6e8 + (51./16.)*i6e6 + (9./5.)*i6e4 + (4./5.)*i6e2 + (
                        -1064887./81920.)*i4e10 + (-28211./3072.)*i4e8 + (-765./128.)*i4e6 + (-27./8.)*i4e4 + (
                        -3./2.)*i4e2 + (85553819./6553600.)*i2e12 + (3194661./327680.)*i2e10 + (28211./4096.)*i2e8 + (
                        2295./512.)*i2e6 + (81./32.)*i2e4 + (9./8.)*i2e2,
            None,
            -1.,
            0.,
            1.
        ),
        (
            '-O',
            'O',
            '2, 1, 1, 0',
            -spin_freq,
            (8192./42567525.)*i14 + (-1024./155925.)*i12e2 + (-1024./467775.)*i12 + (512./4725.)*i10e4 + (
                        256./4725.)*i10e2 + (256./14175.)*i10 + (-64./63.)*i8e6 + (-64./105.)*i8e4 + (
                        -32./105.)*i8e2 + (-32./315.)*i8 + (16./3.)*i6e8 + (32./9.)*i6e6 + (32./15.)*i6e4 + (
                        16./15.)*i6e2 + (16./45.)*i6 + -14.*i4e10 + -10.*i4e8 + (-20./3.)*i4e6 + -4.*i4e4 + -2.*i4e2 + (
                        -2./3.)*i4 + 14.*i2e12 + (21./2.)*i2e10 + (15./2.)*i2e8 + 5.*i2e6 + 3.*i2e4 + (3./2.)*i2e2 + (
                        1./2.)*i2,
            None,
            0.,
            0.,
            1.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 1, 1',
            -spin_freq + orbital_freq,
            (-256./51975.)*i12e2 + (16./175.)*i10e4 + (64./1575.)*i10e2 + (-51./56.)*i8e6 + (-18./35.)*i8e4 + (
                        -8./35.)*i8e2 + (28211./5760.)*i6e8 + (51./16.)*i6e6 + (9./5.)*i6e4 + (4./5.)*i6e2 + (
                        -1064887./81920.)*i4e10 + (-28211./3072.)*i4e8 + (-765./128.)*i4e6 + (-27./8.)*i4e4 + (
                        -3./2.)*i4e2 + (85553819./6553600.)*i2e12 + (3194661./327680.)*i2e10 + (28211./4096.)*i2e8 + (
                        2295./512.)*i2e6 + (81./32.)*i2e4 + (9./8.)*i2e2,
            None,
            1.,
            0.,
            1.
        ),
        (
            '-O + 2*n',
            'O - 2*n',
            '2, 1, 1, 2',
            -spin_freq + 2*orbital_freq,
            (16./175.)*i10e4 + (-4./5.)*i8e6 + (-18./35.)*i8e4 + (1661./360.)*i6e8 + (14./5.)*i6e6 + (9./5.)*i6e4 + (
                        -3919./320.)*i4e10 + (-1661./192.)*i4e8 + (-21./4.)*i4e6 + (-27./8.)*i4e4 + (
                        505601./40960.)*i2e12 + (11757./1280.)*i2e10 + (1661./256.)*i2e8 + (63./16.)*i2e6 + (
                        81./32.)*i2e4,
            None,
            2.,
            0.,
            1.
        ),
        (
            '-O + 3*n',
            'O - 3*n',
            '2, 1, 1, 3',
            -spin_freq + 3*orbital_freq,
            (-2809./2520.)*i8e6 + (6943./1920.)*i6e8 + (2809./720.)*i6e6 + (-2006627./163840.)*i4e10 + (
                        -6943./1024.)*i4e8 + (-2809./384.)*i4e6 + (30355211./2621440.)*i2e12 + (
                        6019881./655360.)*i2e10 + (20829./4096.)*i2e8 + (2809./512.)*i2e6,
            None,
            3.,
            0.,
            1.
        ),
        (
            '-O + 4*n',
            'O - 4*n',
            '2, 1, 1, 4',
            -spin_freq + 4*orbital_freq,
            (5929./720.)*i6e8 + (-3311./640.)*i4e10 + (-5929./384.)*i4e8 + (708871./51200.)*i2e12 + (
                        9933./2560.)*i2e10 + (5929./512.)*i2e8,
            None,
            4.,
            0.,
            1.
        ),
        (
            '-O + 5*n',
            'O - 5*n',
            '2, 1, 1, 5',
            -spin_freq + 5*orbital_freq,
            (-1047843./32768.)*i4e10 + (-2947317./524288.)*i2e12 + (3143529./131072.)*i2e10,
            None,
            5.,
            0.,
            1.
        ),
        (
            '-O + 6*n',
            'O - 6*n',
            '2, 1, 1, 6',
            -spin_freq + 6*orbital_freq,
            (10029889./204800.)*i2e12,
            None,
            6.,
            0.,
            1.
        ),
        (
            '-O - 6*n',
            'O + 6*n',
            '2, 1, 2, -4',
            -spin_freq - 6*orbital_freq,
            (284089./8192.)*i6e8,
            None,
            -6.,
            -2.,
            1.
        ),
        (
            '-O - 5*n',
            'O + 5*n',
            '2, 1, 2, -3',
            -spin_freq - 5*orbital_freq,
            (-714025./147456.)*i8e6 + (-27483625./589824.)*i6e8 + (714025./73728.)*i6e6,
            None,
            -5.,
            -2.,
            1.
        ),
        (
            '-O - 4*n',
            'O + 4*n',
            '2, 1, 2, -2',
            -spin_freq - 4*orbital_freq,
            (2601./10240.)*i10e4 + (1955./384.)*i8e6 + (-289./256.)*i8e4 + (83551./4608.)*i6e8 + (-1955./192.)*i6e6 + (
                        289./128.)*i6e4,
            None,
            -4.,
            -2.,
            1.
        ),
        (
            '-O - 3*n',
            'O + 3*n',
            '2, 1, 2, -1',
            -spin_freq - 3*orbital_freq,
            (-3241./552960.)*i12e2 + (-7749./40960.)*i10e4 + (441./10240.)*i10e2 + (-21975./16384.)*i8e6 + (
                        861./1024.)*i8e4 + (-49./256.)*i8e2 + (-135771./65536.)*i6e8 + (21975./8192.)*i6e6 + (
                        -861./512.)*i6e4 + (49./128.)*i6e2,
            None,
            -3.,
            -2.,
            1.
        ),
        (
            '-O - 2*n',
            'O + 2*n',
            '2, 1, 2, 0',
            -spin_freq - 2*orbital_freq,
            (173./3870720.)*i14 + (463./193536.)*i12e2 + (-463./967680.)*i12 + (567./20480.)*i10e4 + (
                        -9./512.)*i10e2 + (9./2560.)*i10 + (19./288.)*i8e6 + (-63./512.)*i8e4 + (5./64.)*i8e2 + (
                        -1./64.)*i8 + (2065./73728.)*i6e8 + (-19./144.)*i6e6 + (63./256.)*i6e4 + (-5./32.)*i6e2 + (
                        1./32.)*i6,
            None,
            -2.,
            -2.,
            1.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 2, 1',
            -spin_freq - orbital_freq,
            (-463./3870720.)*i12e2 + (-9./40960.)*i10e4 + (9./10240.)*i10e2 + (-13./49152.)*i8e6 + (1./1024.)*i8e4 + (
                        -1./256.)*i8e2 + (305./589824.)*i6e8 + (13./24576.)*i6e6 + (-1./512.)*i6e4 + (1./128.)*i6e2,
            None,
            -1.,
            -2.,
            1.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 2, 3',
            -spin_freq + orbital_freq,
            (-1./147456.)*i8e6 + (11./589824.)*i6e8 + (1./73728.)*i6e6,
            None,
            1.,
            -2.,
            1.
        ),
        (
            '-O + 2*n',
            'O - 2*n',
            '2, 1, 2, 4',
            -spin_freq + 2*orbital_freq,
            (1./18432.)*i6e8,
            None,
            2.,
            -2.,
            1.
        ),
        (
            '-2*O - 5*n',
            '2*O + 5*n',
            '2, 2, 0, -7',
            -2*spin_freq - 5*orbital_freq,
            (244140625./33294385152.)*e14,
            None,
            -5.,
            2.,
            2.
        ),
        (
            '-2*O - 4*n',
            '2*O + 4*n',
            '2, 2, 0, -6',
            -2*spin_freq - 4*orbital_freq,
            (-8./2025.)*i2e12 + (8./2025.)*e14 + (8./2025.)*e12,
            None,
            -4.,
            2.,
            2.
        ),
        (
            '-2*O - 3*n',
            '2*O + 3*n',
            '2, 2, 0, -5',
            -2*spin_freq - 3*orbital_freq,
            (24057./26214400.)*i4e10 + (-6561./2621440.)*i2e12 + (-6561./3276800.)*i2e10 + (
                        4284333./1468006400.)*e14 + (6561./2621440.)*e12 + (6561./3276800.)*e10,
            None,
            -3.,
            2.,
            2.
        ),
        (
            '-2*O - 2*n',
            '2*O + 2*n',
            '2, 2, 0, -4',
            -2*spin_freq - 2*orbital_freq,
            (-23./207360.)*i6e8 + (77./138240.)*i4e10 + (11./27648.)*i4e8 + (-949./691200.)*i2e12 + (
                        -7./5760.)*i2e10 + (-1./1152.)*i2e8 + (2417./1814400.)*e14 + (949./691200.)*e12 + (
                        7./5760.)*e10 + (1./1152.)*e8,
            None,
            -2.,
            2.,
            2.
        ),
        (
            '-2*O - n',
            '2*O + n',
            '2, 2, 0, -3',
            -2*spin_freq - orbital_freq,
            (1957./371589120.)*i8e6 + (-253./6635520.)*i6e8 + (-23./829440.)*i6e6 + (6809./47185920.)*i4e10 + (
                        121./884736.)*i4e8 + (11./110592.)*i4e6 + (-62617./212336640.)*i2e12 + (
                        -619./1966080.)*i2e10 + (-11./36864.)*i2e8 + (-1./4608.)*i2e6 + (
                        31398887./118908518400.)*e14 + (62617./212336640.)*e12 + (619./1966080.)*e10 + (
                        11./36864.)*e8 + (1./4608.)*e6,
            None,
            -1.,
            2.,
            2.
        ),
        (
            '-2*O + n',
            '2*O - n',
            '2, 2, 0, -1',
            -2*spin_freq + orbital_freq,
            (330367./7664025600.)*i12e2 + (12107./116121600.)*i10e4 + (-12107./29030400.)*i10e2 + (
                        25441./123863040.)*i8e6 + (-1957./2580480.)*i8e4 + (1957./645120.)*i8e2 + (
                        -1403./1327104.)*i6e8 + (-299./276480.)*i6e6 + (23./5760.)*i6e4 + (-23./1440.)*i6e2 + (
                        19063./7864320.)*i4e10 + (3355./884736.)*i4e8 + (143./36864.)*i4e6 + (-11./768.)*i4e4 + (
                        11./192.)*i4e2 + (-277229./58982400.)*i2e12 + (-1733./327680.)*i2e10 + (-305./36864.)*i2e8 + (
                        -13./1536.)*i2e6 + (1./32.)*i2e4 + (-1./8.)*i2e2 + (165285343./39636172800.)*e14 + (
                        277229./58982400.)*e12 + (1733./327680.)*e10 + (305./36864.)*e8 + (13./1536.)*e6 + (
                        -1./32.)*e4 + (1./8.)*e2,
            None,
            1.,
            2.,
            2.
        ),
        (
            '-2*O + 2*n',
            '2*O - 2*n',
            '2, 2, 0, 0',
            -2*spin_freq + 2*orbital_freq,
            (-27269./1981324800.)*i14 + (-330367./383201280.)*i12e2 + (330367./1916006400.)*i12 + (
                        -12107./921600.)*i10e4 + (12107./1451520.)*i10e2 + (-12107./7257600.)*i10 + (
                        -154603./2903040.)*i8e6 + (1957./20480.)*i8e4 + (-1957./32256.)*i8e2 + (1957./161280.)*i8 + (
                        -85031./829440.)*i6e8 + (1817./6480.)*i6e6 + (-161./320.)*i6e4 + (23./72.)*i6e2 + (
                        -23./360.)*i6 + (-121451./921600.)*i4e10 + (40667./110592.)*i4e8 + (-869./864.)*i4e6 + (
                        231./128.)*i4e4 + (-55./48.)*i4e2 + (11./48.)*i4 + (-34471./552960.)*i2e12 + (
                        11041./38400.)*i2e10 + (-3697./4608.)*i2e8 + (79./36.)*i2e6 + (-63./16.)*i2e4 + (5./2.)*i2e2 + (
                        -1./2.)*i2 + (-1739939./67737600.)*e14 + (34471./552960.)*e12 + (-11041./38400.)*e10 + (
                        3697./4608.)*e8 + (-79./36.)*e6 + (63./16.)*e4 + (-5./2.)*e2 + (1./2.),
            None,
            2.,
            2.,
            2.
        ),
        (
            '-2*O + 3*n',
            '2*O - 3*n',
            '2, 2, 0, 1',
            -2*spin_freq + 3*orbital_freq,
            (2312569./1094860800.)*i12e2 + (496387./5529600.)*i10e4 + (-84749./4147200.)*i10e2 + (
                        2867005./2752512.)*i8e6 + (-80237./122880.)*i8e4 + (13699./92160.)*i8e2 + (
                        1040911./245760.)*i6e8 + (-33695./6144.)*i6e6 + (6601./1920.)*i6e4 + (-1127./1440.)*i6e2 + (
                        20417133./2621440.)*i4e10 + (-497827./32768.)*i4e8 + (80575./4096.)*i4e6 + (
                        -3157./256.)*i4e4 + (539./192.)*i4e2 + (47982879./6553600.)*i2e12 + (
                        -5568309./327680.)*i2e10 + (135771./4096.)*i2e8 + (-21975./512.)*i2e6 + (861./32.)*i2e4 + (
                        -49./8.)*i2e2 + (534226163./209715200.)*e14 + (-47982879./6553600.)*e12 + (
                        5568309./327680.)*e10 + (-135771./4096.)*e8 + (21975./512.)*e6 + (-861./32.)*e4 + (49./8.)*e2,
            None,
            3.,
            2.,
            2.
        ),
        (
            '-2*O + 4*n',
            '2*O - 4*n',
            '2, 2, 0, 2',
            -2*spin_freq + 4*orbital_freq,
            (-3498923./29030400.)*i10e4 + (-765187./193536.)*i8e6 + (565573./645120.)*i8e4 + (-1921673./51840.)*i6e8 + (
                        8993./432.)*i6e6 + (-6647./1440.)*i6e4 + (-1476299./11520.)*i4e10 + (919061./6912.)*i4e8 + (
                        -21505./288.)*i4e6 + (3179./192.)*i4e4 + (-8421731./46080.)*i2e12 + (134209./480.)*i2e10 + (
                        -83551./288.)*i2e8 + (1955./12.)*i2e6 + (-289./8.)*i2e4 + (-346700573./3628800.)*e14 + (
                        8421731./46080.)*e12 + (-134209./480.)*e10 + (83551./288.)*e8 + (-1955./12.)*e6 + (289./8.)*e4,
            None,
            4.,
            2.,
            2.
        ),
        (
            '-2*O + 5*n',
            '2*O - 5*n',
            '2, 2, 0, 3',
            -2*spin_freq + 5*orbital_freq,
            (279469385./74317824.)*i8e6 + (126424675./1327104.)*i6e8 + (-3284515./165888.)*i6e6 + (
                        6459479125./9437184.)*i4e10 + (-302319875./884736.)*i4e8 + (7854275./110592.)*i4e6 + (
                        72670996375./42467328.)*i2e12 + (-587225375./393216.)*i2e10 + (27483625./36864.)*i2e8 + (
                        -714025./4608.)*i2e6 + (6429415592375./4756340736.)*e14 + (-72670996375./42467328.)*e12 + (
                        587225375./393216.)*e10 + (-27483625./36864.)*e8 + (714025./4608.)*e6,
            None,
            5.,
            2.,
            2.
        ),
        (
            '-2*O + 6*n',
            '2*O - 6*n',
            '2, 2, 0, 4',
            -2*spin_freq + 6*orbital_freq,
            (-6534047./92160.)*i6e8 + (-27022567./20480.)*i4e10 + (3124979./12288.)*i4e8 + (
                        -659870313./102400.)*i2e12 + (7369791./2560.)*i2e10 + (-284089./512.)*i2e8 + (
                        -3055721757./358400.)*e14 + (659870313./102400.)*e12 + (-7369791./2560.)*e10 + (
                        284089./512.)*e8,
            None,
            6.,
            2.,
            2.
        ),
        (
            '-2*O + 7*n',
            '2*O - 7*n',
            '2, 2, 0, 5',
            -2*spin_freq + 7*orbital_freq,
            (573565876499./707788800.)*i4e10 + (140254152605./14155776.)*i2e12 + (-52142352409./29491200.)*i2e10 + (
                        417304029320899./16986931200.)*e14 + (-140254152605./14155776.)*e12 + (
                        52142352409./29491200.)*e10,
            None,
            7.,
            2.,
            2.
        ),
        (
            '-2*O + 8*n',
            '2*O - 8*n',
            '2, 2, 0, 6',
            -2*spin_freq + 8*orbital_freq,
            (-5383010161./1036800.)*i2e12 + (-56914314263./1814400.)*e14 + (5383010161./1036800.)*e12,
            None,
            8.,
            2.,
            2.
        ),
        (
            '-2*O + 9*n',
            '2*O - 9*n',
            '2, 2, 0, 7',
            -2*spin_freq + 9*orbital_freq,
            (147483366698529./10276044800.)*e14,
            None,
            9.,
            2.,
            2.
        ),
        (
            '-2*O - 5*n',
            '2*O + 5*n',
            '2, 2, 1, -5',
            -2*spin_freq - 5*orbital_freq,
            (3143529./524288.)*i4e10,
            None,
            -5.,
            0.,
            2.
        ),
        (
            '-2*O - 4*n',
            '2*O + 4*n',
            '2, 2, 1, -4',
            -2*spin_freq - 4*orbital_freq,
            (-5929./3072.)*i6e8 + (9933./10240.)*i4e10 + (5929./2048.)*i4e8,
            None,
            -4.,
            0.,
            2.
        ),
        (
            '-2*O - 3*n',
            '2*O + 3*n',
            '2, 2, 1, -3',
            -2*spin_freq - 3*orbital_freq,
            (2809./10240.)*i8e6 + (-6943./8192.)*i6e8 + (-2809./3072.)*i6e6 + (6019881./2621440.)*i4e10 + (
                        20829./16384.)*i4e8 + (2809./2048.)*i4e6,
            None,
            -3.,
            0.,
            2.
        ),
        (
            '-2*O - 2*n',
            '2*O + 2*n',
            '2, 2, 1, -2',
            -2*spin_freq - 2*orbital_freq,
            (-51./2240.)*i10e4 + (63./320.)*i8e6 + (81./640.)*i8e4 + (-1661./1536.)*i6e8 + (-21./32.)*i6e6 + (
                        -27./64.)*i6e4 + (11757./5120.)*i4e10 + (1661./1024.)*i4e8 + (63./64.)*i4e6 + (81./128.)*i4e4,
            None,
            -2.,
            0.,
            2.
        ),
        (
            '-2*O - n',
            '2*O + n',
            '2, 2, 1, -1',
            -2*spin_freq - orbital_freq,
            (31./25200.)*i12e2 + (-51./2240.)*i10e4 + (-17./1680.)*i10e2 + (459./2048.)*i8e6 + (81./640.)*i8e4 + (
                        9./160.)*i8e2 + (-28211./24576.)*i6e8 + (-765./1024.)*i6e6 + (-27./64.)*i6e4 + (
                        -3./16.)*i6e2 + (3194661./1310720.)*i4e10 + (28211./16384.)*i4e8 + (2295./2048.)*i4e6 + (
                        81./128.)*i4e4 + (9./32.)*i4e2,
            None,
            -1.,
            0.,
            2.
        ),
        (
            '-2*O',
            '2*O',
            '2, 2, 1, 0',
            -2*spin_freq,
            (-1./20790.)*i14 + (31./18900.)*i12e2 + (31./56700.)*i12 + (-17./630.)*i10e4 + (-17./1260.)*i10e2 + (
                        -17./3780.)*i10 + (1./4.)*i8e6 + (3./20.)*i8e4 + (3./40.)*i8e2 + (1./40.)*i8 + (-5./4.)*i6e8 + (
                        -5./6.)*i6e6 + (-1./2.)*i6e4 + (-1./4.)*i6e2 + (-1./12.)*i6 + (21./8.)*i4e10 + (15./8.)*i4e8 + (
                        5./4.)*i4e6 + (3./4.)*i4e4 + (3./8.)*i4e2 + (1./8.)*i4,
            None,
            0.,
            0.,
            2.
        ),
        (
            '-2*O + n',
            '2*O - n',
            '2, 2, 1, 1',
            -2*spin_freq + orbital_freq,
            (31./25200.)*i12e2 + (-51./2240.)*i10e4 + (-17./1680.)*i10e2 + (459./2048.)*i8e6 + (81./640.)*i8e4 + (
                        9./160.)*i8e2 + (-28211./24576.)*i6e8 + (-765./1024.)*i6e6 + (-27./64.)*i6e4 + (
                        -3./16.)*i6e2 + (3194661./1310720.)*i4e10 + (28211./16384.)*i4e8 + (2295./2048.)*i4e6 + (
                        81./128.)*i4e4 + (9./32.)*i4e2,
            None,
            1.,
            0.,
            2.
        ),
        (
            '-2*O + 2*n',
            '2*O - 2*n',
            '2, 2, 1, 2',
            -2*spin_freq + 2*orbital_freq,
            (-51./2240.)*i10e4 + (63./320.)*i8e6 + (81./640.)*i8e4 + (-1661./1536.)*i6e8 + (-21./32.)*i6e6 + (
                        -27./64.)*i6e4 + (11757./5120.)*i4e10 + (1661./1024.)*i4e8 + (63./64.)*i4e6 + (81./128.)*i4e4,
            None,
            2.,
            0.,
            2.
        ),
        (
            '-2*O + 3*n',
            '2*O - 3*n',
            '2, 2, 1, 3',
            -2*spin_freq + 3*orbital_freq,
            (2809./10240.)*i8e6 + (-6943./8192.)*i6e8 + (-2809./3072.)*i6e6 + (6019881./2621440.)*i4e10 + (
                        20829./16384.)*i4e8 + (2809./2048.)*i4e6,
            None,
            3.,
            0.,
            2.
        ),
        (
            '-2*O + 4*n',
            '2*O - 4*n',
            '2, 2, 1, 4',
            -2*spin_freq + 4*orbital_freq,
            (-5929./3072.)*i6e8 + (9933./10240.)*i4e10 + (5929./2048.)*i4e8,
            None,
            4.,
            0.,
            2.
        ),
        (
            '-2*O + 5*n',
            '2*O - 5*n',
            '2, 2, 1, 5',
            -2*spin_freq + 5*orbital_freq,
            (3143529./524288.)*i4e10,
            None,
            5.,
            0.,
            2.
        ),
        (
            '-2*O - 5*n',
            '2*O + 5*n',
            '2, 2, 2, -3',
            -2*spin_freq - 5*orbital_freq,
            (714025./1179648.)*i8e6,
            None,
            -5.,
            -2.,
            2.
        ),
        (
            '-2*O - 4*n',
            '2*O + 4*n',
            '2, 2, 2, -2',
            -2*spin_freq - 4*orbital_freq,
            (-289./6144.)*i10e4 + (-1955./3072.)*i8e6 + (289./2048.)*i8e4,
            None,
            -4.,
            -2.,
            2.
        ),
        (
            '-2*O - 3*n',
            '2*O + 3*n',
            '2, 2, 2, -1',
            -2*spin_freq - 3*orbital_freq,
            (931./737280.)*i12e2 + (287./8192.)*i10e4 + (-49./6144.)*i10e2 + (21975./131072.)*i8e6 + (
                        -861./8192.)*i8e4 + (49./2048.)*i8e2,
            None,
            -3.,
            -2.,
            2.
        ),
        (
            '-2*O - 2*n',
            '2*O + 2*n',
            '2, 2, 2, 0',
            -2*spin_freq - 2*orbital_freq,
            (-1./96768.)*i14 + (-19./36864.)*i12e2 + (19./184320.)*i12 + (-21./4096.)*i10e4 + (5./1536.)*i10e2 + (
                        -1./1536.)*i10 + (-19./2304.)*i8e6 + (63./4096.)*i8e4 + (-5./512.)*i8e2 + (1./512.)*i8,
            None,
            -2.,
            -2.,
            2.
        ),
        (
            '-2*O - n',
            '2*O + n',
            '2, 2, 2, 1',
            -2*spin_freq - orbital_freq,
            (19./737280.)*i12e2 + (1./24576.)*i10e4 + (-1./6144.)*i10e2 + (13./393216.)*i8e6 + (-1./8192.)*i8e4 + (
                        1./2048.)*i8e2,
            None,
            -1.,
            -2.,
            2.
        ),
        (
            '-2*O + n',
            '2*O - n',
            '2, 2, 2, 3',
            -2*spin_freq + orbital_freq,
            (1./1179648.)*i8e6,
            None,
            1.,
            -2.,
            2.
        )
    )

    return mode_data_output

@njit
def nsr_modes_t16(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray) \
        -> OutputType:
    """ Non-sync Tidal Modes for truncation level 16, tidal order 2

    Parameters
    ----------
    orbital_freq : np.ndarray
        Planet's orbital mean motion [rad s-1]
    spin_freq : np.ndarray
        Planet's spin frequency [rad s-1]
    eccentricity : np.ndarray
        Planet's eccentricity
    inclination : np.ndarray
        Planet's inclination [rads]

    Returns
    -------
    mode_data_tuple : OutputType
        List of mode datas. Each mode data contains information to calculate the heating and torques.
    """
    # Tidal Potential Builder for l-order = 2, NSR = True.
    # Max Eccentricity Order = 16
    # Max Inclination Order = 16
    # Max q = 9.
    # Number of unique modes = 50.
    # Number of unique frequencies = 42.

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

    mode_data_output = (
        (
            '-4*n',
            '4*n',
            '2, 0, 0, -6',
            -4*orbital_freq,
            (1./1350.)*i4e12,
            None,
            -4.,
            2.,
            0.
        ),
        (
            '-3*n',
            '3*n',
            '2, 0, 0, -5',
            -3*orbital_freq,
            (-6561./26214400.)*i6e10 + (19683./41943040.)*i4e12 + (19683./52428800.)*i4e10,
            None,
            -3.,
            2.,
            0.
        ),
        (
            '-2*n',
            '2*n',
            '2, 0, 0, -4',
            -2*orbital_freq,
            (1./30720.)*i8e8 + (-7./46080.)*i6e10 + (-1./9216.)*i6e8 + (949./3686400.)*i4e12 + (7./30720.)*i4e10 + (
                        1./6144.)*i4e8,
            None,
            -2.,
            2.,
            0.
        ),
        (
            '-n',
            'n',
            '2, 0, 0, -3',
            -orbital_freq,
            (-17./11612160.)*i10e6 + (11./983040.)*i8e8 + (1./122880.)*i8e6 + (-619./15728640.)*i6e10 + (
                        -11./294912.)*i6e8 + (-1./36864.)*i6e6 + (62617./1132462080.)*i4e12 + (619./10485760.)*i4e10 + (
                        11./196608.)*i4e8 + (1./24576.)*i4e6,
            None,
            -1.,
            2.,
            0.
        ),
        (
            'n',
            'n',
            '2, 0, 0, -1',
            orbital_freq,
            (-1./110880.)*i14e2 + (-31./1209600.)*i12e4 + (31./302400.)*i12e2 + (-221./3870720.)*i10e6 + (
                        17./80640.)*i10e4 + (-17./20160.)*i10e2 + (61./196608.)*i8e8 + (13./40960.)*i8e6 + (
                        -3./2560.)*i8e4 + (3./640.)*i8e2 + (-1733./2621440.)*i6e10 + (-305./294912.)*i6e8 + (
                        -13./12288.)*i6e6 + (1./256.)*i6e4 + (-1./64.)*i6e2 + (277229./314572800.)*i4e12 + (
                        5199./5242880.)*i4e10 + (305./196608.)*i4e8 + (13./8192.)*i4e6 + (-3./512.)*i4e4 + (
                        3./128.)*i4e2,
            None,
            1.,
            2.,
            0.
        ),
        (
            '2*n',
            '2*n',
            '2, 0, 0, 0',
            2*orbital_freq,
            (5461./2270268000.)*i16 + (1./5544.)*i14e2 + (-1./27720.)*i14 + (31./9600.)*i12e4 + (-31./15120.)*i12e2 + (
                        31./75600.)*i12 + (1343./90720.)*i10e6 + (-17./640.)*i10e4 + (17./1008.)*i10e2 + (
                        -17./5040.)*i10 + (3697./122880.)*i8e8 + (-79./960.)*i8e6 + (189./1280.)*i8e4 + (
                        -3./32.)*i8e2 + (3./160.)*i8 + (11041./307200.)*i6e10 + (-3697./36864.)*i6e8 + (
                        79./288.)*i6e6 + (-63./128.)*i6e4 + (5./16.)*i6e2 + (-1./16.)*i6 + (34471./2949120.)*i4e12 + (
                        -11041./204800.)*i4e10 + (3697./24576.)*i4e8 + (-79./192.)*i4e6 + (189./256.)*i4e4 + (
                        -15./32.)*i4e2 + (3./32.)*i4,
            None,
            2.,
            2.,
            0.
        ),
        (
            '3*n',
            '3*n',
            '2, 0, 0, 1',
            3*orbital_freq,
            (-7./15840.)*i14e2 + (-1271./57600.)*i12e4 + (217./43200.)*i12e2 + (-24905./86016.)*i10e6 + (
                        697./3840.)*i10e4 + (-119./2880.)*i10e2 + (-407313./327680.)*i8e8 + (13185./8192.)*i8e6 + (
                        -2583./2560.)*i8e4 + (147./640.)*i8e2 + (-5568309./2621440.)*i6e10 + (135771./32768.)*i6e8 + (
                        -21975./4096.)*i6e6 + (861./256.)*i6e4 + (-49./64.)*i6e2 + (-143948637./104857600.)*i4e12 + (
                        16704927./5242880.)*i4e10 + (-407313./65536.)*i4e8 + (65925./8192.)*i4e6 + (
                        -2583./512.)*i4e4 + (147./128.)*i4e2,
            None,
            3.,
            2.,
            0.
        ),
        (
            '4*n',
            '4*n',
            '2, 0, 0, 2',
            4*orbital_freq,
            (8959./302400.)*i12e4 + (6647./6048.)*i10e6 + (-4913./20160.)*i10e4 + (83551./7680.)*i8e8 + (
                        -391./64.)*i8e6 + (867./640.)*i8e4 + (134209./3840.)*i6e10 + (-83551./2304.)*i6e8 + (
                        1955./96.)*i6e6 + (-289./64.)*i6e4 + (8421731./245760.)*i4e12 + (-134209./2560.)*i4e10 + (
                        83551./1536.)*i4e8 + (-1955./64.)*i4e6 + (867./128.)*i4e4,
            None,
            4.,
            2.,
            0.
        ),
        (
            '5*n',
            '5*n',
            '2, 0, 0, 3',
            5*orbital_freq,
            (-2427685./2322432.)*i10e6 + (-5496725./196608.)*i8e8 + (142805./24576.)*i8e6 + (
                        -587225375./3145728.)*i6e10 + (27483625./294912.)*i6e8 + (-714025./36864.)*i6e6 + (
                        -72670996375./226492416.)*i4e12 + (587225375./2097152.)*i4e10 + (-27483625./196608.)*i4e8 + (
                        714025./24576.)*i4e6,
            None,
            5.,
            2.,
            0.
        ),
        (
            '6*n',
            '6*n',
            '2, 0, 0, 4',
            6*orbital_freq,
            (852267./40960.)*i8e8 + (7369791./20480.)*i6e10 + (-284089./4096.)*i6e8 + (1979610939./1638400.)*i4e12 + (
                        -22109373./40960.)*i4e10 + (852267./8192.)*i4e8,
            None,
            6.,
            2.,
            0.
        ),
        (
            '7*n',
            '7*n',
            '2, 0, 0, 5',
            7*orbital_freq,
            (-52142352409./235929600.)*i6e10 + (-140254152605./75497472.)*i4e12 + (52142352409./157286400.)*i4e10,
            None,
            7.,
            2.,
            0.
        ),
        (
            '8*n',
            '8*n',
            '2, 0, 0, 6',
            8*orbital_freq,
            (5383010161./5529600.)*i4e12,
            None,
            8.,
            2.,
            0.
        ),
        (
            '-8*n',
            '8*n',
            '2, 0, 1, -8',
            -8*orbital_freq,
            (31801945561./481689600.)*e16,
            None,
            -8.,
            0.,
            0.
        ),
        (
            '-7*n',
            '7*n',
            '2, 0, 1, -7',
            -7*orbital_freq,
            (-186702632281./1887436800.)*i2e14 + (-1328179895713./30198988800.)*e16 + (186702632281./5662310400.)*e14,
            None,
            -7.,
            0.,
            0.
        ),
        (
            '-6*n',
            '6*n',
            '2, 0, 1, -6',
            -6*orbital_freq,
            (130388557./2457600.)*i4e12 + (6926229./179200.)*i2e14 + (-10029889./204800.)*i2e12 + (
                        65845785./3211264.)*e16 + (-2308743./179200.)*e14 + (10029889./614400.)*e12,
            None,
            -6.,
            0.,
            0.
        ),
        (
            '-5*n',
            '5*n',
            '2, 0, 1, -5',
            -5*orbital_freq,
            (-17114769./1310720.)*i6e10 + (-12771707./2097152.)*i4e12 + (13621959./524288.)*i4e10 + (
                        -13524196093./528482304.)*i2e14 + (2947317./524288.)*i2e12 + (-3143529./131072.)*i2e10 + (
                        28491858875./8455716864.)*e16 + (13524196093./1585446912.)*e14 + (-982439./524288.)*e12 + (
                        1047843./131072.)*e10,
            None,
            -5.,
            0.,
            0.
        ),
        (
            '-4*n',
            '4*n',
            '2, 0, 1, -4',
            -4*orbital_freq,
            (163471./92160.)*i8e8 + (-162239./76800.)*i6e10 + (-290521./46080.)*i6e8 + (9215323./614400.)*i4e12 + (
                        43043./10240.)*i4e10 + (77077./6144.)*i4e8 + (-3027367./230400.)*i2e14 + (
                        -708871./51200.)*i2e12 + (-9933./2560.)*i2e10 + (-5929./512.)*i2e8 + (
                        124914751./20643840.)*e16 + (3027367./691200.)*e14 + (708871./153600.)*e12 + (
                        3311./2560.)*e10 + (5929./1536.)*e8,
            None,
            -4.,
            0.,
            0.
        ),
        (
            '-3*n',
            '3*n',
            '2, 0, 1, -3',
            -3*orbital_freq,
            (-2160121./14515200.)*i10e6 + (1339999./1720320.)*i8e8 + (542137./645120.)*i8e6 + (
                        -98324723./19660800.)*i6e10 + (-340207./122880.)*i6e8 + (-137641./46080.)*i6e6 + (
                        394617743./31457280.)*i4e12 + (26086151./2621440.)*i4e10 + (90259./16384.)*i4e8 + (
                        36517./6144.)*i4e6 + (-3165427449./209715200.)*i2e14 + (-30355211./2621440.)*i2e12 + (
                        -6019881./655360.)*i2e10 + (-20829./4096.)*i2e8 + (-2809./512.)*i2e6 + (
                        145747579353./23488102400.)*e16 + (1055142483./209715200.)*e14 + (30355211./7864320.)*e12 + (
                        2006627./655360.)*e10 + (6943./4096.)*e8 + (2809./1536.)*e6,
            None,
            -3.,
            0.,
            0.
        ),
        (
            '-2*n',
            '2*n',
            '2, 0, 1, -2',
            -2*orbital_freq,
            (439./52800.)*i12e4 + (-769./7200.)*i10e6 + (-769./11200.)*i10e4 + (320573./322560.)*i8e8 + (
                        193./320.)*i8e6 + (1737./4480.)*i8e4 + (-192031./38400.)*i6e10 + (-81389./23040.)*i6e8 + (
                        -343./160.)*i6e6 + (-441./320.)*i6e4 + (6572813./491520.)*i4e12 + (50947./5120.)*i4e10 + (
                        21593./3072.)*i4e8 + (273./64.)*i4e6 + (351./128.)*i4e4 + (-25565893./1612800.)*i2e14 + (
                        -505601./40960.)*i2e12 + (-11757./1280.)*i2e10 + (-1661./256.)*i2e8 + (-63./16.)*i2e6 + (
                        -81./32.)*i2e4 + (678544541./103219200.)*e16 + (25565893./4838400.)*e14 + (
                        505601./122880.)*e12 + (3919./1280.)*e10 + (1661./768.)*e8 + (21./16.)*e6 + (27./32.)*e4,
            None,
            -2.,
            0.,
            0.
        ),
        (
            '-n',
            'n',
            '2, 0, 1, -1',
            -orbital_freq,
            (-12289./37837800.)*i14e2 + (439./52800.)*i12e4 + (439./118800.)*i12e2 + (-13073./107520.)*i10e6 + (
                        -769./11200.)*i10e4 + (-769./25200.)*i10e2 + (5444723./5160960.)*i8e8 + (9843./14336.)*i8e6 + (
                        1737./4480.)*i8e4 + (193./1120.)*i8e2 + (-52179463./9830400.)*i6e10 + (
                        -1382339./368640.)*i6e8 + (-2499./1024.)*i6e6 + (-441./320.)*i6e4 + (-49./80.)*i6e2 + (
                        1112199647./78643200.)*i4e12 + (13843531./1310720.)*i4e10 + (366743./49152.)*i4e8 + (
                        9945./2048.)*i4e6 + (351./128.)*i4e4 + (39./32.)*i4e2 + (-24654653741./1468006400.)*i2e14 + (
                        -85553819./6553600.)*i2e12 + (-3194661./327680.)*i2e10 + (-28211./4096.)*i2e8 + (
                        -2295./512.)*i2e6 + (-81./32.)*i2e4 + (-9./8.)*i2e2 + (689312857627./98650030080.)*e16 + (
                        24654653741./4404019200.)*e14 + (85553819./19660800.)*e12 + (1064887./327680.)*e10 + (
                        28211./12288.)*e8 + (765./512.)*e6 + (27./32.)*e4 + (3./8.)*e2,
            None,
            -1.,
            0.,
            0.
        ),
        (
            'n',
            'n',
            '2, 0, 1, 1',
            orbital_freq,
            (-12289./37837800.)*i14e2 + (439./52800.)*i12e4 + (439./118800.)*i12e2 + (-13073./107520.)*i10e6 + (
                        -769./11200.)*i10e4 + (-769./25200.)*i10e2 + (5444723./5160960.)*i8e8 + (9843./14336.)*i8e6 + (
                        1737./4480.)*i8e4 + (193./1120.)*i8e2 + (-52179463./9830400.)*i6e10 + (
                        -1382339./368640.)*i6e8 + (-2499./1024.)*i6e6 + (-441./320.)*i6e4 + (-49./80.)*i6e2 + (
                        1112199647./78643200.)*i4e12 + (13843531./1310720.)*i4e10 + (366743./49152.)*i4e8 + (
                        9945./2048.)*i4e6 + (351./128.)*i4e4 + (39./32.)*i4e2 + (-24654653741./1468006400.)*i2e14 + (
                        -85553819./6553600.)*i2e12 + (-3194661./327680.)*i2e10 + (-28211./4096.)*i2e8 + (
                        -2295./512.)*i2e6 + (-81./32.)*i2e4 + (-9./8.)*i2e2 + (689312857627./98650030080.)*e16 + (
                        24654653741./4404019200.)*e14 + (85553819./19660800.)*e12 + (1064887./327680.)*e10 + (
                        28211./12288.)*e8 + (765./512.)*e6 + (27./32.)*e4 + (3./8.)*e2,
            None,
            1.,
            0.,
            0.
        ),
        (
            '2*n',
            '2*n',
            '2, 0, 1, 2',
            2*orbital_freq,
            (439./52800.)*i12e4 + (-769./7200.)*i10e6 + (-769./11200.)*i10e4 + (320573./322560.)*i8e8 + (
                        193./320.)*i8e6 + (1737./4480.)*i8e4 + (-192031./38400.)*i6e10 + (-81389./23040.)*i6e8 + (
                        -343./160.)*i6e6 + (-441./320.)*i6e4 + (6572813./491520.)*i4e12 + (50947./5120.)*i4e10 + (
                        21593./3072.)*i4e8 + (273./64.)*i4e6 + (351./128.)*i4e4 + (-25565893./1612800.)*i2e14 + (
                        -505601./40960.)*i2e12 + (-11757./1280.)*i2e10 + (-1661./256.)*i2e8 + (-63./16.)*i2e6 + (
                        -81./32.)*i2e4 + (678544541./103219200.)*e16 + (25565893./4838400.)*e14 + (
                        505601./122880.)*e12 + (3919./1280.)*e10 + (1661./768.)*e8 + (21./16.)*e6 + (27./32.)*e4,
            None,
            2.,
            0.,
            0.
        ),
        (
            '3*n',
            '3*n',
            '2, 0, 1, 3',
            3*orbital_freq,
            (-2160121./14515200.)*i10e6 + (1339999./1720320.)*i8e8 + (542137./645120.)*i8e6 + (
                        -98324723./19660800.)*i6e10 + (-340207./122880.)*i6e8 + (-137641./46080.)*i6e6 + (
                        394617743./31457280.)*i4e12 + (26086151./2621440.)*i4e10 + (90259./16384.)*i4e8 + (
                        36517./6144.)*i4e6 + (-3165427449./209715200.)*i2e14 + (-30355211./2621440.)*i2e12 + (
                        -6019881./655360.)*i2e10 + (-20829./4096.)*i2e8 + (-2809./512.)*i2e6 + (
                        145747579353./23488102400.)*e16 + (1055142483./209715200.)*e14 + (30355211./7864320.)*e12 + (
                        2006627./655360.)*e10 + (6943./4096.)*e8 + (2809./1536.)*e6,
            None,
            3.,
            0.,
            0.
        ),
        (
            '4*n',
            '4*n',
            '2, 0, 1, 4',
            4*orbital_freq,
            (163471./92160.)*i8e8 + (-162239./76800.)*i6e10 + (-290521./46080.)*i6e8 + (9215323./614400.)*i4e12 + (
                        43043./10240.)*i4e10 + (77077./6144.)*i4e8 + (-3027367./230400.)*i2e14 + (
                        -708871./51200.)*i2e12 + (-9933./2560.)*i2e10 + (-5929./512.)*i2e8 + (
                        124914751./20643840.)*e16 + (3027367./691200.)*e14 + (708871./153600.)*e12 + (
                        3311./2560.)*e10 + (5929./1536.)*e8,
            None,
            4.,
            0.,
            0.
        ),
        (
            '5*n',
            '5*n',
            '2, 0, 1, 5',
            5*orbital_freq,
            (-17114769./1310720.)*i6e10 + (-12771707./2097152.)*i4e12 + (13621959./524288.)*i4e10 + (
                        -13524196093./528482304.)*i2e14 + (2947317./524288.)*i2e12 + (-3143529./131072.)*i2e10 + (
                        28491858875./8455716864.)*e16 + (13524196093./1585446912.)*e14 + (-982439./524288.)*e12 + (
                        1047843./131072.)*e10,
            None,
            5.,
            0.,
            0.
        ),
        (
            '6*n',
            '6*n',
            '2, 0, 1, 6',
            6*orbital_freq,
            (130388557./2457600.)*i4e12 + (6926229./179200.)*i2e14 + (-10029889./204800.)*i2e12 + (
                        65845785./3211264.)*e16 + (-2308743./179200.)*e14 + (10029889./614400.)*e12,
            None,
            6.,
            0.,
            0.
        ),
        (
            '7*n',
            '7*n',
            '2, 0, 1, 7',
            7*orbital_freq,
            (-186702632281./1887436800.)*i2e14 + (-1328179895713./30198988800.)*e16 + (186702632281./5662310400.)*e14,
            None,
            7.,
            0.,
            0.
        ),
        (
            '8*n',
            '8*n',
            '2, 0, 1, 8',
            8*orbital_freq,
            (31801945561./481689600.)*e16,
            None,
            8.,
            0.,
            0.
        ),
        (
            '-8*n',
            '8*n',
            '2, 0, 2, -6',
            -8*orbital_freq,
            (5383010161./5529600.)*i4e12,
            None,
            -8.,
            -2.,
            0.
        ),
        (
            '-7*n',
            '7*n',
            '2, 0, 2, -5',
            -7*orbital_freq,
            (-52142352409./235929600.)*i6e10 + (-140254152605./75497472.)*i4e12 + (52142352409./157286400.)*i4e10,
            None,
            -7.,
            -2.,
            0.
        ),
        (
            '-6*n',
            '6*n',
            '2, 0, 2, -4',
            -6*orbital_freq,
            (852267./40960.)*i8e8 + (7369791./20480.)*i6e10 + (-284089./4096.)*i6e8 + (1979610939./1638400.)*i4e12 + (
                        -22109373./40960.)*i4e10 + (852267./8192.)*i4e8,
            None,
            -6.,
            -2.,
            0.
        ),
        (
            '-5*n',
            '5*n',
            '2, 0, 2, -3',
            -5*orbital_freq,
            (-2427685./2322432.)*i10e6 + (-5496725./196608.)*i8e8 + (142805./24576.)*i8e6 + (
                        -587225375./3145728.)*i6e10 + (27483625./294912.)*i6e8 + (-714025./36864.)*i6e6 + (
                        -72670996375./226492416.)*i4e12 + (587225375./2097152.)*i4e10 + (-27483625./196608.)*i4e8 + (
                        714025./24576.)*i4e6,
            None,
            -5.,
            -2.,
            0.
        ),
        (
            '-4*n',
            '4*n',
            '2, 0, 2, -2',
            -4*orbital_freq,
            (8959./302400.)*i12e4 + (6647./6048.)*i10e6 + (-4913./20160.)*i10e4 + (83551./7680.)*i8e8 + (
                        -391./64.)*i8e6 + (867./640.)*i8e4 + (134209./3840.)*i6e10 + (-83551./2304.)*i6e8 + (
                        1955./96.)*i6e6 + (-289./64.)*i6e4 + (8421731./245760.)*i4e12 + (-134209./2560.)*i4e10 + (
                        83551./1536.)*i4e8 + (-1955./64.)*i4e6 + (867./128.)*i4e4,
            None,
            -4.,
            -2.,
            0.
        ),
        (
            '-3*n',
            '3*n',
            '2, 0, 2, -1',
            -3*orbital_freq,
            (-7./15840.)*i14e2 + (-1271./57600.)*i12e4 + (217./43200.)*i12e2 + (-24905./86016.)*i10e6 + (
                        697./3840.)*i10e4 + (-119./2880.)*i10e2 + (-407313./327680.)*i8e8 + (13185./8192.)*i8e6 + (
                        -2583./2560.)*i8e4 + (147./640.)*i8e2 + (-5568309./2621440.)*i6e10 + (135771./32768.)*i6e8 + (
                        -21975./4096.)*i6e6 + (861./256.)*i6e4 + (-49./64.)*i6e2 + (-143948637./104857600.)*i4e12 + (
                        16704927./5242880.)*i4e10 + (-407313./65536.)*i4e8 + (65925./8192.)*i4e6 + (
                        -2583./512.)*i4e4 + (147./128.)*i4e2,
            None,
            -3.,
            -2.,
            0.
        ),
        (
            '-2*n',
            '2*n',
            '2, 0, 2, 0',
            -2*orbital_freq,
            (5461./2270268000.)*i16 + (1./5544.)*i14e2 + (-1./27720.)*i14 + (31./9600.)*i12e4 + (-31./15120.)*i12e2 + (
                        31./75600.)*i12 + (323./22680.)*i10e6 + (-17./640.)*i10e4 + (17./1008.)*i10e2 + (
                        -17./5040.)*i10 + (413./24576.)*i8e8 + (-19./240.)*i8e6 + (189./1280.)*i8e4 + (-3./32.)*i8e2 + (
                        3./160.)*i8 + (-4079./307200.)*i6e10 + (-2065./36864.)*i6e8 + (19./72.)*i6e6 + (
                        -63./128.)*i6e4 + (5./16.)*i6e2 + (-1./16.)*i6 + (-2323./983040.)*i4e12 + (
                        4079./204800.)*i4e10 + (2065./24576.)*i4e8 + (-19./48.)*i4e6 + (189./256.)*i4e4 + (
                        -15./32.)*i4e2 + (3./32.)*i4,
            None,
            -2.,
            -2.,
            0.
        ),
        (
            '-n',
            'n',
            '2, 0, 2, 1',
            -orbital_freq,
            (-1./110880.)*i14e2 + (-31./1209600.)*i12e4 + (31./302400.)*i12e2 + (-221./3870720.)*i10e6 + (
                        17./80640.)*i10e4 + (-17./20160.)*i10e2 + (61./196608.)*i8e8 + (13./40960.)*i8e6 + (
                        -3./2560.)*i8e4 + (3./640.)*i8e2 + (-1733./2621440.)*i6e10 + (-305./294912.)*i6e8 + (
                        -13./12288.)*i6e6 + (1./256.)*i6e4 + (-1./64.)*i6e2 + (277229./314572800.)*i4e12 + (
                        5199./5242880.)*i4e10 + (305./196608.)*i4e8 + (13./8192.)*i4e6 + (-3./512.)*i4e4 + (
                        3./128.)*i4e2,
            None,
            -1.,
            -2.,
            0.
        ),
        (
            'n',
            'n',
            '2, 0, 2, 3',
            orbital_freq,
            (-17./11612160.)*i10e6 + (11./983040.)*i8e8 + (1./122880.)*i8e6 + (-619./15728640.)*i6e10 + (
                        -11./294912.)*i6e8 + (-1./36864.)*i6e6 + (62617./1132462080.)*i4e12 + (619./10485760.)*i4e10 + (
                        11./196608.)*i4e8 + (1./24576.)*i4e6,
            None,
            1.,
            -2.,
            0.
        ),
        (
            '2*n',
            '2*n',
            '2, 0, 2, 4',
            2*orbital_freq,
            (1./30720.)*i8e8 + (-7./46080.)*i6e10 + (-1./9216.)*i6e8 + (949./3686400.)*i4e12 + (7./30720.)*i4e10 + (
                        1./6144.)*i4e8,
            None,
            2.,
            -2.,
            0.
        ),
        (
            '3*n',
            '3*n',
            '2, 0, 2, 5',
            3*orbital_freq,
            (-6561./26214400.)*i6e10 + (19683./41943040.)*i4e12 + (19683./52428800.)*i4e10,
            None,
            3.,
            -2.,
            0.
        ),
        (
            '4*n',
            '4*n',
            '2, 0, 2, 6',
            4*orbital_freq,
            (1./1350.)*i4e12,
            None,
            4.,
            -2.,
            0.
        ),
        (
            '-O - 5*n',
            'O + 5*n',
            '2, 1, 0, -7',
            -spin_freq - 5*orbital_freq,
            (244140625./33294385152.)*i2e14,
            None,
            -5.,
            2.,
            1.
        ),
        (
            '-O - 4*n',
            'O + 4*n',
            '2, 1, 0, -6',
            -spin_freq - 4*orbital_freq,
            (-4./1215.)*i4e12 + (8./2025.)*i2e14 + (8./2025.)*i2e12,
            None,
            -4.,
            2.,
            1.
        ),
        (
            '-O - 3*n',
            'O + 3*n',
            '2, 1, 0, -5',
            -spin_freq - 3*orbital_freq,
            (165483./262144000.)*i6e10 + (-2187./1048576.)*i4e12 + (-2187./1310720.)*i4e10 + (
                        4284333./1468006400.)*i2e14 + (6561./2621440.)*i2e12 + (6561./3276800.)*i2e10,
            None,
            -3.,
            2.,
            1.
        ),
        (
            '-O - 2*n',
            'O + 2*n',
            '2, 1, 0, -4',
            -spin_freq - 2*orbital_freq,
            (-145./2322432.)*i8e8 + (1589./4147200.)*i6e10 + (227./829440.)*i6e8 + (-949./829440.)*i4e12 + (
                        -7./6912.)*i4e10 + (-5./6912.)*i4e8 + (2417./1814400.)*i2e14 + (949./691200.)*i2e12 + (
                        7./5760.)*i2e10 + (1./1152.)*i2e8,
            None,
            -2.,
            2.,
            1.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 0, -3',
            -spin_freq - orbital_freq,
            (40277./16721510400.)*i10e6 + (-1595./74317824.)*i8e8 + (-145./9289728.)*i8e6 + (
                        140513./1415577600.)*i6e10 + (2497./26542080.)*i6e8 + (227./3317760.)*i6e6 + (
                        -62617./254803968.)*i4e12 + (-619./2359296.)*i4e10 + (-55./221184.)*i4e8 + (-5./27648.)*i4e6 + (
                        31398887./118908518400.)*i2e14 + (62617./212336640.)*i2e12 + (619./1966080.)*i2e10 + (
                        11./36864.)*i2e8 + (1./4608.)*i2e6,
            None,
            -1.,
            2.,
            1.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 0, -1',
            -spin_freq + orbital_freq,
            (8988527./697426329600.)*i14e2 + (59123./1532805120.)*i12e4 + (-59123./383201280.)*i12e2 + (
                        523601./5573836800.)*i10e6 + (-40277./116121600.)*i10e4 + (40277./29030400.)*i10e2 + (
                        -44225./74317824.)*i8e8 + (-1885./3096576.)*i8e6 + (145./64512.)*i8e4 + (-145./16128.)*i8e2 + (
                        393391./235929600.)*i6e10 + (13847./5308416.)*i6e8 + (2951./1105920.)*i6e6 + (
                        -227./23040.)*i6e4 + (227./5760.)*i6e2 + (-277229./70778880.)*i4e12 + (-1733./393216.)*i4e10 + (
                        -1525./221184.)*i4e8 + (-65./9216.)*i4e6 + (5./192.)*i4e4 + (-5./48.)*i4e2 + (
                        165285343./39636172800.)*i2e14 + (277229./58982400.)*i2e12 + (1733./327680.)*i2e10 + (
                        305./36864.)*i2e8 + (13./1536.)*i2e6 + (-1./32.)*i2e4 + (1./8.)*i2e2,
            None,
            1.,
            2.,
            1.
        ),
        (
            '-O + 2*n',
            'O - 2*n',
            '2, 1, 0, 0',
            -spin_freq + 2*orbital_freq,
            (-3490169./1046139494400.)*i16 + (-8988527./34871316480.)*i14e2 + (8988527./174356582400.)*i14 + (
                        -59123./12165120.)*i12e4 + (59123./19160064.)*i12e2 + (-59123./95800320.)*i12 + (
                        -3181883./130636800.)*i10e6 + (40277./921600.)*i10e4 + (-40277./1451520.)*i10e2 + (
                        40277./7257600.)*i10 + (-536065./9289728.)*i8e8 + (11455./72576.)*i8e6 + (-145./512.)*i8e4 + (
                        725./4032.)*i8e2 + (-145./4032.)*i8 + (-2506307./27648000.)*i6e10 + (839219./3317760.)*i6e8 + (
                        -17933./25920.)*i6e6 + (1589./1280.)*i6e4 + (-227./288.)*i6e2 + (227./1440.)*i6 + (
                        -34471./663552.)*i4e12 + (11041./46080.)*i4e10 + (-18485./27648.)*i4e8 + (395./216.)*i4e6 + (
                        -105./32.)*i4e4 + (25./12.)*i4e2 + (-5./12.)*i4 + (-1739939./67737600.)*i2e14 + (
                        34471./552960.)*i2e12 + (-11041./38400.)*i2e10 + (3697./4608.)*i2e8 + (-79./36.)*i2e6 + (
                        63./16.)*i2e4 + (-5./2.)*i2e2 + (1./2.)*i2,
            None,
            2.,
            2.,
            1.
        ),
        (
            '-O + 3*n',
            'O - 3*n',
            '2, 1, 0, 1',
            -spin_freq + 3*orbital_freq,
            (8988527./14233190400.)*i14e2 + (2424043./72990720.)*i12e4 + (-413861./54743040.)*i12e2 + (
                        11801161./24772608.)*i10e6 + (-1651357./5529600.)*i10e4 + (281939./4147200.)*i10e2 + (
                        6562265./2752512.)*i8e8 + (-1062125./344064.)*i8e6 + (5945./3072.)*i8e4 + (
                        -1015./2304.)*i8e2 + (140445127./26214400.)*i6e10 + (-10273339./983040.)*i6e8 + (
                        332555./24576.)*i6e6 + (-65149./7680.)*i6e4 + (11123./5760.)*i6e2 + (
                        15994293./2621440.)*i4e12 + (-1856103./131072.)*i4e10 + (226285./8192.)*i4e8 + (
                        -36625./1024.)*i4e6 + (1435./64.)*i4e4 + (-245./48.)*i4e2 + (534226163./209715200.)*i2e14 + (
                        -47982879./6553600.)*i2e12 + (5568309./327680.)*i2e10 + (-135771./4096.)*i2e8 + (
                        21975./512.)*i2e6 + (-861./32.)*i2e4 + (49./8.)*i2e2,
            None,
            3.,
            2.,
            1.
        ),
        (
            '-O + 4*n',
            'O - 4*n',
            '2, 1, 0, 2',
            -spin_freq + 4*orbital_freq,
            (-17086547./383201280.)*i12e4 + (-15748307./8709120.)*i10e6 + (11640053./29030400.)*i10e4 + (
                        -12114895./580608.)*i8e8 + (283475./24192.)*i8e6 + (-41905./16128.)*i8e4 + (
                        -30465443./345600.)*i6e10 + (18966077./207360.)*i6e8 + (-88757./1728.)*i6e6 + (
                        65603./5760.)*i6e4 + (-8421731./55296.)*i4e12 + (134209./576.)*i4e10 + (-417755./1728.)*i4e8 + (
                        9775./72.)*i4e6 + (-1445./48.)*i4e4 + (-346700573./3628800.)*i2e14 + (8421731./46080.)*i2e12 + (
                        -134209./480.)*i2e10 + (83551./288.)*i2e8 + (-1955./12.)*i2e6 + (289./8.)*i2e4,
            None,
            4.,
            2.,
            1.
        ),
        (
            '-O + 5*n',
            'O - 5*n',
            '2, 1, 0, 3',
            -spin_freq + 5*orbital_freq,
            (1150351397./668860416.)*i10e6 + (3985125625./74317824.)*i8e8 + (-103533625./9289728.)*i8e6 + (
                        26660032025./56623104.)*i6e10 + (-1247756575./5308416.)*i6e8 + (32416735./663552.)*i6e6 + (
                        363354981875./254803968.)*i4e12 + (-2936126875./2359296.)*i4e10 + (137418125./221184.)*i4e8 + (
                        -3570125./27648.)*i4e6 + (6429415592375./4756340736.)*i2e14 + (
                        -72670996375./42467328.)*i2e12 + (587225375./393216.)*i2e10 + (-27483625./36864.)*i2e8 + (
                        714025./4608.)*i2e6,
            None,
            5.,
            2.,
            1.
        ),
        (
            '-O + 6*n',
            'O - 6*n',
            '2, 1, 0, 4',
            -spin_freq + 6*orbital_freq,
            (-41192905./1032192.)*i8e8 + (-557647519./614400.)*i6e10 + (64488203./368640.)*i6e8 + (
                        -219956771./40960.)*i4e12 + (2456597./1024.)*i4e10 + (-1420445./3072.)*i4e8 + (
                        -3055721757./358400.)*i2e14 + (659870313./102400.)*i2e12 + (-7369791./2560.)*i2e10 + (
                        284089./512.)*i2e8,
            None,
            6.,
            2.,
            1.
        ),
        (
            '-O + 7*n',
            'O - 7*n',
            '2, 1, 0, 5',
            -spin_freq + 7*orbital_freq,
            (11836313996843./21233664000.)*i6e10 + (701270763025./84934656.)*i4e12 + (-52142352409./35389440.)*i4e10 + (
                        417304029320899./16986931200.)*i2e14 + (-140254152605./14155776.)*i2e12 + (
                        52142352409./29491200.)*i2e10,
            None,
            7.,
            2.,
            1.
        ),
        (
            '-O + 8*n',
            'O - 8*n',
            '2, 1, 0, 6',
            -spin_freq + 8*orbital_freq,
            (-5383010161./1244160.)*i4e12 + (-56914314263./1814400.)*i2e14 + (5383010161./1036800.)*i2e12,
            None,
            8.,
            2.,
            1.
        ),
        (
            '-O + 9*n',
            'O - 9*n',
            '2, 1, 0, 7',
            -spin_freq + 9*orbital_freq,
            (147483366698529./10276044800.)*i2e14,
            None,
            9.,
            2.,
            1.
        ),
        (
            '-O - 7*n',
            'O + 7*n',
            '2, 1, 1, -7',
            -spin_freq - 7*orbital_freq,
            (186702632281./1887436800.)*i2e14,
            None,
            -7.,
            0.,
            1.
        ),
        (
            '-O - 6*n',
            'O + 6*n',
            '2, 1, 1, -6',
            -spin_freq - 6*orbital_freq,
            (-10029889./153600.)*i4e12 + (-6926229./179200.)*i2e14 + (10029889./204800.)*i2e12,
            None,
            -6.,
            0.,
            1.
        ),
        (
            '-O - 5*n',
            'O + 5*n',
            '2, 1, 1, -5',
            -spin_freq - 5*orbital_freq,
            (349281./20480.)*i6e10 + (982439./131072.)*i4e12 + (-1047843./32768.)*i4e10 + (
                        13524196093./528482304.)*i2e14 + (-2947317./524288.)*i2e12 + (3143529./131072.)*i2e10,
            None,
            -5.,
            0.,
            1.
        ),
        (
            '-O - 4*n',
            'O + 4*n',
            '2, 1, 1, -4',
            -spin_freq - 4*orbital_freq,
            (-847./360.)*i8e8 + (3311./1200.)*i6e10 + (5929./720.)*i6e8 + (-708871./38400.)*i4e12 + (
                        -3311./640.)*i4e10 + (-5929./384.)*i4e8 + (3027367./230400.)*i2e14 + (708871./51200.)*i2e12 + (
                        9933./2560.)*i2e10 + (5929./512.)*i2e8,
            None,
            -4.,
            0.,
            1.
        ),
        (
            '-O - 3*n',
            'O + 3*n',
            '2, 1, 1, -3',
            -spin_freq - 3*orbital_freq,
            (2809./14175.)*i10e6 + (-6943./6720.)*i8e8 + (-2809./2520.)*i8e6 + (2006627./307200.)*i6e10 + (
                        6943./1920.)*i6e8 + (2809./720.)*i6e6 + (-30355211./1966080.)*i4e12 + (
                        -2006627./163840.)*i4e10 + (-6943./1024.)*i4e8 + (-2809./384.)*i4e6 + (
                        3165427449./209715200.)*i2e14 + (30355211./2621440.)*i2e12 + (6019881./655360.)*i2e10 + (
                        20829./4096.)*i2e8 + (2809./512.)*i2e6,
            None,
            -3.,
            0.,
            1.
        ),
        (
            '-O - 2*n',
            'O + 2*n',
            '2, 1, 1, -2',
            -spin_freq - 2*orbital_freq,
            (-64./5775.)*i12e4 + (32./225.)*i10e6 + (16./175.)*i10e4 + (-1661./1260.)*i8e8 + (-4./5.)*i8e6 + (
                        -18./35.)*i8e4 + (3919./600.)*i6e10 + (1661./360.)*i6e8 + (14./5.)*i6e6 + (9./5.)*i6e4 + (
                        -505601./30720.)*i4e12 + (-3919./320.)*i4e10 + (-1661./192.)*i4e8 + (-21./4.)*i4e6 + (
                        -27./8.)*i4e4 + (25565893./1612800.)*i2e14 + (505601./40960.)*i2e12 + (11757./1280.)*i2e10 + (
                        1661./256.)*i2e8 + (63./16.)*i2e6 + (81./32.)*i2e4,
            None,
            -2.,
            0.,
            1.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 1, -1',
            -spin_freq - orbital_freq,
            (2048./4729725.)*i14e2 + (-64./5775.)*i12e4 + (-256./51975.)*i12e2 + (17./105.)*i10e6 + (16./175.)*i10e4 + (
                        64./1575.)*i10e2 + (-28211./20160.)*i8e8 + (-51./56.)*i8e6 + (-18./35.)*i8e4 + (
                        -8./35.)*i8e2 + (1064887./153600.)*i6e10 + (28211./5760.)*i6e8 + (51./16.)*i6e6 + (
                        9./5.)*i6e4 + (4./5.)*i6e2 + (-85553819./4915200.)*i4e12 + (-1064887./81920.)*i4e10 + (
                        -28211./3072.)*i4e8 + (-765./128.)*i4e6 + (-27./8.)*i4e4 + (-3./2.)*i4e2 + (
                        24654653741./1468006400.)*i2e14 + (85553819./6553600.)*i2e12 + (3194661./327680.)*i2e10 + (
                        28211./4096.)*i2e8 + (2295./512.)*i2e6 + (81./32.)*i2e4 + (9./8.)*i2e2,
            None,
            -1.,
            0.,
            1.
        ),
        (
            '-O',
            'O',
            '2, 1, 1, 0',
            -spin_freq,
            (-8192./638512875.)*i16 + (8192./14189175.)*i14e2 + (8192./42567525.)*i14 + (-2048./155925.)*i12e4 + (
                        -1024./155925.)*i12e2 + (-1024./467775.)*i12 + (512./2835.)*i10e6 + (512./4725.)*i10e4 + (
                        256./4725.)*i10e2 + (256./14175.)*i10 + (-32./21.)*i8e8 + (-64./63.)*i8e6 + (-64./105.)*i8e4 + (
                        -32./105.)*i8e2 + (-32./315.)*i8 + (112./15.)*i6e10 + (16./3.)*i6e8 + (32./9.)*i6e6 + (
                        32./15.)*i6e4 + (16./15.)*i6e2 + (16./45.)*i6 + (-56./3.)*i4e12 + -14.*i4e10 + -10.*i4e8 + (
                        -20./3.)*i4e6 + -4.*i4e4 + -2.*i4e2 + (-2./3.)*i4 + 18.*i2e14 + 14.*i2e12 + (21./2.)*i2e10 + (
                        15./2.)*i2e8 + 5.*i2e6 + 3.*i2e4 + (3./2.)*i2e2 + (1./2.)*i2,
            None,
            0.,
            0.,
            1.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 1, 1',
            -spin_freq + orbital_freq,
            (2048./4729725.)*i14e2 + (-64./5775.)*i12e4 + (-256./51975.)*i12e2 + (17./105.)*i10e6 + (16./175.)*i10e4 + (
                        64./1575.)*i10e2 + (-28211./20160.)*i8e8 + (-51./56.)*i8e6 + (-18./35.)*i8e4 + (
                        -8./35.)*i8e2 + (1064887./153600.)*i6e10 + (28211./5760.)*i6e8 + (51./16.)*i6e6 + (
                        9./5.)*i6e4 + (4./5.)*i6e2 + (-85553819./4915200.)*i4e12 + (-1064887./81920.)*i4e10 + (
                        -28211./3072.)*i4e8 + (-765./128.)*i4e6 + (-27./8.)*i4e4 + (-3./2.)*i4e2 + (
                        24654653741./1468006400.)*i2e14 + (85553819./6553600.)*i2e12 + (3194661./327680.)*i2e10 + (
                        28211./4096.)*i2e8 + (2295./512.)*i2e6 + (81./32.)*i2e4 + (9./8.)*i2e2,
            None,
            1.,
            0.,
            1.
        ),
        (
            '-O + 2*n',
            'O - 2*n',
            '2, 1, 1, 2',
            -spin_freq + 2*orbital_freq,
            (-64./5775.)*i12e4 + (32./225.)*i10e6 + (16./175.)*i10e4 + (-1661./1260.)*i8e8 + (-4./5.)*i8e6 + (
                        -18./35.)*i8e4 + (3919./600.)*i6e10 + (1661./360.)*i6e8 + (14./5.)*i6e6 + (9./5.)*i6e4 + (
                        -505601./30720.)*i4e12 + (-3919./320.)*i4e10 + (-1661./192.)*i4e8 + (-21./4.)*i4e6 + (
                        -27./8.)*i4e4 + (25565893./1612800.)*i2e14 + (505601./40960.)*i2e12 + (11757./1280.)*i2e10 + (
                        1661./256.)*i2e8 + (63./16.)*i2e6 + (81./32.)*i2e4,
            None,
            2.,
            0.,
            1.
        ),
        (
            '-O + 3*n',
            'O - 3*n',
            '2, 1, 1, 3',
            -spin_freq + 3*orbital_freq,
            (2809./14175.)*i10e6 + (-6943./6720.)*i8e8 + (-2809./2520.)*i8e6 + (2006627./307200.)*i6e10 + (
                        6943./1920.)*i6e8 + (2809./720.)*i6e6 + (-30355211./1966080.)*i4e12 + (
                        -2006627./163840.)*i4e10 + (-6943./1024.)*i4e8 + (-2809./384.)*i4e6 + (
                        3165427449./209715200.)*i2e14 + (30355211./2621440.)*i2e12 + (6019881./655360.)*i2e10 + (
                        20829./4096.)*i2e8 + (2809./512.)*i2e6,
            None,
            3.,
            0.,
            1.
        ),
        (
            '-O + 4*n',
            'O - 4*n',
            '2, 1, 1, 4',
            -spin_freq + 4*orbital_freq,
            (-847./360.)*i8e8 + (3311./1200.)*i6e10 + (5929./720.)*i6e8 + (-708871./38400.)*i4e12 + (
                        -3311./640.)*i4e10 + (-5929./384.)*i4e8 + (3027367./230400.)*i2e14 + (708871./51200.)*i2e12 + (
                        9933./2560.)*i2e10 + (5929./512.)*i2e8,
            None,
            4.,
            0.,
            1.
        ),
        (
            '-O + 5*n',
            'O - 5*n',
            '2, 1, 1, 5',
            -spin_freq + 5*orbital_freq,
            (349281./20480.)*i6e10 + (982439./131072.)*i4e12 + (-1047843./32768.)*i4e10 + (
                        13524196093./528482304.)*i2e14 + (-2947317./524288.)*i2e12 + (3143529./131072.)*i2e10,
            None,
            5.,
            0.,
            1.
        ),
        (
            '-O + 6*n',
            'O - 6*n',
            '2, 1, 1, 6',
            -spin_freq + 6*orbital_freq,
            (-10029889./153600.)*i4e12 + (-6926229./179200.)*i2e14 + (10029889./204800.)*i2e12,
            None,
            6.,
            0.,
            1.
        ),
        (
            '-O + 7*n',
            'O - 7*n',
            '2, 1, 1, 7',
            -spin_freq + 7*orbital_freq,
            (186702632281./1887436800.)*i2e14,
            None,
            7.,
            0.,
            1.
        ),
        (
            '-O - 7*n',
            'O + 7*n',
            '2, 1, 2, -5',
            -spin_freq - 7*orbital_freq,
            (52142352409./471859200.)*i6e10,
            None,
            -7.,
            -2.,
            1.
        ),
        (
            '-O - 6*n',
            'O + 6*n',
            '2, 1, 2, -4',
            -spin_freq - 6*orbital_freq,
            (-284089./16384.)*i8e8 + (-7369791./40960.)*i6e10 + (284089./8192.)*i6e8,
            None,
            -6.,
            -2.,
            1.
        ),
        (
            '-O - 5*n',
            'O + 5*n',
            '2, 1, 2, -3',
            -spin_freq - 5*orbital_freq,
            (142805./131072.)*i10e6 + (27483625./1179648.)*i8e8 + (-714025./147456.)*i8e6 + (
                        587225375./6291456.)*i6e10 + (-27483625./589824.)*i6e8 + (714025./73728.)*i6e6,
            None,
            -5.,
            -2.,
            1.
        ),
        (
            '-O - 4*n',
            'O + 4*n',
            '2, 1, 2, -2',
            -spin_freq - 4*orbital_freq,
            (-133807./3870720.)*i12e4 + (-1173./1024.)*i10e6 + (2601./10240.)*i10e4 + (-83551./9216.)*i8e8 + (
                        1955./384.)*i8e6 + (-289./256.)*i8e4 + (-134209./7680.)*i6e10 + (83551./4608.)*i6e8 + (
                        -1955./192.)*i6e6 + (289./128.)*i6e4,
            None,
            -4.,
            -2.,
            1.
        ),
        (
            '-O - 3*n',
            'O + 3*n',
            '2, 1, 2, -1',
            -spin_freq - 3*orbital_freq,
            (1211./2211840.)*i14e2 + (18983./737280.)*i12e4 + (-3241./552960.)*i12e2 + (39555./131072.)*i10e6 + (
                        -7749./40960.)*i10e4 + (441./10240.)*i10e2 + (135771./131072.)*i8e8 + (-21975./16384.)*i8e6 + (
                        861./1024.)*i8e4 + (-49./256.)*i8e2 + (5568309./5242880.)*i6e10 + (-135771./65536.)*i6e8 + (
                        21975./8192.)*i6e6 + (-861./512.)*i6e4 + (49./128.)*i6e2,
            None,
            -3.,
            -2.,
            1.
        ),
        (
            '-O - 2*n',
            'O + 2*n',
            '2, 1, 2, 0',
            -spin_freq - 2*orbital_freq,
            (-437./141926400.)*i16 + (-173./774144.)*i14e2 + (173./3870720.)*i14 + (-463./122880.)*i12e4 + (
                        463./193536.)*i12e2 + (-463./967680.)*i12 + (-19./1280.)*i10e6 + (567./20480.)*i10e4 + (
                        -9./512.)*i10e2 + (9./2560.)*i10 + (-2065./147456.)*i8e8 + (19./288.)*i8e6 + (
                        -63./512.)*i8e4 + (5./64.)*i8e2 + (-1./64.)*i8 + (4079./614400.)*i6e10 + (2065./73728.)*i6e8 + (
                        -19./144.)*i6e6 + (63./256.)*i6e4 + (-5./32.)*i6e2 + (1./32.)*i6,
            None,
            -2.,
            -2.,
            1.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 2, 1',
            -spin_freq - orbital_freq,
            (173./15482880.)*i14e2 + (463./15482880.)*i12e4 + (-463./3870720.)*i12e2 + (39./655360.)*i10e6 + (
                        -9./40960.)*i10e4 + (9./10240.)*i10e2 + (-305./1179648.)*i8e8 + (-13./49152.)*i8e6 + (
                        1./1024.)*i8e4 + (-1./256.)*i8e2 + (1733./5242880.)*i6e10 + (305./589824.)*i6e8 + (
                        13./24576.)*i6e6 + (-1./512.)*i6e4 + (1./128.)*i6e2,
            None,
            -1.,
            -2.,
            1.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 2, 3',
            -spin_freq + orbital_freq,
            (1./655360.)*i10e6 + (-11./1179648.)*i8e8 + (-1./147456.)*i8e6 + (619./31457280.)*i6e10 + (
                        11./589824.)*i6e8 + (1./73728.)*i6e6,
            None,
            1.,
            -2.,
            1.
        ),
        (
            '-O + 2*n',
            'O - 2*n',
            '2, 1, 2, 4',
            -spin_freq + 2*orbital_freq,
            (-1./36864.)*i8e8 + (7./92160.)*i6e10 + (1./18432.)*i6e8,
            None,
            2.,
            -2.,
            1.
        ),
        (
            '-O + 3*n',
            'O - 3*n',
            '2, 1, 2, 5',
            -spin_freq + 3*orbital_freq,
            (6561./52428800.)*i6e10,
            None,
            3.,
            -2.,
            1.
        ),
        (
            '-2*O - 6*n',
            '2*O + 6*n',
            '2, 2, 0, -8',
            -2*spin_freq - 6*orbital_freq,
            (531441./40140800.)*e16,
            None,
            -6.,
            2.,
            2.
        ),
        (
            '-2*O - 5*n',
            '2*O + 5*n',
            '2, 2, 0, -7',
            -2*spin_freq - 5*orbital_freq,
            (-244140625./33294385152.)*i2e14 + (2685546875./532710162432.)*e16 + (244140625./33294385152.)*e14,
            None,
            -5.,
            2.,
            2.
        ),
        (
            '-2*O - 4*n',
            '2*O + 4*n',
            '2, 2, 0, -6',
            -2*spin_freq - 4*orbital_freq,
            (11./6075.)*i4e12 + (-8./2025.)*i2e14 + (-8./2025.)*i2e12 + (23./4725.)*e16 + (8./2025.)*e14 + (
                        8./2025.)*e12,
            None,
            -4.,
            2.,
            2.
        ),
        (
            '-2*O - 3*n',
            '2*O + 3*n',
            '2, 2, 0, -5',
            -2*spin_freq - 3*orbital_freq,
            (-16767./65536000.)*i6e10 + (24057./20971520.)*i4e12 + (24057./26214400.)*i4e10 + (
                        -4284333./1468006400.)*i2e14 + (-6561./2621440.)*i2e12 + (-6561./3276800.)*i2e10 + (
                        9428157./3355443200.)*e16 + (4284333./1468006400.)*e14 + (6561./2621440.)*e12 + (
                        6561./3276800.)*e10,
            None,
            -3.,
            2.,
            2.
        ),
        (
            '-2*O - 2*n',
            '2*O + 2*n',
            '2, 2, 0, -4',
            -2*spin_freq - 2*orbital_freq,
            (1957./92897280.)*i8e8 + (-161./1036800.)*i6e10 + (-23./207360.)*i6e8 + (10439./16588800.)*i4e12 + (
                        77./138240.)*i4e10 + (11./27648.)*i4e8 + (-2417./1814400.)*i2e14 + (-949./691200.)*i2e12 + (
                        -7./5760.)*i2e10 + (-1./1152.)*i2e8 + (22601./18579456.)*e16 + (2417./1814400.)*e14 + (
                        949./691200.)*e12 + (7./5760.)*e10 + (1./1152.)*e8,
            None,
            -2.,
            2.,
            2.
        ),
        (
            '-2*O - n',
            '2*O + n',
            '2, 2, 0, -3',
            -2*spin_freq - orbital_freq,
            (-12107./16721510400.)*i10e6 + (21527./2972712960.)*i8e8 + (1957./371589120.)*i8e6 + (
                        -14237./353894400.)*i6e10 + (-253./6635520.)*i6e8 + (-23./829440.)*i6e6 + (
                        688787./5096079360.)*i4e12 + (6809./47185920.)*i4e10 + (121./884736.)*i4e8 + (
                        11./110592.)*i4e6 + (-31398887./118908518400.)*i2e14 + (-62617./212336640.)*i2e12 + (
                        -619./1966080.)*i2e10 + (-11./36864.)*i2e8 + (-1./4608.)*i2e6 + (
                        147400583./634178764800.)*e16 + (31398887./118908518400.)*e14 + (62617./212336640.)*e12 + (
                        619./1966080.)*e10 + (11./36864.)*e8 + (1./4608.)*e6,
            None,
            -1.,
            2.,
            2.
        ),
        (
            '-2*O + n',
            '2*O - n',
            '2, 2, 0, -1',
            -2*spin_freq + orbital_freq,
            (-27269./7925299200.)*i14e2 + (-330367./30656102400.)*i12e4 + (330367./7664025600.)*i12e2 + (
                        -157391./5573836800.)*i10e6 + (12107./116121600.)*i10e4 + (-12107./29030400.)*i10e2 + (
                        119377./594542592.)*i8e8 + (25441./123863040.)*i8e6 + (-1957./2580480.)*i8e4 + (
                        1957./645120.)*i8e2 + (-39859./58982400.)*i6e10 + (-1403./1327104.)*i6e8 + (
                        -299./276480.)*i6e6 + (23./5760.)*i6e4 + (-23./1440.)*i6e2 + (3049519./1415577600.)*i4e12 + (
                        19063./7864320.)*i4e10 + (3355./884736.)*i4e8 + (143./36864.)*i4e6 + (-11./768.)*i4e4 + (
                        11./192.)*i4e2 + (-165285343./39636172800.)*i2e14 + (-277229./58982400.)*i2e12 + (
                        -1733./327680.)*i2e10 + (-305./36864.)*i2e8 + (-13./1536.)*i2e6 + (1./32.)*i2e4 + (
                        -1./8.)*i2e2 + (49450862117./13317754060800.)*e16 + (165285343./39636172800.)*e14 + (
                        277229./58982400.)*e12 + (1733./327680.)*e10 + (305./36864.)*e8 + (13./1536.)*e6 + (
                        -1./32.)*e4 + (1./8.)*e2,
            None,
            1.,
            2.,
            2.
        ),
        (
            '-2*O + 2*n',
            '2*O - 2*n',
            '2, 2, 0, 0',
            -2*spin_freq + 2*orbital_freq,
            (72518377./83691159552000.)*i16 + (27269./396264960.)*i14e2 + (-27269./1981324800.)*i14 + (
                        330367./243302400.)*i12e4 + (-330367./383201280.)*i12e2 + (330367./1916006400.)*i12 + (
                        956453./130636800.)*i10e6 + (-12107./921600.)*i10e4 + (12107./1451520.)*i10e2 + (
                        -12107./7257600.)*i10 + (7235029./371589120.)*i8e8 + (-154603./2903040.)*i8e6 + (
                        1957./20480.)*i8e4 + (-1957./32256.)*i8e2 + (1957./161280.)*i8 + (253943./6912000.)*i6e10 + (
                        -85031./829440.)*i6e8 + (1817./6480.)*i6e6 + (-161./320.)*i6e4 + (23./72.)*i6e2 + (
                        -23./360.)*i6 + (379181./13271040.)*i4e12 + (-121451./921600.)*i4e10 + (40667./110592.)*i4e8 + (
                        -869./864.)*i4e6 + (231./128.)*i4e4 + (-55./48.)*i4e2 + (11./48.)*i4 + (
                        1739939./67737600.)*i2e14 + (-34471./552960.)*i2e12 + (11041./38400.)*i2e10 + (
                        -3697./4608.)*i2e8 + (79./36.)*i2e6 + (-63./16.)*i2e4 + (5./2.)*i2e2 + (-1./2.)*i2 + (
                        -561889./240844800.)*e16 + (-1739939./67737600.)*e14 + (34471./552960.)*e12 + (
                        -11041./38400.)*e10 + (3697./4608.)*e8 + (-79./36.)*e6 + (63./16.)*e4 + (-5./2.)*e2 + (1./2.),
            None,
            2.,
            2.,
            2.
        ),
        (
            '-2*O + 3*n',
            '2*O - 3*n',
            '2, 2, 0, 1',
            -2*spin_freq + 3*orbital_freq,
            (-27269./161740800.)*i14e2 + (-13545047./1459814400.)*i12e4 + (2312569./1094860800.)*i12e2 + (
                        -3547351./24772608.)*i10e6 + (496387./5529600.)*i10e4 + (-84749./4147200.)*i10e2 + (
                        -88567949./110100480.)*i8e8 + (2867005./2752512.)*i8e6 + (-80237./122880.)*i8e4 + (
                        13699./92160.)*i8e2 + (-14230123./6553600.)*i6e10 + (1040911./245760.)*i6e8 + (
                        -33695./6144.)*i6e6 + (6601./1920.)*i6e4 + (-1127./1440.)*i6e2 + (
                        -175937223./52428800.)*i4e12 + (20417133./2621440.)*i4e10 + (-497827./32768.)*i4e8 + (
                        80575./4096.)*i4e6 + (-3157./256.)*i4e4 + (539./192.)*i4e2 + (-534226163./209715200.)*i2e14 + (
                        47982879./6553600.)*i2e12 + (-5568309./327680.)*i2e10 + (135771./4096.)*i2e8 + (
                        -21975./512.)*i2e6 + (861./32.)*i2e4 + (-49./8.)*i2e2 + (-3973253733./4697620480.)*e16 + (
                        534226163./209715200.)*e14 + (-47982879./6553600.)*e12 + (5568309./327680.)*e10 + (
                        -135771./4096.)*e8 + (21975./512.)*e6 + (-861./32.)*e4 + (49./8.)*e2,
            None,
            3.,
            2.,
            2.
        ),
        (
            '-2*O + 4*n',
            '2*O - 4*n',
            '2, 2, 0, 2',
            -2*spin_freq + 4*orbital_freq,
            (95476063./7664025600.)*i12e4 + (4733837./8709120.)*i10e6 + (-3498923./29030400.)*i10e4 + (
                        163509307./23224320.)*i8e8 + (-765187./193536.)*i8e6 + (565573./645120.)*i8e4 + (
                        3086807./86400.)*i6e10 + (-1921673./51840.)*i6e8 + (8993./432.)*i6e6 + (-6647./1440.)*i6e4 + (
                        92639041./1105920.)*i4e12 + (-1476299./11520.)*i4e10 + (919061./6912.)*i4e8 + (
                        -21505./288.)*i4e6 + (3179./192.)*i4e4 + (346700573./3628800.)*i2e14 + (
                        -8421731./46080.)*i2e12 + (134209./480.)*i2e10 + (-83551./288.)*i2e8 + (1955./12.)*i2e6 + (
                        -289./8.)*i2e4 + (173370469./4147200.)*e16 + (-346700573./3628800.)*e14 + (
                        8421731./46080.)*e12 + (-134209./480.)*e10 + (83551./288.)*e8 + (-1955./12.)*e6 + (289./8.)*e4,
            None,
            4.,
            2.,
            2.
        ),
        (
            '-2*O + 5*n',
            '2*O - 5*n',
            '2, 2, 0, 3',
            -2*spin_freq + 5*orbital_freq,
            (-345788027./668860416.)*i10e6 + (-10757090825./594542592.)*i8e8 + (279469385./74317824.)*i8e6 + (
                        -2701236725./14155776.)*i6e10 + (126424675./1327104.)*i6e8 + (-3284515./165888.)*i6e6 + (
                        -799380960125./1019215872.)*i4e12 + (6459479125./9437184.)*i4e10 + (
                        -302319875./884736.)*i4e8 + (7854275./110592.)*i4e6 + (-6429415592375./4756340736.)*i2e14 + (
                        72670996375./42467328.)*i2e12 + (-587225375./393216.)*i2e10 + (27483625./36864.)*i2e8 + (
                        -714025./4608.)*i2e6 + (-21274828753525./25367150592.)*e16 + (
                        6429415592375./4756340736.)*e14 + (-72670996375./42467328.)*e12 + (587225375./393216.)*e10 + (
                        -27483625./36864.)*e8 + (714025./4608.)*e6,
            None,
            5.,
            2.,
            2.
        ),
        (
            '-2*O + 6*n',
            '2*O - 6*n',
            '2, 2, 0, 4',
            -2*spin_freq + 6*orbital_freq,
            (555962173./41287680.)*i8e8 + (56501731./153600.)*i6e10 + (-6534047./92160.)*i6e8 + (
                        2419524481./819200.)*i4e12 + (-27022567./20480.)*i4e10 + (3124979./12288.)*i4e8 + (
                        3055721757./358400.)*i2e14 + (-659870313./102400.)*i2e12 + (7369791./2560.)*i2e10 + (
                        -284089./512.)*i2e8 + (180544031973./22937600.)*e16 + (-3055721757./358400.)*e14 + (
                        659870313./102400.)*e12 + (-7369791./2560.)*e10 + (284089./512.)*e8,
            None,
            6.,
            2.,
            2.
        ),
        (
            '-2*O + 7*n',
            '2*O - 7*n',
            '2, 2, 0, 5',
            -2*spin_freq + 7*orbital_freq,
            (-1199274105407./5308416000.)*i6e10 + (-1542795678655./339738624.)*i4e12 + (
                        573565876499./707788800.)*i4e10 + (-417304029320899./16986931200.)*i2e14 + (
                        140254152605./14155776.)*i2e12 + (-52142352409./29491200.)*i2e10 + (
                        -9997999669389767./271790899200.)*e16 + (417304029320899./16986931200.)*e14 + (
                        -140254152605./14155776.)*e12 + (52142352409./29491200.)*e10,
            None,
            7.,
            2.,
            2.
        ),
        (
            '-2*O + 8*n',
            '2*O - 8*n',
            '2, 2, 0, 6',
            -2*spin_freq + 8*orbital_freq,
            (59213111771./24883200.)*i4e12 + (56914314263./1814400.)*i2e14 + (-5383010161./1036800.)*i2e12 + (
                        578658802849./6773760.)*e16 + (-56914314263./1814400.)*e14 + (5383010161./1036800.)*e12,
            None,
            8.,
            2.,
            2.
        ),
        (
            '-2*O + 9*n',
            '2*O - 9*n',
            '2, 2, 0, 7',
            -2*spin_freq + 9*orbital_freq,
            (-147483366698529./10276044800.)*i2e14 + (-3064748517300717./32883343360.)*e16 + (
                        147483366698529./10276044800.)*e14,
            None,
            9.,
            2.,
            2.
        ),
        (
            '-2*O + 10*n',
            '2*O - 10*n',
            '2, 2, 0, 8',
            -2*spin_freq + 10*orbital_freq,
            (402063787225./10616832.)*e16,
            None,
            10.,
            2.,
            2.
        ),
        (
            '-2*O - 6*n',
            '2*O + 6*n',
            '2, 2, 1, -6',
            -2*spin_freq - 6*orbital_freq,
            (10029889./819200.)*i4e12,
            None,
            -6.,
            0.,
            2.
        ),
        (
            '-2*O - 5*n',
            '2*O + 5*n',
            '2, 2, 1, -5',
            -2*spin_freq - 5*orbital_freq,
            (-1047843./262144.)*i6e10 + (-2947317./2097152.)*i4e12 + (3143529./524288.)*i4e10,
            None,
            -5.,
            0.,
            2.
        ),
        (
            '-2*O - 4*n',
            '2*O + 4*n',
            '2, 2, 1, -4',
            -2*spin_freq - 4*orbital_freq,
            (5929./10240.)*i8e8 + (-3311./5120.)*i6e10 + (-5929./3072.)*i6e8 + (708871./204800.)*i4e12 + (
                        9933./10240.)*i4e10 + (5929./2048.)*i4e8,
            None,
            -4.,
            0.,
            2.
        ),
        (
            '-2*O - 3*n',
            '2*O + 3*n',
            '2, 2, 1, -3',
            -2*spin_freq - 3*orbital_freq,
            (-47753./967680.)*i10e6 + (20829./81920.)*i8e8 + (2809./10240.)*i8e6 + (-2006627./1310720.)*i6e10 + (
                        -6943./8192.)*i6e8 + (-2809./3072.)*i6e6 + (30355211./10485760.)*i4e12 + (
                        6019881./2621440.)*i4e10 + (20829./16384.)*i4e8 + (2809./2048.)*i4e6,
            None,
            -3.,
            0.,
            2.
        ),
        (
            '-2*O - 2*n',
            '2*O + 2*n',
            '2, 2, 1, -2',
            -2*spin_freq - 2*orbital_freq,
            (31./11200.)*i12e4 + (-17./480.)*i10e6 + (-51./2240.)*i10e4 + (1661./5120.)*i8e8 + (63./320.)*i8e6 + (
                        81./640.)*i8e4 + (-3919./2560.)*i6e10 + (-1661./1536.)*i6e8 + (-21./32.)*i6e6 + (
                        -27./64.)*i6e4 + (505601./163840.)*i4e12 + (11757./5120.)*i4e10 + (1661./1024.)*i4e8 + (
                        63./64.)*i4e6 + (81./128.)*i4e4,
            None,
            -2.,
            0.,
            2.
        ),
        (
            '-2*O - n',
            '2*O + n',
            '2, 2, 1, -1',
            -2*spin_freq - orbital_freq,
            (-1./9240.)*i14e2 + (31./11200.)*i12e4 + (31./25200.)*i12e2 + (-289./7168.)*i10e6 + (-51./2240.)*i10e4 + (
                        -17./1680.)*i10e2 + (28211./81920.)*i8e8 + (459./2048.)*i8e6 + (81./640.)*i8e4 + (
                        9./160.)*i8e2 + (-1064887./655360.)*i6e10 + (-28211./24576.)*i6e8 + (-765./1024.)*i6e6 + (
                        -27./64.)*i6e4 + (-3./16.)*i6e2 + (85553819./26214400.)*i4e12 + (3194661./1310720.)*i4e10 + (
                        28211./16384.)*i4e8 + (2295./2048.)*i4e6 + (81./128.)*i4e4 + (9./32.)*i4e2,
            None,
            -1.,
            0.,
            2.
        ),
        (
            '-2*O',
            '2*O',
            '2, 2, 1, 0',
            -2*spin_freq,
            (5461./1702701000.)*i16 + (-1./6930.)*i14e2 + (-1./20790.)*i14 + (31./9450.)*i12e4 + (31./18900.)*i12e2 + (
                        31./56700.)*i12 + (-17./378.)*i10e6 + (-17./630.)*i10e4 + (-17./1260.)*i10e2 + (
                        -17./3780.)*i10 + (3./8.)*i8e8 + (1./4.)*i8e6 + (3./20.)*i8e4 + (3./40.)*i8e2 + (1./40.)*i8 + (
                        -7./4.)*i6e10 + (-5./4.)*i6e8 + (-5./6.)*i6e6 + (-1./2.)*i6e4 + (-1./4.)*i6e2 + (-1./12.)*i6 + (
                        7./2.)*i4e12 + (21./8.)*i4e10 + (15./8.)*i4e8 + (5./4.)*i4e6 + (3./4.)*i4e4 + (3./8.)*i4e2 + (
                        1./8.)*i4,
            None,
            0.,
            0.,
            2.
        ),
        (
            '-2*O + n',
            '2*O - n',
            '2, 2, 1, 1',
            -2*spin_freq + orbital_freq,
            (-1./9240.)*i14e2 + (31./11200.)*i12e4 + (31./25200.)*i12e2 + (-289./7168.)*i10e6 + (-51./2240.)*i10e4 + (
                        -17./1680.)*i10e2 + (28211./81920.)*i8e8 + (459./2048.)*i8e6 + (81./640.)*i8e4 + (
                        9./160.)*i8e2 + (-1064887./655360.)*i6e10 + (-28211./24576.)*i6e8 + (-765./1024.)*i6e6 + (
                        -27./64.)*i6e4 + (-3./16.)*i6e2 + (85553819./26214400.)*i4e12 + (3194661./1310720.)*i4e10 + (
                        28211./16384.)*i4e8 + (2295./2048.)*i4e6 + (81./128.)*i4e4 + (9./32.)*i4e2,
            None,
            1.,
            0.,
            2.
        ),
        (
            '-2*O + 2*n',
            '2*O - 2*n',
            '2, 2, 1, 2',
            -2*spin_freq + 2*orbital_freq,
            (31./11200.)*i12e4 + (-17./480.)*i10e6 + (-51./2240.)*i10e4 + (1661./5120.)*i8e8 + (63./320.)*i8e6 + (
                        81./640.)*i8e4 + (-3919./2560.)*i6e10 + (-1661./1536.)*i6e8 + (-21./32.)*i6e6 + (
                        -27./64.)*i6e4 + (505601./163840.)*i4e12 + (11757./5120.)*i4e10 + (1661./1024.)*i4e8 + (
                        63./64.)*i4e6 + (81./128.)*i4e4,
            None,
            2.,
            0.,
            2.
        ),
        (
            '-2*O + 3*n',
            '2*O - 3*n',
            '2, 2, 1, 3',
            -2*spin_freq + 3*orbital_freq,
            (-47753./967680.)*i10e6 + (20829./81920.)*i8e8 + (2809./10240.)*i8e6 + (-2006627./1310720.)*i6e10 + (
                        -6943./8192.)*i6e8 + (-2809./3072.)*i6e6 + (30355211./10485760.)*i4e12 + (
                        6019881./2621440.)*i4e10 + (20829./16384.)*i4e8 + (2809./2048.)*i4e6,
            None,
            3.,
            0.,
            2.
        ),
        (
            '-2*O + 4*n',
            '2*O - 4*n',
            '2, 2, 1, 4',
            -2*spin_freq + 4*orbital_freq,
            (5929./10240.)*i8e8 + (-3311./5120.)*i6e10 + (-5929./3072.)*i6e8 + (708871./204800.)*i4e12 + (
                        9933./10240.)*i4e10 + (5929./2048.)*i4e8,
            None,
            4.,
            0.,
            2.
        ),
        (
            '-2*O + 5*n',
            '2*O - 5*n',
            '2, 2, 1, 5',
            -2*spin_freq + 5*orbital_freq,
            (-1047843./262144.)*i6e10 + (-2947317./2097152.)*i4e12 + (3143529./524288.)*i4e10,
            None,
            5.,
            0.,
            2.
        ),
        (
            '-2*O + 6*n',
            '2*O - 6*n',
            '2, 2, 1, 6',
            -2*spin_freq + 6*orbital_freq,
            (10029889./819200.)*i4e12,
            None,
            6.,
            0.,
            2.
        ),
        (
            '-2*O - 6*n',
            '2*O + 6*n',
            '2, 2, 2, -4',
            -2*spin_freq - 6*orbital_freq,
            (284089./131072.)*i8e8,
            None,
            -6.,
            -2.,
            2.
        ),
        (
            '-2*O - 5*n',
            '2*O + 5*n',
            '2, 2, 2, -3',
            -2*spin_freq - 5*orbital_freq,
            (-714025./3538944.)*i10e6 + (-27483625./9437184.)*i8e8 + (714025./1179648.)*i8e6,
            None,
            -5.,
            -2.,
            2.
        ),
        (
            '-2*O - 4*n',
            '2*O + 4*n',
            '2, 2, 2, -2',
            -2*spin_freq - 4*orbital_freq,
            (5491./737280.)*i12e4 + (1955./9216.)*i10e6 + (-289./6144.)*i10e4 + (83551./73728.)*i8e8 + (
                        -1955./3072.)*i8e6 + (289./2048.)*i8e4,
            None,
            -4.,
            -2.,
            2.
        ),
        (
            '-2*O - 3*n',
            '2*O + 3*n',
            '2, 2, 2, -1',
            -2*spin_freq - 3*orbital_freq,
            (-7./55296.)*i14e2 + (-5453./983040.)*i12e4 + (931./737280.)*i12e2 + (-7325./131072.)*i10e6 + (
                        287./8192.)*i10e4 + (-49./6144.)*i10e2 + (-135771./1048576.)*i8e8 + (21975./131072.)*i8e6 + (
                        -861./8192.)*i8e4 + (49./2048.)*i8e2,
            None,
            -3.,
            -2.,
            2.
        ),
        (
            '-2*O - 2*n',
            '2*O + 2*n',
            '2, 2, 2, 0',
            -2*spin_freq - 2*orbital_freq,
            (457./619315200.)*i16 + (5./96768.)*i14e2 + (-1./96768.)*i14 + (133./163840.)*i12e4 + (
                        -19./36864.)*i12e2 + (19./184320.)*i12 + (19./6912.)*i10e6 + (-21./4096.)*i10e4 + (
                        5./1536.)*i10e2 + (-1./1536.)*i10 + (2065./1179648.)*i8e8 + (-19./2304.)*i8e6 + (
                        63./4096.)*i8e4 + (-5./512.)*i8e2 + (1./512.)*i8,
            None,
            -2.,
            -2.,
            2.
        ),
        (
            '-2*O - n',
            '2*O + n',
            '2, 2, 2, 1',
            -2*spin_freq - orbital_freq,
            (-1./387072.)*i14e2 + (-19./2949120.)*i12e4 + (19./737280.)*i12e2 + (-13./1179648.)*i10e6 + (
                        1./24576.)*i10e4 + (-1./6144.)*i10e2 + (305./9437184.)*i8e8 + (13./393216.)*i8e6 + (
                        -1./8192.)*i8e4 + (1./2048.)*i8e2,
            None,
            -1.,
            -2.,
            2.
        ),
        (
            '-2*O + n',
            '2*O - n',
            '2, 2, 2, 3',
            -2*spin_freq + orbital_freq,
            (-1./3538944.)*i10e6 + (11./9437184.)*i8e8 + (1./1179648.)*i8e6,
            None,
            1.,
            -2.,
            2.
        ),
        (
            '-2*O + 2*n',
            '2*O - 2*n',
            '2, 2, 2, 4',
            -2*spin_freq + 2*orbital_freq,
            (1./294912.)*i8e8,
            None,
            2.,
            -2.,
            2.
        )
    )

    return mode_data_output

@njit
def nsr_modes_t18(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray) \
        -> OutputType:
    """ Non-sync Tidal Modes for truncation level 18, tidal order 2

    Parameters
    ----------
    orbital_freq : np.ndarray
        Planet's orbital mean motion [rad s-1]
    spin_freq : np.ndarray
        Planet's spin frequency [rad s-1]
    eccentricity : np.ndarray
        Planet's eccentricity
    inclination : np.ndarray
        Planet's inclination [rads]

    Returns
    -------
    mode_data_tuple : OutputType
        List of mode datas. Each mode data contains information to calculate the heating and torques.
    """
    # Tidal Potential Builder for l-order = 2, NSR = True.
    # Max Eccentricity Order = 18
    # Max Inclination Order = 18
    # Max q = 10.
    # Number of unique modes = 56.
    # Number of unique frequencies = 47.

    i2 = inclination**2
    i4 = inclination**4
    i6 = inclination**6
    i8 = inclination**8
    i10 = inclination**10
    i12 = inclination**12
    i14 = inclination**14
    i16 = inclination**16
    i18 = inclination**18
    e2 = eccentricity**2
    e4 = eccentricity**4
    e6 = eccentricity**6
    e8 = eccentricity**8
    e10 = eccentricity**10
    e12 = eccentricity**12
    e14 = eccentricity**14
    e16 = eccentricity**16
    e18 = eccentricity**18
    i2e2 = i2*e2
    i2e4 = i2*e4
    i2e6 = i2*e6
    i2e8 = i2*e8
    i2e10 = i2*e10
    i2e12 = i2*e12
    i2e14 = i2*e14
    i2e16 = i2*e16
    i4e2 = i4*e2
    i4e4 = i4*e4
    i4e6 = i4*e6
    i4e8 = i4*e8
    i4e10 = i4*e10
    i4e12 = i4*e12
    i4e14 = i4*e14
    i6e2 = i6*e2
    i6e4 = i6*e4
    i6e6 = i6*e6
    i6e8 = i6*e8
    i6e10 = i6*e10
    i6e12 = i6*e12
    i8e2 = i8*e2
    i8e4 = i8*e4
    i8e6 = i8*e6
    i8e8 = i8*e8
    i8e10 = i8*e10
    i10e2 = i10*e2
    i10e4 = i10*e4
    i10e6 = i10*e6
    i10e8 = i10*e8
    i12e2 = i12*e2
    i12e4 = i12*e4
    i12e6 = i12*e6
    i14e2 = i14*e2
    i14e4 = i14*e4
    i16e2 = i16*e2

    mode_data_output = (
        (
            '-5*n',
            '5*n',
            '2, 0, 0, -7',
            -5*orbital_freq,
            (244140625./177570054144.)*i4e14,
            None,
            -5.,
            2.,
            0.
        ),
        (
            '-4*n',
            '4*n',
            '2, 0, 0, -6',
            -4*orbital_freq,
            (-1./2025.)*i6e12 + (1./1350.)*i4e14 + (1./1350.)*i4e12,
            None,
            -4.,
            2.,
            0.
        ),
        (
            '-3*n',
            '3*n',
            '2, 0, 0, -5',
            -3*orbital_freq,
            (19683./262144000.)*i8e10 + (-6561./20971520.)*i6e12 + (-6561./26214400.)*i6e10 + (
                        12852999./23488102400.)*i4e14 + (19683./41943040.)*i4e12 + (19683./52428800.)*i4e10,
            None,
            -3.,
            2.,
            0.
        ),
        (
            '-2*n',
            '2*n',
            '2, 0, 0, -4',
            -2*orbital_freq,
            (-17./2903040.)*i10e8 + (7./153600.)*i8e10 + (1./30720.)*i8e8 + (-949./5529600.)*i6e12 + (
                        -7./46080.)*i6e10 + (-1./9216.)*i6e8 + (2417./9676800.)*i4e14 + (949./3686400.)*i4e12 + (
                        7./30720.)*i4e10 + (1./6144.)*i4e8,
            None,
            -2.,
            2.,
            0.
        ),
        (
            '-n',
            'n',
            '2, 0, 0, -3',
            -orbital_freq,
            (31./174182400.)*i12e6 + (-187./92897280.)*i10e8 + (-17./11612160.)*i10e6 + (619./52428800.)*i8e10 + (
                        11./983040.)*i8e8 + (1./122880.)*i8e6 + (-62617./1698693120.)*i6e12 + (
                        -619./15728640.)*i6e10 + (-11./294912.)*i6e8 + (-1./36864.)*i6e6 + (
                        31398887./634178764800.)*i4e14 + (62617./1132462080.)*i4e12 + (619./10485760.)*i4e10 + (
                        11./196608.)*i4e8 + (1./24576.)*i4e6,
            None,
            -1.,
            2.,
            0.
        ),
        (
            'n',
            'n',
            '2, 0, 0, -1',
            orbital_freq,
            (5461./9081072000.)*i16e2 + (1./443520.)*i14e4 + (-1./110880.)*i14e2 + (403./58060800.)*i12e6 + (
                        -31./1209600.)*i12e4 + (31./302400.)*i12e2 + (-1037./18579456.)*i10e8 + (
                        -221./3870720.)*i10e6 + (17./80640.)*i10e4 + (-17./20160.)*i10e2 + (5199./26214400.)*i8e10 + (
                        61./196608.)*i8e8 + (13./40960.)*i8e6 + (-3./2560.)*i8e4 + (3./640.)*i8e2 + (
                        -277229./471859200.)*i6e12 + (-1733./2621440.)*i6e10 + (-305./294912.)*i6e8 + (
                        -13./12288.)*i6e6 + (1./256.)*i6e4 + (-1./64.)*i6e2 + (165285343./211392921600.)*i4e14 + (
                        277229./314572800.)*i4e12 + (5199./5242880.)*i4e10 + (305./196608.)*i4e8 + (13./8192.)*i4e6 + (
                        -3./512.)*i4e4 + (3./128.)*i4e2,
            None,
            1.,
            2.,
            0.
        ),
        (
            '2*n',
            '2*n',
            '2, 0, 0, 0',
            2*orbital_freq,
            (-257./2043241200.)*i18 + (-5461./454053600.)*i16e2 + (5461./2270268000.)*i16 + (-1./3520.)*i14e4 + (
                        1./5544.)*i14e2 + (-1./27720.)*i14 + (-2449./1360800.)*i12e6 + (31./9600.)*i12e4 + (
                        -31./15120.)*i12e2 + (31./75600.)*i12 + (-62849./11612160.)*i10e8 + (1343./90720.)*i10e6 + (
                        -17./640.)*i10e4 + (17./1008.)*i10e2 + (-17./5040.)*i10 + (-11041./1024000.)*i8e10 + (
                        3697./122880.)*i8e8 + (-79./960.)*i8e6 + (189./1280.)*i8e4 + (-3./32.)*i8e2 + (3./160.)*i8 + (
                        -34471./4423680.)*i6e12 + (11041./307200.)*i6e10 + (-3697./36864.)*i6e8 + (79./288.)*i6e6 + (
                        -63./128.)*i6e4 + (5./16.)*i6e2 + (-1./16.)*i6 + (-1739939./361267200.)*i4e14 + (
                        34471./2949120.)*i4e12 + (-11041./204800.)*i4e10 + (3697./24576.)*i4e8 + (-79./192.)*i4e6 + (
                        189./256.)*i4e4 + (-15./32.)*i4e2 + (3./32.)*i4,
            None,
            2.,
            2.,
            0.
        ),
        (
            '3*n',
            '3*n',
            '2, 0, 0, 1',
            3*orbital_freq,
            (5461./185328000.)*i16e2 + (41./21120.)*i14e4 + (-7./15840.)*i14e2 + (9083./258048.)*i12e6 + (
                        -1271./57600.)*i12e4 + (217./43200.)*i12e2 + (769369./3440640.)*i10e8 + (
                        -24905./86016.)*i10e6 + (697./3840.)*i10e4 + (-119./2880.)*i10e2 + (
                        16704927./26214400.)*i8e10 + (-407313./327680.)*i8e8 + (13185./8192.)*i8e6 + (
                        -2583./2560.)*i8e4 + (147./640.)*i8e2 + (47982879./52428800.)*i6e12 + (
                        -5568309./2621440.)*i6e10 + (135771./32768.)*i6e8 + (-21975./4096.)*i6e6 + (861./256.)*i6e4 + (
                        -49./64.)*i6e2 + (1602678489./3355443200.)*i4e14 + (-143948637./104857600.)*i4e12 + (
                        16704927./5242880.)*i4e10 + (-407313./65536.)*i4e8 + (65925./8192.)*i4e6 + (
                        -2583./512.)*i4e4 + (147./128.)*i4e2,
            None,
            3.,
            2.,
            0.
        ),
        (
            '4*n',
            '4*n',
            '2, 0, 0, 2',
            4*orbital_freq,
            (-289./110880.)*i14e4 + (-12121./90720.)*i12e6 + (8959./302400.)*i12e4 + (-1420367./725760.)*i10e8 + (
                        6647./6048.)*i10e6 + (-4913./20160.)*i10e4 + (-134209./12800.)*i8e10 + (83551./7680.)*i8e8 + (
                        -391./64.)*i8e6 + (867./640.)*i8e4 + (-8421731./368640.)*i6e12 + (134209./3840.)*i6e10 + (
                        -83551./2304.)*i6e8 + (1955./96.)*i6e6 + (-289./64.)*i6e4 + (-346700573./19353600.)*i4e14 + (
                        8421731./245760.)*i4e12 + (-134209./2560.)*i4e10 + (83551./1536.)*i4e8 + (-1955./64.)*i4e6 + (
                        867./128.)*i4e4,
            None,
            4.,
            2.,
            0.
        ),
        (
            '5*n',
            '5*n',
            '2, 0, 0, 3',
            5*orbital_freq,
            (885391./6967296.)*i12e6 + (93444325./18579456.)*i10e8 + (-2427685./2322432.)*i10e6 + (
                        117445075./2097152.)*i8e10 + (-5496725./196608.)*i8e8 + (142805./24576.)*i8e6 + (
                        72670996375./339738624.)*i6e12 + (-587225375./3145728.)*i6e10 + (27483625./294912.)*i6e8 + (
                        -714025./36864.)*i6e6 + (6429415592375./25367150592.)*i4e14 + (
                        -72670996375./226492416.)*i4e12 + (587225375./2097152.)*i4e10 + (-27483625./196608.)*i4e8 + (
                        714025./24576.)*i4e6,
            None,
            5.,
            2.,
            0.
        ),
        (
            '6*n',
            '6*n',
            '2, 0, 0, 4',
            6*orbital_freq,
            (-4829513./1290240.)*i10e8 + (-22109373./204800.)*i8e10 + (852267./40960.)*i8e8 + (
                        -659870313./819200.)*i6e12 + (7369791./20480.)*i6e10 + (-284089./4096.)*i6e8 + (
                        -9167165271./5734400.)*i4e14 + (1979610939./1638400.)*i4e12 + (-22109373./40960.)*i4e10 + (
                        852267./8192.)*i4e8,
            None,
            6.,
            2.,
            0.
        ),
        (
            '7*n',
            '7*n',
            '2, 0, 0, 5',
            7*orbital_freq,
            (52142352409./786432000.)*i8e10 + (140254152605./113246208.)*i6e12 + (-52142352409./235929600.)*i6e10 + (
                        417304029320899./90596966400.)*i4e14 + (-140254152605./75497472.)*i4e12 + (
                        52142352409./157286400.)*i4e10,
            None,
            7.,
            2.,
            0.
        ),
        (
            '8*n',
            '8*n',
            '2, 0, 0, 6',
            8*orbital_freq,
            (-5383010161./8294400.)*i6e12 + (-56914314263./9676800.)*i4e14 + (5383010161./5529600.)*i4e12,
            None,
            8.,
            2.,
            0.
        ),
        (
            '9*n',
            '9*n',
            '2, 0, 0, 7',
            9*orbital_freq,
            (442450100095587./164416716800.)*i4e14,
            None,
            9.,
            2.,
            0.
        ),
        (
            '-9*n',
            '9*n',
            '2, 0, 1, -9',
            -9*orbital_freq,
            (4143587919679849./31568009625600.)*e18,
            None,
            -9.,
            0.,
            0.
        ),
        (
            '-8*n',
            '8*n',
            '2, 0, 1, -8',
            -8*orbital_freq,
            (-31801945561./160563200.)*i2e16 + (-1606626956771./13005619200.)*e18 + (31801945561./481689600.)*e16,
            None,
            -8.,
            0.,
            0.
        ),
        (
            '-7*n',
            '7*n',
            '2, 0, 1, -7',
            -7*orbital_freq,
            (2427134219653./22649241600.)*i4e14 + (1328179895713./10066329600.)*i2e16 + (
                        -186702632281./1887436800.)*i2e14 + (999990833612299./17394617548800.)*e18 + (
                        -1328179895713./30198988800.)*e16 + (186702632281./5662310400.)*e14,
            None,
            -7.,
            0.,
            0.
        ),
        (
            '-6*n',
            '6*n',
            '2, 0, 1, -6',
            -6*orbital_freq,
            (-491464561./18432000.)*i6e12 + (-30013659./716800.)*i4e14 + (130388557./2457600.)*i4e12 + (
                        -197537355./3211264.)*i2e16 + (6926229./179200.)*i2e14 + (-10029889./204800.)*i2e12 + (
                        -941362189./240844800.)*e18 + (65845785./3211264.)*e16 + (-2308743./179200.)*e14 + (
                        10029889./614400.)*e12,
            None,
            -6.,
            0.,
            0.
        ),
        (
            '-5*n',
            '5*n',
            '2, 0, 1, -5',
            -5*orbital_freq,
            (67411233./18350080.)*i8e10 + (48139511./15728640.)*i6e12 + (-17114769./1310720.)*i6e10 + (
                        175814549209./6341787648.)*i4e14 + (-12771707./2097152.)*i4e12 + (13621959./524288.)*i4e10 + (
                        -28491858875./2818572288.)*i2e16 + (-13524196093./528482304.)*i2e14 + (
                        2947317./524288.)*i2e12 + (-3143529./131072.)*i2e10 + (2447539096445./315680096256.)*e18 + (
                        28491858875./8455716864.)*e16 + (13524196093./1585446912.)*e14 + (-982439./524288.)*e12 + (
                        1047843./131072.)*e10,
            None,
            -5.,
            0.,
            0.
        ),
        (
            '-4*n',
            '4*n',
            '2, 0, 1, -4',
            -4*orbital_freq,
            (-651343./2073600.)*i10e8 + (91289./153600.)*i8e10 + (163471./92160.)*i8e8 + (-34734679./4608000.)*i6e12 + (
                        -162239./76800.)*i6e10 + (-290521./46080.)*i6e8 + (39355771./2764800.)*i4e14 + (
                        9215323./614400.)*i4e12 + (43043./10240.)*i4e10 + (77077./6144.)*i4e8 + (
                        -124914751./6881280.)*i2e16 + (-3027367./230400.)*i2e14 + (-708871./51200.)*i2e12 + (
                        -9933./2560.)*i2e10 + (-5929./512.)*i2e8 + (2170376447./309657600.)*e18 + (
                        124914751./20643840.)*e16 + (3027367./691200.)*e14 + (708871./153600.)*e12 + (
                        3311./2560.)*e10 + (5929./1536.)*e8,
            None,
            -4.,
            0.,
            0.
        ),
        (
            '-3*n',
            '3*n',
            '2, 0, 1, -3',
            -3*orbital_freq,
            (1233151./68428800.)*i12e6 + (-5339167./38707200.)*i10e8 + (-2160121./14515200.)*i10e6 + (
                        55325573./39321600.)*i8e10 + (1339999./1720320.)*i8e8 + (542137./645120.)*i8e6 + (
                        -1487405339./235929600.)*i6e12 + (-98324723./19660800.)*i6e10 + (-340207./122880.)*i6e8 + (
                        -137641./46080.)*i6e6 + (13716852279./838860800.)*i4e14 + (394617743./31457280.)*i4e12 + (
                        26086151./2621440.)*i4e10 + (90259./16384.)*i4e8 + (36517./6144.)*i4e6 + (
                        -437242738059./23488102400.)*i2e16 + (-3165427449./209715200.)*i2e14 + (
                        -30355211./2621440.)*i2e12 + (-6019881./655360.)*i2e10 + (-20829./4096.)*i2e8 + (
                        -2809./512.)*i2e6 + (2818313897141./375809638400.)*e18 + (145747579353./23488102400.)*e16 + (
                        1055142483./209715200.)*e14 + (30355211./7864320.)*e12 + (2006627./655360.)*e10 + (
                        6943./4096.)*e8 + (2809./1536.)*e6,
            None,
            -3.,
            0.,
            0.
        ),
        (
            '-2*n',
            '2*n',
            '2, 0, 1, -2',
            -2*orbital_freq,
            (-12289./16816800.)*i14e4 + (3073./237600.)*i12e6 + (439./52800.)*i12e4 + (-1277309./7257600.)*i10e8 + (
                        -769./7200.)*i10e6 + (-769./11200.)*i10e4 + (756367./537600.)*i8e10 + (320573./322560.)*i8e8 + (
                        193./320.)*i8e6 + (1737./4480.)*i8e4 + (-24774449./3686400.)*i6e12 + (-192031./38400.)*i6e10 + (
                        -81389./23040.)*i6e8 + (-343./160.)*i6e6 + (-441./320.)*i6e4 + (332356609./19353600.)*i4e14 + (
                        6572813./491520.)*i4e12 + (50947./5120.)*i4e10 + (21593./3072.)*i4e8 + (273./64.)*i4e6 + (
                        351./128.)*i4e4 + (-678544541./34406400.)*i2e16 + (-25565893./1612800.)*i2e14 + (
                        -505601./40960.)*i2e12 + (-11757./1280.)*i2e10 + (-1661./256.)*i2e8 + (-63./16.)*i2e6 + (
                        -81./32.)*i2e4 + (51883919761./6502809600.)*e18 + (678544541./103219200.)*e16 + (
                        25565893./4838400.)*e14 + (505601./122880.)*e12 + (3919./1280.)*e10 + (1661./768.)*e8 + (
                        21./16.)*e6 + (27./32.)*e4,
            None,
            -2.,
            0.,
            0.
        ),
        (
            '-n',
            'n',
            '2, 0, 1, -1',
            -orbital_freq,
            (3781./174636000.)*i16e2 + (-12289./16816800.)*i14e4 + (-12289./37837800.)*i14e2 + (7463./506880.)*i12e6 + (
                        439./52800.)*i12e4 + (439./118800.)*i12e2 + (-21694259./116121600.)*i10e8 + (
                        -13073./107520.)*i10e6 + (-769./11200.)*i10e4 + (-769./25200.)*i10e2 + (
                        205523191./137625600.)*i8e10 + (5444723./5160960.)*i8e8 + (9843./14336.)*i8e6 + (
                        1737./4480.)*i8e4 + (193./1120.)*i8e2 + (-4192137131./589824000.)*i6e12 + (
                        -52179463./9830400.)*i6e10 + (-1382339./368640.)*i6e8 + (-2499./1024.)*i6e6 + (
                        -441./320.)*i6e4 + (-49./80.)*i6e2 + (320510498633./17616076800.)*i4e14 + (
                        1112199647./78643200.)*i4e12 + (13843531./1310720.)*i4e10 + (366743./49152.)*i4e8 + (
                        9945./2048.)*i4e6 + (351./128.)*i4e4 + (39./32.)*i4e2 + (-689312857627./32883343360.)*i2e16 + (
                        -24654653741./1468006400.)*i2e14 + (-85553819./6553600.)*i2e12 + (-3194661./327680.)*i2e10 + (
                        -28211./4096.)*i2e8 + (-2295./512.)*i2e6 + (-81./32.)*i2e4 + (-9./8.)*i2e2 + (
                        725941889609009./85233625989120.)*e18 + (689312857627./98650030080.)*e16 + (
                        24654653741./4404019200.)*e14 + (85553819./19660800.)*e12 + (1064887./327680.)*e10 + (
                        28211./12288.)*e8 + (765./512.)*e6 + (27./32.)*e4 + (3./8.)*e2,
            None,
            -1.,
            0.,
            0.
        ),
        (
            'n',
            'n',
            '2, 0, 1, 1',
            orbital_freq,
            (3781./174636000.)*i16e2 + (-12289./16816800.)*i14e4 + (-12289./37837800.)*i14e2 + (7463./506880.)*i12e6 + (
                        439./52800.)*i12e4 + (439./118800.)*i12e2 + (-21694259./116121600.)*i10e8 + (
                        -13073./107520.)*i10e6 + (-769./11200.)*i10e4 + (-769./25200.)*i10e2 + (
                        205523191./137625600.)*i8e10 + (5444723./5160960.)*i8e8 + (9843./14336.)*i8e6 + (
                        1737./4480.)*i8e4 + (193./1120.)*i8e2 + (-4192137131./589824000.)*i6e12 + (
                        -52179463./9830400.)*i6e10 + (-1382339./368640.)*i6e8 + (-2499./1024.)*i6e6 + (
                        -441./320.)*i6e4 + (-49./80.)*i6e2 + (320510498633./17616076800.)*i4e14 + (
                        1112199647./78643200.)*i4e12 + (13843531./1310720.)*i4e10 + (366743./49152.)*i4e8 + (
                        9945./2048.)*i4e6 + (351./128.)*i4e4 + (39./32.)*i4e2 + (-689312857627./32883343360.)*i2e16 + (
                        -24654653741./1468006400.)*i2e14 + (-85553819./6553600.)*i2e12 + (-3194661./327680.)*i2e10 + (
                        -28211./4096.)*i2e8 + (-2295./512.)*i2e6 + (-81./32.)*i2e4 + (-9./8.)*i2e2 + (
                        725941889609009./85233625989120.)*e18 + (689312857627./98650030080.)*e16 + (
                        24654653741./4404019200.)*e14 + (85553819./19660800.)*e12 + (1064887./327680.)*e10 + (
                        28211./12288.)*e8 + (765./512.)*e6 + (27./32.)*e4 + (3./8.)*e2,
            None,
            1.,
            0.,
            0.
        ),
        (
            '2*n',
            '2*n',
            '2, 0, 1, 2',
            2*orbital_freq,
            (-12289./16816800.)*i14e4 + (3073./237600.)*i12e6 + (439./52800.)*i12e4 + (-1277309./7257600.)*i10e8 + (
                        -769./7200.)*i10e6 + (-769./11200.)*i10e4 + (756367./537600.)*i8e10 + (320573./322560.)*i8e8 + (
                        193./320.)*i8e6 + (1737./4480.)*i8e4 + (-24774449./3686400.)*i6e12 + (-192031./38400.)*i6e10 + (
                        -81389./23040.)*i6e8 + (-343./160.)*i6e6 + (-441./320.)*i6e4 + (332356609./19353600.)*i4e14 + (
                        6572813./491520.)*i4e12 + (50947./5120.)*i4e10 + (21593./3072.)*i4e8 + (273./64.)*i4e6 + (
                        351./128.)*i4e4 + (-678544541./34406400.)*i2e16 + (-25565893./1612800.)*i2e14 + (
                        -505601./40960.)*i2e12 + (-11757./1280.)*i2e10 + (-1661./256.)*i2e8 + (-63./16.)*i2e6 + (
                        -81./32.)*i2e4 + (51883919761./6502809600.)*e18 + (678544541./103219200.)*e16 + (
                        25565893./4838400.)*e14 + (505601./122880.)*e12 + (3919./1280.)*e10 + (1661./768.)*e8 + (
                        21./16.)*e6 + (27./32.)*e4,
            None,
            2.,
            0.,
            0.
        ),
        (
            '3*n',
            '3*n',
            '2, 0, 1, 3',
            3*orbital_freq,
            (1233151./68428800.)*i12e6 + (-5339167./38707200.)*i10e8 + (-2160121./14515200.)*i10e6 + (
                        55325573./39321600.)*i8e10 + (1339999./1720320.)*i8e8 + (542137./645120.)*i8e6 + (
                        -1487405339./235929600.)*i6e12 + (-98324723./19660800.)*i6e10 + (-340207./122880.)*i6e8 + (
                        -137641./46080.)*i6e6 + (13716852279./838860800.)*i4e14 + (394617743./31457280.)*i4e12 + (
                        26086151./2621440.)*i4e10 + (90259./16384.)*i4e8 + (36517./6144.)*i4e6 + (
                        -437242738059./23488102400.)*i2e16 + (-3165427449./209715200.)*i2e14 + (
                        -30355211./2621440.)*i2e12 + (-6019881./655360.)*i2e10 + (-20829./4096.)*i2e8 + (
                        -2809./512.)*i2e6 + (2818313897141./375809638400.)*e18 + (145747579353./23488102400.)*e16 + (
                        1055142483./209715200.)*e14 + (30355211./7864320.)*e12 + (2006627./655360.)*e10 + (
                        6943./4096.)*e8 + (2809./1536.)*e6,
            None,
            3.,
            0.,
            0.
        ),
        (
            '4*n',
            '4*n',
            '2, 0, 1, 4',
            4*orbital_freq,
            (-651343./2073600.)*i10e8 + (91289./153600.)*i8e10 + (163471./92160.)*i8e8 + (-34734679./4608000.)*i6e12 + (
                        -162239./76800.)*i6e10 + (-290521./46080.)*i6e8 + (39355771./2764800.)*i4e14 + (
                        9215323./614400.)*i4e12 + (43043./10240.)*i4e10 + (77077./6144.)*i4e8 + (
                        -124914751./6881280.)*i2e16 + (-3027367./230400.)*i2e14 + (-708871./51200.)*i2e12 + (
                        -9933./2560.)*i2e10 + (-5929./512.)*i2e8 + (2170376447./309657600.)*e18 + (
                        124914751./20643840.)*e16 + (3027367./691200.)*e14 + (708871./153600.)*e12 + (
                        3311./2560.)*e10 + (5929./1536.)*e8,
            None,
            4.,
            0.,
            0.
        ),
        (
            '5*n',
            '5*n',
            '2, 0, 1, 5',
            5*orbital_freq,
            (67411233./18350080.)*i8e10 + (48139511./15728640.)*i6e12 + (-17114769./1310720.)*i6e10 + (
                        175814549209./6341787648.)*i4e14 + (-12771707./2097152.)*i4e12 + (13621959./524288.)*i4e10 + (
                        -28491858875./2818572288.)*i2e16 + (-13524196093./528482304.)*i2e14 + (
                        2947317./524288.)*i2e12 + (-3143529./131072.)*i2e10 + (2447539096445./315680096256.)*e18 + (
                        28491858875./8455716864.)*e16 + (13524196093./1585446912.)*e14 + (-982439./524288.)*e12 + (
                        1047843./131072.)*e10,
            None,
            5.,
            0.,
            0.
        ),
        (
            '6*n',
            '6*n',
            '2, 0, 1, 6',
            6*orbital_freq,
            (-491464561./18432000.)*i6e12 + (-30013659./716800.)*i4e14 + (130388557./2457600.)*i4e12 + (
                        -197537355./3211264.)*i2e16 + (6926229./179200.)*i2e14 + (-10029889./204800.)*i2e12 + (
                        -941362189./240844800.)*e18 + (65845785./3211264.)*e16 + (-2308743./179200.)*e14 + (
                        10029889./614400.)*e12,
            None,
            6.,
            0.,
            0.
        ),
        (
            '7*n',
            '7*n',
            '2, 0, 1, 7',
            7*orbital_freq,
            (2427134219653./22649241600.)*i4e14 + (1328179895713./10066329600.)*i2e16 + (
                        -186702632281./1887436800.)*i2e14 + (999990833612299./17394617548800.)*e18 + (
                        -1328179895713./30198988800.)*e16 + (186702632281./5662310400.)*e14,
            None,
            7.,
            0.,
            0.
        ),
        (
            '8*n',
            '8*n',
            '2, 0, 1, 8',
            8*orbital_freq,
            (-31801945561./160563200.)*i2e16 + (-1606626956771./13005619200.)*e18 + (31801945561./481689600.)*e16,
            None,
            8.,
            0.,
            0.
        ),
        (
            '9*n',
            '9*n',
            '2, 0, 1, 9',
            9*orbital_freq,
            (4143587919679849./31568009625600.)*e18,
            None,
            9.,
            0.,
            0.
        ),
        (
            '-9*n',
            '9*n',
            '2, 0, 2, -7',
            -9*orbital_freq,
            (442450100095587./164416716800.)*i4e14,
            None,
            -9.,
            -2.,
            0.
        ),
        (
            '-8*n',
            '8*n',
            '2, 0, 2, -6',
            -8*orbital_freq,
            (-5383010161./8294400.)*i6e12 + (-56914314263./9676800.)*i4e14 + (5383010161./5529600.)*i4e12,
            None,
            -8.,
            -2.,
            0.
        ),
        (
            '-7*n',
            '7*n',
            '2, 0, 2, -5',
            -7*orbital_freq,
            (52142352409./786432000.)*i8e10 + (140254152605./113246208.)*i6e12 + (-52142352409./235929600.)*i6e10 + (
                        417304029320899./90596966400.)*i4e14 + (-140254152605./75497472.)*i4e12 + (
                        52142352409./157286400.)*i4e10,
            None,
            -7.,
            -2.,
            0.
        ),
        (
            '-6*n',
            '6*n',
            '2, 0, 2, -4',
            -6*orbital_freq,
            (-4829513./1290240.)*i10e8 + (-22109373./204800.)*i8e10 + (852267./40960.)*i8e8 + (
                        -659870313./819200.)*i6e12 + (7369791./20480.)*i6e10 + (-284089./4096.)*i6e8 + (
                        -9167165271./5734400.)*i4e14 + (1979610939./1638400.)*i4e12 + (-22109373./40960.)*i4e10 + (
                        852267./8192.)*i4e8,
            None,
            -6.,
            -2.,
            0.
        ),
        (
            '-5*n',
            '5*n',
            '2, 0, 2, -3',
            -5*orbital_freq,
            (885391./6967296.)*i12e6 + (93444325./18579456.)*i10e8 + (-2427685./2322432.)*i10e6 + (
                        117445075./2097152.)*i8e10 + (-5496725./196608.)*i8e8 + (142805./24576.)*i8e6 + (
                        72670996375./339738624.)*i6e12 + (-587225375./3145728.)*i6e10 + (27483625./294912.)*i6e8 + (
                        -714025./36864.)*i6e6 + (6429415592375./25367150592.)*i4e14 + (
                        -72670996375./226492416.)*i4e12 + (587225375./2097152.)*i4e10 + (-27483625./196608.)*i4e8 + (
                        714025./24576.)*i4e6,
            None,
            -5.,
            -2.,
            0.
        ),
        (
            '-4*n',
            '4*n',
            '2, 0, 2, -2',
            -4*orbital_freq,
            (-289./110880.)*i14e4 + (-12121./90720.)*i12e6 + (8959./302400.)*i12e4 + (-1420367./725760.)*i10e8 + (
                        6647./6048.)*i10e6 + (-4913./20160.)*i10e4 + (-134209./12800.)*i8e10 + (83551./7680.)*i8e8 + (
                        -391./64.)*i8e6 + (867./640.)*i8e4 + (-8421731./368640.)*i6e12 + (134209./3840.)*i6e10 + (
                        -83551./2304.)*i6e8 + (1955./96.)*i6e6 + (-289./64.)*i6e4 + (-346700573./19353600.)*i4e14 + (
                        8421731./245760.)*i4e12 + (-134209./2560.)*i4e10 + (83551./1536.)*i4e8 + (-1955./64.)*i4e6 + (
                        867./128.)*i4e4,
            None,
            -4.,
            -2.,
            0.
        ),
        (
            '-3*n',
            '3*n',
            '2, 0, 2, -1',
            -3*orbital_freq,
            (5461./185328000.)*i16e2 + (41./21120.)*i14e4 + (-7./15840.)*i14e2 + (9083./258048.)*i12e6 + (
                        -1271./57600.)*i12e4 + (217./43200.)*i12e2 + (769369./3440640.)*i10e8 + (
                        -24905./86016.)*i10e6 + (697./3840.)*i10e4 + (-119./2880.)*i10e2 + (
                        16704927./26214400.)*i8e10 + (-407313./327680.)*i8e8 + (13185./8192.)*i8e6 + (
                        -2583./2560.)*i8e4 + (147./640.)*i8e2 + (47982879./52428800.)*i6e12 + (
                        -5568309./2621440.)*i6e10 + (135771./32768.)*i6e8 + (-21975./4096.)*i6e6 + (861./256.)*i6e4 + (
                        -49./64.)*i6e2 + (1602678489./3355443200.)*i4e14 + (-143948637./104857600.)*i4e12 + (
                        16704927./5242880.)*i4e10 + (-407313./65536.)*i4e8 + (65925./8192.)*i4e6 + (
                        -2583./512.)*i4e4 + (147./128.)*i4e2,
            None,
            -3.,
            -2.,
            0.
        ),
        (
            '-2*n',
            '2*n',
            '2, 0, 2, 0',
            -2*orbital_freq,
            (-257./2043241200.)*i18 + (-5461./454053600.)*i16e2 + (5461./2270268000.)*i16 + (-1./3520.)*i14e4 + (
                        1./5544.)*i14e2 + (-1./27720.)*i14 + (-589./340200.)*i12e6 + (31./9600.)*i12e4 + (
                        -31./15120.)*i12e2 + (31./75600.)*i12 + (-1003./331776.)*i10e8 + (323./22680.)*i10e6 + (
                        -17./640.)*i10e4 + (17./1008.)*i10e2 + (-17./5040.)*i10 + (4079./1024000.)*i8e10 + (
                        413./24576.)*i8e8 + (-19./240.)*i8e6 + (189./1280.)*i8e4 + (-3./32.)*i8e2 + (3./160.)*i8 + (
                        2323./1474560.)*i6e12 + (-4079./307200.)*i6e10 + (-2065./36864.)*i6e8 + (19./72.)*i6e6 + (
                        -63./128.)*i6e4 + (5./16.)*i6e2 + (-1./16.)*i6 + (689797./120422400.)*i4e14 + (
                        -2323./983040.)*i4e12 + (4079./204800.)*i4e10 + (2065./24576.)*i4e8 + (-19./48.)*i4e6 + (
                        189./256.)*i4e4 + (-15./32.)*i4e2 + (3./32.)*i4,
            None,
            -2.,
            -2.,
            0.
        ),
        (
            '-n',
            'n',
            '2, 0, 2, 1',
            -orbital_freq,
            (5461./9081072000.)*i16e2 + (1./443520.)*i14e4 + (-1./110880.)*i14e2 + (403./58060800.)*i12e6 + (
                        -31./1209600.)*i12e4 + (31./302400.)*i12e2 + (-1037./18579456.)*i10e8 + (
                        -221./3870720.)*i10e6 + (17./80640.)*i10e4 + (-17./20160.)*i10e2 + (5199./26214400.)*i8e10 + (
                        61./196608.)*i8e8 + (13./40960.)*i8e6 + (-3./2560.)*i8e4 + (3./640.)*i8e2 + (
                        -277229./471859200.)*i6e12 + (-1733./2621440.)*i6e10 + (-305./294912.)*i6e8 + (
                        -13./12288.)*i6e6 + (1./256.)*i6e4 + (-1./64.)*i6e2 + (165285343./211392921600.)*i4e14 + (
                        277229./314572800.)*i4e12 + (5199./5242880.)*i4e10 + (305./196608.)*i4e8 + (13./8192.)*i4e6 + (
                        -3./512.)*i4e4 + (3./128.)*i4e2,
            None,
            -1.,
            -2.,
            0.
        ),
        (
            'n',
            'n',
            '2, 0, 2, 3',
            orbital_freq,
            (31./174182400.)*i12e6 + (-187./92897280.)*i10e8 + (-17./11612160.)*i10e6 + (619./52428800.)*i8e10 + (
                        11./983040.)*i8e8 + (1./122880.)*i8e6 + (-62617./1698693120.)*i6e12 + (
                        -619./15728640.)*i6e10 + (-11./294912.)*i6e8 + (-1./36864.)*i6e6 + (
                        31398887./634178764800.)*i4e14 + (62617./1132462080.)*i4e12 + (619./10485760.)*i4e10 + (
                        11./196608.)*i4e8 + (1./24576.)*i4e6,
            None,
            1.,
            -2.,
            0.
        ),
        (
            '2*n',
            '2*n',
            '2, 0, 2, 4',
            2*orbital_freq,
            (-17./2903040.)*i10e8 + (7./153600.)*i8e10 + (1./30720.)*i8e8 + (-949./5529600.)*i6e12 + (
                        -7./46080.)*i6e10 + (-1./9216.)*i6e8 + (2417./9676800.)*i4e14 + (949./3686400.)*i4e12 + (
                        7./30720.)*i4e10 + (1./6144.)*i4e8,
            None,
            2.,
            -2.,
            0.
        ),
        (
            '3*n',
            '3*n',
            '2, 0, 2, 5',
            3*orbital_freq,
            (19683./262144000.)*i8e10 + (-6561./20971520.)*i6e12 + (-6561./26214400.)*i6e10 + (
                        12852999./23488102400.)*i4e14 + (19683./41943040.)*i4e12 + (19683./52428800.)*i4e10,
            None,
            3.,
            -2.,
            0.
        ),
        (
            '4*n',
            '4*n',
            '2, 0, 2, 6',
            4*orbital_freq,
            (-1./2025.)*i6e12 + (1./1350.)*i4e14 + (1./1350.)*i4e12,
            None,
            4.,
            -2.,
            0.
        ),
        (
            '5*n',
            '5*n',
            '2, 0, 2, 7',
            5*orbital_freq,
            (244140625./177570054144.)*i4e14,
            None,
            5.,
            -2.,
            0.
        ),
        (
            '-O - 6*n',
            'O + 6*n',
            '2, 1, 0, -8',
            -spin_freq - 6*orbital_freq,
            (531441./40140800.)*i2e16,
            None,
            -6.,
            2.,
            1.
        ),
        (
            '-O - 5*n',
            'O + 5*n',
            '2, 1, 0, -7',
            -spin_freq - 5*orbital_freq,
            (-1220703125./199766310912.)*i4e14 + (2685546875./532710162432.)*i2e16 + (244140625./33294385152.)*i2e14,
            None,
            -5.,
            2.,
            1.
        ),
        (
            '-O - 4*n',
            'O + 4*n',
            '2, 1, 0, -6',
            -spin_freq - 4*orbital_freq,
            (227./182250.)*i6e12 + (-4./1215.)*i4e14 + (-4./1215.)*i4e12 + (23./4725.)*i2e16 + (8./2025.)*i2e14 + (
                        8./2025.)*i2e12,
            None,
            -4.,
            2.,
            1.
        ),
        (
            '-O - 3*n',
            'O + 3*n',
            '2, 1, 0, -5',
            -spin_freq - 3*orbital_freq,
            (-21141./146800640.)*i8e10 + (165483./209715200.)*i6e12 + (165483./262144000.)*i6e10 + (
                        -1428111./587202560.)*i4e14 + (-2187./1048576.)*i4e12 + (-2187./1310720.)*i4e10 + (
                        9428157./3355443200.)*i2e16 + (4284333./1468006400.)*i2e14 + (6561./2621440.)*i2e12 + (
                        6561./3276800.)*i2e10,
            None,
            -3.,
            2.,
            1.
        ),
        (
            '-O - 2*n',
            'O + 2*n',
            '2, 1, 0, -4',
            -spin_freq - 2*orbital_freq,
            (40277./4180377600.)*i10e8 + (-29./331776.)*i8e10 + (-145./2322432.)*i8e8 + (215423./497664000.)*i6e12 + (
                        1589./4147200.)*i6e10 + (227./829440.)*i6e8 + (-2417./2177280.)*i4e14 + (
                        -949./829440.)*i4e12 + (-7./6912.)*i4e10 + (-5./6912.)*i4e8 + (22601./18579456.)*i2e16 + (
                        2417./1814400.)*i2e14 + (949./691200.)*i2e12 + (7./5760.)*i2e10 + (1./1152.)*i2e8,
            None,
            -2.,
            2.,
            1.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 0, -3',
            -spin_freq - orbital_freq,
            (-59123./220723937280.)*i12e6 + (443047./133772083200.)*i10e8 + (40277./16721510400.)*i10e6 + (
                        -17951./792723456.)*i8e10 + (-1595./74317824.)*i8e8 + (-145./9289728.)*i8e6 + (
                        14214059./152882380800.)*i6e12 + (140513./1415577600.)*i6e10 + (2497./26542080.)*i6e8 + (
                        227./3317760.)*i6e6 + (-31398887./142690222080.)*i4e14 + (-62617./254803968.)*i4e12 + (
                        -619./2359296.)*i4e10 + (-55./221184.)*i4e8 + (-5./27648.)*i4e6 + (
                        147400583./634178764800.)*i2e16 + (31398887./118908518400.)*i2e14 + (
                        62617./212336640.)*i2e12 + (619./1966080.)*i2e10 + (11./36864.)*i2e8 + (1./4608.)*i2e6,
            None,
            -1.,
            2.,
            1.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 0, -1',
            -spin_freq + orbital_freq,
            (-3490169./4184557977600.)*i16e2 + (-8988527./2789705318400.)*i14e4 + (8988527./697426329600.)*i14e2 + (
                        -768599./73574645760.)*i12e6 + (59123./1532805120.)*i12e4 + (-59123./383201280.)*i12e2 + (
                        2456897./26754416640.)*i10e8 + (523601./5573836800.)*i10e6 + (-40277./116121600.)*i10e4 + (
                        40277./29030400.)*i10e2 + (-50257./132120576.)*i8e10 + (-44225./74317824.)*i8e8 + (
                        -1885./3096576.)*i8e6 + (145./64512.)*i8e4 + (-145./16128.)*i8e2 + (
                        62930983./42467328000.)*i6e12 + (393391./235929600.)*i6e10 + (13847./5308416.)*i6e8 + (
                        2951./1105920.)*i6e6 + (-227./23040.)*i6e4 + (227./5760.)*i6e2 + (
                        -165285343./47563407360.)*i4e14 + (-277229./70778880.)*i4e12 + (-1733./393216.)*i4e10 + (
                        -1525./221184.)*i4e8 + (-65./9216.)*i4e6 + (5./192.)*i4e4 + (-5./48.)*i4e2 + (
                        49450862117./13317754060800.)*i2e16 + (165285343./39636172800.)*i2e14 + (
                        277229./58982400.)*i2e12 + (1733./327680.)*i2e10 + (305./36864.)*i2e8 + (13./1536.)*i2e6 + (
                        -1./32.)*i2e4 + (1./8.)*i2e2,
            None,
            1.,
            2.,
            1.
        ),
        (
            '-O + 2*n',
            'O - 2*n',
            '2, 1, 0, 0',
            -spin_freq + 2*orbital_freq,
            (2195943977./12804747411456000.)*i18 + (3490169./209227898880.)*i16e2 + (-3490169./1046139494400.)*i16 + (
                        8988527./22140518400.)*i14e4 + (-8988527./34871316480.)*i14e2 + (8988527./174356582400.)*i14 + (
                        4670717./1724405760.)*i12e6 + (-59123./12165120.)*i12e4 + (59123./19160064.)*i12e2 + (
                        -59123./95800320.)*i12 + (148904069./16721510400.)*i10e8 + (-3181883./130636800.)*i10e6 + (
                        40277./921600.)*i10e4 + (-40277./1451520.)*i10e2 + (40277./7257600.)*i10 + (
                        320189./15482880.)*i8e10 + (-536065./9289728.)*i8e8 + (11455./72576.)*i8e6 + (
                        -145./512.)*i8e4 + (725./4032.)*i8e2 + (-145./4032.)*i8 + (7824917./398131200.)*i6e12 + (
                        -2506307./27648000.)*i6e10 + (839219./3317760.)*i6e8 + (-17933./25920.)*i6e6 + (
                        1589./1280.)*i6e4 + (-227./288.)*i6e2 + (227./1440.)*i6 + (1739939./81285120.)*i4e14 + (
                        -34471./663552.)*i4e12 + (11041./46080.)*i4e10 + (-18485./27648.)*i4e8 + (395./216.)*i4e6 + (
                        -105./32.)*i4e4 + (25./12.)*i4e2 + (-5./12.)*i4 + (-561889./240844800.)*i2e16 + (
                        -1739939./67737600.)*i2e14 + (34471./552960.)*i2e12 + (-11041./38400.)*i2e10 + (
                        3697./4608.)*i2e8 + (-79./36.)*i2e6 + (63./16.)*i2e4 + (-5./2.)*i2e2 + (1./2.)*i2,
            None,
            2.,
            2.,
            1.
        ),
        (
            '-O + 3*n',
            'O - 3*n',
            '2, 1, 0, 1',
            -spin_freq + 3*orbital_freq,
            (-3490169./85399142400.)*i16e2 + (-368529607./132843110400.)*i14e4 + (8988527./14233190400.)*i14e2 + (
                        -86615195./1634992128.)*i12e6 + (2424043./72990720.)*i12e4 + (-413861./54743040.)*i12e2 + (
                        -1822816189./4954521600.)*i10e8 + (11801161./24772608.)*i10e6 + (-1651357./5529600.)*i10e4 + (
                        281939./4147200.)*i10e2 + (-17942329./14680064.)*i8e10 + (6562265./2752512.)*i8e8 + (
                        -1062125./344064.)*i8e6 + (5945./3072.)*i8e4 + (-1015./2304.)*i8e2 + (
                        -1210234837./524288000.)*i6e12 + (140445127./26214400.)*i6e10 + (-10273339./983040.)*i6e8 + (
                        332555./24576.)*i6e6 + (-65149./7680.)*i6e4 + (11123./5760.)*i6e2 + (
                        -534226163./251658240.)*i4e14 + (15994293./2621440.)*i4e12 + (-1856103./131072.)*i4e10 + (
                        226285./8192.)*i4e8 + (-36625./1024.)*i4e6 + (1435./64.)*i4e4 + (-245./48.)*i4e2 + (
                        -3973253733./4697620480.)*i2e16 + (534226163./209715200.)*i2e14 + (
                        -47982879./6553600.)*i2e12 + (5568309./327680.)*i2e10 + (-135771./4096.)*i2e8 + (
                        21975./512.)*i2e6 + (-861./32.)*i2e4 + (49./8.)*i2e2,
            None,
            3.,
            2.,
            1.
        ),
        (
            '-O + 4*n',
            'O - 4*n',
            '2, 1, 0, 2',
            -spin_freq + 4*orbital_freq,
            (2597684303./697426329600.)*i14e4 + (23117093./114960384.)*i12e6 + (-17086547./383201280.)*i12e4 + (
                        3365183627./1045094400.)*i10e8 + (-15748307./8709120.)*i10e6 + (11640053./29030400.)*i10e4 + (
                        3892061./193536.)*i8e10 + (-12114895./580608.)*i8e8 + (283475./24192.)*i8e6 + (
                        -41905./16128.)*i8e4 + (1911732937./33177600.)*i6e12 + (-30465443./345600.)*i6e10 + (
                        18966077./207360.)*i6e8 + (-88757./1728.)*i6e6 + (65603./5760.)*i6e4 + (
                        346700573./4354560.)*i4e14 + (-8421731./55296.)*i4e12 + (134209./576.)*i4e10 + (
                        -417755./1728.)*i4e8 + (9775./72.)*i4e6 + (-1445./48.)*i4e4 + (173370469./4147200.)*i2e16 + (
                        -346700573./3628800.)*i2e14 + (8421731./46080.)*i2e12 + (-134209./480.)*i2e10 + (
                        83551./288.)*i2e8 + (-1955./12.)*i2e6 + (289./8.)*i2e4,
            None,
            4.,
            2.,
            1.
        ),
        (
            '-O + 5*n',
            'O - 5*n',
            '2, 1, 0, 3',
            -spin_freq + 5*orbital_freq,
            (-8443060015./44144787456.)*i12e6 + (-44278318565./5350883328.)*i10e8 + (1150351397./668860416.)*i10e6 + (
                        -85147679375./792723456.)*i8e10 + (3985125625./74317824.)*i8e8 + (-103533625./9289728.)*i8e6 + (
                        -3299263235425./6115295232.)*i6e12 + (26660032025./56623104.)*i6e10 + (
                        -1247756575./5308416.)*i6e8 + (32416735./663552.)*i6e6 + (
                        -32147077961875./28538044416.)*i4e14 + (363354981875./254803968.)*i4e12 + (
                        -2936126875./2359296.)*i4e10 + (137418125./221184.)*i4e8 + (-3570125./27648.)*i4e6 + (
                        -21274828753525./25367150592.)*i2e16 + (6429415592375./4756340736.)*i2e14 + (
                        -72670996375./42467328.)*i2e12 + (587225375./393216.)*i2e10 + (-27483625./36864.)*i2e8 + (
                        714025./4608.)*i2e6,
            None,
            5.,
            2.,
            1.
        ),
        (
            '-O + 6*n',
            'O - 6*n',
            '2, 1, 0, 4',
            -spin_freq + 6*orbital_freq,
            (11442252653./1857945600.)*i10e8 + (71241313./344064.)*i8e10 + (-41192905./1032192.)*i8e8 + (
                        49930187017./24576000.)*i6e12 + (-557647519./614400.)*i6e10 + (64488203./368640.)*i6e8 + (
                        1018573919./143360.)*i4e14 + (-219956771./40960.)*i4e12 + (2456597./1024.)*i4e10 + (
                        -1420445./3072.)*i4e8 + (180544031973./22937600.)*i2e16 + (-3055721757./358400.)*i2e14 + (
                        659870313./102400.)*i2e12 + (-7369791./2560.)*i2e10 + (284089./512.)*i2e8,
            None,
            6.,
            2.,
            1.
        ),
        (
            '-O + 7*n',
            'O - 7*n',
            '2, 1, 0, 5',
            -spin_freq + 7*orbital_freq,
            (-216018317123./1698693120.)*i8e10 + (-6367538528267./2038431744.)*i6e12 + (
                        11836313996843./21233664000.)*i6e10 + (-417304029320899./20384317440.)*i4e14 + (
                        701270763025./84934656.)*i4e12 + (-52142352409./35389440.)*i4e10 + (
                        -9997999669389767./271790899200.)*i2e16 + (417304029320899./16986931200.)*i2e14 + (
                        -140254152605./14155776.)*i2e12 + (52142352409./29491200.)*i2e10,
            None,
            7.,
            2.,
            1.
        ),
        (
            '-O + 8*n',
            'O - 8*n',
            '2, 1, 0, 6',
            -spin_freq + 8*orbital_freq,
            (1221943306547./746496000.)*i6e12 + (56914314263./2177280.)*i4e14 + (-5383010161./1244160.)*i4e12 + (
                        578658802849./6773760.)*i2e16 + (-56914314263./1814400.)*i2e14 + (5383010161./1036800.)*i2e12,
            None,
            8.,
            2.,
            1.
        ),
        (
            '-O + 9*n',
            'O - 9*n',
            '2, 1, 0, 7',
            -spin_freq + 9*orbital_freq,
            (-49161122232843./4110417920.)*i4e14 + (-3064748517300717./32883343360.)*i2e16 + (
                        147483366698529./10276044800.)*i2e14,
            None,
            9.,
            2.,
            1.
        ),
        (
            '-O + 10*n',
            'O - 10*n',
            '2, 1, 0, 8',
            -spin_freq + 10*orbital_freq,
            (402063787225./10616832.)*i2e16,
            None,
            10.,
            2.,
            1.
        ),
        (
            '-O - 8*n',
            'O + 8*n',
            '2, 1, 1, -8',
            -spin_freq - 8*orbital_freq,
            (31801945561./160563200.)*i2e16,
            None,
            -8.,
            0.,
            1.
        ),
        (
            '-O - 7*n',
            'O + 7*n',
            '2, 1, 1, -7',
            -spin_freq - 7*orbital_freq,
            (-186702632281./1415577600.)*i4e14 + (-1328179895713./10066329600.)*i2e16 + (
                        186702632281./1887436800.)*i2e14,
            None,
            -7.,
            0.,
            1.
        ),
        (
            '-O - 6*n',
            'O + 6*n',
            '2, 1, 1, -6',
            -spin_freq - 6*orbital_freq,
            (10029889./288000.)*i6e12 + (2308743./44800.)*i4e14 + (-10029889./153600.)*i4e12 + (
                        197537355./3211264.)*i2e16 + (-6926229./179200.)*i2e14 + (10029889./204800.)*i2e12,
            None,
            -6.,
            0.,
            1.
        ),
        (
            '-O - 5*n',
            'O + 5*n',
            '2, 1, 1, -5',
            -spin_freq - 5*orbital_freq,
            (-349281./71680.)*i8e10 + (-982439./245760.)*i6e12 + (349281./20480.)*i6e10 + (
                        -13524196093./396361728.)*i4e14 + (982439./131072.)*i4e12 + (-1047843./32768.)*i4e10 + (
                        28491858875./2818572288.)*i2e16 + (13524196093./528482304.)*i2e14 + (
                        -2947317./524288.)*i2e12 + (3143529./131072.)*i2e10,
            None,
            -5.,
            0.,
            1.
        ),
        (
            '-O - 4*n',
            'O + 4*n',
            '2, 1, 1, -4',
            -spin_freq - 4*orbital_freq,
            (847./2025.)*i10e8 + (-473./600.)*i8e10 + (-847./360.)*i8e8 + (708871./72000.)*i6e12 + (
                        3311./1200.)*i6e10 + (5929./720.)*i6e8 + (-3027367./172800.)*i4e14 + (-708871./38400.)*i4e12 + (
                        -3311./640.)*i4e10 + (-5929./384.)*i4e8 + (124914751./6881280.)*i2e16 + (
                        3027367./230400.)*i2e14 + (708871./51200.)*i2e12 + (9933./2560.)*i2e10 + (5929./512.)*i2e8,
            None,
            -4.,
            0.,
            1.
        ),
        (
            '-O - 3*n',
            'O + 3*n',
            '2, 1, 1, -3',
            -spin_freq - 3*orbital_freq,
            (-11236./467775.)*i12e6 + (6943./37800.)*i10e8 + (2809./14175.)*i10e6 + (-286661./153600.)*i8e10 + (
                        -6943./6720.)*i8e8 + (-2809./2520.)*i8e6 + (30355211./3686400.)*i6e12 + (
                        2006627./307200.)*i6e10 + (6943./1920.)*i6e8 + (2809./720.)*i6e6 + (
                        -1055142483./52428800.)*i4e14 + (-30355211./1966080.)*i4e12 + (-2006627./163840.)*i4e10 + (
                        -6943./1024.)*i4e8 + (-2809./384.)*i4e6 + (437242738059./23488102400.)*i2e16 + (
                        3165427449./209715200.)*i2e14 + (30355211./2621440.)*i2e12 + (6019881./655360.)*i2e10 + (
                        20829./4096.)*i2e8 + (2809./512.)*i2e6,
            None,
            -3.,
            0.,
            1.
        ),
        (
            '-O - 2*n',
            'O + 2*n',
            '2, 1, 1, -2',
            -spin_freq - 2*orbital_freq,
            (512./525525.)*i14e4 + (-128./7425.)*i12e6 + (-64./5775.)*i12e4 + (3322./14175.)*i10e8 + (
                        32./225.)*i10e6 + (16./175.)*i10e4 + (-3919./2100.)*i8e10 + (-1661./1260.)*i8e8 + (
                        -4./5.)*i8e6 + (-18./35.)*i8e4 + (505601./57600.)*i6e12 + (3919./600.)*i6e10 + (
                        1661./360.)*i6e8 + (14./5.)*i6e6 + (9./5.)*i6e4 + (-25565893./1209600.)*i4e14 + (
                        -505601./30720.)*i4e12 + (-3919./320.)*i4e10 + (-1661./192.)*i4e8 + (-21./4.)*i4e6 + (
                        -27./8.)*i4e4 + (678544541./34406400.)*i2e16 + (25565893./1612800.)*i2e14 + (
                        505601./40960.)*i2e12 + (11757./1280.)*i2e10 + (1661./256.)*i2e8 + (63./16.)*i2e6 + (
                        81./32.)*i2e4,
            None,
            -2.,
            0.,
            1.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 1, -1',
            -spin_freq - orbital_freq,
            (-2048./70945875.)*i16e2 + (512./525525.)*i14e4 + (2048./4729725.)*i14e2 + (-68./3465.)*i12e6 + (
                        -64./5775.)*i12e4 + (-256./51975.)*i12e2 + (28211./113400.)*i10e8 + (17./105.)*i10e6 + (
                        16./175.)*i10e4 + (64./1575.)*i10e2 + (-1064887./537600.)*i8e10 + (-28211./20160.)*i8e8 + (
                        -51./56.)*i8e6 + (-18./35.)*i8e4 + (-8./35.)*i8e2 + (85553819./9216000.)*i6e12 + (
                        1064887./153600.)*i6e10 + (28211./5760.)*i6e8 + (51./16.)*i6e6 + (9./5.)*i6e4 + (4./5.)*i6e2 + (
                        -24654653741./1101004800.)*i4e14 + (-85553819./4915200.)*i4e12 + (-1064887./81920.)*i4e10 + (
                        -28211./3072.)*i4e8 + (-765./128.)*i4e6 + (-27./8.)*i4e4 + (-3./2.)*i4e2 + (
                        689312857627./32883343360.)*i2e16 + (24654653741./1468006400.)*i2e14 + (
                        85553819./6553600.)*i2e12 + (3194661./327680.)*i2e10 + (28211./4096.)*i2e8 + (
                        2295./512.)*i2e6 + (81./32.)*i2e4 + (9./8.)*i2e2,
            None,
            -1.,
            0.,
            1.
        ),
        (
            '-O',
            'O',
            '2, 1, 1, 0',
            -spin_freq,
            (65536./97692469875.)*i18 + (-8192./212837625.)*i16e2 + (-8192./638512875.)*i16 + (
                        16384./14189175.)*i14e4 + (8192./14189175.)*i14e2 + (8192./42567525.)*i14 + (
                        -2048./93555.)*i12e6 + (-2048./155925.)*i12e4 + (-1024./155925.)*i12e2 + (
                        -1024./467775.)*i12 + (256./945.)*i10e8 + (512./2835.)*i10e6 + (512./4725.)*i10e4 + (
                        256./4725.)*i10e2 + (256./14175.)*i10 + (-32./15.)*i8e10 + (-32./21.)*i8e8 + (-64./63.)*i8e6 + (
                        -64./105.)*i8e4 + (-32./105.)*i8e2 + (-32./315.)*i8 + (448./45.)*i6e12 + (112./15.)*i6e10 + (
                        16./3.)*i6e8 + (32./9.)*i6e6 + (32./15.)*i6e4 + (16./15.)*i6e2 + (16./45.)*i6 + -24.*i4e14 + (
                        -56./3.)*i4e12 + -14.*i4e10 + -10.*i4e8 + (-20./3.)*i4e6 + -4.*i4e4 + -2.*i4e2 + (-2./3.)*i4 + (
                        45./2.)*i2e16 + 18.*i2e14 + 14.*i2e12 + (21./2.)*i2e10 + (15./2.)*i2e8 + 5.*i2e6 + 3.*i2e4 + (
                        3./2.)*i2e2 + (1./2.)*i2,
            None,
            0.,
            0.,
            1.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 1, 1',
            -spin_freq + orbital_freq,
            (-2048./70945875.)*i16e2 + (512./525525.)*i14e4 + (2048./4729725.)*i14e2 + (-68./3465.)*i12e6 + (
                        -64./5775.)*i12e4 + (-256./51975.)*i12e2 + (28211./113400.)*i10e8 + (17./105.)*i10e6 + (
                        16./175.)*i10e4 + (64./1575.)*i10e2 + (-1064887./537600.)*i8e10 + (-28211./20160.)*i8e8 + (
                        -51./56.)*i8e6 + (-18./35.)*i8e4 + (-8./35.)*i8e2 + (85553819./9216000.)*i6e12 + (
                        1064887./153600.)*i6e10 + (28211./5760.)*i6e8 + (51./16.)*i6e6 + (9./5.)*i6e4 + (4./5.)*i6e2 + (
                        -24654653741./1101004800.)*i4e14 + (-85553819./4915200.)*i4e12 + (-1064887./81920.)*i4e10 + (
                        -28211./3072.)*i4e8 + (-765./128.)*i4e6 + (-27./8.)*i4e4 + (-3./2.)*i4e2 + (
                        689312857627./32883343360.)*i2e16 + (24654653741./1468006400.)*i2e14 + (
                        85553819./6553600.)*i2e12 + (3194661./327680.)*i2e10 + (28211./4096.)*i2e8 + (
                        2295./512.)*i2e6 + (81./32.)*i2e4 + (9./8.)*i2e2,
            None,
            1.,
            0.,
            1.
        ),
        (
            '-O + 2*n',
            'O - 2*n',
            '2, 1, 1, 2',
            -spin_freq + 2*orbital_freq,
            (512./525525.)*i14e4 + (-128./7425.)*i12e6 + (-64./5775.)*i12e4 + (3322./14175.)*i10e8 + (
                        32./225.)*i10e6 + (16./175.)*i10e4 + (-3919./2100.)*i8e10 + (-1661./1260.)*i8e8 + (
                        -4./5.)*i8e6 + (-18./35.)*i8e4 + (505601./57600.)*i6e12 + (3919./600.)*i6e10 + (
                        1661./360.)*i6e8 + (14./5.)*i6e6 + (9./5.)*i6e4 + (-25565893./1209600.)*i4e14 + (
                        -505601./30720.)*i4e12 + (-3919./320.)*i4e10 + (-1661./192.)*i4e8 + (-21./4.)*i4e6 + (
                        -27./8.)*i4e4 + (678544541./34406400.)*i2e16 + (25565893./1612800.)*i2e14 + (
                        505601./40960.)*i2e12 + (11757./1280.)*i2e10 + (1661./256.)*i2e8 + (63./16.)*i2e6 + (
                        81./32.)*i2e4,
            None,
            2.,
            0.,
            1.
        ),
        (
            '-O + 3*n',
            'O - 3*n',
            '2, 1, 1, 3',
            -spin_freq + 3*orbital_freq,
            (-11236./467775.)*i12e6 + (6943./37800.)*i10e8 + (2809./14175.)*i10e6 + (-286661./153600.)*i8e10 + (
                        -6943./6720.)*i8e8 + (-2809./2520.)*i8e6 + (30355211./3686400.)*i6e12 + (
                        2006627./307200.)*i6e10 + (6943./1920.)*i6e8 + (2809./720.)*i6e6 + (
                        -1055142483./52428800.)*i4e14 + (-30355211./1966080.)*i4e12 + (-2006627./163840.)*i4e10 + (
                        -6943./1024.)*i4e8 + (-2809./384.)*i4e6 + (437242738059./23488102400.)*i2e16 + (
                        3165427449./209715200.)*i2e14 + (30355211./2621440.)*i2e12 + (6019881./655360.)*i2e10 + (
                        20829./4096.)*i2e8 + (2809./512.)*i2e6,
            None,
            3.,
            0.,
            1.
        ),
        (
            '-O + 4*n',
            'O - 4*n',
            '2, 1, 1, 4',
            -spin_freq + 4*orbital_freq,
            (847./2025.)*i10e8 + (-473./600.)*i8e10 + (-847./360.)*i8e8 + (708871./72000.)*i6e12 + (
                        3311./1200.)*i6e10 + (5929./720.)*i6e8 + (-3027367./172800.)*i4e14 + (-708871./38400.)*i4e12 + (
                        -3311./640.)*i4e10 + (-5929./384.)*i4e8 + (124914751./6881280.)*i2e16 + (
                        3027367./230400.)*i2e14 + (708871./51200.)*i2e12 + (9933./2560.)*i2e10 + (5929./512.)*i2e8,
            None,
            4.,
            0.,
            1.
        ),
        (
            '-O + 5*n',
            'O - 5*n',
            '2, 1, 1, 5',
            -spin_freq + 5*orbital_freq,
            (-349281./71680.)*i8e10 + (-982439./245760.)*i6e12 + (349281./20480.)*i6e10 + (
                        -13524196093./396361728.)*i4e14 + (982439./131072.)*i4e12 + (-1047843./32768.)*i4e10 + (
                        28491858875./2818572288.)*i2e16 + (13524196093./528482304.)*i2e14 + (
                        -2947317./524288.)*i2e12 + (3143529./131072.)*i2e10,
            None,
            5.,
            0.,
            1.
        ),
        (
            '-O + 6*n',
            'O - 6*n',
            '2, 1, 1, 6',
            -spin_freq + 6*orbital_freq,
            (10029889./288000.)*i6e12 + (2308743./44800.)*i4e14 + (-10029889./153600.)*i4e12 + (
                        197537355./3211264.)*i2e16 + (-6926229./179200.)*i2e14 + (10029889./204800.)*i2e12,
            None,
            6.,
            0.,
            1.
        ),
        (
            '-O + 7*n',
            'O - 7*n',
            '2, 1, 1, 7',
            -spin_freq + 7*orbital_freq,
            (-186702632281./1415577600.)*i4e14 + (-1328179895713./10066329600.)*i2e16 + (
                        186702632281./1887436800.)*i2e14,
            None,
            7.,
            0.,
            1.
        ),
        (
            '-O + 8*n',
            'O - 8*n',
            '2, 1, 1, 8',
            -spin_freq + 8*orbital_freq,
            (31801945561./160563200.)*i2e16,
            None,
            8.,
            0.,
            1.
        ),
        (
            '-O - 8*n',
            'O + 8*n',
            '2, 1, 2, -6',
            -spin_freq - 8*orbital_freq,
            (5383010161./16588800.)*i6e12,
            None,
            -8.,
            -2.,
            1.
        ),
        (
            '-O - 7*n',
            'O + 7*n',
            '2, 1, 2, -5',
            -spin_freq - 7*orbital_freq,
            (-52142352409./943718400.)*i8e10 + (-140254152605./226492416.)*i6e12 + (52142352409./471859200.)*i6e10,
            None,
            -7.,
            -2.,
            1.
        ),
        (
            '-O - 6*n',
            'O + 6*n',
            '2, 1, 2, -4',
            -spin_freq - 6*orbital_freq,
            (2556801./655360.)*i10e8 + (7369791./81920.)*i8e10 + (-284089./16384.)*i8e8 + (
                        659870313./1638400.)*i6e12 + (-7369791./40960.)*i6e10 + (284089./8192.)*i6e8,
            None,
            -6.,
            -2.,
            1.
        ),
        (
            '-O - 5*n',
            'O + 5*n',
            '2, 1, 2, -3',
            -spin_freq - 5*orbital_freq,
            (-66118715./445906944.)*i12e6 + (-5496725./1048576.)*i10e8 + (142805./131072.)*i10e6 + (
                        -587225375./12582912.)*i8e10 + (27483625./1179648.)*i8e8 + (-714025./147456.)*i8e6 + (
                        -72670996375./679477248.)*i6e12 + (587225375./6291456.)*i6e10 + (-27483625./589824.)*i6e8 + (
                        714025./73728.)*i6e6,
            None,
            -5.,
            -2.,
            1.
        ),
        (
            '-O - 4*n',
            'O + 4*n',
            '2, 1, 2, -2',
            -spin_freq - 4*orbital_freq,
            (49997./15482880.)*i14e4 + (181033./1161216.)*i12e6 + (-133807./3870720.)*i12e4 + (83551./40960.)*i10e8 + (
                        -1173./1024.)*i10e6 + (2601./10240.)*i10e4 + (134209./15360.)*i8e10 + (-83551./9216.)*i8e8 + (
                        1955./384.)*i8e6 + (-289./256.)*i8e4 + (8421731./737280.)*i6e12 + (-134209./7680.)*i6e10 + (
                        83551./4608.)*i6e8 + (-1955./192.)*i6e6 + (289./128.)*i6e4,
            None,
            -4.,
            -2.,
            1.
        ),
        (
            '-O - 3*n',
            'O + 3*n',
            '2, 1, 2, -1',
            -spin_freq - 3*orbital_freq,
            (-3059./81100800.)*i16e2 + (-7093./2949120.)*i14e4 + (1211./2211840.)*i14e2 + (-678295./16515072.)*i12e6 + (
                        18983./737280.)*i12e4 + (-3241./552960.)*i12e2 + (-1221939./5242880.)*i10e8 + (
                        39555./131072.)*i10e6 + (-7749./40960.)*i10e4 + (441./10240.)*i10e2 + (
                        -5568309./10485760.)*i8e10 + (135771./131072.)*i8e8 + (-21975./16384.)*i8e6 + (
                        861./1024.)*i8e4 + (-49./256.)*i8e2 + (-47982879./104857600.)*i6e12 + (
                        5568309./5242880.)*i6e10 + (-135771./65536.)*i6e8 + (21975./8192.)*i6e6 + (-861./512.)*i6e4 + (
                        49./128.)*i6e2,
            None,
            -3.,
            -2.,
            1.
        ),
        (
            '-O - 2*n',
            'O + 2*n',
            '2, 1, 2, 0',
            -spin_freq - 2*orbital_freq,
            (2743907./16738231910400.)*i18 + (437./28385280.)*i16e2 + (-437./141926400.)*i16 + (173./491520.)*i14e4 + (
                        -173./774144.)*i14e2 + (173./3870720.)*i14 + (8797./4354560.)*i12e6 + (-463./122880.)*i12e4 + (
                        463./193536.)*i12e2 + (-463./967680.)*i12 + (413./131072.)*i10e8 + (-19./1280.)*i10e6 + (
                        567./20480.)*i10e4 + (-9./512.)*i10e2 + (9./2560.)*i10 + (-4079./1228800.)*i8e10 + (
                        -2065./147456.)*i8e8 + (19./288.)*i8e6 + (-63./512.)*i8e4 + (5./64.)*i8e2 + (-1./64.)*i8 + (
                        -2323./2949120.)*i6e12 + (4079./614400.)*i6e10 + (2065./73728.)*i6e8 + (-19./144.)*i6e6 + (
                        63./256.)*i6e4 + (-5./32.)*i6e2 + (1./32.)*i6,
            None,
            -2.,
            -2.,
            1.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 2, 1',
            -spin_freq - orbital_freq,
            (-437./567705600.)*i16e2 + (-173./61931520.)*i14e4 + (173./15482880.)*i14e2 + (-6019./743178240.)*i12e6 + (
                        463./15482880.)*i12e4 + (-463./3870720.)*i12e2 + (61./1048576.)*i10e8 + (39./655360.)*i10e6 + (
                        -9./40960.)*i10e4 + (9./10240.)*i10e2 + (-1733./10485760.)*i8e10 + (-305./1179648.)*i8e8 + (
                        -13./49152.)*i8e6 + (1./1024.)*i8e4 + (-1./256.)*i8e2 + (277229./943718400.)*i6e12 + (
                        1733./5242880.)*i6e10 + (305./589824.)*i6e8 + (13./24576.)*i6e6 + (-1./512.)*i6e4 + (
                        1./128.)*i6e2,
            None,
            -1.,
            -2.,
            1.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 2, 3',
            -spin_freq + orbital_freq,
            (-463./2229534720.)*i12e6 + (11./5242880.)*i10e8 + (1./655360.)*i10e6 + (-619./62914560.)*i8e10 + (
                        -11./1179648.)*i8e8 + (-1./147456.)*i8e6 + (62617./3397386240.)*i6e12 + (
                        619./31457280.)*i6e10 + (11./589824.)*i6e8 + (1./73728.)*i6e6,
            None,
            1.,
            -2.,
            1.
        ),
        (
            '-O + 2*n',
            'O - 2*n',
            '2, 1, 2, 4',
            -spin_freq + 2*orbital_freq,
            (1./163840.)*i10e8 + (-7./184320.)*i8e10 + (-1./36864.)*i8e8 + (949./11059200.)*i6e12 + (
                        7./92160.)*i6e10 + (1./18432.)*i6e8,
            None,
            2.,
            -2.,
            1.
        ),
        (
            '-O + 3*n',
            'O - 3*n',
            '2, 1, 2, 5',
            -spin_freq + 3*orbital_freq,
            (-6561./104857600.)*i8e10 + (6561./41943040.)*i6e12 + (6561./52428800.)*i6e10,
            None,
            3.,
            -2.,
            1.
        ),
        (
            '-O + 4*n',
            'O - 4*n',
            '2, 1, 2, 6',
            -spin_freq + 4*orbital_freq,
            (1./4050.)*i6e12,
            None,
            4.,
            -2.,
            1.
        ),
        (
            '-2*O - 7*n',
            '2*O + 7*n',
            '2, 2, 0, -9',
            -2*spin_freq - 7*orbital_freq,
            (33232930569601./1408964021452800.)*e18,
            None,
            -7.,
            2.,
            2.
        ),
        (
            '-2*O - 6*n',
            '2*O + 6*n',
            '2, 2, 0, -8',
            -2*spin_freq - 6*orbital_freq,
            (-531441./40140800.)*i2e16 + (177147./40140800.)*e18 + (531441./40140800.)*e16,
            None,
            -6.,
            2.,
            2.
        ),
        (
            '-2*O - 5*n',
            '2*O + 5*n',
            '2, 2, 0, -7',
            -2*spin_freq - 5*orbital_freq,
            (2685546875./799065243648.)*i4e14 + (-2685546875./532710162432.)*i2e16 + (
                        -244140625./33294385152.)*i2e14 + (323974609375./43834436222976.)*e18 + (
                        2685546875./532710162432.)*e16 + (244140625./33294385152.)*e14,
            None,
            -5.,
            2.,
            2.
        ),
        (
            '-2*O - 4*n',
            '2*O + 4*n',
            '2, 2, 0, -6',
            -2*spin_freq - 4*orbital_freq,
            (-46./91125.)*i6e12 + (11./6075.)*i4e14 + (11./6075.)*i4e12 + (-23./4725.)*i2e16 + (-8./2025.)*i2e14 + (
                        -8./2025.)*i2e12 + (68./15309.)*e18 + (23./4725.)*e16 + (8./2025.)*e14 + (8./2025.)*e12,
            None,
            -4.,
            2.,
            2.
        ),
        (
            '-2*O - 3*n',
            '2*O + 3*n',
            '2, 2, 0, -5',
            -2*spin_freq - 3*orbital_freq,
            (1426653./29360128000.)*i8e10 + (-16767./52428800.)*i6e12 + (-16767./65536000.)*i6e10 + (
                        15709221./11744051200.)*i4e14 + (24057./20971520.)*i4e12 + (24057./26214400.)*i4e10 + (
                        -9428157./3355443200.)*i2e16 + (-4284333./1468006400.)*i2e14 + (-6561./2621440.)*i2e12 + (
                        -6561./3276800.)*i2e10 + (270241029./105226698752.)*e18 + (9428157./3355443200.)*e16 + (
                        4284333./1468006400.)*e14 + (6561./2621440.)*e12 + (6561./3276800.)*e10,
            None,
            -3.,
            2.,
            2.
        ),
        (
            '-2*O - 2*n',
            '2*O + 2*n',
            '2, 2, 0, -4',
            -2*spin_freq - 2*orbital_freq,
            (-12107./4180377600.)*i10e8 + (1957./66355200.)*i8e10 + (1957./92897280.)*i8e8 + (
                        -21827./124416000.)*i6e12 + (-161./1036800.)*i6e10 + (-23./207360.)*i6e8 + (
                        26587./43545600.)*i4e14 + (10439./16588800.)*i4e12 + (77./138240.)*i4e10 + (11./27648.)*i4e8 + (
                        -22601./18579456.)*i2e16 + (-2417./1814400.)*i2e14 + (-949./691200.)*i2e12 + (
                        -7./5760.)*i2e10 + (-1./1152.)*i2e8 + (128441./119439360.)*e18 + (22601./18579456.)*e16 + (
                        2417./1814400.)*e14 + (949./691200.)*e12 + (7./5760.)*e10 + (1./1152.)*e8,
            None,
            -2.,
            2.,
            2.
        ),
        (
            '-2*O - n',
            '2*O + n',
            '2, 2, 0, -3',
            -2*spin_freq - orbital_freq,
            (330367./4414478745600.)*i12e6 + (-133177./133772083200.)*i10e8 + (-12107./16721510400.)*i10e6 + (
                        1211383./158544691200.)*i8e10 + (21527./2972712960.)*i8e8 + (1957./371589120.)*i8e6 + (
                        -1440191./38220595200.)*i6e12 + (-14237./353894400.)*i6e10 + (-253./6635520.)*i6e8 + (
                        -23./829440.)*i6e6 + (345387757./2853804441600.)*i4e14 + (688787./5096079360.)*i4e12 + (
                        6809./47185920.)*i4e10 + (121./884736.)*i4e8 + (11./110592.)*i4e6 + (
                        -147400583./634178764800.)*i2e16 + (-31398887./118908518400.)*i2e14 + (
                        -62617./212336640.)*i2e12 + (-619./1966080.)*i2e10 + (-11./36864.)*i2e8 + (-1./4608.)*i2e6 + (
                        167431204877./821895679180800.)*e18 + (147400583./634178764800.)*e16 + (
                        31398887./118908518400.)*e14 + (62617./212336640.)*e12 + (619./1966080.)*e10 + (
                        11./36864.)*e8 + (1./4608.)*e6,
            None,
            -1.,
            2.,
            2.
        ),
        (
            '-2*O + n',
            '2*O - n',
            '2, 2, 0, -1',
            -2*spin_freq + orbital_freq,
            (72518377./334764638208000.)*i16e2 + (27269./31701196800.)*i14e4 + (-27269./7925299200.)*i14e2 + (
                        4294771./1471492915200.)*i12e6 + (-330367./30656102400.)*i12e4 + (330367./7664025600.)*i12e2 + (
                        -738527./26754416640.)*i10e8 + (-157391./5573836800.)*i10e6 + (12107./116121600.)*i10e4 + (
                        -12107./29030400.)*i10e2 + (3391481./26424115200.)*i8e10 + (119377./594542592.)*i8e8 + (
                        25441./123863040.)*i8e6 + (-1957./2580480.)*i8e4 + (1957./645120.)*i8e2 + (
                        -6376267./10616832000.)*i6e12 + (-39859./58982400.)*i6e10 + (-1403./1327104.)*i6e8 + (
                        -299./276480.)*i6e6 + (23./5760.)*i6e4 + (-23./1440.)*i6e2 + (
                        1818138773./951268147200.)*i4e14 + (3049519./1415577600.)*i4e12 + (19063./7864320.)*i4e10 + (
                        3355./884736.)*i4e8 + (143./36864.)*i4e6 + (-11./768.)*i4e4 + (11./192.)*i4e2 + (
                        -49450862117./13317754060800.)*i2e16 + (-165285343./39636172800.)*i2e14 + (
                        -277229./58982400.)*i2e12 + (-1733./327680.)*i2e10 + (-305./36864.)*i2e8 + (-13./1536.)*i2e6 + (
                        1./32.)*i2e4 + (-1./8.)*i2e2 + (12842565048623./3835513169510400.)*e18 + (
                        49450862117./13317754060800.)*e16 + (165285343./39636172800.)*e14 + (277229./58982400.)*e12 + (
                        1733./327680.)*e10 + (305./36864.)*e8 + (13./1536.)*e6 + (-1./32.)*e4 + (1./8.)*e2,
            None,
            1.,
            2.,
            2.
        ),
        (
            '-2*O + 2*n',
            '2*O - 2*n',
            '2, 2, 0, 0',
            -2*spin_freq + 2*orbital_freq,
            (-561142037./12804747411456000.)*i18 + (-72518377./16738231910400.)*i16e2 + (
                        72518377./83691159552000.)*i16 + (-27269./251596800.)*i14e4 + (27269./396264960.)*i14e2 + (
                        -27269./1981324800.)*i14 + (-26098993./34488115200.)*i12e6 + (330367./243302400.)*i12e4 + (
                        -330367./383201280.)*i12e2 + (330367./1916006400.)*i12 + (-44759579./16721510400.)*i10e8 + (
                        956453./130636800.)*i10e6 + (-12107./921600.)*i10e4 + (12107./1451520.)*i10e2 + (
                        -12107./7257600.)*i10 + (-21607237./3096576000.)*i8e10 + (7235029./371589120.)*i8e8 + (
                        -154603./2903040.)*i8e6 + (1957./20480.)*i8e4 + (-1957./32256.)*i8e2 + (1957./161280.)*i8 + (
                        -792833./99532800.)*i6e12 + (253943./6912000.)*i6e10 + (-85031./829440.)*i6e8 + (
                        1817./6480.)*i6e6 + (-161./320.)*i6e4 + (23./72.)*i6e2 + (-23./360.)*i6 + (
                        -19139329./1625702400.)*i4e14 + (379181./13271040.)*i4e12 + (-121451./921600.)*i4e10 + (
                        40667./110592.)*i4e8 + (-869./864.)*i4e6 + (231./128.)*i4e4 + (-55./48.)*i4e2 + (11./48.)*i4 + (
                        561889./240844800.)*i2e16 + (1739939./67737600.)*i2e14 + (-34471./552960.)*i2e12 + (
                        11041./38400.)*i2e10 + (-3697./4608.)*i2e8 + (79./36.)*i2e6 + (-63./16.)*i2e4 + (5./2.)*i2e2 + (
                        -1./2.)*i2 + (-459927151./75246796800.)*e18 + (-561889./240844800.)*e16 + (
                        -1739939./67737600.)*e14 + (34471./552960.)*e12 + (-11041./38400.)*e10 + (3697./4608.)*e8 + (
                        -79./36.)*e6 + (63./16.)*e4 + (-5./2.)*e2 + (1./2.),
            None,
            2.,
            2.,
            2.
        ),
        (
            '-2*O + 3*n',
            '2*O - 3*n',
            '2, 2, 0, 1',
            -2*spin_freq + 3*orbital_freq,
            (72518377./6831931392000.)*i16e2 + (1118029./1509580800.)*i14e4 + (-27269./161740800.)*i14e2 + (
                        96797531./6539968512.)*i12e6 + (-13545047./1459814400.)*i12e4 + (2312569./1094860800.)*i12e2 + (
                        547926499./4954521600.)*i10e8 + (-3547351./24772608.)*i10e6 + (496387./5529600.)*i10e4 + (
                        -84749./4147200.)*i10e2 + (1210797857./2936012800.)*i8e10 + (-88567949./110100480.)*i8e8 + (
                        2867005./2752512.)*i8e6 + (-80237./122880.)*i8e4 + (13699./92160.)*i8e2 + (
                        122622913./131072000.)*i6e12 + (-14230123./6553600.)*i6e10 + (1040911./245760.)*i6e8 + (
                        -33695./6144.)*i6e6 + (6601./1920.)*i6e4 + (-1127./1440.)*i6e2 + (
                        5876487793./5033164800.)*i4e14 + (-175937223./52428800.)*i4e12 + (20417133./2621440.)*i4e10 + (
                        -497827./32768.)*i4e8 + (80575./4096.)*i4e6 + (-3157./256.)*i4e4 + (539./192.)*i4e2 + (
                        3973253733./4697620480.)*i2e16 + (-534226163./209715200.)*i2e14 + (47982879./6553600.)*i2e12 + (
                        -5568309./327680.)*i2e10 + (135771./4096.)*i2e8 + (-21975./512.)*i2e6 + (861./32.)*i2e4 + (
                        -49./8.)*i2e2 + (6258845529./30064771072.)*e18 + (-3973253733./4697620480.)*e16 + (
                        534226163./209715200.)*e14 + (-47982879./6553600.)*e12 + (5568309./327680.)*e10 + (
                        -135771./4096.)*e8 + (21975./512.)*e6 + (-861./32.)*e4 + (49./8.)*e2,
            None,
            3.,
            2.,
            2.
        ),
        (
            '-2*O + 4*n',
            '2*O - 4*n',
            '2, 2, 0, 2',
            -2*spin_freq + 4*orbital_freq,
            (-7880741./7925299200.)*i14e4 + (-129173497./2299207680.)*i12e6 + (95476063./7664025600.)*i12e4 + (
                        -1011551957./1045094400.)*i10e8 + (4733837./8709120.)*i10e6 + (-3498923./29030400.)*i10e4 + (
                        -262647013./38707200.)*i8e10 + (163509307./23224320.)*i8e8 + (-765187./193536.)*i8e6 + (
                        565573./645120.)*i8e4 + (-193699813./8294400.)*i6e12 + (3086807./86400.)*i6e10 + (
                        -1921673./51840.)*i6e8 + (8993./432.)*i6e6 + (-6647./1440.)*i6e4 + (
                        -3813706303./87091200.)*i4e14 + (92639041./1105920.)*i4e12 + (-1476299./11520.)*i4e10 + (
                        919061./6912.)*i4e8 + (-21505./288.)*i4e6 + (3179./192.)*i4e4 + (-173370469./4147200.)*i2e16 + (
                        346700573./3628800.)*i2e14 + (-8421731./46080.)*i2e12 + (134209./480.)*i2e10 + (
                        -83551./288.)*i2e8 + (1955./12.)*i2e6 + (-289./8.)*i2e4 + (-2949969133./182891520.)*e18 + (
                        173370469./4147200.)*e16 + (-346700573./3628800.)*e14 + (8421731./46080.)*e12 + (
                        -134209./480.)*e10 + (83551./288.)*e8 + (-1955./12.)*e6 + (289./8.)*e4,
            None,
            4.,
            2.,
            2.
        ),
        (
            '-2*O + 5*n',
            '2*O - 5*n',
            '2, 2, 0, 3',
            -2*spin_freq + 5*orbital_freq,
            (9435611887./176579149824.)*i12e6 + (13309769915./5350883328.)*i10e8 + (-345788027./668860416.)*i10e6 + (
                        229840011775./6341787648.)*i8e10 + (-10757090825./594542592.)*i8e8 + (
                        279469385./74317824.)*i8e6 + (334286583325./1528823808.)*i6e12 + (
                        -2701236725./14155776.)*i6e10 + (126424675./1327104.)*i6e8 + (-3284515./165888.)*i6e6 + (
                        70723571516125./114152177664.)*i4e14 + (-799380960125./1019215872.)*i4e12 + (
                        6459479125./9437184.)*i4e10 + (-302319875./884736.)*i4e8 + (7854275./110592.)*i4e6 + (
                        21274828753525./25367150592.)*i2e16 + (-6429415592375./4756340736.)*i2e14 + (
                        72670996375./42467328.)*i2e12 + (-587225375./393216.)*i2e10 + (27483625./36864.)*i2e8 + (
                        -714025./4608.)*i2e6 + (2044426346565875./4696546738176.)*e18 + (
                        -21274828753525./25367150592.)*e16 + (6429415592375./4756340736.)*e14 + (
                        -72670996375./42467328.)*e12 + (587225375./393216.)*e10 + (-27483625./36864.)*e8 + (
                        714025./4608.)*e6,
            None,
            5.,
            2.,
            2.
        ),
        (
            '-2*O + 6*n',
            '2*O - 6*n',
            '2, 2, 0, 4',
            -2*spin_freq + 6*orbital_freq,
            (-3439465523./1857945600.)*i10e8 + (-4807560329./68812800.)*i8e10 + (555962173./41287680.)*i8e8 + (
                        -5059005733./6144000.)*i6e12 + (56501731./153600.)*i6e10 + (-6534047./92160.)*i6e8 + (
                        -11204313109./2867200.)*i4e14 + (2419524481./819200.)*i4e12 + (-27022567./20480.)*i4e10 + (
                        3124979./12288.)*i4e8 + (-180544031973./22937600.)*i2e16 + (3055721757./358400.)*i2e14 + (
                        -659870313./102400.)*i2e12 + (7369791./2560.)*i2e10 + (-284089./512.)*i2e8 + (
                        -25998653133./4587520.)*e18 + (180544031973./22937600.)*e16 + (-3055721757./358400.)*e14 + (
                        659870313./102400.)*e12 + (-7369791./2560.)*e10 + (284089./512.)*e8,
            None,
            6.,
            2.,
            2.
        ),
        (
            '-2*O + 7*n',
            '2*O - 7*n',
            '2, 2, 0, 5',
            -2*spin_freq + 7*orbital_freq,
            (14577511952059./339738624000.)*i8e10 + (645169101983./509607936.)*i6e12 + (
                        -1199274105407./5308416000.)*i6e10 + (4590344322529889./407686348800.)*i4e14 + (
                        -1542795678655./339738624.)*i4e12 + (573565876499./707788800.)*i4e10 + (
                        9997999669389767./271790899200.)*i2e16 + (-417304029320899./16986931200.)*i2e14 + (
                        140254152605./14155776.)*i2e12 + (-52142352409./29491200.)*i2e10 + (
                        1518065224694500853./39137889484800.)*e18 + (-9997999669389767./271790899200.)*e16 + (
                        417304029320899./16986931200.)*e14 + (-140254152605./14155776.)*e12 + (
                        52142352409./29491200.)*e10,
            None,
            7.,
            2.,
            2.
        ),
        (
            '-2*O + 8*n',
            '2*O - 8*n',
            '2, 2, 0, 6',
            -2*spin_freq + 8*orbital_freq,
            (-123809233703./186624000.)*i6e12 + (-626057456893./43545600.)*i4e14 + (59213111771./24883200.)*i4e12 + (
                        -578658802849./6773760.)*i2e16 + (56914314263./1814400.)*i2e14 + (
                        -5383010161./1036800.)*i2e12 + (-391340609035087./2743372800.)*e18 + (
                        578658802849./6773760.)*e16 + (-56914314263./1814400.)*e14 + (5383010161./1036800.)*e12,
            None,
            8.,
            2.,
            2.
        ),
        (
            '-2*O + 9*n',
            '2*O - 9*n',
            '2, 2, 0, 7',
            -2*spin_freq + 9*orbital_freq,
            (540772344561273./82208358400.)*i4e14 + (3064748517300717./32883343360.)*i2e16 + (
                        -147483366698529./10276044800.)*i2e14 + (415947859083950607./1503238553600.)*e18 + (
                        -3064748517300717./32883343360.)*e16 + (147483366698529./10276044800.)*e14,
            None,
            9.,
            2.,
            2.
        ),
        (
            '-2*O + 10*n',
            '2*O - 10*n',
            '2, 2, 0, 8',
            -2*spin_freq + 10*orbital_freq,
            (-402063787225./10616832.)*i2e16 + (-176187983600875./668860416.)*e18 + (402063787225./10616832.)*e16,
            None,
            10.,
            2.,
            2.
        ),
        (
            '-2*O + 11*n',
            '2*O - 11*n',
            '2, 2, 0, 9',
            -2*spin_freq + 11*orbital_freq,
            (6648821549377771726369./69039237051187200.)*e18,
            None,
            11.,
            2.,
            2.
        ),
        (
            '-2*O - 7*n',
            '2*O + 7*n',
            '2, 2, 1, -7',
            -2*spin_freq - 7*orbital_freq,
            (186702632281./7549747200.)*i4e14,
            None,
            -7.,
            0.,
            2.
        ),
        (
            '-2*O - 6*n',
            '2*O + 6*n',
            '2, 2, 1, -6',
            -2*spin_freq - 6*orbital_freq,
            (-10029889./1228800.)*i6e12 + (-6926229./716800.)*i4e14 + (10029889./819200.)*i4e12,
            None,
            -6.,
            0.,
            2.
        ),
        (
            '-2*O - 5*n',
            '2*O + 5*n',
            '2, 2, 1, -5',
            -2*spin_freq - 5*orbital_freq,
            (3143529./2621440.)*i8e10 + (982439./1048576.)*i6e12 + (-1047843./262144.)*i6e10 + (
                        13524196093./2113929216.)*i4e14 + (-2947317./2097152.)*i4e12 + (3143529./524288.)*i4e10,
            None,
            -5.,
            0.,
            2.
        ),
        (
            '-2*O - 4*n',
            '2*O + 4*n',
            '2, 2, 1, -4',
            -2*spin_freq - 4*orbital_freq,
            (-14399./138240.)*i10e8 + (9933./51200.)*i8e10 + (5929./10240.)*i8e8 + (-708871./307200.)*i6e12 + (
                        -3311./5120.)*i6e10 + (-5929./3072.)*i6e8 + (3027367./921600.)*i4e14 + (
                        708871./204800.)*i4e12 + (9933./10240.)*i4e10 + (5929./2048.)*i4e8,
            None,
            -4.,
            0.,
            2.
        ),
        (
            '-2*O - 3*n',
            '2*O + 3*n',
            '2, 2, 1, -3',
            -2*spin_freq - 3*orbital_freq,
            (87079./14515200.)*i12e6 + (-118031./2580480.)*i10e8 + (-47753./967680.)*i10e6 + (
                        6019881./13107200.)*i8e10 + (20829./81920.)*i8e8 + (2809./10240.)*i8e6 + (
                        -30355211./15728640.)*i6e12 + (-2006627./1310720.)*i6e10 + (-6943./8192.)*i6e8 + (
                        -2809./3072.)*i6e6 + (3165427449./838860800.)*i4e14 + (30355211./10485760.)*i4e12 + (
                        6019881./2621440.)*i4e10 + (20829./16384.)*i4e8 + (2809./2048.)*i4e6,
            None,
            -3.,
            0.,
            2.
        ),
        (
            '-2*O - 2*n',
            '2*O + 2*n',
            '2, 2, 1, -2',
            -2*spin_freq - 2*orbital_freq,
            (-3./12320.)*i14e4 + (31./7200.)*i12e6 + (31./11200.)*i12e4 + (-28237./483840.)*i10e8 + (
                        -17./480.)*i10e6 + (-51./2240.)*i10e4 + (11757./25600.)*i8e10 + (1661./5120.)*i8e8 + (
                        63./320.)*i8e6 + (81./640.)*i8e4 + (-505601./245760.)*i6e12 + (-3919./2560.)*i6e10 + (
                        -1661./1536.)*i6e8 + (-21./32.)*i6e6 + (-27./64.)*i6e4 + (25565893./6451200.)*i4e14 + (
                        505601./163840.)*i4e12 + (11757./5120.)*i4e10 + (1661./1024.)*i4e8 + (63./64.)*i4e6 + (
                        81./128.)*i4e4,
            None,
            -2.,
            0.,
            2.
        ),
        (
            '-2*O - n',
            '2*O + n',
            '2, 2, 1, -1',
            -2*spin_freq - orbital_freq,
            (5461./756756000.)*i16e2 + (-3./12320.)*i14e4 + (-1./9240.)*i14e2 + (527./107520.)*i12e6 + (
                        31./11200.)*i12e4 + (31./25200.)*i12e2 + (-479587./7741440.)*i10e8 + (-289./7168.)*i10e6 + (
                        -51./2240.)*i10e4 + (-17./1680.)*i10e2 + (3194661./6553600.)*i8e10 + (28211./81920.)*i8e8 + (
                        459./2048.)*i8e6 + (81./640.)*i8e4 + (9./160.)*i8e2 + (-85553819./39321600.)*i6e12 + (
                        -1064887./655360.)*i6e10 + (-28211./24576.)*i6e8 + (-765./1024.)*i6e6 + (-27./64.)*i6e4 + (
                        -3./16.)*i6e2 + (24654653741./5872025600.)*i4e14 + (85553819./26214400.)*i4e12 + (
                        3194661./1310720.)*i4e10 + (28211./16384.)*i4e8 + (2295./2048.)*i4e6 + (81./128.)*i4e4 + (
                        9./32.)*i4e2,
            None,
            -1.,
            0.,
            2.
        ),
        (
            '-2*O',
            '2*O',
            '2, 2, 1, 0',
            -2*spin_freq,
            (-257./1532430900.)*i18 + (5461./567567000.)*i16e2 + (5461./1702701000.)*i16 + (-1./3465.)*i14e4 + (
                        -1./6930.)*i14e2 + (-1./20790.)*i14 + (31./5670.)*i12e6 + (31./9450.)*i12e4 + (
                        31./18900.)*i12e2 + (31./56700.)*i12 + (-17./252.)*i10e8 + (-17./378.)*i10e6 + (
                        -17./630.)*i10e4 + (-17./1260.)*i10e2 + (-17./3780.)*i10 + (21./40.)*i8e10 + (3./8.)*i8e8 + (
                        1./4.)*i8e6 + (3./20.)*i8e4 + (3./40.)*i8e2 + (1./40.)*i8 + (-7./3.)*i6e12 + (-7./4.)*i6e10 + (
                        -5./4.)*i6e8 + (-5./6.)*i6e6 + (-1./2.)*i6e4 + (-1./4.)*i6e2 + (-1./12.)*i6 + (9./2.)*i4e14 + (
                        7./2.)*i4e12 + (21./8.)*i4e10 + (15./8.)*i4e8 + (5./4.)*i4e6 + (3./4.)*i4e4 + (3./8.)*i4e2 + (
                        1./8.)*i4,
            None,
            0.,
            0.,
            2.
        ),
        (
            '-2*O + n',
            '2*O - n',
            '2, 2, 1, 1',
            -2*spin_freq + orbital_freq,
            (5461./756756000.)*i16e2 + (-3./12320.)*i14e4 + (-1./9240.)*i14e2 + (527./107520.)*i12e6 + (
                        31./11200.)*i12e4 + (31./25200.)*i12e2 + (-479587./7741440.)*i10e8 + (-289./7168.)*i10e6 + (
                        -51./2240.)*i10e4 + (-17./1680.)*i10e2 + (3194661./6553600.)*i8e10 + (28211./81920.)*i8e8 + (
                        459./2048.)*i8e6 + (81./640.)*i8e4 + (9./160.)*i8e2 + (-85553819./39321600.)*i6e12 + (
                        -1064887./655360.)*i6e10 + (-28211./24576.)*i6e8 + (-765./1024.)*i6e6 + (-27./64.)*i6e4 + (
                        -3./16.)*i6e2 + (24654653741./5872025600.)*i4e14 + (85553819./26214400.)*i4e12 + (
                        3194661./1310720.)*i4e10 + (28211./16384.)*i4e8 + (2295./2048.)*i4e6 + (81./128.)*i4e4 + (
                        9./32.)*i4e2,
            None,
            1.,
            0.,
            2.
        ),
        (
            '-2*O + 2*n',
            '2*O - 2*n',
            '2, 2, 1, 2',
            -2*spin_freq + 2*orbital_freq,
            (-3./12320.)*i14e4 + (31./7200.)*i12e6 + (31./11200.)*i12e4 + (-28237./483840.)*i10e8 + (
                        -17./480.)*i10e6 + (-51./2240.)*i10e4 + (11757./25600.)*i8e10 + (1661./5120.)*i8e8 + (
                        63./320.)*i8e6 + (81./640.)*i8e4 + (-505601./245760.)*i6e12 + (-3919./2560.)*i6e10 + (
                        -1661./1536.)*i6e8 + (-21./32.)*i6e6 + (-27./64.)*i6e4 + (25565893./6451200.)*i4e14 + (
                        505601./163840.)*i4e12 + (11757./5120.)*i4e10 + (1661./1024.)*i4e8 + (63./64.)*i4e6 + (
                        81./128.)*i4e4,
            None,
            2.,
            0.,
            2.
        ),
        (
            '-2*O + 3*n',
            '2*O - 3*n',
            '2, 2, 1, 3',
            -2*spin_freq + 3*orbital_freq,
            (87079./14515200.)*i12e6 + (-118031./2580480.)*i10e8 + (-47753./967680.)*i10e6 + (
                        6019881./13107200.)*i8e10 + (20829./81920.)*i8e8 + (2809./10240.)*i8e6 + (
                        -30355211./15728640.)*i6e12 + (-2006627./1310720.)*i6e10 + (-6943./8192.)*i6e8 + (
                        -2809./3072.)*i6e6 + (3165427449./838860800.)*i4e14 + (30355211./10485760.)*i4e12 + (
                        6019881./2621440.)*i4e10 + (20829./16384.)*i4e8 + (2809./2048.)*i4e6,
            None,
            3.,
            0.,
            2.
        ),
        (
            '-2*O + 4*n',
            '2*O - 4*n',
            '2, 2, 1, 4',
            -2*spin_freq + 4*orbital_freq,
            (-14399./138240.)*i10e8 + (9933./51200.)*i8e10 + (5929./10240.)*i8e8 + (-708871./307200.)*i6e12 + (
                        -3311./5120.)*i6e10 + (-5929./3072.)*i6e8 + (3027367./921600.)*i4e14 + (
                        708871./204800.)*i4e12 + (9933./10240.)*i4e10 + (5929./2048.)*i4e8,
            None,
            4.,
            0.,
            2.
        ),
        (
            '-2*O + 5*n',
            '2*O - 5*n',
            '2, 2, 1, 5',
            -2*spin_freq + 5*orbital_freq,
            (3143529./2621440.)*i8e10 + (982439./1048576.)*i6e12 + (-1047843./262144.)*i6e10 + (
                        13524196093./2113929216.)*i4e14 + (-2947317./2097152.)*i4e12 + (3143529./524288.)*i4e10,
            None,
            5.,
            0.,
            2.
        ),
        (
            '-2*O + 6*n',
            '2*O - 6*n',
            '2, 2, 1, 6',
            -2*spin_freq + 6*orbital_freq,
            (-10029889./1228800.)*i6e12 + (-6926229./716800.)*i4e14 + (10029889./819200.)*i4e12,
            None,
            6.,
            0.,
            2.
        ),
        (
            '-2*O + 7*n',
            '2*O - 7*n',
            '2, 2, 1, 7',
            -2*spin_freq + 7*orbital_freq,
            (186702632281./7549747200.)*i4e14,
            None,
            7.,
            0.,
            2.
        ),
        (
            '-2*O - 7*n',
            '2*O + 7*n',
            '2, 2, 2, -5',
            -2*spin_freq - 7*orbital_freq,
            (52142352409./7549747200.)*i8e10,
            None,
            -7.,
            -2.,
            2.
        ),
        (
            '-2*O - 6*n',
            '2*O + 6*n',
            '2, 2, 2, -4',
            -2*spin_freq - 6*orbital_freq,
            (-284089./393216.)*i10e8 + (-7369791./655360.)*i8e10 + (284089./131072.)*i8e8,
            None,
            -6.,
            -2.,
            2.
        ),
        (
            '-2*O - 5*n',
            '2*O + 5*n',
            '2, 2, 2, -3',
            -2*spin_freq - 5*orbital_freq,
            (2713295./84934656.)*i12e6 + (27483625./28311552.)*i10e8 + (-714025./3538944.)*i10e6 + (
                        587225375./100663296.)*i8e10 + (-27483625./9437184.)*i8e8 + (714025./1179648.)*i8e6,
            None,
            -5.,
            -2.,
            2.
        ),
        (
            '-2*O - 4*n',
            '2*O + 4*n',
            '2, 2, 2, -2',
            -2*spin_freq - 4*orbital_freq,
            (-289./387072.)*i14e4 + (-7429./221184.)*i12e6 + (5491./737280.)*i12e4 + (-83551./221184.)*i10e8 + (
                        1955./9216.)*i10e6 + (-289./6144.)*i10e4 + (-134209./122880.)*i8e10 + (83551./73728.)*i8e8 + (
                        -1955./3072.)*i8e6 + (289./2048.)*i8e4,
            None,
            -4.,
            -2.,
            2.
        ),
        (
            '-2*O - 3*n',
            '2*O + 3*n',
            '2, 2, 2, -1',
            -2*spin_freq - 3*orbital_freq,
            (3199./353894400.)*i16e2 + (41./73728.)*i14e4 + (-7./55296.)*i14e2 + (27835./3145728.)*i12e6 + (
                        -5453./983040.)*i12e4 + (931./737280.)*i12e2 + (45257./1048576.)*i10e8 + (
                        -7325./131072.)*i10e6 + (287./8192.)*i10e4 + (-49./6144.)*i10e2 + (5568309./83886080.)*i8e10 + (
                        -135771./1048576.)*i8e8 + (21975./131072.)*i8e6 + (-861./8192.)*i8e4 + (49./2048.)*i8e2,
            None,
            -3.,
            -2.,
            2.
        ),
        (
            '-2*O - 2*n',
            '2*O + 2*n',
            '2, 2, 2, 0',
            -2*spin_freq - 2*orbital_freq,
            (-491./12262440960.)*i18 + (-457./123863040.)*i16e2 + (457./619315200.)*i16 + (-1./12288.)*i14e4 + (
                        5./96768.)*i14e2 + (-1./96768.)*i14 + (-361./829440.)*i12e6 + (133./163840.)*i12e4 + (
                        -19./36864.)*i12e2 + (19./184320.)*i12 + (-2065./3538944.)*i10e8 + (19./6912.)*i10e6 + (
                        -21./4096.)*i10e4 + (5./1536.)*i10e2 + (-1./1536.)*i10 + (4079./9830400.)*i8e10 + (
                        2065./1179648.)*i8e8 + (-19./2304.)*i8e6 + (63./4096.)*i8e4 + (-5./512.)*i8e2 + (1./512.)*i8,
            None,
            -2.,
            -2.,
            2.
        ),
        (
            '-2*O - n',
            '2*O + n',
            '2, 2, 2, 1',
            -2*spin_freq - orbital_freq,
            (457./2477260800.)*i16e2 + (1./1548288.)*i14e4 + (-1./387072.)*i14e2 + (247./141557760.)*i12e6 + (
                        -19./2949120.)*i12e4 + (19./737280.)*i12e2 + (-305./28311552.)*i10e8 + (-13./1179648.)*i10e6 + (
                        1./24576.)*i10e4 + (-1./6144.)*i10e2 + (1733./83886080.)*i8e10 + (305./9437184.)*i8e8 + (
                        13./393216.)*i8e6 + (-1./8192.)*i8e4 + (1./2048.)*i8e2,
            None,
            -1.,
            -2.,
            2.
        ),
        (
            '-2*O + n',
            '2*O - n',
            '2, 2, 2, 3',
            -2*spin_freq + orbital_freq,
            (19./424673280.)*i12e6 + (-11./28311552.)*i10e8 + (-1./3538944.)*i10e6 + (619./503316480.)*i8e10 + (
                        11./9437184.)*i8e8 + (1./1179648.)*i8e6,
            None,
            1.,
            -2.,
            2.
        ),
        (
            '-2*O + 2*n',
            '2*O - 2*n',
            '2, 2, 2, 4',
            -2*spin_freq + 2*orbital_freq,
            (-1./884736.)*i10e8 + (7./1474560.)*i8e10 + (1./294912.)*i8e8,
            None,
            2.,
            -2.,
            2.
        ),
        (
            '-2*O + 3*n',
            '2*O - 3*n',
            '2, 2, 2, 5',
            -2*spin_freq + 3*orbital_freq,
            (6561./838860800.)*i8e10,
            None,
            3.,
            -2.,
            2.
        )
    )

    return mode_data_output

@njit
def nsr_modes_t20(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray) \
        -> OutputType:
    """ Non-sync Tidal Modes for truncation level 20, tidal order 2

    Parameters
    ----------
    orbital_freq : np.ndarray
        Planet's orbital mean motion [rad s-1]
    spin_freq : np.ndarray
        Planet's spin frequency [rad s-1]
    eccentricity : np.ndarray
        Planet's eccentricity
    inclination : np.ndarray
        Planet's inclination [rads]

    Returns
    -------
    mode_data_tuple : OutputType
        List of mode datas. Each mode data contains information to calculate the heating and torques.
    """
    # Tidal Potential Builder for l-order = 2, NSR = True.
    # Max Eccentricity Order = 20
    # Max Inclination Order = 20
    # Max q = 11.
    # Number of unique modes = 62.
    # Number of unique frequencies = 52.

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

    mode_data_output = (
        (
            '-6*n',
            '6*n',
            '2, 0, 0, -8',
            -6*orbital_freq,
            (1594323./642252800.)*i4e16,
            None,
            -6.,
            2.,
            0.
        ),
        (
            '-5*n',
            '5*n',
            '2, 0, 0, -7',
            -5*orbital_freq,
            (-244140625./266355081216.)*i6e14 + (2685546875./2841120866304.)*i4e16 + (244140625./177570054144.)*i4e14,
            None,
            -5.,
            2.,
            0.
        ),
        (
            '-4*n',
            '4*n',
            '2, 0, 0, -6',
            -4*orbital_freq,
            (1./6750.)*i8e12 + (-1./2025.)*i6e14 + (-1./2025.)*i6e12 + (23./25200.)*i4e16 + (1./1350.)*i4e14 + (
                        1./1350.)*i4e12,
            None,
            -4.,
            2.,
            0.
        ),
        (
            '-3*n',
            '3*n',
            '2, 0, 0, -5',
            -3*orbital_freq,
            (-12393./917504000.)*i10e10 + (19683./209715200.)*i8e12 + (19683./262144000.)*i8e10 + (
                        -4284333./11744051200.)*i6e14 + (-6561./20971520.)*i6e12 + (-6561./26214400.)*i6e10 + (
                        28284471./53687091200.)*i4e16 + (12852999./23488102400.)*i4e14 + (19683./41943040.)*i4e12 + (
                        19683./52428800.)*i4e10,
            None,
            -3.,
            2.,
            0.
        ),
        (
            '-2*n',
            '2*n',
            '2, 0, 0, -4',
            -2*orbital_freq,
            (31./43545600.)*i12e8 + (-17./2073600.)*i10e10 + (-17./2903040.)*i10e8 + (949./18432000.)*i8e12 + (
                        7./153600.)*i8e10 + (1./30720.)*i8e8 + (-2417./14515200.)*i6e14 + (-949./5529600.)*i6e12 + (
                        -7./46080.)*i6e10 + (-1./9216.)*i6e8 + (22601./99090432.)*i4e16 + (2417./9676800.)*i4e14 + (
                        949./3686400.)*i4e12 + (7./30720.)*i4e10 + (1./6144.)*i4e8,
            None,
            -2.,
            2.,
            0.
        ),
        (
            '-n',
            'n',
            '2, 0, 0, -3',
            -orbital_freq,
            (-1./63866880.)*i14e6 + (341./1393459200.)*i12e8 + (31./174182400.)*i12e6 + (-10523./4954521600.)*i10e10 + (
                        -187./92897280.)*i10e8 + (-17./11612160.)*i10e6 + (62617./5662310400.)*i8e12 + (
                        619./52428800.)*i8e10 + (11./983040.)*i8e8 + (1./122880.)*i8e6 + (
                        -31398887./951268147200.)*i6e14 + (-62617./1698693120.)*i6e12 + (-619./15728640.)*i6e10 + (
                        -11./294912.)*i6e8 + (-1./36864.)*i6e6 + (147400583./3382286745600.)*i4e16 + (
                        31398887./634178764800.)*i4e14 + (62617./1132462080.)*i4e12 + (619./10485760.)*i4e10 + (
                        11./196608.)*i4e8 + (1./24576.)*i4e6,
            None,
            -1.,
            2.,
            0.
        ),
        (
            'n',
            'n',
            '2, 0, 0, -1',
            orbital_freq,
            (-257./8172964800.)*i18e2 + (-5461./36324288000.)*i16e4 + (5461./9081072000.)*i16e2 + (
                        -13./21288960.)*i14e6 + (1./443520.)*i14e4 + (-1./110880.)*i14e2 + (1891./278691840.)*i12e8 + (
                        403./58060800.)*i12e6 + (-31./1209600.)*i12e4 + (31./302400.)*i12e2 + (
                        -29461./825753600.)*i10e10 + (-1037./18579456.)*i10e8 + (-221./3870720.)*i10e6 + (
                        17./80640.)*i10e4 + (-17./20160.)*i10e2 + (277229./1572864000.)*i8e12 + (
                        5199./26214400.)*i8e10 + (61./196608.)*i8e8 + (13./40960.)*i8e6 + (-3./2560.)*i8e4 + (
                        3./640.)*i8e2 + (-165285343./317089382400.)*i6e14 + (-277229./471859200.)*i6e12 + (
                        -1733./2621440.)*i6e10 + (-305./294912.)*i6e8 + (-13./12288.)*i6e6 + (1./256.)*i6e4 + (
                        -1./64.)*i6e2 + (49450862117./71028021657600.)*i4e16 + (165285343./211392921600.)*i4e14 + (
                        277229./314572800.)*i4e12 + (5199./5242880.)*i4e10 + (305./196608.)*i4e8 + (13./8192.)*i4e6 + (
                        -3./512.)*i4e4 + (3./128.)*i4e2,
            None,
            1.,
            2.,
            0.
        ),
        (
            '2*n',
            '2*n',
            '2, 0, 0, 0',
            2*orbital_freq,
            (73./13783770000.)*i20 + (257./408648240.)*i18e2 + (-257./2043241200.)*i18 + (5461./288288000.)*i16e4 + (
                        -5461./454053600.)*i16e2 + (5461./2270268000.)*i16 + (79./498960.)*i14e6 + (-1./3520.)*i14e4 + (
                        1./5544.)*i14e2 + (-1./27720.)*i14 + (114607./174182400.)*i12e8 + (-2449./1360800.)*i12e6 + (
                        31./9600.)*i12e4 + (-31./15120.)*i12e2 + (31./75600.)*i12 + (187697./96768000.)*i10e10 + (
                        -62849./11612160.)*i10e8 + (1343./90720.)*i10e6 + (-17./640.)*i10e4 + (17./1008.)*i10e2 + (
                        -17./5040.)*i10 + (34471./14745600.)*i8e12 + (-11041./1024000.)*i8e10 + (3697./122880.)*i8e8 + (
                        -79./960.)*i8e6 + (189./1280.)*i8e4 + (-3./32.)*i8e2 + (3./160.)*i8 + (
                        1739939./541900800.)*i6e14 + (-34471./4423680.)*i6e12 + (11041./307200.)*i6e10 + (
                        -3697./36864.)*i6e8 + (79./288.)*i6e6 + (-63./128.)*i6e4 + (5./16.)*i6e2 + (-1./16.)*i6 + (
                        -561889./1284505600.)*i4e16 + (-1739939./361267200.)*i4e14 + (34471./2949120.)*i4e12 + (
                        -11041./204800.)*i4e10 + (3697./24576.)*i4e8 + (-79./192.)*i4e6 + (189./256.)*i4e4 + (
                        -15./32.)*i4e2 + (3./32.)*i4,
            None,
            2.,
            2.,
            0.
        ),
        (
            '3*n',
            '3*n',
            '2, 0, 0, 1',
            3*orbital_freq,
            (-257./166795200.)*i18e2 + (-223901./1729728000.)*i16e4 + (5461./185328000.)*i16e2 + (
                        -1465./473088.)*i14e6 + (41./21120.)*i14e4 + (-7./15840.)*i14e2 + (
                        -1402967./51609600.)*i12e8 + (9083./258048.)*i12e6 + (-1271./57600.)*i12e4 + (
                        217./43200.)*i12e2 + (-10517917./91750400.)*i10e10 + (769369./3440640.)*i10e8 + (
                        -24905./86016.)*i10e6 + (697./3840.)*i10e4 + (-119./2880.)*i10e2 + (
                        -143948637./524288000.)*i8e12 + (16704927./26214400.)*i8e10 + (-407313./327680.)*i8e8 + (
                        13185./8192.)*i8e6 + (-2583./2560.)*i8e4 + (147./640.)*i8e2 + (
                        -534226163./1677721600.)*i6e14 + (47982879./52428800.)*i6e12 + (-5568309./2621440.)*i6e10 + (
                        135771./32768.)*i6e8 + (-21975./4096.)*i6e6 + (861./256.)*i6e4 + (-49./64.)*i6e2 + (
                        -11919761199./75161927680.)*i4e16 + (1602678489./3355443200.)*i4e14 + (
                        -143948637./104857600.)*i4e12 + (16704927./5242880.)*i4e10 + (-407313./65536.)*i4e8 + (
                        65925./8192.)*i4e6 + (-2583./512.)*i4e4 + (147./128.)*i4e2,
            None,
            3.,
            2.,
            0.
        ),
        (
            '4*n',
            '4*n',
            '2, 0, 0, 2',
            4*orbital_freq,
            (1578229./9081072000.)*i16e4 + (391./33264.)*i14e6 + (-289./110880.)*i14e4 + (2590081./10886400.)*i12e8 + (
                        -12121./90720.)*i12e6 + (8959./302400.)*i12e4 + (2281553./1209600.)*i10e10 + (
                        -1420367./725760.)*i10e8 + (6647./6048.)*i10e6 + (-4913./20160.)*i10e4 + (
                        8421731./1228800.)*i8e12 + (-134209./12800.)*i8e10 + (83551./7680.)*i8e8 + (-391./64.)*i8e6 + (
                        867./640.)*i8e4 + (346700573./29030400.)*i6e14 + (-8421731./368640.)*i6e12 + (
                        134209./3840.)*i6e10 + (-83551./2304.)*i6e8 + (1955./96.)*i6e6 + (-289./64.)*i6e4 + (
                        173370469./22118400.)*i4e16 + (-346700573./19353600.)*i4e14 + (8421731./245760.)*i4e12 + (
                        -134209./2560.)*i4e10 + (83551./1536.)*i4e8 + (-1955./64.)*i4e6 + (867./128.)*i4e4,
            None,
            4.,
            2.,
            0.
        ),
        (
            '5*n',
            '5*n',
            '2, 0, 0, 3',
            5*orbital_freq,
            (-142805./12773376.)*i14e6 + (-34079695./55738368.)*i12e8 + (885391./6967296.)*i12e6 + (
                        -1996566275./198180864.)*i10e10 + (93444325./18579456.)*i10e8 + (-2427685./2322432.)*i10e6 + (
                        -14534199275./226492416.)*i8e12 + (117445075./2097152.)*i8e10 + (-5496725./196608.)*i8e8 + (
                        142805./24576.)*i8e6 + (-6429415592375./38050725888.)*i6e14 + (
                        72670996375./339738624.)*i6e12 + (-587225375./3145728.)*i6e10 + (27483625./294912.)*i6e8 + (
                        -714025./36864.)*i6e6 + (-21274828753525./135291469824.)*i4e16 + (
                        6429415592375./25367150592.)*i4e14 + (-72670996375./226492416.)*i4e12 + (
                        587225375./2097152.)*i4e10 + (-27483625./196608.)*i4e8 + (714025./24576.)*i4e6,
            None,
            5.,
            2.,
            0.
        ),
        (
            '6*n',
            '6*n',
            '2, 0, 0, 4',
            6*orbital_freq,
            (8806759./19353600.)*i12e8 + (41762149./2150400.)*i10e10 + (-4829513./1290240.)*i10e8 + (
                        1979610939./8192000.)*i8e12 + (-22109373./204800.)*i8e10 + (852267./40960.)*i8e8 + (
                        3055721757./2867200.)*i6e14 + (-659870313./819200.)*i6e12 + (7369791./20480.)*i6e10 + (
                        -284089./4096.)*i6e8 + (541632095919./367001600.)*i4e16 + (-9167165271./5734400.)*i4e14 + (
                        1979610939./1638400.)*i4e12 + (-22109373./40960.)*i4e10 + (852267./8192.)*i4e8,
            None,
            6.,
            2.,
            0.
        ),
        (
            '7*n',
            '7*n',
            '2, 0, 0, 5',
            7*orbital_freq,
            (-126631427279./10616832000.)*i10e10 + (-28050830521./75497472.)*i8e12 + (52142352409./786432000.)*i8e10 + (
                        -417304029320899./135895449600.)*i6e14 + (140254152605./113246208.)*i6e12 + (
                        -52142352409./235929600.)*i6e10 + (-9997999669389767./1449551462400.)*i4e16 + (
                        417304029320899./90596966400.)*i4e14 + (-140254152605./75497472.)*i4e12 + (
                        52142352409./157286400.)*i4e10,
            None,
            7.,
            2.,
            0.
        ),
        (
            '8*n',
            '8*n',
            '2, 0, 0, 6',
            8*orbital_freq,
            (5383010161./27648000.)*i8e12 + (56914314263./14515200.)*i6e14 + (-5383010161./8294400.)*i6e12 + (
                        578658802849./36126720.)*i4e16 + (-56914314263./9676800.)*i4e14 + (5383010161./5529600.)*i4e12,
            None,
            8.,
            2.,
            0.
        ),
        (
            '9*n',
            '9*n',
            '2, 0, 0, 7',
            9*orbital_freq,
            (-147483366698529./82208358400.)*i6e14 + (-9194245551902151./526133493760.)*i4e16 + (
                        442450100095587./164416716800.)*i4e14,
            None,
            9.,
            2.,
            0.
        ),
        (
            '10*n',
            '10*n',
            '2, 0, 0, 8',
            10*orbital_freq,
            (402063787225./56623104.)*i4e16,
            None,
            10.,
            2.,
            0.
        ),
        (
            '-10*n',
            '10*n',
            '2, 0, 1, -10',
            -10*orbital_freq,
            (58301303109841./224737099776.)*e20,
            None,
            -10.,
            0.,
            0.
        ),
        (
            '-9*n',
            '9*n',
            '2, 0, 1, -9',
            -9*orbital_freq,
            (-4143587919679849./10522669875200.)*i2e18 + (-9482058568573459./30064771072000.)*e20 + (
                        4143587919679849./31568009625600.)*e18,
            None,
            -9.,
            0.,
            0.
        ),
        (
            '-8*n',
            '8*n',
            '2, 0, 1, -8',
            -8*orbital_freq,
            (413425292293./1926758400.)*i4e16 + (1606626956771./4335206400.)*i2e18 + (
                        -31801945561./160563200.)*i2e16 + (583590180249631./3511517184000.)*e20 + (
                        -1606626956771./13005619200.)*e18 + (31801945561./481689600.)*e16,
            None,
            -8.,
            0.,
            0.
        ),
        (
            '-7*n',
            '7*n',
            '2, 0, 1, -7',
            -7*orbital_freq,
            (-9148428981769./169869312000.)*i6e14 + (-17266338644269./120795955200.)*i4e16 + (
                        2427134219653./22649241600.)*i4e14 + (-999990833612299./5798205849600.)*i2e18 + (
                        1328179895713./10066329600.)*i2e16 + (-186702632281./1887436800.)*i2e14 + (
                        -105829510207034821./3131031158784000.)*e20 + (999990833612299./17394617548800.)*e18 + (
                        -1328179895713./30198988800.)*e16 + (186702632281./5662310400.)*e14,
            None,
            -7.,
            0.,
            0.
        ),
        (
            '-6*n',
            '6*n',
            '2, 0, 1, -6',
            -6*orbital_freq,
            (1935768577./258048000.)*i8e12 + (5387067./256000.)*i6e14 + (-491464561./18432000.)*i6e12 + (
                        855995205./12845056.)*i4e16 + (-30013659./716800.)*i4e14 + (130388557./2457600.)*i4e12 + (
                        941362189./80281600.)*i2e18 + (-197537355./3211264.)*i2e16 + (6926229./179200.)*i2e14 + (
                        -10029889./204800.)*i2e12 + (168060374177./12845056000.)*e20 + (-941362189./240844800.)*e18 + (
                        65845785./3211264.)*e16 + (-2308743./179200.)*e14 + (10029889./614400.)*e12,
            None,
            -6.,
            0.,
            0.
        ),
        (
            '-5*n',
            '5*n',
            '2, 0, 1, -5',
            -5*orbital_freq,
            (-29844121./45875200.)*i10e10 + (-189610727./220200960.)*i8e12 + (67411233./18350080.)*i8e10 + (
                        -94669372651./6794772480.)*i6e14 + (48139511./15728640.)*i6e12 + (-17114769./1310720.)*i6e10 + (
                        370394165375./33822867456.)*i4e16 + (175814549209./6341787648.)*i4e14 + (
                        -12771707./2097152.)*i4e12 + (13621959./524288.)*i4e10 + (
                        -2447539096445./105226698752.)*i2e18 + (-28491858875./2818572288.)*i2e16 + (
                        -13524196093./528482304.)*i2e14 + (2947317./524288.)*i2e12 + (-3143529./131072.)*i2e10 + (
                        321698310085813./43834436222976.)*e20 + (2447539096445./315680096256.)*e18 + (
                        28491858875./8455716864.)*e16 + (13524196093./1585446912.)*e14 + (-982439./524288.)*e12 + (
                        1047843./131072.)*e10,
            None,
            -5.,
            0.,
            0.
        ),
        (
            '-4*n',
            '4*n',
            '2, 0, 1, -4',
            -4*orbital_freq,
            (236621./6220800.)*i12e8 + (-363737./3456000.)*i10e10 + (-651343./2073600.)*i10e8 + (
                        136812103./64512000.)*i8e12 + (91289./153600.)*i8e10 + (163471./92160.)*i8e8 + (
                        -148340983./20736000.)*i6e14 + (-34734679./4608000.)*i6e12 + (-162239./76800.)*i6e10 + (
                        -290521./46080.)*i6e8 + (1623891763./82575360.)*i4e16 + (39355771./2764800.)*i4e14 + (
                        9215323./614400.)*i4e12 + (43043./10240.)*i4e10 + (77077./6144.)*i4e8 + (
                        -2170376447./103219200.)*i2e18 + (-124914751./6881280.)*i2e16 + (-3027367./230400.)*i2e14 + (
                        -708871./51200.)*i2e12 + (-9933./2560.)*i2e10 + (-5929./512.)*i2e8 + (
                        4872964014199./585252864000.)*e20 + (2170376447./309657600.)*e18 + (
                        124914751./20643840.)*e16 + (3027367./691200.)*e14 + (708871./153600.)*e12 + (
                        3311./2560.)*e10 + (5929./1536.)*e8,
            None,
            -4.,
            0.,
            0.
        ),
        (
            '-3*n',
            '3*n',
            '2, 0, 1, -3',
            -3*orbital_freq,
            (-34519801./21794572800.)*i14e6 + (3047977./182476800.)*i12e8 + (1233151./68428800.)*i12e6 + (
                        -220442309./884736000.)*i10e10 + (-5339167./38707200.)*i10e8 + (-2160121./14515200.)*i10e6 + (
                        5858555723./3303014400.)*i8e12 + (55325573./39321600.)*i8e10 + (1339999./1720320.)*i8e8 + (
                        542137./645120.)*i8e6 + (-17233993889./2097152000.)*i6e14 + (-1487405339./235929600.)*i6e12 + (
                        -98324723./19660800.)*i6e10 + (-340207./122880.)*i6e8 + (-137641./46080.)*i6e6 + (
                        1894718531589./93952409600.)*i4e16 + (13716852279./838860800.)*i4e14 + (
                        394617743./31457280.)*i4e12 + (26086151./2621440.)*i4e10 + (90259./16384.)*i4e8 + (
                        36517./6144.)*i4e6 + (-8454941691423./375809638400.)*i2e18 + (
                        -437242738059./23488102400.)*i2e16 + (-3165427449./209715200.)*i2e14 + (
                        -30355211./2621440.)*i2e12 + (-6019881./655360.)*i2e10 + (-20829./4096.)*i2e8 + (
                        -2809./512.)*i2e6 + (466859091912363./52613349376000.)*e20 + (
                        2818313897141./375809638400.)*e18 + (145747579353./23488102400.)*e16 + (
                        1055142483./209715200.)*e14 + (30355211./7864320.)*e12 + (2006627./655360.)*e10 + (
                        6943./4096.)*e8 + (2809./1536.)*e6,
            None,
            -3.,
            0.,
            0.
        ),
        (
            '-2*n',
            '2*n',
            '2, 0, 1, -2',
            -2*orbital_freq,
            (3781./77616000.)*i16e4 + (-12289./10810800.)*i14e6 + (-12289./16816800.)*i14e4 + (
                        66289./3110400.)*i12e8 + (3073./237600.)*i12e6 + (439./52800.)*i12e4 + (
                        -3013711./12096000.)*i10e10 + (-1277309./7257600.)*i10e8 + (-769./7200.)*i10e6 + (
                        -769./11200.)*i10e4 + (97580993./51609600.)*i8e12 + (756367./537600.)*i8e10 + (
                        320573./322560.)*i8e8 + (193./320.)*i8e6 + (1737./4480.)*i8e4 + (
                        -178961251./20736000.)*i6e14 + (-24774449./3686400.)*i6e12 + (-192031./38400.)*i6e10 + (
                        -81389./23040.)*i6e8 + (-343./160.)*i6e6 + (-441./320.)*i6e4 + (
                        8821079033./412876800.)*i4e16 + (332356609./19353600.)*i4e14 + (6572813./491520.)*i4e12 + (
                        50947./5120.)*i4e10 + (21593./3072.)*i4e8 + (273./64.)*i4e6 + (351./128.)*i4e4 + (
                        -51883919761./2167603200.)*i2e18 + (-678544541./34406400.)*i2e16 + (
                        -25565893./1612800.)*i2e14 + (-505601./40960.)*i2e12 + (-11757./1280.)*i2e10 + (
                        -1661./256.)*i2e8 + (-63./16.)*i2e6 + (-81./32.)*i2e4 + (6351400670087./668860416000.)*e20 + (
                        51883919761./6502809600.)*e18 + (678544541./103219200.)*e16 + (25565893./4838400.)*e14 + (
                        505601./122880.)*e12 + (3919./1280.)*e10 + (1661./768.)*e8 + (21./16.)*e6 + (27./32.)*e4,
            None,
            -2.,
            0.,
            0.
        ),
        (
            '-n',
            'n',
            '2, 0, 1, -1',
            -orbital_freq,
            (-28087./24810786000.)*i18e2 + (3781./77616000.)*i16e4 + (3781./174636000.)*i16e2 + (
                        -208913./161441280.)*i14e6 + (-12289./16816800.)*i14e4 + (-12289./37837800.)*i14e2 + (
                        12384629./547430400.)*i12e8 + (7463./506880.)*i12e6 + (439./52800.)*i12e4 + (
                        439./118800.)*i12e2 + (-818898103./3096576000.)*i10e10 + (-21694259./116121600.)*i10e8 + (
                        -13073./107520.)*i10e6 + (-769./11200.)*i10e4 + (-769./25200.)*i10e2 + (
                        16511887067./8257536000.)*i8e12 + (205523191./137625600.)*i8e10 + (5444723./5160960.)*i8e8 + (
                        9843./14336.)*i8e6 + (1737./4480.)*i8e4 + (193./1120.)*i8e2 + (
                        -172582576187./18874368000.)*i6e14 + (-4192137131./589824000.)*i6e12 + (
                        -52179463./9830400.)*i6e10 + (-1382339./368640.)*i6e8 + (-2499./1024.)*i6e6 + (
                        -441./320.)*i6e4 + (-49./80.)*i6e2 + (8961067149151./394600120320.)*i4e16 + (
                        320510498633./17616076800.)*i4e14 + (1112199647./78643200.)*i4e12 + (
                        13843531./1310720.)*i4e10 + (366743./49152.)*i4e8 + (9945./2048.)*i4e6 + (351./128.)*i4e4 + (
                        39./32.)*i4e2 + (-725941889609009./28411208663040.)*i2e18 + (
                        -689312857627./32883343360.)*i2e16 + (-24654653741./1468006400.)*i2e14 + (
                        -85553819./6553600.)*i2e12 + (-3194661./327680.)*i2e10 + (-28211./4096.)*i2e8 + (
                        -2295./512.)*i2e6 + (-81./32.)*i2e4 + (-9./8.)*i2e2 + (
                        2343932511383522191./230130790170624000.)*e20 + (725941889609009./85233625989120.)*e18 + (
                        689312857627./98650030080.)*e16 + (24654653741./4404019200.)*e14 + (85553819./19660800.)*e12 + (
                        1064887./327680.)*e10 + (28211./12288.)*e8 + (765./512.)*e6 + (27./32.)*e4 + (3./8.)*e2,
            None,
            -1.,
            0.,
            0.
        ),
        (
            'n',
            'n',
            '2, 0, 1, 1',
            orbital_freq,
            (-28087./24810786000.)*i18e2 + (3781./77616000.)*i16e4 + (3781./174636000.)*i16e2 + (
                        -208913./161441280.)*i14e6 + (-12289./16816800.)*i14e4 + (-12289./37837800.)*i14e2 + (
                        12384629./547430400.)*i12e8 + (7463./506880.)*i12e6 + (439./52800.)*i12e4 + (
                        439./118800.)*i12e2 + (-818898103./3096576000.)*i10e10 + (-21694259./116121600.)*i10e8 + (
                        -13073./107520.)*i10e6 + (-769./11200.)*i10e4 + (-769./25200.)*i10e2 + (
                        16511887067./8257536000.)*i8e12 + (205523191./137625600.)*i8e10 + (5444723./5160960.)*i8e8 + (
                        9843./14336.)*i8e6 + (1737./4480.)*i8e4 + (193./1120.)*i8e2 + (
                        -172582576187./18874368000.)*i6e14 + (-4192137131./589824000.)*i6e12 + (
                        -52179463./9830400.)*i6e10 + (-1382339./368640.)*i6e8 + (-2499./1024.)*i6e6 + (
                        -441./320.)*i6e4 + (-49./80.)*i6e2 + (8961067149151./394600120320.)*i4e16 + (
                        320510498633./17616076800.)*i4e14 + (1112199647./78643200.)*i4e12 + (
                        13843531./1310720.)*i4e10 + (366743./49152.)*i4e8 + (9945./2048.)*i4e6 + (351./128.)*i4e4 + (
                        39./32.)*i4e2 + (-725941889609009./28411208663040.)*i2e18 + (
                        -689312857627./32883343360.)*i2e16 + (-24654653741./1468006400.)*i2e14 + (
                        -85553819./6553600.)*i2e12 + (-3194661./327680.)*i2e10 + (-28211./4096.)*i2e8 + (
                        -2295./512.)*i2e6 + (-81./32.)*i2e4 + (-9./8.)*i2e2 + (
                        2343932511383522191./230130790170624000.)*e20 + (725941889609009./85233625989120.)*e18 + (
                        689312857627./98650030080.)*e16 + (24654653741./4404019200.)*e14 + (85553819./19660800.)*e12 + (
                        1064887./327680.)*e10 + (28211./12288.)*e8 + (765./512.)*e6 + (27./32.)*e4 + (3./8.)*e2,
            None,
            1.,
            0.,
            0.
        ),
        (
            '2*n',
            '2*n',
            '2, 0, 1, 2',
            2*orbital_freq,
            (3781./77616000.)*i16e4 + (-12289./10810800.)*i14e6 + (-12289./16816800.)*i14e4 + (
                        66289./3110400.)*i12e8 + (3073./237600.)*i12e6 + (439./52800.)*i12e4 + (
                        -3013711./12096000.)*i10e10 + (-1277309./7257600.)*i10e8 + (-769./7200.)*i10e6 + (
                        -769./11200.)*i10e4 + (97580993./51609600.)*i8e12 + (756367./537600.)*i8e10 + (
                        320573./322560.)*i8e8 + (193./320.)*i8e6 + (1737./4480.)*i8e4 + (
                        -178961251./20736000.)*i6e14 + (-24774449./3686400.)*i6e12 + (-192031./38400.)*i6e10 + (
                        -81389./23040.)*i6e8 + (-343./160.)*i6e6 + (-441./320.)*i6e4 + (
                        8821079033./412876800.)*i4e16 + (332356609./19353600.)*i4e14 + (6572813./491520.)*i4e12 + (
                        50947./5120.)*i4e10 + (21593./3072.)*i4e8 + (273./64.)*i4e6 + (351./128.)*i4e4 + (
                        -51883919761./2167603200.)*i2e18 + (-678544541./34406400.)*i2e16 + (
                        -25565893./1612800.)*i2e14 + (-505601./40960.)*i2e12 + (-11757./1280.)*i2e10 + (
                        -1661./256.)*i2e8 + (-63./16.)*i2e6 + (-81./32.)*i2e4 + (6351400670087./668860416000.)*e20 + (
                        51883919761./6502809600.)*e18 + (678544541./103219200.)*e16 + (25565893./4838400.)*e14 + (
                        505601./122880.)*e12 + (3919./1280.)*e10 + (1661./768.)*e8 + (21./16.)*e6 + (27./32.)*e4,
            None,
            2.,
            0.,
            0.
        ),
        (
            '3*n',
            '3*n',
            '2, 0, 1, 3',
            3*orbital_freq,
            (-34519801./21794572800.)*i14e6 + (3047977./182476800.)*i12e8 + (1233151./68428800.)*i12e6 + (
                        -220442309./884736000.)*i10e10 + (-5339167./38707200.)*i10e8 + (-2160121./14515200.)*i10e6 + (
                        5858555723./3303014400.)*i8e12 + (55325573./39321600.)*i8e10 + (1339999./1720320.)*i8e8 + (
                        542137./645120.)*i8e6 + (-17233993889./2097152000.)*i6e14 + (-1487405339./235929600.)*i6e12 + (
                        -98324723./19660800.)*i6e10 + (-340207./122880.)*i6e8 + (-137641./46080.)*i6e6 + (
                        1894718531589./93952409600.)*i4e16 + (13716852279./838860800.)*i4e14 + (
                        394617743./31457280.)*i4e12 + (26086151./2621440.)*i4e10 + (90259./16384.)*i4e8 + (
                        36517./6144.)*i4e6 + (-8454941691423./375809638400.)*i2e18 + (
                        -437242738059./23488102400.)*i2e16 + (-3165427449./209715200.)*i2e14 + (
                        -30355211./2621440.)*i2e12 + (-6019881./655360.)*i2e10 + (-20829./4096.)*i2e8 + (
                        -2809./512.)*i2e6 + (466859091912363./52613349376000.)*e20 + (
                        2818313897141./375809638400.)*e18 + (145747579353./23488102400.)*e16 + (
                        1055142483./209715200.)*e14 + (30355211./7864320.)*e12 + (2006627./655360.)*e10 + (
                        6943./4096.)*e8 + (2809./1536.)*e6,
            None,
            3.,
            0.,
            0.
        ),
        (
            '4*n',
            '4*n',
            '2, 0, 1, 4',
            4*orbital_freq,
            (236621./6220800.)*i12e8 + (-363737./3456000.)*i10e10 + (-651343./2073600.)*i10e8 + (
                        136812103./64512000.)*i8e12 + (91289./153600.)*i8e10 + (163471./92160.)*i8e8 + (
                        -148340983./20736000.)*i6e14 + (-34734679./4608000.)*i6e12 + (-162239./76800.)*i6e10 + (
                        -290521./46080.)*i6e8 + (1623891763./82575360.)*i4e16 + (39355771./2764800.)*i4e14 + (
                        9215323./614400.)*i4e12 + (43043./10240.)*i4e10 + (77077./6144.)*i4e8 + (
                        -2170376447./103219200.)*i2e18 + (-124914751./6881280.)*i2e16 + (-3027367./230400.)*i2e14 + (
                        -708871./51200.)*i2e12 + (-9933./2560.)*i2e10 + (-5929./512.)*i2e8 + (
                        4872964014199./585252864000.)*e20 + (2170376447./309657600.)*e18 + (
                        124914751./20643840.)*e16 + (3027367./691200.)*e14 + (708871./153600.)*e12 + (
                        3311./2560.)*e10 + (5929./1536.)*e8,
            None,
            4.,
            0.,
            0.
        ),
        (
            '5*n',
            '5*n',
            '2, 0, 1, 5',
            5*orbital_freq,
            (-29844121./45875200.)*i10e10 + (-189610727./220200960.)*i8e12 + (67411233./18350080.)*i8e10 + (
                        -94669372651./6794772480.)*i6e14 + (48139511./15728640.)*i6e12 + (-17114769./1310720.)*i6e10 + (
                        370394165375./33822867456.)*i4e16 + (175814549209./6341787648.)*i4e14 + (
                        -12771707./2097152.)*i4e12 + (13621959./524288.)*i4e10 + (
                        -2447539096445./105226698752.)*i2e18 + (-28491858875./2818572288.)*i2e16 + (
                        -13524196093./528482304.)*i2e14 + (2947317./524288.)*i2e12 + (-3143529./131072.)*i2e10 + (
                        321698310085813./43834436222976.)*e20 + (2447539096445./315680096256.)*e18 + (
                        28491858875./8455716864.)*e16 + (13524196093./1585446912.)*e14 + (-982439./524288.)*e12 + (
                        1047843./131072.)*e10,
            None,
            5.,
            0.,
            0.
        ),
        (
            '6*n',
            '6*n',
            '2, 0, 1, 6',
            6*orbital_freq,
            (1935768577./258048000.)*i8e12 + (5387067./256000.)*i6e14 + (-491464561./18432000.)*i6e12 + (
                        855995205./12845056.)*i4e16 + (-30013659./716800.)*i4e14 + (130388557./2457600.)*i4e12 + (
                        941362189./80281600.)*i2e18 + (-197537355./3211264.)*i2e16 + (6926229./179200.)*i2e14 + (
                        -10029889./204800.)*i2e12 + (168060374177./12845056000.)*e20 + (-941362189./240844800.)*e18 + (
                        65845785./3211264.)*e16 + (-2308743./179200.)*e14 + (10029889./614400.)*e12,
            None,
            6.,
            0.,
            0.
        ),
        (
            '7*n',
            '7*n',
            '2, 0, 1, 7',
            7*orbital_freq,
            (-9148428981769./169869312000.)*i6e14 + (-17266338644269./120795955200.)*i4e16 + (
                        2427134219653./22649241600.)*i4e14 + (-999990833612299./5798205849600.)*i2e18 + (
                        1328179895713./10066329600.)*i2e16 + (-186702632281./1887436800.)*i2e14 + (
                        -105829510207034821./3131031158784000.)*e20 + (999990833612299./17394617548800.)*e18 + (
                        -1328179895713./30198988800.)*e16 + (186702632281./5662310400.)*e14,
            None,
            7.,
            0.,
            0.
        ),
        (
            '8*n',
            '8*n',
            '2, 0, 1, 8',
            8*orbital_freq,
            (413425292293./1926758400.)*i4e16 + (1606626956771./4335206400.)*i2e18 + (
                        -31801945561./160563200.)*i2e16 + (583590180249631./3511517184000.)*e20 + (
                        -1606626956771./13005619200.)*e18 + (31801945561./481689600.)*e16,
            None,
            8.,
            0.,
            0.
        ),
        (
            '9*n',
            '9*n',
            '2, 0, 1, 9',
            9*orbital_freq,
            (-4143587919679849./10522669875200.)*i2e18 + (-9482058568573459./30064771072000.)*e20 + (
                        4143587919679849./31568009625600.)*e18,
            None,
            9.,
            0.,
            0.
        ),
        (
            '10*n',
            '10*n',
            '2, 0, 1, 10',
            10*orbital_freq,
            (58301303109841./224737099776.)*e20,
            None,
            10.,
            0.,
            0.
        ),
        (
            '-10*n',
            '10*n',
            '2, 0, 2, -8',
            -10*orbital_freq,
            (402063787225./56623104.)*i4e16,
            None,
            -10.,
            -2.,
            0.
        ),
        (
            '-9*n',
            '9*n',
            '2, 0, 2, -7',
            -9*orbital_freq,
            (-147483366698529./82208358400.)*i6e14 + (-9194245551902151./526133493760.)*i4e16 + (
                        442450100095587./164416716800.)*i4e14,
            None,
            -9.,
            -2.,
            0.
        ),
        (
            '-8*n',
            '8*n',
            '2, 0, 2, -6',
            -8*orbital_freq,
            (5383010161./27648000.)*i8e12 + (56914314263./14515200.)*i6e14 + (-5383010161./8294400.)*i6e12 + (
                        578658802849./36126720.)*i4e16 + (-56914314263./9676800.)*i4e14 + (5383010161./5529600.)*i4e12,
            None,
            -8.,
            -2.,
            0.
        ),
        (
            '-7*n',
            '7*n',
            '2, 0, 2, -5',
            -7*orbital_freq,
            (-126631427279./10616832000.)*i10e10 + (-28050830521./75497472.)*i8e12 + (52142352409./786432000.)*i8e10 + (
                        -417304029320899./135895449600.)*i6e14 + (140254152605./113246208.)*i6e12 + (
                        -52142352409./235929600.)*i6e10 + (-9997999669389767./1449551462400.)*i4e16 + (
                        417304029320899./90596966400.)*i4e14 + (-140254152605./75497472.)*i4e12 + (
                        52142352409./157286400.)*i4e10,
            None,
            -7.,
            -2.,
            0.
        ),
        (
            '-6*n',
            '6*n',
            '2, 0, 2, -4',
            -6*orbital_freq,
            (8806759./19353600.)*i12e8 + (41762149./2150400.)*i10e10 + (-4829513./1290240.)*i10e8 + (
                        1979610939./8192000.)*i8e12 + (-22109373./204800.)*i8e10 + (852267./40960.)*i8e8 + (
                        3055721757./2867200.)*i6e14 + (-659870313./819200.)*i6e12 + (7369791./20480.)*i6e10 + (
                        -284089./4096.)*i6e8 + (541632095919./367001600.)*i4e16 + (-9167165271./5734400.)*i4e14 + (
                        1979610939./1638400.)*i4e12 + (-22109373./40960.)*i4e10 + (852267./8192.)*i4e8,
            None,
            -6.,
            -2.,
            0.
        ),
        (
            '-5*n',
            '5*n',
            '2, 0, 2, -3',
            -5*orbital_freq,
            (-142805./12773376.)*i14e6 + (-34079695./55738368.)*i12e8 + (885391./6967296.)*i12e6 + (
                        -1996566275./198180864.)*i10e10 + (93444325./18579456.)*i10e8 + (-2427685./2322432.)*i10e6 + (
                        -14534199275./226492416.)*i8e12 + (117445075./2097152.)*i8e10 + (-5496725./196608.)*i8e8 + (
                        142805./24576.)*i8e6 + (-6429415592375./38050725888.)*i6e14 + (
                        72670996375./339738624.)*i6e12 + (-587225375./3145728.)*i6e10 + (27483625./294912.)*i6e8 + (
                        -714025./36864.)*i6e6 + (-21274828753525./135291469824.)*i4e16 + (
                        6429415592375./25367150592.)*i4e14 + (-72670996375./226492416.)*i4e12 + (
                        587225375./2097152.)*i4e10 + (-27483625./196608.)*i4e8 + (714025./24576.)*i4e6,
            None,
            -5.,
            -2.,
            0.
        ),
        (
            '-4*n',
            '4*n',
            '2, 0, 2, -2',
            -4*orbital_freq,
            (1578229./9081072000.)*i16e4 + (391./33264.)*i14e6 + (-289./110880.)*i14e4 + (2590081./10886400.)*i12e8 + (
                        -12121./90720.)*i12e6 + (8959./302400.)*i12e4 + (2281553./1209600.)*i10e10 + (
                        -1420367./725760.)*i10e8 + (6647./6048.)*i10e6 + (-4913./20160.)*i10e4 + (
                        8421731./1228800.)*i8e12 + (-134209./12800.)*i8e10 + (83551./7680.)*i8e8 + (-391./64.)*i8e6 + (
                        867./640.)*i8e4 + (346700573./29030400.)*i6e14 + (-8421731./368640.)*i6e12 + (
                        134209./3840.)*i6e10 + (-83551./2304.)*i6e8 + (1955./96.)*i6e6 + (-289./64.)*i6e4 + (
                        173370469./22118400.)*i4e16 + (-346700573./19353600.)*i4e14 + (8421731./245760.)*i4e12 + (
                        -134209./2560.)*i4e10 + (83551./1536.)*i4e8 + (-1955./64.)*i4e6 + (867./128.)*i4e4,
            None,
            -4.,
            -2.,
            0.
        ),
        (
            '-3*n',
            '3*n',
            '2, 0, 2, -1',
            -3*orbital_freq,
            (-257./166795200.)*i18e2 + (-223901./1729728000.)*i16e4 + (5461./185328000.)*i16e2 + (
                        -1465./473088.)*i14e6 + (41./21120.)*i14e4 + (-7./15840.)*i14e2 + (
                        -1402967./51609600.)*i12e8 + (9083./258048.)*i12e6 + (-1271./57600.)*i12e4 + (
                        217./43200.)*i12e2 + (-10517917./91750400.)*i10e10 + (769369./3440640.)*i10e8 + (
                        -24905./86016.)*i10e6 + (697./3840.)*i10e4 + (-119./2880.)*i10e2 + (
                        -143948637./524288000.)*i8e12 + (16704927./26214400.)*i8e10 + (-407313./327680.)*i8e8 + (
                        13185./8192.)*i8e6 + (-2583./2560.)*i8e4 + (147./640.)*i8e2 + (
                        -534226163./1677721600.)*i6e14 + (47982879./52428800.)*i6e12 + (-5568309./2621440.)*i6e10 + (
                        135771./32768.)*i6e8 + (-21975./4096.)*i6e6 + (861./256.)*i6e4 + (-49./64.)*i6e2 + (
                        -11919761199./75161927680.)*i4e16 + (1602678489./3355443200.)*i4e14 + (
                        -143948637./104857600.)*i4e12 + (16704927./5242880.)*i4e10 + (-407313./65536.)*i4e8 + (
                        65925./8192.)*i4e6 + (-2583./512.)*i4e4 + (147./128.)*i4e2,
            None,
            -3.,
            -2.,
            0.
        ),
        (
            '-2*n',
            '2*n',
            '2, 0, 2, 0',
            -2*orbital_freq,
            (73./13783770000.)*i20 + (257./408648240.)*i18e2 + (-257./2043241200.)*i18 + (5461./288288000.)*i16e4 + (
                        -5461./454053600.)*i16e2 + (5461./2270268000.)*i16 + (19./124740.)*i14e6 + (-1./3520.)*i14e4 + (
                        1./5544.)*i14e2 + (-1./27720.)*i14 + (1829./4976640.)*i12e8 + (-589./340200.)*i12e6 + (
                        31./9600.)*i12e4 + (-31./15120.)*i12e2 + (31./75600.)*i12 + (-69343./96768000.)*i10e10 + (
                        -1003./331776.)*i10e8 + (323./22680.)*i10e6 + (-17./640.)*i10e4 + (17./1008.)*i10e2 + (
                        -17./5040.)*i10 + (-2323./4915200.)*i8e12 + (4079./1024000.)*i8e10 + (413./24576.)*i8e8 + (
                        -19./240.)*i8e6 + (189./1280.)*i8e4 + (-3./32.)*i8e2 + (3./160.)*i8 + (
                        -689797./180633600.)*i6e14 + (2323./1474560.)*i6e12 + (-4079./307200.)*i6e10 + (
                        -2065./36864.)*i6e8 + (19./72.)*i6e6 + (-63./128.)*i6e4 + (5./16.)*i6e2 + (-1./16.)*i6 + (
                        18956717./3853516800.)*i4e16 + (689797./120422400.)*i4e14 + (-2323./983040.)*i4e12 + (
                        4079./204800.)*i4e10 + (2065./24576.)*i4e8 + (-19./48.)*i4e6 + (189./256.)*i4e4 + (
                        -15./32.)*i4e2 + (3./32.)*i4,
            None,
            -2.,
            -2.,
            0.
        ),
        (
            '-n',
            'n',
            '2, 0, 2, 1',
            -orbital_freq,
            (-257./8172964800.)*i18e2 + (-5461./36324288000.)*i16e4 + (5461./9081072000.)*i16e2 + (
                        -13./21288960.)*i14e6 + (1./443520.)*i14e4 + (-1./110880.)*i14e2 + (1891./278691840.)*i12e8 + (
                        403./58060800.)*i12e6 + (-31./1209600.)*i12e4 + (31./302400.)*i12e2 + (
                        -29461./825753600.)*i10e10 + (-1037./18579456.)*i10e8 + (-221./3870720.)*i10e6 + (
                        17./80640.)*i10e4 + (-17./20160.)*i10e2 + (277229./1572864000.)*i8e12 + (
                        5199./26214400.)*i8e10 + (61./196608.)*i8e8 + (13./40960.)*i8e6 + (-3./2560.)*i8e4 + (
                        3./640.)*i8e2 + (-165285343./317089382400.)*i6e14 + (-277229./471859200.)*i6e12 + (
                        -1733./2621440.)*i6e10 + (-305./294912.)*i6e8 + (-13./12288.)*i6e6 + (1./256.)*i6e4 + (
                        -1./64.)*i6e2 + (49450862117./71028021657600.)*i4e16 + (165285343./211392921600.)*i4e14 + (
                        277229./314572800.)*i4e12 + (5199./5242880.)*i4e10 + (305./196608.)*i4e8 + (13./8192.)*i4e6 + (
                        -3./512.)*i4e4 + (3./128.)*i4e2,
            None,
            -1.,
            -2.,
            0.
        ),
        (
            'n',
            'n',
            '2, 0, 2, 3',
            orbital_freq,
            (-1./63866880.)*i14e6 + (341./1393459200.)*i12e8 + (31./174182400.)*i12e6 + (-10523./4954521600.)*i10e10 + (
                        -187./92897280.)*i10e8 + (-17./11612160.)*i10e6 + (62617./5662310400.)*i8e12 + (
                        619./52428800.)*i8e10 + (11./983040.)*i8e8 + (1./122880.)*i8e6 + (
                        -31398887./951268147200.)*i6e14 + (-62617./1698693120.)*i6e12 + (-619./15728640.)*i6e10 + (
                        -11./294912.)*i6e8 + (-1./36864.)*i6e6 + (147400583./3382286745600.)*i4e16 + (
                        31398887./634178764800.)*i4e14 + (62617./1132462080.)*i4e12 + (619./10485760.)*i4e10 + (
                        11./196608.)*i4e8 + (1./24576.)*i4e6,
            None,
            1.,
            -2.,
            0.
        ),
        (
            '2*n',
            '2*n',
            '2, 0, 2, 4',
            2*orbital_freq,
            (31./43545600.)*i12e8 + (-17./2073600.)*i10e10 + (-17./2903040.)*i10e8 + (949./18432000.)*i8e12 + (
                        7./153600.)*i8e10 + (1./30720.)*i8e8 + (-2417./14515200.)*i6e14 + (-949./5529600.)*i6e12 + (
                        -7./46080.)*i6e10 + (-1./9216.)*i6e8 + (22601./99090432.)*i4e16 + (2417./9676800.)*i4e14 + (
                        949./3686400.)*i4e12 + (7./30720.)*i4e10 + (1./6144.)*i4e8,
            None,
            2.,
            -2.,
            0.
        ),
        (
            '3*n',
            '3*n',
            '2, 0, 2, 5',
            3*orbital_freq,
            (-12393./917504000.)*i10e10 + (19683./209715200.)*i8e12 + (19683./262144000.)*i8e10 + (
                        -4284333./11744051200.)*i6e14 + (-6561./20971520.)*i6e12 + (-6561./26214400.)*i6e10 + (
                        28284471./53687091200.)*i4e16 + (12852999./23488102400.)*i4e14 + (19683./41943040.)*i4e12 + (
                        19683./52428800.)*i4e10,
            None,
            3.,
            -2.,
            0.
        ),
        (
            '4*n',
            '4*n',
            '2, 0, 2, 6',
            4*orbital_freq,
            (1./6750.)*i8e12 + (-1./2025.)*i6e14 + (-1./2025.)*i6e12 + (23./25200.)*i4e16 + (1./1350.)*i4e14 + (
                        1./1350.)*i4e12,
            None,
            4.,
            -2.,
            0.
        ),
        (
            '5*n',
            '5*n',
            '2, 0, 2, 7',
            5*orbital_freq,
            (-244140625./266355081216.)*i6e14 + (2685546875./2841120866304.)*i4e16 + (244140625./177570054144.)*i4e14,
            None,
            5.,
            -2.,
            0.
        ),
        (
            '6*n',
            '6*n',
            '2, 0, 2, 8',
            6*orbital_freq,
            (1594323./642252800.)*i4e16,
            None,
            6.,
            -2.,
            0.
        ),
        (
            '-O - 7*n',
            'O + 7*n',
            '2, 1, 0, -9',
            -spin_freq - 7*orbital_freq,
            (33232930569601./1408964021452800.)*i2e18,
            None,
            -7.,
            2.,
            1.
        ),
        (
            '-O - 6*n',
            'O + 6*n',
            '2, 1, 0, -8',
            -spin_freq - 6*orbital_freq,
            (-177147./16056320.)*i4e16 + (177147./40140800.)*i2e18 + (531441./40140800.)*i2e16,
            None,
            -6.,
            2.,
            1.
        ),
        (
            '-O - 5*n',
            'O + 5*n',
            '2, 1, 0, -7',
            -spin_freq - 5*orbital_freq,
            (11083984375./4794391461888.)*i6e14 + (-13427734375./3196260974592.)*i4e16 + (
                        -1220703125./199766310912.)*i4e14 + (323974609375./43834436222976.)*i2e18 + (
                        2685546875./532710162432.)*i2e16 + (244140625./33294385152.)*i2e14,
            None,
            -5.,
            2.,
            1.
        ),
        (
            '-O - 4*n',
            'O + 4*n',
            '2, 1, 0, -6',
            -spin_freq - 4*orbital_freq,
            (-29./102060.)*i8e12 + (227./182250.)*i6e14 + (227./182250.)*i6e12 + (-23./5670.)*i4e16 + (
                        -4./1215.)*i4e14 + (-4./1215.)*i4e12 + (68./15309.)*i2e18 + (23./4725.)*i2e16 + (
                        8./2025.)*i2e14 + (8./2025.)*i2e12,
            None,
            -4.,
            2.,
            1.
        ),
        (
            '-O - 3*n',
            'O + 3*n',
            '2, 1, 0, -5',
            -spin_freq - 3*orbital_freq,
            (3262437./146800640000.)*i10e10 + (-21141./117440512.)*i8e12 + (-21141./146800640.)*i8e10 + (
                        108060399./117440512000.)*i6e14 + (165483./209715200.)*i6e12 + (165483./262144000.)*i6e10 + (
                        -3142719./1342177280.)*i4e16 + (-1428111./587202560.)*i4e14 + (-2187./1048576.)*i4e12 + (
                        -2187./1310720.)*i4e10 + (270241029./105226698752.)*i2e18 + (9428157./3355443200.)*i2e16 + (
                        4284333./1468006400.)*i2e14 + (6561./2621440.)*i2e12 + (6561./3276800.)*i2e10,
            None,
            -3.,
            2.,
            1.
        ),
        (
            '-O - 2*n',
            'O + 2*n',
            '2, 1, 0, -4',
            -spin_freq - 2*orbital_freq,
            (-59123./55180984320.)*i12e8 + (40277./2985984000.)*i10e10 + (40277./4180377600.)*i10e8 + (
                        -27521./278691840.)*i8e12 + (-29./331776.)*i8e10 + (-145./2322432.)*i8e8 + (
                        548659./1306368000.)*i6e14 + (215423./497664000.)*i6e12 + (1589./4147200.)*i6e10 + (
                        227./829440.)*i6e8 + (-113005./111476736.)*i4e16 + (-2417./2177280.)*i4e14 + (
                        -949./829440.)*i4e12 + (-7./6912.)*i4e10 + (-5./6912.)*i4e8 + (128441./119439360.)*i2e18 + (
                        22601./18579456.)*i2e16 + (2417./1814400.)*i2e14 + (949./691200.)*i2e12 + (7./5760.)*i2e10 + (
                        1./1152.)*i2e8,
            None,
            -2.,
            2.,
            1.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 0, -3',
            -spin_freq - orbital_freq,
            (8988527./401717565849600.)*i14e6 + (-59123./160526499840.)*i12e8 + (-59123./220723937280.)*i12e6 + (
                        24931463./7134511104000.)*i10e10 + (443047./133772083200.)*i10e8 + (
                        40277./16721510400.)*i10e6 + (-1815893./85614133248.)*i8e12 + (-17951./792723456.)*i8e10 + (
                        -1595./74317824.)*i8e8 + (-145./9289728.)*i8e6 + (7127547349./85614133248000.)*i6e14 + (
                        14214059./152882380800.)*i6e12 + (140513./1415577600.)*i6e10 + (2497./26542080.)*i6e8 + (
                        227./3317760.)*i6e6 + (-147400583./761014517760.)*i4e16 + (-31398887./142690222080.)*i4e14 + (
                        -62617./254803968.)*i4e12 + (-619./2359296.)*i4e10 + (-55./221184.)*i4e8 + (-5./27648.)*i4e6 + (
                        167431204877./821895679180800.)*i2e18 + (147400583./634178764800.)*i2e16 + (
                        31398887./118908518400.)*i2e14 + (62617./212336640.)*i2e12 + (619./1966080.)*i2e10 + (
                        11./36864.)*i2e8 + (1./4608.)*i2e6,
            None,
            -1.,
            2.,
            1.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 0, -1',
            -spin_freq + orbital_freq,
            (2195943977./51218989645824000.)*i18e2 + (3490169./16738231910400.)*i16e4 + (
                        -3490169./4184557977600.)*i16e2 + (8988527./10300450406400.)*i14e6 + (
                        -8988527./2789705318400.)*i14e4 + (8988527./697426329600.)*i14e2 + (
                        -3606503./353158299648.)*i12e8 + (-768599./73574645760.)*i12e6 + (59123./1532805120.)*i12e4 + (
                        -59123./383201280.)*i12e2 + (69800041./1189085184000.)*i10e10 + (
                        2456897./26754416640.)*i10e8 + (523601./5573836800.)*i10e6 + (-40277./116121600.)*i10e4 + (
                        40277./29030400.)*i10e2 + (-8039641./23781703680.)*i8e12 + (-50257./132120576.)*i8e10 + (
                        -44225./74317824.)*i8e8 + (-1885./3096576.)*i8e6 + (145./64512.)*i8e4 + (-145./16128.)*i8e2 + (
                        37519772861./28538044416000.)*i6e14 + (62930983./42467328000.)*i6e12 + (
                        393391./235929600.)*i6e10 + (13847./5308416.)*i6e8 + (2951./1105920.)*i6e6 + (
                        -227./23040.)*i6e4 + (227./5760.)*i6e2 + (-49450862117./15981304872960.)*i4e16 + (
                        -165285343./47563407360.)*i4e14 + (-277229./70778880.)*i4e12 + (-1733./393216.)*i4e10 + (
                        -1525./221184.)*i4e8 + (-65./9216.)*i4e6 + (5./192.)*i4e4 + (-5./48.)*i4e2 + (
                        12842565048623./3835513169510400.)*i2e18 + (49450862117./13317754060800.)*i2e16 + (
                        165285343./39636172800.)*i2e14 + (277229./58982400.)*i2e12 + (1733./327680.)*i2e10 + (
                        305./36864.)*i2e8 + (13./1536.)*i2e6 + (-1./32.)*i2e4 + (1./8.)*i2e2,
            None,
            1.,
            2.,
            1.
        ),
        (
            '-O + 2*n',
            'O - 2*n',
            '2, 1, 0, 0',
            -spin_freq + 2*orbital_freq,
            (-3479571749./486580401635328000.)*i20 + (-2195943977./2560949482291200.)*i18e2 + (
                        2195943977./12804747411456000.)*i18 + (-3490169./132843110400.)*i16e4 + (
                        3490169./209227898880.)*i16e2 + (-3490169./1046139494400.)*i16 + (
                        -710093633./3138418483200.)*i14e6 + (8988527./22140518400.)*i14e4 + (
                        -8988527./34871316480.)*i14e2 + (8988527./174356582400.)*i14 + (
                        -218577731./220723937280.)*i12e8 + (4670717./1724405760.)*i12e6 + (-59123./12165120.)*i12e4 + (
                        59123./19160064.)*i12e2 + (-59123./95800320.)*i12 + (-444698357./139345920000.)*i10e10 + (
                        148904069./16721510400.)*i10e8 + (-3181883./130636800.)*i10e6 + (40277./921600.)*i10e4 + (
                        -40277./1451520.)*i10e2 + (40277./7257600.)*i10 + (-999659./222953472.)*i8e12 + (
                        320189./15482880.)*i8e10 + (-536065./9289728.)*i8e8 + (11455./72576.)*i8e6 + (
                        -145./512.)*i8e4 + (725./4032.)*i8e2 + (-145./4032.)*i8 + (-394966153./48771072000.)*i6e14 + (
                        7824917./398131200.)*i6e12 + (-2506307./27648000.)*i6e10 + (839219./3317760.)*i6e8 + (
                        -17933./25920.)*i6e6 + (1589./1280.)*i6e4 + (-227./288.)*i6e2 + (227./1440.)*i6 + (
                        561889./289013760.)*i4e16 + (1739939./81285120.)*i4e14 + (-34471./663552.)*i4e12 + (
                        11041./46080.)*i4e10 + (-18485./27648.)*i4e8 + (395./216.)*i4e6 + (-105./32.)*i4e4 + (
                        25./12.)*i4e2 + (-5./12.)*i4 + (-459927151./75246796800.)*i2e18 + (
                        -561889./240844800.)*i2e16 + (-1739939./67737600.)*i2e14 + (34471./552960.)*i2e12 + (
                        -11041./38400.)*i2e10 + (3697./4608.)*i2e8 + (-79./36.)*i2e6 + (63./16.)*i2e4 + (
                        -5./2.)*i2e2 + (1./2.)*i2,
            None,
            2.,
            2.,
            1.
        ),
        (
            '-O + 3*n',
            'O - 3*n',
            '2, 1, 0, 1',
            -spin_freq + 3*orbital_freq,
            (2195943977./1045285502976000.)*i18e2 + (143096929./797058662400.)*i16e4 + (
                        -3490169./85399142400.)*i16e2 + (2633638411./595137134592.)*i14e6 + (
                        -368529607./132843110400.)*i14e4 + (8988527./14233190400.)*i14e2 + (
                        2675729611./65399685120.)*i12e8 + (-86615195./1634992128.)*i12e6 + (
                        2424043./72990720.)*i12e4 + (-413861./54743040.)*i12e2 + (24919420177./132120576000.)*i10e10 + (
                        -1822816189./4954521600.)*i10e8 + (11801161./24772608.)*i10e6 + (-1651357./5529600.)*i10e4 + (
                        281939./4147200.)*i10e2 + (22087357./41943040.)*i8e12 + (-17942329./14680064.)*i8e10 + (
                        6562265./2752512.)*i8e8 + (-1062125./344064.)*i8e6 + (5945./3072.)*i8e4 + (
                        -1015./2304.)*i8e2 + (121269339001./150994944000.)*i6e14 + (-1210234837./524288000.)*i6e12 + (
                        140445127./26214400.)*i6e10 + (-10273339./983040.)*i6e8 + (332555./24576.)*i6e6 + (
                        -65149./7680.)*i6e4 + (11123./5760.)*i6e2 + (1324417911./1879048192.)*i4e16 + (
                        -534226163./251658240.)*i4e14 + (15994293./2621440.)*i4e12 + (-1856103./131072.)*i4e10 + (
                        226285./8192.)*i4e8 + (-36625./1024.)*i4e6 + (1435./64.)*i4e4 + (-245./48.)*i4e2 + (
                        6258845529./30064771072.)*i2e18 + (-3973253733./4697620480.)*i2e16 + (
                        534226163./209715200.)*i2e14 + (-47982879./6553600.)*i2e12 + (5568309./327680.)*i2e10 + (
                        -135771./4096.)*i2e8 + (21975./512.)*i2e6 + (-861./32.)*i2e4 + (49./8.)*i2e2,
            None,
            3.,
            2.,
            1.
        ),
        (
            '-O + 4*n',
            'O - 4*n',
            '2, 1, 0, 2',
            -spin_freq + 4*orbital_freq,
            (-1008658841./4184557977600.)*i16e4 + (-3514514057./209227898880.)*i14e6 + (
                        2597684303./697426329600.)*i14e4 + (-4939785773./13795246080.)*i12e8 + (
                        23117093./114960384.)*i12e6 + (-17086547./383201280.)*i12e4 + (
                        -5405535893./1741824000.)*i10e10 + (3365183627./1045094400.)*i10e8 + (
                        -15748307./8709120.)*i10e6 + (11640053./29030400.)*i10e4 + (-244230199./18579456.)*i8e12 + (
                        3892061./193536.)*i8e10 + (-12114895./580608.)*i8e8 + (283475./24192.)*i8e6 + (
                        -41905./16128.)*i8e4 + (-78701030071./2612736000.)*i6e14 + (1911732937./33177600.)*i6e12 + (
                        -30465443./345600.)*i6e10 + (18966077./207360.)*i6e8 + (-88757./1728.)*i6e6 + (
                        65603./5760.)*i6e4 + (-173370469./4976640.)*i4e16 + (346700573./4354560.)*i4e14 + (
                        -8421731./55296.)*i4e12 + (134209./576.)*i4e10 + (-417755./1728.)*i4e8 + (9775./72.)*i4e6 + (
                        -1445./48.)*i4e4 + (-2949969133./182891520.)*i2e18 + (173370469./4147200.)*i2e16 + (
                        -346700573./3628800.)*i2e14 + (8421731./46080.)*i2e12 + (-134209./480.)*i2e10 + (
                        83551./288.)*i2e8 + (-1955./12.)*i2e6 + (289./8.)*i2e4,
            None,
            4.,
            2.,
            1.
        ),
        (
            '-O + 5*n',
            'O - 5*n',
            '2, 1, 0, 3',
            -spin_freq + 5*orbital_freq,
            (19747793819./1236054048768.)*i14e6 + (324982872175./353158299648.)*i12e8 + (
                        -8443060015./44144787456.)*i12e6 + (946067057155./57076088832.)*i10e10 + (
                        -44278318565./5350883328.)*i10e8 + (1150351397./668860416.)*i10e6 + (
                        10537294474375./85614133248.)*i8e12 + (-85147679375./792723456.)*i8e10 + (
                        3985125625./74317824.)*i8e8 + (-103533625./9289728.)*i8e6 + (
                        291895467893825./684913065984.)*i6e14 + (-3299263235425./6115295232.)*i6e12 + (
                        26660032025./56623104.)*i6e10 + (-1247756575./5308416.)*i6e8 + (32416735./663552.)*i6e6 + (
                        106374143767625./152202903552.)*i4e16 + (-32147077961875./28538044416.)*i4e14 + (
                        363354981875./254803968.)*i4e12 + (-2936126875./2359296.)*i4e10 + (137418125./221184.)*i4e8 + (
                        -3570125./27648.)*i4e6 + (2044426346565875./4696546738176.)*i2e18 + (
                        -21274828753525./25367150592.)*i2e16 + (6429415592375./4756340736.)*i2e14 + (
                        -72670996375./42467328.)*i2e12 + (587225375./393216.)*i2e10 + (-27483625./36864.)*i2e8 + (
                        714025./4608.)*i2e6,
            None,
            5.,
            2.,
            1.
        ),
        (
            '-O + 6*n',
            'O - 6*n',
            '2, 1, 0, 4',
            -spin_freq + 6*orbital_freq,
            (-16796193947./24524881920.)*i12e8 + (-98944357369./3096576000.)*i10e10 + (
                        11442252653./1857945600.)*i10e8 + (-6378746359./13762560.)*i8e12 + (71241313./344064.)*i8e10 + (
                        -41192905./1032192.)*i8e8 + (-231216279613./86016000.)*i6e14 + (
                        49930187017./24576000.)*i6e12 + (-557647519./614400.)*i6e10 + (64488203./368640.)*i6e8 + (
                        -60181343991./9175040.)*i4e16 + (1018573919./143360.)*i4e14 + (-219956771./40960.)*i4e12 + (
                        2456597./1024.)*i4e10 + (-1420445./3072.)*i4e8 + (-25998653133./4587520.)*i2e18 + (
                        180544031973./22937600.)*i2e16 + (-3055721757./358400.)*i2e14 + (659870313./102400.)*i2e12 + (
                        -7369791./2560.)*i2e10 + (284089./512.)*i2e8,
            None,
            6.,
            2.,
            1.
        ),
        (
            '-O + 7*n',
            'O - 7*n',
            '2, 1, 0, 5',
            -spin_freq + 7*orbital_freq,
            (300019646853899./15288238080000.)*i10e10 + (2905264589675./4076863488.)*i8e12 + (
                        -216018317123./1698693120.)*i8e10 + (94728014655844073./12230590464000.)*i6e14 + (
                        -6367538528267./2038431744.)*i6e12 + (11836313996843./21233664000.)*i6e10 + (
                        9997999669389767./326149079040.)*i4e16 + (-417304029320899./20384317440.)*i4e14 + (
                        701270763025./84934656.)*i4e12 + (-52142352409./35389440.)*i4e10 + (
                        1518065224694500853./39137889484800.)*i2e18 + (-9997999669389767./271790899200.)*i2e16 + (
                        417304029320899./16986931200.)*i2e14 + (-140254152605./14155776.)*i2e12 + (
                        52142352409./29491200.)*i2e10,
            None,
            7.,
            2.,
            1.
        ),
        (
            '-O + 8*n',
            'O - 8*n',
            '2, 1, 0, 6',
            -spin_freq + 8*orbital_freq,
            (-156107294669./418037760.)*i8e12 + (-12919549337701./1306368000.)*i6e14 + (
                        1221943306547./746496000.)*i6e12 + (-578658802849./8128512.)*i4e16 + (
                        56914314263./2177280.)*i4e14 + (-5383010161./1244160.)*i4e12 + (
                        -391340609035087./2743372800.)*i2e18 + (578658802849./6773760.)*i2e16 + (
                        -56914314263./1814400.)*i2e14 + (5383010161./1036800.)*i2e12,
            None,
            8.,
            2.,
            1.
        ),
        (
            '-O + 9*n',
            'O - 9*n',
            '2, 1, 0, 7',
            -spin_freq + 9*orbital_freq,
            (3719858248951787./822083584000.)*i6e14 + (1021582839100239./13153337344.)*i4e16 + (
                        -49161122232843./4110417920.)*i4e14 + (415947859083950607./1503238553600.)*i2e18 + (
                        -3064748517300717./32883343360.)*i2e16 + (147483366698529./10276044800.)*i2e14,
            None,
            9.,
            2.,
            1.
        ),
        (
            '-O + 10*n',
            'O - 10*n',
            '2, 1, 0, 8',
            -spin_freq + 10*orbital_freq,
            (-2010318936125./63700992.)*i4e16 + (-176187983600875./668860416.)*i2e18 + (402063787225./10616832.)*i2e16,
            None,
            10.,
            2.,
            1.
        ),
        (
            '-O + 11*n',
            'O - 11*n',
            '2, 1, 0, 9',
            -spin_freq + 11*orbital_freq,
            (6648821549377771726369./69039237051187200.)*i2e18,
            None,
            11.,
            2.,
            1.
        ),
        (
            '-O - 9*n',
            'O + 9*n',
            '2, 1, 1, -9',
            -spin_freq - 9*orbital_freq,
            (4143587919679849./10522669875200.)*i2e18,
            None,
            -9.,
            0.,
            1.
        ),
        (
            '-O - 8*n',
            'O + 8*n',
            '2, 1, 1, -8',
            -spin_freq - 8*orbital_freq,
            (-31801945561./120422400.)*i4e16 + (-1606626956771./4335206400.)*i2e18 + (31801945561./160563200.)*i2e16,
            None,
            -8.,
            0.,
            1.
        ),
        (
            '-O - 7*n',
            'O + 7*n',
            '2, 1, 1, -7',
            -spin_freq - 7*orbital_freq,
            (186702632281./2654208000.)*i6e14 + (1328179895713./7549747200.)*i4e16 + (
                        -186702632281./1415577600.)*i4e14 + (999990833612299./5798205849600.)*i2e18 + (
                        -1328179895713./10066329600.)*i2e16 + (186702632281./1887436800.)*i2e14,
            None,
            -7.,
            0.,
            1.
        ),
        (
            '-O - 6*n',
            'O + 6*n',
            '2, 1, 1, -6',
            -spin_freq - 6*orbital_freq,
            (-10029889./1008000.)*i8e12 + (-769581./28000.)*i6e14 + (10029889./288000.)*i6e12 + (
                        -65845785./802816.)*i4e16 + (2308743./44800.)*i4e14 + (-10029889./153600.)*i4e12 + (
                        -941362189./80281600.)*i2e18 + (197537355./3211264.)*i2e16 + (-6926229./179200.)*i2e14 + (
                        10029889./204800.)*i2e12,
            None,
            -6.,
            0.,
            1.
        ),
        (
            '-O - 5*n',
            'O + 5*n',
            '2, 1, 1, -5',
            -spin_freq - 5*orbital_freq,
            (38809./44800.)*i10e10 + (982439./860160.)*i8e12 + (-349281./71680.)*i8e10 + (
                        13524196093./743178240.)*i6e14 + (-982439./245760.)*i6e12 + (349281./20480.)*i6e10 + (
                        -28491858875./2113929216.)*i4e16 + (-13524196093./396361728.)*i4e14 + (
                        982439./131072.)*i4e12 + (-1047843./32768.)*i4e10 + (2447539096445./105226698752.)*i2e18 + (
                        28491858875./2818572288.)*i2e16 + (13524196093./528482304.)*i2e14 + (
                        -2947317./524288.)*i2e12 + (3143529./131072.)*i2e10,
            None,
            -5.,
            0.,
            1.
        ),
        (
            '-O - 4*n',
            'O + 4*n',
            '2, 1, 1, -4',
            -spin_freq - 4*orbital_freq,
            (-308./6075.)*i12e8 + (473./3375.)*i10e10 + (847./2025.)*i10e8 + (-708871./252000.)*i8e12 + (
                        -473./600.)*i8e10 + (-847./360.)*i8e8 + (3027367./324000.)*i6e14 + (708871./72000.)*i6e12 + (
                        3311./1200.)*i6e10 + (5929./720.)*i6e8 + (-124914751./5160960.)*i4e16 + (
                        -3027367./172800.)*i4e14 + (-708871./38400.)*i4e12 + (-3311./640.)*i4e10 + (
                        -5929./384.)*i4e8 + (2170376447./103219200.)*i2e18 + (124914751./6881280.)*i2e16 + (
                        3027367./230400.)*i2e14 + (708871./51200.)*i2e12 + (9933./2560.)*i2e10 + (5929./512.)*i2e8,
            None,
            -4.,
            0.,
            1.
        ),
        (
            '-O - 3*n',
            'O + 3*n',
            '2, 1, 1, -3',
            -spin_freq - 3*orbital_freq,
            (89888./42567525.)*i14e6 + (-6943./311850.)*i12e8 + (-11236./467775.)*i12e6 + (286661./864000.)*i10e10 + (
                        6943./37800.)*i10e8 + (2809./14175.)*i10e6 + (-30355211./12902400.)*i8e12 + (
                        -286661./153600.)*i8e10 + (-6943./6720.)*i8e8 + (-2809./2520.)*i8e6 + (
                        351714161./32768000.)*i6e14 + (30355211./3686400.)*i6e12 + (2006627./307200.)*i6e10 + (
                        6943./1920.)*i6e8 + (2809./720.)*i6e6 + (-145747579353./5872025600.)*i4e16 + (
                        -1055142483./52428800.)*i4e14 + (-30355211./1966080.)*i4e12 + (-2006627./163840.)*i4e10 + (
                        -6943./1024.)*i4e8 + (-2809./384.)*i4e6 + (8454941691423./375809638400.)*i2e18 + (
                        437242738059./23488102400.)*i2e16 + (3165427449./209715200.)*i2e14 + (
                        30355211./2621440.)*i2e12 + (6019881./655360.)*i2e10 + (20829./4096.)*i2e8 + (2809./512.)*i2e6,
            None,
            -3.,
            0.,
            1.
        ),
        (
            '-O - 2*n',
            'O + 2*n',
            '2, 1, 1, -2',
            -spin_freq - 2*orbital_freq,
            (-512./7882875.)*i16e4 + (1024./675675.)*i14e6 + (512./525525.)*i14e4 + (-1208./42525.)*i12e8 + (
                        -128./7425.)*i12e6 + (-64./5775.)*i12e4 + (7838./23625.)*i10e10 + (3322./14175.)*i10e8 + (
                        32./225.)*i10e6 + (16./175.)*i10e4 + (-505601./201600.)*i8e12 + (-3919./2100.)*i8e10 + (
                        -1661./1260.)*i8e8 + (-4./5.)*i8e6 + (-18./35.)*i8e4 + (25565893./2268000.)*i6e14 + (
                        505601./57600.)*i6e12 + (3919./600.)*i6e10 + (1661./360.)*i6e8 + (14./5.)*i6e6 + (
                        9./5.)*i6e4 + (-678544541./25804800.)*i4e16 + (-25565893./1209600.)*i4e14 + (
                        -505601./30720.)*i4e12 + (-3919./320.)*i4e10 + (-1661./192.)*i4e8 + (-21./4.)*i4e6 + (
                        -27./8.)*i4e4 + (51883919761./2167603200.)*i2e18 + (678544541./34406400.)*i2e16 + (
                        25565893./1612800.)*i2e14 + (505601./40960.)*i2e12 + (11757./1280.)*i2e10 + (
                        1661./256.)*i2e8 + (63./16.)*i2e6 + (81./32.)*i2e4,
            None,
            -2.,
            0.,
            1.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 1, -1',
            -spin_freq - orbital_freq,
            (16384./10854718875.)*i18e2 + (-512./7882875.)*i16e4 + (-2048./70945875.)*i16e2 + (544./315315.)*i14e6 + (
                        512./525525.)*i14e4 + (2048./4729725.)*i14e2 + (-28211./935550.)*i12e8 + (-68./3465.)*i12e6 + (
                        -64./5775.)*i12e4 + (-256./51975.)*i12e2 + (1064887./3024000.)*i10e10 + (
                        28211./113400.)*i10e8 + (17./105.)*i10e6 + (16./175.)*i10e4 + (64./1575.)*i10e2 + (
                        -85553819./32256000.)*i8e12 + (-1064887./537600.)*i8e10 + (-28211./20160.)*i8e8 + (
                        -51./56.)*i8e6 + (-18./35.)*i8e4 + (-8./35.)*i8e2 + (24654653741./2064384000.)*i6e14 + (
                        85553819./9216000.)*i6e12 + (1064887./153600.)*i6e10 + (28211./5760.)*i6e8 + (51./16.)*i6e6 + (
                        9./5.)*i6e4 + (4./5.)*i6e2 + (-689312857627./24662507520.)*i4e16 + (
                        -24654653741./1101004800.)*i4e14 + (-85553819./4915200.)*i4e12 + (-1064887./81920.)*i4e10 + (
                        -28211./3072.)*i4e8 + (-765./128.)*i4e6 + (-27./8.)*i4e4 + (-3./2.)*i4e2 + (
                        725941889609009./28411208663040.)*i2e18 + (689312857627./32883343360.)*i2e16 + (
                        24654653741./1468006400.)*i2e14 + (85553819./6553600.)*i2e12 + (3194661./327680.)*i2e10 + (
                        28211./4096.)*i2e8 + (2295./512.)*i2e6 + (81./32.)*i2e4 + (9./8.)*i2e2,
            None,
            -1.,
            0.,
            1.
        ),
        (
            '-O',
            'O',
            '2, 1, 1, 0',
            -spin_freq,
            (-262144./9280784638125.)*i20 + (65536./32564156625.)*i18e2 + (65536./97692469875.)*i18 + (
                        -16384./212837625.)*i16e4 + (-8192./212837625.)*i16e2 + (-8192./638512875.)*i16 + (
                        16384./8513505.)*i14e6 + (16384./14189175.)*i14e4 + (8192./14189175.)*i14e2 + (
                        8192./42567525.)*i14 + (-1024./31185.)*i12e8 + (-2048./93555.)*i12e6 + (
                        -2048./155925.)*i12e4 + (-1024./155925.)*i12e2 + (-1024./467775.)*i12 + (256./675.)*i10e10 + (
                        256./945.)*i10e8 + (512./2835.)*i10e6 + (512./4725.)*i10e4 + (256./4725.)*i10e2 + (
                        256./14175.)*i10 + (-128./45.)*i8e12 + (-32./15.)*i8e10 + (-32./21.)*i8e8 + (-64./63.)*i8e6 + (
                        -64./105.)*i8e4 + (-32./105.)*i8e2 + (-32./315.)*i8 + (64./5.)*i6e14 + (448./45.)*i6e12 + (
                        112./15.)*i6e10 + (16./3.)*i6e8 + (32./9.)*i6e6 + (32./15.)*i6e4 + (16./15.)*i6e2 + (
                        16./45.)*i6 + -30.*i4e16 + -24.*i4e14 + (-56./3.)*i4e12 + -14.*i4e10 + -10.*i4e8 + (
                        -20./3.)*i4e6 + -4.*i4e4 + -2.*i4e2 + (-2./3.)*i4 + (55./2.)*i2e18 + (
                        45./2.)*i2e16 + 18.*i2e14 + 14.*i2e12 + (21./2.)*i2e10 + (15./2.)*i2e8 + 5.*i2e6 + 3.*i2e4 + (
                        3./2.)*i2e2 + (1./2.)*i2,
            None,
            0.,
            0.,
            1.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 1, 1',
            -spin_freq + orbital_freq,
            (16384./10854718875.)*i18e2 + (-512./7882875.)*i16e4 + (-2048./70945875.)*i16e2 + (544./315315.)*i14e6 + (
                        512./525525.)*i14e4 + (2048./4729725.)*i14e2 + (-28211./935550.)*i12e8 + (-68./3465.)*i12e6 + (
                        -64./5775.)*i12e4 + (-256./51975.)*i12e2 + (1064887./3024000.)*i10e10 + (
                        28211./113400.)*i10e8 + (17./105.)*i10e6 + (16./175.)*i10e4 + (64./1575.)*i10e2 + (
                        -85553819./32256000.)*i8e12 + (-1064887./537600.)*i8e10 + (-28211./20160.)*i8e8 + (
                        -51./56.)*i8e6 + (-18./35.)*i8e4 + (-8./35.)*i8e2 + (24654653741./2064384000.)*i6e14 + (
                        85553819./9216000.)*i6e12 + (1064887./153600.)*i6e10 + (28211./5760.)*i6e8 + (51./16.)*i6e6 + (
                        9./5.)*i6e4 + (4./5.)*i6e2 + (-689312857627./24662507520.)*i4e16 + (
                        -24654653741./1101004800.)*i4e14 + (-85553819./4915200.)*i4e12 + (-1064887./81920.)*i4e10 + (
                        -28211./3072.)*i4e8 + (-765./128.)*i4e6 + (-27./8.)*i4e4 + (-3./2.)*i4e2 + (
                        725941889609009./28411208663040.)*i2e18 + (689312857627./32883343360.)*i2e16 + (
                        24654653741./1468006400.)*i2e14 + (85553819./6553600.)*i2e12 + (3194661./327680.)*i2e10 + (
                        28211./4096.)*i2e8 + (2295./512.)*i2e6 + (81./32.)*i2e4 + (9./8.)*i2e2,
            None,
            1.,
            0.,
            1.
        ),
        (
            '-O + 2*n',
            'O - 2*n',
            '2, 1, 1, 2',
            -spin_freq + 2*orbital_freq,
            (-512./7882875.)*i16e4 + (1024./675675.)*i14e6 + (512./525525.)*i14e4 + (-1208./42525.)*i12e8 + (
                        -128./7425.)*i12e6 + (-64./5775.)*i12e4 + (7838./23625.)*i10e10 + (3322./14175.)*i10e8 + (
                        32./225.)*i10e6 + (16./175.)*i10e4 + (-505601./201600.)*i8e12 + (-3919./2100.)*i8e10 + (
                        -1661./1260.)*i8e8 + (-4./5.)*i8e6 + (-18./35.)*i8e4 + (25565893./2268000.)*i6e14 + (
                        505601./57600.)*i6e12 + (3919./600.)*i6e10 + (1661./360.)*i6e8 + (14./5.)*i6e6 + (
                        9./5.)*i6e4 + (-678544541./25804800.)*i4e16 + (-25565893./1209600.)*i4e14 + (
                        -505601./30720.)*i4e12 + (-3919./320.)*i4e10 + (-1661./192.)*i4e8 + (-21./4.)*i4e6 + (
                        -27./8.)*i4e4 + (51883919761./2167603200.)*i2e18 + (678544541./34406400.)*i2e16 + (
                        25565893./1612800.)*i2e14 + (505601./40960.)*i2e12 + (11757./1280.)*i2e10 + (
                        1661./256.)*i2e8 + (63./16.)*i2e6 + (81./32.)*i2e4,
            None,
            2.,
            0.,
            1.
        ),
        (
            '-O + 3*n',
            'O - 3*n',
            '2, 1, 1, 3',
            -spin_freq + 3*orbital_freq,
            (89888./42567525.)*i14e6 + (-6943./311850.)*i12e8 + (-11236./467775.)*i12e6 + (286661./864000.)*i10e10 + (
                        6943./37800.)*i10e8 + (2809./14175.)*i10e6 + (-30355211./12902400.)*i8e12 + (
                        -286661./153600.)*i8e10 + (-6943./6720.)*i8e8 + (-2809./2520.)*i8e6 + (
                        351714161./32768000.)*i6e14 + (30355211./3686400.)*i6e12 + (2006627./307200.)*i6e10 + (
                        6943./1920.)*i6e8 + (2809./720.)*i6e6 + (-145747579353./5872025600.)*i4e16 + (
                        -1055142483./52428800.)*i4e14 + (-30355211./1966080.)*i4e12 + (-2006627./163840.)*i4e10 + (
                        -6943./1024.)*i4e8 + (-2809./384.)*i4e6 + (8454941691423./375809638400.)*i2e18 + (
                        437242738059./23488102400.)*i2e16 + (3165427449./209715200.)*i2e14 + (
                        30355211./2621440.)*i2e12 + (6019881./655360.)*i2e10 + (20829./4096.)*i2e8 + (2809./512.)*i2e6,
            None,
            3.,
            0.,
            1.
        ),
        (
            '-O + 4*n',
            'O - 4*n',
            '2, 1, 1, 4',
            -spin_freq + 4*orbital_freq,
            (-308./6075.)*i12e8 + (473./3375.)*i10e10 + (847./2025.)*i10e8 + (-708871./252000.)*i8e12 + (
                        -473./600.)*i8e10 + (-847./360.)*i8e8 + (3027367./324000.)*i6e14 + (708871./72000.)*i6e12 + (
                        3311./1200.)*i6e10 + (5929./720.)*i6e8 + (-124914751./5160960.)*i4e16 + (
                        -3027367./172800.)*i4e14 + (-708871./38400.)*i4e12 + (-3311./640.)*i4e10 + (
                        -5929./384.)*i4e8 + (2170376447./103219200.)*i2e18 + (124914751./6881280.)*i2e16 + (
                        3027367./230400.)*i2e14 + (708871./51200.)*i2e12 + (9933./2560.)*i2e10 + (5929./512.)*i2e8,
            None,
            4.,
            0.,
            1.
        ),
        (
            '-O + 5*n',
            'O - 5*n',
            '2, 1, 1, 5',
            -spin_freq + 5*orbital_freq,
            (38809./44800.)*i10e10 + (982439./860160.)*i8e12 + (-349281./71680.)*i8e10 + (
                        13524196093./743178240.)*i6e14 + (-982439./245760.)*i6e12 + (349281./20480.)*i6e10 + (
                        -28491858875./2113929216.)*i4e16 + (-13524196093./396361728.)*i4e14 + (
                        982439./131072.)*i4e12 + (-1047843./32768.)*i4e10 + (2447539096445./105226698752.)*i2e18 + (
                        28491858875./2818572288.)*i2e16 + (13524196093./528482304.)*i2e14 + (
                        -2947317./524288.)*i2e12 + (3143529./131072.)*i2e10,
            None,
            5.,
            0.,
            1.
        ),
        (
            '-O + 6*n',
            'O - 6*n',
            '2, 1, 1, 6',
            -spin_freq + 6*orbital_freq,
            (-10029889./1008000.)*i8e12 + (-769581./28000.)*i6e14 + (10029889./288000.)*i6e12 + (
                        -65845785./802816.)*i4e16 + (2308743./44800.)*i4e14 + (-10029889./153600.)*i4e12 + (
                        -941362189./80281600.)*i2e18 + (197537355./3211264.)*i2e16 + (-6926229./179200.)*i2e14 + (
                        10029889./204800.)*i2e12,
            None,
            6.,
            0.,
            1.
        ),
        (
            '-O + 7*n',
            'O - 7*n',
            '2, 1, 1, 7',
            -spin_freq + 7*orbital_freq,
            (186702632281./2654208000.)*i6e14 + (1328179895713./7549747200.)*i4e16 + (
                        -186702632281./1415577600.)*i4e14 + (999990833612299./5798205849600.)*i2e18 + (
                        -1328179895713./10066329600.)*i2e16 + (186702632281./1887436800.)*i2e14,
            None,
            7.,
            0.,
            1.
        ),
        (
            '-O + 8*n',
            'O - 8*n',
            '2, 1, 1, 8',
            -spin_freq + 8*orbital_freq,
            (-31801945561./120422400.)*i4e16 + (-1606626956771./4335206400.)*i2e18 + (31801945561./160563200.)*i2e16,
            None,
            8.,
            0.,
            1.
        ),
        (
            '-O + 9*n',
            'O - 9*n',
            '2, 1, 1, 9',
            -spin_freq + 9*orbital_freq,
            (4143587919679849./10522669875200.)*i2e18,
            None,
            9.,
            0.,
            1.
        ),
        (
            '-O - 9*n',
            'O + 9*n',
            '2, 1, 2, -7',
            -spin_freq - 9*orbital_freq,
            (147483366698529./164416716800.)*i6e14,
            None,
            -9.,
            -2.,
            1.
        ),
        (
            '-O - 8*n',
            'O + 8*n',
            '2, 1, 2, -6',
            -spin_freq - 8*orbital_freq,
            (-5383010161./33177600.)*i8e12 + (-56914314263./29030400.)*i6e14 + (5383010161./16588800.)*i6e12,
            None,
            -8.,
            -2.,
            1.
        ),
        (
            '-O - 7*n',
            'O + 7*n',
            '2, 1, 2, -5',
            -spin_freq - 7*orbital_freq,
            (52142352409./4194304000.)*i10e10 + (140254152605./452984832.)*i8e12 + (-52142352409./943718400.)*i8e10 + (
                        417304029320899./271790899200.)*i6e14 + (-140254152605./226492416.)*i6e12 + (
                        52142352409./471859200.)*i6e10,
            None,
            -7.,
            -2.,
            1.
        ),
        (
            '-O - 6*n',
            'O + 6*n',
            '2, 1, 2, -4',
            -spin_freq - 6*orbital_freq,
            (-131533207./247726080.)*i12e8 + (-66328119./3276800.)*i10e10 + (2556801./655360.)*i10e8 + (
                        -659870313./3276800.)*i8e12 + (7369791./81920.)*i8e10 + (-284089./16384.)*i8e8 + (
                        -3055721757./5734400.)*i6e14 + (659870313./1638400.)*i6e12 + (-7369791./40960.)*i6e10 + (
                        284089./8192.)*i6e8,
            None,
            -6.,
            -2.,
            1.
        ),
        (
            '-O - 5*n',
            'O + 5*n',
            '2, 1, 2, -3',
            -spin_freq - 5*orbital_freq,
            (24705265./1783627776.)*i14e6 + (2544983675./3567255552.)*i12e8 + (-66118715./445906944.)*i12e6 + (
                        352335225./33554432.)*i10e10 + (-5496725./1048576.)*i10e8 + (142805./131072.)*i10e6 + (
                        72670996375./1358954496.)*i8e12 + (-587225375./12582912.)*i8e10 + (27483625./1179648.)*i8e8 + (
                        -714025./147456.)*i8e6 + (6429415592375./76101451776.)*i6e14 + (
                        -72670996375./679477248.)*i6e12 + (587225375./6291456.)*i6e10 + (-27483625./589824.)*i6e8 + (
                        714025./73728.)*i6e6,
            None,
            -5.,
            -2.,
            1.
        ),
        (
            '-O - 4*n',
            'O + 4*n',
            '2, 1, 2, -2',
            -spin_freq - 4*orbital_freq,
            (-126293./567705600.)*i16e4 + (-67643./4644864.)*i14e6 + (49997./15482880.)*i14e4 + (
                        -38684113./139345920.)*i12e8 + (181033./1161216.)*i12e6 + (-133807./3870720.)*i12e4 + (
                        -402627./204800.)*i10e10 + (83551./40960.)*i10e8 + (-1173./1024.)*i10e6 + (
                        2601./10240.)*i10e4 + (-8421731./1474560.)*i8e12 + (134209./15360.)*i8e10 + (
                        -83551./9216.)*i8e8 + (1955./384.)*i8e6 + (-289./256.)*i8e4 + (-346700573./58060800.)*i6e14 + (
                        8421731./737280.)*i6e12 + (-134209./7680.)*i6e10 + (83551./4608.)*i6e8 + (-1955./192.)*i6e6 + (
                        289./128.)*i6e4,
            None,
            -4.,
            -2.,
            1.
        ),
        (
            '-O - 3*n',
            'O + 3*n',
            '2, 1, 2, -1',
            -spin_freq - 3*orbital_freq,
            (2743907./1366386278400.)*i18e2 + (17917./108134400.)*i16e4 + (-3059./81100800.)*i16e2 + (
                        253445./66060288.)*i14e6 + (-7093./2949120.)*i14e4 + (1211./2211840.)*i14e2 + (
                        20953991./660602880.)*i12e8 + (-678295./16515072.)*i12e6 + (18983./737280.)*i12e4 + (
                        -3241./552960.)*i12e2 + (50114781./419430400.)*i10e10 + (-1221939./5242880.)*i10e8 + (
                        39555./131072.)*i10e6 + (-7749./40960.)*i10e4 + (441./10240.)*i10e2 + (
                        47982879./209715200.)*i8e12 + (-5568309./10485760.)*i8e10 + (135771./131072.)*i8e8 + (
                        -21975./16384.)*i8e6 + (861./1024.)*i8e4 + (-49./256.)*i8e2 + (534226163./3355443200.)*i6e14 + (
                        -47982879./104857600.)*i6e12 + (5568309./5242880.)*i6e10 + (-135771./65536.)*i6e8 + (
                        21975./8192.)*i6e6 + (-861./512.)*i6e4 + (49./128.)*i6e2,
            None,
            -3.,
            -2.,
            1.
        ),
        (
            '-O - 2*n',
            'O + 2*n',
            '2, 1, 2, 0',
            -spin_freq - 2*orbital_freq,
            (-166711./23911759872000.)*i20 + (-2743907./3347646382080.)*i18e2 + (2743907./16738231910400.)*i18 + (
                        -437./18022400.)*i16e4 + (437./28385280.)*i16e2 + (-437./141926400.)*i16 + (
                        -3287./17418240.)*i14e6 + (173./491520.)*i14e4 + (-173./774144.)*i14e2 + (173./3870720.)*i14 + (
                        -27317./63700992.)*i12e8 + (8797./4354560.)*i12e6 + (-463./122880.)*i12e4 + (
                        463./193536.)*i12e2 + (-463./967680.)*i12 + (12237./16384000.)*i10e10 + (413./131072.)*i10e8 + (
                        -19./1280.)*i10e6 + (567./20480.)*i10e4 + (-9./512.)*i10e2 + (9./2560.)*i10 + (
                        2323./5898240.)*i8e12 + (-4079./1228800.)*i8e10 + (-2065./147456.)*i8e8 + (19./288.)*i8e6 + (
                        -63./512.)*i8e4 + (5./64.)*i8e2 + (-1./64.)*i8 + (689797./361267200.)*i6e14 + (
                        -2323./2949120.)*i6e12 + (4079./614400.)*i6e10 + (2065./73728.)*i6e8 + (-19./144.)*i6e6 + (
                        63./256.)*i6e4 + (-5./32.)*i6e2 + (1./32.)*i6,
            None,
            -2.,
            -2.,
            1.
        ),
        (
            '-O - n',
            'O + n',
            '2, 1, 2, 1',
            -spin_freq - orbital_freq,
            (2743907./66952927641600.)*i18e2 + (437./2270822400.)*i16e4 + (-437./567705600.)*i16e2 + (
                        2249./2972712960.)*i14e6 + (-173./61931520.)*i14e4 + (173./15482880.)*i14e2 + (
                        -28243./3567255552.)*i12e8 + (-6019./743178240.)*i12e6 + (463./15482880.)*i12e4 + (
                        -463./3870720.)*i12e2 + (15597./419430400.)*i10e10 + (61./1048576.)*i10e8 + (
                        39./655360.)*i10e6 + (-9./40960.)*i10e4 + (9./10240.)*i10e2 + (-277229./1887436800.)*i8e12 + (
                        -1733./10485760.)*i8e10 + (-305./1179648.)*i8e8 + (-13./49152.)*i8e6 + (1./1024.)*i8e4 + (
                        -1./256.)*i8e2 + (165285343./634178764800.)*i6e14 + (277229./943718400.)*i6e12 + (
                        1733./5242880.)*i6e10 + (305./589824.)*i6e8 + (13./24576.)*i6e6 + (-1./512.)*i6e4 + (
                        1./128.)*i6e2,
            None,
            -1.,
            -2.,
            1.
        ),
        (
            '-O + n',
            'O - n',
            '2, 1, 2, 3',
            -spin_freq + orbital_freq,
            (173./8918138880.)*i14e6 + (-5093./17836277760.)*i12e8 + (-463./2229534720.)*i12e6 + (
                        1857./838860800.)*i10e10 + (11./5242880.)*i10e8 + (1./655360.)*i10e6 + (
                        -62617./6794772480.)*i8e12 + (-619./62914560.)*i8e10 + (-11./1179648.)*i8e8 + (
                        -1./147456.)*i8e6 + (31398887./1902536294400.)*i6e14 + (62617./3397386240.)*i6e12 + (
                        619./31457280.)*i6e10 + (11./589824.)*i6e8 + (1./73728.)*i6e6,
            None,
            1.,
            -2.,
            1.
        ),
        (
            '-O + 2*n',
            'O - 2*n',
            '2, 1, 2, 4',
            -spin_freq + 2*orbital_freq,
            (-463./557383680.)*i12e8 + (7./819200.)*i10e10 + (1./163840.)*i10e8 + (-949./22118400.)*i8e12 + (
                        -7./184320.)*i8e10 + (-1./36864.)*i8e8 + (2417./29030400.)*i6e14 + (949./11059200.)*i6e12 + (
                        7./92160.)*i6e10 + (1./18432.)*i6e8,
            None,
            2.,
            -2.,
            1.
        ),
        (
            '-O + 3*n',
            'O - 3*n',
            '2, 1, 2, 5',
            -spin_freq + 3*orbital_freq,
            (59049./4194304000.)*i10e10 + (-6561./83886080.)*i8e12 + (-6561./104857600.)*i8e10 + (
                        4284333./23488102400.)*i6e14 + (6561./41943040.)*i6e12 + (6561./52428800.)*i6e10,
            None,
            3.,
            -2.,
            1.
        ),
        (
            '-O + 4*n',
            'O - 4*n',
            '2, 1, 2, 6',
            -spin_freq + 4*orbital_freq,
            (-1./8100.)*i8e12 + (1./4050.)*i6e14 + (1./4050.)*i6e12,
            None,
            4.,
            -2.,
            1.
        ),
        (
            '-O + 5*n',
            'O - 5*n',
            '2, 1, 2, 7',
            -spin_freq + 5*orbital_freq,
            (244140625./532710162432.)*i6e14,
            None,
            5.,
            -2.,
            1.
        ),
        (
            '-2*O - 8*n',
            '2*O + 8*n',
            '2, 2, 0, -10',
            -2*spin_freq - 8*orbital_freq,
            (8388608./200930625.)*e20,
            None,
            -8.,
            2.,
            2.
        ),
        (
            '-2*O - 7*n',
            '2*O + 7*n',
            '2, 2, 0, -9',
            -2*spin_freq - 7*orbital_freq,
            (-33232930569601./1408964021452800.)*i2e18 + (-33232930569601./28179280429056000.)*e20 + (
                        33232930569601./1408964021452800.)*e18,
            None,
            -7.,
            2.,
            2.
        ),
        (
            '-2*O - 6*n',
            '2*O + 6*n',
            '2, 2, 0, -8',
            -2*spin_freq - 6*orbital_freq,
            (1948617./321126400.)*i4e16 + (-177147./40140800.)*i2e18 + (-531441./40140800.)*i2e16 + (
                        18128043./1605632000.)*e20 + (177147./40140800.)*e18 + (531441./40140800.)*e16,
            None,
            -6.,
            2.,
            2.
        ),
        (
            '-2*O - 5*n',
            '2*O + 5*n',
            '2, 2, 0, -7',
            -2*spin_freq - 5*orbital_freq,
            (-1123046875./1198597865472.)*i6e14 + (29541015625./12785043898368.)*i4e16 + (
                        2685546875./799065243648.)*i4e14 + (-323974609375./43834436222976.)*i2e18 + (
                        -2685546875./532710162432.)*i2e16 + (-244140625./33294385152.)*i2e14 + (
                        21589111328125./3682092642729984.)*e20 + (323974609375./43834436222976.)*e18 + (
                        2685546875./532710162432.)*e16 + (244140625./33294385152.)*e14,
            None,
            -5.,
            2.,
            2.
        ),
        (
            '-2*O - 4*n',
            '2*O + 4*n',
            '2, 2, 0, -6',
            -2*spin_freq - 4*orbital_freq,
            (1957./20412000.)*i8e12 + (-46./91125.)*i6e14 + (-46./91125.)*i6e12 + (253./113400.)*i4e16 + (
                        11./6075.)*i4e14 + (11./6075.)*i4e12 + (-68./15309.)*i2e18 + (-23./4725.)*i2e16 + (
                        -8./2025.)*i2e14 + (-8./2025.)*i2e12 + (1733531./428652000.)*e20 + (68./15309.)*e18 + (
                        23./4725.)*e16 + (8./2025.)*e14 + (8./2025.)*e12,
            None,
            -4.,
            2.,
            2.
        ),
        (
            '-2*O - 3*n',
            '2*O + 3*n',
            '2, 2, 0, -5',
            -2*spin_freq - 3*orbital_freq,
            (-980667./146800640000.)*i10e10 + (1426653./23488102400.)*i8e12 + (1426653./29360128000.)*i8e10 + (
                        -10948851./29360128000.)*i6e14 + (-16767./52428800.)*i6e12 + (-16767./65536000.)*i6e10 + (
                        34569909./26843545600.)*i4e16 + (15709221./11744051200.)*i4e14 + (24057./20971520.)*i4e12 + (
                        24057./26214400.)*i4e10 + (-270241029./105226698752.)*i2e18 + (-9428157./3355443200.)*i2e16 + (
                        -4284333./1468006400.)*i2e14 + (-6561./2621440.)*i2e12 + (-6561./3276800.)*i2e10 + (
                        591397458471./263066746880000.)*e20 + (270241029./105226698752.)*e18 + (
                        9428157./3355443200.)*e16 + (4284333./1468006400.)*e14 + (6561./2621440.)*e12 + (
                        6561./3276800.)*e10,
            None,
            -3.,
            2.,
            2.
        ),
        (
            '-2*O - 2*n',
            '2*O + 2*n',
            '2, 2, 0, -4',
            -2*spin_freq - 2*orbital_freq,
            (330367./1103619686400.)*i12e8 + (-12107./2985984000.)*i10e10 + (-12107./4180377600.)*i10e8 + (
                        1857193./55738368000.)*i8e12 + (1957./66355200.)*i8e10 + (1957./92897280.)*i8e8 + (
                        -55591./326592000.)*i6e14 + (-21827./124416000.)*i6e12 + (-161./1036800.)*i6e10 + (
                        -23./207360.)*i6e8 + (248611./445906944.)*i4e16 + (26587./43545600.)*i4e14 + (
                        10439./16588800.)*i4e12 + (77./138240.)*i4e10 + (11./27648.)*i4e8 + (
                        -128441./119439360.)*i2e18 + (-22601./18579456.)*i2e16 + (-2417./1814400.)*i2e14 + (
                        -949./691200.)*i2e12 + (-7./5760.)*i2e10 + (-1./1152.)*i2e8 + (
                        3291434567./3511517184000.)*e20 + (128441./119439360.)*e18 + (22601./18579456.)*e16 + (
                        2417./1814400.)*e14 + (949./691200.)*e12 + (7./5760.)*e10 + (1./1152.)*e8,
            None,
            -2.,
            2.,
            2.
        ),
        (
            '-2*O - n',
            '2*O + n',
            '2, 2, 0, -3',
            -2*spin_freq - orbital_freq,
            (-27269./4564972339200.)*i14e6 + (330367./3210529996800.)*i12e8 + (330367./4414478745600.)*i12e6 + (
                        -7494233./7134511104000.)*i10e10 + (-133177./133772083200.)*i10e8 + (
                        -12107./16721510400.)*i10e6 + (122541469./17122826649600.)*i8e12 + (
                        1211383./158544691200.)*i8e10 + (21527./2972712960.)*i8e8 + (1957./371589120.)*i8e6 + (
                        -722174401./21403533312000.)*i6e14 + (-1440191./38220595200.)*i6e12 + (
                        -14237./353894400.)*i6e10 + (-253./6635520.)*i6e8 + (-23./829440.)*i6e6 + (
                        1621406413./15220290355200.)*i4e16 + (345387757./2853804441600.)*i4e14 + (
                        688787./5096079360.)*i4e12 + (6809./47185920.)*i4e10 + (121./884736.)*i4e8 + (
                        11./110592.)*i4e6 + (-167431204877./821895679180800.)*i2e18 + (
                        -147400583./634178764800.)*i2e16 + (-31398887./118908518400.)*i2e14 + (
                        -62617./212336640.)*i2e12 + (-619./1966080.)*i2e10 + (-11./36864.)*i2e8 + (-1./4608.)*i2e6 + (
                        20572630185001./115065395085312000.)*e20 + (167431204877./821895679180800.)*e18 + (
                        147400583./634178764800.)*e16 + (31398887./118908518400.)*e14 + (62617./212336640.)*e12 + (
                        619./1966080.)*e10 + (11./36864.)*e8 + (1./4608.)*e6,
            None,
            -1.,
            2.,
            2.
        ),
        (
            '-2*O + n',
            '2*O - n',
            '2, 2, 0, -1',
            -2*spin_freq + orbital_freq,
            (-561142037./51218989645824000.)*i18e2 + (-72518377./1339058552832000.)*i16e4 + (
                        72518377./334764638208000.)*i16e2 + (-27269./117050572800.)*i14e6 + (
                        27269./31701196800.)*i14e4 + (-27269./7925299200.)*i14e2 + (20152387./7063165992960.)*i12e8 + (
                        4294771./1471492915200.)*i12e6 + (-330367./30656102400.)*i12e4 + (330367./7664025600.)*i12e2 + (
                        -20981431./1189085184000.)*i10e10 + (-738527./26754416640.)*i10e8 + (
                        -157391./5573836800.)*i10e6 + (12107./116121600.)*i10e4 + (-12107./29030400.)*i10e2 + (
                        542537153./4756340736000.)*i8e12 + (3391481./26424115200.)*i8e10 + (119377./594542592.)*i8e8 + (
                        25441./123863040.)*i8e6 + (-1957./2580480.)*i8e4 + (1957./645120.)*i8e2 + (
                        -3801562889./7134511104000.)*i6e14 + (-6376267./10616832000.)*i6e12 + (
                        -39859./58982400.)*i6e10 + (-1403./1327104.)*i6e8 + (-299./276480.)*i6e6 + (23./5760.)*i6e4 + (
                        -23./1440.)*i6e2 + (543959483287./319626097459200.)*i4e16 + (
                        1818138773./951268147200.)*i4e14 + (3049519./1415577600.)*i4e12 + (19063./7864320.)*i4e10 + (
                        3355./884736.)*i4e8 + (143./36864.)*i4e6 + (-11./768.)*i4e4 + (11./192.)*i4e2 + (
                        -12842565048623./3835513169510400.)*i2e18 + (-49450862117./13317754060800.)*i2e16 + (
                        -165285343./39636172800.)*i2e14 + (-277229./58982400.)*i2e12 + (-1733./327680.)*i2e10 + (
                        -305./36864.)*i2e8 + (-13./1536.)*i2e6 + (1./32.)*i2e4 + (-1./8.)*i2e2 + (
                        2104490540764777./690392370511872000.)*e20 + (12842565048623./3835513169510400.)*e18 + (
                        49450862117./13317754060800.)*e16 + (165285343./39636172800.)*e14 + (277229./58982400.)*e12 + (
                        1733./327680.)*e10 + (305./36864.)*e8 + (13./1536.)*e6 + (-1./32.)*e4 + (1./8.)*e2,
            None,
            1.,
            2.,
            2.
        ),
        (
            '-2*O + 2*n',
            '2*O - 2*n',
            '2, 2, 0, 0',
            -2*spin_freq + 2*orbital_freq,
            (17616175987./9731608032706560000.)*i20 + (561142037./2560949482291200.)*i18e2 + (
                        -561142037./12804747411456000.)*i18 + (72518377./10627448832000.)*i16e4 + (
                        -72518377./16738231910400.)*i16e2 + (72518377./83691159552000.)*i16 + (
                        2154251./35663846400.)*i14e6 + (-27269./251596800.)*i14e4 + (27269./396264960.)*i14e2 + (
                        -27269./1981324800.)*i14 + (1221366799./4414478745600.)*i12e8 + (
                        -26098993./34488115200.)*i12e6 + (330367./243302400.)*i12e4 + (-330367./383201280.)*i12e2 + (
                        330367./1916006400.)*i12 + (133673387./139345920000.)*i10e10 + (
                        -44759579./16721510400.)*i10e8 + (956453./130636800.)*i10e6 + (-12107./921600.)*i10e4 + (
                        12107./1451520.)*i10e2 + (-12107./7257600.)*i10 + (67459747./44590694400.)*i8e12 + (
                        -21607237./3096576000.)*i8e10 + (7235029./371589120.)*i8e8 + (-154603./2903040.)*i8e6 + (
                        1957./20480.)*i8e4 + (-1957./32256.)*i8e2 + (1957./161280.)*i8 + (
                        40018597./12192768000.)*i6e14 + (-792833./99532800.)*i6e12 + (253943./6912000.)*i6e10 + (
                        -85031./829440.)*i6e8 + (1817./6480.)*i6e6 + (-161./320.)*i6e4 + (23./72.)*i6e2 + (
                        -23./360.)*i6 + (-6180779./5780275200.)*i4e16 + (-19139329./1625702400.)*i4e14 + (
                        379181./13271040.)*i4e12 + (-121451./921600.)*i4e10 + (40667./110592.)*i4e8 + (
                        -869./864.)*i4e6 + (231./128.)*i4e4 + (-55./48.)*i4e2 + (11./48.)*i4 + (
                        459927151./75246796800.)*i2e18 + (561889./240844800.)*i2e16 + (1739939./67737600.)*i2e14 + (
                        -34471./552960.)*i2e12 + (11041./38400.)*i2e10 + (-3697./4608.)*i2e8 + (79./36.)*i2e6 + (
                        -63./16.)*i2e4 + (5./2.)*i2e2 + (-1./2.)*i2 + (-475843828001./105345515520000.)*e20 + (
                        -459927151./75246796800.)*e18 + (-561889./240844800.)*e16 + (-1739939./67737600.)*e14 + (
                        34471./552960.)*e12 + (-11041./38400.)*e10 + (3697./4608.)*e8 + (-79./36.)*e6 + (63./16.)*e4 + (
                        -5./2.)*e2 + (1./2.),
            None,
            2.,
            2.,
            2.
        ),
        (
            '-2*O + 3*n',
            '2*O - 3*n',
            '2, 2, 0, 1',
            -2*spin_freq + 3*orbital_freq,
            (-561142037./1045285502976000.)*i18e2 + (-2973253457./63764692992000.)*i16e4 + (
                        72518377./6831931392000.)*i16e2 + (-7989817./6762921984.)*i14e6 + (
                        1118029./1509580800.)*i14e4 + (-27269./161740800.)*i14e2 + (
                        -14951419319./1307993702400.)*i12e8 + (96797531./6539968512.)*i12e6 + (
                        -13545047./1459814400.)*i12e4 + (2312569./1094860800.)*i12e2 + (
                        -7490613007./132120576000.)*i10e10 + (547926499./4954521600.)*i10e8 + (
                        -3547351./24772608.)*i10e6 + (496387./5529600.)*i10e4 + (-84749./4147200.)*i10e2 + (
                        -1490515781./8388608000.)*i8e12 + (1210797857./2936012800.)*i8e10 + (
                        -88567949./110100480.)*i8e8 + (2867005./2752512.)*i8e6 + (-80237./122880.)*i8e4 + (
                        13699./92160.)*i8e2 + (-12287201749./37748736000.)*i6e14 + (122622913./131072000.)*i6e12 + (
                        -14230123./6553600.)*i6e10 + (1040911./245760.)*i6e8 + (-33695./6144.)*i6e6 + (
                        6601./1920.)*i6e4 + (-1127./1440.)*i6e2 + (-14568597021./37580963840.)*i4e16 + (
                        5876487793./5033164800.)*i4e14 + (-175937223./52428800.)*i4e12 + (20417133./2621440.)*i4e10 + (
                        -497827./32768.)*i4e8 + (80575./4096.)*i4e6 + (-3157./256.)*i4e4 + (539./192.)*i4e2 + (
                        -6258845529./30064771072.)*i2e18 + (3973253733./4697620480.)*i2e16 + (
                        -534226163./209715200.)*i2e14 + (47982879./6553600.)*i2e12 + (-5568309./327680.)*i2e10 + (
                        135771./4096.)*i2e8 + (-21975./512.)*i2e6 + (861./32.)*i2e4 + (-49./8.)*i2e2 + (
                        -8388292638491./105226698752000.)*e20 + (6258845529./30064771072.)*e18 + (
                        -3973253733./4697620480.)*e16 + (534226163./209715200.)*e14 + (-47982879./6553600.)*e12 + (
                        5568309./327680.)*e10 + (-135771./4096.)*e8 + (21975./512.)*e6 + (-861./32.)*e4 + (49./8.)*e2,
            None,
            3.,
            2.,
            2.
        ),
        (
            '-2*O + 4*n',
            '2*O - 4*n',
            '2, 2, 0, 2',
            -2*spin_freq + 4*orbital_freq,
            (20957810953./334764638208000.)*i16e4 + (10662179./2377589760.)*i14e6 + (-7880741./7925299200.)*i14e4 + (
                        27602493217./275904921600.)*i12e8 + (-129173497./2299207680.)*i12e6 + (
                        95476063./7664025600.)*i12e4 + (1624868363./1741824000.)*i10e10 + (
                        -1011551957./1045094400.)*i10e8 + (4733837./8709120.)*i10e6 + (-3498923./29030400.)*i10e4 + (
                        16481327567./3715891200.)*i8e12 + (-262647013./38707200.)*i8e10 + (
                        163509307./23224320.)*i8e8 + (-765187./193536.)*i8e6 + (565573./645120.)*i8e4 + (
                        7974113179./653184000.)*i6e14 + (-193699813./8294400.)*i6e12 + (3086807./86400.)*i6e10 + (
                        -1921673./51840.)*i6e8 + (8993./432.)*i6e6 + (-6647./1440.)*i6e4 + (
                        1907075159./99532800.)*i4e16 + (-3813706303./87091200.)*i4e14 + (92639041./1105920.)*i4e12 + (
                        -1476299./11520.)*i4e10 + (919061./6912.)*i4e8 + (-21505./288.)*i4e6 + (3179./192.)*i4e4 + (
                        2949969133./182891520.)*i2e18 + (-173370469./4147200.)*i2e16 + (346700573./3628800.)*i2e14 + (
                        -8421731./46080.)*i2e12 + (134209./480.)*i2e10 + (-83551./288.)*i2e8 + (1955./12.)*i2e6 + (
                        -289./8.)*i2e4 + (38506861007131./7023034368000.)*e20 + (-2949969133./182891520.)*e18 + (
                        173370469./4147200.)*e16 + (-346700573./3628800.)*e14 + (8421731./46080.)*e12 + (
                        -134209./480.)*e10 + (83551./288.)*e8 + (-1955./12.)*e6 + (289./8.)*e4,
            None,
            4.,
            2.,
            2.
        ),
        (
            '-2*O + 5*n',
            '2*O - 5*n',
            '2, 2, 0, 3',
            -2*spin_freq + 5*orbital_freq,
            (-59909993./14046068736.)*i14e6 + (-363187309615./1412633198592.)*i12e8 + (
                        9435611887./176579149824.)*i12e6 + (-284381504605./57076088832.)*i10e10 + (
                        13309769915./5350883328.)*i10e8 + (-345788027./668860416.)*i10e6 + (
                        -28443427981175./684913065984.)*i8e12 + (229840011775./6341787648.)*i8e10 + (
                        -10757090825./594542592.)*i8e8 + (279469385./74317824.)*i8e6 + (
                        -29575311724925./171228266496.)*i6e14 + (334286583325./1528823808.)*i6e12 + (
                        -2701236725./14155776.)*i6e10 + (126424675./1327104.)*i6e8 + (-3284515./165888.)*i6e6 + (
                        -234023116288775./608811614208.)*i4e16 + (70723571516125./114152177664.)*i4e14 + (
                        -799380960125./1019215872.)*i4e12 + (6459479125./9437184.)*i4e10 + (
                        -302319875./884736.)*i4e8 + (7854275./110592.)*i4e6 + (
                        -2044426346565875./4696546738176.)*i2e18 + (21274828753525./25367150592.)*i2e16 + (
                        -6429415592375./4756340736.)*i2e14 + (72670996375./42467328.)*i2e12 + (
                        -587225375./393216.)*i2e10 + (27483625./36864.)*i2e8 + (-714025./4608.)*i2e6 + (
                        -180939012603859375./920523160682496.)*e20 + (2044426346565875./4696546738176.)*e18 + (
                        -21274828753525./25367150592.)*e16 + (6429415592375./4756340736.)*e14 + (
                        -72670996375./42467328.)*e12 + (587225375./393216.)*e10 + (-27483625./36864.)*e8 + (
                        714025./4608.)*e6,
            None,
            5.,
            2.,
            2.
        ),
        (
            '-2*O + 6*n',
            '2*O - 6*n',
            '2, 2, 0, 4',
            -2*spin_freq + 6*orbital_freq,
            (93853630663./490497638400.)*i12e8 + (29742019879./3096576000.)*i10e10 + (
                        -3439465523./1857945600.)*i10e8 + (430455400847./2752512000.)*i8e12 + (
                        -4807560329./68812800.)*i8e10 + (555962173./41287680.)*i8e8 + (23427200137./21504000.)*i6e14 + (
                        -5059005733./6144000.)*i6e12 + (56501731./153600.)*i6e10 + (-6534047./92160.)*i6e8 + (
                        661994783901./183500800.)*i4e16 + (-11204313109./2867200.)*i4e14 + (
                        2419524481./819200.)*i4e12 + (-27022567./20480.)*i4e10 + (3124979./12288.)*i4e8 + (
                        25998653133./4587520.)*i2e18 + (-180544031973./22937600.)*i2e16 + (
                        3055721757./358400.)*i2e14 + (-659870313./102400.)*i2e12 + (7369791./2560.)*i2e10 + (
                        -284089./512.)*i2e8 + (21810011108497./6422528000.)*e20 + (-25998653133./4587520.)*e18 + (
                        180544031973./22937600.)*e16 + (-3055721757./358400.)*e14 + (659870313./102400.)*e12 + (
                        -7369791./2560.)*e10 + (284089./512.)*e8,
            None,
            6.,
            2.,
            2.
        ),
        (
            '-2*O + 7*n',
            '2*O - 7*n',
            '2, 2, 0, 5',
            -2*spin_freq + 7*orbital_freq,
            (-90183922945109./15288238080000.)*i10e10 + (-7842210761371./32614907904.)*i8e12 + (
                        14577511952059./339738624000.)*i8e10 + (-9597992674380677./3057647616000.)*i6e14 + (
                        645169101983./509607936.)*i6e12 + (-1199274105407./5308416000.)*i6e10 + (
                        -109977996363287437./6522981580800.)*i4e16 + (4590344322529889./407686348800.)*i4e14 + (
                        -1542795678655./339738624.)*i4e12 + (573565876499./707788800.)*i4e10 + (
                        -1518065224694500853./39137889484800.)*i2e18 + (9997999669389767./271790899200.)*i2e16 + (
                        -417304029320899./16986931200.)*i2e14 + (140254152605./14155776.)*i2e12 + (
                        -52142352409./29491200.)*i2e10 + (-373553969127047054087./11741366845440000.)*e20 + (
                        1518065224694500853./39137889484800.)*e18 + (-9997999669389767./271790899200.)*e16 + (
                        417304029320899./16986931200.)*e14 + (-140254152605./14155776.)*e12 + (
                        52142352409./29491200.)*e10,
            None,
            7.,
            2.,
            2.
        ),
        (
            '-2*O + 8*n',
            '2*O - 8*n',
            '2, 2, 0, 6',
            -2*spin_freq + 8*orbital_freq,
            (10534550885077./83607552000.)*i8e12 + (1309029228049./326592000.)*i6e14 + (
                        -123809233703./186624000.)*i6e12 + (6365246831339./162570240.)*i4e16 + (
                        -626057456893./43545600.)*i4e14 + (59213111771./24883200.)*i4e12 + (
                        391340609035087./2743372800.)*i2e18 + (-578658802849./6773760.)*i2e16 + (
                        56914314263./1814400.)*i2e14 + (-5383010161./1036800.)*i2e12 + (
                        74152056390168773./438939648000.)*e20 + (-391340609035087./2743372800.)*e18 + (
                        578658802849./6773760.)*e16 + (-56914314263./1814400.)*e14 + (5383010161./1036800.)*e12,
            None,
            8.,
            2.,
            2.
        ),
        (
            '-2*O + 9*n',
            '2*O - 9*n',
            '2, 2, 0, 7',
            -2*spin_freq + 9*orbital_freq,
            (-376901937118463./205520896000.)*i6e14 + (-11237411230102629./263066746880.)*i4e16 + (
                        540772344561273./82208358400.)*i4e14 + (-415947859083950607./1503238553600.)*i2e18 + (
                        3064748517300717./32883343360.)*i2e16 + (-147483366698529./10276044800.)*i2e14 + (
                        -107205587887490347401./210453397504000.)*e20 + (415947859083950607./1503238553600.)*e18 + (
                        -3064748517300717./32883343360.)*e16 + (147483366698529./10276044800.)*e14,
            None,
            9.,
            2.,
            2.
        ),
        (
            '-2*O + 10*n',
            '2*O - 10*n',
            '2, 2, 0, 8',
            -2*spin_freq + 10*orbital_freq,
            (4422701659475./254803968.)*i4e16 + (176187983600875./668860416.)*i2e18 + (
                        -402063787225./10616832.)*i2e16 + (285521789807747375./337105649664.)*e20 + (
                        -176187983600875./668860416.)*e18 + (402063787225./10616832.)*e16,
            None,
            10.,
            2.,
            2.
        ),
        (
            '-2*O + 11*n',
            '2*O - 11*n',
            '2, 2, 0, 9',
            -2*spin_freq + 11*orbital_freq,
            (-6648821549377771726369./69039237051187200.)*i2e18 + (
                        -987278781529197450645517./1380784741023744000.)*e20 + (
                        6648821549377771726369./69039237051187200.)*e18,
            None,
            11.,
            2.,
            2.
        ),
        (
            '-2*O + 12*n',
            '2*O - 12*n',
            '2, 2, 0, 10',
            -2*spin_freq + 12*orbital_freq,
            (3816001995797209./16056320000.)*e20,
            None,
            12.,
            2.,
            2.
        ),
        (
            '-2*O - 8*n',
            '2*O + 8*n',
            '2, 2, 1, -8',
            -2*spin_freq - 8*orbital_freq,
            (31801945561./642252800.)*i4e16,
            None,
            -8.,
            0.,
            2.
        ),
        (
            '-2*O - 7*n',
            '2*O + 7*n',
            '2, 2, 1, -7',
            -2*spin_freq - 7*orbital_freq,
            (-186702632281./11324620800.)*i6e14 + (-1328179895713./40265318400.)*i4e16 + (
                        186702632281./7549747200.)*i4e14,
            None,
            -7.,
            0.,
            2.
        ),
        (
            '-2*O - 6*n',
            '2*O + 6*n',
            '2, 2, 1, -6',
            -2*spin_freq - 6*orbital_freq,
            (10029889./4096000.)*i8e12 + (2308743./358400.)*i6e14 + (-10029889./1228800.)*i6e12 + (
                        197537355./12845056.)*i4e16 + (-6926229./716800.)*i4e14 + (10029889./819200.)*i4e12,
            None,
            -6.,
            0.,
            2.
        ),
        (
            '-2*O - 5*n',
            '2*O + 5*n',
            '2, 2, 1, -5',
            -2*spin_freq - 5*orbital_freq,
            (-1979259./9175040.)*i10e10 + (-2947317./10485760.)*i8e12 + (3143529./2621440.)*i8e10 + (
                        -13524196093./3170893824.)*i6e14 + (982439./1048576.)*i6e12 + (-1047843./262144.)*i6e10 + (
                        28491858875./11274289152.)*i4e16 + (13524196093./2113929216.)*i4e14 + (
                        -2947317./2097152.)*i4e12 + (3143529./524288.)*i4e10,
            None,
            -5.,
            0.,
            2.
        ),
        (
            '-2*O - 4*n',
            '2*O + 4*n',
            '2, 2, 1, -4',
            -2*spin_freq - 4*orbital_freq,
            (26257./2073600.)*i12e8 + (-8041./230400.)*i10e10 + (-14399./138240.)*i10e8 + (708871./1024000.)*i8e12 + (
                        9933./51200.)*i8e10 + (5929./10240.)*i8e8 + (-3027367./1382400.)*i6e14 + (
                        -708871./307200.)*i6e12 + (-3311./5120.)*i6e10 + (-5929./3072.)*i6e8 + (
                        124914751./27525120.)*i4e16 + (3027367./921600.)*i4e14 + (708871./204800.)*i4e12 + (
                        9933./10240.)*i4e10 + (5929./2048.)*i4e8,
            None,
            -4.,
            0.,
            2.
        ),
        (
            '-2*O - 3*n',
            '2*O + 3*n',
            '2, 2, 1, -3',
            -2*spin_freq - 3*orbital_freq,
            (-2809./5322240.)*i14e6 + (215233./38707200.)*i12e8 + (87079./14515200.)*i12e6 + (
                        -4873237./58982400.)*i10e10 + (-118031./2580480.)*i10e8 + (-47753./967680.)*i10e6 + (
                        30355211./52428800.)*i8e12 + (6019881./13107200.)*i8e10 + (20829./81920.)*i8e8 + (
                        2809./10240.)*i8e6 + (-1055142483./419430400.)*i6e14 + (-30355211./15728640.)*i6e12 + (
                        -2006627./1310720.)*i6e10 + (-6943./8192.)*i6e8 + (-2809./3072.)*i6e6 + (
                        437242738059./93952409600.)*i4e16 + (3165427449./838860800.)*i4e14 + (
                        30355211./10485760.)*i4e12 + (6019881./2621440.)*i4e10 + (20829./16384.)*i4e8 + (
                        2809./2048.)*i4e6,
            None,
            -3.,
            0.,
            2.
        ),
        (
            '-2*O - 2*n',
            '2*O + 2*n',
            '2, 2, 1, -2',
            -2*spin_freq - 2*orbital_freq,
            (5461./336336000.)*i16e4 + (-1./2640.)*i14e6 + (-3./12320.)*i14e4 + (51491./7257600.)*i12e8 + (
                        31./7200.)*i12e6 + (31./11200.)*i12e4 + (-66623./806400.)*i10e10 + (-28237./483840.)*i10e8 + (
                        -17./480.)*i10e6 + (-51./2240.)*i10e4 + (505601./819200.)*i8e12 + (11757./25600.)*i8e10 + (
                        1661./5120.)*i8e8 + (63./320.)*i8e6 + (81./640.)*i8e4 + (-25565893./9676800.)*i6e14 + (
                        -505601./245760.)*i6e12 + (-3919./2560.)*i6e10 + (-1661./1536.)*i6e8 + (-21./32.)*i6e6 + (
                        -27./64.)*i6e4 + (678544541./137625600.)*i4e16 + (25565893./6451200.)*i4e14 + (
                        505601./163840.)*i4e12 + (11757./5120.)*i4e10 + (1661./1024.)*i4e8 + (63./64.)*i4e6 + (
                        81./128.)*i4e4,
            None,
            -2.,
            0.,
            2.
        ),
        (
            '-2*O - n',
            '2*O + n',
            '2, 2, 1, -1',
            -2*spin_freq - orbital_freq,
            (-257./681080400.)*i18e2 + (5461./336336000.)*i16e4 + (5461./756756000.)*i16e2 + (-17./39424.)*i14e6 + (
                        -3./12320.)*i14e4 + (-1./9240.)*i14e2 + (874541./116121600.)*i12e8 + (527./107520.)*i12e6 + (
                        31./11200.)*i12e4 + (31./25200.)*i12e2 + (-18103079./206438400.)*i10e10 + (
                        -479587./7741440.)*i10e8 + (-289./7168.)*i10e6 + (-51./2240.)*i10e4 + (-17./1680.)*i10e2 + (
                        85553819./131072000.)*i8e12 + (3194661./6553600.)*i8e10 + (28211./81920.)*i8e8 + (
                        459./2048.)*i8e6 + (81./640.)*i8e4 + (9./160.)*i8e2 + (-24654653741./8808038400.)*i6e14 + (
                        -85553819./39321600.)*i6e12 + (-1064887./655360.)*i6e10 + (-28211./24576.)*i6e8 + (
                        -765./1024.)*i6e6 + (-27./64.)*i6e4 + (-3./16.)*i6e2 + (689312857627./131533373440.)*i4e16 + (
                        24654653741./5872025600.)*i4e14 + (85553819./26214400.)*i4e12 + (3194661./1310720.)*i4e10 + (
                        28211./16384.)*i4e8 + (2295./2048.)*i4e6 + (81./128.)*i4e4 + (9./32.)*i4e2,
            None,
            -1.,
            0.,
            2.
        ),
        (
            '-2*O',
            '2*O',
            '2, 2, 1, 0',
            -2*spin_freq,
            (73./10337827500.)*i20 + (-257./510810300.)*i18e2 + (-257./1532430900.)*i18 + (5461./283783500.)*i16e4 + (
                        5461./567567000.)*i16e2 + (5461./1702701000.)*i16 + (-1./2079.)*i14e6 + (-1./3465.)*i14e4 + (
                        -1./6930.)*i14e2 + (-1./20790.)*i14 + (31./3780.)*i12e8 + (31./5670.)*i12e6 + (
                        31./9450.)*i12e4 + (31./18900.)*i12e2 + (31./56700.)*i12 + (-17./180.)*i10e10 + (
                        -17./252.)*i10e8 + (-17./378.)*i10e6 + (-17./630.)*i10e4 + (-17./1260.)*i10e2 + (
                        -17./3780.)*i10 + (7./10.)*i8e12 + (21./40.)*i8e10 + (3./8.)*i8e8 + (1./4.)*i8e6 + (
                        3./20.)*i8e4 + (3./40.)*i8e2 + (1./40.)*i8 + -3.*i6e14 + (-7./3.)*i6e12 + (-7./4.)*i6e10 + (
                        -5./4.)*i6e8 + (-5./6.)*i6e6 + (-1./2.)*i6e4 + (-1./4.)*i6e2 + (-1./12.)*i6 + (45./8.)*i4e16 + (
                        9./2.)*i4e14 + (7./2.)*i4e12 + (21./8.)*i4e10 + (15./8.)*i4e8 + (5./4.)*i4e6 + (3./4.)*i4e4 + (
                        3./8.)*i4e2 + (1./8.)*i4,
            None,
            0.,
            0.,
            2.
        ),
        (
            '-2*O + n',
            '2*O - n',
            '2, 2, 1, 1',
            -2*spin_freq + orbital_freq,
            (-257./681080400.)*i18e2 + (5461./336336000.)*i16e4 + (5461./756756000.)*i16e2 + (-17./39424.)*i14e6 + (
                        -3./12320.)*i14e4 + (-1./9240.)*i14e2 + (874541./116121600.)*i12e8 + (527./107520.)*i12e6 + (
                        31./11200.)*i12e4 + (31./25200.)*i12e2 + (-18103079./206438400.)*i10e10 + (
                        -479587./7741440.)*i10e8 + (-289./7168.)*i10e6 + (-51./2240.)*i10e4 + (-17./1680.)*i10e2 + (
                        85553819./131072000.)*i8e12 + (3194661./6553600.)*i8e10 + (28211./81920.)*i8e8 + (
                        459./2048.)*i8e6 + (81./640.)*i8e4 + (9./160.)*i8e2 + (-24654653741./8808038400.)*i6e14 + (
                        -85553819./39321600.)*i6e12 + (-1064887./655360.)*i6e10 + (-28211./24576.)*i6e8 + (
                        -765./1024.)*i6e6 + (-27./64.)*i6e4 + (-3./16.)*i6e2 + (689312857627./131533373440.)*i4e16 + (
                        24654653741./5872025600.)*i4e14 + (85553819./26214400.)*i4e12 + (3194661./1310720.)*i4e10 + (
                        28211./16384.)*i4e8 + (2295./2048.)*i4e6 + (81./128.)*i4e4 + (9./32.)*i4e2,
            None,
            1.,
            0.,
            2.
        ),
        (
            '-2*O + 2*n',
            '2*O - 2*n',
            '2, 2, 1, 2',
            -2*spin_freq + 2*orbital_freq,
            (5461./336336000.)*i16e4 + (-1./2640.)*i14e6 + (-3./12320.)*i14e4 + (51491./7257600.)*i12e8 + (
                        31./7200.)*i12e6 + (31./11200.)*i12e4 + (-66623./806400.)*i10e10 + (-28237./483840.)*i10e8 + (
                        -17./480.)*i10e6 + (-51./2240.)*i10e4 + (505601./819200.)*i8e12 + (11757./25600.)*i8e10 + (
                        1661./5120.)*i8e8 + (63./320.)*i8e6 + (81./640.)*i8e4 + (-25565893./9676800.)*i6e14 + (
                        -505601./245760.)*i6e12 + (-3919./2560.)*i6e10 + (-1661./1536.)*i6e8 + (-21./32.)*i6e6 + (
                        -27./64.)*i6e4 + (678544541./137625600.)*i4e16 + (25565893./6451200.)*i4e14 + (
                        505601./163840.)*i4e12 + (11757./5120.)*i4e10 + (1661./1024.)*i4e8 + (63./64.)*i4e6 + (
                        81./128.)*i4e4,
            None,
            2.,
            0.,
            2.
        ),
        (
            '-2*O + 3*n',
            '2*O - 3*n',
            '2, 2, 1, 3',
            -2*spin_freq + 3*orbital_freq,
            (-2809./5322240.)*i14e6 + (215233./38707200.)*i12e8 + (87079./14515200.)*i12e6 + (
                        -4873237./58982400.)*i10e10 + (-118031./2580480.)*i10e8 + (-47753./967680.)*i10e6 + (
                        30355211./52428800.)*i8e12 + (6019881./13107200.)*i8e10 + (20829./81920.)*i8e8 + (
                        2809./10240.)*i8e6 + (-1055142483./419430400.)*i6e14 + (-30355211./15728640.)*i6e12 + (
                        -2006627./1310720.)*i6e10 + (-6943./8192.)*i6e8 + (-2809./3072.)*i6e6 + (
                        437242738059./93952409600.)*i4e16 + (3165427449./838860800.)*i4e14 + (
                        30355211./10485760.)*i4e12 + (6019881./2621440.)*i4e10 + (20829./16384.)*i4e8 + (
                        2809./2048.)*i4e6,
            None,
            3.,
            0.,
            2.
        ),
        (
            '-2*O + 4*n',
            '2*O - 4*n',
            '2, 2, 1, 4',
            -2*spin_freq + 4*orbital_freq,
            (26257./2073600.)*i12e8 + (-8041./230400.)*i10e10 + (-14399./138240.)*i10e8 + (708871./1024000.)*i8e12 + (
                        9933./51200.)*i8e10 + (5929./10240.)*i8e8 + (-3027367./1382400.)*i6e14 + (
                        -708871./307200.)*i6e12 + (-3311./5120.)*i6e10 + (-5929./3072.)*i6e8 + (
                        124914751./27525120.)*i4e16 + (3027367./921600.)*i4e14 + (708871./204800.)*i4e12 + (
                        9933./10240.)*i4e10 + (5929./2048.)*i4e8,
            None,
            4.,
            0.,
            2.
        ),
        (
            '-2*O + 5*n',
            '2*O - 5*n',
            '2, 2, 1, 5',
            -2*spin_freq + 5*orbital_freq,
            (-1979259./9175040.)*i10e10 + (-2947317./10485760.)*i8e12 + (3143529./2621440.)*i8e10 + (
                        -13524196093./3170893824.)*i6e14 + (982439./1048576.)*i6e12 + (-1047843./262144.)*i6e10 + (
                        28491858875./11274289152.)*i4e16 + (13524196093./2113929216.)*i4e14 + (
                        -2947317./2097152.)*i4e12 + (3143529./524288.)*i4e10,
            None,
            5.,
            0.,
            2.
        ),
        (
            '-2*O + 6*n',
            '2*O - 6*n',
            '2, 2, 1, 6',
            -2*spin_freq + 6*orbital_freq,
            (10029889./4096000.)*i8e12 + (2308743./358400.)*i6e14 + (-10029889./1228800.)*i6e12 + (
                        197537355./12845056.)*i4e16 + (-6926229./716800.)*i4e14 + (10029889./819200.)*i4e12,
            None,
            6.,
            0.,
            2.
        ),
        (
            '-2*O + 7*n',
            '2*O - 7*n',
            '2, 2, 1, 7',
            -2*spin_freq + 7*orbital_freq,
            (-186702632281./11324620800.)*i6e14 + (-1328179895713./40265318400.)*i4e16 + (
                        186702632281./7549747200.)*i4e14,
            None,
            7.,
            0.,
            2.
        ),
        (
            '-2*O + 8*n',
            '2*O - 8*n',
            '2, 2, 1, 8',
            -2*spin_freq + 8*orbital_freq,
            (31801945561./642252800.)*i4e16,
            None,
            8.,
            0.,
            2.
        ),
        (
            '-2*O - 8*n',
            '2*O + 8*n',
            '2, 2, 2, -6',
            -2*spin_freq - 8*orbital_freq,
            (5383010161./265420800.)*i8e12,
            None,
            -8.,
            -2.,
            2.
        ),
        (
            '-2*O - 7*n',
            '2*O + 7*n',
            '2, 2, 2, -5',
            -2*spin_freq - 7*orbital_freq,
            (-52142352409./22649241600.)*i10e10 + (-140254152605./3623878656.)*i8e12 + (52142352409./7549747200.)*i8e10,
            None,
            -7.,
            -2.,
            2.
        ),
        (
            '-2*O - 6*n',
            '2*O + 6*n',
            '2, 2, 2, -4',
            -2*spin_freq - 6*orbital_freq,
            (5397691./47185920.)*i12e8 + (2456597./655360.)*i10e10 + (-284089./393216.)*i10e8 + (
                        659870313./26214400.)*i8e12 + (-7369791./655360.)*i8e10 + (284089./131072.)*i8e8,
            None,
            -6.,
            -2.,
            2.
        ),
        (
            '-2*O - 5*n',
            '2*O + 5*n',
            '2, 2, 2, -3',
            -2*spin_freq - 5*orbital_freq,
            (-714025./222953472.)*i14e6 + (-104437775./679477248.)*i12e8 + (2713295./84934656.)*i12e6 + (
                        -587225375./301989888.)*i10e10 + (27483625./28311552.)*i10e8 + (-714025./3538944.)*i10e6 + (
                        -72670996375./10871635968.)*i8e12 + (587225375./100663296.)*i8e10 + (
                        -27483625./9437184.)*i8e8 + (714025./1179648.)*i8e6,
            None,
            -5.,
            -2.,
            2.
        ),
        (
            '-2*O - 4*n',
            '2*O + 4*n',
            '2, 2, 2, -2',
            -2*spin_freq - 4*orbital_freq,
            (132073./2477260800.)*i16e4 + (1955./580608.)*i14e6 + (-289./387072.)*i14e4 + (1587469./26542080.)*i12e8 + (
                        -7429./221184.)*i12e6 + (5491./737280.)*i12e4 + (134209./368640.)*i10e10 + (
                        -83551./221184.)*i10e8 + (1955./9216.)*i10e6 + (-289./6144.)*i10e4 + (
                        8421731./11796480.)*i8e12 + (-134209./122880.)*i8e10 + (83551./73728.)*i8e8 + (
                        -1955./3072.)*i8e6 + (289./2048.)*i8e4,
            None,
            -4.,
            -2.,
            2.
        ),
        (
            '-2*O - 3*n',
            '2*O + 3*n',
            '2, 2, 2, -1',
            -2*spin_freq - 3*orbital_freq,
            (-3437./7007109120.)*i18e2 + (-18737./471859200.)*i16e4 + (3199./353894400.)*i16e2 + (
                        -7325./8257536.)*i14e6 + (41./73728.)*i14e4 + (-7./55296.)*i14e2 + (
                        -859883./125829120.)*i12e8 + (27835./3145728.)*i12e6 + (-5453./983040.)*i12e4 + (
                        931./737280.)*i12e2 + (-1856103./83886080.)*i10e10 + (45257./1048576.)*i10e8 + (
                        -7325./131072.)*i10e6 + (287./8192.)*i10e4 + (-49./6144.)*i10e2 + (
                        -47982879./1677721600.)*i8e12 + (5568309./83886080.)*i8e10 + (-135771./1048576.)*i8e8 + (
                        21975./131072.)*i8e6 + (-861./8192.)*i8e4 + (49./2048.)*i8e2,
            None,
            -3.,
            -2.,
            2.
        ),
        (
            '-2*O - 2*n',
            '2*O + 2*n',
            '2, 2, 2, 0',
            -2*spin_freq - 2*orbital_freq,
            (164573./95647039488000.)*i20 + (491./2452488192.)*i18e2 + (-491./12262440960.)*i18 + (
                        457./78643200.)*i16e4 + (-457./123863040.)*i16e2 + (457./619315200.)*i16 + (
                        19./435456.)*i14e6 + (-1./12288.)*i14e4 + (5./96768.)*i14e2 + (-1./96768.)*i14 + (
                        7847./84934656.)*i12e8 + (-361./829440.)*i12e6 + (133./163840.)*i12e4 + (-19./36864.)*i12e2 + (
                        19./184320.)*i12 + (-4079./29491200.)*i10e10 + (-2065./3538944.)*i10e8 + (19./6912.)*i10e6 + (
                        -21./4096.)*i10e4 + (5./1536.)*i10e2 + (-1./1536.)*i10 + (-2323./47185920.)*i8e12 + (
                        4079./9830400.)*i8e10 + (2065./1179648.)*i8e8 + (-19./2304.)*i8e6 + (63./4096.)*i8e4 + (
                        -5./512.)*i8e2 + (1./512.)*i8,
            None,
            -2.,
            -2.,
            2.
        ),
        (
            '-2*O - n',
            '2*O + n',
            '2, 2, 2, 1',
            -2*spin_freq - orbital_freq,
            (-491./49049763840.)*i18e2 + (-457./9909043200.)*i16e4 + (457./2477260800.)*i16e2 + (
                        -13./74317824.)*i14e6 + (1./1548288.)*i14e4 + (-1./387072.)*i14e2 + (1159./679477248.)*i12e8 + (
                        247./141557760.)*i12e6 + (-19./2949120.)*i12e4 + (19./737280.)*i12e2 + (
                        -1733./251658240.)*i10e10 + (-305./28311552.)*i10e8 + (-13./1179648.)*i10e6 + (
                        1./24576.)*i10e4 + (-1./6144.)*i10e2 + (277229./15099494400.)*i8e12 + (
                        1733./83886080.)*i8e10 + (305./9437184.)*i8e8 + (13./393216.)*i8e6 + (-1./8192.)*i8e4 + (
                        1./2048.)*i8e2,
            None,
            -1.,
            -2.,
            2.
        ),
        (
            '-2*O + n',
            '2*O - n',
            '2, 2, 2, 3',
            -2*spin_freq + orbital_freq,
            (-1./222953472.)*i14e6 + (209./3397386240.)*i12e8 + (19./424673280.)*i12e6 + (-619./1509949440.)*i10e10 + (
                        -11./28311552.)*i10e8 + (-1./3538944.)*i10e6 + (62617./54358179840.)*i8e12 + (
                        619./503316480.)*i8e10 + (11./9437184.)*i8e8 + (1./1179648.)*i8e6,
            None,
            1.,
            -2.,
            2.
        ),
        (
            '-2*O + 2*n',
            '2*O - 2*n',
            '2, 2, 2, 4',
            -2*spin_freq + 2*orbital_freq,
            (19./106168320.)*i12e8 + (-7./4423680.)*i10e10 + (-1./884736.)*i10e8 + (949./176947200.)*i8e12 + (
                        7./1474560.)*i8e10 + (1./294912.)*i8e8,
            None,
            2.,
            -2.,
            2.
        ),
        (
            '-2*O + 3*n',
            '2*O - 3*n',
            '2, 2, 2, 5',
            -2*spin_freq + 3*orbital_freq,
            (-2187./838860800.)*i10e10 + (6561./671088640.)*i8e12 + (6561./838860800.)*i8e10,
            None,
            3.,
            -2.,
            2.
        ),
        (
            '-2*O + 4*n',
            '2*O - 4*n',
            '2, 2, 2, 6',
            -2*spin_freq + 4*orbital_freq,
            (1./64800.)*i8e12,
            None,
            4.,
            -2.,
            2.
        )
    )

    return mode_data_output