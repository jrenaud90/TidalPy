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