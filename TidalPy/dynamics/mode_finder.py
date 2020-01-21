import numpy as np

from typing import Callable

from ..performance import njit
from .nsr_modes_l2_NEW import OutputType as NsrOutputType
from .sync_modes_l2_NEW import OutputType as SyncOutputType
from . import MODE_ZERO_TOL

NsrModeFuncType = Callable[[np.ndarray, np.ndarray, np.ndarray, np.ndarray], NsrOutputType]
SyncModeFuncType = Callable[[np.ndarray, np.ndarray, np.ndarray, np.ndarray], SyncOutputType]


@njit
def nsr_mode_finder(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray,
                    nsr_mode_function: NsrModeFuncType):
    """ Find the coefficients for NSR tidal dissipation using a provide mode function.

    The mode function sets the limits on eccentricity and inclination truncations.

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
    nsr_mode_function : NsrModeFuncType
        Mode function from the nsr_modes_l2 or nsr_modes_l3 packages.

    Returns
    -------
    mode_names : List[str]
        Mode symbols. aka '2n - 2o'
    freq_names : List[str]
        Mode's frequency symbols (abs of mode). aka '2n' for the tidal mode '-2n'
    modes : List[np.ndarray]
        List of calculated mode values [rad s-1]
    freqs : List[np.ndarray]
        List of calculated frequency values (abs of mode) [rad s-1]
    heating_subterms : List[np.ndarray]
        List of reduced tidal heating subterms [rad s-1]
    ztorque_subterms : List[np.ndarray]
        List of reduced tidal ztorque subterms [unitless] for i=0, [rads] for i!=0
    dudm_subterms : List[np.ndarray]
        List of reduced derivatives of tidal potential wrt mean anomaly [unitless] for i=0, [rads] for i!=0
    dudw_subterms : List[np.ndarray]
        List of reduced derivatives of tidal potential wrt periapsis [unitless] for i=0, [rads] for i!=0
    dudo_subterms : List[np.ndarray]
        List of reduced derivatives of tidal potential wrt node [unitless] for i=0, [rads] for i!=0

    """

    # Find data for each tidal mode (main calculation step)
    mode_datas = nsr_mode_function(orbital_freq, spin_freq, eccentricity, inclination)

    # Store results for subsequent calculations
    heating_subterms = list()
    ztorque_subterms = list()
    dudm_subterms = list()
    dudw_subterms = list()
    dudo_subterms = list()
    freqs = list()
    modes = list()
    freq_names = list()
    mode_names = list()
    for mode_name, freq_name, lmpq_str, mode, ei_coeff, _, dudm_coeff, dudw_coeff, dudo_coeff in mode_datas:
        sgn = np.sign(mode)
        freq = np.abs(mode)

        # Calculate subterm for each parameter
        heating_subterm = ei_coeff * freq
        ztorque_subterm = dudo_coeff * sgn * ei_coeff
        dudm_subterm = dudm_coeff * sgn * ei_coeff
        dudw_subterm = dudw_coeff * sgn * ei_coeff
        dudo_subterm = dudo_coeff * sgn * ei_coeff

        # Check for very small frequencies (Im[k](0) == 0 --> all subterms at this freq == 0)
        freq_too_low_index = freq < MODE_ZERO_TOL
        if np.any(freq_too_low_index):
            freq[freq_too_low_index] = 0.
            mode[freq_too_low_index] = 0.
            heating_subterm[freq_too_low_index] = 0.
            ztorque_subterm[freq_too_low_index] = 0.
            dudm_subterm[freq_too_low_index] = 0.
            dudw_subterm[freq_too_low_index] = 0.
            dudo_subterm[freq_too_low_index] = 0.

        # Save results
        freqs.append(freq)
        modes.append(mode)
        freq_names.append(freq_name)
        mode_names.append(mode_name)
        heating_subterms.append(heating_subterm)
        ztorque_subterms.append(ztorque_subterm)
        dudm_subterms.append(dudm_subterm)
        dudw_subterms.append(dudw_subterm)
        dudo_subterms.append(dudo_subterm)

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return mode_names, freq_names, modes, freqs, \
           heating_subterms, ztorque_subterms, dudm_subterms, dudw_subterms, dudo_subterms


@njit
def sync_mode_finder(orbital_freq: np.ndarray, spin_freq: np.ndarray, eccentricity: np.ndarray, inclination: np.ndarray,
                     sync_mode_function: SyncModeFuncType):
    """ Find the coefficients for spin-sync tidal dissipation using a provide mode function.

    The mode function sets the limits on eccentricity and inclination truncations.

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
    sync_mode_function : SyncModeFuncType
        Mode function from the sync_modes_l2 or sync_modes_l3 packages.

    Returns
    -------

    """

    # Find data for each tidal mode (main calculation step)
    mode_datas = sync_mode_function(orbital_freq, spin_freq, eccentricity, inclination)

    # Store results for subsequent calculations
    heating_subterms = list()
    ztorque_subterms = list()
    dudm_subterms = list()
    dudw_subterms = list()
    dudo_subterms = list()
    freqs = list()
    modes = list()
    freq_names = list()
    mode_names = list()
    for mode_name, freq_name, lmpq_str, mode, _, heating_subterm, dudm_subterm, dudw_subterm, dudo_subterm in mode_datas:
        freq = np.abs(mode)

        ztorque_subterm = np.copy(dudo_subterm)

        # Check for very small frequencies (Im[k](0) == 0 --> all subterms at this freq == 0)
        freq_too_low_index = freq < MODE_ZERO_TOL
        if np.any(freq_too_low_index):
            freq[freq_too_low_index] = 0.
            mode[freq_too_low_index] = 0.
            heating_subterm[freq_too_low_index] = 0.
            ztorque_subterm[freq_too_low_index] = 0.
            dudm_subterm[freq_too_low_index] = 0.
            dudw_subterm[freq_too_low_index] = 0.
            dudo_subterm[freq_too_low_index] = 0.

        # Save results
        freqs.append(freq)
        modes.append(mode)
        freq_names.append(freq_name)
        mode_names.append(mode_name)
        heating_subterms.append(heating_subterm)
        ztorque_subterms.append(ztorque_subterm)
        dudm_subterms.append(dudm_subterm)
        dudw_subterms.append(dudw_subterm)
        dudo_subterms.append(dudo_subterm)

    # As of Numba June 2019, tuple(List) is not supported. If that support comes then these should be uncommented.
    # modes = tuple(modes)
    # freqs = tuple(freqs)

    return mode_names, modes, freqs, heating_subterms, ztorque_subterms, dudm_subterms, dudw_subterms, dudo_subterms