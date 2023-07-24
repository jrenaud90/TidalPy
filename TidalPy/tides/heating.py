""" Functions used to estimate tidal heating within a solid planet or layer. """

import numpy as np

from ..utilities.performance import njit

@njit(cacheable=True)
def calculate_volumetric_heating(stress: np.ndarray, strain: np.ndarray) -> np.ndarray:
    """ Calculates the tidal heating rate per unit volume based on the tidal stresses and strains.

    Parameters
    ----------
    stress : np.ndarray
        Tidal stress tensor (complex np.ndarray) [Pa]
    strain : np.ndarray
        Tidal strain tensor (complex np.ndarray) [unitless]

    Returns
    -------
    volumetric_heating : np.ndarray
        Tidal heating rate per unit volume [W m-3]
    """

    # Find real and imaginary components of the stress and strain tensor.
    stress_real = np.real(stress)
    stress_imag = np.imag(stress)
    strain_real = np.real(strain)
    strain_imag = np.imag(strain)

    # Find heating rate per unit volume.
    # Heating is equal to imag[o] * real[s] - real[o] * imag[s] but we need to multiply by two for the cross terms
    #  since it is part of a symmetric matrix but only one side of the matrix is calculated.
    volumetric_heating = (
        # Im[s_rr] Re[e_rr] - Re[s_rr] Im[e_rr]
        stress_imag[0] * strain_real[0] - stress_real[0] * strain_imag[0] +
        # Im[s_thth] Re[e_thth] - Re[s_thth] Im[e_thth]
        stress_imag[1] * strain_real[1] - stress_real[1] * strain_imag[1] +
        # Im[s_phiphi] Re[e_phiphi] - Re[s_phiphi] Im[e_phiphi]
        stress_imag[2] * strain_real[2] - stress_real[2] * strain_imag[2] +
        # 2 Im[s_rth] Re[e_rth] - Re[s_rth] Im[e_rth]
        2. * (stress_imag[3] * strain_real[3] - stress_real[3] * strain_imag[3]) +
        # 2 Im[s_rphi] Re[e_rphi] - Re[s_rphi] Im[e_rphi]
        2. * (stress_imag[4] * strain_real[4] - stress_real[4] * strain_imag[4]) +
        # 2 Im[s_thphi] Re[e_thphi] - Re[s_thphi] Im[e_thphi]
        2. * (stress_imag[5] * strain_real[5] - stress_real[5] * strain_imag[5])
    )

    # TODO: Without this abs term the resulting heating maps are very blotchy around
    #    Europa book does have an abs at Equation 42, Page 102
    volumetric_heating = np.abs(volumetric_heating)

    return volumetric_heating
