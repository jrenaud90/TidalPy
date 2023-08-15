""" Calculation displacement based on tidal solutions and tidal potential.

References
----------
SVC16 : Sabadini, Vermeerson, & Cambiotti (2016, DOI: 10.1007/978-94-017-7552-6)
HH14  : Henning & Hurford (2014, DOI: 10.1088/0004-637X/789/1/30)
TB05  : Tobie et al. (2005, DOI: 10.1016/j.icarus.2005.04.006)
B13   : Beuthe (2013, DOI: 10.1016/j.icarus.2012.11.020)
"""

from typing import Tuple, TYPE_CHECKING

import numpy as np

from TidalPy.utilities.performance import njit

if TYPE_CHECKING:
    from TidalPy.utilities.types import FloatArray


@njit(cacheable=True)
def calculate_displacements(
    tidal_potential: 'FloatArray',
    tidal_potential_partial_theta: 'FloatArray', tidal_potential_partial_phi: 'FloatArray',
    tidal_solution_y: np.ndarray, colatitude: 'FloatArray'
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """ Calculate tidal displacements using the tidal potential and its partial derivatives as well as the y-solution
    vector.

    Parameters
    ----------
    tidal_potential : FloatArray
        The tidal potential
    tidal_potential_partial_theta : FloatArray
        The partial derivative of the tidal potential with respect to the colatitude
    tidal_potential_partial_phi : FloatArray
        The partial derivative of the tidal potential with respect to the longitude
    tidal_solution_y : np.ndarray
        Vector of tidal solutions often denoted as a lower-case "y"
    colatitude : FloatArray
        Co-latitude where the calculations are conducted

    Returns
    -------
    radial_displacement : FloatArray
        The tidal displacement in the radial direction (r)
    polar_displacement : FloatArray
        The tidal displacement in the polar direction (theta)
    azimuthal_displacement : FloatArray
        The tidal displacement in the azimuthal direction (phi)
    """

    # Shortcuts
    radius_n = tidal_solution_y.shape[1]
    y1 = tidal_solution_y[0, :]
    y3 = tidal_solution_y[2, :]
    shape = (radius_n, *tidal_potential.shape)
    potential_dphi_over_sin = tidal_potential_partial_phi / np.sin(colatitude)

    # Build output arrays
    radial_displacement = np.empty(shape, dtype=np.complex128)
    polar_displacement = np.empty(shape, dtype=np.complex128)
    azimuthal_displacement = np.empty(shape, dtype=np.complex128)

    for ri in range(radius_n):

        # Pull out radius dependent values
        y1_ = y1[ri]
        y3_ = y3[ri]

        # Calculate displacement
        radial_displacement[ri, :, :, :] = y1_ * tidal_potential
        polar_displacement[ri, :, :, :] = y3_ * tidal_potential_partial_theta
        azimuthal_displacement[ri, :, :, :] = y3_ * potential_dphi_over_sin

    return radial_displacement, polar_displacement, azimuthal_displacement
