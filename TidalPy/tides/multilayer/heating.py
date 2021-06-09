""" Functions related to the calculation of tidal heating using method described in TB05

References
----------
SVC16 : Sabadini, Vermeerson, & Cambiotti (2016, DOI: 10.1007/978-94-017-7552-6)
HH14  : Henning & Hurford (2014, DOI: 10.1088/0004-637X/789/1/30)
TB05  : Tobie et al. (2005, DOI: 10.1016/j.icarus.2005.04.006)
ID    : IcyDwarf Code by Marc Neveu (https://github.com/MarcNeveu/IcyDwarf/blob/master/IcyDwarf/Thermal.h)
"""

import numpy as np

from ...utilities.types import FloatArray
from ...utilities.performance import njit
from ...constants import G

@njit(cacheable=True)
def calc_radial_tidal_heating(eccentricity: FloatArray, orbital_frequency: FloatArray, semi_major_axis: FloatArray,
                              tidal_host_mass: float,
                              radius_array: np.ndarray, radial_sensitivity_to_shear: np.ndarray,
                              complex_shear_modulus: np.ndarray, order_l: int = 2):
    """ Calculate tidal heating as a function of radius.

    Assumptions
    -----------
    - World is in synchronous rotation with no obliquity.
    - Eccentricity is low enough to only include the e^2 truncation.

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital eccentricity of the world.
    orbital_frequency : FloatArray
        Orbital mean motion of the world [rad s-1].
    semi_major_axis : FloatArray
        Orbital separation of the world and the tidal host [m].
    tidal_host_mass : float
        Mass of the tidal host [kg]
    radius_array : np.ndarray
        Radius array throughout the target world [m].
    radial_sensitivity_to_shear : np.ndarray
        Radial displacement's sensitivity of to shear stress.
        See definition and calculation in decompose.py
    complex_shear_modulus : np.ndarray
        Complex shear modulus at each radius as calculated by the world's layer rheology [Pa]
    order_l : int = 2
        Tidal harmonic order, l.

    Returns
    -------
    radial_tidal_heating : np.ndarray
        Volumetric tidal heating as a function of radius [W m-3]

    """

    world_radius = radius_array[-1]

    # TODO: This term appears to contain all the information that could be upgraded in the future to eliminate some or
    #    all of the assumptions mentioned in the doc strings.
    portion_to_be_upgraded = (7. * eccentricity**2 * orbital_frequency)

    # The below is modified from TB05. That reference assumed l=2 and that the tidal host's mass was much larger than
    #    the target planet's mass. We have removed these restrictions in the below. The derivation was done by changing
    #    Eq. 2 in TB05 to the general form and then proceeding with their derivations (Eqs. 35--37) except that we did
    #    not set l=2.

    # To keep things consistent with other parts of TidalPy and the references that it is built off of, we are defining
    #    a few extra terms here. Some of their components are redundant and will be divided out shortly. We are choosing
    #    consistency at a (very) slight performance hit.
    tidal_susceptibility = (3. / 2.) * G * tidal_host_mass**2 * world_radius**5 / semi_major_axis**6

    radial_tidal_heating = (tidal_susceptibility / world_radius) * \
        (G * radial_sensitivity_to_shear * np.imag(complex_shear_modulus) / ((2. * order_l + 1.) * radius_array**2)) *\
        portion_to_be_upgraded

    # TODO: Finding that the heating rate as a function of depth can go negative. Setting those to zero for now.
    radial_tidal_heating[radial_tidal_heating < 0.] = 0.

    return radial_tidal_heating

