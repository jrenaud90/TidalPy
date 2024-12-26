# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

""" Functions related to the calculation of tidal heating using method described in TB05

References
----------
SVC16 : Sabadini, Vermeerson, & Cambiotti (2016, DOI: 10.1007/978-94-017-7552-6)
HH14  : Henning & Hurford (2014, DOI: 10.1088/0004-637X/789/1/30)
TB05  : Tobie et al. (2005, DOI: 10.1016/j.icarus.2005.04.006)
ID    : IcyDwarf Code by Marc Neveu (https://github.com/MarcNeveu/IcyDwarf/blob/master/IcyDwarf/Thermal.h)
"""

from typing import TYPE_CHECKING

import numpy as np

from TidalPy.constants cimport d_G

cdef void cf_calc_radial_tidal_heating(
        double* tidal_heating_arr_ptr,
        size_t total_slices,
        double eccentricity,
        double orbital_frequency,
        double semi_major_axis,
        double host_mass,
        double* radius_arr_ptr,
        double* radial_sensitivity_to_shear_arr_ptr,
        double complex* complex_shear_modulus_arr_ptr,
        double* radial_sensitivity_to_bulk_arr_ptr,
        double complex* complex_bulk_modulus_arr_ptr,
        int degree_l,
        double G_to_use
        ) noexcept nogil:
    """ Calculate tidal heating as a function of radius.

    Assumptions
    -----------
    - World is in synchronous rotation with no obliquity.
    - Eccentricity is low enough to only include the e^2 truncation.

    """

    cdef double world_radius = radius_arr_ptr[total_slices - 1]

    # TODO: This term appears to contain all the information that could be upgraded in the future to eliminate some or
    #    all of the assumptions mentioned in the doc strings.
    cdef double portion_to_be_upgraded = (7. * eccentricity**2 * orbital_frequency)

    # The below is modified from TB05. That reference assumed l=2 and that the tidal host's mass was much larger than
    #    the target planet's mass. We have removed these restrictions in the below. The derivation was done by changing
    #    Eq. 2 in TB05 to the general form and then proceeding with their derivations (Eqs. 35--37) except that we did
    #    not set l=2.

    # To keep things consistent with other parts of TidalPy and the references that it is built off of, we are defining
    #    a few extra terms here. Some of their components are redundant and will be divided out shortly. We are choosing
    #    consistency at a (very) slight performance hit.
    cdef double tidal_susceptibility_overR = (3. / 2.) * G_to_use * host_mass**2 * world_radius**(5 - 1) / semi_major_axis**6


    radial_tidal_heating = \
        (tidal_susceptibility / world_radius) * (G / ((2. * order_l + 1.) * radius_array**2)) * portion_to_be_upgraded * \
        (radial_sensitivity_to_shear * np.imag(complex_shear_modulus) + 
         radial_sensitivity_to_bulk * np.imag(complex_bulk_modulus))

    # TODO: Finding that the heating rate as a function of depth can go negative. Setting those to zero for now.
    radial_tidal_heating[radial_tidal_heating < 0.] = 0.

    return radial_tidal_heating
