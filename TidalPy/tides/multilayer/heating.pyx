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
from libc.math cimport isnan
from libcpp cimport bool as cpp_bool
from libcpp.vector cimport vector

import numpy as np
cimport numpy as cnp
cnp.import_array()

from TidalPy.exceptions import ArgumentException
from TidalPy.constants cimport d_G, d_NAN_DBL
from TidalPy.RadialSolver cimport RadialSolverSolution
from TidalPy.Material.eos.eos_solution cimport EOSSolutionCC
from TidalPy.tides.multilayer.sensitivity cimport cf_calc_sensitivity_to_shear, cf_calc_sensitivity_to_bulk


cdef void cf_calc_radial_volumetric_tidal_heating(
        double* volumetric_tidal_heating_arr_ptr,
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

    References
    ----------
    Eq. 37 in Tobie+ (2005)

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
    cdef double tidal_susceptibility_overR = \
        (3. / 2.) * G_to_use * host_mass**2 * world_radius**(5 - 1) / semi_major_axis**6
    cdef double G_2lp1 = G_to_use / (2.0 * degree_l + 1.0)

    cdef double r
    cdef Py_ssize_t slice_i
    cdef Py_ssize_t total_slices_ssize = <Py_ssize_t>total_slices
    cdef double volumetric_tidal_heating
    for slice_i in range(total_slices_ssize):
        r = radius_arr_ptr[slice_i]
        if r == 0.0:
            volumetric_tidal_heating_arr_ptr[slice_i] = 0.0
        else:
            volumetric_tidal_heating = \
                tidal_susceptibility_overR * G_2lp1 / (r**2) * portion_to_be_upgraded * (
                    radial_sensitivity_to_bulk_arr_ptr[slice_i] * complex_bulk_modulus_arr_ptr[slice_i].imag +
                    radial_sensitivity_to_shear_arr_ptr[slice_i] * complex_shear_modulus_arr_ptr[slice_i].imag
                )
            if (volumetric_tidal_heating < 0.0) or isnan(volumetric_tidal_heating):
                # TODO: Finding that the heating rate as a function of depth can go negative. Setting those to zero for now.
                # TODO: This check was from a older version of TidalPy (Pre v0.5.0) so it may no longer be needed.
                volumetric_tidal_heating_arr_ptr[slice_i] = 0.0
            else:
                volumetric_tidal_heating_arr_ptr[slice_i] = volumetric_tidal_heating


def calc_radial_volumetric_tidal_heating(
        double eccentricity,
        double orbital_frequency,
        double semi_major_axis,
        double host_mass,
        double[::1] radius_arr,
        double[::1] radial_sensitivity_to_shear_arr,
        double complex[::1] complex_shear_modulus_arr,
        double[::1] radial_sensitivity_to_bulk_arr,
        double complex[::1] complex_bulk_modulus_arr,
        int degree_l = 2,
        double G_to_use = d_G,
        cpp_bool perform_checks = True
        ):

    # Get size of arrays, perform checks
    cdef size_t total_slices = radius_arr.size
    if perform_checks:
        if radial_sensitivity_to_shear_arr.size != total_slices:
            raise ArgumentException("Unexpected size found for `radial_sensitivity_to_shear_arr`.")
        if complex_shear_modulus_arr.size != total_slices:
            raise ArgumentException("Unexpected size found for `complex_shear_modulus_arr`.")
        if radial_sensitivity_to_bulk_arr.size != total_slices:
            raise ArgumentException("Unexpected size found for `radial_sensitivity_to_bulk_arr`.")
        if complex_bulk_modulus_arr.size != total_slices:
            raise ArgumentException("Unexpected size found for `complex_bulk_modulus_arr`.")
    
    # Setup array pointers
    cdef double* radius_arr_ptr                        = &radius_arr[0]
    cdef double* radial_sensitivity_to_shear_arr_ptr   = &radial_sensitivity_to_shear_arr[0]
    cdef double complex* complex_shear_modulus_arr_ptr = &complex_shear_modulus_arr[0]
    cdef double* radial_sensitivity_to_bulk_arr_ptr    = &radial_sensitivity_to_bulk_arr[0]
    cdef double complex* complex_bulk_modulus_arr_ptr  = &complex_bulk_modulus_arr[0]
    
    # Create output array, initialize with zero instead of empty / nan
    cdef cnp.ndarray[cnp.float64_t, ndim=1] volumetric_tidal_heating_arr = np.zeros(total_slices, dtype=np.float64, order='C')
    cdef double[::1] volumetric_tidal_heating_arr_view = volumetric_tidal_heating_arr
    cdef double* volumetric_tidal_heating_arr_ptr      = &volumetric_tidal_heating_arr_view[0]

    # Call Cythonized function
    cf_calc_radial_volumetric_tidal_heating(
        volumetric_tidal_heating_arr_ptr,
        total_slices,
        eccentricity,
        orbital_frequency,
        semi_major_axis,
        host_mass,
        radius_arr_ptr,
        radial_sensitivity_to_shear_arr_ptr,
        complex_shear_modulus_arr_ptr,
        radial_sensitivity_to_bulk_arr_ptr,
        complex_bulk_modulus_arr_ptr,
        degree_l,
        G_to_use
        )
    
    return volumetric_tidal_heating_arr


def calc_radial_volumetric_tidal_heating_from_rs_solution(
        double eccentricity,
        double orbital_frequency,
        double semi_major_axis,
        double host_mass,
        RadialSolverSolution rs_solution,
        cpp_bool perform_checks = True
        ):

    # Build output arrays
    cdef size_t num_ytypes = rs_solution.num_ytypes
    if num_ytypes != 1:
        raise NotImplementedError("`calc_radial_volumetric_tidal_heating_from_rs_solution` currently only supports radial solver solutions with one ytype (e.g., just 'tidal' or just 'loading').")
    cdef size_t total_slices = rs_solution.radius_array_size
    cdef cnp.ndarray[cnp.float64_t, ndim=1] radial_volumetric_tidal_heating_arr = np.zeros(total_slices, dtype=np.float64, order='C')
    cdef double[::1] radial_volumetric_tidal_heating_view = radial_volumetric_tidal_heating_arr
    cdef double* radial_volumetric_tidal_heating_ptr      = &radial_volumetric_tidal_heating_view[0]

    # Get other metadata
    cdef int degree_l = rs_solution.degree_l

    # Build intermediate arrays
    cdef vector[double] radial_sensitivity_to_bulk_vec = vector[double]()
    radial_sensitivity_to_bulk_vec.resize(total_slices)
    cdef double* radial_sensitivity_to_bulk_ptr = &radial_sensitivity_to_bulk_vec[0]

    cdef vector[double] radial_sensitivity_to_shear_vec = vector[double]()
    radial_sensitivity_to_shear_vec.resize(total_slices)
    cdef double* radial_sensitivity_to_shear_ptr = &radial_sensitivity_to_shear_vec[0]

    # Get physical state pointers
    cdef EOSSolutionCC* eos_solution_ptr     = rs_solution.solution_storage_ptr.get_eos_solution_ptr()
    cdef double* radius_array_ptr            = eos_solution_ptr.radius_array_vec.data()
    cdef double* shear_modulus_array_dbl_ptr = eos_solution_ptr.complex_shear_array_vec.data()
    cdef double* bulk_modulus_array_dbl_ptr  = eos_solution_ptr.complex_bulk_array_vec.data()
    cdef double* rs_radial_solution_dbl_ptr  = rs_solution.solution_storage_ptr.full_solution_vec.data()

    # Convert the float version of the shear/bulk and radial solutions to double complex
    cdef double complex* shear_modulus_array_ptr = <double complex*>shear_modulus_array_dbl_ptr
    cdef double complex* bulk_modulus_array_ptr  = <double complex*>bulk_modulus_array_dbl_ptr
    cdef double complex* rs_radial_solution_ptr  = <double complex*>rs_radial_solution_dbl_ptr
    
    # Calculate the radial sensitivity to shear and bulk
    cf_calc_sensitivity_to_bulk(
        radial_sensitivity_to_bulk_ptr,
        rs_radial_solution_ptr, 
        radius_array_ptr,
        shear_modulus_array_ptr,
        bulk_modulus_array_ptr,
        total_slices,
        num_ytypes,
        degree_l
        )
    
    cf_calc_sensitivity_to_shear(
        radial_sensitivity_to_shear_ptr,
        rs_radial_solution_ptr, 
        radius_array_ptr,
        shear_modulus_array_ptr,
        bulk_modulus_array_ptr,
        total_slices,
        num_ytypes,
        degree_l
        )

    # Now calculate the volumetric heating
    cf_calc_radial_volumetric_tidal_heating(
        radial_volumetric_tidal_heating_ptr,
        total_slices,
        eccentricity,
        orbital_frequency,
        semi_major_axis,
        host_mass,
        radius_array_ptr,
        radial_sensitivity_to_shear_ptr,
        shear_modulus_array_ptr,
        radial_sensitivity_to_bulk_ptr,
        bulk_modulus_array_ptr,
        degree_l,
        d_G
        )
    
    return radial_volumetric_tidal_heating_arr
