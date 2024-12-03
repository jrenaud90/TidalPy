# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from CyRK.array.interp cimport interpj_ptr, interp_ptr, interp_complex_ptr

from TidalPy.Material.eos.ode cimport EOSOutput
from TidalPy.utilities.math.complex cimport cmplx_NAN

cdef void preeval_interpolate(
        # Values that will be updated by the function
        char* preeval_output,
        # Input that is used by the pre-eval
        double radius,
        double* radial_solutions,
        char* preeval_input
        ) noexcept nogil:

    # Cast input to the proper structure for this function
    cdef EOS_ODEInput* ode_args        = <EOS_ODEInput*>preeval_input
    cdef InterpolateEOSInput* eos_data = <InterpolateEOSInput*>ode_args.eos_input_ptr

    # Cast output to the proper structure
    cdef EOSOutput* output = <EOSOutput*>preeval_output

    # Set state variables based on an interpolation using the provided radius.
    # The first interpolation will be the slowest as it must find the closest index.
    # We will use this index in the other interpolations.
    cdef (double, Py_ssize_t) interp_out = interpj_ptr(
        radius,
        eos_data.radius_array_ptr,
        eos_data.density_array_ptr,
        eos_data.num_slices,
        provided_j=-2
        )

    cdef Py_ssize_t index_j
    output.density = interp_out[0]
    index_j        = interp_out[1]

    # TODO: Interpolate the static bulk and shear mod; interpolate the bulk and shear viscosity, then apply rheology here.
    cdef double complex complex_bulk
    if ode_args.update_bulk:
        complex_bulk = interp_complex_ptr(
            radius,
            eos_data.radius_array_ptr,
            eos_data.bulk_modulus_array_ptr,
            eos_data.num_slices,
            provided_j=index_j
            )
        
        output.bulk_modulus = complex_bulk
    else:
        output.bulk_modulus = cmplx_NAN

    cdef double complex complex_shear
    if ode_args.update_shear:
        complex_shear = interp_complex_ptr(
            radius,
            eos_data.radius_array_ptr,
            eos_data.shear_modulus_array_ptr,
            eos_data.num_slices,
            provided_j=index_j
            )

        # Apply shear rheology
        output.shear_modulus = complex_shear
    else:
        output.shear_modulus = cmplx_NAN

    # Done
