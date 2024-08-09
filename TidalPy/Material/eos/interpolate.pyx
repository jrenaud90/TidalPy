# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from CyRK.array.interp cimport interpj_ptr, interp_ptr

from TidalPy.Material.eos.common cimport EOSOutput
from TidalPy.utilities.math.complex cimport cmplx_NAN
from TidalPy.utilities.constants_x cimport NAN_DBL

cdef void preeval_interpolate(
        # Values that will be updated by the function
        void* preeval_output,
        # Input that is used by the pre-eval
        double radius,
        double* radial_solutions,
        void* preeval_input
        ) noexcept nogil:

    # Cast input to the proper structure for this function
    cdef ODEInput* ode_args = <ODEInput*>preeval_input
    cdef interpolate_input* eos_data = <interpolate_input*>ode_args.eos_input_ptr

    # Cast output to the proper structure
    cdef EOSOutput* output = <EOSOutput*>preeval_output

    # Set state variables based on an interpolation using the provided radius.
    cdef (double, Py_ssize_t) interp_out
    cdef Py_ssize_t index_j

    # The first interpolation will be the slowest as it must find the closest index.
    # We will use this index in the other interpolations.
    interp_out = interpj_ptr(
        radius,
        eos_data.radius_array_ptr,
        eos_data.density_array_ptr,
        eos_data.num_slices,
        provided_j=-2
        )

    output.density = interp_out[0]
    index_j        = interp_out[1]

    # Rheology C++ class only works with double pointers so we need to convert the double complex output to double*
    cdef double* complex_bulk_ptr = <double*>&output.bulk_modulus
    cdef double static_bulk
    cdef double bulk_viscosity
    if ode_args.update_bulk:
        static_bulk = interp_ptr(
            radius,
            eos_data.radius_array_ptr,
            eos_data.bulk_modulus_array_ptr,
            eos_data.num_slices,
            provided_j=index_j
            )
        
        bulk_viscosity = interp_ptr(
            radius,
            eos_data.radius_array_ptr,
            eos_data.bulk_viscosity_array_ptr,
            eos_data.num_slices,
            provided_j=index_j
            )
        
        # Apply bulk rheology
        output.static_bulk_modulus = static_bulk
        output.bulk_viscosity      = bulk_viscosity
        eos_data.bulk_rheology.get().call(complex_bulk_ptr, ode_args.frequency, static_bulk, bulk_viscosity)
    else:
        output.static_bulk_modulus = NAN_DBL
        output.bulk_viscosity      = NAN_DBL
        output.bulk_modulus        = cmplx_NAN

    cdef double* complex_shear_ptr = <double*>&output.shear_modulus
    cdef double static_shear
    cdef double shear_viscosity
    if ode_args.update_shear:
        static_shear = interp_ptr(
            radius,
            eos_data.radius_array_ptr,
            eos_data.shear_modulus_array_ptr,
            eos_data.num_slices,
            provided_j=index_j
            )

        shear_viscosity = interp_ptr(
            radius,
            eos_data.radius_array_ptr,
            eos_data.shear_viscosity_array_ptr,
            eos_data.num_slices,
            provided_j=index_j
            )

        # Apply shear rheology
        output.static_shear_modulus = static_shear
        output.shear_viscosity      = shear_viscosity
        eos_data.shear_rheology.get().call(complex_shear_ptr, ode_args.frequency, static_shear, shear_viscosity)
    else:
        output.static_shear_modulus = NAN_DBL
        output.shear_viscosity      = NAN_DBL
        output.shear_modulus        = cmplx_NAN
