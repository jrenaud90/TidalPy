#pragma once

#include <complex>
#include <limits>
#include <utility>
#include "interp_.hpp"
#include "ode_.hpp"
#include "constants_.hpp"


/// Input data needed for interpolation-based EOS evaluation.
///
/// Stores pointers to arrays of radius, density, and complex moduli data
/// that are used by the interpolation pre-evaluation function.
struct c_InterpolateEOSInput
{
    size_t num_slices = 0;
    double* radius_array_ptr  = nullptr;
    double* density_array_ptr = nullptr;
    std::complex<double>* bulk_modulus_array_ptr  = nullptr;
    std::complex<double>* shear_modulus_array_ptr = nullptr;
};

// The function signature matches CyRK's PreEvalFunc:
//   void(char* preeval_output, double radius, double* radial_solutions, char* preeval_input)

void c_preeval_interpolate(
        // Values that will be updated by the function
        char* preeval_output,
        // Input that is used by the pre-eval
        double radius,
        double* radial_solutions,
        char* preeval_input
        ) noexcept
{
    // Cast input to the proper structure for this function
    c_EOS_ODEInput* ode_args = reinterpret_cast<c_EOS_ODEInput*>(preeval_input);
    c_InterpolateEOSInput* eos_data = reinterpret_cast<c_InterpolateEOSInput*>(ode_args->eos_input_ptr);

    // Cast output to the proper structure
    c_EOSOutput* output = reinterpret_cast<c_EOSOutput*>(preeval_output);

    // Find the shared index_j to use across all three interpolations.
    // We do this explicitly because the provided cf_interp functions read the pointer but don't output the new j.
    int b_search_code = 0;
    
    // Formulate the initial guess to pass to the binary search
    const double left_x = eos_data->radius_array_ptr[0];
    const double right_x = eos_data->radius_array_ptr[eos_data->num_slices - 1];

    // Find initial index guess
    size_t j_guess = static_cast<size_t>(
        eos_data->num_slices * std::fabs(std::floor(radius / (right_x - left_x)))
    );
    j_guess = std::max<size_t>(std::min<size_t>(j_guess, eos_data->num_slices), 0);
    size_t index_j = cf_binary_search_with_guess(
        radius, 
        eos_data->radius_array_ptr, 
        eos_data->num_slices, 
        j_guess, 
        &b_search_code
    );

    // Interpolate Density
    double density_result = 0.0;
    cf_interp(
        &radius,
        eos_data->radius_array_ptr,
        eos_data->density_array_ptr,
        eos_data->num_slices,
        &index_j,
        &density_result
    );
    output->density = density_result;


    // Interpolate Bulk Modulus
    if (ode_args->update_bulk)
    {
        double bulk_result[2] = {0.0, 0.0};
        
        // Cast the complex array to a double array (interleaved real/imag) for cf_interp_complex
        auto* bulk_ptr = reinterpret_cast<double*>(eos_data->bulk_modulus_array_ptr);
        
        cf_interp_complex(
            radius,                      // Note: cf_interp_complex takes double
            eos_data->radius_array_ptr,
            bulk_ptr,
            eos_data->num_slices,
            &index_j,
            bulk_result
        );
        output->bulk_modulus = std::complex<double>(bulk_result[0], bulk_result[1]);
    }
    else
    {
        output->bulk_modulus = {TidalPyConstants::d_NAN, 0.0};
    }


    // Interpolate Shear Modulus
    if (ode_args->update_shear)
    {
        double shear_result[2] = {TidalPyConstants::d_NAN, 0.0};
        
        auto* shear_ptr = reinterpret_cast<double*>(eos_data->shear_modulus_array_ptr);
        
        cf_interp_complex(
            radius,
            eos_data->radius_array_ptr,
            shear_ptr,
            eos_data->num_slices,
            &index_j,
            shear_result
        );
        
        // Apply shear rheology
        output->shear_modulus = std::complex<double>(shear_result[0], shear_result[1]);
    }
    else
    {
        output->shear_modulus = {TidalPyConstants::d_NAN, 0.0};
    }
}
