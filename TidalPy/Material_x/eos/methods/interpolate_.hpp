#pragma once

#include <complex>
#include <limits>
#include <utility>

#include "../ode_.hpp"
#include "../../../constants_.hpp"
#include "../../../utilities/arrays/interp_.hpp"


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

    // Slice index the last pre-evaluation landed on, carried forward as the seed for the next
    // binary search.
    size_t last_slice_index = 0;
};

// The function signature matches CyRK's PreEvalFunc:
//   void(char* preeval_output, double radius, double* radial_solutions, char* preeval_input)

inline void c_preeval_interpolate(
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
    // We do this explicitly because the provided c_interp functions read the pointer but don't output the new j.
    int b_search_code = 0;

    // Seed the search with the slice the previous call used (c_InterpolateEOSInput::last_slice_index).
    // The estimate this replaced divided the radius by the layer's radius span and took the floor of the
    // result, which is zero for every radius inside the layer, so every call bisected the whole array.
    size_t index_j = c_binary_search_with_guess(
        radius,
        eos_data->radius_array_ptr,
        eos_data->num_slices,
        eos_data->last_slice_index,
        &b_search_code
    );

    // Carry the interval forward. The search returns num_slices for a radius past the top of the array
    // and 0 with b_search_code == -1 below its bottom, so clamp to a seedable index.
    eos_data->last_slice_index =
        (eos_data->num_slices > 0 && index_j >= eos_data->num_slices) ? eos_data->num_slices - 1 : index_j;

    // Interpolate Density
    double density_result = 0.0;
    c_interp(
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
        
        // Cast the complex array to a double array (interleaved real/imag) for c_interp_complex
        auto* bulk_ptr = reinterpret_cast<double*>(eos_data->bulk_modulus_array_ptr);
        
        c_interp_complex(
            radius,                      // Note: c_interp_complex takes double
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
        
        c_interp_complex(
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
