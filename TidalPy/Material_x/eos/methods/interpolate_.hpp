#pragma once

#include <complex>
#include <limits>
#include <utility>

#include "../ode_.hpp"
#include "../../../constants_.hpp"
#include "../../../utilities/arrays/interp_.hpp"


/// Non-owning table pointers for the interpolation EOS pre-evaluation.
struct c_InterpolateEOSInput
{
    size_t num_slices = 0;
    double* radius_array_ptr  = nullptr;
    double* density_array_ptr = nullptr;
    std::complex<double>* bulk_modulus_array_ptr  = nullptr;
    std::complex<double>* shear_modulus_array_ptr = nullptr;

    // Slice index of the last pre-evaluation, the seed for the next binary search.
    size_t last_slice_index = 0;
};

/// Interpolation EOS pre-evaluation (CyRK PreEvalFunc signature).
inline void c_preeval_interpolate(
        char* preeval_output,
        double radius,
        double* radial_solutions,
        char* preeval_input
        ) noexcept
{
    c_EOS_ODEInput* ode_args = reinterpret_cast<c_EOS_ODEInput*>(preeval_input);
    c_InterpolateEOSInput* eos_data = reinterpret_cast<c_InterpolateEOSInput*>(ode_args->eos_input_ptr);
    c_EOSOutput* output = reinterpret_cast<c_EOSOutput*>(preeval_output);

    // One search shared by the three interpolations (c_interp reads the index but does not update it).
    int b_search_code = 0;
    size_t index_j = c_binary_search_with_guess(
        radius,
        eos_data->radius_array_ptr,
        eos_data->num_slices,
        eos_data->last_slice_index,
        &b_search_code
    );

    // The search returns num_slices past the top of the array, so clamp to a seedable index.
    eos_data->last_slice_index =
        (eos_data->num_slices > 0 && index_j >= eos_data->num_slices) ? eos_data->num_slices - 1 : index_j;

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


    if (ode_args->update_bulk)
    {
        double bulk_result[2] = {0.0, 0.0};

        // c_interp_complex takes the complex array as interleaved real/imag doubles.
        auto* bulk_ptr = reinterpret_cast<double*>(eos_data->bulk_modulus_array_ptr);

        c_interp_complex(
            radius,
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

        output->shear_modulus = std::complex<double>(shear_result[0], shear_result[1]);
    }
    else
    {
        output->shear_modulus = {TidalPyConstants::d_NAN, 0.0};
    }
}
