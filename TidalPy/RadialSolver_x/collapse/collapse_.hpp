// collapse_.hpp - Collapse multiple shooting solutions into single solution
// Ported from TidalPy/RadialSolver/collapse/collapse.pyx
#pragma once

#include <complex>


// Collapse multiple independent shooting solutions into a single solution using the constant vector.
//
// Parameters
// ----------
// solution_ptr : complex*
//     Output array for final computed solutions.
// constant_vector_ptr : complex*
//     Constants used to scale each independent solution.
// storage_by_solution : complex**
//     Array of pointers to precomputed intermediate solutions for each slice.
// layer_radius_ptr : double*
//     Radii for each slice.
// layer_density_ptr : double*
//     Density for each slice.
// layer_gravity_ptr : double*
//     Gravity for each slice.
// frequency_to_use : double
//     Frequency parameter for y_3 calculation in dynamic liquid layers.
// layer_start_index : size_t
//     Starting index for this layer's slices.
// num_layer_slices : size_t
//     Number of slices in this layer.
// num_sols : size_t
//     Number of independent solutions.
// max_num_y : size_t
//     Maximum number of y-values per solution.
// num_ys : size_t
//     Number of y-values defined for this layer type.
// num_output_ys : size_t
//     Number of y-values in the output array.
// ytype_i : size_t
//     Y-type index.
// layer_type : int
//     0=solid, 1=liquid.
// layer_is_static : bool
// layer_is_incomp : bool
inline void c_collapse_layer_solution(
        std::complex<double>* solution_ptr,
        std::complex<double>* constant_vector_ptr,
        std::complex<double>** storage_by_solution,
        double* layer_radius_ptr,
        double* layer_density_ptr,
        double* layer_gravity_ptr,
        double frequency_to_use,
        size_t layer_start_index,
        size_t num_layer_slices,
        size_t num_sols,
        size_t max_num_y,
        size_t num_ys,
        size_t num_output_ys,
        size_t ytype_i,
        int layer_type,
        bool layer_is_static,
        bool layer_is_incomp) noexcept
{
    size_t solution_i;
    size_t y_i;
    size_t y_rhs_i;
    size_t lhs_y_index;
    size_t slice_i;
    size_t slice_i_shifted;
    size_t slice_end;
    bool calculate_y3;

    slice_end    = layer_start_index + num_layer_slices;
    calculate_y3 = false;

    if ((layer_type != 0) && (!layer_is_static))
    {
        // For dynamic liquid layers, calculate y3 from other ys after the main loop.
        calculate_y3 = true;
    }

    y_rhs_i = 0;
    y_i = 0;
    while (max_num_y > y_i)
    {
        lhs_y_index = ytype_i * max_num_y + y_i;
        // Bail out early for ys that the layer type is not defined for.
        if (layer_type == 0)
        {
            // Solid layers have same number of ys as max_num_y.
            y_rhs_i = y_i;
        } else
        {
            // Liquid layers have fewer ys than max_num_y.
            if (layer_is_static) {
                // Static liquid only has y5 (stored at index 0).
                if (y_i == 4) {
                    y_rhs_i = 0;
                } else {
                    ++y_i;
                    continue;
                }
            } else
            {
                // Dynamic liquid has y1, y2, y5, y6 (indices 0, 1, 2, 3)
                if (y_i < 2) {
                    y_rhs_i = y_i;
                } else if ((y_i > 3) && (y_i < 6)) {
                    y_rhs_i = y_i - 2;
                } else {
                    ++y_i;
                    continue;
                }
            }
        }

        // Collapse solutions together for each slice.
        solution_i = 0;
        while (num_sols > solution_i)
        {
            slice_i = 0;
            slice_i_shifted = layer_start_index;
            while (slice_end > slice_i_shifted) {
                if (solution_i == 0)
                {
                    solution_ptr[slice_i_shifted * num_output_ys + lhs_y_index] =
                        (constant_vector_ptr[solution_i] *
                        storage_by_solution[solution_i][slice_i * num_ys + y_rhs_i]);
                } else
                {
                    solution_ptr[slice_i_shifted * num_output_ys + lhs_y_index] +=
                        (constant_vector_ptr[solution_i] *
                        storage_by_solution[solution_i][slice_i * num_ys + y_rhs_i]);
                }

                ++slice_i;
                ++slice_i_shifted;
            }
            ++solution_i;
        }
        ++y_i;
    }

    if (calculate_y3)
    {
        // Handle y3 for dynamic liquid layers.
        lhs_y_index = ytype_i * max_num_y;
        slice_i = 0;
        slice_i_shifted = layer_start_index;
        while (slice_end > slice_i_shifted)
        {
            solution_ptr[slice_i_shifted * num_output_ys + (lhs_y_index + 2)] =
                (1.0 / (frequency_to_use * frequency_to_use * layer_radius_ptr[slice_i])) *
                (solution_ptr[slice_i_shifted * num_output_ys + (lhs_y_index + 0)] * layer_gravity_ptr[slice_i] -
                solution_ptr[slice_i_shifted * num_output_ys + (lhs_y_index + 1)] / layer_density_ptr[slice_i] -
                solution_ptr[slice_i_shifted * num_output_ys + (lhs_y_index + 4)]);

            ++slice_i;
            ++slice_i_shifted;
        }
    }
}
