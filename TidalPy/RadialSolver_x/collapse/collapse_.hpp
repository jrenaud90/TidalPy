// collapse_.hpp - Collapse the independent shooting solutions into a single radial solution.
#pragma once

#include <vector>
#include <complex>


// Collapse the independent shooting solutions into one, weighted by the constant vector.
//
// Parameters
// ----------
// solution_ptr : complex*
//     Output solution array.
// constant_vector_ptr : complex*
//     Scale factor for each independent solution.
// storage_by_solution : vector of vector of complex
//     Per-solution intermediate y values, slice-major.
// layer_radius_ptr, layer_density_ptr, layer_gravity_ptr : double*
//     Radius [m], density [kg m-3], and gravity [m s-2] at each slice of this layer.
// frequency_to_use : double
//     Forcing frequency [rad s-1]; only used for y_3 in dynamic liquid layers.
// num_ys, num_output_ys : size_t
//     Y count defined for this layer type, and in the output array.
// layer_type : int
//     0=solid, 1=liquid.
inline void c_collapse_layer_solution(
        std::complex<double>* solution_ptr,
        std::complex<double>* constant_vector_ptr,
        std::vector<std::vector<std::complex<double>>>& storage_by_solution,
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
