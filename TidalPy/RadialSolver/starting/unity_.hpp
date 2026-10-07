// unity_.hpp: unit-vector starting conditions.
//
// Each independent solution starts as a unit vector: no assumption is made about the layer at the start, so the
// start holds singular content as well as regular. Integrating outward, the singular part decays relative to the
// regular part as (r0 / r)^(2l - 1) in a solid and (r0 / r)^(2l + 1) in a liquid, so the error it leaves at the surface
// falls as the starting radius r0 shrinks.
// The components are the leading ones of Martens' (2016) free-constant solutions: y1, y4, y6 in a solid; y1, y6 in a
// dynamic liquid; y7 in a static liquid. A comparison of every combination of unit vectors on homogeneous bodies
// found this choice among the most accurate; sets with y3, or y2 with y5, can leave the regular solutions degenerate.
#pragma once

#include <complex>
#include <cstddef>


// Fill starting_conditions_ptr (num_ys per solution) with unit vectors for the layer type (0 = solid, 1 = liquid).
inline void c_unity_starting_conditions(
        const int layer_type,
        const bool is_static,
        const size_t num_ys,
        std::complex<double>* starting_conditions_ptr) noexcept
{
    const bool is_liquid = (layer_type != 0);
    if (!is_liquid)
    {
        // Solid y1..y6: y1, y4, y6.
        constexpr size_t components[3] = {0, 3, 5};
        for (size_t solution_i = 0; solution_i < 3; ++solution_i)
        {
            for (size_t y_i = 0; y_i < 6; ++y_i)
            {
                starting_conditions_ptr[solution_i * num_ys + y_i] = 0.0;
            }
            starting_conditions_ptr[solution_i * num_ys + components[solution_i]] = 1.0;
        }
    }
    else if (is_static)
    {
        // Static liquid (y5, y7): y7.
        starting_conditions_ptr[0] = 0.0;
        starting_conditions_ptr[1] = 1.0;
    }
    else
    {
        // Dynamic liquid (y1, y2, y5, y6): y1, y6.
        constexpr size_t components[2] = {0, 3};
        for (size_t solution_i = 0; solution_i < 2; ++solution_i)
        {
            for (size_t y_i = 0; y_i < 4; ++y_i)
            {
                starting_conditions_ptr[solution_i * num_ys + y_i] = 0.0;
            }
            starting_conditions_ptr[solution_i * num_ys + components[solution_i]] = 1.0;
        }
    }
}
