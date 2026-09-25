#pragma once
/* Split an ascending radius array into per-layer runs.
 *
 * A layered profile repeats each interface radius, once as the top of the layer below and once as the base
 * of the layer above, so the two sides can carry their own density and moduli. Every consumer that indexes
 * such an array by layer must agree on where one run ends, and they stop agreeing the moment the rule is
 * written out twice. This is that rule, written once.
 *
 * A run starts at the first unclaimed slice and ends at the first copy of its upper radius; the second copy
 * starts the layer above. A slice past the upper radius with no copy at it also ends the run, which is what
 * happens when a caller's boundaries do not fall exactly on its grid.
 */

#include <cstddef>
#include <vector>

#include "../math_x/numerics_.hpp"  // c_isclose

namespace tidalpy {

/// Interface copies are written from the same value so they match far inside this; the tolerance is here
/// for a boundary that arrived through a conversion.
inline constexpr double d_LAYER_BOUNDARY_RTOL = 1.0e-9;

/// Writes the first index and slice count of each layer. A layer that finds no slices gets a count of
/// zero, which the caller checks against whatever minimum it needs: five for the shooting method, two for
/// an interpolated material.
inline void c_partition_radius_by_layer(
        const double* radius_ptr,
        const std::size_t num_slices,
        const double* upper_radius_bylayer_ptr,
        const std::size_t num_layers,
        std::vector<std::size_t>& first_slice_out,
        std::vector<std::size_t>& num_slices_out) noexcept
{
    first_slice_out.assign(num_layers, 0);
    num_slices_out.assign(num_layers, 0);
    if (radius_ptr == nullptr || upper_radius_bylayer_ptr == nullptr || num_slices == 0)
    {
        return;
    }

    std::size_t next_first = 0;
    for (std::size_t layer_i = 0; layer_i < num_layers; ++layer_i)
    {
        const double layer_upper = upper_radius_bylayer_ptr[layer_i];
        first_slice_out[layer_i] = next_first;

        std::size_t count          = 0;
        std::size_t interface_hits = 0;
        for (std::size_t slice_i = next_first; slice_i < num_slices; ++slice_i)
        {
            const double radius_here = radius_ptr[slice_i];
            if (c_isclose(radius_here, layer_upper, d_LAYER_BOUNDARY_RTOL, 0.0))
            {
                if (++interface_hits > 1)
                {
                    break;
                }
            }
            else if (radius_here > layer_upper)
            {
                break;
            }
            ++count;
        }
        num_slices_out[layer_i] = count;
        next_first += count;
    }
}

} // namespace tidalpy
