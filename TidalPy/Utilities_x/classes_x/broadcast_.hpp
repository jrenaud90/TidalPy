#pragma once
/* Broadcasting for the physics models' vectorized calls. Each input of an element-wise call holds either one value
 * per point or a single value used at every point, so a constant input is never copied out to the full length.
 */

#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <string>

namespace tidalpy {

// Number of points of an element-wise call over inputs of the given sizes. Every size must be 1 or the one common
// size; any other combination throws std::invalid_argument naming the caller.
inline std::size_t c_broadcast_length(std::initializer_list<std::size_t> input_sizes, const char* caller) {
    std::size_t num_points = 1;
    for (const std::size_t input_size : input_sizes) {
        if (input_size == 1) { continue; }
        if (num_points == 1) {
            num_points = input_size;
        }
        else if (input_size != num_points) {
            throw std::invalid_argument(
                std::string("TidalPy::") + caller + ": every input must hold one value or the same number of "
                "values, but inputs hold " + std::to_string(num_points) + " and " + std::to_string(input_size));
        }
    }
    return num_points;
}

// Index step through an input of the given size: 0 for a single value used at every point, else 1.
inline std::size_t c_broadcast_stride(std::size_t input_size) noexcept {
    return (input_size == 1) ? 0 : 1;
}

} // namespace tidalpy
