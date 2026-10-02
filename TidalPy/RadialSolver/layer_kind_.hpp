#pragma once
/* How many independent solutions the shooting method carries through a layer, and how many radial functions (ys) each
 * holds, by the kind of layer. A solid carries 3 solutions of y1 to y6, a dynamic liquid 2 of y1, y2, y5, y6, and a
 * static liquid 1 of y5 and y7 (Takeuchi and Saito 1972; Saito 1974). The starting conditions, the interface and
 * surface conditions, and the wrappers' buffer checks all read these.
 */

#include <cstddef>

// layer_type: 0 = solid, anything else a liquid.
inline constexpr std::size_t c_layer_num_solutions(int layer_type, bool is_static) noexcept {
    return (layer_type == 0) ? 3 : (is_static ? 1 : 2);
}
inline constexpr std::size_t c_layer_num_ys(int layer_type, bool is_static) noexcept {
    return (layer_type == 0) ? 6 : (is_static ? 2 : 4);
}
