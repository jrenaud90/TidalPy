#pragma once
/* How the shooting method carries each kind of layer: its independent solutions, the radial functions (ys) each holds
 * and in which order, the stored ys the surface conditions constrain, and the integration tolerance each y takes. A
 * solid carries 3 solutions of y1 to y6; a dynamic liquid 2 of (y1, y2, y5, y6); a static liquid 1 of (y5, y7)
 * (Takeuchi and Saito 1972; Saito 1974). The starting conditions, the surface conditions and their conditioning
 * estimate, the solution's collapse, the tolerance scaling, and the wrappers' buffer checks all read this one table.
 */

#include <array>
#include <cstddef>

// A radial function by its full index: 0 to 5 for y1 to y6, and 6 for a static liquid's y7.
inline constexpr std::size_t C_NUM_FULL_YS = 7;
// A full y the layer kind does not store.
inline constexpr std::size_t C_Y_NOT_STORED = static_cast<std::size_t>(-1);

struct c_LayerKindLayout {
    std::size_t num_solutions;
    std::size_t num_ys;
    // The full y each stored slot holds, slots 0 to num_ys - 1.
    std::array<std::size_t, 6> full_y_of_slot;
    // The stored slots the surface boundary conditions constrain, one per solution: y2, y4, y6 for a solid; y2, y6
    // for a dynamic liquid; y7 for a static liquid.
    std::array<std::size_t, 3> surface_condition_slots;
    // The factor on the integration rtol of each stored slot: tighter on the ys that drive the instability (y2 and y3
    // of a solid, y2 of a dynamic liquid) when the solve scales its tolerances by layer (Issue #44 tests the values).
    std::array<double, 6> rtol_scale;

    // The stored slot of a full y, or C_Y_NOT_STORED.
    constexpr std::size_t slot_of_full_y(std::size_t full_y) const noexcept {
        for (std::size_t slot = 0; slot < this->num_ys; ++slot) {
            if (this->full_y_of_slot[slot] == full_y) { return slot; }
        }
        return C_Y_NOT_STORED;
    }
};

inline constexpr c_LayerKindLayout C_SOLID_LAYOUT = {
    3, 6, {0, 1, 2, 3, 4, 5}, {1, 3, 5}, {1.0, 0.1, 0.1, 1.0, 1.0, 1.0}};
inline constexpr c_LayerKindLayout C_DYNAMIC_LIQUID_LAYOUT = {
    2, 4, {0, 1, 4, 5, C_Y_NOT_STORED, C_Y_NOT_STORED}, {1, 3, C_Y_NOT_STORED}, {1.0, 0.01, 1.0, 1.0, 1.0, 1.0}};
inline constexpr c_LayerKindLayout C_STATIC_LIQUID_LAYOUT = {
    1, 2, {4, 6, C_Y_NOT_STORED, C_Y_NOT_STORED, C_Y_NOT_STORED, C_Y_NOT_STORED},
    {1, C_Y_NOT_STORED, C_Y_NOT_STORED}, {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}};

// layer_type: 0 = solid, anything else a liquid.
inline constexpr const c_LayerKindLayout& c_layer_layout(int layer_type, bool is_static) noexcept {
    return (layer_type == 0) ? C_SOLID_LAYOUT : (is_static ? C_STATIC_LIQUID_LAYOUT : C_DYNAMIC_LIQUID_LAYOUT);
}

inline constexpr std::size_t c_layer_num_solutions(int layer_type, bool is_static) noexcept {
    return c_layer_layout(layer_type, is_static).num_solutions;
}
inline constexpr std::size_t c_layer_num_ys(int layer_type, bool is_static) noexcept {
    return c_layer_layout(layer_type, is_static).num_ys;
}
