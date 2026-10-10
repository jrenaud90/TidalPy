#pragma once
/* How the shooting method carries each kind of layer: its independent solutions, the radial functions (ys) each holds
 * and in which order, the stored ys the surface conditions constrain, and the integration tolerance each y takes. A
 * solid carries 3 solutions of y1 to y6; a dynamic liquid 2 of (y1, y2, y5, y6); a static liquid 1 of (y5, y7)
 * (Takeuchi and Saito 1972; Saito 1974). The starting conditions, the surface conditions and their conditioning
 * estimate, the solution's collapse, the tolerance scaling, and the wrappers' buffer checks all read this one table.
 * A dynamic liquid is integrated with P = y2 - rho g y1 + rho y5 in y2's slot (derivatives/odes_.hpp); everything
 * outside its integration reads y2.
 */

#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>

#include "../constants_.hpp"   // TidalPyConstants::d_PI

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
    // A degree-1 loading solve fixes its reference frame with y5 = 1 at the surface (k' = 0, the frame of the body's
    // own center of mass), which takes the place of the last surface condition (y6 for a solid or a dynamic liquid,
    // y7 for a static liquid): the slots it solves, and the slot of the condition it drops.
    std::array<std::size_t, 3> degree1_frame_slots;
    std::size_t degree1_dropped_slot;
    // The factor on the integration rtol of each stored slot when the solve scales its tolerances by layer: tighter on
    // the solid's y2 and y3, which drive its instability (Issue #44). A dynamic liquid's P slot takes the plain rtol,
    // since P carries no such instability.
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
    3, 6, {0, 1, 2, 3, 4, 5}, {1, 3, 5}, {1, 3, 4}, 5, {1.0, 0.1, 0.1, 1.0, 1.0, 1.0}};
inline constexpr c_LayerKindLayout C_DYNAMIC_LIQUID_LAYOUT = {
    2, 4, {0, 1, 4, 5, C_Y_NOT_STORED, C_Y_NOT_STORED}, {1, 3, C_Y_NOT_STORED}, {1, 2, C_Y_NOT_STORED}, 3,
    {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}};
inline constexpr c_LayerKindLayout C_STATIC_LIQUID_LAYOUT = {
    1, 2, {4, 6, C_Y_NOT_STORED, C_Y_NOT_STORED, C_Y_NOT_STORED, C_Y_NOT_STORED},
    {1, C_Y_NOT_STORED, C_Y_NOT_STORED}, {0, C_Y_NOT_STORED, C_Y_NOT_STORED}, 1, {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}};

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

// The weight of each stored y of a layer kind, one over its characteristic size, in a planet of radius R and bulk
// density rho_bar, with gravitational constant G, in whatever units these are given. The sizes are
// 1 / (pi G rho_bar R) for the displacements y1 and y3, rho_bar for the stresses y2 and y4 (and a dynamic liquid's P),
// 1 for y5, and 1 / R for y6 and y7. They are TidalPy's non-dimensional units, so every weight is 1 in a
// non-dimensional solve. Unused slots hold 1.
inline std::array<double, 6> c_layer_y_weights(
        int layer_type,
        bool is_static,
        double radius,
        double bulk_density,
        double G) noexcept {
    const double displacement = TidalPyConstants::d_PI * G * bulk_density * radius;
    const double stress       = 1.0 / bulk_density;
    const std::array<double, C_NUM_FULL_YS> full_weights = {
        displacement, stress, displacement, stress, 1.0, radius, radius};
    const c_LayerKindLayout& layout = c_layer_layout(layer_type, is_static);
    std::array<double, 6> slot_weights = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    for (std::size_t slot = 0; slot < layout.num_ys; ++slot) {
        slot_weights[slot] = full_weights[layout.full_y_of_slot[slot]];
    }
    return slot_weights;
}

// The unit a dynamic liquid's P is stored in, as a fraction of a stress: omega^2 / (pi G rho_bar) below the frequency
// unit sqrt(pi G rho_bar), else 1 (any consistent units). P = -rho omega^2 r y3 is of that order, so the stored P, of
// order rho r y3 in a stress unit, holds y3 to the integration tolerance and, through a re-orthonormalization, to
// working precision however long the period. In a stress unit the integrator's absolute tolerance and the stress
// weight of the independence would see only P's roundoff.
inline double c_pressure_unit(double frequency, double bulk_density, double G) noexcept {
    return std::min(1.0, frequency * frequency / (TidalPyConstants::d_PI * G * bulk_density));
}

// A dynamic liquid's stored P = (y2 - rho g y1 + rho y5) / pressure_unit from TS72's y2, and back, at a point of
// gravity g and density rho.
inline std::complex<double> c_pressure_variable_from_y2(
        const std::complex<double>& y1,
        const std::complex<double>& y2,
        const std::complex<double>& y5,
        double gravity,
        double density,
        double pressure_unit) noexcept {
    return (y2 - density * gravity * y1 + density * y5) / pressure_unit;
}

inline std::complex<double> c_y2_from_pressure_variable(
        const std::complex<double>& y1,
        const std::complex<double>& pressure_variable,
        const std::complex<double>& y5,
        double gravity,
        double density,
        double pressure_unit) noexcept {
    return pressure_unit * pressure_variable + density * gravity * y1 - density * y5;
}
