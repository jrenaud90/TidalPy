#pragma once
/* Where a radial solve starts, shared by the shooting method (shooting_.hpp) and the propagation matrix (matrix_.hpp) so
 * the two agree on the automatic start and on which layer a start radius lies in.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "../constants_.hpp"   // tidalpy_config_ptr

// The integration cannot start at r = 0 (a singularity); starting higher is more stable but skips more of the planet.
// A requested radius of 0 takes the automatic choice after Martens' thesis and the LoadDef manual,
// R tol^(1 / l), capped by [numerical] max_start_radius_fraction.
inline double c_auto_starting_radius(
        double requested_radius, double planet_radius, double start_radius_tol, double degree_l) noexcept {
    if (requested_radius != 0.0) { return requested_radius; }
    const double automatic = planet_radius * std::pow(start_radius_tol, 1.0 / degree_l);
    return std::fmin(automatic, tidalpy_config_ptr->d_MAX_START_RADIUS_FRAC * planet_radius);
}

// The layer holding the starting radius, [lower, upper), so a start on an interface begins in the layer above it.
// False for a radius at or below zero, at or above the surface, or NaN, which lies in no layer.
inline bool c_find_start_layer(
        const std::vector<double>& upper_radius_bylayer,
        std::size_t num_layers,
        double starting_radius,
        std::size_t& layer_out) noexcept {
    for (std::size_t layer_i = 0; layer_i < num_layers; ++layer_i) {
        const double lower = (layer_i == 0) ? 0.0 : upper_radius_bylayer[layer_i - 1];
        if ((starting_radius > 0.0) && (lower <= starting_radius) && (starting_radius < upper_radius_bylayer[layer_i])) {
            layer_out = layer_i;
            return true;
        }
    }
    return false;
}
