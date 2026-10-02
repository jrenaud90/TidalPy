#pragma once

#include <cstddef>


/// Maximum number of y-values (radial functions y1--y6) per solution.
constexpr size_t C_MAX_NUM_Y = 6;

/// Maximum number of real values per solution (2x complex: real + imaginary for each y).
constexpr size_t C_MAX_NUM_Y_REAL = 2 * C_MAX_NUM_Y;

/// Maximum number of independent solutions per layer (solid layers have 3).
constexpr size_t C_MAX_NUM_SOL = 3;

/// Maximum number of surface boundary conditions ("ytypes") one solve may produce
constexpr size_t C_MAX_NUM_YTYPES = 5;

/// Surface conditions per ytype (on y2, y4, y6), and the buffer that holds them for every ytype.
constexpr size_t C_NUM_SURFACE_CONDITIONS = 3;
constexpr size_t C_MAX_SURFACE_CONDITIONS = C_MAX_NUM_YTYPES * C_NUM_SURFACE_CONDITIONS;

/// Fewest radial slices a layer of a supplied profile, or of the propagation matrix's grid, may have.
constexpr size_t C_RS_MIN_SLICES_PER_LAYER = 5;
