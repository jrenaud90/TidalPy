#pragma once

#include <cstddef>


/// Maximum number of y-values (radial functions y1--y6) per solution.
constexpr size_t C_MAX_NUM_Y = 6;

/// Maximum number of real values per solution (2x complex: real + imaginary for each y).
constexpr size_t C_MAX_NUM_Y_REAL = 2 * C_MAX_NUM_Y;

/// Maximum number of independent solutions per layer (solid layers have 3).
constexpr size_t C_MAX_NUM_SOL = 3;
