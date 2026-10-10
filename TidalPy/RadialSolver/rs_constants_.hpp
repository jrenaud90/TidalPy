#pragma once

#include <array>
#include <complex>
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

/// A change of basis of a layer's solutions (orthonormalize_.hpp), row-major [row * C_MAX_NUM_SOL + column]: the
/// solutions before it are those after it times it. The upper-triangular R of a re-orthonormalization, or a dynamic
/// liquid's split of its starting solutions.
using c_BasisChange = std::array<std::complex<double>, C_MAX_NUM_SOL * C_MAX_NUM_SOL>;

inline c_BasisChange c_identity_basis_change() noexcept
{
    c_BasisChange identity{};
    for (size_t solution_i = 0; solution_i < C_MAX_NUM_SOL; ++solution_i)
    {
        identity[solution_i * C_MAX_NUM_SOL + solution_i] = std::complex<double>(1.0, 0.0);
    }
    return identity;
}
