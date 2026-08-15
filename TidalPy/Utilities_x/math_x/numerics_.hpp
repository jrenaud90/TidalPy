#pragma once
/*
 * numerics_.hpp — small shared numerical helpers for all TidalPy C++ code.
 *
 * Header-only and dependency-free (uses only the C++ standard library) so it can be
 * included from any extension, including the lowest-level utility modules. The NaN
 * returned by the guarded helpers is the same quiet NaN as TidalPyConstants::d_NAN.
 */

#include <cmath>
#include <algorithm>
#include <limits>


// =====================================================================================================================
// Floating-point comparison
// =====================================================================================================================

/// Check whether two doubles are approximately equal (NaN-safe; mirrors math.isclose).
inline bool c_isclose(
        double value_a,
        double value_b,
        double rtol = 1e-9,
        double atol = 0.0)
{
    // Check for nans
    if (std::isnan(value_a) || std::isnan(value_b))
    {
        return false;
    }

    // Check for pure equivalence
    if (value_a == value_b)
    {
        return true;
    }

    // Check for closeness
    const double lhs = std::abs(value_a - value_b);
    const double rhs = std::max(rtol * std::max(std::abs(value_a), std::abs(value_b)), atol);

    return lhs <= rhs;
}


// =====================================================================================================================
// Guarded growth functions
//
// Scientific formulas with exponential or power-law growth can silently overflow to inf
// (or produce NaN from a domain error). These wrappers return a quiet NaN instead so a
// bad parameter regime is visible to the caller rather than propagating garbage.
// =====================================================================================================================

/// std::pow that returns NaN instead of inf (or a domain-error NaN) when the result is not finite.
inline double c_safe_pow(double base, double exponent)
{
    const double result = std::pow(base, exponent);
    if (!std::isfinite(result))
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return result;
}

/// std::exp that returns NaN instead of inf when the result overflows.
inline double c_safe_exp(double exponent)
{
    const double result = std::exp(exponent);
    if (!std::isfinite(result))
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return result;
}
