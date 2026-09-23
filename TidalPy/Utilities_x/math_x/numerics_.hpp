#pragma once
/* Small shared numerical helpers.
 *
 * Dependency-free (standard library only) so even the lowest-level utility modules can include it.
 */

#include <cmath>
#include <algorithm>
#include <limits>


/// NaN-safe; mirrors math.isclose.
inline bool c_isclose(
        double value_a,
        double value_b,
        double rtol = 1e-9,
        double atol = 0.0)
{
    if (std::isnan(value_a) || std::isnan(value_b))
    {
        return false;
    }

    if (value_a == value_b)
    {
        return true;
    }

    const double lhs = std::abs(value_a - value_b);
    const double rhs = std::max(rtol * std::max(std::abs(value_a), std::abs(value_b)), atol);

    return lhs <= rhs;
}


// Formulas with exponential or power-law growth silently overflow to inf (or a domain-error NaN). These
// wrappers return a quiet NaN instead, so a bad parameter regime is visible rather than propagating.

inline double c_safe_pow(double base, double exponent)
{
    const double result = std::pow(base, exponent);
    if (!std::isfinite(result))
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return result;
}

inline double c_safe_exp(double exponent)
{
    const double result = std::exp(exponent);
    if (!std::isfinite(result))
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return result;
}
