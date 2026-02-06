#pragma once

#include <cmath>
#include <algorithm>


bool c_isclose(
        double a,
        double b,
        double rtol = 1e-9,
        double atol = 0.0)
{
    // Check for nans
    if (std::isnan(a) || std::isnan(b))
    {
        return false;
    }

    // Check for pure equivalence
    if (a == b)
    {
        return true;
    }

    // Check for closeness
    double lhs = std::abs(a - b);
    double rhs = std::max(rtol * std::max(std::abs(a), std::abs(b)), atol);

    return lhs <= rhs;
}
