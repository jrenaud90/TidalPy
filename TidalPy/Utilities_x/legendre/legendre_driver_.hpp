#pragma once
/*
 * legendre_driver_.hpp - Dispatch to the precomputed associated-Legendre tables (l = 2..10) and a
 * generic (any l, m) fallback. Returns P_lm(cos theta) with its first/second colatitude derivatives.
 * All angles in radians; colatitude theta in [0, pi]. Uses the Condon-Shortley phase convention.
 */

#include <cmath>
#include <cstdint>

#include "constants_.hpp"
#include "legendre_common_.hpp"
#include "legendre_l2_.hpp"
#include "legendre_l3_.hpp"
#include "legendre_l4_.hpp"
#include "legendre_l5_.hpp"
#include "legendre_l6_.hpp"
#include "legendre_l7_.hpp"
#include "legendre_l8_.hpp"
#include "legendre_l9_.hpp"
#include "legendre_l10_.hpp"

namespace tidalpy {

constexpr int C_LEGENDRE_MAX_TABLE_DEGREE = 10;

// Precomputed associated Legendre triple for supported degrees (l = 2..10). Returns NaN triple for an
// unsupported degree or an out-of-range order (m < 0 or m > l).
inline c_LegendreValue c_legendre(int degree_l, int order_m, double colatitude)
{
    const double cos_t = std::cos(colatitude);
    const double sin_t = std::sin(colatitude);
    if (order_m < 0 || order_m > degree_l)
    {
        return c_LegendreValue {TidalPyConstants::d_NAN, TidalPyConstants::d_NAN, TidalPyConstants::d_NAN};
    }
    switch (degree_l)
    {
    case 2: return c_legendre_l2(order_m, cos_t, sin_t);
    case 3: return c_legendre_l3(order_m, cos_t, sin_t);
    case 4: return c_legendre_l4(order_m, cos_t, sin_t);
    case 5: return c_legendre_l5(order_m, cos_t, sin_t);
    case 6: return c_legendre_l6(order_m, cos_t, sin_t);
    case 7: return c_legendre_l7(order_m, cos_t, sin_t);
    case 8: return c_legendre_l8(order_m, cos_t, sin_t);
    case 9: return c_legendre_l9(order_m, cos_t, sin_t);
    case 10: return c_legendre_l10(order_m, cos_t, sin_t);
    default:
        return c_LegendreValue {TidalPyConstants::d_NAN, TidalPyConstants::d_NAN, TidalPyConstants::d_NAN};
    }
}

} // namespace tidalpy
