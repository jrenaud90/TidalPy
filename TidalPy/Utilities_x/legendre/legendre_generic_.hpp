#pragma once
/*
 * legendre_generic_.hpp - Generic (any degree l, order m) associated Legendre evaluator built on the
 * vendored xsf special-function library. Returns P_lm(cos theta) with its first and second colatitude
 * derivatives, matching the Condon-Shortley phase convention of the precomputed tables
 * (legendre_l{2..10}_.hpp).
 *
 * The precomputed tables (c_legendre in legendre_driver_.hpp) cover l = 2..10 and are the fast path;
 * this generic routine is a fallback for degrees outside that range (not needed by the current tidal
 * pipeline, but kept beside the tables so future arbitrary-degree work single-sources the convention).
 *
 * xsf evaluates P_lm(x) with x = cos(theta); the derivatives come from xsf's dual (auto-diff) numbers
 * wrt x, converted to colatitude by the chain rule:
 *   dP/dtheta   = -sin(theta) * dP/dx
 *   d2P/dtheta2 =  sin^2(theta) * d2P/dx2 - cos(theta) * dP/dx
 */

#include <cmath>
#include <cstdint>

#include "xsf/legendre.h"
#include "xsf/dual.h"

#include "constants_.hpp"
#include "legendre_common_.hpp"

namespace tidalpy {

// Generic unnormalized associated Legendre triple (Condon-Shortley phase) at colatitude theta [rad].
// Returns a NaN triple for an out-of-range order (m < 0 or m > l). Valid for any degree l >= 0; use the
// precomputed c_legendre for l = 2..10 (faster).
inline c_LegendreValue c_legendre_generic(int degree_l, int order_m, double colatitude)
{
    if (order_m < 0 || order_m > degree_l)
    {
        return c_LegendreValue {TidalPyConstants::d_NAN, TidalPyConstants::d_NAN, TidalPyConstants::d_NAN};
    }

    const double cos_t = std::cos(colatitude);
    const double sin_t = std::sin(colatitude);

    // Second-order dual in x = cos(theta): value, dP/dx, d2P/dx2. branch_cut = 2 selects the
    // real on-interval (Ferrers) associated Legendre functions, matching scipy / the tables.
    xsf::dual<double, 2> x_dual({cos_t, 1.0, 0.0});
    xsf::dual<double, 2> p_dual =
        xsf::assoc_legendre_p(xsf::assoc_legendre_unnorm_policy{}, degree_l, order_m, x_dual, 2);

    const double p_val   = p_dual[0];
    const double dp_dx   = p_dual[1];
    const double d2p_dx2 = p_dual[2];   // xsf dual stores the true derivative (not a Taylor coefficient)

    return c_LegendreValue {
        p_val,
        -sin_t * dp_dx,
        sin_t * sin_t * d2p_dx2 - cos_t * dp_dx};
}

} // namespace tidalpy
