#pragma once
/* Associated Legendre evaluator for any degree l and order m, built on the vendored xsf library. The
 * precomputed tables (c_legendre in legendre_driver_.hpp) cover l = 2..10 and are the fast path; this
 * covers everything else under the same Condon-Shortley phase convention.
 *
 * xsf evaluates P_lm(x) with x = cos(theta); the derivatives come from its dual (auto-diff) numbers wrt
 * x, converted to colatitude by the chain rule:
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

// Unnormalized associated Legendre triple at colatitude theta [rad], for any degree l >= 0. NaN triple
// for an out-of-range order.
inline c_LegendreValue c_legendre_generic(int degree_l, int order_m, double colatitude)
{
    if (order_m < 0 || order_m > degree_l)
    {
        return c_LegendreValue {TidalPyConstants::d_NAN, TidalPyConstants::d_NAN, TidalPyConstants::d_NAN};
    }

    const double cos_t = std::cos(colatitude);
    const double sin_t = std::sin(colatitude);

    // Second-order dual in x = cos(theta): value, dP/dx, d2P/dx2. branch_cut = 2 selects the real
    // on-interval (Ferrers) functions, matching scipy and the tables.
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
