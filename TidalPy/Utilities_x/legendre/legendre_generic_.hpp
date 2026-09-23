#pragma once
/* Associated Legendre evaluator for any degree l and order m, built on the vendored xsf library. The
 * precomputed tables (c_legendre in legendre_driver_.hpp) cover l = 2..10 and are the fast path; this
 * covers everything else under the same Condon-Shortley phase convention.
 *
 * xsf evaluates P_l^m(x) with x = cos(theta). The colatitude derivatives come from the order recurrence
 *   dP_l^k/dtheta = (1/2) [P_l^(k+1) - (l+k)(l-k+1) P_l^(k-1)]        (Condon-Shortley phase)
 * applied once and twice, with P_l^(-k) = (-1)^k (l-k)!/(l+k)! P_l^k and P_l^k = 0 for |k| > l. Every term is a
 * smooth function of theta, so the poles are exact. (The chain rule through x, dP/dtheta = -sin(theta) dP/dx,
 * multiplies a derivative that is infinite at the poles for odd m by a zero there.)
 */

#include <cmath>
#include <cstdint>

#include "xsf/legendre.h"

#include "constants_.hpp"
#include "legendre_common_.hpp"

namespace tidalpy {

namespace detail {

// P_l^k(cos theta) for any integer order k. branch_cut = 2 selects the real on-interval (Ferrers) functions,
// matching scipy and the tables.
inline double legendre_any_order(int degree_l, int order_k, double cos_t)
{
    if (order_k > degree_l || order_k < -degree_l) { return 0.0; }
    const int abs_k = (order_k < 0) ? -order_k : order_k;
    const double value = xsf::assoc_legendre_p(xsf::assoc_legendre_unnorm_policy{}, degree_l, abs_k, cos_t, 2);
    if (order_k >= 0) { return value; }
    // P_l^(-k) = (-1)^k (l-k)!/(l+k)! P_l^k.
    double ratio = 1.0;
    for (int factor = degree_l - abs_k + 1; factor <= degree_l + abs_k; ++factor) {
        ratio /= static_cast<double>(factor);
    }
    return ((abs_k % 2 == 0) ? ratio : -ratio) * value;
}

// dP_l^k/dtheta for any integer order k, from the order recurrence.
inline double legendre_dtheta_any_order(int degree_l, int order_k, double cos_t)
{
    const double ladder = static_cast<double>(degree_l + order_k) * static_cast<double>(degree_l - order_k + 1);
    return 0.5 * (legendre_any_order(degree_l, order_k + 1, cos_t)
                  - ladder * legendre_any_order(degree_l, order_k - 1, cos_t));
}

}  // namespace detail

// Unnormalized associated Legendre triple at colatitude theta [rad], for any degree l >= 0. NaN triple
// for an out-of-range order.
inline c_LegendreValue c_legendre_generic(int degree_l, int order_m, double colatitude)
{
    if (order_m < 0 || order_m > degree_l)
    {
        return c_LegendreValue {TidalPyConstants::d_NAN, TidalPyConstants::d_NAN, TidalPyConstants::d_NAN};
    }

    const double cos_t  = std::cos(colatitude);
    const double ladder = static_cast<double>(degree_l + order_m) * static_cast<double>(degree_l - order_m + 1);

    return c_LegendreValue {
        detail::legendre_any_order(degree_l, order_m, cos_t),
        detail::legendre_dtheta_any_order(degree_l, order_m, cos_t),
        0.5 * (detail::legendre_dtheta_any_order(degree_l, order_m + 1, cos_t)
               - ladder * detail::legendre_dtheta_any_order(degree_l, order_m - 1, cos_t))};
}

} // namespace tidalpy
