#pragma once
/* Associated Legendre evaluator for any degree l and order m. The precomputed tables (c_legendre in
 * legendre_driver_.hpp) cover l = 2..10 and are the fast path; this covers everything else under the same
 * unnormalized, Condon-Shortley phase convention (scipy's assoc_legendre_p with branch_cut = 2).
 *
 * The values come from the standard upward recurrence in degree at fixed order, written in cos(theta) and
 * sin(theta) themselves:
 *   P_k^k = (-1)^k (2k - 1)!! sin^k(theta),   P_(k+1)^k = (2k + 1) cos(theta) P_k^k,
 *   (l - k) P_l^k = (2l - 1) cos(theta) P_(l-1)^k - (l + k - 1) P_(l-2)^k.
 * An evaluator in x = cos(theta) alone rebuilds sin(theta) as sqrt(1 - x^2), which keeps only a few digits near
 * the poles (4e-5 relative at theta = 1e-6).
 *
 * The colatitude derivatives come from the order recurrence
 *   dP_l^k/dtheta = (1/2) [P_l^(k+1) - (l+k)(l-k+1) P_l^(k-1)]
 * applied once and twice, with P_l^(-k) = (-1)^k (l-k)!/(l+k)! P_l^k and P_l^k = 0 for |k| > l. Every term is a
 * smooth function of theta, so the poles are exact. (The chain rule through x, dP/dtheta = -sin(theta) dP/dx,
 * multiplies a derivative that is infinite at the poles for odd m by a zero there.)
 */

#include <cmath>
#include <cstdint>

#include "constants_.hpp"
#include "legendre_common_.hpp"

namespace tidalpy {

namespace detail {

// P_l^k(cos theta) for 0 <= k <= l, from the recurrence in degree.
inline double legendre_nonnegative_order(int degree_l, int order_k, double cos_t, double sin_t)
{
    double p_kk = 1.0;
    for (int i = 1; i <= order_k; ++i) {
        p_kk *= -static_cast<double>(2 * i - 1) * sin_t;
    }
    if (degree_l == order_k) { return p_kk; }
    double p_previous = p_kk;
    double p_current  = static_cast<double>(2 * order_k + 1) * cos_t * p_kk;
    for (int degree = order_k + 2; degree <= degree_l; ++degree) {
        const double p_next = (static_cast<double>(2 * degree - 1) * cos_t * p_current
                               - static_cast<double>(degree + order_k - 1) * p_previous)
                            / static_cast<double>(degree - order_k);
        p_previous = p_current;
        p_current  = p_next;
    }
    return p_current;
}

// P_l^k(cos theta) for any integer order k.
inline double legendre_any_order(int degree_l, int order_k, double cos_t, double sin_t)
{
    if (order_k > degree_l || order_k < -degree_l) { return 0.0; }
    const int abs_k = (order_k < 0) ? -order_k : order_k;
    const double value = legendre_nonnegative_order(degree_l, abs_k, cos_t, sin_t);
    if (order_k >= 0) { return value; }
    // P_l^(-k) = (-1)^k (l-k)!/(l+k)! P_l^k.
    double ratio = 1.0;
    for (int factor = degree_l - abs_k + 1; factor <= degree_l + abs_k; ++factor) {
        ratio /= static_cast<double>(factor);
    }
    return ((abs_k % 2 == 0) ? ratio : -ratio) * value;
}

// dP_l^k/dtheta for any integer order k, from the order recurrence.
inline double legendre_dtheta_any_order(int degree_l, int order_k, double cos_t, double sin_t)
{
    const double ladder = static_cast<double>(degree_l + order_k) * static_cast<double>(degree_l - order_k + 1);
    return 0.5 * (legendre_any_order(degree_l, order_k + 1, cos_t, sin_t)
                  - ladder * legendre_any_order(degree_l, order_k - 1, cos_t, sin_t));
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
    const double sin_t  = std::sin(colatitude);
    const double ladder = static_cast<double>(degree_l + order_m) * static_cast<double>(degree_l - order_m + 1);

    return c_LegendreValue {
        detail::legendre_any_order(degree_l, order_m, cos_t, sin_t),
        detail::legendre_dtheta_any_order(degree_l, order_m, cos_t, sin_t),
        0.5 * (detail::legendre_dtheta_any_order(degree_l, order_m + 1, cos_t, sin_t)
               - ladder * detail::legendre_dtheta_any_order(degree_l, order_m - 1, cos_t, sin_t))};
}

} // namespace tidalpy
