#pragma once
/*
 * love_.hpp - c_LoveNumbers: container for the three complex tidal Love numbers.
 *
 *   k - potential Love number (tidal modification of the gravity field)
 *   h - radial displacement Love number (radial surface deformation)
 *   l - tangential displacement Love number (horizontal surface deformation)
 *
 * All three are dimensionless complex numbers whose imaginary part carries the dissipation at the
 * tidal forcing frequency. This header is the data container only; the radial solver computes them.
 */

#include <cmath>
#include <complex>

namespace tidalpy {

struct c_LoveNumbers {
    std::complex<double> k = {0.0, 0.0};  // potential Love number           [dimensionless]
    std::complex<double> h = {0.0, 0.0};  // radial displacement Love number [dimensionless]
    std::complex<double> l = {0.0, 0.0};  // tangential displacement Love number [dimensionless]

    c_LoveNumbers() noexcept = default;

    c_LoveNumbers(std::complex<double> k_in,
                  std::complex<double> h_in,
                  std::complex<double> l_in) noexcept
        : k(k_in), h(h_in), l(l_in) {}

    bool operator==(const c_LoveNumbers& o) const noexcept {
        return k == o.k && h == o.h && l == o.l;
    }
    bool operator!=(const c_LoveNumbers& o) const noexcept { return !(*this == o); }
};

// The factors the analytic lag models multiply a Love number by, shared by the tide models (tide_.hpp) and the
// quasi-homogeneous Love methods (love_method_.hpp). Each depends on the forcing frequency's magnitude only: the global
// collapse passes |omega| and carries a mode's sign in its coefficients, and a direct call with a negative frequency
// gets the same lag. The callers decide what a zero or unset parameter means.
//
// Constant phase lag: 1 - i / Q.
inline std::complex<double> c_constant_phase_lag(double fixed_q) noexcept {
    return {1.0, -1.0 / fixed_q};
}
// Constant time lag: 1 - i |omega| dt, divided by Q as well for the time lag with a quality factor (Q = 1 without).
inline std::complex<double> c_constant_time_lag(double frequency, double fixed_dt, double fixed_q = 1.0) noexcept {
    return {1.0, -std::abs(frequency) * fixed_dt / fixed_q};
}

} // namespace tidalpy
