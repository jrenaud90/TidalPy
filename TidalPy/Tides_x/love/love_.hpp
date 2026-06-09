#pragma once
/*
 * love_.hpp — c_LoveNumbers: container for the three complex tidal Love numbers.
 *
 * The three Love numbers k, h, l describe how a body deforms under an external
 * tidal potential:
 *   k — potential Love number (tidal modification of the gravity field)
 *   h — radial displacement Love number (radial surface deformation)
 *   l — tangential displacement Love number (horizontal surface deformation)
 *
 * All three are dimensionless complex numbers; the imaginary part represents
 * dissipation at the tidal forcing frequency.
 *
 * This header provides only the data container.  Computation of Love numbers
 * from material properties is handled by the radial solver (Phase 10+).
 */

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

} // namespace tidalpy
