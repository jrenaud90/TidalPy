#pragma once
/* Orbital-element conversions by Kepler's third law for a two-body orbit, n^2 a^3 = mu with mu = G (M_host + M_target).
 *
 * The caller forms mu, reading G at call time, so a reinitialized config is honored. Neither function checks its
 * inputs.
 */

#include <cmath>

namespace tidalpy {

// Mean motion n = sqrt(mu / a^3) [rad s-1] from a semi-major axis [m] and gravitational parameter mu [m3 s-2].
inline double c_semi_a2orbital_motion(double semi_major_axis, double gravitational_parameter) noexcept {
    return std::sqrt(gravitational_parameter / (semi_major_axis * semi_major_axis * semi_major_axis));
}

// Semi-major axis a = (mu / n^2)^(1/3) [m] from a mean motion [rad s-1] and gravitational parameter mu [m3 s-2].
inline double c_orbital_motion2semi_a(double orbital_motion, double gravitational_parameter) noexcept {
    return std::cbrt(gravitational_parameter / (orbital_motion * orbital_motion));
}

} // namespace tidalpy
