#pragma once
/* dy1/dr of a solid layer (Takeuchi and Saito 1972), the one place both the radial ODEs (odes_.hpp) and the 3D strain
 * at a point (Tides/multilayer/strain_radial_.hpp) take it from, so the two cannot drift apart. y1_y3_term is
 * 2 y1 - l(l+1) y3 and r_inverse is 1 / r. Light on includes so the strain code can use it without the ODE machinery.
 */

#include <complex>

// Compressible: (1 / (lambda + 2 mu)) [ y2 - (lambda / r)(2 y1 - l(l+1) y3) ].
inline std::complex<double> c_solid_dy1_compressible(
        const std::complex<double>& y1_y3_term,
        const std::complex<double>& y2,
        const std::complex<double>& lame,
        const std::complex<double>& lame_2mu,
        double r_inverse) noexcept {
    // In the radial ODEs' own order of operations, so their results do not move.
    return (1.0 / lame_2mu) * (y1_y3_term * -lame * r_inverse + y2);
}

// Incompressible: -(2 y1 - l(l+1) y3) / r.
inline std::complex<double> c_solid_dy1_incompressible(const std::complex<double>& y1_y3_term, double r_inverse) noexcept {
    return -y1_y3_term * r_inverse;
}
