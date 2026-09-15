#pragma once
/*
 * potential_point_.hpp - c_PotentialPointC: the tidal potential's angular factor and its derivatives at one point
 * for one mode, as complex phasor amplitudes. A small dependency-free value type shared by the 3D potential engine
 * (potential_3d_.hpp, which produces it) and the 3D strain/stress kernel (multilayer/kernel_.hpp, which consumes
 * it). Kept standalone so the kernel does not have to pull in the potential engine.
 */

#include <complex>

namespace tidalpy {

// The tidal potential angular factor U and its first and second theta (colatitude) / phi (longitude) derivatives
// at one point and one mode, carried as COMPLEX amplitudes with the mode's time factor e^{i omega t} pulled out
// (U(t) = Re[U_c e^{i omega t}]). U already carries the radial factor (its r^2 coefficient, evaluated at the surface
// radius) and the orbital amplitude, so it multiplies the radial strain coefficients directly in the kernel. With
// complex amplitudes the cycle average is exact, h_bar = (omega/2) Im(sigma_c : conj(eps_c)), needing no time grid,
// and the 90-degree phase between U (cos) and its phi derivatives (which bring a factor i*m) is carried by the
// imaginary parts.
struct c_PotentialPointC {
    c_PotentialPointC() {}
    c_PotentialPointC(
            std::complex<double> U_,
            std::complex<double> dU_dtheta_,
            std::complex<double> dU_dphi_,
            std::complex<double> d2U_dtheta2_,
            std::complex<double> d2U_dphi2_,
            std::complex<double> d2U_dtheta_dphi_) :
        U(U_),
        dU_dtheta(dU_dtheta_),
        dU_dphi(dU_dphi_),
        d2U_dtheta2(d2U_dtheta2_),
        d2U_dphi2(d2U_dphi2_),
        d2U_dtheta_dphi(d2U_dtheta_dphi_) {}

    std::complex<double> U               {0.0, 0.0};
    std::complex<double> dU_dtheta       {0.0, 0.0};
    std::complex<double> dU_dphi         {0.0, 0.0};
    std::complex<double> d2U_dtheta2     {0.0, 0.0};
    std::complex<double> d2U_dphi2       {0.0, 0.0};
    std::complex<double> d2U_dtheta_dphi {0.0, 0.0};
};

}  // namespace tidalpy
