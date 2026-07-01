#pragma once
/*
 * potential_point_.hpp - c_PotentialPoint: the tidal potential's angular factor and its derivatives at
 * one point and one mode. A tiny dependency-free value type shared by the tidal-potential model classes
 * (tidal_potential_base_.hpp, which produce it) and the 3D strain/stress kernel (multilayer/kernel_.hpp,
 * which consumes it). Kept standalone so the kernel does not have to pull in c_PhysicsBase.
 */

namespace tidalpy {

// The tidal potential's angular factor U and its first/second theta(colatitude)/phi(longitude)
// derivatives at one point and one mode (all real - the potential is cos/sin in the mode phase). U
// already carries the radial factor (its r^2 coefficient, evaluated at the surface radius) and the
// orbital amplitude so it can be multiplied directly by the radial strain coefficients in the kernel.
struct c_PotentialPoint {
    c_PotentialPoint() {}
    c_PotentialPoint(
            double U_,
            double dU_dtheta_,
            double dU_dphi_,
            double d2U_dtheta2_,
            double d2U_dphi2_,
            double d2U_dtheta_dphi_) :
        U(U_),
        dU_dtheta(dU_dtheta_),
        dU_dphi(dU_dphi_),
        d2U_dtheta2(d2U_dtheta2_),
        d2U_dphi2(d2U_dphi2_),
        d2U_dtheta_dphi(d2U_dtheta_dphi_) {}

    double U               {0.0};
    double dU_dtheta       {0.0};
    double dU_dphi         {0.0};
    double d2U_dtheta2     {0.0};
    double d2U_dphi2       {0.0};
    double d2U_dtheta_dphi {0.0};
};

}  // namespace tidalpy
