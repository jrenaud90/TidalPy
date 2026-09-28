#pragma once
/* Orbital rate equations from tidal dissipation. Rates only; the System class integrates them and sums
 * the two bodies of a dual-body dissipation system.
 *
 * The potential derivatives become disturbing-function derivatives through the reduced mass:
 *   dR/dX = -((M_target + M_host) / M_target) * dU/dX ,   X in {mean anomaly M, arg pericenter w}.
 * (The 1/M_host that removes one host-mass power relative to the heating is already carried in the
 * dU/dX values from the global tidal collapse.) Then, following Boue & Efroimsky (2019, CMDA) Eqs.
 * 116-117:
 *   da/dt = (2 / (n a)) dR/dM
 *   de/dt = (sqrt(1-e^2) / (n a^2 e)) ( sqrt(1-e^2) dR/dM - dR/dw )
 *   dn/dt = -(3/2)(n / a) da/dt               // Kepler's third law differentiated
 *
 * At small e the bracket of de/dt is O(e^2) while each term is O(1) (with a non-synchronous spin), so it is
 * evaluated as
 *   sqrt(1-e^2) dR/dM - dR/dw = -(e^2 / (1 + sqrt(1-e^2))) dR/dM + d(R)/d(M - w),
 * with d(R)/d(M - w) from the per-mode sum dU_dM_minus_dw of the collapse, which holds no cancellation.
 */

#include <cmath>

#include "constants_.hpp"

namespace tidalpy {

// The orbital and mass state the rate equations need.
struct c_OrbitState {
    double orbital_frequency = 0.0;  // mean motion n           [rad s-1]
    double semi_major_axis   = 0.0;  // a                       [m]
    double eccentricity      = 0.0;  // e                       [dimensionless]
    double target_mass       = 0.0;  // the dissipating body    [kg]
    double host_mass         = 0.0;  // the companion (host)    [kg]
};

struct c_OrbitDerivatives {
    double da_dt = 0.0;  // [m s-1]
    double de_dt = 0.0;  // [s-1]
    double dn_dt = 0.0;  // [rad s-2]
};

class c_OrbitSolver {
public:
    c_OrbitSolver() noexcept = default;

    // da/dt = (2 / (n a)) dR/dM.
    double calc_da_dt(const c_OrbitState& state, double dU_dM) const noexcept {
        const double n_a = state.orbital_frequency * state.semi_major_axis;
        if (std::abs(n_a) <= TidalPyConstants::d_EPS) {
            return TidalPyConstants::d_NAN;
        }
        const double dR_dM = this->calc_dR(state, dU_dM);
        return 2.0 * dR_dM / n_a;
    }

    // de/dt = (sqrt(1-e^2) / (n a^2 e)) ( sqrt(1-e^2) dR/dM - dR/dw ). Zero for a circular orbit, where
    // the 1/e term is indeterminate. dU_dM_minus_dw is the per-mode sum of dU_dM - dU_dw (see the header); NaN
    // forms it from the two sums instead, which loses precision as e^2 approaches the rounding of either.
    double calc_de_dt(
            const c_OrbitState& state,
            double dU_dM,
            double dU_dw,
            double dU_dM_minus_dw = TidalPyConstants::d_NAN) const noexcept {
        const double denom = state.orbital_frequency * state.semi_major_axis
                           * state.semi_major_axis * state.eccentricity;
        if (std::abs(denom) <= TidalPyConstants::d_EPS) {
            return 0.0;
        }
        const double e2       = state.eccentricity * state.eccentricity;
        const double ecc_term = std::sqrt(1.0 - e2);
        const double dR_dM    = this->calc_dR(state, dU_dM);
        const double dR_dMw   = std::isfinite(dU_dM_minus_dw)
            ? this->calc_dR(state, dU_dM_minus_dw)
            : (dR_dM - this->calc_dR(state, dU_dw));
        // sqrt(1-e^2) - 1 = -e^2 / (1 + sqrt(1-e^2)), with no cancellation.
        return (ecc_term / denom) * (-(e2 / (1.0 + ecc_term)) * dR_dM + dR_dMw);
    }

    // dn/dt = -(3/2)(n / a) da/dt, from Kepler's third law.
    double calc_dn_dt(
            double orbital_frequency,
            double semi_major_axis,
            double da_dt) const noexcept {
        if (std::abs(semi_major_axis) <= TidalPyConstants::d_EPS) {
            return TidalPyConstants::d_NAN;
        }
        return -1.5 * (orbital_frequency / semi_major_axis) * da_dt;
    }

    c_OrbitDerivatives calc_derivatives(
            const c_OrbitState& state,
            double dU_dM,
            double dU_dw,
            double dU_dM_minus_dw = TidalPyConstants::d_NAN) const noexcept {
        c_OrbitDerivatives out;
        out.da_dt = this->calc_da_dt(state, dU_dM);
        out.de_dt = this->calc_de_dt(state, dU_dM, dU_dw, dU_dM_minus_dw);
        out.dn_dt = this->calc_dn_dt(state.orbital_frequency, state.semi_major_axis, out.da_dt);
        return out;
    }

private:
    // dR/dX = -((M_target + M_host)/M_target) dU/dX, from the dissipating body's potential derivative.
    double calc_dR(const c_OrbitState& state, double dU_dX) const noexcept {
        if (std::abs(state.target_mass) <= TidalPyConstants::d_EPS) {
            return TidalPyConstants::d_NAN;
        }
        return -((state.target_mass + state.host_mass) / state.target_mass) * dU_dX;
    }
};

}  // namespace tidalpy
