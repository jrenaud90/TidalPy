#pragma once
/*
 * orbit_solver_.hpp - c_OrbitSolver: orbital rate equations from tidal dissipation.
 *
 * Given the tidal-potential derivatives of a dissipating body (dU/dM, dU/dw, produced by
 * world.calc_tides) and the orbital state, this returns the instantaneous rates of change of the
 * semi-major axis (da/dt), eccentricity (de/dt), and mean motion (dn/dt). Rates only, no time
 * integration. For a dual-body dissipation system each body's rates are additive in the disturbing-function
 * derivatives; the System class sums them.
 *
 * The tidal-potential derivatives are converted to derivatives of the disturbing function R via the
 * reduced mass: for the dissipating (target) body with companion (host),
 *   dR/dX = -((M_target + M_host) / M_target) * dU/dX ,   X in {mean anomaly M, arg pericenter w}.
 * (The 1/M_host that removes one host-mass power relative to the heating is already carried in the
 * dU/dX values from the global tidal collapse.) Then, following Boue & Efroimsky (2019, CMDA) Eqs.
 * 116-117:
 *   da/dt = (2 / (n a)) dR/dM
 *   de/dt = (sqrt(1-e^2) / (n a^2 e)) ( sqrt(1-e^2) dR/dM - dR/dw )
 *   dn/dt = -(3/2)(n / a) da/dt               # (Kepler's third law differentiated)
 *
 * All quantities MKS: n [rad s-1], a [m], e [-], masses [kg], dU/dX [J kg-1 rad-1], da/dt [m s-1],
 * de/dt [s-1], dn/dt [rad s-2].
 */

#include <cmath>

#include "constants_.hpp"   // TidalPyConstants::d_EPS, d_NAN

namespace tidalpy {

// -------------------------------------------------------------------------------
// c_OrbitState - the orbital + mass state the rate equations need.
// -------------------------------------------------------------------------------
struct c_OrbitState {
    double orbital_frequency = 0.0;  // mean motion n           [rad s-1]
    double semi_major_axis   = 0.0;  // a                       [m]
    double eccentricity      = 0.0;  // e                       [dimensionless]
    double target_mass       = 0.0;  // the dissipating body    [kg]
    double host_mass         = 0.0;  // the companion (host)    [kg]
};

// -------------------------------------------------------------------------------
// c_OrbitDerivatives - the instantaneous orbital rates.
// -------------------------------------------------------------------------------
struct c_OrbitDerivatives {
    double da_dt = 0.0;  // [m s-1]
    double de_dt = 0.0;  // [s-1]
    double dn_dt = 0.0;  // [rad s-2]
};

// -------------------------------------------------------------------------------
// c_OrbitSolver - orbital rate calculator (rates only).
// -------------------------------------------------------------------------------
class c_OrbitSolver {
public:
    c_OrbitSolver() noexcept = default;

    // Semi-major-axis rate [m s-1]: da/dt = (2 / (n a)) dR/dM.
    double calc_da_dt(const c_OrbitState& state, double dU_dM) const noexcept {
        const double n_a = state.orbital_frequency * state.semi_major_axis;
        if (std::abs(n_a) <= TidalPyConstants::d_EPS) {
            return TidalPyConstants::d_NAN;
        }
        const double dR_dM = this->calc_dR(state, dU_dM);
        return 2.0 * dR_dM / n_a;
    }

    // Eccentricity rate [s-1]: de/dt = (sqrt(1-e^2) / (n a^2 e)) ( sqrt(1-e^2) dR/dM - dR/dw ).
    // Returns 0 for a circular (or degenerate) orbit, where the 1/e term is indeterminate.
    double calc_de_dt(const c_OrbitState& state, double dU_dM, double dU_dw) const noexcept {
        const double denom = state.orbital_frequency * state.semi_major_axis
                           * state.semi_major_axis * state.eccentricity;
        if (std::abs(denom) <= TidalPyConstants::d_EPS) {
            return 0.0;
        }
        const double dR_dM = this->calc_dR(state, dU_dM);
        const double dR_dw = this->calc_dR(state, dU_dw);
        const double ecc_term = std::sqrt(1.0 - state.eccentricity * state.eccentricity);
        return (ecc_term / denom) * (ecc_term * dR_dM - dR_dw);
    }

    // Mean-motion rate [rad s-2]: dn/dt = -(3/2)(n / a) da/dt (from Kepler's third law).
    double calc_dn_dt(
            double orbital_frequency,
            double semi_major_axis,
            double da_dt) const noexcept {
        if (std::abs(semi_major_axis) <= TidalPyConstants::d_EPS) {
            return TidalPyConstants::d_NAN;
        }
        return -1.5 * (orbital_frequency / semi_major_axis) * da_dt;
    }

    // All three rates in one call.
    c_OrbitDerivatives calc_derivatives(
            const c_OrbitState& state,
            double dU_dM,
            double dU_dw) const noexcept {
        c_OrbitDerivatives out;
        out.da_dt = this->calc_da_dt(state, dU_dM);
        out.de_dt = this->calc_de_dt(state, dU_dM, dU_dw);
        out.dn_dt = this->calc_dn_dt(state.orbital_frequency, state.semi_major_axis, out.da_dt);
        return out;
    }

private:
    // Disturbing-function derivative dR/dX = -((M_target + M_host)/M_target) dU/dX from a potential
    // derivative dU/dX of the dissipating body.
    double calc_dR(const c_OrbitState& state, double dU_dX) const noexcept {
        if (std::abs(state.target_mass) <= TidalPyConstants::d_EPS) {
            return TidalPyConstants::d_NAN;
        }
        return -((state.target_mass + state.host_mass) / state.target_mass) * dU_dX;
    }
};

}  // namespace tidalpy
