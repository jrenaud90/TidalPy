#pragma once
/*
 * spin_.hpp - c_Spin: rotational (spin) dynamics of a tidally interacting body.
 *
 * Provides the body's moment of inertia [kg m2], the tidal spin-rate change dspin/dt [rad s-2] from
 * the tidal potential derivative dU/dO [J kg-1 rad-1], and the synchronous spin rate [rad s-1].
 * Rates only; the System class integrates them. dU/dO comes from the global tidal solve
 * (c_GlobalTideResult.dU_dO) and is passed in as a plain scalar so this module does not depend on
 * the Tides_x headers.
 *
 * The moment of inertia is I = f M R^2 with the conventional dimensionless factor f = C / (M R^2).
 * Reference: Ferraz-Mello et al. (2008) for the spin-rate torque.
 */

#include <cmath>
#include <stdexcept>
#include <string>

#include "constants_.hpp"   // TidalPyConstants::d_EPS, d_NAN

namespace tidalpy {

// Largest physical moment-of-inertia factor C / (M R^2): all of the mass in a thin surface shell.
constexpr double d_SPIN_MOI_FACTOR_MAX = 2.0 / 3.0;

// -------------------------------------------------------------------------------
// c_SpinConfig - the (few) configurable properties of the spin model.
// -------------------------------------------------------------------------------
struct c_SpinConfig {
    // Conventional dimensionless moment-of-inertia factor f = C / (M R^2): 0.4 for a uniform sphere, smaller for
    // a centrally condensed body (0.3307 for the Earth), and at most 2/3 (a thin hollow shell).
    double moment_of_inertia_factor = 0.4;
};

// -------------------------------------------------------------------------------
// c_validate_moi_factor - reject an unphysical moment-of-inertia factor.
//
// A factor must be finite and within (0, 2/3]. The upper bound is the thin hollow shell, the most outwardly
// concentrated body with non-negative density. A value above it, such as 1.0 from the ratio-to-a-uniform-sphere
// convention, is rejected rather than silently producing a moment of inertia 2.5 times too large.
// Throws std::invalid_argument (surfaced as ValueError in Cython).
// -------------------------------------------------------------------------------
inline void c_validate_moi_factor(double factor) {
    if (!std::isfinite(factor) || factor <= 0.0 || factor > d_SPIN_MOI_FACTOR_MAX) {
        throw std::invalid_argument(
            "TidalPy: moment_of_inertia_factor must be the conventional C / (M R^2), within (0, 2/3]; got "
            + std::to_string(factor) + " (a uniform sphere is 0.4)");
    }
}

// -------------------------------------------------------------------------------
// c_Spin - spin-dynamics calculator (rates only).
// -------------------------------------------------------------------------------
class c_Spin {
public:
    c_Spin() noexcept = default;

    // Throws std::invalid_argument if the configured moment-of-inertia factor is unphysical.
    explicit c_Spin(const c_SpinConfig& config)
        : p_config(config)
    {
        c_validate_moi_factor(config.moment_of_inertia_factor);
    }

    const c_SpinConfig& get_config() const noexcept { return this->p_config; }

    // Moment of inertia [kg m2] from the conventional structure factor: I = f M R^2.
    double calc_moment_of_inertia(
            double mass,
            double radius) const noexcept {
        return this->p_config.moment_of_inertia_factor * mass * radius * radius;
    }

    // Tidal spin-rate change [rad s-2]: dspin/dt = M_host * dU/dO / I (Ferraz-Mello et al. 2008), where
    // the tidal polar torque is M_host * dU/dO. dU_dO is the c_GlobalTideResult value [J kg-1 rad-1].
    // Returns NaN for a non-positive moment of inertia.
    //
    // Assumptions
    // -----------
    // The spin axis is aligned with the orbit normal (the polar torque drives spin along that axis).
    double calc_dspin_dt(
            double host_mass,
            double dU_dO,
            double moment_of_inertia) const noexcept {
        if (std::abs(moment_of_inertia) <= TidalPyConstants::d_EPS) {
            return TidalPyConstants::d_NAN;
        }
        return host_mass * dU_dO / moment_of_inertia;
    }

    // Synchronous spin rate [rad s-1]: a tidally locked body spins at the orbital mean motion.
    double calc_synchronous_spin(double orbital_frequency) const noexcept {
        return orbital_frequency;
    }

private:
    c_SpinConfig p_config {};
};

}  // namespace tidalpy
