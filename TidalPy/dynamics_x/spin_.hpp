#pragma once
/* Rotational (spin) dynamics of a tidally interacting body. Rates only; the System class integrates them.
 *
 * dU/dO comes from the global tidal solve (c_GlobalTideResult.dU_dO) but is taken as a plain scalar, so
 * this module does not depend on the Tides_x headers.
 *
 * Reference: Ferraz-Mello et al. (2008) for the spin-rate torque.
 */

#include <cmath>
#include <stdexcept>
#include <string>

#include "constants_.hpp"

namespace tidalpy {

// Largest physical moment-of-inertia factor C / (M R^2): all of the mass in a thin surface shell.
constexpr double d_SPIN_MOI_FACTOR_MAX = 2.0 / 3.0;

struct c_SpinConfig {
    // Conventional f = C / (M R^2): 0.4 for a uniform sphere, less for a centrally condensed body
    // (0.3307 for the Earth), at most 2/3.
    double moment_of_inertia_factor = 0.4;
};

// The upper bound is the thin hollow shell, the most outwardly concentrated body with non-negative
// density. Rejecting anything above it catches a 1.0 from the ratio-to-a-uniform-sphere convention,
// which would otherwise give a moment of inertia 2.5 times too large.
inline void c_validate_moi_factor(double factor) {
    if (!std::isfinite(factor) || factor <= 0.0 || factor > d_SPIN_MOI_FACTOR_MAX) {
        throw std::invalid_argument(
            "TidalPy: moment_of_inertia_factor must be the conventional C / (M R^2), within (0, 2/3]; got "
            + std::to_string(factor) + " (a uniform sphere is 0.4)");
    }
}

class c_Spin {
public:
    c_Spin() noexcept = default;

    explicit c_Spin(const c_SpinConfig& config)
        : p_config(config)
    {
        c_validate_moi_factor(config.moment_of_inertia_factor);
    }

    const c_SpinConfig& get_config() const noexcept { return this->p_config; }

    // I = f M R^2.
    double calc_moment_of_inertia(
            double mass,
            double radius) const noexcept {
        return this->p_config.moment_of_inertia_factor * mass * radius * radius;
    }

    // dspin/dt = M_host * dU/dO / I (Ferraz-Mello et al. 2008); the tidal polar torque is M_host * dU/dO.
    // Assumes the spin axis is aligned with the orbit normal, the axis the polar torque drives.
    double calc_dspin_dt(
            double host_mass,
            double dU_dO,
            double moment_of_inertia) const noexcept {
        if (std::abs(moment_of_inertia) <= TidalPyConstants::d_EPS) {
            return TidalPyConstants::d_NAN;
        }
        return host_mass * dU_dO / moment_of_inertia;
    }

    // A tidally locked body spins at the orbital mean motion.
    double calc_synchronous_spin(double orbital_frequency) const noexcept {
        return orbital_frequency;
    }

private:
    c_SpinConfig p_config {};
};

}  // namespace tidalpy
