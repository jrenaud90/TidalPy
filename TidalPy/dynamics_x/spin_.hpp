#pragma once
/*
 * spin_.hpp - c_Spin: rotational (spin) dynamics of a tidally interacting body.
 *
 * A world's spin evolves under the tidal torque exerted by its host. This class provides the
 * rate quantities only (no time integration; the System class integrates them):
 *   - the moment of inertia of the body,
 *   - the tidal spin-rate change dspin/dt from the tidal potential derivative dU/dO, and
 *   - the synchronous spin rate.
 *
 * The tidal potential derivative dU/dO is produced by the global tidal solve
 * (c_GlobalTideResult.dU_dO, from world.calc_tides); it is passed in here as a plain scalar so this
 * module does not depend on the Tides_x headers.
 *
 * References: Ferraz-Mello et al. (2008) for the spin-rate torque; the moment-of-inertia form is the
 * standard uniform-density spherical-shell result scaled by a dimensionless structure factor.
 *
 * All quantities MKS: masses in kg, radii in m, frequencies in rad s-1, moment of inertia in kg m2,
 * dU/dO in J kg-1 rad-1, dspin/dt in rad s-2.
 */

#include <cmath>

#include "constants_.hpp"   // TidalPyConstants::d_EPS, d_NAN

namespace tidalpy {

// -------------------------------------------------------------------------------
// c_SpinConfig - the (few) configurable properties of the spin model.
// -------------------------------------------------------------------------------
struct c_SpinConfig {
    // Dimensionless moment-of-inertia factor C / (M R^2) relative to the uniform-density value: 1 for a
    // uniform sphere, smaller for a centrally condensed body (~0.33 for the Earth). Scales the ideal MoI.
    double moment_of_inertia_factor = 1.0;
};

// -------------------------------------------------------------------------------
// c_Spin - spin-dynamics calculator (rates only).
// -------------------------------------------------------------------------------
class c_Spin {
public:
    c_Spin() noexcept = default;
    explicit c_Spin(const c_SpinConfig& config) noexcept
        : p_config(config) {}

    const c_SpinConfig& get_config() const noexcept { return this->p_config; }

    // Moment of inertia [kg m2] of a uniform-density spherical shell (a solid sphere when radius_inner = 0),
    // scaled by the configured structure factor:
    //   I = factor * (2/5) M (R_outer^5 - R_inner^5) / (R_outer^3 - R_inner^3).
    // Returns NaN for a degenerate (zero-thickness) shell.
    //
    // Assumptions
    // -----------
    // Uniform density within the shell; the structure factor absorbs any real radial density variation.
    double calc_moment_of_inertia(
            double mass,
            double radius_outer,
            double radius_inner = 0.0) const noexcept {
        const double r_outer3 = radius_outer * radius_outer * radius_outer;
        const double r_inner3 = radius_inner * radius_inner * radius_inner;
        const double volume_term = r_outer3 - r_inner3;
        if (std::abs(volume_term) <= TidalPyConstants::d_EPS) {
            return TidalPyConstants::d_NAN;
        }
        const double r_outer5 = r_outer3 * radius_outer * radius_outer;
        const double r_inner5 = r_inner3 * radius_inner * radius_inner;
        return this->p_config.moment_of_inertia_factor * 0.4 * mass * (r_outer5 - r_inner5) / volume_term;
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

    // Synchronous spin rate [rad s-1]: a synchronously (tidally locked) rotating body spins at the
    // orbital mean motion, so the synchronous spin equals the orbital frequency.
    double calc_synchronous_spin(double orbital_frequency) const noexcept {
        return orbital_frequency;
    }

private:
    c_SpinConfig p_config {};
};

}  // namespace tidalpy
