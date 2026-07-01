#pragma once
/*
 * tidal_potential_base_.hpp - c_TidalPotentialBase: abstract base for TidalPy's tidal-potential
 * truncation models (the 2D colatitude/longitude/time potential used by the 3D stress/strain kernel).
 *
 * Inherits c_PhysicsBase (Utilities_x/classes_x/physics_base_.hpp), mirroring c_TideBase. A tidal
 * potential model evaluates, at a point (radius, colatitude, longitude, time) and a given orbital/spin
 * state, the active tidal modes: each mode's signed forcing frequency and the potential angular factor
 * U together with its first/second colatitude(theta)/longitude(phi) derivatives. The frequency response
 * (complex moduli + radial y-solution at the mode frequency) is applied downstream by the dissipation
 * model (c_RheologyTide), which sums each mode's stress/strain.
 *
 * The concrete models live in tidal_potential_.hpp:
 *   c_SyncLowEPotential  (alias "sync_low_e") - synchronous rotation, low eccentricity, no obliquity (1 mode).
 *   c_NSRModesPotential  (alias "nsr_modes")  - moderate eccentricity, non-synchronous rotation, no obliquity
 *                                               (up to 9 modes n, 2n, 3n, 2o+-kn).
 *
 * All quantities MKS; frequencies in rad s-1; angles in radians. calc_modes is const. Universal
 * constants (G, the spin-orbit mode threshold) are read from the TidalPy config singleton, not passed in.
 */

#include <array>
#include <cmath>
#include <string>
#include <cstdint>

#include "physics_base_.hpp"     // tidalpy::c_PhysicsBase
#include "constants_.hpp"        // TidalPyConstants
#include "potential_point_.hpp"  // tidalpy::c_PotentialPoint

namespace tidalpy {

// Orbital/spin state the tidal potential depends on. Filled per call by the world's tide solve (the world
// stays stateless w.r.t. the orbit), matching how c_TideSolveConfig carries the 1D state.
struct c_TidalPotentialState {
    double orbital_frequency = 0.0;   // n   [rad s-1]
    double spin_frequency    = 0.0;   // o   [rad s-1]
    double eccentricity      = 0.0;
    double obliquity         = 0.0;   // [rad] (unused by the no-obliquity truncations)
    double host_mass         = 0.0;   // [kg]
    double semi_major_axis   = 0.0;   // [m]
};

// Maximum number of modes any currently-supported l = 2 truncation produces (extend when obliquity modes
// are added). The mode set is a fixed-capacity value type so it can be returned without heap allocation.
constexpr int C_MAX_TIDAL_POTENTIAL_MODES = 9;

// The active tidal modes at a point: signed frequency + potential angular factor per mode.
struct c_TidalPotentialModeSet {
    int num_modes = 0;
    double mode_frequency[C_MAX_TIDAL_POTENTIAL_MODES] = {0.0};   // signed [rad s-1]
    c_PotentialPoint potential[C_MAX_TIDAL_POTENTIAL_MODES];
};

class c_TidalPotentialBase : public c_PhysicsBase {
public:
    c_TidalPotentialBase() = default;

    explicit c_TidalPotentialBase(const std::string& model_name) : c_PhysicsBase(model_name) {}

    ~c_TidalPotentialBase() override = default;

    // Evaluate the active tidal modes (signed frequency + potential angular factor) at one point.
    //
    // Assumptions
    // -----------
    // - The potential's r^2 coefficient uses `radius` directly; callers evaluating the 3D kernel pass the
    //   surface radius (the radial dependence lives in the y-functions), matching the legacy collapse.
    // - Universal constants (G, the spin-orbit mode threshold) are read from the TidalPy config singleton.
    virtual c_TidalPotentialModeSet calc_modes(
            const c_TidalPotentialState& state,
            double radius,
            double colatitude,
            double longitude,
            double time) const = 0;

    // Number of modes this truncation produces (independent of position/time).
    virtual int num_modes() const = 0;
};

} // namespace tidalpy
