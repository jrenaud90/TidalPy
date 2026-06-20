#pragma once
/*
 * tide_result_.hpp — lightweight result + config structs for the global (1D) tidal solve.
 *
 * These are split out from tide_collapse_.hpp (which pulls in the heavy global-potential
 * engine + eccentricity/obliquity tables) so the world class can store a tide config and
 * result without compiling those tables into every translation unit that includes the
 * world header. The orchestration that actually runs the global potential + collapse lives
 * in structures_x/worlds/world_tides_.hpp and is compiled into one extension only.
 *
 * All quantities MKS; frequencies in rad s-1; angles in radians.
 *
 * These structs live in the global namespace (matching the global-potential + collapse
 * code they pair with); the world class in namespace tidalpy references them unqualified.
 */

// -------------------------------------------------------------------------------
// c_TideConfig — the world's stored [tides] configuration (truncation + degree range).
// The dissipation model itself is held separately on the world (c_TideBase).
// -------------------------------------------------------------------------------
struct c_TideConfig {
    int min_degree_l            = 2;    // lowest tidal harmonic degree (>= 2)
    int max_degree_l            = 2;    // highest tidal harmonic degree (<= 10)
    int eccentricity_truncation = 6;    // eccentricity-function truncation level
    int obliquity_truncation    = 10;   // obliquity-function truncation (0=off, 2, 4, 10=general)
};

// -------------------------------------------------------------------------------
// c_TideSolveConfig — the per-call orbital/spin state for calc_tides. Passed in each
// solve (later filled by the System class); the world stays stateless w.r.t. the orbit.
// -------------------------------------------------------------------------------
struct c_TideSolveConfig {
    double orbital_frequency = 0.0;   // orbital mean motion n          [rad s-1]
    double spin_frequency    = 0.0;   // spin rate of the deformed body [rad s-1]
    double eccentricity      = 0.0;   // orbital eccentricity           [dimensionless]
    double obliquity         = 0.0;   // axial tilt                     [radians]
    double semi_major_axis   = 0.0;   // orbital semi-major axis        [m]
    double host_mass         = 0.0;   // mass of the tidal host         [kg]
};

// -------------------------------------------------------------------------------
// c_GlobalTideResult — the collapsed global tidal solution.
// -------------------------------------------------------------------------------
struct c_GlobalTideResult {
    double tidal_heating = 0.0;  // total global tidal heating                         [W]
    double dU_dM         = 0.0;  // potential derivative wrt mean anomaly              [J kg-1 rad-1]
    double dU_dw         = 0.0;  // potential derivative wrt argument of pericenter    [J kg-1 rad-1]
    double dU_dO         = 0.0;  // potential derivative wrt longitude of node         [J kg-1 rad-1]
    int num_modes        = 0;    // number of active (nonzero-frequency) modes summed
    int error_code       = 0;    // propagated from the potential solve
};
