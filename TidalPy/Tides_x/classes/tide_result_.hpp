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

#include <cstddef>
#include <vector>

// -------------------------------------------------------------------------------
// c_TideConfig — the world's stored [tides] configuration (truncation + degree range).
// The dissipation model itself is held separately on the world (c_TideBase).
// -------------------------------------------------------------------------------
struct c_TideConfig {
    int min_degree_l            = 2;    // lowest tidal harmonic degree (>= 2)
    int max_degree_l            = 2;    // highest tidal harmonic degree (<= 10)
    int eccentricity_truncation = 6;    // eccentricity-function truncation level
    int obliquity_truncation    = 10;   // obliquity-function truncation (0=off, 2, 4, 10=general)
    // Width [decades] of the log-Gaussian bell used by the tidal_timescale layer scale method
    // (scale = exp(-0.5*(log10(maxwell_time/forcing_period)/width)^2)).
    double tidal_timescale_width_decades = 1.0;
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

// -------------------------------------------------------------------------------
// 3D tidal-heating collapse (flavor flags + result).
//
// The 3D tidal heating can be produced as a full grid over (radius, colatitude, longitude[, time]) or
// reduced (integrated) along any of the spatial dimensions:
//   * orbit_averaged = true  -> the secular (cycle-averaged) volumetric heating density h_bar [W m-3].
//     It is longitude- and time-independent (the per-mode e^{i m phi} cancels), so the longitude axis
//     is constant and the longitude integral is analytic (a factor 2*pi). No time axis.
//   * orbit_averaged = false -> the instantaneous mechanical power density sigma_ij(t) eps_dot_ij(t)
//     [W m-3], evaluated at each user-supplied time. It depends on longitude and time and orbit-averages
//     back to h_bar. A time axis is present.
//
// Reduction convention (marginal densities): if any spatial axis is summed, each surviving spatial axis
// carries its Jacobian (r^2 for radius, sin theta for colatitude, 1 for longitude), and each SUMMED axis
// is integrated with its Jacobian and quadrature (colatitude: Gauss-Legendre in cos theta, which absorbs
// the sin theta weight; radius: per-layer trapezoid; longitude: 2*pi analytic when averaged, trapezoid
// over [0, 2*pi) when instantaneous). A plain integral over the surviving axes then recovers the total.
// If NO axis is summed the output is the raw density (it matches the point/scalar path). Whole-planet and
// per-layer totals are produced only when all three spatial axes are summed (per-time when instantaneous).
//
// Non-summed spatial axes use the user-supplied radii / colatitudes / longitudes arrays; summed axes use
// internal integration grids. The time axis (instantaneous only) always uses the user-supplied times.
// -------------------------------------------------------------------------------
struct c_Heating3DCollapseConfig {
    bool orbit_averaged   = true;   // true: secular density; false: instantaneous sigma:eps_dot vs time
    bool latitude_summed  = false;  // integrate over colatitude (Gauss-Legendre, sin theta weight)
    bool longitude_summed = false;  // integrate over longitude (2*pi analytic when averaged; else trapezoid)
    bool radial_summed    = false;  // integrate over radius (per-layer trapezoid, r^2 weight)
    // Integration resolutions. latitude_nodes / radial_slices tuned empirically: for a homogeneous
    // degree-2 body the collapsed total is within ~1% of the 1D global heating by ~4 Gauss-Legendre
    // colatitude nodes and ~8 radial slices per layer; these defaults add margin for higher degree l and
    // layered bodies (Gauss-Legendre is cheap). Raise them for many-layer or high-degree configs.
    int  latitude_nodes   = 16;     // Gauss-Legendre order for the colatitude integral
    int  longitude_nodes  = 64;     // trapezoid nodes for the instantaneous longitude integral
    int  radial_slices    = 16;     // trapezoid slices per layer for the radial integral
};

struct c_Heating3DCollapsed {
    // Primary output, flattened row-major over the surviving axes in the fixed order
    // [radius, colatitude, longitude, time]; a spatial axis is dropped when summed, and the time axis is
    // present only when not orbit-averaged. `shape` lists the surviving-axis lengths in that order.
    std::vector<double> values;
    std::vector<std::size_t> shape;
    std::vector<double> radii;          // radius axis actually used         [m]
    std::vector<double> colatitudes;    // colatitude axis actually used     [rad]
    std::vector<double> longitudes;     // longitude axis actually used      [rad]
    std::vector<double> times;          // time axis (empty when orbit-averaged) [s]
    // Per-layer totals [W], only when all three spatial axes are summed; flattened [n_layers, n_times]
    // (n_times = 1 when orbit-averaged).
    std::vector<double> layer_totals;
    std::size_t n_layers = 0;
    std::size_t n_times  = 1;
    bool all_spatial_summed = false;    // whole-planet + per-layer totals are populated
};
