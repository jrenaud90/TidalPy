#pragma once
/*
 * tide_result_.hpp - result and config structs for the global (1D) tidal solve.
 *
 * Split out from tide_collapse_.hpp (which pulls in the global-potential engine and the
 * eccentricity/obliquity tables) so the world class can store a tide config and result without
 * compiling those tables into every translation unit that includes the world header. The
 * orchestration that runs the potential and the collapse lives in structures_x/worlds/world_tides_.hpp.
 *
 * All quantities MKS; frequencies in rad s-1; angles in radians. These structs live in the global
 * namespace, matching the collapse code they pair with.
 */

#include <cstddef>
#include <limits>
#include <vector>

#include "constants_.hpp"   // TidalPyConstants::d_PI (colatitude band default)

// c_TideConfig: the world's stored [tides] configuration. The dissipation model itself is held
// separately on the world (c_TideBase).
struct c_TideConfig {
    int min_degree_l            = 2;    // lowest tidal harmonic degree (>= 2)
    int max_degree_l            = 2;    // highest tidal harmonic degree (<= 10)
    int eccentricity_truncation = 3;    // eccentricity-function truncation level (e^3 in G, e^6 in G^2)
    int obliquity_truncation    = 10;   // obliquity-function truncation (0=off, 2, 4, 10=general)
    // Width [decades] of the log-Gaussian bell used by the tidal_timescale layer scale method
    // (scale = exp(-0.5*(log10(maxwell_time/forcing_period)/width)^2)).
    double tidal_timescale_width_decades = 1.0;
    // How the world obtains its Love numbers when the tide model asks for them (c_LoveMethod as an int:
    // 0 radial_solver, 1 propagation_matrix, 2 homogeneous, 3 cpl, 4 ctl, 5 laterally_inhomogeneous).
    int love_method = 0;
    // Fixed quality factor / time lag [s] for the cpl / ctl Love methods (NaN: take them from the tide model).
    double love_fixed_q  = std::numeric_limits<double>::quiet_NaN();
    double love_fixed_dt = std::numeric_limits<double>::quiet_NaN();
};

// c_TideSolveConfig: the per-call orbital and spin state for calc_tides. The world stays stateless
// with respect to the orbit.
struct c_TideSolveConfig {
    double orbital_frequency = 0.0;   // orbital mean motion n          [rad s-1]
    double spin_frequency    = 0.0;   // spin rate of the deformed body [rad s-1]
    double eccentricity      = 0.0;   // orbital eccentricity           [dimensionless]
    double obliquity         = 0.0;   // axial tilt                     [radians]
    double semi_major_axis   = 0.0;   // orbital semi-major axis        [m]
    double host_mass         = 0.0;   // mass of the tidal host         [kg]
};

// c_GlobalTideResult: the collapsed global tidal solution.
struct c_GlobalTideResult {
    double tidal_heating = 0.0;  // total global tidal heating                         [W]
    double dU_dM         = 0.0;  // potential derivative wrt mean anomaly              [J kg-1 rad-1]
    double dU_dw         = 0.0;  // potential derivative wrt argument of pericenter    [J kg-1 rad-1]
    double dU_dO         = 0.0;  // potential derivative wrt longitude of node         [J kg-1 rad-1]
    int num_modes        = 0;    // number of active (nonzero-frequency) modes summed
    int error_code       = 0;    // propagated from the potential solve
};

// -------------------------------------------------------------------------------
// 3D tidal-heating collapse (flavor flags and result).
//
// orbit_averaged = true gives the secular volumetric heating density h_bar [W m-3]: the time average
// of the instantaneous power at each point. It carries no time axis and still depends on longitude
// wherever waves at one frequency have different azimuthal structure, as they do for a synchronously
// rotating body; the longitude integral is analytic (2*pi times the longitude mean, which drops those
// cross terms). orbit_averaged = false gives the instantaneous power density sigma_ij eps_dot_ij
// [W m-3] at each supplied time, which time-averages to h_bar.
//
// Reduction convention (marginal densities): when any spatial axis is summed, each surviving spatial
// axis carries its Jacobian (r^2 for radius, sin theta for colatitude, 1 for longitude) and each
// summed axis is integrated with its Jacobian and quadrature (colatitude: Gauss-Legendre in cos theta,
// which absorbs the sin theta weight; radius: per-layer trapezoid; longitude: 2*pi analytic when
// averaged, trapezoid over [0, 2*pi) when instantaneous), so a plain integral over the surviving axes
// recovers the total. With no axis summed the output is the raw density. Whole-planet and per-layer
// totals appear only when all three spatial axes are summed.
//
// Non-summed spatial axes use the supplied radii, colatitudes, and longitudes; summed axes use
// internal integration grids. The time axis always uses the supplied times.
// -------------------------------------------------------------------------------
struct c_Heating3DCollapseConfig {
    bool orbit_averaged   = true;   // true: secular density; false: instantaneous sigma:eps_dot vs time
    bool latitude_summed  = false;  // integrate over colatitude (Gauss-Legendre, sin theta weight)
    bool longitude_summed = false;  // integrate over longitude (2*pi analytic when averaged; else trapezoid)
    bool radial_summed    = false;  // integrate over radius (Gauss-Legendre inside each layer, r^2 weight)
    // Integration resolutions. Both integrals use Gauss-Legendre nodes and the radial nodes stay inside
    // each layer, so the collapsed total converges quickly to the 1D global heating. Raise them if a
    // refined call still moves the total.
    int  latitude_nodes   = 16;     // Gauss-Legendre order for the colatitude integral
    int  longitude_nodes  = 64;     // trapezoid nodes for the instantaneous longitude integral
    int  radial_slices    = 16;     // Gauss-Legendre nodes per layer for the radial integral
    int  num_threads      = 1;      // threads for the per-point evaluation after the radial solves
    // When latitude_summed for the secular heating, do the colatitude integral with the precomputed
    // analytic angular Gram table (exact, no theta grid) instead of the quadrature above.
    bool latitude_analytic = true;
    // Colatitude band [rad] for the latitude integral. The analytic Gram table is full-sphere only, so a
    // band narrower than [0, pi] always falls back on the Gauss-Legendre quadrature.
    double colatitude_min = 0.0;
    double colatitude_max = TidalPyConstants::d_PI;
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
