#pragma once
/*
 * eccentricity_accuracy_.hpp - how far each eccentricity truncation level can be trusted, and which to choose.
 *
 * The limits are measured (2026-09-24 for degrees 2 and 3, 2026-09-29 for degrees 4 to 10) with TidalPy's tables
 * against the exact heating, from exact Hansen coefficients summed over |q| <= 4000: for each level, degree, and
 * tolerance, the largest eccentricity (on a 0.005 grid) at which that degree's heating stays within the tolerance,
 * worst case over constant-phase-lag, constant-time-lag, and Maxwell tides (relaxation times of 1/3 and 1/10 of the
 * inverse mean motion) at spin rates of 0.5, 1, and 2.3 times the mean motion. Every level errs low. Higher degrees
 * generally lose accuracy at a lower eccentricity, so a solve through max_degree_l takes the tightest limit of degrees
 * 2 to max_degree_l. Level 50 is also capped where its mode sum cancels by more than about 6e5, which multiplies any
 * error in the Love numbers: 0.78 at degree 2, 0.75 at degree 3, falling to 0.555 at degree 8 (the truncation limits
 * are lower at degrees 9 and 10). See Documentation/Tides/Eccentricity.md.
 *
 * Standard-library only, so any extension can include it.
 */

#include "../truncation_accuracy_.hpp"

// The truncation code of the exact eccentricity functions: Hansen coefficients from the exact Kepler orbit, the mode
// range chosen by a heating tolerance (eccentricity_exact_.hpp).
inline constexpr int C_ECCENTRICITY_EXACT = -1;

// Default tail tolerance of the exact eccentricity functions; the configuration's [tides]
// eccentricity_exact_tolerance supplies the value a built world uses.
inline constexpr double C_ECCENTRICITY_EXACT_TOLERANCE = 1.0e-4;

inline constexpr c_TruncationAccuracyRow C_ECCENTRICITY_ACCURACY[7] = {
    // One row per degree (l = 2 to 10), one column per tolerance (1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 1e-1).
    { 2, {{0.0,   0.0,   0.0,   0.005, 0.02,  0.075},   // l = 2
          {0.0,   0.0,   0.0,   0.005, 0.015, 0.06 },   // l = 3
          {0.0,   0.0,   0.0,   0.0,   0.015, 0.05 },   // l = 4
          {0.0,   0.0,   0.0,   0.0,   0.01,  0.04 },   // l = 5
          {0.0,   0.0,   0.0,   0.0,   0.01,  0.035},   // l = 6
          {0.0,   0.0,   0.0,   0.0,   0.01,  0.03 },   // l = 7
          {0.0,   0.0,   0.0,   0.0,   0.005, 0.03 },   // l = 8
          {0.0,   0.0,   0.0,   0.0,   0.005, 0.025},   // l = 9
          {0.0,   0.0,   0.0,   0.0,   0.005, 0.025}}}, // l = 10
    { 4, {{0.0,   0.005, 0.03,  0.05,  0.095, 0.19 },   // l = 2
          {0.0,   0.005, 0.02,  0.04,  0.08,  0.155},   // l = 3
          {0.0,   0.005, 0.02,  0.035, 0.065, 0.13 },   // l = 4
          {0.0,   0.005, 0.015, 0.03,  0.055, 0.115},   // l = 5
          {0.0,   0.005, 0.015, 0.025, 0.05,  0.1  },   // l = 6
          {0.0,   0.0,   0.01,  0.025, 0.045, 0.09 },   // l = 7
          {0.0,   0.0,   0.01,  0.02,  0.04,  0.08 },   // l = 8
          {0.0,   0.0,   0.01,  0.02,  0.035, 0.075},   // l = 9
          {0.0,   0.0,   0.01,  0.015, 0.035, 0.07 }}}, // l = 10
    { 6, {{0.015, 0.035, 0.075, 0.115, 0.175, 0.285},   // l = 2
          {0.01,  0.025, 0.06,  0.095, 0.145, 0.24 },   // l = 3
          {0.01,  0.025, 0.055, 0.08,  0.125, 0.21 },   // l = 4
          {0.01,  0.02,  0.045, 0.07,  0.11,  0.185},   // l = 5
          {0.005, 0.015, 0.04,  0.06,  0.095, 0.165},   // l = 6
          {0.005, 0.015, 0.035, 0.055, 0.085, 0.15 },   // l = 7
          {0.005, 0.015, 0.035, 0.05,  0.08,  0.135},   // l = 8
          {0.005, 0.01,  0.03,  0.045, 0.075, 0.125},   // l = 9
          {0.005, 0.01,  0.025, 0.04,  0.065, 0.115}}}, // l = 10
    { 8, {{0.035, 0.07,  0.125, 0.175, 0.245, 0.36 },   // l = 2
          {0.03,  0.06,  0.105, 0.15,  0.21,  0.31 },   // l = 3
          {0.025, 0.05,  0.095, 0.13,  0.18,  0.275},   // l = 4
          {0.025, 0.045, 0.08,  0.115, 0.16,  0.245},   // l = 5
          {0.02,  0.04,  0.075, 0.1,   0.145, 0.22 },   // l = 6
          {0.02,  0.035, 0.065, 0.09,  0.13,  0.2  },   // l = 7
          {0.015, 0.03,  0.06,  0.085, 0.12,  0.185},   // l = 8
          {0.015, 0.03,  0.055, 0.075, 0.11,  0.17 },   // l = 9
          {0.015, 0.025, 0.05,  0.07,  0.105, 0.16 }}}, // l = 10
    {10, {{0.065, 0.11,  0.18,  0.235, 0.31,  0.405},   // l = 2
          {0.055, 0.09,  0.155, 0.2,   0.265, 0.36 },   // l = 3
          {0.05,  0.08,  0.135, 0.175, 0.235, 0.325},   // l = 4
          {0.045, 0.07,  0.12,  0.155, 0.21,  0.3  },   // l = 5
          {0.04,  0.065, 0.105, 0.14,  0.19,  0.275},   // l = 6
          {0.035, 0.06,  0.095, 0.13,  0.175, 0.25 },   // l = 7
          {0.03,  0.055, 0.09,  0.12,  0.16,  0.23 },   // l = 8
          {0.03,  0.05,  0.08,  0.11,  0.15,  0.215},   // l = 9
          {0.025, 0.045, 0.075, 0.1,   0.14,  0.2  }}}, // l = 10
    {20, {{0.23,  0.295, 0.375, 0.425, 0.49,  0.595},   // l = 2
          {0.205, 0.265, 0.34,  0.38,  0.44,  0.535},   // l = 3
          {0.185, 0.24,  0.315, 0.35,  0.4,   0.485},   // l = 4
          {0.165, 0.22,  0.275, 0.31,  0.355, 0.415},   // l = 5
          {0.145, 0.185, 0.235, 0.265, 0.3,   0.34 },   // l = 6
          {0.135, 0.17,  0.215, 0.245, 0.275, 0.315},   // l = 7
          {0.13,  0.17,  0.215, 0.25,  0.285, 0.335},   // l = 8
          {0.125, 0.165, 0.215, 0.245, 0.285, 0.33 },   // l = 9
          {0.115, 0.155, 0.205, 0.24,  0.28,  0.32 }}}, // l = 10
    {50, {{0.515, 0.58,  0.65,  0.695, 0.745, 0.78 },   // l = 2
          {0.485, 0.545, 0.615, 0.66,  0.71,  0.75 },   // l = 3
          {0.455, 0.515, 0.585, 0.625, 0.675, 0.7  },   // l = 4
          {0.43,  0.485, 0.555, 0.6,   0.65,  0.655},   // l = 5
          {0.405, 0.46,  0.53,  0.57,  0.615, 0.615},   // l = 6
          {0.385, 0.44,  0.505, 0.545, 0.585, 0.585},   // l = 7
          {0.36,  0.405, 0.46,  0.495, 0.535, 0.555},   // l = 8
          {0.315, 0.35,  0.395, 0.42,  0.45,  0.485},   // l = 9
          {0.285, 0.315, 0.35,  0.37,  0.395, 0.42 }}}, // l = 10
};

// The largest eccentricity at which the level's heating stays within `tolerance` of the exact value at every degree
// from 2 to `max_degree_l`, using the tabulated tolerance at or below the requested one (so the answer never promises
// more than was measured); 0.0 for a tolerance below the smallest tabulated one, and NaN for an untabulated level (the
// exact option has no such limit).
inline double c_eccentricity_accuracy_limit(int eccentricity_truncation, double tolerance, int max_degree_l) noexcept {
    return c_truncation_accuracy_limit(C_ECCENTRICITY_ACCURACY, eccentricity_truncation, tolerance, max_degree_l);
}

// The eccentricity above which a level's heating can be 10% or more below the exact value; the solve warns past it.
inline double c_eccentricity_truncation_limit(int eccentricity_truncation, int max_degree_l) noexcept {
    return c_eccentricity_accuracy_limit(eccentricity_truncation, 1.0e-1, max_degree_l);
}

// The lowest tabulated level whose heating stays within `tolerance` at `eccentricity`, or C_ECCENTRICITY_EXACT when
// none does (e past about 0.75 at degree 3 and lower at higher degrees, or a tolerance too tight for the tables).
inline int c_recommend_eccentricity_truncation(double eccentricity, double tolerance, int max_degree_l) noexcept {
    return c_recommend_truncation(C_ECCENTRICITY_ACCURACY, eccentricity, tolerance, max_degree_l, C_ECCENTRICITY_EXACT);
}
