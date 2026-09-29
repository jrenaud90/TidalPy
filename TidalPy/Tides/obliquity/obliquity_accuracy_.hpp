#pragma once
/*
 * obliquity_accuracy_.hpp - how far each obliquity truncation level can be trusted, and which to choose.
 *
 * The limits are measured (2026-09-24 for degrees 2 and 3, 2026-09-29 for degrees 4 to 10) with TidalPy's tables
 * against the general (exact) obliquity functions: for each level, degree, and tolerance, the largest obliquity [rad]
 * (on a 0.005 rad grid) at which that degree's heating stays within the tolerance, worst case over constant-phase-lag,
 * constant-time-lag, and Maxwell tides (relaxation times of 1/3 and 1/10 of the inverse mean motion) at spin rates of
 * 0.5, 1, and 2.3 times the mean motion and eccentricities of 0 and 0.2 (the limits are the same at both). Each higher
 * degree loses accuracy at a lower obliquity, so a solve through max_degree_l takes the tightest limit of degrees 2 to
 * max_degree_l. Level 0 (off) holds only at zero obliquity. See Documentation/Tides/Obliquity.md.
 *
 * Standard-library only, so any extension can include it.
 */

#include <cmath>

#include "../truncation_accuracy_.hpp"

// Obliquity off: the functions at I = 0 (level 0 of the tables).
inline constexpr int C_OBLIQUITY_OFF = 0;

// The truncation code of the general obliquity functions: the exact half-angle form, any obliquity.
inline constexpr int C_OBLIQUITY_GENERAL = -1;

// Limits [rad].
inline constexpr c_TruncationAccuracyRow C_OBLIQUITY_ACCURACY[3] = {
    // One row per degree (l = 2 to 10), one column per tolerance (1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 1e-1).
    {0, {}},                                           // every degree: zero obliquity only
    {2, {{0.0,   0.0,   0.01,  0.045, 0.145, 0.465},   // l = 2
         {0.0,   0.0,   0.005, 0.03,  0.095, 0.31 },   // l = 3
         {0.0,   0.0,   0.005, 0.02,  0.075, 0.235},   // l = 4
         {0.0,   0.0,   0.005, 0.015, 0.06,  0.19 },   // l = 5
         {0.0,   0.0,   0.005, 0.015, 0.05,  0.155},   // l = 6
         {0.0,   0.0,   0.0,   0.01,  0.04,  0.135},   // l = 7
         {0.0,   0.0,   0.0,   0.01,  0.035, 0.12 },   // l = 8
         {0.0,   0.0,   0.0,   0.01,  0.03,  0.105},   // l = 9
         {0.0,   0.0,   0.0,   0.005, 0.03,  0.095}}}, // l = 10
    {4, {{0.015, 0.045, 0.15,  0.265, 0.47,  0.83 },   // l = 2
         {0.01,  0.03,  0.1,   0.18,  0.32,  0.565},   // l = 3
         {0.005, 0.025, 0.075, 0.14,  0.245, 0.43 },   // l = 4
         {0.005, 0.02,  0.06,  0.11,  0.2,   0.35 },   // l = 5
         {0.005, 0.015, 0.05,  0.095, 0.17,  0.295},   // l = 6
         {0.0,   0.01,  0.045, 0.08,  0.145, 0.255},   // l = 7
         {0.0,   0.01,  0.04,  0.07,  0.125, 0.225},   // l = 8
         {0.0,   0.01,  0.035, 0.065, 0.115, 0.2  },   // l = 9
         {0.0,   0.01,  0.03,  0.055, 0.1,   0.18 }}}, // l = 10
};

// The largest obliquity [rad] at which the level's heating stays within `tolerance` of the general value at every
// degree from 2 to `max_degree_l`, using the tabulated tolerance at or below the requested one (so the answer never
// promises more than was measured); 0.0 for a tolerance below the smallest tabulated one, and NaN for an untabulated
// level (the general functions have no such limit).
inline double c_obliquity_accuracy_limit(int obliquity_truncation, double tolerance, int max_degree_l) noexcept {
    return c_truncation_accuracy_limit(C_OBLIQUITY_ACCURACY, obliquity_truncation, tolerance, max_degree_l);
}

// The obliquity [rad] above which a level's heating can be 10% or more off the general value; the solve warns past it.
inline double c_obliquity_truncation_limit(int obliquity_truncation, int max_degree_l) noexcept {
    return c_obliquity_accuracy_limit(obliquity_truncation, 1.0e-1, max_degree_l);
}

// The lowest level whose heating stays within `tolerance` at `obliquity` [rad] (off only at zero obliquity), or
// C_OBLIQUITY_GENERAL when none does.
inline int c_recommend_obliquity_truncation(double obliquity, double tolerance, int max_degree_l) noexcept {
    return c_recommend_truncation(
        C_OBLIQUITY_ACCURACY, std::abs(obliquity), tolerance, max_degree_l, C_OBLIQUITY_GENERAL);
}
