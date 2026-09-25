#pragma once
/*
 * obliquity_accuracy_.hpp - how far each obliquity truncation level can be trusted, and which to choose.
 *
 * The limits are measured (2026-09-24) with TidalPy's tables against the general (exact) obliquity functions: for each
 * level, degree, and tolerance, the largest obliquity [rad] (on a 0.005 rad grid) at which the heating stays within the
 * tolerance, worst case over constant-phase-lag, constant-time-lag, and Maxwell tides (relaxation times of 1/3 and 1/10
 * of the orbital period) at spin rates of 0.5, 1, and 2.3 times the mean motion and eccentricities of 0 and 0.2 (the
 * limits are the same at both). Degree 3 loses accuracy at a lower obliquity than degree 2, so a solve that includes
 * degree 3 or higher uses the degree-3 limits. Level 0 (off) holds only at zero obliquity. See
 * Documentation/Tides_x/obliquity.md.
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
    {0, {0.0,   0.0,   0.0,   0.0,   0.0,   0.0  }, {0.0,   0.0,   0.0,   0.0,   0.0,   0.0  }},
    {2, {0.0,   0.0,   0.01,  0.045, 0.145, 0.465}, {0.0,   0.0,   0.005, 0.03,  0.095, 0.31 }},
    {4, {0.015, 0.045, 0.15,  0.265, 0.47,  0.83 }, {0.01,  0.03,  0.1,   0.18,  0.32,  0.565}},
};

// The largest obliquity [rad] at which the level's heating stays within `tolerance` of the general value, using the
// tabulated tolerance at or below the requested one (so the answer never promises more than was measured); 0.0 for
// a tolerance below the smallest tabulated one, and NaN for an untabulated level (the general functions have no such
// limit).
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
