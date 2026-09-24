#pragma once
/*
 * eccentricity_accuracy_.hpp - how far each eccentricity truncation level can be trusted, and which to choose.
 *
 * The limits are measured (2026-09-24) with TidalPy's tables against the exact heating, from exact Hansen coefficients
 * summed over |q| <= 4000: for each level, degree, and tolerance, the largest eccentricity (on a 0.005 grid) at which the
 * heating stays within the tolerance, worst case over constant-phase-lag, constant-time-lag, and Maxwell tides
 * (relaxation times of 1/3 and 1/10 of the orbital period) at spin rates of 0.5, 1, and 2.3 times the mean motion.
 * Every level errs low. Degree 3 loses accuracy at a lower eccentricity than degree 2, so a solve that includes degree 3
 * or higher uses the degree-3 limits. Level 50 is capped at 0.78 (0.75 at degree 3): past that its mode sum cancels by
 * more than about 5e5, which multiplies any error in the Love numbers. See Documentation/Tides_x/eccentricity.md.
 *
 * Standard-library only, so any extension can include it.
 */

#include <cmath>
#include <limits>

// The truncation code of the exact eccentricity functions: Hansen coefficients from the exact Kepler orbit, the mode
// range chosen by a heating tolerance (eccentricity_exact_.hpp).
inline constexpr int C_ECCENTRICITY_EXACT = -1;

// Default tail tolerance of the exact eccentricity functions; the configuration's [tides]
// eccentricity_exact_tolerance supplies the value a built world uses.
inline constexpr double C_ECCENTRICITY_EXACT_TOLERANCE = 1.0e-4;

inline constexpr int C_ECCENTRICITY_ACCURACY_NUM_TOLERANCES = 6;
inline constexpr double C_ECCENTRICITY_ACCURACY_TOLERANCES[C_ECCENTRICITY_ACCURACY_NUM_TOLERANCES] = {
    1.0e-8, 1.0e-6, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1};

struct c_EccentricityAccuracyRow {
    int truncation;
    double degree_two[C_ECCENTRICITY_ACCURACY_NUM_TOLERANCES];     // limits for max_degree_l == 2
    double degree_three[C_ECCENTRICITY_ACCURACY_NUM_TOLERANCES];   // limits for max_degree_l >= 3
};

inline constexpr int C_ECCENTRICITY_ACCURACY_NUM_LEVELS = 7;
inline constexpr c_EccentricityAccuracyRow C_ECCENTRICITY_ACCURACY[C_ECCENTRICITY_ACCURACY_NUM_LEVELS] = {
    { 2, {0.0,   0.0,   0.0,   0.005, 0.02,  0.075}, {0.0,   0.0,   0.0,   0.005, 0.015, 0.06 }},
    { 4, {0.0,   0.005, 0.03,  0.05,  0.095, 0.19 }, {0.0,   0.005, 0.02,  0.04,  0.08,  0.155}},
    { 6, {0.015, 0.035, 0.075, 0.115, 0.175, 0.285}, {0.01,  0.025, 0.06,  0.095, 0.145, 0.24 }},
    { 8, {0.035, 0.07,  0.125, 0.175, 0.245, 0.36 }, {0.03,  0.06,  0.105, 0.15,  0.21,  0.31 }},
    {10, {0.065, 0.11,  0.18,  0.235, 0.31,  0.405}, {0.055, 0.09,  0.155, 0.2,   0.265, 0.36 }},
    {20, {0.23,  0.295, 0.375, 0.425, 0.49,  0.595}, {0.205, 0.265, 0.34,  0.38,  0.44,  0.535}},
    {50, {0.515, 0.58,  0.65,  0.695, 0.745, 0.78 }, {0.485, 0.545, 0.615, 0.66,  0.71,  0.75 }},
};

// The largest eccentricity at which the level's heating stays within `tolerance` of the exact value, using the
// tabulated tolerance at or below the requested one (so the answer never promises more than was measured); 0.0 for
// a tolerance below the smallest tabulated one, and NaN for an untabulated level (the exact option has no such limit).
inline double c_eccentricity_accuracy_limit(int eccentricity_truncation, double tolerance, int max_degree_l) noexcept {
    int column = -1;
    for (int i = 0; i < C_ECCENTRICITY_ACCURACY_NUM_TOLERANCES; ++i) {
        if (C_ECCENTRICITY_ACCURACY_TOLERANCES[i] <= tolerance) { column = i; }
    }
    for (int level = 0; level < C_ECCENTRICITY_ACCURACY_NUM_LEVELS; ++level) {
        const c_EccentricityAccuracyRow& row = C_ECCENTRICITY_ACCURACY[level];
        if (row.truncation != eccentricity_truncation) { continue; }
        if (column < 0) { return 0.0; }
        return (max_degree_l <= 2) ? row.degree_two[column] : row.degree_three[column];
    }
    return std::numeric_limits<double>::quiet_NaN();
}

// The eccentricity above which a level's heating can be 10% or more below the exact value; the solve warns past it.
inline double c_eccentricity_truncation_limit(int eccentricity_truncation, int max_degree_l) noexcept {
    return c_eccentricity_accuracy_limit(eccentricity_truncation, 1.0e-1, max_degree_l);
}

// The lowest tabulated level whose heating stays within `tolerance` at `eccentricity`, or C_ECCENTRICITY_EXACT when
// none does (e past about 0.75, or a tolerance too tight for the tables).
inline int c_recommend_eccentricity_truncation(double eccentricity, double tolerance, int max_degree_l) noexcept {
    for (int level = 0; level < C_ECCENTRICITY_ACCURACY_NUM_LEVELS; ++level) {
        const int truncation = C_ECCENTRICITY_ACCURACY[level].truncation;
        if (eccentricity <= c_eccentricity_accuracy_limit(truncation, tolerance, max_degree_l)) {
            return truncation;
        }
    }
    return C_ECCENTRICITY_EXACT;
}
