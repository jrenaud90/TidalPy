#pragma once
/*
 * truncation_accuracy_.hpp - the measured accuracy table format shared by the eccentricity and obliquity truncation
 * levels (eccentricity/eccentricity_accuracy_.hpp, obliquity/obliquity_accuracy_.hpp), and its lookups.
 *
 * Standard-library only, so any extension can include it.
 */

#include <algorithm>
#include <cmath>
#include <limits>

// The spin bands the limits are measured in: |spin rate / mean motion| in [0, 1.5] (synchronous rotation and the 3:2
// resonance), [1.5, 5], and [5, 30], each including its edges. A table per band, then a table over every measured spin
// rate, which a ratio past the last edge (not measured) or not finite takes. A band holds both signs of the ratio;
// retrograde spins were measured at -0.5, -1, -2, -5, and -10 times the mean motion.
inline constexpr int C_TRUNCATION_NUM_SPIN_BANDS = 3;
inline constexpr double C_TRUNCATION_SPIN_BAND_EDGES[C_TRUNCATION_NUM_SPIN_BANDS] = {1.5, 5.0, 30.0};
inline constexpr int C_TRUNCATION_NUM_SPIN_TABLES = C_TRUNCATION_NUM_SPIN_BANDS + 1;

// The table of a spin ratio (spin rate / mean motion, either sign): its band, or the any-spin table
// (C_TRUNCATION_NUM_SPIN_BANDS) for a ratio past the last edge or NaN.
inline int c_truncation_spin_table(double spin_ratio) noexcept {
    const double magnitude = std::abs(spin_ratio);
    for (int band = 0; band < C_TRUNCATION_NUM_SPIN_BANDS; ++band) {
        if (magnitude <= C_TRUNCATION_SPIN_BAND_EDGES[band]) { return band; }
    }
    return C_TRUNCATION_NUM_SPIN_BANDS;
}

inline constexpr int C_TRUNCATION_ACCURACY_NUM_TOLERANCES = 6;
inline constexpr double C_TRUNCATION_ACCURACY_TOLERANCES[C_TRUNCATION_ACCURACY_NUM_TOLERANCES] = {
    1.0e-8, 1.0e-6, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1};

// The tabulated degrees: l = 2 to 10, each measured on its own.
inline constexpr int C_TRUNCATION_ACCURACY_MIN_DEGREE = 2;
inline constexpr int C_TRUNCATION_ACCURACY_MAX_DEGREE = 10;
inline constexpr int C_TRUNCATION_ACCURACY_NUM_DEGREES =
    C_TRUNCATION_ACCURACY_MAX_DEGREE - C_TRUNCATION_ACCURACY_MIN_DEGREE + 1;

// One truncation level's limits: one row per degree (l = 2 first), one value per tabulated tolerance.
struct c_TruncationAccuracyRow {
    int truncation;
    double by_degree[C_TRUNCATION_ACCURACY_NUM_DEGREES][C_TRUNCATION_ACCURACY_NUM_TOLERANCES];
};

// The limit of a level at `tolerance` for a solve that includes degrees 2 to `max_degree_l`: the tightest of those
// degrees' limits (a degree past the table uses its last degree's). It uses the tabulated tolerance at or below the
// requested one (so the answer never promises more than was measured); 0.0 for a tolerance below the smallest
// tabulated one, and NaN for a level the table does not hold.
template <int NumLevels>
inline double c_truncation_accuracy_limit(
        const c_TruncationAccuracyRow (&rows)[NumLevels],
        int truncation,
        double tolerance,
        int max_degree_l) noexcept {
    int column = -1;
    for (int i = 0; i < C_TRUNCATION_ACCURACY_NUM_TOLERANCES; ++i) {
        if (C_TRUNCATION_ACCURACY_TOLERANCES[i] <= tolerance) { column = i; }
    }
    const int last_degree =
        std::clamp(max_degree_l, C_TRUNCATION_ACCURACY_MIN_DEGREE, C_TRUNCATION_ACCURACY_MAX_DEGREE);
    for (int level = 0; level < NumLevels; ++level) {
        const c_TruncationAccuracyRow& row = rows[level];
        if (row.truncation != truncation) { continue; }
        if (column < 0) { return 0.0; }
        double limit = row.by_degree[0][column];
        for (int degree_l = C_TRUNCATION_ACCURACY_MIN_DEGREE + 1; degree_l <= last_degree; ++degree_l) {
            limit = std::min(limit, row.by_degree[degree_l - C_TRUNCATION_ACCURACY_MIN_DEGREE][column]);
        }
        return limit;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

// The lowest level whose limit at `tolerance` reaches `magnitude`, or `fallback` when none does.
template <int NumLevels>
inline int c_recommend_truncation(
        const c_TruncationAccuracyRow (&rows)[NumLevels],
        double magnitude,
        double tolerance,
        int max_degree_l,
        int fallback) noexcept {
    for (int level = 0; level < NumLevels; ++level) {
        const int truncation = rows[level].truncation;
        if (magnitude <= c_truncation_accuracy_limit(rows, truncation, tolerance, max_degree_l)) {
            return truncation;
        }
    }
    return fallback;
}
