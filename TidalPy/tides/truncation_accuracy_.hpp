#pragma once
/*
 * truncation_accuracy_.hpp - the measured accuracy table format shared by the eccentricity and obliquity truncation
 * levels (eccentricity/eccentricity_accuracy_.hpp, obliquity/obliquity_accuracy_.hpp), and its lookups.
 *
 * Standard-library only, so any extension can include it.
 */

#include <limits>

inline constexpr int C_TRUNCATION_ACCURACY_NUM_TOLERANCES = 6;
inline constexpr double C_TRUNCATION_ACCURACY_TOLERANCES[C_TRUNCATION_ACCURACY_NUM_TOLERANCES] = {
    1.0e-8, 1.0e-6, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1};

// One truncation level's limits, one per tabulated tolerance.
struct c_TruncationAccuracyRow {
    int truncation;
    double degree_two[C_TRUNCATION_ACCURACY_NUM_TOLERANCES];     // limits for max_degree_l == 2
    double degree_three[C_TRUNCATION_ACCURACY_NUM_TOLERANCES];   // limits for max_degree_l >= 3
};

// The limit of a level at `tolerance`, using the tabulated tolerance at or below the requested one (so the answer
// never promises more than was measured); 0.0 for a tolerance below the smallest tabulated one, and NaN for a level
// the table does not hold.
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
    for (int level = 0; level < NumLevels; ++level) {
        const c_TruncationAccuracyRow& row = rows[level];
        if (row.truncation != truncation) { continue; }
        if (column < 0) { return 0.0; }
        return (max_degree_l <= 2) ? row.degree_two[column] : row.degree_three[column];
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
