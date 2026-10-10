#pragma once
/*
 * obliquity_accuracy_.hpp - how far each obliquity truncation level can be trusted, and which to choose.
 *
 * The limits were measured on 2026-10-02 (with a script kept outside the repo) from TidalPy's tables against the
 * general (exact) obliquity functions: for each level, degree, tolerance, and spin band, the largest obliquity [rad]
 * of the grid up to which that degree's heating stays within the tolerance. The grid is log-spaced below 0.005 rad
 * (1, 1.25, 1.6, 2, 2.5, 3.2, 4, 5, 6.3, and 8 times each power of ten from 1e-5) and steps by 0.005 rad above it.
 * Each limit is the worst case over constant-phase-lag and constant-time-lag tides, Maxwell tides with relaxation
 * times of 1e-3 to 1e3 inverse mean motions (ten per decade; one relaxation peak at every time bounds any sum of
 * them), and Andrade tides, at the spin rates of its band (truncation_accuracy_.hpp) and eccentricities of 0, 0.001,
 * 0.01, 0.05, and 0.2 (eccentricity level 20). The limits of the bands that hold synchronous rotation are set by it,
 * at low eccentricity with a long relaxation time (-Im k falling as 1 / frequency, as for a stiff Maxwell body), where
 * the obliquity tide carries the heating; level 2 errs high there and level 4 low. Each higher degree loses accuracy
 * at a lower obliquity, so a solve through max_degree_l takes the tightest limit of degrees 2 to max_degree_l. Level
 * 0 (off) holds only at zero obliquity. See Documentation/Tides/Obliquity.md.
 *
 * Standard-library only, so any extension can include it.
 */

#include <cmath>
#include <limits>

#include "../truncation_accuracy_.hpp"

// Obliquity off: the functions at I = 0 (level 0 of the tables).
inline constexpr int C_OBLIQUITY_OFF = 0;

// The truncation code of the general obliquity functions: the exact half-angle form, any obliquity.
inline constexpr int C_OBLIQUITY_GENERAL = -1;

// One table per spin band (truncation_accuracy_.hpp), then the any-spin table; in each, one row per degree (l = 2 to
// 10) and one column per tolerance (1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 1e-1).
inline constexpr c_TruncationAccuracyRow C_OBLIQUITY_ACCURACY[C_TRUNCATION_NUM_SPIN_TABLES][3] = {
    {   // |spin / n| <= 1.5
        {0, {}},   // every degree: zero obliquity only
        {2, {{1e-4,   0.001,  0.01,   0.035,  0.11,   0.35  },   // l = 2
             {6.3e-5, 6.3e-4, 0.005,  0.02,   0.075,  0.24  },   // l = 3
             {5e-5,   5e-4,   0.005,  0.015,  0.06,   0.185 },   // l = 4
             {4e-5,   4e-4,   0.004,  0.015,  0.045,  0.15  },   // l = 5
             {4e-5,   4e-4,   0.004,  0.01,   0.04,   0.125 },   // l = 6
             {3.2e-5, 3.2e-4, 0.0032, 0.01,   0.035,  0.11  },   // l = 7
             {2.5e-5, 2.5e-4, 0.0025, 0.005,  0.03,   0.095 },   // l = 8
             {2.5e-5, 2.5e-4, 0.0025, 0.005,  0.025,  0.085 },   // l = 9
             {2.5e-5, 2.5e-4, 0.0025, 0.005,  0.025,  0.075 }}}, // l = 10
        {4, {{0.01,   0.04,   0.13,   0.23,   0.405,  0.695 },   // l = 2
             {0.005,  0.025,  0.09,   0.16,   0.28,   0.48  },   // l = 3
             {0.005,  0.02,   0.065,  0.12,   0.215,  0.37  },   // l = 4
             {0.005,  0.015,  0.055,  0.1,    0.175,  0.3   },   // l = 5
             {0.004,  0.015,  0.045,  0.08,   0.145,  0.25  },   // l = 6
             {0.004,  0.01,   0.04,   0.07,   0.125,  0.215 },   // l = 7
             {0.0032, 0.01,   0.035,  0.06,   0.11,   0.19  },   // l = 8
             {0.0032, 0.01,   0.03,   0.055,  0.1,    0.17  },   // l = 9
             {0.0025, 0.005,  0.025,  0.05,   0.09,   0.155 }}}, // l = 10
    },
    {   // 1.5 <= |spin / n| <= 5
        {0, {}},   // every degree: zero obliquity only
        {2, {{0.005,  0.025,  0.09,   0.165,  0.3,    0.58  },   // l = 2
             {0.005,  0.02,   0.07,   0.13,   0.245,  0.48  },   // l = 3
             {0.005,  0.015,  0.055,  0.1,    0.185,  0.35  },   // l = 4
             {0.005,  0.015,  0.05,   0.09,   0.17,   0.33  },   // l = 5
             {0.004,  0.01,   0.04,   0.07,   0.135,  0.255 },   // l = 6
             {0.004,  0.01,   0.04,   0.07,   0.135,  0.255 },   // l = 7
             {0.0032, 0.01,   0.03,   0.055,  0.105,  0.205 },   // l = 8
             {0.0032, 0.01,   0.035,  0.06,   0.11,   0.21  },   // l = 9
             {0.0025, 0.005,  0.025,  0.045,  0.09,   0.17  }}}, // l = 10
        {4, {{0.045,  0.105,  0.23,   0.345,  0.525,  0.795 },   // l = 2
             {0.035,  0.08,   0.18,   0.27,   0.4,    0.595 },   // l = 3
             {0.025,  0.06,   0.135,  0.2,    0.305,  0.465 },   // l = 4
             {0.025,  0.055,  0.125,  0.185,  0.27,   0.41  },   // l = 5
             {0.02,   0.045,  0.095,  0.145,  0.22,   0.335 },   // l = 6
             {0.02,   0.04,   0.09,   0.135,  0.2,    0.305 },   // l = 7
             {0.015,  0.035,  0.075,  0.115,  0.17,   0.26  },   // l = 8
             {0.015,  0.03,   0.075,  0.11,   0.16,   0.245 },   // l = 9
             {0.01,   0.025,  0.06,   0.095,  0.14,   0.215 }}}, // l = 10
    },
    {   // 5 <= |spin / n| <= 30
        {0, {}},   // every degree: zero obliquity only
        {2, {{0.005,  0.025,  0.09,   0.165,  0.3,    0.565 },   // l = 2
             {0.005,  0.02,   0.075,  0.135,  0.24,   0.44  },   // l = 3
             {0.005,  0.015,  0.06,   0.105,  0.195,  0.37  },   // l = 4
             {0.005,  0.015,  0.05,   0.095,  0.17,   0.32  },   // l = 5
             {0.004,  0.01,   0.045,  0.08,   0.15,   0.29  },   // l = 6
             {0.004,  0.01,   0.04,   0.075,  0.135,  0.26  },   // l = 7
             {0.0032, 0.01,   0.035,  0.065,  0.12,   0.24  },   // l = 8
             {0.0032, 0.01,   0.035,  0.06,   0.115,  0.22  },   // l = 9
             {0.0032, 0.01,   0.03,   0.055,  0.105,  0.205 }}}, // l = 10
        {4, {{0.045,  0.105,  0.225,  0.335,  0.505,  0.76  },   // l = 2
             {0.035,  0.08,   0.17,   0.255,  0.37,   0.555 },   // l = 3
             {0.025,  0.06,   0.135,  0.2,    0.3,    0.455 },   // l = 4
             {0.025,  0.05,   0.115,  0.17,   0.255,  0.38  },   // l = 5
             {0.02,   0.045,  0.1,    0.145,  0.22,   0.335 },   // l = 6
             {0.015,  0.04,   0.09,   0.13,   0.195,  0.295 },   // l = 7
             {0.015,  0.035,  0.08,   0.115,  0.175,  0.265 },   // l = 8
             {0.015,  0.03,   0.07,   0.105,  0.16,   0.245 },   // l = 9
             {0.01,   0.03,   0.065,  0.1,    0.145,  0.225 }}}, // l = 10
    },
    {   // any spin rate
        {0, {}},   // every degree: zero obliquity only
        {2, {{1e-4,   0.001,  0.01,   0.035,  0.11,   0.35  },   // l = 2
             {6.3e-5, 6.3e-4, 0.005,  0.02,   0.075,  0.24  },   // l = 3
             {5e-5,   5e-4,   0.005,  0.015,  0.06,   0.185 },   // l = 4
             {4e-5,   4e-4,   0.004,  0.015,  0.045,  0.15  },   // l = 5
             {4e-5,   4e-4,   0.004,  0.01,   0.04,   0.125 },   // l = 6
             {3.2e-5, 3.2e-4, 0.0032, 0.01,   0.035,  0.11  },   // l = 7
             {2.5e-5, 2.5e-4, 0.0025, 0.005,  0.03,   0.095 },   // l = 8
             {2.5e-5, 2.5e-4, 0.0025, 0.005,  0.025,  0.085 },   // l = 9
             {2.5e-5, 2.5e-4, 0.0025, 0.005,  0.025,  0.075 }}}, // l = 10
        {4, {{0.01,   0.04,   0.13,   0.23,   0.405,  0.695 },   // l = 2
             {0.005,  0.025,  0.09,   0.16,   0.28,   0.48  },   // l = 3
             {0.005,  0.02,   0.065,  0.12,   0.215,  0.37  },   // l = 4
             {0.005,  0.015,  0.055,  0.1,    0.175,  0.3   },   // l = 5
             {0.004,  0.015,  0.045,  0.08,   0.145,  0.25  },   // l = 6
             {0.004,  0.01,   0.04,   0.07,   0.125,  0.215 },   // l = 7
             {0.0032, 0.01,   0.035,  0.06,   0.11,   0.19  },   // l = 8
             {0.0032, 0.01,   0.03,   0.055,  0.1,    0.17  },   // l = 9
             {0.0025, 0.005,  0.025,  0.05,   0.09,   0.155 }}}, // l = 10
    },
};

// The largest obliquity [rad] at which the level's heating stays within `tolerance` of the general value at every
// degree from 2 to `max_degree_l` and every spin rate of the band of `spin_ratio` (spin rate / mean motion; NaN, the
// default, for any spin rate), using the tabulated tolerance at or below the requested one (so the answer never
// promises more than was measured); 0.0 for a tolerance below the smallest tabulated one, and NaN for an untabulated
// level (the general functions have no such limit).
inline double c_obliquity_accuracy_limit(
        int obliquity_truncation,
        double tolerance,
        int max_degree_l,
        double spin_ratio = std::numeric_limits<double>::quiet_NaN()) noexcept {
    return c_truncation_accuracy_limit(
        C_OBLIQUITY_ACCURACY[c_truncation_spin_table(spin_ratio)], obliquity_truncation, tolerance, max_degree_l);
}

// The obliquity [rad] above which a level's heating can be 10% or more off the general value; the solve warns past it.
inline double c_obliquity_truncation_limit(
        int obliquity_truncation,
        int max_degree_l,
        double spin_ratio = std::numeric_limits<double>::quiet_NaN()) noexcept {
    return c_obliquity_accuracy_limit(obliquity_truncation, 1.0e-1, max_degree_l, spin_ratio);
}

// The lowest level whose heating stays within `tolerance` at `obliquity` [rad] and the spin ratio's band (off only at
// zero obliquity), or C_OBLIQUITY_GENERAL when none does.
inline int c_recommend_obliquity_truncation(
        double obliquity,
        double tolerance,
        int max_degree_l,
        double spin_ratio = std::numeric_limits<double>::quiet_NaN()) noexcept {
    return c_recommend_truncation(
        C_OBLIQUITY_ACCURACY[c_truncation_spin_table(spin_ratio)], std::abs(obliquity), tolerance, max_degree_l,
        C_OBLIQUITY_GENERAL);
}
