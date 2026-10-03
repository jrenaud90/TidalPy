#pragma once
/*
 * eccentricity_accuracy_.hpp - how far each eccentricity truncation level can be trusted, and which to choose.
 *
 * The limits were measured on 2026-10-02 (with a script kept outside the repo) from TidalPy's tables against the exact
 * functions at a tail tolerance of 1e-13: for each level, degree, and tolerance, the largest eccentricity of the grid
 * up to which that degree's heating stays within the tolerance. The grid is log-spaced below 0.005 (1, 1.25, 1.6, 2,
 * 2.5, 3.2, 4, 5, 6.3, and 8 times each power of ten from 1e-5) and steps by 0.005 above it, so 0.0 means the level
 * misses the tolerance already at e = 1e-5. Each limit is the worst case over constant-phase-lag and constant-time-lag
 * tides, Maxwell tides with relaxation times of 1e-3 to 1e3 inverse mean motions (ten per decade; one relaxation peak
 * at every time bounds any sum of them), and Andrade tides, at spin rates of 0 to 20 times the mean motion in steps of
 * 0.25 and at 25, 30, -0.5, -1, -2, -5, and -10 times it. The spin rate matters most at the high levels: a body
 * spinning several times faster than its mean motion with a strongly frequency-dependent rheology weighs the high-|q|
 * modes unevenly, and level 50 then holds 10% only to e = 0.57 at degree 2, where a synchronous body holds it to about
 * 0.8. Levels can err either way. Higher degrees generally lose accuracy at a lower eccentricity, so a solve through
 * max_degree_l takes the tightest limit of degrees 2 to max_degree_l. Within every limit the mode sum cancels by less
 * than 6e5, the cap that set level 50's earlier limits. See Documentation/Tides/Eccentricity.md.
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
    { 2, {{2e-5,    2e-4,    0.002,   0.005,   0.02,    0.075  },   // l = 2
          {1.6e-5,  1.6e-4,  0.0016,  0.005,   0.015,   0.06   },   // l = 3
          {1.25e-5, 1.25e-4, 0.00125, 0.004,   0.015,   0.05   },   // l = 4
          {1.25e-5, 1.25e-4, 0.00125, 0.004,   0.01,    0.04   },   // l = 5
          {1e-5,    1e-4,    0.001,   0.0032,  0.01,    0.035  },   // l = 6
          {1e-5,    1e-4,    0.001,   0.0032,  0.01,    0.03   },   // l = 7
          {0.0,     8e-5,    8e-4,    0.0025,  0.005,   0.03   },   // l = 8
          {0.0,     8e-5,    8e-4,    0.0025,  0.005,   0.025  },   // l = 9
          {0.0,     6.3e-5,  6.3e-4,  0.002,   0.005,   0.025  }}}, // l = 10
    { 4, {{0.0025,  0.005,   0.03,    0.05,    0.095,   0.18   },   // l = 2
          {0.002,   0.005,   0.02,    0.04,    0.08,    0.145  },   // l = 3
          {0.002,   0.005,   0.02,    0.035,   0.065,   0.125  },   // l = 4
          {0.0016,  0.005,   0.015,   0.03,    0.055,   0.105  },   // l = 5
          {0.00125, 0.005,   0.015,   0.025,   0.05,    0.095  },   // l = 6
          {0.00125, 0.004,   0.01,    0.025,   0.045,   0.085  },   // l = 7
          {0.00125, 0.004,   0.01,    0.02,    0.04,    0.075  },   // l = 8
          {0.001,   0.0032,  0.01,    0.02,    0.035,   0.07   },   // l = 9
          {0.001,   0.0032,  0.01,    0.015,   0.035,   0.065  }}}, // l = 10
    { 6, {{0.015,   0.035,   0.075,   0.11,    0.17,    0.265  },   // l = 2
          {0.01,    0.025,   0.06,    0.095,   0.14,    0.22   },   // l = 3
          {0.01,    0.025,   0.055,   0.08,    0.12,    0.19   },   // l = 4
          {0.01,    0.02,    0.045,   0.07,    0.105,   0.165  },   // l = 5
          {0.005,   0.015,   0.04,    0.06,    0.09,    0.145  },   // l = 6
          {0.005,   0.015,   0.035,   0.055,   0.085,   0.135  },   // l = 7
          {0.005,   0.015,   0.035,   0.05,    0.075,   0.12   },   // l = 8
          {0.005,   0.01,    0.03,    0.045,   0.07,    0.11   },   // l = 9
          {0.005,   0.01,    0.025,   0.04,    0.065,   0.1    }}}, // l = 10
    { 8, {{0.035,   0.07,    0.125,   0.17,    0.235,   0.335  },   // l = 2
          {0.03,    0.06,    0.105,   0.145,   0.195,   0.28   },   // l = 3
          {0.025,   0.05,    0.09,    0.125,   0.17,    0.245  },   // l = 4
          {0.025,   0.045,   0.08,    0.11,    0.15,    0.22   },   // l = 5
          {0.02,    0.04,    0.07,    0.095,   0.135,   0.195  },   // l = 6
          {0.02,    0.035,   0.065,   0.085,   0.125,   0.175  },   // l = 7
          {0.015,   0.03,    0.06,    0.08,    0.11,    0.16   },   // l = 8
          {0.015,   0.03,    0.055,   0.075,   0.1,     0.145  },   // l = 9
          {0.015,   0.025,   0.05,    0.065,   0.095,   0.135  }}}, // l = 10
    {10, {{0.065,   0.11,    0.175,   0.225,   0.29,    0.395  },   // l = 2
          {0.055,   0.09,    0.15,    0.19,    0.25,    0.335  },   // l = 3
          {0.05,    0.08,    0.13,    0.165,   0.22,    0.29   },   // l = 4
          {0.045,   0.07,    0.115,   0.15,    0.195,   0.255  },   // l = 5
          {0.04,    0.065,   0.1,     0.135,   0.175,   0.23   },   // l = 6
          {0.035,   0.06,    0.095,   0.12,    0.155,   0.215  },   // l = 7
          {0.03,    0.055,   0.085,   0.11,    0.14,    0.175  },   // l = 8
          {0.03,    0.05,    0.08,    0.1,     0.135,   0.19   },   // l = 9
          {0.025,   0.045,   0.07,    0.095,   0.12,    0.16   }}}, // l = 10
    {20, {{0.215,   0.27,    0.345,   0.385,   0.44,    0.5    },   // l = 2
          {0.185,   0.235,   0.295,   0.335,   0.38,    0.435  },   // l = 3
          {0.165,   0.205,   0.26,    0.295,   0.335,   0.38   },   // l = 4
          {0.15,    0.19,    0.24,    0.275,   0.31,    0.36   },   // l = 5
          {0.13,    0.165,   0.21,    0.235,   0.27,    0.31   },   // l = 6
          {0.12,    0.15,    0.195,   0.22,    0.25,    0.285  },   // l = 7
          {0.11,    0.135,   0.175,   0.195,   0.225,   0.255  },   // l = 8
          {0.1,     0.125,   0.165,   0.185,   0.21,    0.245  },   // l = 9
          {0.095,   0.125,   0.16,    0.18,    0.21,    0.24   }}}, // l = 10
    {50, {{0.405,   0.445,   0.49,    0.515,   0.545,   0.57   },   // l = 2
          {0.37,    0.41,    0.45,    0.475,   0.5,     0.525  },   // l = 3
          {0.34,    0.375,   0.415,   0.44,    0.46,    0.485  },   // l = 4
          {0.32,    0.35,    0.39,    0.41,    0.43,    0.455  },   // l = 5
          {0.295,   0.33,    0.365,   0.385,   0.405,   0.425  },   // l = 6
          {0.28,    0.31,    0.34,    0.36,    0.38,    0.4    },   // l = 7
          {0.265,   0.29,    0.325,   0.34,    0.36,    0.38   },   // l = 8
          {0.25,    0.28,    0.31,    0.325,   0.345,   0.365  },   // l = 9
          {0.235,   0.26,    0.29,    0.305,   0.325,   0.345  }}}, // l = 10
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
// none does (e past 0.57 at degree 2 and lower at higher degrees, or a tolerance too tight for the tables).
inline int c_recommend_eccentricity_truncation(double eccentricity, double tolerance, int max_degree_l) noexcept {
    return c_recommend_truncation(C_ECCENTRICITY_ACCURACY, eccentricity, tolerance, max_degree_l, C_ECCENTRICITY_EXACT);
}
