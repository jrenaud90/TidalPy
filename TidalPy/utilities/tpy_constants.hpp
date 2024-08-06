#include <limits>

// Numerical constants
static const double NAN_DBL = std::numeric_limits<double>::quiet_NaN();
static const double INF_DBL = std::numeric_limits<double>::infinity();
static const long double PI_LDBL = 3.141592653589793238462643383279502884L;
static const double PI_DBL = (double)PI_LDBL;

// Astrophysical Constants
static const double G = 6.67430e-11;

// Extremes
// Forcing Frequency Extremes
// Assume that any forcing period larger than a Gyr leads to a zero in frequency.
// Converting to frequency is roughly 1.0e-17 rads s-1
static const double MIN_FREQUENCY = 1.0e-17;
// Assume max frequency is for a forcing period of 1 micro-second
static const double MAX_FREQUENCY = 1.0e8;

/// Shear/Bulk Modulus Extremes
static const double MIN_MODULUS = 1.0e-3;