#include <limits>

// Numerical constants
static const double NAN_DBL = std::numeric_limits<double>::quiet_NaN();
static const double INF_DBL = std::numeric_limits<double>::infinity();
static const long double PI_LDBL = 3.141592653589793238462643383279502884L;
static const double PI_DBL = (double)PI_LDBL;

// Astrophysical Constants
static const double au = 149597870700.0; // Au m-1
static const double G = 6.67430e-11;     // m^3 kg^-1 s^-2

// Sun
static const double mass_solar       = 1.988435e30;   // [kg]
static const double radius_solar     = 6.957e8;       // [m]
static const double luminosity_solar = 3.848e26;      // [Watts]

// Earth
static const double mass_earth   = 5.9721986e24;  // [kg]
static const double radius_earth = 6.371008e6;    // [m]

// Jupiter
static const double mass_jupiter   = 1.89813e27;  // [kg]
static const double radius_jupiter = 6.9911e7;    // [m]

// Pluto
static const double mass_pluto   = 1.309e22;  // [kg]
static const double radius_pluto = 1.1899e6;  // [m]

// Io
static const double mass_io   = 8.9298e22;  // [kg]
static const double radius_io = 1.82149e6;  // [m]

// Other Physics
static const double molar_gas_constant        = 8.3144598;           // J mol^-1 K^-1
static const double stefan_boltzmann_constant = 5.6703744191844e-8;  // W m^-2 K^-4

// Mathematics
static const double ppm = 1.e-6;
static const double ppb = 1.e-9;

// SI prefixes
static const double quetta = 1e30;
static const double ronna = 1e27;
static const double yotta = 1e24;
static const double zetta = 1e21;
static const double exa = 1e18;
static const double peta = 1e15;
static const double tera = 1e12;
static const double giga = 1e9;
static const double mega = 1e6;
static const double kilo = 1e3;
static const double hecto = 1e2;
static const double deka = 1e1;
static const double deci = 1e-1;
static const double centi = 1e-2;
static const double milli = 1e-3;
static const double micro = 1e-6;
static const double nano = 1e-9;
static const double pico = 1e-12;
static const double femto = 1e-15;
static const double atto = 1e-18;
static const double zepto = 1e-21;
static const double yocto = 1e-24;
static const double ronto = 1e-27;
static const double quecto = 1e-30;

// Extremes
// Forcing Frequency Extremes
// Assume that any forcing period larger than a Gyr leads to a zero in frequency.
// Converting to frequency is roughly 1.0e-17 rads s-1
static const double MIN_FREQUENCY = 1.0e-17;
// Assume max frequency is for a forcing period of 1 micro-second
static const double MAX_FREQUENCY = 1.0e8;

/// Shear/Bulk Modulus Extremes
static const double MIN_MODULUS = 1.0e-3;
