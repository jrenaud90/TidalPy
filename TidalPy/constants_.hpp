#pragma once

#include <limits>
#include <numbers>


// TidalPy constants and runtime parameter struct
struct TidalPyConstants
{
    // Static members will be loaded in via TidalPy's configuration file.
    // True constants are left as const

    // Mathematics
    static constexpr double d_ppm = 1.0e-6;
    static constexpr double d_ppb = 1.0e-9;
    static constexpr double d_INF = std::numeric_limits<double>::infinity();
    static constexpr double d_PI = std::numbers::pi;
    static constexpr double d_NAN = std::numeric_limits<double>::quiet_NaN();

    // Computational
    static constexpr double d_DBL_MAX = std::numeric_limits<double>::max();
    static constexpr double d_DBL_MIN = std::numeric_limits<double>::min();
    static constexpr double d_DBL_MANT_DIGITS = std::numeric_limits<double>::digits;
    static constexpr double d_EPS = std::numeric_limits<double>::epsilon();
    static constexpr double d_EPS_10 = 10 * std::numeric_limits<double>::epsilon();
    static constexpr double d_EPS_100 = 100 * std::numeric_limits<double>::epsilon();

    // Sun
    static constexpr double d_MASS_SOLAR = 1.988435e30;
    static constexpr double d_RADIUS_SOLAR = 6.957e8;
    static constexpr double d_LUMINOSITY_SOLAR = 3.848e26;

    // TRAPPIST-1: So we have small M-dwarf parameter at hand. Set with Agol+ 2021 data
    static constexpr double d_MASS_TRAP1 = 0.0898 * 1.988435e30;
    static constexpr double d_RADIUS_TRAP1 = 0.1192 * 6.957e8;
    static constexpr double d_LUMINOSITY_TRAP1 = 0.000553 * 3.848e26;

    // Earth
    static constexpr double d_MASS_EARTH = 5.9721986e24;
    static constexpr double d_RADIUS_EARTH = 6.371008e6;

    // Jupiter
    static constexpr double d_MASS_JUPITER = 1.898125e27;
    static constexpr double d_RADIUS_JUPITER = 69.8860e6;

    // Pluto
    static constexpr double d_MASS_PLUTO = 1.309e22;
    static constexpr double d_RADIUS_PLUTO = 1.1899e6;

    // Io
    static constexpr double d_MASS_IO  = 8.9298e22;
    static constexpr double d_RADIUS_IO = 1.82149e6;
};

// Runtime Configurable Parameters (Mutable Shared State)
struct TidalPyConfig
{
    // Forcing Frequency Extremes
    double d_MIN_FREQUENCY; // Updated from TidalPy.config['tides']['modes']['minimum_frequency']
    double d_MAX_FREQUENCY; // Updated from TidalPy.config['tides']['modes']['maximum_frequency']
    double d_MIN_SPIN_ORBIT_DIFF;// Updated from TidalPy.config['tides']['modes']['min_spin_orbit_diff']

    // Material Extremes    
    double d_MIN_VISCOSITY;// Updated from TidalPy.config['physics']['materials']['minimum_viscosity']
    double d_MIN_MODULUS;// Updated from TidalPy.config['physics']['materials']['minimum_modulus']

    // Planet Extremes
    double d_MIN_THICKNESS; // Updated from TidalPy.config['layers']['minimum_layer_thickness']

    // Astro / Physics Constants
    // The below are updated from SciPy
    double d_G;
    double d_AU;
    double d_SBC;
    double d_R;
    double d_K_BOLTZMAN;
    
    double d_TEST_CONST;

    // Constructor to set defaults (NaNs)
    TidalPyConfig() {
        double nan = std::numeric_limits<double>::quiet_NaN();
        d_MIN_FREQUENCY = nan;
        d_MAX_FREQUENCY = nan;
        d_MIN_SPIN_ORBIT_DIFF = nan;
        d_MIN_VISCOSITY = nan;
        d_MIN_MODULUS = nan;
        d_MIN_THICKNESS = nan;
        d_G = nan;
        d_AU = nan;
        d_SBC = nan;
        d_R = nan;
        d_K_BOLTZMAN = nan;
        d_TEST_CONST = nan;
    }
};

inline TidalPyConfig* tidalpy_config_ptr;