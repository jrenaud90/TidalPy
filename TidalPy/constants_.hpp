#pragma once

#include <limits>
#include <numbers>
#include <stdexcept>


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
    // Natural log of one half. The decay constant for a half life t_half is d_LN_HALF / t_half,
    // so this is negative: ln(1/2) = -ln(2).
    static constexpr double d_LN_HALF = -std::numbers::ln2;

    // Time.
    static constexpr double d_SECONDS_PER_MYR = 1.0e6 * 365.25 * 86400.0;

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
    static constexpr double d_LUMINOSITY_SOLAR = 3.828e26;  // IAU 2015 nominal [W]

    // TRAPPIST-1: So we have small M-dwarf parameter at hand. Set with Agol+ 2021 data
    static constexpr double d_MASS_TRAP1 = 0.0898 * d_MASS_SOLAR;
    static constexpr double d_RADIUS_TRAP1 = 0.1192 * d_RADIUS_SOLAR;
    static constexpr double d_LUMINOSITY_TRAP1 = 0.000553 * d_LUMINOSITY_SOLAR;

    // Earth
    static constexpr double d_MASS_EARTH = 5.9721986e24;
    static constexpr double d_RADIUS_EARTH = 6.371008e6;

    // Jupiter
    static constexpr double d_MASS_JUPITER = 1.898125e27;
    static constexpr double d_RADIUS_JUPITER = 69.911e6;  // IAU nominal mean radius [m]

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

    // Smallest magnitude a denominator may take before a guard substitutes it, shared by every module
    // that divides by a quantity which can reach zero (rheology, cooling, radiogenics).
    double d_NUMERICAL_FLOOR; // Updated from TidalPy.config_x['numerical']['numerical_floor']

    // Relative tolerance on layer-boundary continuity: how far a layer's inner radius may sit from the
    // previous layer's outer radius before the geometry is rejected.
    double d_LAYER_CONTINUITY_RTOL; // Updated from TidalPy.config_x['numerical']['layer_continuity_rtol']

    // Largest fraction of the planet radius a radial-solver integration may start from. Caps the
    // solver's own automatic choice and rejects a caller's starting radius above it.
    double d_MAX_START_RADIUS_FRAC; // Updated from TidalPy.config_x['numerical']['max_start_radius_fraction']

    // Whole-planet EOS solve defaults, from TidalPy.config_x['eos_solver']. Read by every EOS solve that is
    // not handed an explicit value. The method is CyRK's ODEMethod enum as an int (-1 until the config is
    // loaded).
    int    d_EOS_SOLVER_METHOD;
    double d_EOS_SOLVER_RTOL;
    double d_EOS_SOLVER_ATOL;
    double d_EOS_SOLVER_PRESSURE_TOL;
    int    d_EOS_SOLVER_MAX_ITERS;
    int    d_EOS_SOLVER_SLICES_PER_LAYER;
    bool   d_EOS_SOLVER_NONDIMENSIONALIZE;

    // Radial (Love number) solve defaults, from TidalPy.config_x['radial_solver']. Read by the world Love
    // solves, the tide paths that build their own, and the standalone radial_solver.
    int    d_RADIAL_SOLVER_METHOD;
    double d_RADIAL_SOLVER_RTOL;
    double d_RADIAL_SOLVER_ATOL;
    bool   d_RADIAL_SOLVER_USE_KAMATA;
    double d_RADIAL_SOLVER_START_RADIUS_TOL;
    bool   d_RADIAL_SOLVER_SCALE_RTOLS;
    int    d_RADIAL_SOLVER_MAX_NUM_STEPS;
    int    d_RADIAL_SOLVER_EXPECTED_SIZE;
    int    d_RADIAL_SOLVER_MAX_RAM_MB;
    bool   d_RADIAL_SOLVER_NONDIMENSIONALIZE;

    // Astro / Physics Constants
    // The below are updated from SciPy
    double d_G;
    double d_AU;
    double d_SBC;
    double d_R;
    double d_K_BOLTZMANN;
    
    double d_TEST_CONST;

    TidalPyConfig() {
        double nan = std::numeric_limits<double>::quiet_NaN();
        d_MIN_FREQUENCY = nan;
        d_MAX_FREQUENCY = nan;
        d_MIN_SPIN_ORBIT_DIFF = nan;
        d_MIN_VISCOSITY = nan;
        d_MIN_MODULUS = nan;
        d_MIN_THICKNESS = nan;
        d_NUMERICAL_FLOOR = nan;
        d_LAYER_CONTINUITY_RTOL = nan;
        d_MAX_START_RADIUS_FRAC = nan;
        d_EOS_SOLVER_METHOD = -1;
        d_EOS_SOLVER_RTOL = nan;
        d_EOS_SOLVER_ATOL = nan;
        d_EOS_SOLVER_PRESSURE_TOL = nan;
        d_EOS_SOLVER_MAX_ITERS = -1;
        d_EOS_SOLVER_SLICES_PER_LAYER = -1;
        d_EOS_SOLVER_NONDIMENSIONALIZE = true;
        d_RADIAL_SOLVER_METHOD = -1;
        d_RADIAL_SOLVER_RTOL = nan;
        d_RADIAL_SOLVER_ATOL = nan;
        d_RADIAL_SOLVER_USE_KAMATA = false;
        d_RADIAL_SOLVER_START_RADIUS_TOL = nan;
        d_RADIAL_SOLVER_SCALE_RTOLS = false;
        d_RADIAL_SOLVER_MAX_NUM_STEPS = -1;
        d_RADIAL_SOLVER_EXPECTED_SIZE = -1;
        d_RADIAL_SOLVER_MAX_RAM_MB = -1;
        d_RADIAL_SOLVER_NONDIMENSIONALIZE = true;
        d_G = nan;
        d_AU = nan;
        d_SBC = nan;
        d_R = nan;
        d_K_BOLTZMANN = nan;
        d_TEST_CONST = nan;
    }
};

// C++17 inline variable guarantees only one instance per extension/DLL
inline TidalPyConfig* tidalpy_config_ptr = nullptr;

// Cython will call this to inject the correct memory address
inline void set_tidalpy_config_ptr(TidalPyConfig* ptr)
{
    tidalpy_config_ptr = ptr;
}

// Newton's constant from the shared runtime config (populated from SciPy at package initialization). Returns NaN
// when the pointer was never wired so a missing initialization shows up in results instead of being masked by a
// hard-coded literal.
inline double c_get_G() noexcept
{
    return (tidalpy_config_ptr != nullptr) ? tidalpy_config_ptr->d_G : TidalPyConstants::d_NAN;
}
