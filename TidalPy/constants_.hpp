#pragma once

#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

// Asks the compiler to inline a small per-point kernel into the sweep that calls it. MSVC otherwise leaves some of
// these as calls inside the larger vectorized loops, which costs a sweep about a sixth of its time.
#if defined(_MSC_VER)
#define TIDALPY_FORCE_INLINE __forceinline
#else
#define TIDALPY_FORCE_INLINE inline __attribute__((always_inline))
#endif


// TidalPy constants and runtime parameter struct
struct TidalPyConstants
{
    // Static members will be loaded in via TidalPy's configuration file.
    // True constants are left as const

    // Mathematics
    static constexpr double d_ppm = 1.0e-6;
    static constexpr double d_ppb = 1.0e-9;
    static constexpr double d_INF = std::numeric_limits<double>::infinity();
    static constexpr double d_PI  = std::numbers::pi;
    static constexpr double d_NAN = std::numeric_limits<double>::quiet_NaN();
    // Natural log of one half. The decay constant for a half life t_half is d_LN_HALF / t_half,
    // so this is negative: ln(1/2) = -ln(2).
    static constexpr double d_LN_HALF = -std::numbers::ln2;

    // Time.
    static constexpr double d_SECONDS_PER_MYR = 1.0e6 * 365.25 * 86400.0;

    // Computational
    static constexpr double d_DBL_MAX         = std::numeric_limits<double>::max();
    static constexpr double d_DBL_MIN         = std::numeric_limits<double>::min();
    static constexpr double d_DBL_MANT_DIGITS = std::numeric_limits<double>::digits;
    static constexpr double d_EPS             = std::numeric_limits<double>::epsilon();
    static constexpr double d_EPS_10          = 10 * std::numeric_limits<double>::epsilon();
    static constexpr double d_EPS_100         = 100 * std::numeric_limits<double>::epsilon();

    // Sun
    static constexpr double d_MASS_SOLAR       = 1.988435e30;
    static constexpr double d_RADIUS_SOLAR     = 6.957e8;
    static constexpr double d_LUMINOSITY_SOLAR = 3.828e26;  // IAU 2015 nominal [W]

    // TRAPPIST-1: So we have small M-dwarf parameter at hand. Set with Agol+ 2021 data
    static constexpr double d_MASS_TRAP1       = 0.0898 * d_MASS_SOLAR;
    static constexpr double d_RADIUS_TRAP1     = 0.1192 * d_RADIUS_SOLAR;
    static constexpr double d_LUMINOSITY_TRAP1 = 0.000553 * d_LUMINOSITY_SOLAR;

    // Earth
    static constexpr double d_MASS_EARTH   = 5.9721986e24;
    static constexpr double d_RADIUS_EARTH = 6.371008e6;

    // Jupiter
    static constexpr double d_MASS_JUPITER   = 1.898125e27;
    static constexpr double d_RADIUS_JUPITER = 69.911e6;  // IAU nominal mean radius [m]

    // Pluto
    static constexpr double d_MASS_PLUTO   = 1.309e22;
    static constexpr double d_RADIUS_PLUTO = 1.1899e6;

    // Io
    static constexpr double d_MASS_IO   = 8.9298e22;
    static constexpr double d_RADIUS_IO = 1.82149e6;
};

// Runtime Configurable Parameters (Mutable Shared State)
struct TidalPyConfig
{
    // Forcing Frequency Extremes
    double d_MIN_FREQUENCY; // Updated from TidalPy.config['numerical']['minimum_frequency']
    double d_MAX_FREQUENCY; // Updated from TidalPy.config['numerical']['maximum_frequency']

    // Material Extremes
    double d_MIN_MODULUS; // Updated from TidalPy.config['numerical']['minimum_modulus']
    // Rigidity mu / (rho g R) below which a melt-weakened solid is solved as a static liquid.
    double d_MIN_SOLID_RIGIDITY; // Updated from TidalPy.config['numerical']['minimum_solid_rigidity']

    // Planet Extremes
    double d_MIN_THICKNESS; // Updated from TidalPy.config['numerical']['minimum_layer_thickness']

    // Smallest magnitude a denominator may take before a guard substitutes it, shared by every module
    // that divides by a quantity which can reach zero (rheology, cooling, radiogenics).
    double d_NUMERICAL_FLOOR; // Updated from TidalPy.config['numerical']['numerical_floor']

    // Relative tolerance on layer-boundary continuity: how far a layer's inner radius may sit from the
    // previous layer's outer radius before the geometry is rejected.
    double d_LAYER_CONTINUITY_RTOL; // Updated from TidalPy.config['numerical']['layer_continuity_rtol']

    // Largest fraction of the planet radius a radial-solver integration may start from. Caps the
    // solver's own automatic choice and rejects a caller's starting radius above it.
    double d_MAX_START_RADIUS_FRAC; // Updated from TidalPy.config['numerical']['max_start_radius_fraction']

    // Smallest equilibrated reciprocal condition number the radial solver's surface boundary-condition system
    // may have; below it the system is singular to working precision and the solve fails.
    double d_MIN_SURFACE_RCOND; // Updated from TidalPy.config['numerical']['minimum_surface_rcond']

    // Relative tolerance within which two tidal-mode frequencies are the same one, and a frequency is zero.
    double d_FREQUENCY_MATCH_RTOL; // Updated from TidalPy.config['numerical']['frequency_match_rtol']

    // Smallest Nusselt number the convection cooling model reports.
    double d_MIN_NUSSELT; // Updated from TidalPy.config['numerical']['minimum_nusselt']

    // Largest factor by which a solved world's enclosed mass may differ from its stated mass before the EOS solve
    // fails as having no hydrostatic structure near that mass.
    double d_MAX_EOS_MASS_RATIO; // Updated from TidalPy.config['numerical']['maximum_eos_mass_ratio']

    // Density-from-pressure inversion of the compressible material EOS models (Birch-Murnaghan, Vinet): the
    // relative convergence tolerance and the iteration cap. Read when a model is built without its own values.
    double d_EOS_INVERT_RTOL;      // Updated from TidalPy.config['numerical']['eos_invert_rtol']
    int    d_EOS_INVERT_MAX_ITERS; // Updated from TidalPy.config['numerical']['eos_invert_max_iters']

    // Quadrature resolutions of the 3D tidal heating integrals (calc_3d_tides), -1 until the config is loaded:
    // the Gauss-Legendre order of the colatitude integral, the trapezoid nodes of the instantaneous longitude
    // integral, and the Gauss-Legendre nodes per layer of the radial integral.
    int d_TIDES_3D_LATITUDE_NODES;       // Updated from TidalPy.config['numerical']['tides_3d_latitude_nodes']
    int d_TIDES_3D_LONGITUDE_NODES;      // Updated from TidalPy.config['numerical']['tides_3d_longitude_nodes']
    int d_TIDES_3D_RADIAL_SLICES;        // Updated from TidalPy.config['numerical']['tides_3d_radial_slices']
    int d_TIDES_3D_MIN_RADII_PER_THREAD; // Updated from TidalPy.config['numerical']['tides_3d_min_radii_per_thread']
    int d_LOVE_SOLVE_THREADS;            // Updated from TidalPy.config['numerical']['love_solve_threads']
    int d_LOVE_SOLVE_MIN_PARALLEL;       // Updated from TidalPy.config['numerical']['love_solve_min_parallel']

    // Whole-planet EOS solve defaults, from TidalPy.config['eos_solver']. Read by every EOS solve that is
    // not handed an explicit value. The method is CyRK's ODEMethod enum as an int (-1 until the config is
    // loaded).
    int    d_EOS_SOLVER_METHOD;
    double d_EOS_SOLVER_RTOL;
    double d_EOS_SOLVER_ATOL;
    double d_EOS_SOLVER_PRESSURE_TOL;
    int    d_EOS_SOLVER_MAX_ITERS;
    int    d_EOS_SOLVER_SLICES_PER_LAYER;
    bool   d_EOS_SOLVER_NONDIMENSIONALIZE;
    bool   d_EOS_SOLVER_SOLVE_TEMPERATURE;// Updated from config['eos_solver']['solve_temperature']

    // Radial (Love number) solve defaults, from TidalPy.config['radial_solver']. Read by the world Love
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
        d_MIN_MODULUS = nan;
        d_MIN_SOLID_RIGIDITY = nan;
        d_MIN_THICKNESS = nan;
        d_NUMERICAL_FLOOR = nan;
        d_LAYER_CONTINUITY_RTOL = nan;
        d_MAX_START_RADIUS_FRAC = nan;
        d_MIN_SURFACE_RCOND = nan;
        d_FREQUENCY_MATCH_RTOL = nan;
        d_MIN_NUSSELT = nan;
        d_MAX_EOS_MASS_RATIO = nan;
        d_EOS_INVERT_RTOL = nan;
        d_EOS_INVERT_MAX_ITERS = -1;
        d_TIDES_3D_LATITUDE_NODES = -1;
        d_TIDES_3D_LONGITUDE_NODES = -1;
        d_TIDES_3D_RADIAL_SLICES = -1;
        d_TIDES_3D_MIN_RADII_PER_THREAD = -1;
        d_LOVE_SOLVE_THREADS = 1;
        d_LOVE_SOLVE_MIN_PARALLEL = -1;
        d_EOS_SOLVER_METHOD = -1;
        d_EOS_SOLVER_RTOL = nan;
        d_EOS_SOLVER_ATOL = nan;
        d_EOS_SOLVER_PRESSURE_TOL = nan;
        d_EOS_SOLVER_MAX_ITERS = -1;
        d_EOS_SOLVER_SLICES_PER_LAYER = -1;
        d_EOS_SOLVER_NONDIMENSIONALIZE = true;
        d_EOS_SOLVER_SOLVE_TEMPERATURE = true;
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

// Guards a denominator that may approach zero (a zero forcing frequency, say): a magnitude below the config's
// numerical floor is replaced by the floor, keeping the sign. Shared by every module that divides by such a quantity
// (rheology, cooling, radiogenics). An unwired config pointer leaves the value unchanged, as an unloaded (NaN) floor
// does.
// Thermal expansivity relative to its value at a reference density, alpha / alpha0, for an Anderson-Gruneisen
// parameter that itself falls with compression, delta_T = delta_T0 (rho0 / rho)^kappa (Anderson 1967; Chopelas and
// Boehler 1992). Integrating d ln(alpha) = delta_T d ln(V) gives
//     alpha / alpha0 = exp[(delta_T0 / kappa) ((rho0 / rho)^kappa - 1)],
// and the kappa -> 0 limit, a constant delta_T0, is (rho0 / rho)^delta_T0. 1 (a constant alpha) for delta_T0 = 0 or a
// density or reference density that is not positive.
inline double c_anderson_gruneisen_factor(
        double reference_density,
        double density,
        double delta_t0,
        double kappa) noexcept
{
    if ((delta_t0 == 0.0) || !(reference_density > 0.0) || !(density > 0.0)) {
        return 1.0;
    }
    const double expansion = reference_density / density;
    if (kappa == 0.0) {
        return std::pow(expansion, delta_t0);
    }
    return std::exp((delta_t0 / kappa) * (std::pow(expansion, kappa) - 1.0));
}

inline double c_guard_denominator(double value) noexcept
{
    if (tidalpy_config_ptr == nullptr) {
        return value;
    }
    const double floor_value = tidalpy_config_ptr->d_NUMERICAL_FLOOR;
    if (std::abs(value) < floor_value) {
        return (value < 0.0) ? -floor_value : floor_value;
    }
    return value;
}
