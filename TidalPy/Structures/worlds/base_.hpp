#pragma once
/*
 * base_.hpp: c_BaseWorld, the world class every TidalPy world type derives from (extends c_StructureBase).
 *
 * A world owns an ordered stack of layers, inner to outer with index 0 innermost, which may be empty: a world with
 * no layers still dissipates tidally through the analytic tide models. It holds the world-level identification and
 * the orbital and thermal scalars (albedo, emissivity, obliquity, spin frequency), bulk geometry and the
 * equilibrium temperature, the spin model, the tide model and its configuration, the whole-planet aggregates and
 * geometry validation, and the whole-planet EOS and radial (Love number) solves, which walk every layer. Every
 * world owns a call lock (c_WorldCallLock) that its solves, its tide setters, and the reads of its results take.
 * All MKS: radius [m], mass [kg], angles [rad], frequency [rad/s].
 *
 * Binary payload: radius and mass (c_StructureBase), the name, world type, albedo, emissivity, obliquity, and spin
 * frequency, then the tide section: min_degree_l, max_degree_l, eccentricity_truncation, obliquity_truncation, and
 * love_method (int32_t each), layer_tidal_heating (uint8_t), love_fixed_q, love_fixed_dt, and
 * eccentricity_exact_tolerance, and the tide model behind a presence flag; then the layer count, the spin model's
 * moment-of-inertia factor, and the pinned [eos_solver] and [radial_solver] settings (each key a presence flag, then
 * its value when set), then each layer's own complete binary record, in index order. Solved state (the EOS profile,
 * Love numbers, tide results) is not saved.
 */

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <istream>
#include <iomanip>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "structure_base_.hpp"
#include "../layers/call_lock_.hpp"   // c_WorldCallLock
#include "../layers/factory_.hpp"

#include "constants_.hpp"
#include "../../Dynamics/spin_.hpp"
#include "../../Material/eos/solver_.hpp"
#include "../../Material/eos/methods/material_preeval_.hpp"
#include "thermal_layout_.hpp"

// The RadialSolver sub-modules compile into this translation unit so the shooting CyRK integration runs in
// the extension that owns the CySolverResult objects, with no cross-extension call().
#include "../../Utilities/math/numerics_.hpp"
#include "../../Utilities/dimensions/nondimensional_.hpp"
#include "../../Utilities/arrays/interp_.hpp"
#include "../../RadialSolver/rs_constants_.hpp"
#include "../../RadialSolver/rs_solution_.hpp"
#include "../../RadialSolver/love_.hpp"
#include "../../RadialSolver/shooting_.hpp"
#include "../../RadialSolver/world_radial_solver_.hpp"

// The tide model hierarchy and its config and result structs are light, with no global-potential tables.
// calc_tides, which drives the heavy engine and collapse, is defined out-of-line in world_tides_.hpp, so the
// eccentricity and obliquity tables compile only into the extensions that call it.
#include "../../Tides/classes/tide_base_.hpp"     // tidalpy::c_TideBase, c_LoveNumbers
#include "../../Tides/classes/tide_result_.hpp"   // c_TideConfig, c_TideSolveConfig, c_GlobalTideResult
// tide_.hpp gives c_Grid3DAxes for the 3D grid declarations and c_tide_from_binary for the tide model in the world
// record. It adds no includes this header does not already have and forward-declares c_BaseWorld, so the
// dependency stays one-directional.
#include "../../Tides/classes/tide_.hpp"
#include "../../Tides/potential/truncation_warnings_.hpp"   // c_warn_tide_truncations
#include "../../Tides/love/love_method_.hpp"
// Relative paths, not bare names, so every extension that includes base_.hpp resolves these without
// needing Utilities/lookups on its include path.
#include "../../Utilities/lookups/keys_.hpp"
#include "../../Utilities/lookups/intmap_.hpp"

namespace tidalpy {

// Construction parameters for c_BaseWorld and its subclasses.
struct c_WorldConfig {
    std::string name;
    std::string world_type_str = "world";  // "star", "gasgiant", "terrestrial", ...
    double      radius     = 0.0;      // [m]
    double      mass       = 0.0;      // [kg]
    double      albedo     = 0.3;      // [dimensionless]
    double      emissivity = 1.0;      // [dimensionless]
    double      obliquity  = 0.0;      // [rad]
    double      spin_frequency = 0.0;  // [rad/s]
};

// Every EOS and Love-number solve starts from the [eos_solver] and [radial_solver] sections of the shared
// config, and a caller overrides only the fields it passes. The member initializers of the two structs below
// stand in only when the config has not been loaded, which does not happen after a normal import.
inline ODEMethod c_ode_method_from_config(int method_int, ODEMethod fallback) noexcept {
    if (method_int > static_cast<int>(ODEMethod::RK_BASE_METHOD)
        && method_int <= static_cast<int>(ODEMethod::RADAU)) {
        return static_cast<ODEMethod>(method_int);
    }
    return fallback;
}

// True once update_constants has filled the solver sections of the shared config.
inline bool c_solver_config_loaded() noexcept {
    return (tidalpy_config_ptr != nullptr) && (tidalpy_config_ptr->d_EOS_SOLVER_METHOD >= 0);
}

// Grouped so c_BaseWorld::solve_eos takes a single argument.
struct c_WorldEOSSolveConfig {
    double    surface_pressure    = 0.0;                 // [Pa]
    size_t    slices_per_layer    = 100;                 // radial samples per layer (>= 2)
    double    G_to_use            = -1.0;                // [m^3 kg^-1 s^-2]; < 0 -> TidalPy config G
    ODEMethod integration_method  = ODEMethod::DOP853;
    double    rtol                = 1.0e-10;
    double    atol                = 1.0e-14;
    double    pressure_tol        = 1.0e-8;              // relative to the central-pressure scale
    size_t    max_iters           = 100;
    bool      nondimensionalize   = true;                // integrate in non-dimensional units
    // Overrides the temperature each layer carries; NaN uses each layer's own.
    double    temperature         = TidalPyConstants::d_NAN;
    // Carry temperature and heat flow through the solve. Without a temperature contrast there is no profile
    // to integrate, so the solve keeps its four structure variables whatever this says.
    bool      solve_temperature   = true;
    // Temperature [K] the outermost layer radiates to. NaN leaves no flow through the surface.
    double    surface_temperature = TidalPyConstants::d_NAN;
    // On the clock the radiogenics models share; NaN takes each model's own reference time.
    double    time                = TidalPyConstants::d_NAN;
    // Cap on the passes that relax the boundary layers and interface flows against the structure, and the
    // relative change in the interface temperatures and flows that ends them ([eos_solver] max_thermal_passes and
    // thermal_tol).
    size_t    max_thermal_passes  = 12;
    double    thermal_tol         = 1.0e-8;
    // Whether this solve redefines the mass every layer holding its mass holds (the mass inside its current radii);
    // otherwise its configured mass, or the world's first solve, sets it.
    bool      reset_layer_masses  = false;
    bool      verbose             = false;

    c_WorldEOSSolveConfig() {
        if (!c_solver_config_loaded()) { return; }
        const TidalPyConfig& config = *tidalpy_config_ptr;
        this->slices_per_layer   = static_cast<size_t>(config.d_EOS_SOLVER_SLICES_PER_LAYER);
        this->integration_method = c_ode_method_from_config(config.d_EOS_SOLVER_METHOD, this->integration_method);
        this->rtol               = config.d_EOS_SOLVER_RTOL;
        this->atol               = config.d_EOS_SOLVER_ATOL;
        this->pressure_tol       = config.d_EOS_SOLVER_PRESSURE_TOL;
        this->max_iters          = static_cast<size_t>(config.d_EOS_SOLVER_MAX_ITERS);
        this->nondimensionalize  = config.d_EOS_SOLVER_NONDIMENSIONALIZE;
        this->solve_temperature  = config.d_EOS_SOLVER_SOLVE_TEMPERATURE;
        this->max_thermal_passes = static_cast<size_t>(config.d_EOS_SOLVER_MAX_THERMAL_PASSES);
        this->thermal_tol        = config.d_EOS_SOLVER_THERMAL_TOL;
    }
};

// Parameters for the whole-planet Love-number solve.
struct c_LoveSolveConfig {
    double    frequency = 1.0e-5;   // [rad/s]; tidal forcing frequency
    int       degree_l  = 2;        // harmonic degree
    // In order; 1 = tidal, 2 = loading, 0 = free.
    std::vector<int> bc_models = {1};
    // c_LoveMethod as int: 0 radial_solver, 1 propagation_matrix,
    // 2 homogeneous, 3 cpl, 4 ctl, 5 laterally_inhomogeneous
    int       love_method = 0;
    double    fixed_q            = TidalPyConstants::d_NAN;   // cpl quality factor (NaN: from the tide model)
    double    fixed_dt           = TidalPyConstants::d_NAN;   // ctl time lag [s] (NaN: from the tide model)
    int       core_model         = 0;                         // propagation-matrix core starting condition (0-4)
    int       starting_method    = 0;                         // c_StartingMethod (shooting method only)
    int       degree1_frame      = 0;                         // c_Degree1Frame (degree-1 loading only)
    bool      nondimensionalize  = true;
    double    starting_radius    = 0.0;                       // [m]; 0 -> auto
    double    start_radius_tol   = 1.0e-5;
    ODEMethod integration_method = ODEMethod::DOP853;
    double    rtol               = 3.0e-8;
    double    atol               = 3.0e-8;
    bool      scale_rtols        = false;
    size_t    max_num_steps      = 500000;
    size_t    expected_size      = 128;
    size_t    max_ram_MB         = 500;
    double    max_step           = 0.0;
    bool      verbose            = false;
    bool      warnings           = true;
    // Shooting method: keep only what the Love numbers need, so the radial functions below the surface are not
    // available afterwards. The internal tide solves leave it off.
    bool      love_only          = false;

    c_LoveSolveConfig() {
        if (!c_solver_config_loaded()) { return; }
        const TidalPyConfig& config = *tidalpy_config_ptr;
        this->starting_method    = config.d_RADIAL_SOLVER_STARTING_METHOD;
        this->degree1_frame      = config.d_RADIAL_SOLVER_DEGREE1_FRAME;
        this->nondimensionalize  = config.d_RADIAL_SOLVER_NONDIMENSIONALIZE;
        this->start_radius_tol   = config.d_RADIAL_SOLVER_START_RADIUS_TOL;
        this->integration_method = c_ode_method_from_config(config.d_RADIAL_SOLVER_METHOD, this->integration_method);
        this->rtol               = config.d_RADIAL_SOLVER_RTOL;
        this->atol               = config.d_RADIAL_SOLVER_ATOL;
        this->scale_rtols        = config.d_RADIAL_SOLVER_SCALE_RTOLS;
        this->max_num_steps      = static_cast<size_t>(config.d_RADIAL_SOLVER_MAX_NUM_STEPS);
        this->expected_size      = static_cast<size_t>(config.d_RADIAL_SOLVER_EXPECTED_SIZE);
        this->max_ram_MB         = static_cast<size_t>(config.d_RADIAL_SOLVER_MAX_RAM_MB);
    }

    // From a plain array, for the Cython wrappers.
    void set_bc_models(const int* models_ptr, size_t num_models) {
        if (models_ptr == nullptr || num_models == 0) { return; }
        this->bc_models.assign(models_ptr, models_ptr + num_models);
    }
};

// The [eos_solver] and [radial_solver] keys a world file may pin (c_BaseWorld::set_eos_solver_overrides and
// set_radial_solver_overrides), each a row of its section's table: the TidalPy_Configs.toml key, how the value is kept
// and written in a binary record, and how it applies to the section's solve configuration. A pinned key replaces the
// configuration's value for every solve the world runs; a call's own argument still wins over it, and a key left out
// keeps following the configuration. The table is the one place a key is listed in C++: apply, the binary record, and
// the Python set_solver_defaults and get_solver_defaults all iterate it, and its keys must match schema.py's
// _SOLVER_KEY_RULES (a test checks).
enum class c_SolverSettingKind : uint8_t {
    Method = 0,   // an ODEMethod (int32 in a binary record)
    Real   = 1,   // a double
    Count  = 2,   // a non-negative integer (uint64 in a binary record)
    Flag   = 3,   // a switch (uint8 in a binary record)
    StartingMethod = 4,   // a c_StartingMethod (int32 in a binary record)
    Degree1Frame   = 5,   // a c_Degree1Frame (int32 in a binary record)
};

template <class Config>
struct c_SolverSettingSpec {
    const char*         key;
    c_SolverSettingKind kind;
    void (*apply)(Config& cfg, double value);
};

struct c_EOSSolverSection {
    using Config = c_WorldEOSSolveConfig;
    static const std::vector<c_SolverSettingSpec<Config>>& specs() {
        static const std::vector<c_SolverSettingSpec<Config>> table = {
            {"integration_method", c_SolverSettingKind::Method,
             [](Config& cfg, double value) {
                 cfg.integration_method = static_cast<ODEMethod>(static_cast<int>(value));
             }},
            {"rtol", c_SolverSettingKind::Real, [](Config& cfg, double value) { cfg.rtol = value; }},
            {"atol", c_SolverSettingKind::Real, [](Config& cfg, double value) { cfg.atol = value; }},
            {"pressure_tol", c_SolverSettingKind::Real, [](Config& cfg, double value) { cfg.pressure_tol = value; }},
            {"max_iters", c_SolverSettingKind::Count,
             [](Config& cfg, double value) { cfg.max_iters = static_cast<size_t>(value); }},
            {"slices_per_layer", c_SolverSettingKind::Count,
             [](Config& cfg, double value) { cfg.slices_per_layer = static_cast<size_t>(value); }},
            {"nondimensionalize", c_SolverSettingKind::Flag,
             [](Config& cfg, double value) { cfg.nondimensionalize = (value != 0.0); }},
            {"solve_temperature", c_SolverSettingKind::Flag,
             [](Config& cfg, double value) { cfg.solve_temperature = (value != 0.0); }},
            {"max_thermal_passes", c_SolverSettingKind::Count,
             [](Config& cfg, double value) { cfg.max_thermal_passes = static_cast<size_t>(value); }},
            {"thermal_tol", c_SolverSettingKind::Real, [](Config& cfg, double value) { cfg.thermal_tol = value; }},
        };
        return table;
    }
};

struct c_RadialSolverSection {
    using Config = c_LoveSolveConfig;
    static const std::vector<c_SolverSettingSpec<Config>>& specs() {
        static const std::vector<c_SolverSettingSpec<Config>> table = {
            {"integration_method", c_SolverSettingKind::Method,
             [](Config& cfg, double value) {
                 cfg.integration_method = static_cast<ODEMethod>(static_cast<int>(value));
             }},
            {"rtol", c_SolverSettingKind::Real, [](Config& cfg, double value) { cfg.rtol = value; }},
            {"atol", c_SolverSettingKind::Real, [](Config& cfg, double value) { cfg.atol = value; }},
            {"starting_method", c_SolverSettingKind::StartingMethod,
             [](Config& cfg, double value) { cfg.starting_method = static_cast<int>(value); }},
            {"degree1_frame", c_SolverSettingKind::Degree1Frame,
             [](Config& cfg, double value) { cfg.degree1_frame = static_cast<int>(value); }},
            {"start_radius_tolerance", c_SolverSettingKind::Real,
             [](Config& cfg, double value) { cfg.start_radius_tol = value; }},
            {"scale_rtols", c_SolverSettingKind::Flag,
             [](Config& cfg, double value) { cfg.scale_rtols = (value != 0.0); }},
            {"max_num_steps", c_SolverSettingKind::Count,
             [](Config& cfg, double value) { cfg.max_num_steps = static_cast<size_t>(value); }},
            {"expected_size", c_SolverSettingKind::Count,
             [](Config& cfg, double value) { cfg.expected_size = static_cast<size_t>(value); }},
            {"max_ram_mb", c_SolverSettingKind::Count,
             [](Config& cfg, double value) { cfg.max_ram_MB = static_cast<size_t>(value); }},
            {"nondimensionalize", c_SolverSettingKind::Flag,
             [](Config& cfg, double value) { cfg.nondimensionalize = (value != 0.0); }},
        };
        return table;
    }
};

// The keys of one section a world pins, held by row of the section's table. A value is kept as a double whatever its
// kind (an ODEMethod's enum value, a count, 0 or 1 for a switch), exactly for every value a setting can take.
template <class Section>
class c_SolverOverrides {
public:
    using Config = typename Section::Config;

    void apply(Config& cfg) const noexcept {
        const auto& specs = Section::specs();
        for (std::size_t row = 0; row < specs.size(); ++row) {
            if (this->p_values[row]) { specs[row].apply(cfg, *this->p_values[row]); }
        }
    }

    // Pins a key; throws std::invalid_argument for a key the section does not have.
    void set(const std::string& key, double value) { this->p_values[p_row(key)] = value; }
    bool has(const std::string& key) const { return this->p_values[p_row(key)].has_value(); }
    // The pinned value; NaN when the key is not pinned.
    double get(const std::string& key) const {
        const auto& value = this->p_values[p_row(key)];
        return value ? *value : TidalPyConstants::d_NAN;
    }

    static std::vector<std::string> keys() {
        std::vector<std::string> out;
        for (const auto& spec : Section::specs()) { out.emplace_back(spec.key); }
        return out;
    }
    // The kind of a key (c_SolverSettingKind as an int); throws for an unknown key.
    static int kind(const std::string& key) { return static_cast<int>(Section::specs()[p_row(key)].kind); }

    // Each row in table order: a presence byte, then the value at its kind's fixed width.
    void write(std::ostream& out) const {
        const auto& specs = Section::specs();
        for (std::size_t row = 0; row < specs.size(); ++row) {
            const std::optional<double>& value = this->p_values[row];
            const uint8_t present = value.has_value() ? 1 : 0;
            out.write(reinterpret_cast<const char*>(&present), sizeof(uint8_t));
            if (!present) { continue; }
            switch (specs[row].kind) {
                case c_SolverSettingKind::Method: p_write_as<int32_t>(out, *value);  break;
                case c_SolverSettingKind::Real:   p_write_as<double>(out, *value);   break;
                case c_SolverSettingKind::Count:  p_write_as<uint64_t>(out, *value); break;
                case c_SolverSettingKind::Flag:   p_write_as<uint8_t>(out, *value);  break;
                case c_SolverSettingKind::StartingMethod: p_write_as<int32_t>(out, *value); break;
                case c_SolverSettingKind::Degree1Frame:   p_write_as<int32_t>(out, *value); break;
            }
        }
    }
    void read(std::istream& in) {
        const auto& specs = Section::specs();
        for (std::size_t row = 0; row < specs.size(); ++row) {
            uint8_t present = 0;
            in.read(reinterpret_cast<char*>(&present), sizeof(uint8_t));
            if (!present) { this->p_values[row].reset(); continue; }
            switch (specs[row].kind) {
                case c_SolverSettingKind::Method: this->p_values[row] = p_read_as<int32_t>(in);  break;
                case c_SolverSettingKind::Real:   this->p_values[row] = p_read_as<double>(in);   break;
                case c_SolverSettingKind::Count:  this->p_values[row] = p_read_as<uint64_t>(in); break;
                case c_SolverSettingKind::Flag:   this->p_values[row] = p_read_as<uint8_t>(in);  break;
                case c_SolverSettingKind::StartingMethod: this->p_values[row] = p_read_as<int32_t>(in); break;
                case c_SolverSettingKind::Degree1Frame:   this->p_values[row] = p_read_as<int32_t>(in); break;
            }
        }
    }

private:
    static std::size_t p_row(const std::string& key) {
        const auto& specs = Section::specs();
        for (std::size_t row = 0; row < specs.size(); ++row) {
            if (key == specs[row].key) { return row; }
        }
        throw std::invalid_argument("TidalPy: unknown solver setting '" + key + "'.");
    }
    template <typename Stored>
    static void p_write_as(std::ostream& out, double value) {
        const Stored stored = static_cast<Stored>(value);
        out.write(reinterpret_cast<const char*>(&stored), sizeof(Stored));
    }
    template <typename Stored>
    static double p_read_as(std::istream& in) {
        Stored stored {};
        in.read(reinterpret_cast<char*>(&stored), sizeof(Stored));
        return static_cast<double>(stored);
    }

    std::vector<std::optional<double>> p_values = std::vector<std::optional<double>>(Section::specs().size());
};

using c_EOSSolverOverrides    = c_SolverOverrides<c_EOSSolverSection>;
using c_RadialSolverOverrides = c_SolverOverrides<c_RadialSolverSection>;

// c_WorldCallLock, the world's call lock, is defined in layers/call_lock_.hpp: a world's layers take it too. The mutex
// it locks is c_BaseWorld::p_call_mutex.

// A copy of everything a solve_eos result reports, taken under the world's call lock so that another thread's
// solve_eos cannot replace the solution while it is read. The profile arrays are empty when no solve has populated
// the world; the per-layer entries are filled only after a completed solve (get_eos_solved).
struct c_WorldEOSReport {
    bool        solved             = false;
    bool        success            = false;
    std::string message;
    int         iterations         = 0;
    // Integrations of the whole structure over every thermal pass, the repeats that keep a converged pass included.
    int         structure_integrations = 0;
    bool        max_iters_hit      = false;
    double      pressure_error     = TidalPyConstants::d_NAN;  // [Pa]
    double      surface_gravity    = TidalPyConstants::d_NAN;  // [m s-2]
    double      surface_pressure   = TidalPyConstants::d_NAN;  // [Pa]
    double      central_pressure   = TidalPyConstants::d_NAN;  // [Pa]
    double      planet_mass        = TidalPyConstants::d_NAN;  // [kg]
    double      planet_moi         = TidalPyConstants::d_NAN;  // [kg m2]
    std::size_t thermal_passes     = 0;
    bool        thermal_converged  = false;

    // Radial profile of the solve (SI), one entry per reported sample.
    std::vector<double> radius;
    std::vector<double> gravity;
    std::vector<double> pressure;
    std::vector<double> mass;
    std::vector<double> moi;
    std::vector<double> density;
    std::vector<double> temperature;
    std::vector<double> heat_flow;

    // One entry per layer.
    std::vector<c_LayerThermal> layer_thermal;
    std::vector<double>         layer_radius_outer;      // [m]

    // The solid and liquid zones the solve found (c_BaseWorld::get_zones).
    std::vector<c_EOSZone> zones;
};

// What one EOS solve evaluates its materials with: every layer's material (shared, and immutable, so holding it costs
// nothing), the per-layer inputs the structure ODE reaches them through, and the heat sources of a thermal solve. The
// solution co-owns it, so a retained or exported solution keeps answering exactly as solved, whatever later happens
// to the layers.
struct c_EOSSolveState {
    std::vector<std::shared_ptr<const c_Material>> materials;
    std::vector<c_MaterialPreevalInput>            inputs;
    c_Heating                                      heating;
};

// A radial solve calc_tides has already run for one (degree, |omega|) group of its tidal modes. It lends these to the
// per-layer heating integral, which needs the same solves, so they are not run twice. Non-owning.
struct c_RetainedRadialSolve {
    int                              degree_l  = 0;
    double                           frequency = 0.0;       // |omega| [rad s-1]
    const ::c_RadialSolutionStorage* storage   = nullptr;
};

// Everything one global tidal solve gathers (c_BaseWorld::calc_tides), before any of it reaches the world: calc_tides
// commits it only once the whole solve has succeeded.
struct c_TideSolveOutcome {
    c_GlobalTideResult tide_result;
    // Per-mode Love numbers (k, h, l) keyed by the tidal mode (l, m, p, q); empty for the analytic models.
    c_IntMap<c_Key4, c_LoveNumbers> tide_love;
    // Each layer's share of the heating [W] (see calc_tides), NaN where none is resolved.
    std::vector<double> layer_heating;
    // What the tidal heat source spreads in later solves.
    c_TidalHeatingRecord heating_record;
};

// One tidal layer's part of a quasi-homogeneous Love solve (the homogeneous, cpl, and ctl methods): the Love numbers
// of a homogeneous planet made of the layer's averaged material, and the tidal scale the world weighs them by.
struct c_LayerLove {
    std::size_t          layer_index = 0;
    double               tidal_scale = 0.0;                               // dimensionless
    c_LoveNumbers        love;                                            // the layer's own Love numbers
    std::complex<double> shear_modulus = {TidalPyConstants::d_NAN, 0.0};  // its complex shear modulus [Pa]
    double               volume = 0.0;                                    // [m3]
};

// Everything one Love-number solve produces: the radial solver and its solution storage for the radial methods, or
// the Love numbers of the quasi-homogeneous methods. The world keeps one for solve_love_numbers and the getters that
// report it; calc_tides and the 3D paths solve into workspaces of their own, so they neither overwrite the world's
// last solve nor share solver state with it.
struct c_LoveWorkspace {
    bool         solved      = false;
    c_LoveMethod method_last = c_LoveMethod::RadialSolver;

    // The radial methods: the cached solver, the world layer each of its layers (the world's zones) belongs to, and
    // whether that layer takes the solid equations (1) or the liquid ones (0).
    std::unique_ptr<::c_WorldRadialSolver> radial_solver;
    std::vector<std::size_t>               radial_world_layer;
    std::vector<char>                      radial_solid_layer;

    // The quasi-homogeneous methods: the world's Love numbers (the tidal-scale-weighted sum of the layers') and
    // each tidal layer's part.
    bool                     analytic_success    = false;
    int                      analytic_error_code = -100;
    std::string              analytic_message    = "No love-number solve has been run.";
    c_LoveNumbers            analytic;
    std::vector<c_LayerLove> analytic_layers;

    bool is_analytic() const noexcept { return c_love_method_is_homogeneous(this->method_last); }

    // Forget the analytic results, so nothing reports a value from an earlier solve.
    void reset_analytic() noexcept {
        this->analytic_success    = false;
        this->analytic_error_code = -100;
        this->analytic            = c_LoveNumbers();
        this->analytic_layers.clear();
    }

    // Everything solved describes a structure that is about to change.
    void invalidate() noexcept {
        this->solved = false;
        this->reset_analytic();
        if (this->radial_solver) { this->radial_solver->invalidate(); }
    }

    // Non-owning; null until a radial solve builds the cache.
    const ::c_RadialSolutionStorage* get_storage() const noexcept {
        return this->radial_solver ? this->radial_solver->get_storage() : nullptr;
    }

    bool get_success() const noexcept {
        if (this->is_analytic()) { return this->analytic_success; }
        const auto* storage = this->get_storage();
        return (storage && this->solved) ? storage->success : false;
    }
    // The storage while it still describes the structure (c_WorldRadialSolver::get_storage_current), else null, so
    // the diagnostics of a solve on an earlier structure never read as the current structure's.
    const ::c_RadialSolutionStorage* get_current_storage() const noexcept {
        return (this->radial_solver && this->radial_solver->get_storage_current()) ? this->get_storage() : nullptr;
    }
    int get_error_code() const noexcept {
        if (this->is_analytic()) { return this->analytic_error_code; }
        const auto* storage = this->get_current_storage();
        return storage ? storage->error_code : -100;
    }
    const std::string& get_message() const noexcept {
        static const std::string no_message = "No love-number solve has been run.";
        if (this->is_analytic()) { return this->analytic_message; }
        const auto* storage = this->get_current_storage();
        return storage ? storage->message : no_message;
    }
    std::size_t get_num_ytypes() const noexcept {
        if (this->is_analytic()) { return this->analytic_success ? 1 : 0; }
        const auto* storage = this->get_current_storage();
        return storage ? storage->num_ytypes : 0;
    }
    double get_surface_amplification() const noexcept {
        if (this->is_analytic()) { return 0.0; }
        const auto* storage = this->get_current_storage();
        return storage ? storage->surface_amplification : 0.0;
    }
    double get_surface_rcond() const noexcept {
        if (this->is_analytic()) { return TidalPyConstants::d_NAN; }
        const auto* storage = this->get_current_storage();
        return storage ? storage->surface_rcond : TidalPyConstants::d_NAN;
    }
    double get_surface_frame_residual() const noexcept {
        if (this->is_analytic()) { return TidalPyConstants::d_NAN; }
        const auto* storage = this->get_current_storage();
        return storage ? storage->surface_frame_residual : TidalPyConstants::d_NAN;
    }

    // The Love numbers for a boundary-condition ytype; the analytic methods hold one tidal set at index 0. NaN when
    // no solve describes the current structure.
    c_LoveNumbers get_love(std::size_t ytype_idx = 0) const noexcept {
        const std::complex<double> nan_value(TidalPyConstants::d_NAN, 0.0);
        const c_LoveNumbers nan_love(nan_value, nan_value, nan_value);
        if (this->is_analytic()) {
            return (this->analytic_success && ytype_idx == 0) ? this->analytic : nan_love;
        }
        const auto* storage = this->get_storage();
        if (!storage || !this->solved || ytype_idx >= storage->complex_love_vec.size()) { return nan_love; }
        const ::c_LoveNumbers& solved_love = storage->complex_love_vec[ytype_idx];
        return c_LoveNumbers(solved_love.k, solved_love.h, solved_love.l);
    }

    // The radial functions (SI) at a radius [m] for a ytype and y-index (0..5 -> y1..y6); NaN when unsolved, after
    // an analytic solve (which has no radial functions), or where the solution does not reach.
    std::complex<double> get_radial_y(double radius, std::size_t ytype_idx, std::size_t y_idx) const noexcept {
        const std::complex<double> nan_value(TidalPyConstants::d_NAN, 0.0);
        if (this->is_analytic()) { return nan_value; }
        const auto* storage = this->get_storage();
        if (!storage || !this->solved || !storage->success || y_idx >= C_MAX_NUM_Y) { return nan_value; }
        std::complex<double> y_at_r[C_MAX_NUM_Y];
        if (!storage->get_radial_solution(radius, ytype_idx, y_at_r)) { return nan_value; }
        return y_at_r[y_idx];
    }
    // Whether the latest solve kept only its Love numbers and surface values (c_LoveSolveConfig::love_only).
    bool get_love_only() const noexcept {
        if (this->is_analytic()) { return false; }
        const auto* storage = this->get_storage();
        return storage && this->solved && storage->get_love_only();
    }
    std::complex<double> get_surface_y(std::size_t ytype_idx, std::size_t y_idx) const noexcept {
        const std::complex<double> nan_value(TidalPyConstants::d_NAN, 0.0);
        if (this->is_analytic()) { return nan_value; }
        const auto* storage = this->get_storage();
        if (!storage || !this->solved || !storage->success || y_idx >= C_MAX_NUM_Y) { return nan_value; }
        std::complex<double> surface_y[C_MAX_NUM_Y];
        if (!storage->get_surface_y(ytype_idx, surface_y)) { return nan_value; }
        return surface_y[y_idx];
    }
};

// Half-width of the central difference across a zone boundary that gives the latent heat it carries, as a fraction of
// the world radius: wide enough that the difference of two dense reads (good to the solve's rtol) resolves the slope
// of T - T_m, narrow enough that the curvature of the profile does not enter at that accuracy.
inline constexpr double C_LATENT_DIFFERENCE_FRACTION = 1.0e-5;

// The evaluation-layout entries (eos_layout_.hpp) a world reads from its own EOS solution rather than from
// the dense output of the layer that holds the radius.
inline bool c_is_world_eos_field(std::size_t field_index) noexcept {
    return (field_index == C_EOS_TEMPERATURE_INDEX) || (field_index == C_EOS_HEAT_FLOW_INDEX);
}

class c_BaseWorld : public c_StructureBase, public c_LayerOwner {
public:
    // Absolute gap allowed between a layer's inner radius and the previous layer's outer radius.
    static double layer_continuity_tol(double previous_outer_radius) noexcept {
        const double scale = (previous_outer_radius > 1.0) ? previous_outer_radius : 1.0;
        return tidalpy_config_ptr->d_LAYER_CONTINUITY_RTOL * scale;
    }

    c_BaseWorld() = default;

    explicit c_BaseWorld(const c_WorldConfig& cfg)
        : c_StructureBase(cfg.radius, cfg.mass),
          p_name(cfg.name),
          p_world_type(cfg.world_type_str),
          p_albedo(cfg.albedo),
          p_emissivity(cfg.emissivity),
          p_obliquity(cfg.obliquity),
          p_spin_frequency(cfg.spin_frequency)
    {}

    ~c_BaseWorld() override = default;

    // Getters (const, MKS)
    const std::string& get_name()             const noexcept { return this->p_name; }
    const std::string& get_world_type()       const noexcept { return this->p_world_type; }
    double             get_albedo()           const noexcept { return this->p_albedo; }
    double             get_emissivity()       const noexcept { return this->p_emissivity; }
    double             get_obliquity()        const noexcept { return this->p_obliquity; }
    double             get_spin_frequency()   const noexcept { return this->p_spin_frequency; }

    // Bulk geometry (const, MKS) from the world's own stored radius and mass.
    double calc_surface_gravity() const noexcept {
        return this->c_StructureBase::calc_surface_gravity(this->p_mass, this->p_radius);
    }
    double calc_escape_velocity() const noexcept {
        return this->c_StructureBase::calc_escape_velocity(this->p_mass, this->p_radius);
    }
    double calc_mean_density() const noexcept {
        return this->c_StructureBase::calc_mean_density(
            this->p_mass, this->calc_volume_sphere(this->p_radius));
    }

    // Fast-rotator radiative equilibrium temperature [K] over a uniform-temperature surface:
    //   T_eq = [(1 - A) * F / (4 * eps * sigma)]^(1/4)
    // with F the incident insolation flux [W/m^2], A the bond albedo, eps the emissivity, and sigma the
    // Stefan-Boltzmann constant. Returns 0.0 for non-positive flux or when the config pointer is null.
    double calc_equilibrium_temperature(double insolation_flux) const noexcept {
        if (insolation_flux <= 0.0 || tidalpy_config_ptr == nullptr) { return 0.0; }
        const double sigma = tidalpy_config_ptr->d_SBC;
        const double eps   = (this->p_emissivity > 0.0) ? this->p_emissivity : 1.0;
        const double absorbed = (1.0 - this->p_albedo) * insolation_flux;
        return std::pow(absorbed / (4.0 * eps * sigma), 0.25);
    }

    // Mutators
    void set_name(const std::string& name)      { this->p_name = name; }

    // The system this world belongs to, as a source of the orbital state its tides are raised in (non-owning;
    // null outside a system), and this world's index there. The system sets both when the world is added and
    // clears them when it goes away.
    void set_tide_state_provider(const c_TideStateProvider* provider_ptr, std::size_t world_index) noexcept {
        this->p_tide_state_provider_ptr = provider_ptr;
        this->p_tide_state_index        = world_index;
    }
    const c_TideStateProvider* get_tide_state_provider() const noexcept { return this->p_tide_state_provider_ptr; }
    std::size_t get_tide_state_index() const noexcept { return this->p_tide_state_index; }

    // The tidal state the world's system gives it; false outside a system, or with no tidal host or usable orbit.
    bool get_tide_state(c_TideSolveConfig& state_out) const {
        if (this->p_tide_state_provider_ptr == nullptr) { return false; }
        return this->p_tide_state_provider_ptr->get_tide_state(this->p_tide_state_index, state_out);
    }
    void set_spin_frequency(double freq) noexcept { this->p_spin_frequency = freq; }
    void set_obliquity(double obliq)        noexcept { this->p_obliquity = obliq; }

    // Global (1D) tidal dissipation. Attach a tide model and a tide config, then call calc_tides(orbital state) to
    // collapse the global tidal modes into the total heating and the three orbital potential derivatives. The
    // analytic models (cpl, ctl, ctl_q) work on any world, layers or not; the rheology model needs a solved EOS.
    // Both setters take the call lock, so neither frees or changes what a tidal solve on another thread is using.
    void set_tide_model(std::shared_ptr<const c_TideBase> tide) noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        this->p_tide         = std::move(tide);
        this->p_tides_solved = false;
    }
    // Through the generic model handle the Python wrappers hold; throws std::invalid_argument for another family.
    void set_tide_model_handle(const std::shared_ptr<c_PhysicsBase>& model) {
        this->set_tide_model(c_share_as<c_TideBase>(model, "a world's tide model"));
    }
    bool get_tide_model_set() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_tide != nullptr;
    }
    const c_TideBase* get_tide_model() const noexcept { return this->p_tide.get(); }
    // The attached model as the generic handle the Python wrappers hold; empty without one.
    std::shared_ptr<c_PhysicsBase> share_tide_model() const {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return c_share_physics_of(this->p_tide);
    }

    // Throws std::invalid_argument for a degree range outside 2 <= min <= max <= 10 (the tabulated degrees), a
    // non-positive Love-method Q or negative time lag (NaN leaves either unset), or an exact eccentricity tolerance
    // outside (0, 1).
    void set_tide_config(const c_TideConfig& cfg) {
        validate_tide_config(cfg);
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        this->p_tide_config  = cfg;
        this->p_tides_solved = false;
    }
    static void validate_tide_config(const c_TideConfig& cfg) {
        if (!((cfg.min_degree_l >= 2) && (cfg.min_degree_l <= cfg.max_degree_l) && (cfg.max_degree_l <= 10))) {
            throw std::invalid_argument(
                "TidalPy: the tidal degree range must satisfy 2 <= min_degree_l <= max_degree_l <= 10; got " +
                std::to_string(cfg.min_degree_l) + " to " + std::to_string(cfg.max_degree_l) + ".");
        }
        if (cfg.love_fixed_q <= 0.0) {
            throw std::invalid_argument("TidalPy: love_fixed_q must be positive (NaN leaves it unset).");
        }
        if (cfg.love_fixed_dt < 0.0) {
            throw std::invalid_argument("TidalPy: love_fixed_dt must not be negative (NaN leaves it unset).");
        }
        if (!((cfg.eccentricity_exact_tolerance > 0.0) && (cfg.eccentricity_exact_tolerance < 1.0))) {
            throw std::invalid_argument("TidalPy: eccentricity_exact_tolerance must be in (0, 1).");
        }
    }
    const c_TideConfig& get_tide_config() const noexcept { return this->p_tide_config; }

    // The results of the most recent calc_tides, each read under the call lock so it never mixes two solves.
    bool get_tides_solved() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_tides_solved;
    }
    double get_tidal_heating() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_tides_solved ? this->p_tide_result.tidal_heating : TidalPyConstants::d_NAN;
    }
    // The tidal heating spread over the world's surface [W m-2], heating / (4 pi R^2); NaN before a calc_tides.
    double get_tidal_heat_flux() const noexcept {
        return this->get_tidal_heating() / this->calc_surface_area(this->get_radius());
    }
    double get_tidal_dU_dM() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_tides_solved ? this->p_tide_result.dU_dM : TidalPyConstants::d_NAN;
    }
    double get_tidal_dU_dw() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_tides_solved ? this->p_tide_result.dU_dw : TidalPyConstants::d_NAN;
    }
    double get_tidal_dU_dO() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_tides_solved ? this->p_tide_result.dU_dO : TidalPyConstants::d_NAN;
    }
    // The per-mode sum of dU/dM - dU/dw, which keeps de/dt exact at small eccentricity where the two separate sums
    // nearly cancel.
    double get_tidal_dU_dM_minus_dw() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_tides_solved ? this->p_tide_result.dU_dM_minus_dw : TidalPyConstants::d_NAN;
    }
    int get_num_tidal_modes() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_tides_solved ? this->p_tide_result.num_modes : 0;
    }

    // The whole collapsed global tidal result (heating, the three potential derivatives, mode and error codes)
    // from the most recent calc_tides. Check get_tides_solved() first: these fields keep their unsolved
    // defaults of zero, where the scalar getters above return NaN. The reference is read unlocked, so a caller that
    // shares the world with other threads holds the call lock (get_call_mutex) from its calc_tides through the read.
    const c_GlobalTideResult& get_tide_result() const noexcept { return this->p_tide_result; }

    // Complex potential Love number k_l for the tidal mode (l, m, p, q) from the most recent
    // rheology calc_tides. NaN for the analytic models (no radial solution) or an inactive mode.
    std::complex<double> get_tidal_love_k(int degree_l, int m, int p, int q) const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        bool found = false;
        c_Key4 lmpq_key(static_cast<int16_t>(degree_l), static_cast<int16_t>(m),
                        static_cast<int16_t>(p), static_cast<int16_t>(q));
        c_LoveNumbers love = this->p_tide_solver_love.get(found, lmpq_key);
        if (!found) {
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        }
        return love.k;
    }

    // The world's call lock (c_WorldCallLock), for a caller that must hold it across several calls, such as a
    // calc_tides and the read of its result (get_tide_result). Non-owning; the world owns it.
    std::recursive_mutex* get_call_mutex() const noexcept { return this->p_call_mutex.get(); }

    // Add a layer, inner to outer; its inner radius must match the current outermost radius, 0 for the
    // first layer, within the continuity tolerance.
    void add_layer(std::unique_ptr<c_Layer> layer) {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        if (!layer) {
            throw std::invalid_argument("TidalPy: cannot add a null layer to a world");
        }
        const std::string rejection = this->layer_rejection_reason(*layer);
        if (!rejection.empty()) { throw std::invalid_argument(rejection); }
        // The layer's profile reads now take turns with this world's calls, and its setters reach this world.
        layer->set_owner(this, this->p_call_mutex.get());
        this->p_layers.push_back(std::move(layer));
        // The solved structure describes the old stack, and so do its tides.
        this->p_reset_solved_state();
        this->p_warm_start_central_pressure = TidalPyConstants::d_NAN;   // a different planet now
        this->p_tidal_heating_record = c_TidalHeatingRecord();
    }

    // Why add_layer would refuse `layer`, or an empty string when it would take it: the layer's index must be its
    // place in the stack, the layer must continue the stack (its inner radius at the previous layer's outer radius),
    // stay inside the world radius, and carry a name no other layer has, since layers are reached by name.
    std::string layer_rejection_reason(const c_Layer& layer) const {
        const std::size_t position = this->p_layers.size();
        if (layer.get_layer_index() != static_cast<int>(position)) {
            return "TidalPy: layer '" + layer.get_name() + "' has layer_index " +
                   std::to_string(layer.get_layer_index()) + ", but it would be layer " + std::to_string(position) +
                   " of world '" + this->get_name() + "' (indices count from 0, inner to outer); give layer_index " +
                   std::to_string(position) + " or leave it out.";
        }
        const double prev_outer = this->p_layers.empty() ? 0.0
                                : this->p_layers.back()->get_radius_outer();
        if (std::abs(layer.get_radius_inner() - prev_outer) > layer_continuity_tol(prev_outer)) {
            return "TidalPy: layer geometry is not continuous. Layer '" + layer.get_name() + "' starts at " +
                   std::to_string(layer.get_radius_inner()) + " m, but the stack of world '" + this->get_name() +
                   "' ends at " + std::to_string(prev_outer) + " m (add layers inner-to-outer, innermost starting at "
                   "radius 0; leave radius_inner out to start at the top of the stack).";
        }
        if (layer.get_radius_outer() > this->p_radius + layer_continuity_tol(this->p_radius)) {
            return "TidalPy: layer '" + layer.get_name() + "' reaches " + std::to_string(layer.get_radius_outer()) +
                   " m, past the radius of world '" + this->get_name() + "' (" + std::to_string(this->p_radius) +
                   " m).";
        }
        for (const auto& existing : this->p_layers) {
            if (existing->get_name() == layer.get_name()) {
                return "TidalPy: world '" + this->get_name() + "' already has a layer named '" + layer.get_name() +
                       "'; give each layer its own name.";
            }
        }
        return std::string();
    }

    // Non-owning observer pointer.
    c_Layer* get_layer(std::size_t index) const {
        if (index >= this->p_layers.size()) {
            throw std::out_of_range("TidalPy: layer index out of range");
        }
        return this->p_layers[index].get();
    }

    std::size_t get_num_layers() const noexcept { return this->p_layers.size(); }

    // Whole-planet EOS profile queries, MKS. Each reads the dense output of the layer containing r, clamped at the
    // surface. NaN when no layer contains r or the EOS has not been solved.
    //
    // Every profile read holds the world's call lock (c_WorldCallLock) for the whole call, so it takes turns with
    // solve_eos, load_binary, and the other locked calls, which replace the profile it reads; get_eos_fields and
    // calc_complex_moduli take it once for a whole array of radii. A lock that cannot be taken (a broken process)
    // terminates, since these are noexcept, rather than returning a value read without it.
    double get_density(double radius) const noexcept {
        return this->p_eos_field(radius, C_EOS_DENSITY_INDEX);
    }

    // Every solved quantity at a radius in one dense evaluation, in the layout of eos_layout_.hpp, through the
    // layer the single-quantity getters use; NaN throughout for a world with no layers.
    void get_eos_state(double radius, double* y_out) const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        this->p_layer_eos_state(radius, y_out);
    }

    double get_gravity(double radius) const noexcept {
        return this->p_eos_field(radius, C_EOS_GRAVITY_INDEX);
    }

    double get_pressure(double radius) const noexcept {
        return this->p_eos_field(radius, C_EOS_PRESSURE_INDEX);
    }

    // Temperature [K] and the heat flowing outward through the sphere of that radius [W]. A solve with no
    // temperature contrast reports each layer's own temperature and no flow. Read from the world's solution
    // itself, so NaN where no layer spans the radius (no clamp at the surface).
    double get_temperature(double radius) const noexcept {
        return this->p_eos_field(radius, C_EOS_TEMPERATURE_INDEX);
    }

    double get_heat_flow(double radius) const noexcept {
        return this->p_eos_field(radius, C_EOS_HEAT_FLOW_INDEX);
    }

    // Vectorized profile read: entries field_indices[0 .. num_fields) of the evaluation layout (eos_layout_.hpp) at
    // each of radii[0 .. num_radii) [m], written field-major: values_out[field_i * num_radii + radius_i]. Each
    // entry is what its single getter returns (get_temperature and get_heat_flow for those two indices), from one
    // dense evaluation per radius and source; an index outside the layout gives NaN. The world's call lock is taken
    // once for the whole call, so every value comes from one solve and a loop over radii pays for one lock.
    //
    // Assumes values_out holds num_fields * num_radii doubles.
    void get_eos_fields(
            const std::size_t* field_indices,
            std::size_t num_fields,
            const double* radii,
            std::size_t num_radii,
            double* values_out) const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        // Temperature and heat flow come from the world's solution, the rest from the layer's dense output.
        bool needs_layer_state = false;
        bool needs_world_state = false;
        for (std::size_t field_i = 0; field_i < num_fields; ++field_i) {
            if (c_is_world_eos_field(field_indices[field_i])) { needs_world_state = true; }
            else                                              { needs_layer_state = true; }
        }
        double layer_state[C_EOS_DY_VALUES];
        double world_state[C_EOS_DY_VALUES];
        for (std::size_t radius_i = 0; radius_i < num_radii; ++radius_i) {
            const double radius = radii[radius_i];
            if (needs_layer_state) { this->p_layer_eos_state(radius, layer_state); }
            if (needs_world_state) { this->p_world_eos_state(radius, world_state); }
            for (std::size_t field_i = 0; field_i < num_fields; ++field_i) {
                const std::size_t field_index = field_indices[field_i];
                double value = TidalPyConstants::d_NAN;
                if (field_index < C_EOS_DY_VALUES) {
                    value = c_is_world_eos_field(field_index) ? world_state[field_index] : layer_state[field_index];
                }
                values_out[field_i * num_radii + radius_i] = value;
            }
        }
    }

    // Per-layer results of the last solve (the solve_eos result reports each of them per layer): temperatures,
    // heat flows, boundary layers, and the Rayleigh and Nusselt numbers of a convecting layer.
    const std::vector<c_LayerThermal>& get_layer_thermal() const noexcept { return this->p_layer_thermal; }
    size_t get_thermal_passes()    const noexcept { return this->p_thermal_passes; }
    bool   get_thermal_converged() const noexcept { return this->p_thermal_converged; }

    // Rate of change of a layer's temperature [K s-1] from the heat entering, leaving, and generated in it:
    //   (C + C_latent) dT/dt = L_in - L_out + H,
    // with L_in and L_out the heat flows of the last solve, H every heat source of a layer with use_heating, C the heat
    // its profile stores per kelvin of its temperature (calc_layer_thermal_capacity), and C_latent the latent heat
    // the boundaries between its solid and liquid zones absorb per kelvin where its material melts at one temperature
    // (a Stefan condition; calc_layer_latent_capacity). The radiogenic and prescribed heat are the last solve's;
    // the tidal heat is the latest calc_tides (the tidal heat source), so a step of an evolution is solve_eos,
    // calc_tides, then this rate, with no second solve. NaN for a layer with no heat capacity (one without a
    // material). Read under the call lock, since solve_eos replaces the thermal state.
    double calc_layer_temperature_rate(std::size_t layer_index) const {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        if (layer_index >= this->p_layer_thermal.size()) { return TidalPyConstants::d_NAN; }
        const c_LayerThermal& thermal = this->p_layer_thermal[layer_index];
        const c_Layer* layer = this->p_layers[layer_index].get();
        const double mass = layer->get_mass();
        this->p_ensure_capacities();
        const double thermal_capacity = this->p_layer_thermal_capacity[layer_index];
        const double latent_capacity  = this->p_layer_latent_capacity[layer_index];
        if (!(thermal_capacity > TidalPyConstants::d_EPS)) { return TidalPyConstants::d_NAN; }
        double heating = 0.0;
        if (this->p_solve_state) {
            const c_Heating& sources = this->p_solve_state->heating;
            heating += sources.calc_layer_power(c_HeatSourceKind::Radiogenic, layer_index, mass);
            heating += sources.calc_layer_power(c_HeatSourceKind::Prescribed, layer_index, mass);
        }
        const std::vector<double>& tidal_power = this->p_tidal_heating_record.layer_power;
        if (layer->get_use_heating() && (layer_index < tidal_power.size()) && std::isfinite(tidal_power[layer_index])) {
            heating += tidal_power[layer_index];
        }
        return (thermal.heat_flow_in - thermal.heat_flow_out + heating) / (thermal_capacity + latent_capacity);
    }

    // The heat a layer's profile stores per kelvin of its temperature [J K-1] (c_layer_thermal_capacity): the integral
    // of rho c_p over the layer, weighted by how far each point moves with the layer's temperature, with the latent heat
    // of a melting range included. Computed from the last solve on the first call after it, then cached. NaN before a
    // successful solve or for a layer index past the last layer.
    double calc_layer_thermal_capacity(std::size_t layer_index) const {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        if (layer_index >= this->p_layer_thermal.size()) { return TidalPyConstants::d_NAN; }
        this->p_ensure_capacities();
        return this->p_layer_thermal_capacity[layer_index];
    }

    // The latent heat [J K-1] the boundaries between a layer's solid and liquid zones absorb per kelvin of its
    // temperature, where its material melts at one temperature (c_zone_boundary_latent_capacity); zero for a layer
    // without such a boundary. Computed and cached as calc_layer_thermal_capacity is. NaN before a successful solve or
    // for a layer index past the last layer.
    double calc_layer_latent_capacity(std::size_t layer_index) const {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        if (layer_index >= this->p_layer_thermal.size()) { return TidalPyConstants::d_NAN; }
        this->p_ensure_capacities();
        return this->p_layer_latent_capacity[layer_index];
    }

    // Heat sources (heating_.hpp). Every source acts in the layers with use_heating, through a solve that carries
    // temperature.
    //
    // The tidal source holds the last calc_tides heating of each layer, and every later solve_eos spreads it, so an
    // evolution alternates calc_tides and solve_eos. It is kept, not invalidated, when the structure changes;
    // clear_tidal_heating forgets it. NaN entries are layers calc_tides gave no heating.
    std::vector<double> get_tidal_heat_source() const {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_tidal_heating_record.layer_power;
    }
    void clear_tidal_heating() noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        this->p_tidal_heating_record = c_TidalHeatingRecord();
    }

    // A layer's prescribed heating: a power [W] spread over the layer by mass, or a specific rate [W kg-1]; NaN for
    // both clears it. Both finite throws std::invalid_argument, as do an infinite value and a layer index out of
    // range. The EOS solve
    // reads it, so the world forgets its solved structure.
    void set_prescribed_heating(std::size_t layer_index, double power, double specific_rate) {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        if (layer_index >= this->p_layers.size()) {
            throw std::invalid_argument("TidalPy: no layer at index " + std::to_string(layer_index) + ".");
        }
        if (std::isfinite(power) && std::isfinite(specific_rate)) {
            throw std::invalid_argument(
                "TidalPy: prescribe a layer's heating as a power or as a specific rate, not both.");
        }
        if (std::isinf(power) || std::isinf(specific_rate)) {
            throw std::invalid_argument("TidalPy: a prescribed heating must be finite (NaN clears it).");
        }
        if (this->p_prescribed_heating.size() < this->p_layers.size()) {
            this->p_prescribed_heating.resize(this->p_layers.size());
        }
        this->p_prescribed_heating[layer_index].power         = std::isfinite(power) ? power : TidalPyConstants::d_NAN;
        this->p_prescribed_heating[layer_index].specific_rate =
            std::isfinite(specific_rate) ? specific_rate : TidalPyConstants::d_NAN;
        this->p_reset_solved_state();
    }
    c_PrescribedLayerHeating get_prescribed_heating(std::size_t layer_index) const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return (layer_index < this->p_prescribed_heating.size())
            ? this->p_prescribed_heating[layer_index] : c_PrescribedLayerHeating();
    }

    // The volumetric heating [W m-3] the last solve's sources give at a radius [m], every source summed (zero in a
    // layer without use_heating); NaN before a successful solve or where no layer holds the radius. A solve with
    // solve_temperature off carries no heat flow, so the heating reported here did not act on it.
    double get_heating(double radius) const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        if (!this->p_eos_solved || !this->p_solve_state || this->p_layers.empty()) { return TidalPyConstants::d_NAN; }
        std::size_t layer_index = this->p_layers.size() - 1;
        for (std::size_t layer_i = 0; layer_i < this->p_layers.size(); ++layer_i) {
            if (radius <= this->p_layers[layer_i]->get_radius_outer()) { layer_index = layer_i; break; }
        }
        if ((radius < 0.0) || (radius > this->p_layers.back()->get_radius_outer())) { return TidalPyConstants::d_NAN; }
        double state[C_EOS_DY_VALUES];
        this->p_layers[layer_index]->p_eos_state(radius, state);
        return this->p_solve_state->heating.calc_heating(layer_index, radius, state[C_EOS_DENSITY_INDEX]);
    }

    // Viscoelastic profile queries, post-melt with pre-melt variants. NaN when no layer contains r or the EOS has
    // not been solved.
    double get_shear_modulus(double radius) const noexcept {
        return this->p_eos_field(radius, C_EOS_SHEAR_MODULUS_INDEX);
    }
    double get_bulk_modulus(double radius) const noexcept {
        return this->p_eos_field(radius, C_EOS_BULK_MODULUS_INDEX);
    }
    double get_shear_viscosity(double radius) const noexcept {
        return this->p_eos_field(radius, C_EOS_SHEAR_VISCOSITY_INDEX);
    }
    double get_bulk_viscosity(double radius) const noexcept {
        return this->p_eos_field(radius, C_EOS_BULK_VISCOSITY_INDEX);
    }
    double get_melt_fraction(double radius) const noexcept {
        return this->p_eos_field(radius, C_EOS_MELT_FRACTION_INDEX);
    }

    // The only per-frequency step: find the layer and apply its rheology to the stored post-melt static
    // modulus and viscosity. The 3D tide paths call these with the call lock already held, which the recursive
    // lock allows.
    std::complex<double> calc_complex_shear_modulus(double radius, double frequency) const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_complex_modulus(true, radius, frequency);
    }
    std::complex<double> calc_complex_bulk_modulus(double radius, double frequency) const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_complex_modulus(false, radius, frequency);
    }

    // Vectorized form of the two above at one frequency [rad s-1]: the shear (is_shear) or bulk complex modulus [Pa]
    // at each of radii[0 .. num_radii) [m], taking the call lock once for the whole call.
    //
    // Assumes moduli_out holds num_radii values.
    void calc_complex_moduli(
            bool is_shear,
            const double* radii,
            std::size_t num_radii,
            double frequency,
            std::complex<double>* moduli_out) const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        for (std::size_t radius_i = 0; radius_i < num_radii; ++radius_i) {
            moduli_out[radius_i] = this->p_complex_modulus(is_shear, radii[radius_i], frequency);
        }
    }

    // True after a successful EOS solve of the current layer stack, until something invalidates it.
    bool get_eos_solved() const noexcept { return this->p_eos_solved; }

    // A precondition for the world-level EOS solve.
    bool get_all_materials_set() const noexcept {
        if (this->p_layers.empty()) { return false; }
        for (const auto& layer : this->p_layers) {
            if (!layer->get_material_set()) { return false; }
        }
        return true;
    }

    // Whole-planet equation-of-state solve. Integrates gravity, pressure, enclosed mass, and moment of
    // inertia from center to surface with each layer's material (with the layer's switches) as the local density
    // source, then populates every layer's c_LayerEOSData and mass. All MKS.
    //
    // A layer that holds its mass (is_volume_fixed off) ends where it encloses that mass, found by the integration
    // itself, and carries every layer above it: the layers above keep their volumes and the world's radius is the top
    // of the outermost. A layer whose material can change state is split into solid and liquid zones where its
    // rigidity margin changes sign (c_EOSLayerBounds); get_zones lists them, and the Love solve integrates each zone
    // as a layer of its own.
    //
    // Assumes spherical symmetry, and that each layer's density comes from its material: from the pressure (and the
    // temperature, with thermal expansion) for the analytic laws, the radius for a tabulated one.
    void solve_eos(const c_WorldEOSSolveConfig& cfg) {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        const std::size_t n_layers = this->p_layers.size();
        if (n_layers == 0) {
            throw std::invalid_argument("TidalPy: cannot solve the EOS for a world with no layers");
        }
        if (!this->get_all_materials_set()) {
            throw std::invalid_argument(
                "TidalPy: every layer must have a material (Layer.material) before the world EOS can be solved");
        }
        if (cfg.slices_per_layer < 2) {
            throw std::invalid_argument("TidalPy: slices_per_layer must be at least 2");
        }

        double G_to_use = cfg.G_to_use;
        if (G_to_use < 0.0) {
            G_to_use = c_get_G();
        }

        const std::size_t slices = cfg.slices_per_layer;

        // The bulk density guess as the volume-weighted EOS density at the surface pressure.
        double total_volume  = 0.0;
        double mass_estimate = 0.0;
        // What a layer holding its mass holds on to: its mass, from its configuration or from the first solve. A solve
        // that fails or throws puts these back as they were, so a mass taken from a rejected structure never sticks.
        const std::vector<double> reference_mass_before_solve = this->p_reference_mass;
        if (cfg.reset_layer_masses) { this->p_reference_mass.clear(); }
        this->p_reference_mass.resize(n_layers, TidalPyConstants::d_NAN);
        for (std::size_t i = 0; i < n_layers; ++i) {
            c_Layer* layer   = this->p_layers[i].get();
            const double r_inner = layer->get_radius_inner();
            const double r_outer = layer->get_radius_outer();
            if (!layer->get_is_volume_fixed() && !cfg.reset_layer_masses
                && !std::isfinite(this->p_reference_mass[i])
                && (layer->get_mass() > TidalPyConstants::d_EPS)) {
                this->p_reference_mass[i] = layer->get_mass();
            }
            c_ThermoPoint guess_point;
            guess_point.pressure    = cfg.surface_pressure;
            guess_point.temperature = TidalPyConstants::d_NAN;
            guess_point.radius      = 0.5 * (r_inner + r_outer);
            const double rho_mid    = layer->get_material()->calc_density(guess_point, layer->get_switches());
            const double shell_vol = (4.0 / 3.0) * TidalPyConstants::d_PI
                                   * (r_outer * r_outer * r_outer - r_inner * r_inner * r_inner);
            total_volume  += shell_vol;
            mass_estimate += rho_mid * shell_vol;
        }
        const double planet_bulk_density = (total_volume > TidalPyConstants::d_EPS) ? (mass_estimate / total_volume) : 3500.0;
        const double planet_radius       = this->p_layers.back()->get_radius_outer();
        // Tides, the potential, and the homogeneous Love methods use the world radius, so the layers must end there.
        if (std::abs(planet_radius - this->p_radius) > layer_continuity_tol(this->p_radius)) {
            throw std::invalid_argument(
                "TidalPy: world '" + this->get_name() + "' has radius " + std::to_string(this->p_radius) +
                " m but its outermost layer ends at " + std::to_string(planet_radius) +
                " m; the layers must fill the world exactly.");
        }

        // The default non-dimensional solve integrates in the radius, the bulk density, and 1/sqrt(pi G rho),
        // so the tolerances mean the same thing for every planet and the central pressure is of order one.
        // An SI solve keeps every scale at one.
        double length_scale  = 1.0;
        double density_scale = 1.0;
        double pascal_scale  = 1.0;
        double mass_scale    = 1.0;
        double second2_scale = 1.0;
        std::unique_ptr<c_NonDimensionalScales> scales_uptr;
        if (cfg.nondimensionalize) {
            scales_uptr   = std::make_unique<c_NonDimensionalScales>(planet_radius, planet_bulk_density);
            length_scale  = scales_uptr->length_conversion;
            density_scale = scales_uptr->density_conversion;
            pascal_scale  = scales_uptr->pascal_conversion;
            mass_scale    = scales_uptr->mass_conversion;
            second2_scale = scales_uptr->second2_conversion;
        }
        const double G_solve = G_to_use / (length_scale * length_scale * length_scale / (mass_scale * second2_scale));

        // The shear modulus [Pa] at or below which a material that can change state is a liquid: minimum_solid_rigidity
        // times rho g R of the world's stated mass and radius, fixed for the whole solve so its state events have one
        // threshold, which the Love solve and the thermal network share.
        const double liquid_shear = c_liquid_shear_threshold(this->p_mass, this->p_radius, G_to_use);

        // Per-layer pre-eval functions and inputs. The solution object is rebuilt once per thermal pass.
        std::vector<PreEvalFunc>    eos_function_vec;
        std::vector<c_EOS_ODEInput> eos_input_vec;
        eos_function_vec.reserve(n_layers);
        eos_input_vec.reserve(n_layers);

        // Everything this solve evaluates its materials with is its own: each layer's material (shared and
        // immutable, so no copy is made) with the layer's switches, the per-layer inputs the structure ODE reaches
        // it through (the solver keeps a copy of each c_EOS_ODEInput, whose eos_input_ptr points into these), and
        // the heat sources. The solution co-owns it, so a retained or exported solution keeps reading exactly what
        // was solved, whatever later happens to the layers or another solve.
        auto solve_state = std::make_shared<c_EOSSolveState>();
        solve_state->materials.reserve(n_layers);
        solve_state->inputs.assign(n_layers, c_MaterialPreevalInput());

        const double gravity_scale = length_scale / second2_scale;
        c_EOS_ODEInput ode_input;
        ode_input.G_to_use      = G_solve;
        ode_input.full_state    = false;
        ode_input.thermal_state = false;
        ode_input.length_scale  = length_scale;
        ode_input.gravity_scale = gravity_scale;
        // Whether every layer's density is set by its radius and temperature alone, which makes the surface pressure
        // of a solve at fixed temperatures fall one for one with the central pressure.
        bool density_independent_of_pressure = true;
        for (std::size_t i = 0; i < n_layers; ++i) {
            solve_state->materials.push_back(this->p_layers[i]->share_material());
            c_MaterialPreevalInput& material_input = solve_state->inputs[i];
            material_input.material_ptr  = solve_state->materials[i].get();
            material_input.switches      = this->p_layers[i]->get_switches();
            if ((material_input.material_ptr == nullptr)
                    || material_input.material_ptr->get_density_depends_on_pressure(material_input.switches)) {
                density_independent_of_pressure = false;
            }
            material_input.length_scale  = length_scale;
            material_input.pascal_scale  = pascal_scale;
            material_input.density_scale = density_scale;
            ode_input.eos_input_ptr = reinterpret_cast<char*>(&material_input);
            eos_function_vec.push_back(c_preeval_material);
            eos_input_vec.push_back(ode_input);
        }

        // Thermal layout. Each layer carries its own temperature and its cooling model says how heat moves
        // inside it; a uniform override replaces every layer's value.
        std::vector<c_LayerThermal> layer_thermal;
        c_init_layer_thermal(this->p_layers, layer_thermal, cfg.temperature);
        // A layer with no temperature of its own (the 0 K default) takes no part in the thermal network, and its
        // material sits at the cold limit of its laws, rigid where its viscosity follows the temperature. That is
        // seldom meant, so such a layer is named once, until its temperature is set.
        for (std::size_t i = 0; i < n_layers; ++i) {
            c_Layer* layer = this->p_layers[i].get();
            if (layer_thermal[i].in_network || layer->get_cold_warning_issued() || !layer->get_material_set()) {
                continue;
            }
            c_MaterialState reference_state;
            layer->calc_reference_state(reference_state);
            if (!std::isinf(reference_state.shear_viscosity)) { continue; }
            layer->set_cold_warning_issued();
            TIDALPY_LOG_WARN(
                "TidalPy: layer '{}' of world '{}' has no temperature of its own ({} K), so its material's viscosity "
                "is at its cold, rigid limit and the layer takes no part in the thermal network. Set the layer's "
                "temperature (temperature_k in a world file).",
                layer->get_name(), this->get_name(), layer->get_temperature());
        }
        // Heat sources, prepared once from the world's state; after each pass they take the layers' radii and masses
        // from it (a power is spread by mass). They need the heat flow of a thermal solve to act through, so a solve
        // with temperature switched off leaves them out.
        c_WorldState world_state;
        world_state.time           = cfg.time;
        world_state.layers_ptr     = &this->p_layers;
        world_state.tidal_ptr      = &this->p_tidal_heating_record;
        world_state.prescribed_ptr = &this->p_prescribed_heating;
        c_Heating& heating = solve_state->heating;
        heating.update_sources(world_state, length_scale, density_scale);
        const bool heating_active = heating.get_is_active();
        // An isothermal solve has no heat flow for the heating to act through. A world that pins one (the bundled
        // worlds do, to keep their fitted temperatures) while its layers carry heating for thermal solves is expected,
        // so this is a debug note rather than a warning.
        if (heating_active && !cfg.solve_temperature) {
            TIDALPY_LOG_DEBUG(
                "TidalPy: world '{}' has layers with use_heating set, but solve_temperature is off, so this EOS "
                "solve carries no heat flow and the heating is ignored.", this->get_name());
        }
        // A temperature contrast or a heated layer gives the solve a profile to integrate.
        const bool thermal_contrast = cfg.solve_temperature
            && (c_thermal_contrast_present(layer_thermal, cfg.surface_temperature) || heating_active);
        for (std::size_t i = 0; i < n_layers; ++i) {
            eos_input_vec[i].heating_ptr = thermal_contrast ? &heating : nullptr;
            eos_input_vec[i].layer_index = i;
        }

        // Pass 0 is isothermal at each layer's own temperature, which is the whole solve for a world with no
        // temperature contrast. Each later pass integrates the profile, then relaxes the boundary layers,
        // interface temperatures, and heat flows against the structure it produced.
        std::vector<c_EOSSegment> segment_vec;
        std::vector<c_EOSLayerBounds> layer_bounds;
        std::shared_ptr<c_EOSSolution> solution;
        // The central pressure the secant iteration starts from, in solve units: the world's last converged solve
        // (kept when a layer change forgets the solved state), then the pass before. A re-solve after a small change
        // then converges in a pass or two. NaN, on a world that has never been solved, starts from a uniform sphere.
        double central_pressure_guess = this->p_warm_start_central_pressure / pascal_scale;
        const std::size_t last_pass = thermal_contrast ? cfg.max_thermal_passes : 0;
        std::size_t thermal_passes = 0;
        bool thermal_converged     = !thermal_contrast;
        int structure_integrations = 0;
        // The first layer that holds its mass in the last pass, which moved itself and every layer above it.
        std::size_t first_moved = n_layers;
        try {
            for (std::size_t pass = 0; pass <= last_pass; ++pass) {
                const bool integrate_temperature = thermal_contrast && (pass > 0);
                c_build_thermal_segments(layer_thermal, integrate_temperature, length_scale, segment_vec);
                first_moved = this->p_layer_bounds(layer_thermal, length_scale, mass_scale, liquid_shear, layer_bounds);
                // The pass before this one (none on the first) places the segment tops of a layer holding its mass.
                p_set_segment_mass_fractions(solution.get(), layer_thermal, layer_bounds, length_scale, segment_vec);
                // The structure derivatives vanish above the planet radius, which a layer holding its mass leaves open.
                const double ode_planet_radius = (first_moved < n_layers)
                    ? TidalPyConstants::d_INF : layer_thermal.back().radius_outer / length_scale;
                for (std::size_t i = 0; i < n_layers; ++i) {
                    eos_input_vec[i].planet_radius = ode_planet_radius;
                    // The material always sees the temperature (its viscosities and melting do); its density sees
                    // it only when the layer uses thermal expansion.
                    solve_state->inputs[i].temperature           = layer_thermal[i].temperature;
                    solve_state->inputs[i].use_state_temperature = integrate_temperature;
                }

                c_EOSSolverSettings settings;
                settings.planet_bulk_density    = planet_bulk_density / density_scale;
                settings.surface_pressure       = cfg.surface_pressure / pascal_scale;
                settings.G_to_use               = G_solve;
                settings.integration_method     = cfg.integration_method;
                settings.rtol                   = cfg.rtol;
                settings.atol                   = cfg.atol;
                settings.pressure_tol           = cfg.pressure_tol;
                settings.max_iters              = cfg.max_iters;
                settings.verbose                = cfg.verbose;
                settings.integrate_temperature  = integrate_temperature;
                settings.central_pressure_guess = central_pressure_guess;
                settings.slices_per_layer       = slices;
                settings.length_scale           = length_scale;
                settings.pascal_scale           = pascal_scale;
                // An integrated temperature profile could carry a pressure dependence into the density.
                settings.density_independent_of_pressure = density_independent_of_pressure && !integrate_temperature;

                solution = std::make_shared<c_EOSSolution>(n_layers);
                solution->input_keepalive        = solve_state;
                solution->liquid_shear_threshold = liquid_shear;
                c_solve_eos(solution.get(), eos_function_vec, eos_input_vec, layer_bounds, segment_vec, settings);
                structure_integrations += solution->structure_integrations;

                // Return the solution to SI: the arrays, layer radii, and pressure error are scaled in place, and
                // every later evaluation of the retained integrators (call_si) converts on the way in and out.
                if (cfg.nondimensionalize) {
                    solution->dimensionalize_data(scales_uptr.get(), true);
                }
                if (!solution->success) { break; }
                thermal_passes = pass;
                central_pressure_guess = solution->central_pressure / pascal_scale;

                // The layers a layer holding its mass moved are where this pass put them, and a layer holding its
                // mass with none set takes what this pass found inside its radii.
                for (std::size_t i = 0; i < n_layers; ++i) {
                    if (i >= first_moved) {
                        if (i > first_moved) { layer_thermal[i].radius_inner = solution->get_layer_radius_inner(i); }
                        layer_thermal[i].radius_outer = solution->upper_radius_bylayer_vec[i];
                    }
                    if (!this->p_layers[i]->get_is_volume_fixed() && !std::isfinite(this->p_reference_mass[i])) {
                        double structure[C_EOS_Y_VALUES];
                        solution->call_y_si(i, layer_thermal[i].radius_outer, structure);
                        const double mass_outer = structure[C_EOS_MASS_INDEX];
                        solution->call_y_si(i, layer_thermal[i].radius_inner, structure);
                        this->p_reference_mass[i] = mass_outer - structure[C_EOS_MASS_INDEX];
                    }
                }

                if (heating_active) {
                    // A source that spreads a power by mass spreads it over the layers as this pass left them.
                    std::vector<c_HeatedLayer> heated_layers(n_layers);
                    double structure[C_EOS_Y_VALUES];
                    for (std::size_t i = 0; i < n_layers; ++i) {
                        heated_layers[i].radius_inner = layer_thermal[i].radius_inner;
                        heated_layers[i].radius_outer = layer_thermal[i].radius_outer;
                        solution->call_y_si(i, layer_thermal[i].radius_outer, structure);
                        const double mass_outer = structure[C_EOS_MASS_INDEX];
                        solution->call_y_si(i, layer_thermal[i].radius_inner, structure);
                        heated_layers[i].mass = mass_outer - structure[C_EOS_MASS_INDEX];
                    }
                    heating.update_layers(heated_layers);
                }
                if (thermal_contrast) {
                    const double thermal_change = c_update_layer_thermal(
                        *solution, this->p_layers, cfg.surface_temperature, layer_thermal, &heating);
                    thermal_converged = integrate_temperature && (thermal_change < cfg.thermal_tol);
                }
                if (thermal_converged) { break; }
            }
        } catch (...) {
            // Nothing half-solved survives a throw: the world is unsolved, its layers where they were.
            this->p_reset_solved_state();
            this->p_reference_mass = reference_mass_before_solve;
            throw;
        }

        // Commit. Whatever was solved on top of the previous structure no longer describes the world.
        this->mark_structure_dirty();
        this->p_thermal_passes     = thermal_passes;
        this->p_thermal_converged  = thermal_converged;

        // A structure far from the world's stated mass has no hydrostatic solution near it: the only surface-pressure
        // root the solve could find lies on a collapsed branch, at an absurd central pressure. That is a failed solve.
        if (solution->success) {
            const double mass_limit  = (tidalpy_config_ptr != nullptr)
                ? tidalpy_config_ptr->d_MAX_EOS_MASS_RATIO : TidalPyConstants::d_NAN;
            const double stated_mass = this->get_mass();
            if (std::isfinite(mass_limit) && (mass_limit > 1.0) && (stated_mass > 0.0)) {
                const double mass_ratio = solution->mass / stated_mass;
                if (!(mass_ratio <= mass_limit && mass_ratio >= 1.0 / mass_limit)) {
                    std::ostringstream message;
                    message << std::setprecision(4)
                        << "TidalPy: the solved structure of world '" << this->get_name() << "' holds " << mass_ratio
                        << " times its stated mass of " << stated_mass << " kg (central pressure "
                        << solution->central_pressure << " Pa), outside the factor of " << mass_limit
                        << " that [numerical] maximum_eos_mass_ratio allows. These layers have no hydrostatic "
                           "structure near the stated mass; check their materials (densities and bulk moduli) and "
                           "radii.";
                    solution->success = false;
                    solution->message = message.str();
                }
            }
        }
        // A layer held at an end of its material's pressure law: failed in tension, warned about in compression.
        if (solution->success) {
            this->p_check_pressure_law_ranges(*solution, *solve_state, slices);
        }

        // Store scalar results.
        this->p_eos_success          = solution->success;
        this->p_eos_message          = solution->message;
        this->p_eos_iterations       = solution->iterations;
        this->p_eos_structure_integrations = structure_integrations;
        this->p_eos_max_iters_hit    = solution->max_iters_hit;
        this->p_eos_pressure_error   = solution->pressure_error;
        this->p_surface_gravity_eos  = solution->surface_gravity;
        this->p_surface_pressure_eos = solution->surface_pressure;
        this->p_central_pressure     = solution->central_pressure;
        this->p_planet_mass_eos      = solution->mass;
        this->p_planet_moi_eos       = solution->moi;
        this->p_eos_solved           = solution->success && solution->other_vecs_set;

        if (!this->p_eos_solved) {
            // A failed solve leaves the world unsolved rather than holding the previous structure: no layer profile,
            // thermal state, or radial setup from before survives, and the layers stay where they were. The failed
            // solution stays for its diagnostics only.
            this->p_layer_thermal.clear();
            this->p_solve_state.reset();
            this->p_love = c_LoveWorkspace();
            for (const auto& layer_uptr : this->p_layers) { layer_uptr->clear_eos_data(); }
            this->p_reference_mass = reference_mass_before_solve;
            this->p_eos_solution = std::move(solution);
            return;
        }
        this->p_warm_start_central_pressure = this->p_central_pressure;
        this->p_solve_state   = solve_state;

        // The layers a layer holding its mass moved go where the solve put them; the world's radius is the top of the
        // outermost. A layer below them keeps its own radii exactly.
        if (first_moved < n_layers) {
            for (std::size_t i = first_moved; i < n_layers; ++i) {
                c_Layer* layer = this->p_layers[i].get();
                const double radius_inner = (i > first_moved)
                    ? solution->get_layer_radius_inner(i) : layer->get_radius_inner();
                layer->set_radii(radius_inner, solution->upper_radius_bylayer_vec[i]);
                layer_thermal[i].radius_inner = radius_inner;
                layer_thermal[i].radius_outer = solution->upper_radius_bylayer_vec[i];
            }
            this->p_radius = this->p_layers.back()->get_radius_outer();
        }
        this->p_layer_thermal = std::move(layer_thermal);
        if (!thermal_converged) {
            TIDALPY_LOG_WARN(
                "TidalPy: world '{}' EOS solve used all {} thermal passes before its interface temperatures and heat "
                "flows settled to thermal_tol = {:.1e}; the temperature profile is from the last pass.",
                this->get_name(), cfg.max_thermal_passes, cfg.thermal_tol);
        }
        for (std::size_t layer_index = 0; layer_index < n_layers; ++layer_index) {
            if (!this->p_layer_thermal[layer_index].boundary_fallback) { continue; }
            TIDALPY_LOG_WARN(
                "TidalPy: layer '{}' of world '{}' convects, but its cooling model gave no boundary-layer thickness "
                "(Rayleigh number {:.3e}); each boundary layer takes {} of the layer's thickness instead. Check the "
                "layer's viscosity at its temperature.",
                this->p_layers[layer_index]->get_name(), this->get_name(),
                this->p_layer_thermal[layer_index].rayleigh_number, d_MAX_BOUNDARY_FRACTION);
        }

        // Populate each layer's structure and viscoelastic profile from its slice of the full arrays.
        {
            const std::size_t total_slices = solution->radius_array_size;
            for (std::size_t layer_index = 0; layer_index < n_layers; ++layer_index) {
                const std::size_t slice_start = layer_index * slices;
                const std::size_t slice_end   = slice_start + slices;
                if (slice_end > total_slices) { break; }
                c_Layer* layer = this->p_layers[layer_index].get();

                // The enclosed mass gained across its slices. Adjacent layers share the interface slice, so
                // the layer masses sum exactly to the planet mass.
                layer->set_mass(
                    solution->mass_array_vec[slice_end - 1] - solution->mass_array_vec[slice_start]);

                c_LayerEOSData eos_data;
                eos_data.populate(
                    std::vector<double>(solution->radius_array_vec.begin()   + slice_start, solution->radius_array_vec.begin()   + slice_end),
                    std::vector<double>(solution->density_array_vec.begin()  + slice_start, solution->density_array_vec.begin()  + slice_end),
                    std::vector<double>(solution->gravity_array_vec.begin()  + slice_start, solution->gravity_array_vec.begin()  + slice_end),
                    std::vector<double>(solution->pressure_array_vec.begin() + slice_start, solution->pressure_array_vec.begin() + slice_end));

                // The dense-output evaluator: the solution's SI-radius call, which runs the layer's retained
                // CySolverResult and its EOS function. The captured shared_ptr co-owns the whole solution, so
                // the dense data and solver stay callable post-solve, and through it the solve state its EOS
                // arguments point at. The lambda compiles in this,
                // CyRK-owning, extension, so the CySolverResult is only ever called by the CyRK copy that
                // built it. The slice arrays above are the fallback for a manual update_eos_data.
                if (layer_index < solution->num_layers
                    && !solution->cysolver_results_uptr_vec.empty()) {
                    std::shared_ptr<c_EOSSolution> solution_owner = solution;
                    eos_data.set_dense_eval(
                        [solution_owner, layer_index](double radius, double* y_out) {
                            solution_owner->call_si(layer_index, radius, y_out);
                        });
                }
                layer->update_eos_data(eos_data);
            }
        }

        this->p_eos_solution = std::move(solution);
        this->p_update_zones();
        // The heat capacities wait for their first reader (p_ensure_capacities): their quadratures cost more than an
        // isothermal solve, and only a temperature rate or a capacity getter needs them.
        this->p_capacities_current = false;
    }

    // The solid and liquid zones of the last successful EOS solve, inner to outer: each layer is one zone unless its
    // material changed state inside it (a layer that can change state, c_Layer::get_can_change_state), where its
    // post-melt rigidity mu / (rho g R) crossed the config's minimum_solid_rigidity. The integration found each
    // boundary as a root of that margin. A zone thinner than [numerical] minimum_zone_fraction of the world's radius
    // takes the state of its thicker neighbor in the layer, since the radial solver cannot grid it. Radii [m] and
    // enclosed masses [kg]; a liquid zone is the state the solve found, whatever the layer's flags now say (the Love
    // solve reads those afresh, get_radial_zones).
    const std::vector<c_EOSZone>& get_zones() const noexcept { return this->p_zones; }

    // The zones the radial solver integrates, under the layers' current flags: a layer that can change state
    // (c_Layer::get_can_change_state, which the solve found zones for; a change to it forgets the solve) takes the
    // solved zones, and any other layer is one zone in the state its flags give it.
    std::vector<c_EOSZone> get_radial_zones() const {
        std::vector<c_EOSZone> zones;
        std::size_t zone_i = 0;
        for (std::size_t layer_i = 0; layer_i < this->p_layers.size(); ++layer_i) {
            const c_Layer* layer = this->p_layers[layer_i].get();
            std::size_t layer_zone_end = zone_i;
            while ((layer_zone_end < this->p_zones.size()) && (this->p_zones[layer_zone_end].layer_index == layer_i)) {
                ++layer_zone_end;
            }
            if (!layer->get_can_change_state() || (layer_zone_end == zone_i)) {
                c_EOSZone zone;
                zone.layer_index  = layer_i;
                zone.radius_inner = layer->get_radius_inner();
                zone.radius_outer = layer->get_radius_outer();
                zone.liquid       = layer->get_is_liquid();
                if (layer_zone_end > zone_i) {
                    zone.mass_inner = this->p_zones[zone_i].mass_inner;
                    zone.mass_outer = this->p_zones[layer_zone_end - 1].mass_outer;
                }
                zones.push_back(zone);
            } else {
                zones.insert(zones.end(), this->p_zones.begin() + static_cast<std::ptrdiff_t>(zone_i),
                             this->p_zones.begin() + static_cast<std::ptrdiff_t>(layer_zone_end));
            }
            zone_i = layer_zone_end;
        }
        return zones;
    }

    // Whether the radial solver takes a radius [m] as liquid: inside a liquid zone (get_radial_zones). A radius on a
    // zone boundary takes the zone below it.
    bool get_is_liquid_at(double radius) const {
        const std::vector<c_EOSZone> zones = this->get_radial_zones();
        for (const c_EOSZone& zone : zones) {
            if (radius <= zone.radius_outer) { return zone.liquid; }
        }
        return zones.empty() ? false : zones.back().liquid;
    }

    // The liquid zones of the layers that are not liquid throughout, the stretches melting made liquid. A copy taken
    // under the call lock, since solve_eos replaces the zones.
    std::vector<c_EOSZone> get_molten_regions() const {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        std::vector<c_EOSZone> regions;
        for (const c_EOSZone& zone : this->get_radial_zones()) {
            if (!zone.liquid) { continue; }
            if (!this->p_layers[zone.layer_index]->get_is_liquid()) { regions.push_back(zone); }
        }
        return regions;
    }

    // A copy of the solved zones, taken under the call lock.
    std::vector<c_EOSZone> get_zones_copy() const {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_zones;
    }

    // Everything solved on top of the structure describes the structure it was solved with, the cached
    // radial-solver setup included: a re-solve changes the structure and moduli even when the grid size does
    // not. solve_eos calls this before it replaces the structure, so nothing can be read against the new one.
    void mark_structure_dirty() noexcept {
        this->p_love.invalidate();
        this->p_tides_solved          = false;
        this->p_tide_solver_love.clear();
        this->p_layer_tidal_heating.clear();
        for (const auto& layer_uptr : this->p_layers) {
            layer_uptr->set_tidal_heating(TidalPyConstants::d_NAN);
        }
        this->p_zones.clear();
        this->p_capacities_current = false;
    }

    // A layer this world owns changed something its EOS solve reads: its radii (moved through a view), its material
    // or a model the material holds, its temperature or thermal switches, its cooling or radiogenics, or whether it
    // holds its volume (c_LayerOwner; the layer's setters call this). The solved structure no longer describes the
    // world, so it is forgotten with everything built on it, as before the first solve_eos: a solve that needs the
    // EOS raises until solve_eos runs again. The central pressure it converged to stays as the next solve's starting
    // guess.
    void update_after_layer_change() override {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        this->p_reset_solved_state();
    }

protected:
    // Checks each layer of a just-finished solve against the range its material's pressure law represents
    // (c_PressureLawRange; the Birch-Murnaghan and Vinet laws are the ones with limits), on the material's base phase.
    // Past either end the law holds the density at that end with a bulk modulus near zero. The law sees the cold
    // pressure: the pressure less the thermal pressure alpha0 K0 (T - T_ref) when the layer uses thermal expansion,
    // at the layer's temperature or, in a solve that integrates temperature, the solved one.
    //
    // Past the tension end the material cannot hold together at its temperature and pressure, so what the layer
    // reports there is the pinned state rather than a structure: the solve fails with a message naming the layer.
    // Past the compression end (Birch-Murnaghan with K0' below 4 turns over at a finite compression) only the
    // deepest slices are held, so the solve warns and stands.
    void p_check_pressure_law_ranges(c_EOSSolution& solution, const c_EOSSolveState& solve_state, std::size_t slices) {
        const std::size_t n_layers = std::min(this->p_layers.size(), solve_state.materials.size());
        for (std::size_t layer_index = 0; layer_index < n_layers; ++layer_index) {
            if (!solve_state.materials[layer_index] || (layer_index >= solve_state.inputs.size())) { continue; }
            const c_EOSBase& eos = solve_state.materials[layer_index]->get_base_phase().get_eos();
            const c_PressureLawRange law_range = eos.get_pressure_law_range();
            if (!std::isfinite(law_range.pressure_min) && !std::isfinite(law_range.pressure_max)) { continue; }

            // The coldest and the most compressed slice of the layer, in the pressure the law sees.
            const c_MaterialPreevalInput& input = solve_state.inputs[layer_index];
            const bool thermal                  = input.switches.use_thermal_expansion;
            const std::size_t slice_start       = layer_index * slices;
            const std::size_t slice_end = std::min(slice_start + slices, solution.pressure_array_vec.size());
            double cold_pressure_min    = TidalPyConstants::d_INF;
            double cold_pressure_max    = -TidalPyConstants::d_INF;
            double thermal_pressure_min = 0.0;   // at the slice of cold_pressure_min [Pa]
            for (std::size_t slice_i = slice_start; slice_i < slice_end; ++slice_i) {
                const double temperature =
                    (input.use_state_temperature && (slice_i < solution.temperature_array_vec.size()))
                    ? solution.temperature_array_vec[slice_i] : input.temperature;
                const double thermal_pressure = eos.calc_thermal_pressure(temperature, thermal);
                const double cold_pressure = solution.pressure_array_vec[slice_i] - thermal_pressure;
                if (cold_pressure < cold_pressure_min) {
                    cold_pressure_min    = cold_pressure;
                    thermal_pressure_min = thermal_pressure;
                }
                cold_pressure_max = std::max(cold_pressure_max, cold_pressure);
            }

            if (cold_pressure_min < law_range.pressure_min) {
                std::ostringstream message;
                message << std::setprecision(4)
                    << "TidalPy: layer '" << this->p_layers[layer_index]->get_name() << "' of world '"
                    << this->get_name() << "' is in tension past what its material's pressure law represents: the "
                       "pressure the law sees (the pressure less a thermal pressure alpha0 K0 (T - T_ref) of "
                    << thermal_pressure_min << " Pa) falls to " << cold_pressure_min << " Pa, below the law's tension "
                       "limit of " << law_range.pressure_min << " Pa, where its density is held at the law's smallest "
                       "compression with a bulk modulus near zero. Lower the layer's temperature or its material's "
                       "thermal expansivity, or check the law's reference bulk modulus and its derivative.";
                solution.success = false;
                solution.message = message.str();
                return;
            }
            if (cold_pressure_max > law_range.pressure_max) {
                TIDALPY_LOG_WARN(
                    "TidalPy: layer '{}' of world '{}' reaches {:.4e} Pa in the pressure its material's law sees, past "
                    "the {:.4e} Pa the law represents; its density is held at the law's largest compression there. "
                    "Check the law's bulk modulus derivative (a Birch-Murnaghan K0' below 4 turns over).",
                    this->p_layers[layer_index]->get_name(), this->get_name(), cold_pressure_max,
                    law_range.pressure_max);
            }
        }
    }

    // Forget everything solved: the EOS solution, the layer profiles, the thermal state, and every result built
    // on them. Called when the layers no longer match what was solved (a layer added, a binary load) and when a
    // solve throws, so no reader can mistake an old structure for the current one.
    // A solved mass more than this fraction away from the stated mass is warned about (p_warn_if_mass_differs).
    static constexpr double mass_mismatch_warn_fraction = 1.0e-2;

    // Warns, once per world, when a successful solve's mass differs from the stated mass by more than
    // mass_mismatch_warn_fraction: the two then describe different planets.
    void p_warn_if_mass_differs() {
        const double stated_mass = this->get_mass();
        if (this->p_mass_mismatch_warned || !(stated_mass > 0.0) || !std::isfinite(this->p_planet_mass_eos)) {
            return;
        }
        const double mass_difference = std::abs(this->p_planet_mass_eos / stated_mass - 1.0);
        if (!(mass_difference > mass_mismatch_warn_fraction)) { return; }
        this->p_mass_mismatch_warned = true;
        TIDALPY_LOG_WARN(
            "TidalPy: world '{}' solved to a mass of {:.6e} kg, {:.1f} percent off its stated mass of {:.6e} kg. Its "
            "Love numbers, tides, and moment of inertia follow the solved structure, while its orbit in a System, "
            "calc_surface_gravity, and calc_mean_density use the stated mass. Adjust the layers' materials or radii, "
            "or the stated mass, so the two agree. Shown once per world.",
            this->get_name(), this->p_planet_mass_eos, 100.0 * mass_difference, stated_mass);
    }
    bool p_mass_mismatch_warned = false;

    void p_reset_solved_state() {
        this->mark_structure_dirty();
        this->p_eos_solved  = false;
        this->p_eos_success = false;
        this->p_eos_message = "EOS not yet solved.";
        this->p_eos_solution.reset();
        this->p_solve_state.reset();
        this->p_layer_thermal.clear();
        this->p_love = c_LoveWorkspace();
        for (const auto& layer_uptr : this->p_layers) { layer_uptr->clear_eos_data(); }
    }

public:
    // Solve the EOS and copy its result out before any other thread can start a solve on this world. A world's own
    // solve also checks its solved mass against its stated one (p_warn_if_mass_differs); the temporary profile world
    // of the standalone radial solver, whose stated mass only sets its scales, calls solve_eos and skips the check.
    c_WorldEOSReport solve_eos_report(const c_WorldEOSSolveConfig& cfg) {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        this->solve_eos(cfg);
        if (this->p_eos_solved) { this->p_warn_if_mass_differs(); }
        return this->get_eos_report();
    }

    // A consistent copy of the last solve's result (see c_WorldEOSReport).
    c_WorldEOSReport get_eos_report() const {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        c_WorldEOSReport report;
        report.solved             = this->p_eos_solved;
        report.success            = this->p_eos_success;
        report.message            = this->p_eos_message;
        report.iterations         = this->p_eos_iterations;
        report.structure_integrations = this->p_eos_structure_integrations;
        report.max_iters_hit      = this->p_eos_max_iters_hit;
        report.pressure_error     = this->p_eos_pressure_error;
        report.surface_gravity    = this->p_surface_gravity_eos;
        report.surface_pressure   = this->p_surface_pressure_eos;
        report.central_pressure   = this->p_central_pressure;
        report.planet_mass        = this->p_planet_mass_eos;
        report.planet_moi         = this->p_planet_moi_eos;
        report.thermal_passes     = this->p_thermal_passes;
        report.thermal_converged  = this->p_thermal_converged;
        report.zones              = this->p_zones;

        const c_EOSSolution* solution = this->p_eos_solution.get();
        if (solution == nullptr || !this->p_eos_solved) { return report; }
        const std::size_t num_points = solution->radius_array_size;
        const auto copy_profile = [num_points](const std::vector<double>& source, std::vector<double>& target) {
            const std::size_t count = std::min(num_points, source.size());
            target.assign(source.begin(), source.begin() + static_cast<std::ptrdiff_t>(count));
            target.resize(num_points, TidalPyConstants::d_NAN);
        };
        copy_profile(solution->radius_array_vec,      report.radius);
        copy_profile(solution->gravity_array_vec,     report.gravity);
        copy_profile(solution->pressure_array_vec,    report.pressure);
        copy_profile(solution->mass_array_vec,        report.mass);
        copy_profile(solution->moi_array_vec,         report.moi);
        copy_profile(solution->density_array_vec,     report.density);
        copy_profile(solution->temperature_array_vec, report.temperature);
        copy_profile(solution->heat_flow_array_vec,   report.heat_flow);

        const std::size_t num_layers = std::min(this->p_layers.size(), this->p_layer_thermal.size());
        report.layer_thermal.assign(this->p_layer_thermal.begin(),
                                    this->p_layer_thermal.begin() + static_cast<std::ptrdiff_t>(num_layers));
        report.layer_radius_outer.reserve(num_layers);
        for (std::size_t layer_i = 0; layer_i < num_layers; ++layer_i) {
            report.layer_radius_outer.push_back(this->p_layers[layer_i]->get_radius_outer());
        }
        return report;
    }

    // Valid after solve_eos; NaN or empty otherwise. The message is a copy taken under the call lock.
    bool               get_eos_success()          const noexcept { return this->p_eos_success; }
    std::string get_eos_message() const {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_eos_message;
    }
    int                get_eos_iterations()       const noexcept { return this->p_eos_iterations; }
    // The central-pressure iteration reached max_iters. A solve that then still misses pressure_tol is a failure
    // (get_eos_success false) whose solution is kept for its diagnostics only.
    bool               get_eos_max_iters_hit()    const noexcept { return this->p_eos_max_iters_hit; }
    double             get_eos_pressure_error()   const noexcept { return this->p_eos_pressure_error; }
    double             get_surface_gravity_eos()  const noexcept { return this->p_surface_gravity_eos; }
    double             get_surface_pressure_eos() const noexcept { return this->p_surface_pressure_eos; }
    double             get_central_pressure()     const noexcept { return this->p_central_pressure; }
    double             get_planet_mass_eos()      const noexcept { return this->p_planet_mass_eos; }
    double             get_planet_moi_eos()       const noexcept { return this->p_planet_moi_eos; }

    // Rates only. The world drives its c_Spin model with its own EOS-based moment of inertia, so the
    // spin-rate change uses the structure-resolved value rather than the uniform-density one. The setter takes the
    // call lock, so an evolution call on another thread never reads a half-replaced model.
    void set_spin_model(const c_Spin& spin) noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        this->p_spin = spin;
    }
    const c_Spin&  get_spin_model() const noexcept { return this->p_spin; }

    // C / (M R^2) [dimensionless] from the solved moment of inertia and the solved mass; NaN before a successful
    // EOS solve (the spin model's factor, the fallback get_moment_of_inertia uses then, is a separate setting).
    double get_moment_of_inertia_factor() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        if (!this->p_eos_solved) { return TidalPyConstants::d_NAN; }
        const double radius = this->get_radius();
        const double denominator = this->p_planet_mass_eos * radius * radius;
        if (!(denominator > TidalPyConstants::d_EPS)) { return TidalPyConstants::d_NAN; }
        return this->p_planet_moi_eos / denominator;
    }

    // The EOS-solved value once the EOS has been solved, else the spin model's factor * M R^2 estimate.
    double get_moment_of_inertia() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        if (this->p_eos_solved && std::isfinite(this->p_planet_moi_eos)) {
            return this->p_planet_moi_eos;
        }
        return this->p_spin.calc_moment_of_inertia(this->get_mass(), this->get_radius());
    }

    // M_host * dU/dO / I, from the last calc_tides' dU/dO and the world's moment of inertia.
    double calc_spin_derivative(double host_mass) const {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        if (!this->get_tides_solved()) {
            throw std::runtime_error(
                "TidalPy: spin derivative needs a tidal solve first: call calc_tides()");
        }
        return this->p_spin.calc_dspin_dt(host_mass, this->get_tidal_dU_dO(), this->get_moment_of_inertia());
    }

    // The orbital mean motion.
    double calc_synchronous_spin(double orbital_frequency) const noexcept {
        return this->p_spin.calc_synchronous_spin(orbital_frequency);
    }

    // Non-owning; the source of the radial profile arrays, null before solve_eos.
    const c_EOSSolution* get_eos_solution() const noexcept { return this->p_eos_solution.get(); }

    // Whole-planet Love-number solve. Uses the world's EOS arrays and the layers' rheology models to build
    // the frequency-dependent complex moduli, then calls the shooting solver directly. The world class is
    // the complete interface; c_RadialSolutionStorage is an internal detail of the cached radial solver.
    // Spherical symmetry, all MKS.
    //
    // This first step rebuilds the workspace's cached setup when it no longer matches the current EOS and config.
    // False, with an error stamped on the solver storage, for an invalid per-layer slice partitioning.
    bool ensure_radial_cache(const c_LoveSolveConfig& cfg, c_LoveWorkspace& workspace) const {
        if (!this->p_eos_solved || !this->p_eos_solution)
            throw std::invalid_argument("TidalPy: solve_eos must succeed before solve_love_numbers");
        if (tidalpy_config_ptr == nullptr)
            throw std::runtime_error("TidalPy: config not initialized - call initialize_tidalpy_config() first");

        const std::size_t n_layers     = this->p_layers.size();
        const std::size_t total_slices = this->p_eos_solution->radius_array_size;
        if (n_layers == 0)
            throw std::invalid_argument("TidalPy: world has no layers");
        if (total_slices < 5)
            throw std::invalid_argument("TidalPy: EOS solution has too few radial slices (< 5)");

        if (!workspace.radial_solver)
            workspace.radial_solver = std::make_unique<::c_WorldRadialSolver>();
        ::c_WorldRadialSolver* solver = workspace.radial_solver.get();

        // The solver's layers are the world's zones (get_radial_zones): every layer once, except that a layer the
        // solve split into solid and liquid zones gives one solver layer per zone. A zone takes its layer's static
        // and compressibility flags, and the liquid equations where it is liquid. Gathered before the cache check
        // because the flags are user-mutable without an EOS re-solve, so a cache hit is only valid when they still
        // match.
        const c_EOSSolution* world_eos = this->p_eos_solution.get();
        const std::size_t slices = total_slices / n_layers;
        std::vector<int>         layer_types;
        std::vector<char>        is_static_flags;
        std::vector<char>        is_incompressible_flags;
        std::vector<double>      upper_radii;
        std::vector<std::size_t> world_layer_of;
        std::vector<double>      radius_grid;
        radius_grid.reserve(total_slices);
        const auto add_solver_layer = [&](std::size_t world_layer, int layer_type, bool is_static,
                                          bool is_incompressible, double radius_inner, double radius_outer) {
            layer_types.push_back(layer_type);
            is_static_flags.push_back(is_static);
            is_incompressible_flags.push_back(is_incompressible);
            upper_radii.push_back(radius_outer);
            world_layer_of.push_back(world_layer);
            const c_Layer* layer = this->p_layers[world_layer].get();
            const double layer_thickness = layer->get_radius_outer() - layer->get_radius_inner();
            const double* layer_radii = world_eos->radius_array_vec.data() + world_layer * slices;
            if ((radius_inner == layer->get_radius_inner()) && (radius_outer == layer->get_radius_outer())) {
                // A whole layer keeps the EOS slices, so a world with nothing molten solves on its own grid.
                radius_grid.insert(radius_grid.end(), layer_radii, layer_radii + slices);
                return;
            }
            // A stretch takes its share of the layer's slices, at least as many as the shooting method needs.
            const double share = (layer_thickness > 0.0) ? (radius_outer - radius_inner) / layer_thickness : 1.0;
            const std::size_t count = std::max(
                C_RS_MIN_SLICES_PER_LAYER, static_cast<std::size_t>(std::ceil(share * static_cast<double>(slices))));
            for (std::size_t slice_i = 0; slice_i < count; ++slice_i) {
                const double fraction = static_cast<double>(slice_i) / static_cast<double>(count - 1);
                radius_grid.push_back((slice_i + 1 == count)
                    ? radius_outer : radius_inner + fraction * (radius_outer - radius_inner));
            }
        };

        for (const c_EOSZone& zone : this->get_radial_zones()) {
            const c_Layer* layer = this->p_layers[zone.layer_index].get();
            add_solver_layer(zone.layer_index, zone.liquid ? 1 : 0, layer->get_is_static(),
                             layer->get_is_incompressible(), zone.radius_inner, zone.radius_outer);
        }

        const std::size_t n_solver_layers = layer_types.size();
        auto layer_types_arr = std::make_unique<int[]>(n_solver_layers);
        auto is_static_arr   = std::make_unique<bool[]>(n_solver_layers);
        auto is_incomp_arr   = std::make_unique<bool[]>(n_solver_layers);
        for (std::size_t solver_layer_i = 0; solver_layer_i < n_solver_layers; ++solver_layer_i) {
            layer_types_arr[solver_layer_i] = layer_types[solver_layer_i];
            is_static_arr[solver_layer_i]   = static_cast<bool>(is_static_flags[solver_layer_i]);
            is_incomp_arr[solver_layer_i]   = static_cast<bool>(is_incompressible_flags[solver_layer_i]);
        }
        workspace.radial_world_layer = world_layer_of;
        workspace.radial_solid_layer.assign(n_solver_layers, 0);
        for (std::size_t solver_layer_i = 0; solver_layer_i < n_solver_layers; ++solver_layer_i) {
            workspace.radial_solid_layer[solver_layer_i] = (layer_types[solver_layer_i] == 0) ? 1 : 0;
        }

        const std::size_t num_ytypes = cfg.bc_models.empty() ? 1 : cfg.bc_models.size();
        if (solver->cache_matches(
                n_solver_layers, radius_grid.size(), cfg.degree_l, cfg.nondimensionalize, num_ytypes)
            && solver->layer_flags_match(
                layer_types_arr.get(), is_static_arr.get(), is_incomp_arr.get(), upper_radii.data(), n_solver_layers))
            return true;

        const double r_planet = world_eos->radius;
        const double vol      = (4.0 / 3.0) * TidalPyConstants::d_PI * r_planet * r_planet * r_planet;
        const double bulk_rho = (vol > TidalPyConstants::d_EPS) ? this->p_planet_mass_eos / vol : 3500.0;

        return solver->build_cache(
            radius_grid,
            upper_radii,
            layer_types_arr.get(),
            is_static_arr.get(),
            is_incomp_arr.get(),
            n_solver_layers,
            r_planet,
            bulk_rho,
            cfg.degree_l,
            cfg.nondimensionalize,
            *world_eos,
            num_ytypes
        );
    }

    c_LoveSolveRuntimeConfig make_runtime_config(const c_LoveSolveConfig& cfg) const {
        c_LoveSolveRuntimeConfig rt;
        rt.frequency          = cfg.frequency;
        rt.bc_models          = cfg.bc_models;
        rt.use_prop_matrix    = (c_love_method_from_int(cfg.love_method) == c_LoveMethod::PropagationMatrix);
        rt.core_model         = cfg.core_model;
        rt.starting_method    = cfg.starting_method;
        rt.degree1_frame      = cfg.degree1_frame;
        rt.starting_radius    = cfg.starting_radius;
        rt.start_radius_tol   = cfg.start_radius_tol;
        rt.integration_method = cfg.integration_method;
        rt.rtol               = cfg.rtol;
        rt.atol               = cfg.atol;
        rt.scale_rtols        = cfg.scale_rtols;
        rt.max_num_steps      = cfg.max_num_steps;
        rt.expected_size      = cfg.expected_size;
        rt.max_ram_MB         = cfg.max_ram_MB;
        rt.max_step           = cfg.max_step;
        rt.verbose            = cfg.verbose;
        rt.love_only          = cfg.love_only;
        return rt;
    }

    // The world's own solve: its results are what the Love getters report and release_radial_storage hands out.
    void solve_love_numbers(const c_LoveSolveConfig& cfg) {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        if (cfg.warnings && c_love_method_uses_radial_solver(c_love_method_from_int(cfg.love_method))) {
            this->warn_if_dynamic_liquid_unstable(cfg.degree_l, cfg.frequency, cfg.rtol);
        }
        this->solve_love_numbers(cfg, nullptr, this->p_love);
    }

    // Reject a Love-solve configuration that no method can answer before anything is solved. A degree below 2 has
    // no tidal Love number (degree 1 is a translation of the body), and a frequency or starting-radius tolerance
    // out of range otherwise surfaces as a failed integration, or as a start deep in the planet that reports success
    // with a wrong answer.
    static void validate_love_config(const c_LoveSolveConfig& cfg) {
        // Under a tidal or free-surface condition degree 1 is a translation of the whole body, with no Love number
        // to find. A surface load of degree 1 does deform the body (its Love numbers depend on the reference frame;
        // Farrell 1972), so a solve for loading alone may ask for it. Boundary models: 0 free, 1 tidal, 2 loading.
        bool loading_only = !cfg.bc_models.empty();
        for (const int bc_model : cfg.bc_models) {
            if (bc_model != 2) { loading_only = false; }
        }
        const int min_degree = loading_only ? 1 : 2;
        if (cfg.degree_l < min_degree) {
            throw std::invalid_argument(
                "TidalPy: degree_l must be 2 or more for a tidal or free-surface Love-number solve (degree 1 is a "
                "translation of the body), and 1 or more for a loading solve; got " + std::to_string(cfg.degree_l)
                + ".");
        }
        // The homogeneous-sphere formulas (homogeneous, cpl, ctl) give the tidal Love numbers only.
        if (c_love_method_is_homogeneous(c_love_method_from_int(cfg.love_method))) {
            for (const int bc_model : cfg.bc_models) {
                if (bc_model != 1) {
                    throw std::invalid_argument(
                        std::string("TidalPy: the ") + c_love_method_name(c_love_method_from_int(cfg.love_method)) +
                        " Love-number method gives tidal Love numbers only; solve for loading or a free surface "
                        "with love_method 'radial_solver'.");
                }
            }
        }
        if (!std::isfinite(cfg.frequency) || !(cfg.frequency > 0.0)) {
            throw std::invalid_argument(
                "TidalPy: the forcing frequency of a Love-number solve must be finite and positive [rad s-1]; the "
                "Love numbers at -omega are the complex conjugates of those at omega.");
        }
        if (tidalpy_config_ptr != nullptr) {
            if ((cfg.frequency < tidalpy_config_ptr->d_MIN_FREQUENCY)
                || (cfg.frequency > tidalpy_config_ptr->d_MAX_FREQUENCY)) {
                throw std::invalid_argument(
                    "TidalPy: the forcing frequency " + std::to_string(cfg.frequency) + " rad s-1 is outside the "
                    "range set by [numerical] minimum_frequency and maximum_frequency in the TidalPy configuration "
                    "(is it in rad s-1?).");
            }
        }
        if (!(cfg.start_radius_tol > 0.0) || !(cfg.start_radius_tol < 1.0)) {
            throw std::invalid_argument(
                "TidalPy: start_radius_tol must be between 0 and 1 (exclusive); got "
                + std::to_string(cfg.start_radius_tol) + ".");
        }
        if (!(cfg.starting_radius >= 0.0) || !(cfg.max_step >= 0.0)) {
            throw std::invalid_argument("TidalPy: starting_radius and max_step must be zero (automatic) or positive.");
        }
    }

    // Relative error above which a dynamic liquid layer's Love solve is flagged (warn_if_dynamic_liquid_unstable).
    static constexpr double dynamic_liquid_warn_error = 1.0e-2;

    // Largest estimated relative error that this world's dynamic liquid layers bring to a radial Love solve at
    // (degree_l, frequency, rtol), with the layer it comes from (npos when there is none) and its growth exponent.
    // The dynamic liquid equations carry no density-gradient term, so a liquid whose density does not follow its
    // bulk modulus (d rho / dr = -rho^2 g / K; a constant-density liquid is the common case) is stratified,
    // N^2 = -g (rho g / K + (d rho / dr) / rho), and where N^2 < 0 its solutions grow as exp(E), with
    // E = sqrt(l (l + 1)) / omega * integral of sqrt(-N^2) / r dr, amplifying the integration error: the estimate is
    // rtol exp(E). An incompressible liquid has no K term in N^2. Against static-liquid solves of the bundled worlds
    // and of neutral and constant-density test bodies it flagged every breakdown measured and no neutral liquid.
    double estimate_dynamic_liquid_error(
            int degree_l,
            double frequency,
            double rtol,
            std::size_t& worst_layer,
            double& worst_growth) const {
        worst_layer  = static_cast<std::size_t>(-1);
        worst_growth = 0.0;
        if (!this->p_eos_solved || !this->p_eos_solution || !(frequency > 0.0) || (degree_l < 1)) { return 0.0; }
        const double degree_l_dbl = static_cast<double>(degree_l);
        const double sqrt_llp1    = std::sqrt(degree_l_dbl * (degree_l_dbl + 1.0));
        constexpr std::size_t num_samples = 256;
        std::vector<double> radius(num_samples), density(num_samples), gravity(num_samples), bulk(num_samples);
        double state[C_EOS_DY_VALUES];
        double worst_error = 0.0;
        // Every liquid the radial solver integrates with the dynamic equations: a liquid layer, or a liquid zone of a
        // layer that melted.
        for (const c_EOSZone& zone : this->get_radial_zones()) {
            const std::size_t layer_i = zone.layer_index;
            const c_Layer* phys = this->p_layers[layer_i].get();
            if (!zone.liquid || phys->get_is_static()) { continue; }
            const bool incompressible = phys->get_is_incompressible();
            // Short of the centre, where 1 / r is singular.
            const double radius_outer = zone.radius_outer;
            const double radius_lower = std::max(zone.radius_inner, 1.0e-3 * radius_outer);
            if (!(radius_outer > radius_lower)) { continue; }
            for (std::size_t i = 0; i < num_samples; ++i) {
                radius[i] = radius_lower
                    + (radius_outer - radius_lower) * static_cast<double>(i) / static_cast<double>(num_samples - 1);
                this->p_eos_solution->call_si(layer_i, radius[i], state, c_ZoneState::Liquid);
                density[i] = state[C_EOS_DENSITY_INDEX];
                gravity[i] = state[C_EOS_GRAVITY_INDEX];
                bulk[i]    = state[C_EOS_BULK_MODULUS_INDEX];
            }
            double integral = 0.0;
            double previous = 0.0;
            for (std::size_t i = 0; i < num_samples; ++i) {
                const std::size_t below = (i == 0) ? 0 : i - 1;
                const std::size_t above = (i + 1 == num_samples) ? i : i + 1;
                const double density_slope = (density[above] - density[below]) / (radius[above] - radius[below]);
                double minus_n2 = 0.0;
                if ((gravity[i] > 0.0) && (density[i] > 0.0)) {
                    minus_n2 = gravity[i] * density_slope / density[i];
                    if (!incompressible && (bulk[i] > 0.0)) {
                        minus_n2 += density[i] * gravity[i] * gravity[i] / bulk[i];
                    }
                }
                const double integrand = (std::isfinite(minus_n2) && (minus_n2 > 0.0))
                    ? std::sqrt(minus_n2) / radius[i] : 0.0;
                if (i > 0) { integral += 0.5 * (integrand + previous) * (radius[i] - radius[i - 1]); }
                previous = integrand;
            }
            const double growth = sqrt_llp1 * integral / frequency;
            const double error  = rtol * std::exp(std::min(growth, 700.0));
            if (error > worst_error) {
                worst_error  = error;
                worst_layer  = layer_i;
                worst_growth = growth;
            }
        }
        return worst_error;
    }

    // True the first time it is called on this world: the heating integrals log the radial nodes they leave out below
    // the radial solver's starting radius once per world, not once per call of a sweep or an evolution.
    bool claim_center_nodes_warning() const noexcept {
        if (this->p_center_nodes_warned) { return false; }
        this->p_center_nodes_warned = true;
        return true;
    }

    // Warns, once per world, when estimate_dynamic_liquid_error passes dynamic_liquid_warn_error. The callers pass
    // the solve's degree and frequency, or for a set of solves the largest degree and smallest frequency among them,
    // which bounds every solve in the set. Runs on the calling thread, before any Love-solve workers start.
    void warn_if_dynamic_liquid_unstable(int degree_l, double frequency, double rtol) const {
        if (this->p_dynamic_liquid_warned) { return; }
        std::size_t layer_i = 0;
        double growth       = 0.0;
        const double error = this->estimate_dynamic_liquid_error(degree_l, frequency, rtol, layer_i, growth);
        if (!(error > dynamic_liquid_warn_error) || (layer_i >= this->p_layers.size())) { return; }
        this->p_dynamic_liquid_warned = true;
        TIDALPY_LOG_WARN(
            "TidalPy: world '{}' solves a liquid in layer '{}' with the dynamic equations at a forcing period of "
            "{:.3g} days (degree {}), where its solutions grow by about exp({:.1f}): the estimated relative error is "
            "{:.1e} at rtol {:.0e}. A dynamic liquid whose density does not follow its bulk modulus, such as a "
            "constant-density one, grows unstable at long periods: make the layer incompressible, give it a "
            "pressure-dependent EOS, or treat it as static. Shown once per world.",
            this->get_name(), this->p_layers[layer_i]->get_name(),
            2.0 * TidalPyConstants::d_PI / frequency / 86400.0, degree_l, growth, error, rtol);
    }

    // Composite Simpson intervals per tidal layer for the quasi-homogeneous averages; even, so 129 nodes.
    static constexpr std::size_t homogeneous_quadrature_intervals = 128;

    // One tidal layer as the quasi-homogeneous Love methods see it: its volume-averaged post-melt shear and bulk
    // moduli, its log-volume-averaged post-melt viscosities, and its tidal scale. Frequency independent.
    struct c_HomogeneousLayer {
        const c_Layer* layer = nullptr;
        std::size_t layer_index  = 0;
        double shear_modulus   = TidalPyConstants::d_NAN;   // [Pa]
        double shear_viscosity = TidalPyConstants::d_NAN;   // [Pa s]
        double tidal_scale     = 0.0;                       // dimensionless
        double volume          = 0.0;                       // [m3]
    };

    // Reusable averages for a run of quasi-homogeneous Love solves against one unchanged interior: every tidal
    // layer is read once, and each frequency then costs one rheology evaluation per layer.
    struct c_HomogeneousLoveCache {
        bool built = false;
        std::vector<c_HomogeneousLayer> layers;
    };

    // A Love-number solve into `workspace`, which holds everything it produces. Nothing on the world changes, so the
    // tide paths solve into workspaces of their own without touching the world's last solve.
    void solve_love_numbers(
            const c_LoveSolveConfig& cfg,
            c_HomogeneousLoveCache* cache,
            c_LoveWorkspace& workspace) const {
        validate_love_config(cfg);
        const c_LoveMethod method = c_love_method_from_int(cfg.love_method);
        workspace.solved      = false;
        workspace.method_last = method;
        // The analytic results describe the last analytic solve only; a radial solve clears them.
        workspace.reset_analytic();
        if (c_love_method_is_homogeneous(method)) {
            this->solve_love_numbers_homogeneous(cfg, method, cache, workspace);
            return;
        }
        if (method == c_LoveMethod::LaterallyInhomogeneous) {
            throw std::logic_error(
                "TidalPy: the laterally_inhomogeneous Love-number method is reserved for the 3D Love solver and is "
                "not implemented.");
        }
        if (!this->ensure_radial_cache(cfg, workspace)) { workspace.solved = false; return; }
        ::c_WorldRadialSolver* solver = workspace.radial_solver.get();

        // Hand the solver a provider that evaluates each layer's material state at the radius the integrator
        // asks for. Nothing is sampled onto the slice grid, so the Love numbers carry no first-order slice
        // error, and the layer the provider is asked about is the one the solver is integrating: an interface
        // radius belongs to two layers, and a lookup by radius alone would give both copies the lower layer's
        // material.
        //
        // The provider stays on the storage after the solve, so the solution keeps answering at any radius:
        // the y3 of a dynamic liquid layer is rebuilt from the density and gravity it reads, and an exported
        // solution reports the complex moduli this solve used.
        solver->set_material_eval(
            this->make_material_eval(cfg.frequency, workspace.radial_world_layer, workspace.radial_solid_layer));

        c_LoveSolveRuntimeConfig rt = this->make_runtime_config(cfg);
        solver->solve(rt);
        workspace.solved = solver->get_solved();
    }

    // Solve from externally supplied complex moduli, the standalone array API path, not the layer rheology.
    void solve_love_numbers_supplied(
            const c_LoveSolveConfig& cfg,
            const std::complex<double>* shear_in,
            const std::complex<double>* bulk_in,
            const double* radius_in,
            std::size_t n_in) {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        validate_love_config(cfg);
        const c_LoveMethod method = c_love_method_from_int(cfg.love_method);
        this->p_love.solved = false;
        this->p_love.reset_analytic();
        if (!c_love_method_uses_radial_solver(method)) {
            throw std::invalid_argument(
                "TidalPy: solve_love_numbers_supplied supports only the radial_solver and propagation_matrix "
                "Love-number methods.");
        }
        this->p_love.method_last = method;
        if (!this->ensure_radial_cache(cfg, this->p_love)) { this->p_love.solved = false; return; }
        ::c_WorldRadialSolver* solver = this->p_love.radial_solver.get();

        // The run of supplied points bounding each layer: the last point at its base through the first at its
        // top. The provider then interpolates inside that run, so a repeated interface radius gives the upper
        // copy to the layer above and the lower copy to the layer below; interpolating across all layers at
        // once would give a lower layer's top the upper layer's value. The bounds match within the same
        // tolerance the cache uses, so copies differing by rounding still count as the interface.
        const double* radius_in_end = radius_in + n_in;
        std::vector<std::size_t> in_first_by_layer(this->p_layers.size(), 0);
        std::vector<std::size_t> in_count_by_layer(this->p_layers.size(), 0);
        for (std::size_t layer_i = 0; layer_i < this->p_layers.size(); ++layer_i) {
            const double r_inner   = this->p_layers[layer_i]->get_radius_inner();
            const double r_outer   = this->p_layers[layer_i]->get_radius_outer();
            const double tolerance = 1.0e-9 * std::fabs(r_outer);
            const std::size_t past_base = static_cast<std::size_t>(
                std::upper_bound(radius_in, radius_in_end, r_inner + tolerance) - radius_in);
            const std::size_t at_top = static_cast<std::size_t>(
                std::lower_bound(radius_in, radius_in_end, r_outer - tolerance) - radius_in);
            const std::size_t in_first = (past_base > 0) ? past_base - 1 : 0;
            const std::size_t in_last  = (at_top < n_in) ? at_top : n_in - 1;
            in_first_by_layer[layer_i] = in_first;
            in_count_by_layer[layer_i] = in_last - in_first + 1;
        }
        // The provider owns copies of the supplied profile rather than borrowing the caller's arrays: an
        // exported solution reads it, so it has to outlive this call. Same for the layer runs.
        std::vector<double>               radius_copy(radius_in, radius_in + n_in);
        std::vector<std::complex<double>> shear_copy(shear_in, shear_in + n_in);
        std::vector<std::complex<double>> bulk_copy(bulk_in, bulk_in + n_in);
        std::shared_ptr<const c_EOSSolution> eos_solution = this->p_eos_solution;
        const std::vector<std::size_t> world_layer_of = this->p_love.radial_world_layer;
        const std::vector<char> solid_of = this->p_love.radial_solid_layer;
        solver->set_material_eval(
            [eos_solution, radius_copy, shear_copy, bulk_copy, in_first_by_layer, in_count_by_layer, world_layer_of,
             solid_of](
                    std::size_t solver_layer_index,
                    double radius_si,
                    double* state_out,
                    std::complex<double>& shear_out,
                    std::complex<double>& bulk_out) {
                // The solver counts the zones of a split layer as layers of their own, and each reads its own side of
                // the boundary it shares with the other.
                const std::size_t layer_index = world_layer_of[solver_layer_index];
                const bool zone_is_solid = (solver_layer_index < solid_of.size()) && solid_of[solver_layer_index];
                // One dense read gives the structure and the density; the supplied profile gives the moduli.
                eos_solution->call_si(
                    layer_index, radius_si, state_out, zone_is_solid ? c_ZoneState::Solid : c_ZoneState::Liquid);
                const std::size_t in_first = in_first_by_layer[layer_index];
                const std::size_t in_count = in_count_by_layer[layer_index];
                // The supplied grid is usually uniform within a layer, so seed the search by fraction.
                const double* layer_radius_ptr = radius_copy.data() + in_first;
                const double layer_span        =
                    (in_count > 1) ? layer_radius_ptr[in_count - 1] - layer_radius_ptr[0] : 0.0;
                const std::size_t guess = (layer_span > 0.0 && radius_si > layer_radius_ptr[0])
                    ? static_cast<std::size_t>(
                        static_cast<double>(in_count - 1) * (radius_si - layer_radius_ptr[0]) / layer_span)
                    : 0;
                shear_out = c_interp_complex(
                    radius_si, layer_radius_ptr, shear_copy.data() + in_first, in_count, guess);
                bulk_out = c_interp_complex(
                    radius_si, layer_radius_ptr, bulk_copy.data() + in_first, in_count, guess);
                // Supplied moduli carry no separate unrelaxed value, so the real part stands in, and the
                // profile said nothing about viscosity: it gave the response directly.
                state_out[C_EOS_SHEAR_MODULUS_INDEX]   = shear_out.real();
                state_out[C_EOS_BULK_MODULUS_INDEX]    = bulk_out.real();
                state_out[C_EOS_SHEAR_VISCOSITY_INDEX] = TidalPyConstants::d_NAN;
                state_out[C_EOS_BULK_VISCOSITY_INDEX]  = TidalPyConstants::d_NAN;
            });

        c_LoveSolveRuntimeConfig rt = this->make_runtime_config(cfg);
        rt.redim_eos_arrays = true;
        // The returned solution's result grid is the caller's radius array, as they gave it.
        rt.sample_radius_si = radius_in;
        rt.num_sample_radii = n_in;
        solver->solve(rt);

        // The world's own EOS solution is dimensional, so the released storage reports SI scalars.
        if (solver->get_storage() != nullptr) {
            this->copy_eos_scalars_si(solver->get_storage()->get_eos_solution_ptr());
        }
        this->p_love.solved = solver->get_solved();
    }

    // The provider the radial solver reads at each integration radius: one dense EOS call for the
    // frequency-independent state, then the layer's rheology, the only part that knows the frequency. In a solid zone
    // the complex shear modulus is floored at the config's minimum_complex_rigidity times rho g R
    // (c_floor_complex_shear): a viscously relaxed solid at a near-static frequency would otherwise leave the solid
    // equations dividing by about omega eta.
    //
    // The callable co-owns the solved EOS (which keeps the layers' materials alive through its input_keepalive) and
    // the rheologies, and resolves the layers here, once, rather than at every radius. It holds nothing of this
    // world, so it stays valid after the world is gone.
    c_EOSSolution::MaterialEval make_material_eval(
            double frequency,
            const std::vector<std::size_t>& world_layer_of,
            const std::vector<char>& solid_of) const {
        const std::size_t n_layers = this->p_layers.size();
        std::vector<std::shared_ptr<const c_RheologyBase>> shear_bylayer(n_layers);
        std::vector<std::shared_ptr<const c_RheologyBase>> bulk_bylayer(n_layers);
        for (std::size_t layer_i = 0; layer_i < n_layers; ++layer_i) {
            // None for a layer with use_tides off: it responds elastically and dissipates nothing.
            shear_bylayer[layer_i] = this->p_layers[layer_i]->share_tidal_shear_rheology();
            bulk_bylayer[layer_i]  = this->p_layers[layer_i]->share_tidal_bulk_rheology();
        }
        std::shared_ptr<const c_EOSSolution> eos_solution = this->p_eos_solution;
        const double shear_floor = c_complex_shear_floor(this->p_mass, this->p_radius, c_get_G());
        return [eos_solution, shear_bylayer, bulk_bylayer, world_layer_of, solid_of, shear_floor, frequency](
                std::size_t solver_layer_index,
                double radius_si,
                double* state_out,
                std::complex<double>& shear_out,
                std::complex<double>& bulk_out) {
            // The solver counts the zones of a split layer as layers of their own, and each reads its own side of the
            // boundary it shares with the other.
            const std::size_t layer_index = world_layer_of[solver_layer_index];
            const bool zone_is_solid = (solver_layer_index < solid_of.size()) && solid_of[solver_layer_index];
            eos_solution->call_si(
                layer_index, radius_si, state_out, zone_is_solid ? c_ZoneState::Solid : c_ZoneState::Liquid);
            if (layer_index >= shear_bylayer.size()) { return; }
            const double static_shear = state_out[C_EOS_SHEAR_MODULUS_INDEX];
            const double static_bulk  = state_out[C_EOS_BULK_MODULUS_INDEX];
            // Purely real (no dissipation) where no rheology is attached.
            shear_out = shear_bylayer[layer_index]
                ? shear_bylayer[layer_index]->calc_complex_modulus(
                    static_shear, state_out[C_EOS_SHEAR_VISCOSITY_INDEX], frequency)
                : std::complex<double>(static_shear, 0.0);
            if (zone_is_solid) {
                shear_out = c_floor_complex_shear(shear_out, shear_floor);
            }
            bulk_out = bulk_bylayer[layer_index]
                ? bulk_bylayer[layer_index]->calc_complex_modulus(
                    static_bulk, state_out[C_EOS_BULK_VISCOSITY_INDEX], frequency)
                : std::complex<double>(static_bulk, 0.0);
        };
    }

    // One-shot export to a RadialSolverSolution.
    std::unique_ptr<::c_RadialSolutionStorage> release_radial_storage() {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        // Only the latest radial solve on the current structure: a later analytic solve or a structure change
        // means the stored solution no longer describes this world.
        if (!this->p_love.radial_solver || this->p_love.is_analytic()
            || !this->p_love.radial_solver->get_storage_current()) { return nullptr; }
        std::unique_ptr<::c_RadialSolutionStorage> storage = this->p_love.radial_solver->release_storage();
        // The released solution's result grid is filled on the solve's EOS grid, which a plain solve leaves empty.
        if (storage && storage->success && (storage->get_sample_radii_si().empty())) { storage->sample_onto_grid(); }
        // The cached solve keeps its EOS scalars in solve units for the next frequency; a released solution
        // reports them, so it takes the world's SI values. Nothing solves on the storage after this.
        if (storage) { this->copy_eos_scalars_si(storage->get_eos_solution_ptr()); }
        return storage;
    }

    // The planet scalars of the world's own (dimensional) EOS solution, copied onto a storage's EOS solution
    // that is about to be reported to a caller.
    void copy_eos_scalars_si(c_EOSSolution* dst) const {
        const c_EOSSolution* src = this->p_eos_solution.get();
        if (dst == nullptr || src == nullptr) { return; }
        dst->radius           = src->radius;
        dst->mass             = src->mass;
        dst->moi              = src->moi;
        dst->surface_gravity  = src->surface_gravity;
        dst->surface_pressure = src->surface_pressure;
        dst->central_pressure = src->central_pressure;
        dst->pressure_error   = src->pressure_error;
    }

    // Quasi-homogeneous Love numbers (the homogeneous, cpl, and ctl methods). Each tidal layer is treated as a
    // homogeneous incompressible planet made of its own averaged material: its volume-averaged post-melt shear
    // modulus and its log-volume-averaged post-melt viscosity (composite Simpson in radius with the r^2 weight),
    // with the planet's radius, EOS bulk density, and EOS surface gravity. The homogeneous method applies the
    // layer's rheology to that average at the forcing frequency; cpl and ctl take the static average and impose
    // the constant phase lag (1 - i/Q) or time lag (1 - i omega dt), with Q or dt from the solve config or, unset,
    // from the tides config or the attached tide model. The world's Love numbers are the sum over the tidal layers
    // of each layer's tidal scale (its volume fraction unless set) times its Love numbers, so a one-layer planet
    // gets exactly the homogeneous value and a small, weak layer cannot dominate the planet's dissipation.
    //
    // A layer takes part with the moduli its material gives, so a fluid layer adds its fluid Love number; a layer
    // that is not tidal, or one with a zero tidal scale, takes no part.
    void build_homogeneous_layers(c_HomogeneousLoveCache& cache) const {
        const std::size_t n_intervals = homogeneous_quadrature_intervals;
        const double planet_radius = this->get_radius();
        const double planet_volume =
            (4.0 / 3.0) * TidalPyConstants::d_PI * planet_radius * planet_radius * planet_radius;
        cache.layers.clear();
        for (std::size_t layer_i = 0; layer_i < this->p_layers.size(); ++layer_i) {
            const c_Layer* layer = this->p_layers[layer_i].get();
            const double tidal_scale = layer->calc_tidal_scale(planet_volume);
            if (!(tidal_scale > 0.0)) { continue; }
            const double r_inner = layer->get_radius_inner();
            const double r_outer = layer->get_radius_outer();
            if (!(r_outer > r_inner)) { continue; }

            const double dr = (r_outer - r_inner) / static_cast<double>(n_intervals);
            double weight_sum          = 0.0;
            double shear_sum           = 0.0;
            double log_shear_visc_sum  = 0.0;
            // One dense evaluation per node gives all four quantities.
            double state[C_EOS_DY_VALUES];
            for (std::size_t i = 0; i <= n_intervals; ++i) {
                const double r = (i == n_intervals) ? r_outer : r_inner + static_cast<double>(i) * dr;
                const double simpson = (i == 0 || i == n_intervals) ? 1.0 : ((i % 2 == 1) ? 4.0 : 2.0);
                const double weight = simpson * r * r;
                layer->get_eos_state(r, state);
                weight_sum         += weight;
                shear_sum          += weight * state[C_EOS_SHEAR_MODULUS_INDEX];     // post-melt
                log_shear_visc_sum += weight * std::log10(state[C_EOS_SHEAR_VISCOSITY_INDEX]);
            }
            c_HomogeneousLayer averaged;
            averaged.layer           = layer;
            averaged.layer_index     = layer_i;
            averaged.shear_modulus   = shear_sum / weight_sum;
            averaged.shear_viscosity = std::pow(10.0, log_shear_visc_sum / weight_sum);
            averaged.tidal_scale     = tidal_scale;
            averaged.volume          = layer->get_volume();
            cache.layers.push_back(averaged);
        }
        cache.built = true;
    }

    void solve_love_numbers_homogeneous(
            const c_LoveSolveConfig& cfg,
            c_LoveMethod method,
            c_HomogeneousLoveCache* cache,
            c_LoveWorkspace& workspace) const {
        workspace.solved = false;
        workspace.reset_analytic();
        if (!this->p_eos_solved || !this->p_eos_solution) {
            throw std::invalid_argument("TidalPy: solve_eos() must be called before solve_love_numbers().");
        }
        const bool use_static = (method != c_LoveMethod::Homogeneous);

        // Without a caller-owned cache the averages go into a local one that dies with this solve.
        c_HomogeneousLoveCache local_cache;
        c_HomogeneousLoveCache& active = (cache != nullptr) ? *cache : local_cache;
        if (!active.built) {
            this->build_homogeneous_layers(active);
        }
        if (active.layers.empty()) {
            workspace.analytic_error_code = -40;
            workspace.analytic_message =
                "TidalPy: the homogeneous Love-number methods need at least one tidal layer (use_tides, with a nonzero "
                "tidal_scale) with a shear modulus.";
            return;
        }

        // The bulk density of the solved structure, which the EOS surface gravity also comes from; the declared
        // mass can differ from what the layers hold.
        const double radius = this->get_radius();
        const double density_bulk =
            this->p_planet_mass_eos / ((4.0 / 3.0) * TidalPyConstants::d_PI * radius * radius * radius);
        const double gravity = this->p_surface_gravity_eos;
        if (!(gravity > 0.0) || !std::isfinite(gravity)) {
            workspace.analytic_error_code = -42;
            workspace.analytic_message =
                "TidalPy: the EOS surface gravity is not finite and positive; re-run solve_eos() before the "
                "homogeneous Love-number methods.";
            return;
        }

        // Precedence for Q and dt: the solve config, then the [tides] config, then the attached tide model.
        double fixed_q  = cfg.fixed_q;
        double fixed_dt = cfg.fixed_dt;
        if (method == c_LoveMethod::HomogeneousCPL) {
            if (!std::isfinite(fixed_q)) { fixed_q = this->get_tide_config().love_fixed_q; }
            if (!std::isfinite(fixed_q) && this->p_tide) { fixed_q = this->p_tide->get_fixed_q(cfg.degree_l); }
            if (!(fixed_q > 0.0) || !std::isfinite(fixed_q)) {
                throw std::invalid_argument(
                    "TidalPy: the cpl Love-number method needs a positive fixed_q for degree "
                    + std::to_string(cfg.degree_l) + " (pass fixed_q, set it in the tides config, or attach a tide "
                    "model that carries a fixed Q).");
            }
        } else if (method == c_LoveMethod::HomogeneousCTL) {
            if (!std::isfinite(fixed_dt)) { fixed_dt = this->get_tide_config().love_fixed_dt; }
            if (!std::isfinite(fixed_dt) && this->p_tide) { fixed_dt = this->p_tide->get_fixed_dt(cfg.degree_l); }
            if (!(fixed_dt >= 0.0) || !std::isfinite(fixed_dt)) {
                throw std::invalid_argument(
                    "TidalPy: the ctl Love-number method needs a non-negative fixed_dt for degree "
                    + std::to_string(cfg.degree_l) + " (pass fixed_dt, set love_fixed_dt_s in the tides config, or "
                    "attach a tide model that carries a fixed time lag).");
            }
        }

        c_LoveNumbers world_love;
        std::vector<c_LayerLove> layer_loves;
        layer_loves.reserve(active.layers.size());
        for (const c_HomogeneousLayer& averaged : active.layers) {
            const c_RheologyBase* rheology = averaged.layer->get_shear_rheology();
            const std::complex<double> shear = (use_static || rheology == nullptr)
                ? std::complex<double>(averaged.shear_modulus, 0.0)
                : rheology->calc_complex_modulus(averaged.shear_modulus, averaged.shear_viscosity, cfg.frequency);
            if (!std::isfinite(shear.real()) || !std::isfinite(shear.imag())) {
                workspace.analytic_error_code = -41;
                workspace.analytic_message =
                    "TidalPy: layer '" + averaged.layer->get_name() + "' has a non-finite averaged shear modulus "
                    "(viscosity or rheology not set?); the homogeneous Love-number methods need finite moduli in every "
                    "tidal layer.";
                return;
            }
            c_LoveNumbers layer_love = c_calc_homogeneous_love_numbers(
                shear,
                density_bulk,
                gravity,
                radius,
                cfg.degree_l);
            if (method == c_LoveMethod::HomogeneousCPL) {
                layer_love = c_apply_fixed_q(layer_love, fixed_q);
            } else if (method == c_LoveMethod::HomogeneousCTL) {
                layer_love = c_apply_fixed_dt(layer_love, cfg.frequency, fixed_dt);
            }
            world_love.k += averaged.tidal_scale * layer_love.k;
            world_love.h += averaged.tidal_scale * layer_love.h;
            world_love.l += averaged.tidal_scale * layer_love.l;
            c_LayerLove part;
            part.layer_index   = averaged.layer_index;
            part.tidal_scale   = averaged.tidal_scale;
            part.love          = layer_love;
            part.shear_modulus = shear;
            part.volume        = averaged.volume;
            layer_loves.push_back(part);
        }

        workspace.analytic            = world_love;
        workspace.analytic_layers     = std::move(layer_loves);
        workspace.analytic_success    = true;
        workspace.analytic_error_code = 0;
        workspace.analytic_message    =
            std::string("Quasi-homogeneous Love numbers (") + c_love_method_name(method) + ").";
        workspace.solved = true;
    }

protected:
    // The [eos_solver] and [radial_solver] settings this world's file pinned. Every solve starts from the
    // TidalPy configuration with these applied on top.
    c_EOSSolverOverrides    p_eos_solver_overrides;
    c_RadialSolverOverrides p_radial_solver_overrides;

public:
    // The settings are read by solves on other threads (calc_tides builds its Love-solve config from them), so each
    // setter, getter, and config builder below holds the call lock; the getters return copies.
    void set_eos_solver_overrides(const c_EOSSolverOverrides& overrides) noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        this->p_eos_solver_overrides = overrides;
    }
    void set_radial_solver_overrides(const c_RadialSolverOverrides& overrides) noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        this->p_radial_solver_overrides = overrides;
    }
    c_EOSSolverOverrides get_eos_solver_overrides() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_eos_solver_overrides;
    }
    c_RadialSolverOverrides get_radial_solver_overrides() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_radial_solver_overrides;
    }

    // The [eos_solver] section of the TidalPy configuration with the world's pinned keys on top.
    c_WorldEOSSolveConfig make_eos_solve_config() const {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        c_WorldEOSSolveConfig cfg;
        this->p_eos_solver_overrides.apply(cfg);
        return cfg;
    }

    // The world's configured method and its cpl / ctl parameters from the [tides] config, plus its pinned
    // [radial_solver] keys. The tide paths start from this, so the configured method drives every solve. Called on
    // the thread that holds the call lock, never from the Love-solve workers.
    c_LoveSolveConfig make_love_solve_config() const {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        c_LoveSolveConfig cfg;
        this->p_radial_solver_overrides.apply(cfg);
        const c_TideConfig& tide_cfg = this->get_tide_config();
        cfg.love_method = tide_cfg.love_method;
        cfg.fixed_q     = tide_cfg.love_fixed_q;
        cfg.fixed_dt    = tide_cfg.love_fixed_dt;
        return cfg;
    }

    // Same, for paths that need the depth-resolved radial solution: the analytic methods have no radial
    // y-functions, so they are rejected with an explanatory error.
    c_LoveSolveConfig make_radial_love_solve_config() const {
        c_LoveSolveConfig cfg = this->make_love_solve_config();
        const c_LoveMethod method = c_love_method_from_int(cfg.love_method);
        if (!c_love_method_uses_radial_solver(method)) {
            throw std::runtime_error(
                std::string("TidalPy: 3D tidal stress/strain/heating needs a depth-resolved radial solution, but the "
                            "world's Love-number method is '") + c_love_method_name(method)
                + "'. Use radial_solver or propagation_matrix (set_tide_config(love_method=...)).");
        }
        return cfg;
    }

    int  get_love_method_last_int() const noexcept { return static_cast<int>(this->p_love.method_last); }

    // Diagnostics of the last quasi-homogeneous solve, NaN after a radial-solver solve: the tidal-scale-weighted mean
    // of the layers' complex shear moduli [Pa], and the volume of the layers that took part [m3]. Like the Love
    // getters below, each holds the call lock, since the Love solves replace what they read.
    std::complex<double> get_love_analytic_shear() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        if (!this->p_love.analytic_success) { return std::complex<double>(TidalPyConstants::d_NAN, 0.0); }
        std::complex<double> weighted(0.0, 0.0);
        double scale_sum = 0.0;
        for (const c_LayerLove& part : this->p_love.analytic_layers) {
            weighted  += part.tidal_scale * part.shear_modulus;
            scale_sum += part.tidal_scale;
        }
        return (scale_sum > 0.0) ? weighted / scale_sum : std::complex<double>(TidalPyConstants::d_NAN, 0.0);
    }
    double get_love_analytic_tidal_volume() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        if (!this->p_love.analytic_success) { return TidalPyConstants::d_NAN; }
        double volume = 0.0;
        for (const c_LayerLove& part : this->p_love.analytic_layers) { volume += part.volume; }
        return volume;
    }
    // Each tidal layer's part of the last quasi-homogeneous solve; empty otherwise. A copy, so a caller can iterate
    // it while another thread's Love solve replaces the world's own.
    std::vector<c_LayerLove> get_love_layer_parts() const {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_love.analytic_layers;
    }

    // The Love results read the solver storage that solve_eos and the Love solves replace, so each read holds the
    // call lock; the message is returned as a copy for the same reason.
    bool get_love_solved() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_love.solved;
    }
    bool get_love_success() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_love.get_success();
    }
    int get_love_error_code() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_love.get_error_code();
    }
    std::string get_love_message() const {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_love.get_message();
    }
    std::size_t get_love_num_ytypes() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_love.get_num_ytypes();
    }
    // Worst-case error amplification of the surface boundary condition solve; 0 before a solve and for the
    // analytic methods. See c_estimate_surface_amplification in RadialSolver/boundaries/boundaries_.hpp.
    double get_love_surface_amplification() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_love.get_surface_amplification();
    }
    // Reciprocal condition number of the surface boundary condition system; NaN before a shooting-method solve and
    // for the other methods. See c_estimate_surface_rcond in RadialSolver/boundaries/boundaries_.hpp.
    double get_love_surface_rcond() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_love.get_surface_rcond();
    }
    // A degree-1 loading solve's frame residual: how far the surface condition its frame row replaced is from met;
    // NaN for any other solve. See c_apply_surface_bc in RadialSolver/boundaries/boundaries_.hpp.
    double get_love_surface_frame_residual() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_love.get_surface_frame_residual();
    }
    // For the given boundary-condition ytype index; the analytic methods hold a single tidal set at index 0.
    // NaN when no solve describes the current structure: never solved, failed, or followed by a solve_eos.
    std::complex<double> get_love_number_k(std::size_t ytype_idx = 0) const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_love.get_love(ytype_idx).k;
    }
    std::complex<double> get_love_number_h(std::size_t ytype_idx = 0) const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_love.get_love(ytype_idx).h;
    }
    std::complex<double> get_love_number_l(std::size_t ytype_idx = 0) const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_love.get_love(ytype_idx).l;
    }
    // The quality factor -s |k| / Im(k) and phase lag arctan2(-s Im(k), |Re(k)|) [rad] (s the sign of Re(k)) of the
    // last solve's k for a boundary-condition ytype, by the definitions of ::c_LoveNumbers (RadialSolver/love_.hpp)
    // that the standalone radial solver reports; NaN when no solve describes the structure (get_love_number_k NaN).
    double get_love_q_k(std::size_t ytype_idx = 0) const noexcept {
        const std::complex<double> love_k = this->get_love_number_k(ytype_idx);
        if (!std::isfinite(love_k.real())) { return TidalPyConstants::d_NAN; }
        return ::c_LoveNumbers(love_k, love_k, love_k).get_Q_k();
    }
    double get_love_lag_k(std::size_t ytype_idx = 0) const noexcept {
        const std::complex<double> love_k = this->get_love_number_k(ytype_idx);
        if (!std::isfinite(love_k.real())) { return TidalPyConstants::d_NAN; }
        return ::c_LoveNumbers(love_k, love_k, love_k).get_lag_k();
    }
    // Whether the latest Love solve kept only its Love numbers and surface values (love_only).
    bool get_love_only() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_love.get_love_only();
    }
    // Surface y-value (SI) for a ytype and y-index (0..5 -> y1..y6). NaN if unsolved or after an analytic
    // solve, which has no radial functions.
    std::complex<double> get_love_surface_y(
            std::size_t ytype_idx, std::size_t y_idx) const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_love.get_surface_y(ytype_idx, y_idx);
    }
    // The same at an arbitrary radius [m]. The shooting method evaluates its dense per-layer interpolants,
    // accurate anywhere including between EOS grid slices; the matrix method interpolates its grid linearly.
    // NaN if unsolved, analytic, out of range, or below the solver's starting radius.
    std::complex<double> get_radial_solution_y(
            double radius,
            std::size_t ytype_idx,
            std::size_t y_idx) const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_love.get_radial_y(radius, ytype_idx, y_idx);
    }

    // Run the global tidal solve for the supplied orbital and spin state, and share the heating among the layers.
    // The analytic models need no layers; the rheology model runs the world's Love solve at each unique tidal
    // frequency, so it needs a solved EOS. Defined out-of-line in world_tides_.hpp, which carries the heavy
    // global-potential engine.
    void calc_tides(const c_TideSolveConfig& state);

    // On-demand 3D tidal stress, strain, and heating. The tidal potential is built from the world's [tides]
    // truncation config; there is no potential-model object. The orchestration lives on c_RheologyTide,
    // out-of-line in world_tides_.hpp, which carries the kernel and potential-engine headers.
    //
    // The secular (cycle and orbit-averaged) 3D volumetric heating [W m-3] is the time-averaged power
    // density, independent of longitude and time, and its volume integral equals the 1D global heating.
    double get_3d_tidal_heating(
            const c_TideSolveConfig& state,
            double radius,
            double colatitude);

    // Batch form over paired (radii[i], colatitudes[i]) points. Same physics as the scalar form, but the
    // radial solve is amortized across the points, one per unique (l, frequency), which makes it the
    // efficient way to build a map. num_threads applies to the radial solves and the per-point evaluation.
    void get_3d_tidal_heating_array(
            const c_TideSolveConfig& state,
            const double* radii,
            const double* colatitudes,
            size_t num_points,
            double* out_heating,
            int num_threads = 0);

    // Instantaneous tidal displacements [m] on the (radius, colatitude, longitude, time) grid; see
    // c_RheologyTide::calc_3d_displacements_grid. Same preconditions as get_3d_tidal_heating.
    void get_3d_displacements_grid(
            const c_TideSolveConfig& state,
            const c_Grid3DAxes& axes,
            double* out_disp,
            int num_threads = 0);

    // Instantaneous stress [Pa] and strain on the same grid; see c_RheologyTide::calc_3d_stress_strain_grid.
    void get_3d_stress_strain_grid(
            const c_TideSolveConfig& state,
            const c_Grid3DAxes& axes,
            double* out_stress,
            double* out_strain,
            int num_threads = 0);

    // Collapsed secular 3D tidal heating: the radial profile, colatitude profile, per-layer totals, or
    // whole-planet total, per the flags in c_Heating3DCollapseConfig.
    c_Heating3DCollapsed calc_3d_tides(
            const c_TideSolveConfig& state,
            const double* radii,
            size_t num_radii,
            const double* colatitudes,
            size_t num_colatitudes,
            const double* longitudes,
            size_t num_longitudes,
            const double* times,
            size_t num_times,
            const c_Heating3DCollapseConfig& cfg);

    // The axes and output shape calc_3d_tides produces for these inputs, values and layer_totals left empty,
    // so a caller can allocate the buffers calc_3d_tides_into fills. Needs only the layer geometry.
    c_Heating3DCollapsed calc_3d_tides_layout(
            const double* radii,
            size_t num_radii,
            const double* colatitudes,
            size_t num_colatitudes,
            const double* longitudes,
            size_t num_longitudes,
            const double* times,
            size_t num_times,
            const c_Heating3DCollapseConfig& cfg);

    // calc_3d_tides written into caller buffers. The radial solves and the per-point evaluation run on up to
    // cfg.num_threads threads; the result is identical for any thread count.
    void calc_3d_tides_into(
            const c_TideSolveConfig& state,
            const double* radii,
            size_t num_radii,
            const double* colatitudes,
            size_t num_colatitudes,
            const double* longitudes,
            size_t num_longitudes,
            const double* times,
            size_t num_times,
            const c_Heating3DCollapseConfig& cfg,
            double* out_values,
            double* out_layer_totals);

    // Each layer's heating [W] from a radial-solver tide solve: the volume integral of the radial solution's
    // orbit-averaged heating density over the layer, scaled so the layers sum to `total_heating`, and, when
    // `record_out` is given, the radial profile it integrated (defined in world_tides_.hpp).
    void calc_layer_tidal_heating_radial(
            const c_TideSolveConfig& state,
            double total_heating,
            std::vector<double>& out,
            const std::vector<c_RetainedRadialSolve>* retained_solves = nullptr,
            c_TidalHeatingRecord* record_out = nullptr);

    // The per-layer heating integral leaves out the radial nodes with no radial solution (non-finite shell power),
    // which lie below the radial solver's starting radius. They are the innermost, where a shell's power falls as r^2
    // at a nearly constant heating density, so the share of the integral they would hold is estimated from the
    // innermost solved node's, and a share above the radial solver's rtol is logged once per world.
    void p_warn_of_center_nodes(
            const std::vector<double>& radii,
            const std::vector<double>& node_weights,
            const std::vector<std::size_t>& node_layer,
            const std::vector<double>& shell_power,
            double integral_sum) const {
        if (this->p_center_nodes_warned || !(std::abs(integral_sum) > 0.0)) { return; }
        std::size_t num_missing = 0;
        double missing_radius = 0.0;
        double inner_power = TidalPyConstants::d_NAN;
        double inner_radius = TidalPyConstants::d_NAN;
        for (std::size_t node_i = 0; node_i < radii.size(); ++node_i) {
            if (!std::isfinite(shell_power[node_i])) {
                ++num_missing;
                missing_radius = std::max(missing_radius, radii[node_i]);
            } else if (!std::isfinite(inner_power)) {
                inner_power  = std::abs(shell_power[node_i]);
                inner_radius = radii[node_i];
            }
        }
        if ((num_missing == 0) || !(inner_radius > 0.0)) { return; }
        double missing_estimate = 0.0;
        for (std::size_t node_i = 0; node_i < radii.size(); ++node_i) {
            if (std::isfinite(shell_power[node_i]) || !this->p_layers[node_layer[node_i]]->get_use_tides()) {
                continue;
            }
            const double ratio = radii[node_i] / inner_radius;
            missing_estimate += node_weights[node_i] * inner_power * ratio * ratio;
        }
        const double share = missing_estimate / std::abs(integral_sum);
        const c_LoveSolveConfig love_cfg = this->make_radial_love_solve_config();
        if (!love_cfg.warnings || !(share > love_cfg.rtol) || !this->claim_center_nodes_warning()) { return; }
        TIDALPY_LOG_WARN(
            "TidalPy: world '{}': {} of the {} radial nodes of the per-layer tidal heating integral, up to r = {:.4e} "
            "m, have no radial solution (they lie below the radial solver's starting radius) and are left out. By the "
            "heating at the innermost solved radius they would hold about {:.1e} of the integral. Lower the starting "
            "radius (start_radius_tolerance or starting_radius) to include them. Shown once per world.",
            this->get_name(), num_missing, radii.size(), missing_radius, share);
    }

    // The radial solves calc_tides lends its per-layer heating integral; null outside that call (see
    // c_RetainedRadialSolve). The 3D radial-group solve takes a matching one instead of solving again.
    const std::vector<c_RetainedRadialSolve>* get_retained_radial_solves() const noexcept {
        return this->p_retained_radial_solves;
    }

    // The tidal scale of each layer in the quasi-homogeneous Love methods (c_Layer::calc_tidal_scale).
    double get_layer_tidal_scale(std::size_t index) const {
        const double planet_radius = this->get_radius();
        const double planet_volume =
            (4.0 / 3.0) * TidalPyConstants::d_PI * planet_radius * planet_radius * planet_radius;
        return this->get_layer(index)->calc_tidal_scale(planet_volume);
    }

    // The heating [W] the last calc_tides put in a layer; NaN before one, or when calc_tides was not asked for it.
    double get_layer_tidal_heating(std::size_t index) const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        if (!this->p_tides_solved || index >= this->p_layer_tidal_heating.size()) {
            return TidalPyConstants::d_NAN;
        }
        return this->p_layer_tidal_heating[index];
    }

    // Total mass [kg]: the sum of the layer masses.
    double calc_total_mass() const noexcept {
        double total = 0.0;
        for (const auto& layer : this->p_layers) { total += layer->get_mass(); }
        return total;
    }

    // Radiogenic heating [W] of every layer with a radiogenics model, at a time [s].
    double calc_internal_heating(double time) const noexcept {
        double total = 0.0;
        for (const auto& layer : this->p_layers) {
            total += layer->calc_radiogenic_heating(time, layer->get_mass());
        }
        return total;
    }

    // True when every layer boundary is continuous and the innermost layer starts at radius 0.
    bool validate_layers() const noexcept {
        double prev_outer = 0.0;
        for (const auto& layer : this->p_layers) {
            const double inner = layer->get_radius_inner();
            const double tol   = layer_continuity_tol(prev_outer);
            if (std::abs(inner - prev_outer) > tol) { return false; }
            prev_outer = layer->get_radius_outer();
        }
        return true;
    }

    uint32_t get_binary_class_id() const override { return static_cast<uint32_t>(BinaryClassID::BaseWorld); }

    // What load_binary reads a file into first, so a bad file never reaches this world, its layers, or its solved
    // state (c_TidalPyBaseClass::make_binary_scratch).
    std::unique_ptr<c_TidalPyBaseClass> make_binary_scratch() const override {
        return std::make_unique<c_BaseWorld>();
    }

protected:
    void p_write_payload(std::ostream& out) const override {
        c_StructureBase::p_write_payload(out);
        write_binary_string(out, this->p_name);
        write_binary_string(out, this->p_world_type);
        const double scalars[4] = {this->p_albedo, this->p_emissivity, this->p_obliquity, this->p_spin_frequency};
        out.write(reinterpret_cast<const char*>(scalars), sizeof(scalars));
        this->write_tide_section(out);
        const auto n_layers = static_cast<uint64_t>(this->p_layers.size());
        out.write(reinterpret_cast<const char*>(&n_layers), sizeof(uint64_t));
        const double moi_factor = this->p_spin.get_config().moment_of_inertia_factor;
        out.write(reinterpret_cast<const char*>(&moi_factor), sizeof(double));
        this->write_solver_overrides(out);
        // Each layer writes its own complete record, its models included.
        for (const auto& layer : this->p_layers) { layer->write_binary(out); }
        // The prescribed heating, one (power, specific rate) pair per layer (NaN for none), as it changes the solve.
        for (std::size_t layer_i = 0; layer_i < this->p_layers.size(); ++layer_i) {
            const c_PrescribedLayerHeating prescribed = (layer_i < this->p_prescribed_heating.size())
                ? this->p_prescribed_heating[layer_i] : c_PrescribedLayerHeating();
            const double pair[2] = {prescribed.power, prescribed.specific_rate};
            out.write(reinterpret_cast<const char*>(pair), sizeof(pair));
        }
    }

    // The layers are read whole before they replace the world's, so a corrupt record leaves the stack whole.
    void p_read_payload(std::istream& in, bool force) override {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        try {
            c_StructureBase::p_read_payload(in, force);
            this->p_name       = read_binary_string(in);
            this->p_world_type = read_binary_string(in);
            double scalars[4] = {0.0, 0.0, 0.0, 0.0};
            in.read(reinterpret_cast<char*>(scalars), sizeof(scalars));
            this->p_albedo         = scalars[0];
            this->p_emissivity     = scalars[1];
            this->p_obliquity      = scalars[2];
            this->p_spin_frequency = scalars[3];
            this->read_tide_section(in, force);
            uint64_t n_layers = 0;
            in.read(reinterpret_cast<char*>(&n_layers), sizeof(uint64_t));
            c_SpinConfig spin_config;
            in.read(reinterpret_cast<char*>(&spin_config.moment_of_inertia_factor), sizeof(double));
            if (!in) {
                throw std::runtime_error("TidalPy: failed to read world binary data");
            }
            try {
                this->p_spin = c_Spin(spin_config);
            } catch (const std::invalid_argument& error) {
                throw std::runtime_error(std::string("TidalPy: corrupt world binary data: ") + error.what());
            }
            this->read_solver_overrides(in);
            check_binary_count(in, n_layers, TIDALPY_BINARY_HEADER_BYTES, "layer");
            std::vector<std::unique_ptr<c_Layer>> loaded_layers;
            loaded_layers.reserve(n_layers);
            for (uint64_t i = 0; i < n_layers; ++i) {
                loaded_layers.push_back(c_layer_from_binary(in, force));
                loaded_layers.back()->set_owner(this, this->p_call_mutex.get());
            }
            std::vector<c_PrescribedLayerHeating> loaded_prescribed(static_cast<std::size_t>(n_layers));
            for (c_PrescribedLayerHeating& prescribed : loaded_prescribed) {
                double pair[2] = {TidalPyConstants::d_NAN, TidalPyConstants::d_NAN};
                in.read(reinterpret_cast<char*>(pair), sizeof(pair));
                prescribed.power         = pair[0];
                prescribed.specific_rate = pair[1];
            }
            if (!in) {
                throw std::runtime_error("TidalPy: failed to read world binary data (prescribed heating)");
            }
            this->p_layers = std::move(loaded_layers);
            this->p_prescribed_heating = std::move(loaded_prescribed);
        } catch (...) {
            // The world fields may already hold the new record's values, so nothing solved describes this world.
            this->p_reset_solved_state();
            throw;
        }
        // Nothing solved describes the loaded layers, and the masses floating layers held and the tides belong to the
        // old ones.
        this->p_reference_mass.clear();
        this->p_tidal_heating_record = c_TidalHeatingRecord();
        this->p_reset_solved_state();
        this->p_warm_start_central_pressure = TidalPyConstants::d_NAN;
    }

    // The [eos_solver] then the [radial_solver] keys, each in its table's order (c_SolverOverrides::write).
    void write_solver_overrides(std::ostream& out) const {
        this->p_eos_solver_overrides.write(out);
        this->p_radial_solver_overrides.write(out);
    }

    // Reads into locals and commits only when the whole block was read.
    void read_solver_overrides(std::istream& in) {
        c_EOSSolverOverrides eos;
        eos.read(in);
        c_RadialSolverOverrides radial;
        radial.read(in);
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read the world solver settings binary data");
        }
        this->p_eos_solver_overrides    = eos;
        this->p_radial_solver_overrides = radial;
    }

    // The bounds of every layer for one pass of solve_eos (c_EOSLayerBounds), on the radii the pass lays its segments
    // out on (layer_thermal): a layer holding its mass, with a mass to hold, ends where it encloses that mass, and a
    // layer whose material can change state carries its rigidity margin against liquid_shear [Pa]. Returns the first
    // layer that holds its mass, or the number of layers.
    std::size_t p_layer_bounds(
            const std::vector<c_LayerThermal>& layer_thermal,
            double length_scale,
            double mass_scale,
            double liquid_shear,
            std::vector<c_EOSLayerBounds>& bounds_out) const {
        const std::size_t n_layers = this->p_layers.size();
        bounds_out.assign(n_layers, c_EOSLayerBounds());
        std::size_t first_moved = n_layers;
        for (std::size_t layer_i = 0; layer_i < n_layers; ++layer_i) {
            const c_Layer* layer = this->p_layers[layer_i].get();
            c_EOSLayerBounds& bounds = bounds_out[layer_i];
            bounds.radius_inner_si = layer_thermal[layer_i].radius_inner;
            bounds.radius_outer_si = layer_thermal[layer_i].radius_outer;
            bounds.radius_inner    = bounds.radius_inner_si / length_scale;
            bounds.radius_outer    = bounds.radius_outer_si / length_scale;
            // The volume the layer was given, which it keeps over a base a layer below has moved.
            const double given_inner = layer->get_radius_inner() / length_scale;
            const double given_outer = layer->get_radius_outer() / length_scale;
            bounds.volume_term = given_outer * given_outer * given_outer - given_inner * given_inner * given_inner;
            bounds.holds_mass  = !layer->get_is_volume_fixed() && std::isfinite(this->p_reference_mass[layer_i]);
            if (bounds.holds_mass) {
                if (!(this->p_reference_mass[layer_i] > 0.0)) {
                    throw std::invalid_argument(
                        "TidalPy: layer '" + layer->get_name() + "' of world '" + this->get_name() +
                        "' holds its mass (is_volume_fixed is off), but that mass is not positive.");
                }
                bounds.target_mass = this->p_reference_mass[layer_i] / mass_scale;
                if (first_moved == n_layers) { first_moved = layer_i; }
            }
            bounds.liquid = layer->get_is_liquid();
            if (layer->get_can_change_state()) {
                bounds.state_event           = c_preeval_rigidity_margin;
                bounds.minimum_shear_modulus = liquid_shear;
            }
        }
        return first_moved;
    }

    // A layer holding its mass ends each of its segments where it encloses the share of that mass the segment's top
    // enclosed in the last pass, so the segments move with the layer as it expands or contracts (a Lagrangian
    // coordinate). Its last segment ends at the whole mass, and so does every segment of a first pass.
    static void p_set_segment_mass_fractions(
            const c_EOSSolution* last_solution,
            const std::vector<c_LayerThermal>& layer_thermal,
            const std::vector<c_EOSLayerBounds>& layer_bounds,
            double length_scale,
            std::vector<c_EOSSegment>& segment_vec) {
        double structure[C_EOS_Y_VALUES];
        for (std::size_t segment_i = 0; segment_i < segment_vec.size(); ++segment_i) {
            c_EOSSegment& segment = segment_vec[segment_i];
            const std::size_t layer_i = segment.layer_index;
            segment.upper_mass_fraction = 1.0;
            const bool last_of_layer = (segment_i + 1 == segment_vec.size())
                || (segment_vec[segment_i + 1].layer_index != layer_i);
            if (last_of_layer || (last_solution == nullptr) || !layer_bounds[layer_i].holds_mass) { continue; }
            last_solution->call_y_si(layer_i, layer_thermal[layer_i].radius_inner, structure);
            const double mass_inner = structure[C_EOS_MASS_INDEX];
            last_solution->call_y_si(layer_i, layer_thermal[layer_i].radius_outer, structure);
            const double mass_outer = structure[C_EOS_MASS_INDEX];
            last_solution->call_y_si(layer_i, segment.upper_radius * length_scale, structure);
            const double fraction = (structure[C_EOS_MASS_INDEX] - mass_inner) / (mass_outer - mass_inner);
            if (std::isfinite(fraction)) { segment.upper_mass_fraction = std::clamp(fraction, 0.0, 1.0); }
        }
    }

    // The zones of the last successful solve (get_zones): the solution's runs of pieces in one state, each layer's
    // ends on its own radii, a zone thinner than [numerical] minimum_zone_fraction of the world's radius given the
    // state of its thicker neighbor in the layer (the radial solver cannot grid it), and neighbors left in one state
    // merged.
    void p_update_zones() {
        this->p_zones.clear();
        const c_EOSSolution* solution = this->p_eos_solution.get();
        if (!this->p_eos_solved || (solution == nullptr)) { return; }
        const std::vector<c_EOSZone> solved_zones = solution->build_zones();
        const double min_fraction = (tidalpy_config_ptr != nullptr)
            ? tidalpy_config_ptr->d_MIN_ZONE_FRACTION : TidalPyConstants::d_NAN;
        const double min_thickness = (std::isfinite(min_fraction) && (min_fraction > 0.0))
            ? min_fraction * this->p_radius : 0.0;

        std::size_t zone_i = 0;
        while (zone_i < solved_zones.size()) {
            const std::size_t layer_i = solved_zones[zone_i].layer_index;
            std::size_t layer_zone_end = zone_i;
            while ((layer_zone_end < solved_zones.size()) && (solved_zones[layer_zone_end].layer_index == layer_i)) {
                ++layer_zone_end;
            }
            const c_Layer* layer = this->p_layers[layer_i].get();
            std::vector<c_EOSZone> zones(solved_zones.begin() + static_cast<std::ptrdiff_t>(zone_i),
                                         solved_zones.begin() + static_cast<std::ptrdiff_t>(layer_zone_end));
            zones.front().radius_inner = layer->get_radius_inner();
            zones.back().radius_outer  = layer->get_radius_outer();

            const auto thickness_of = [&zones](std::size_t i) { return zones[i].radius_outer - zones[i].radius_inner; };
            for (std::size_t thin_i = 0; (zones.size() > 1) && (thin_i < zones.size()); ++thin_i) {
                if (thickness_of(thin_i) >= min_thickness) { continue; }
                const std::size_t neighbor =
                    (thin_i == 0) ? 1
                    : (thin_i + 1 == zones.size()) ? thin_i - 1
                    : (thickness_of(thin_i - 1) >= thickness_of(thin_i + 1)) ? thin_i - 1
                    : thin_i + 1;
                if (zones[neighbor].liquid == zones[thin_i].liquid) { continue; }
                TIDALPY_LOG_INFO(
                    "TidalPy: layer '{}' of world '{}' has a {} zone only {:.3e} m thick between {:.6e} and {:.6e} m, "
                    "under [numerical] minimum_zone_fraction of the world's radius; the radial solver treats it as {}, "
                    "like its neighbor.",
                    layer->get_name(), this->get_name(), zones[thin_i].liquid ? "liquid" : "solid",
                    thickness_of(thin_i), zones[thin_i].radius_inner, zones[thin_i].radius_outer,
                    zones[neighbor].liquid ? "liquid" : "solid");
                zones[thin_i].liquid = zones[neighbor].liquid;
            }

            // An absorbed zone can leave two neighbors in the same state; they are one zone.
            for (const c_EOSZone& zone : zones) {
                const bool continues = (this->p_zones.size() > 0) && (this->p_zones.back().layer_index == layer_i)
                    && (this->p_zones.back().liquid == zone.liquid);
                if (continues) {
                    this->p_zones.back().radius_outer = zone.radius_outer;
                    this->p_zones.back().mass_outer   = zone.mass_outer;
                } else {
                    this->p_zones.push_back(zone);
                }
            }
            for (std::size_t kept_i = this->p_zones.size(); kept_i-- > 0;) {
                const c_EOSZone& zone = this->p_zones[kept_i];
                if (zone.layer_index != layer_i) { break; }
                if (zone.liquid && !layer->get_is_liquid()) {
                    TIDALPY_LOG_INFO(
                        "TidalPy: layer '{}' of world '{}' is liquid between {:.6e} and {:.6e} m; the radial solver "
                        "treats that zone as a liquid.",
                        layer->get_name(), this->get_name(), zone.radius_inner, zone.radius_outer);
                }
            }
            zone_i = layer_zone_end;
        }
    }

    // Fills the per-layer heat capacities of the last solve (calc_layer_thermal_capacity, calc_layer_latent_capacity)
    // on the first read after it. Every change to the solved state (mark_structure_dirty, a new solve) marks them stale.
    // The caller holds the call lock, so the cache is filled once even when several threads read it.
    void p_ensure_capacities() const {
        if (this->p_capacities_current) { return; }
        this->p_update_latent_capacity();
        this->p_update_thermal_capacity();
        this->p_capacities_current = true;
    }

    // The latent heat each layer's zone boundaries absorb per kelvin of its temperature (c_zone_boundary_latent_
    // capacity), from the zones and thermal state of the last solve. The central difference across a boundary spans
    // C_LATENT_DIFFERENCE_FRACTION of the world's radius, or a quarter of the thinner zone beside it, so both of its
    // points lie inside the two zones.
    void p_update_latent_capacity() const {
        this->p_layer_latent_capacity.assign(this->p_layer_thermal.size(), 0.0);
        const c_EOSSolution* solution = this->p_eos_solution.get();
        if ((solution == nullptr) || !this->p_solve_state) { return; }
        for (std::size_t zone_i = 1; zone_i < this->p_zones.size(); ++zone_i) {
            const c_EOSZone& below = this->p_zones[zone_i - 1];
            const c_EOSZone& above = this->p_zones[zone_i];
            const std::size_t layer_i = above.layer_index;
            if ((below.layer_index != layer_i) || (below.liquid == above.liquid)) { continue; }
            if ((layer_i >= this->p_layer_thermal.size()) || (layer_i >= this->p_solve_state->materials.size())) {
                continue;
            }
            c_ZoneBoundary boundary;
            boundary.layer_index = layer_i;
            boundary.radius      = above.radius_inner;
            boundary.half_width  = std::min(
                C_LATENT_DIFFERENCE_FRACTION * this->p_radius,
                0.25 * std::min(below.radius_outer - below.radius_inner, above.radius_outer - above.radius_inner));
            boundary.solid_below = !below.liquid;
            this->p_layer_latent_capacity[layer_i] += c_zone_boundary_latent_capacity(
                *solution,
                *this->p_solve_state->materials[layer_i],
                this->p_solve_state->inputs[layer_i].switches,
                this->p_layer_thermal,
                boundary);
        }
    }

    // The heat each layer's profile stores per kelvin of its temperature (c_layer_thermal_capacity), from the profile
    // of the last solve.
    void p_update_thermal_capacity() const {
        this->p_layer_thermal_capacity.assign(this->p_layer_thermal.size(), TidalPyConstants::d_NAN);
        const c_EOSSolution* solution = this->p_eos_solution.get();
        const std::size_t num_layers = std::min(this->p_layers.size(), this->p_layer_thermal.size());
        for (std::size_t layer_i = 0; layer_i < num_layers; ++layer_i) {
            this->p_layer_thermal_capacity[layer_i] = (solution == nullptr) ? TidalPyConstants::d_NAN
                : c_layer_thermal_capacity(*solution, *this->p_layers[layer_i], this->p_layer_thermal, layer_i);
        }
    }

    // One entry of the evaluation layout at a radius, as get_eos_fields reports it, under the call lock.
    double p_eos_field(double radius, std::size_t field_index) const noexcept {
        double value = TidalPyConstants::d_NAN;
        this->get_eos_fields(&field_index, 1, &radius, 1, &value);
        return value;
    }

    // The evaluation layout at a radius from the dense output of the layer that holds it (clamped at the surface);
    // NaN throughout for a world with no layers. The caller holds the call lock, which is also the layer's, so the
    // layer is read through its unlocked helper and an array read takes the lock once.
    void p_layer_eos_state(double radius, double* state_out) const noexcept {
        const c_Layer* layer = this->find_layer_for_radius(radius);
        if (layer == nullptr) {
            for (std::size_t value_i = 0; value_i < C_EOS_DY_VALUES; ++value_i) {
                state_out[value_i] = TidalPyConstants::d_NAN;
            }
            return;
        }
        layer->p_eos_state(radius, state_out);
    }

    // The evaluation layout at a radius from the world's own solution, through the layer whose span holds it; NaN
    // throughout before a successful solve or where no layer spans the radius. The caller holds the call lock.
    void p_world_eos_state(double radius, double* state_out) const noexcept {
        for (std::size_t value_i = 0; value_i < C_EOS_DY_VALUES; ++value_i) {
            state_out[value_i] = TidalPyConstants::d_NAN;
        }
        if (!this->p_eos_solved || !this->p_eos_solution) { return; }
        const std::size_t n_solved = std::min(this->p_layers.size(), this->p_eos_solution->num_layers);
        for (std::size_t layer_i = 0; layer_i < n_solved; ++layer_i) {
            const c_Layer* layer = this->p_layers[layer_i].get();
            if (radius >= layer->get_radius_inner() && radius <= layer->get_radius_outer()) {
                this->p_eos_solution->call_si(layer_i, radius, state_out);
                return;
            }
        }
    }

    // The shear (is_shear) or bulk complex modulus [Pa] at a radius [m] and frequency [rad s-1] from the rheology
    // of the layer that holds the radius; NaN only for a world with no layers. The caller holds the call lock, which
    // is also the layer's.
    std::complex<double> p_complex_modulus(bool is_shear, double radius, double frequency) const noexcept {
        const c_Layer* layer = this->find_layer_for_radius(radius);
        if (layer == nullptr) { return std::complex<double>(TidalPyConstants::d_NAN, 0.0); }
        return layer->p_complex_modulus(is_shear, radius, frequency);
    }

    // The layer whose radial span contains radius [m], non-owning. Radii beyond the surface clamp to the
    // outermost layer, radii below the innermost inner radius fall in the innermost one, and the result is
    // null only for a world with no layers. Public because it is a const observer query like
    // get_density(radius); the 3D tidal-heating path reads the layer's solid and incompressible flags.
public:
    const c_Layer* find_layer_for_radius(double radius) const noexcept {
        if (this->p_layers.empty()) { return nullptr; }
        for (const auto& layer : this->p_layers) {
            if (radius <= layer->get_radius_outer()) { return layer.get(); }
        }
        return this->p_layers.back().get();
    }

protected:
    // The tide section: the [tides] configuration, then the tide model behind a presence flag, so a loaded world
    // dissipates as the saved one did.
    void write_tide_section(std::ostream& out) const {
        const c_TideConfig& cfg = this->p_tide_config;
        const int32_t ints[5] = {
            cfg.min_degree_l, cfg.max_degree_l, cfg.eccentricity_truncation, cfg.obliquity_truncation,
            cfg.love_method};
        out.write(reinterpret_cast<const char*>(ints), sizeof(ints));
        const uint8_t layer_heating_byte = cfg.layer_tidal_heating ? 1 : 0;
        out.write(reinterpret_cast<const char*>(&layer_heating_byte), sizeof(uint8_t));
        out.write(reinterpret_cast<const char*>(&cfg.love_fixed_q),  sizeof(double));
        out.write(reinterpret_cast<const char*>(&cfg.love_fixed_dt), sizeof(double));
        out.write(reinterpret_cast<const char*>(&cfg.eccentricity_exact_tolerance), sizeof(double));
        write_optional_binary(out, this->p_tide);
    }

    // Reads into locals and commits only once the whole section is read and valid, so a corrupt record throws
    // std::runtime_error without replacing the world's tide model or configuration.
    void read_tide_section(std::istream& in, bool force) {
        int32_t ints[5] = {0, 0, 0, 0, 0};
        in.read(reinterpret_cast<char*>(ints), sizeof(ints));
        uint8_t layer_heating_byte = 1;
        in.read(reinterpret_cast<char*>(&layer_heating_byte), sizeof(uint8_t));
        c_TideConfig cfg;
        cfg.min_degree_l            = ints[0];
        cfg.max_degree_l            = ints[1];
        cfg.eccentricity_truncation = ints[2];
        cfg.obliquity_truncation    = ints[3];
        cfg.love_method             = ints[4];
        cfg.layer_tidal_heating     = (layer_heating_byte != 0);
        in.read(reinterpret_cast<char*>(&cfg.love_fixed_q),  sizeof(double));
        in.read(reinterpret_cast<char*>(&cfg.love_fixed_dt), sizeof(double));
        in.read(reinterpret_cast<char*>(&cfg.eccentricity_exact_tolerance), sizeof(double));
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read world tide configuration binary data");
        }
        try {
            validate_tide_config(cfg);
        } catch (const std::invalid_argument& error) {
            throw std::runtime_error(std::string("TidalPy: corrupt world binary data: ") + error.what());
        }
        std::shared_ptr<const c_TideBase> tide = read_optional_binary<c_TideBase>(in, force, c_tide_from_binary);

        this->p_tide_config  = cfg;
        this->p_tide         = std::move(tide);
        this->p_tides_solved = false;
        this->p_tide_result  = c_GlobalTideResult();
        this->p_tide_solver_love.clear();
    }

    std::string p_name;
    const c_TideStateProvider* p_tide_state_provider_ptr = nullptr;
    std::size_t p_tide_state_index = 0;
    std::string p_world_type = "world";
    double      p_albedo     = 0.3;       // [dimensionless]
    double      p_emissivity = 1.0;       // [dimensionless]
    double      p_obliquity  = 0.0;       // [rad]
    double      p_spin_frequency = 0.0;   // [rad/s]

    // The global tidal solve of calc_tides (world_tides_.hpp). Fills `outcome` as it goes, so a failed global
    // potential leaves its error code in outcome.tide_result when it throws.
    void p_solve_tides(const c_TideSolveConfig& state, c_TideSolveOutcome& outcome);

    // Checks the orbital state a tidal solve is about to use: throws std::invalid_argument for an eccentricity
    // outside [0, 1), a semi-major axis that is not positive, or an orbital frequency that is not finite, and warns
    // once per world (and truncation level) when the truncations misstate the tides (c_warn_tide_truncations), and
    // once per world for a nearly synchronous spin.
    void p_check_tide_state(const c_TideSolveConfig& state) const {
        if (!((state.eccentricity >= 0.0) && (state.eccentricity < 1.0))) {
            throw std::invalid_argument(
                "TidalPy: world '" + this->get_name() + "' tides need an eccentricity in [0, 1); got " +
                std::to_string(state.eccentricity) + ".");
        }
        if (!(state.semi_major_axis > 0.0)) {
            throw std::invalid_argument(
                "TidalPy: world '" + this->get_name() + "' tides need a positive semi-major axis; got " +
                std::to_string(state.semi_major_axis) + " m.");
        }
        if (!std::isfinite(state.orbital_frequency)) {
            throw std::invalid_argument(
                "TidalPy: world '" + this->get_name() + "' tides need a finite orbital frequency; got " +
                std::to_string(state.orbital_frequency) + " rad s-1.");
        }
        c_warn_tide_truncations(
            "world '" + this->get_name() + "'", "its [tides] table or set_tide_config", state.eccentricity,
            state.obliquity, this->p_tide_config.eccentricity_truncation, this->p_tide_config.obliquity_truncation,
            this->p_tide_config.max_degree_l, state.spin_frequency / state.orbital_frequency,
            this->p_truncation_warnings_shown);
        this->p_warn_if_near_synchronous(state);
    }
    mutable c_TruncationWarningsShown p_truncation_warnings_shown;

    // A spin within this fraction of the mean motion, but not equal to it, is taken for a synchronous rotation the
    // caller did not quite set (p_warn_if_near_synchronous).
    static constexpr double near_synchronous_warn_fraction = 1.0e-3;

    // Warns, once per world, when the spin is within near_synchronous_warn_fraction of the orbital mean motion but
    // not equal to it: the slow forcing term at (spin - n) then adds heating that can differ from the synchronous
    // value by orders of magnitude, which is rarely what was meant.
    void p_warn_if_near_synchronous(const c_TideSolveConfig& state) const {
        if (this->p_near_synchronous_warned || !(std::abs(state.orbital_frequency) > TidalPyConstants::d_EPS)) {
            return;
        }
        const double spin_offset = std::abs(state.spin_frequency / state.orbital_frequency - 1.0);
        if (!((spin_offset > 0.0) && (spin_offset < near_synchronous_warn_fraction))) { return; }
        this->p_near_synchronous_warned = true;
        TIDALPY_LOG_WARN(
            "TidalPy: world '{}' spins at {:.6e} rad s-1, a fraction {:.2e} away from its orbital mean motion of "
            "{:.6e} rad s-1. A slightly non-synchronous spin adds a slow forcing term that can change the tidal "
            "heating by orders of magnitude. A synchronous world should spin at exactly its orbital frequency "
            "(System.set_synchronous_rotation, or spin_frequency = orbital_frequency). Shown once per world.",
            this->get_name(), state.spin_frequency, spin_offset, state.orbital_frequency);
    }
    mutable bool p_near_synchronous_warned = false;

    // Global (1D) tidal dissipation state. The configuration and model are serialized (the tide section); the
    // results are not (recompute with calc_tides).
    c_TideConfig                         p_tide_config;
    std::shared_ptr<const c_TideBase>    p_tide;
    c_GlobalTideResult                   p_tide_result;
    bool                                 p_tides_solved = false;
    // Per-mode radial-solver Love numbers (k, h, l) keyed by the tidal mode (l, m, p, q),
    // retained from the most recent rheology calc_tides (empty for the analytic models).
    c_IntMap<c_Key4, c_LoveNumbers>      p_tide_solver_love;

    // Taken by the calls c_WorldCallLock lists; held through a pointer so the world stays movable.
    std::unique_ptr<std::recursive_mutex> p_call_mutex = std::make_unique<std::recursive_mutex>();

    std::vector<std::unique_ptr<c_Layer>> p_layers;

    // Set by solve_eos and not serialized, so a load repopulates them by re-solving.
    bool        p_eos_success          = false;
    bool        p_eos_solved           = false;
    std::string p_eos_message          = "EOS not yet solved.";
    int         p_eos_iterations       = -1;
    int         p_eos_structure_integrations = 0;
    bool        p_eos_max_iters_hit    = false;
    double      p_eos_pressure_error   = std::numeric_limits<double>::quiet_NaN();
    double      p_surface_gravity_eos  = std::numeric_limits<double>::quiet_NaN();
    double      p_surface_pressure_eos = std::numeric_limits<double>::quiet_NaN();
    double      p_central_pressure     = std::numeric_limits<double>::quiet_NaN();
    // The last successful solve's central pressure [Pa]: the next solve's starting guess, which outlives a reset.
    double      p_warm_start_central_pressure = std::numeric_limits<double>::quiet_NaN();
    double      p_planet_mass_eos      = std::numeric_limits<double>::quiet_NaN();
    double      p_planet_moi_eos       = std::numeric_limits<double>::quiet_NaN();
    c_Spin      p_spin {};              // spin-dynamics model (uses the world's EOS moment of inertia)
    std::shared_ptr<c_EOSSolution> p_eos_solution;  // retained full-planet solution (co-owned by layer dense evaluators)
    // The materials, inputs, and heat sources of the last successful solve, co-owned by its solution.
    std::shared_ptr<c_EOSSolveState> p_solve_state;
    // Thermal description of every layer from the last successful solve, and how the thermal passes ended.
    std::vector<c_LayerThermal> p_layer_thermal;
    // Each layer's heat capacities from the last successful solve [J K-1], filled on first read (p_ensure_capacities)
    // under the call lock; p_capacities_current is false until then and after anything replaces the solved state.
    mutable std::vector<double> p_layer_thermal_capacity;
    mutable std::vector<double> p_layer_latent_capacity;
    mutable bool                p_capacities_current = false;
    // The mass each layer holding its mass holds on to [kg]; NaN for a layer that holds its volume instead.
    std::vector<double> p_reference_mass;
    size_t p_thermal_passes    = 0;
    bool   p_thermal_converged = true;

    // The world's own Love solve (solve_love_numbers): the cached radial solver, rebuilt only when the EOS grid or
    // the layer assumptions change, or the quasi-homogeneous results. Not serialized. The call lock that guards it
    // is c_BaseWorld::p_call_mutex.
    c_LoveWorkspace p_love;
    // Set by warn_if_dynamic_liquid_unstable, which only the calling thread runs.
    mutable bool p_dynamic_liquid_warned = false;
    // Set by claim_center_nodes_warning, under the call lock.
    mutable bool p_center_nodes_warned = false;
    // Set only inside calc_layer_tidal_heating_radial, under the call lock (get_retained_radial_solves).
    const std::vector<c_RetainedRadialSolve>* p_retained_radial_solves = nullptr;
    // The solid and liquid zones of the last successful EOS solve (get_zones).
    std::vector<c_EOSZone> p_zones;

    // Each layer's share of the last calc_tides heating (not serialized).
    std::vector<double>                  p_layer_tidal_heating;
    // The tidal heat source: the last calc_tides heating, which every later solve_eos spreads over the heated
    // layers; it outlives the solves (an evolution alternates the two) until cleared, a layer is added, or a binary
    // is loaded. Not serialized.
    c_TidalHeatingRecord p_tidal_heating_record;
    // The caller's per-layer power or specific rate, by layer index. Not serialized.
    std::vector<c_PrescribedLayerHeating> p_prescribed_heating;
};

} // namespace tidalpy
