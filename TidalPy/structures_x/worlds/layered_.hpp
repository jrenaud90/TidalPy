#pragma once
/* A world built from an ordered stack of layers.
 *
 * Owns its layers inner to outer, index 0 innermost, and provides the whole-planet aggregates and geometry
 * validation. The whole-planet EOS and radial (Love number) solves, which walk every layer, are methods here.
 *
 * Binary payload: the c_BaseWorld fields and tide section (tide configuration and model), the layer count, the
 * spin model's moment-of-inertia factor, and the pinned [eos_solver] and [radial_solver] settings (each key a
 * presence flag, then its value when set), followed by each layer's own complete binary record, in index order, as
 * separate appended records. Solved state (the EOS profile, Love numbers, tide results) is not saved.
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

#include "base_.hpp"
#include "../layers/factory_.hpp"

#include "constants_.hpp"
#include "../../dynamics_x/spin_.hpp"
#include "../../Material_x/eos/solver_.hpp"
#include "material_.hpp"
#include "thermal_layout_.hpp"

// The RadialSolver sub-modules compile into this translation unit so the shooting CyRK integration runs in
// the extension that owns the CySolverResult objects, with no cross-extension call().
#include "../../Utilities_x/math_x/numerics_.hpp"
#include "../../Utilities_x/dimensions/nondimensional_.hpp"
#include "../../Utilities_x/arrays/interp_.hpp"
#include "../../RadialSolver_x/rs_constants_.hpp"
#include "../../RadialSolver_x/rs_solution_.hpp"
#include "../../RadialSolver_x/love_.hpp"
#include "../../RadialSolver_x/shooting_.hpp"
#include "../../RadialSolver_x/world_radial_solver_.hpp"

// The tide model hierarchy and its config and result structs are light, with no global-potential tables.
// The calc_tides orchestration that drives the heavy engine and collapse is defined out-of-line in
// world_tides_.hpp, so the eccentricity and obliquity tables compile into that one extension alone.
#include "../../Tides_x/classes/tide_base_.hpp"
#include "../../Tides_x/classes/tide_result_.hpp"
// Only the 3D grid declarations below need c_Grid3DAxes, so an extension that includes this header without
// world_tides_.hpp still compiles. tide_.hpp adds no includes this header does not already have and
// forward-declares c_LayeredWorld, so the dependency stays one-directional.
#include "../../Tides_x/classes/tide_.hpp"          // c_Grid3DAxes
#include "../../Tides_x/love/love_method_.hpp"
// Relative paths, not bare names, so every extension that includes layered_.hpp resolves these without
// needing Utilities_x/lookups on its include path.
#include "../../Utilities_x/lookups/keys_.hpp"
#include "../../Utilities_x/lookups/intmap_.hpp"

namespace tidalpy {

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

// True once update_constants_x has filled the solver sections of the shared config.
inline bool c_solver_config_loaded() noexcept {
    return (tidalpy_config_ptr != nullptr) && (tidalpy_config_ptr->d_EOS_SOLVER_METHOD >= 0);
}

// Grouped so c_LayeredWorld::solve_eos takes a single argument.
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
    // relative change in the interface temperatures and flows that ends them.
    size_t    max_thermal_passes  = 12;
    double    thermal_tol         = 1.0e-8;
    // Relative change in a floating layer's radius that ends the passes, and whether this solve redefines
    // the mass every floating layer holds; otherwise the world's first solve sets it.
    double    radius_tol          = 1.0e-8;
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
    }
};

// Parameters for the whole-planet Love-number solve.
struct c_LoveSolveConfig {
    double    frequency = 1.0e-5;            // [rad/s]; tidal forcing frequency
    int       degree_l  = 2;                  // harmonic degree
    // In order; 1 = tidal, 2 = loading, 0 = free.
    std::vector<int> bc_models = {1};
    int       love_method = 0;                  // c_LoveMethod as int: 0 radial_solver, 1 propagation_matrix,
                                                       // 2 homogeneous, 3 cpl, 4 ctl, 5 laterally_inhomogeneous
    double    fixed_q            = TidalPyConstants::d_NAN;   // cpl quality factor (NaN: from the tide model)
    double    fixed_dt           = TidalPyConstants::d_NAN;   // ctl time lag [s] (NaN: from the tide model)
    int       core_model         = 0;                  // propagation-matrix core starting condition (0-4)
    bool      use_kamata         = false;
    bool      nondimensionalize  = true;
    double    starting_radius    = 0.0;                // [m]; 0 -> auto
    double    start_radius_tol   = 1.0e-5;
    ODEMethod integration_method = ODEMethod::DOP853;
    double    rtol               = 1.0e-6;
    double    atol               = 1.0e-10;
    bool      scale_rtols        = false;
    size_t    max_num_steps      = 500000;
    size_t    expected_size      = 128;
    size_t    max_ram_MB         = 500;
    double    max_step           = 0.0;
    bool      verbose            = false;
    bool      warnings           = true;

    c_LoveSolveConfig() {
        if (!c_solver_config_loaded()) { return; }
        const TidalPyConfig& config = *tidalpy_config_ptr;
        this->use_kamata         = config.d_RADIAL_SOLVER_USE_KAMATA;
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

// The [eos_solver] keys a world file may pin so the file reproduces its run on its own. A key left unset
// falls through to the TidalPy configuration at every solve, so a configuration changed after the world was
// built still reaches it; a call's own arguments win over both.
struct c_EOSSolverOverrides {
    std::optional<ODEMethod> integration_method;
    std::optional<double>    rtol;
    std::optional<double>    atol;
    std::optional<double>    pressure_tol;
    std::optional<size_t>    max_iters;
    std::optional<size_t>    slices_per_layer;
    std::optional<bool>      nondimensionalize;
    std::optional<bool>      solve_temperature;

    void apply(c_WorldEOSSolveConfig& cfg) const noexcept {
        if (this->integration_method) { cfg.integration_method = *this->integration_method; }
        if (this->rtol)               { cfg.rtol               = *this->rtol; }
        if (this->atol)               { cfg.atol               = *this->atol; }
        if (this->pressure_tol)       { cfg.pressure_tol       = *this->pressure_tol; }
        if (this->max_iters)          { cfg.max_iters          = *this->max_iters; }
        if (this->slices_per_layer)   { cfg.slices_per_layer   = *this->slices_per_layer; }
        if (this->nondimensionalize)  { cfg.nondimensionalize  = *this->nondimensionalize; }
        if (this->solve_temperature)  { cfg.solve_temperature  = *this->solve_temperature; }
    }
};

// The [radial_solver] keys a world file may pin; the same rules as c_EOSSolverOverrides.
struct c_RadialSolverOverrides {
    std::optional<ODEMethod> integration_method;
    std::optional<double>    rtol;
    std::optional<double>    atol;
    std::optional<bool>      use_kamata;
    std::optional<double>    start_radius_tol;
    std::optional<bool>      scale_rtols;
    std::optional<size_t>    max_num_steps;
    std::optional<size_t>    expected_size;
    std::optional<size_t>    max_ram_MB;
    std::optional<bool>      nondimensionalize;

    void apply(c_LoveSolveConfig& cfg) const noexcept {
        if (this->integration_method) { cfg.integration_method = *this->integration_method; }
        if (this->rtol)               { cfg.rtol               = *this->rtol; }
        if (this->atol)               { cfg.atol               = *this->atol; }
        if (this->use_kamata)         { cfg.use_kamata         = *this->use_kamata; }
        if (this->start_radius_tol)   { cfg.start_radius_tol   = *this->start_radius_tol; }
        if (this->scale_rtols)        { cfg.scale_rtols        = *this->scale_rtols; }
        if (this->max_num_steps)      { cfg.max_num_steps      = *this->max_num_steps; }
        if (this->expected_size)      { cfg.expected_size      = *this->expected_size; }
        if (this->max_ram_MB)         { cfg.max_ram_MB         = *this->max_ram_MB; }
        if (this->nondimensionalize)  { cfg.nondimensionalize  = *this->nondimensionalize; }
    }
};

// c_WorldCallLock, the world's call lock, is defined in layers/base_.hpp: a world's layers take it too.

// A copy of everything a solve_eos result reports, taken under the world's call lock so that another thread's
// solve_eos cannot replace the solution while it is read. The profile arrays are empty when no solve has populated
// the world; the per-layer entries are filled only after a completed solve (get_eos_solved).
struct c_WorldEOSReport {
    bool        solved             = false;
    bool        success            = false;
    std::string message;
    int         iterations         = 0;
    bool        max_iters_hit      = false;
    double      pressure_error     = TidalPyConstants::d_NAN;  // [Pa]
    double      surface_gravity    = TidalPyConstants::d_NAN;  // [m s-2]
    double      surface_pressure   = TidalPyConstants::d_NAN;  // [Pa]
    double      central_pressure   = TidalPyConstants::d_NAN;  // [Pa]
    double      planet_mass        = TidalPyConstants::d_NAN;  // [kg]
    double      planet_moi         = TidalPyConstants::d_NAN;  // [kg m2]
    std::size_t thermal_passes     = 0;
    bool        thermal_converged  = false;
    bool        geometry_converged = false;

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
    std::vector<double>         layer_temperature_rate;  // [K s-1]
    std::vector<double>         layer_radius_outer;      // [m]
};

// What one EOS solve evaluates its materials with: a copy of every layer's material model, the per-layer inputs
// the structure ODE reaches them through, and the heat sources of a thermal solve. The solution co-owns it, so a
// retained or exported solution keeps answering exactly as solved.
struct c_EOSSolveState {
    std::vector<std::unique_ptr<c_MaterialEOSBase>> materials;
    std::vector<c_MaterialEOSInput>                 inputs;
    c_Heating                                       heating;
};

// A radial solve calc_tides has already run for one (degree, |omega|) group of its tidal modes. It lends these to the
// per-layer heating integral, which needs the same solves, so they are not run twice. Non-owning.
struct c_RetainedRadialSolve {
    int                              degree_l  = 0;
    double                           frequency = 0.0;       // |omega| [rad s-1]
    const ::c_RadialSolutionStorage* storage   = nullptr;
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

    // The radial methods: the cached solver, and the world layer each of its layers belongs to.
    std::unique_ptr<::c_WorldRadialSolver> radial_solver;
    std::vector<std::size_t>               radial_world_layer;

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
    int get_error_code() const noexcept {
        if (this->is_analytic()) { return this->analytic_error_code; }
        const auto* storage = this->get_storage();
        return storage ? storage->error_code : -100;
    }
    const std::string& get_message() const noexcept {
        static const std::string no_message = "No love-number solve has been run.";
        if (this->is_analytic()) { return this->analytic_message; }
        const auto* storage = this->get_storage();
        return storage ? storage->message : no_message;
    }
    std::size_t get_num_ytypes() const noexcept {
        if (this->is_analytic()) { return this->analytic_success ? 1 : 0; }
        const auto* storage = this->get_storage();
        return storage ? storage->num_ytypes : 0;
    }
    double get_surface_amplification() const noexcept {
        if (this->is_analytic()) { return 0.0; }
        const auto* storage = this->get_storage();
        return storage ? storage->surface_amplification : 0.0;
    }
    double get_surface_rcond() const noexcept {
        if (this->is_analytic()) { return TidalPyConstants::d_NAN; }
        const auto* storage = this->get_storage();
        return storage ? storage->surface_rcond : TidalPyConstants::d_NAN;
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

// One stretch of a layer as the radial solver integrates it. A layer is one stretch unless its partial-melt model
// makes part of it molten; it is then split at the edges of each molten stretch
// (c_LayeredWorld::update_radial_segments).
struct c_RadialSegment {
    std::size_t world_layer  = 0;
    double      radius_inner = 0.0;     // [m]
    double      radius_outer = 0.0;     // [m]
    bool        molten       = false;   // at the liquid_shear floor or below the minimum solid rigidity
};

// Round-off allowance on the liquid-shear floor. A molten point's post-melt shear modulus is the partial-melt
// model's liquid_shear returned through the solve's unit conversions, while a solid lies orders of magnitude above.
inline constexpr double d_MOLTEN_SHEAR_RTOL = 1.0e-9;

// The evaluation-layout entries (eos_layout_.hpp) a layered world reads from its own EOS solution rather than from
// the dense output of the layer that holds the radius.
inline bool c_is_world_eos_field(std::size_t field_index) noexcept {
    return (field_index == C_EOS_TEMPERATURE_INDEX) || (field_index == C_EOS_HEAT_FLOW_INDEX);
}

class c_LayeredWorld : public c_BaseWorld {
public:
    // Absolute gap allowed between a layer's inner radius and the previous layer's outer radius.
    static double layer_continuity_tol(double previous_outer_radius) noexcept {
        const double scale = (previous_outer_radius > 1.0) ? previous_outer_radius : 1.0;
        return tidalpy_config_ptr->d_LAYER_CONTINUITY_RTOL * scale;
    }

    c_LayeredWorld() = default;

    explicit c_LayeredWorld(const c_WorldConfig& cfg) : c_BaseWorld(cfg) {}

    ~c_LayeredWorld() override = default;

    // Add a layer, inner to outer; its inner radius must match the current outermost radius, 0 for the
    // first layer, within the continuity tolerance.
    void add_layer(std::unique_ptr<c_BaseLayer> layer) {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        if (!layer) {
            throw std::invalid_argument("TidalPy: cannot add a null layer to a world");
        }
        const std::string rejection = this->layer_rejection_reason(*layer);
        if (!rejection.empty()) { throw std::invalid_argument(rejection); }
        // The layer's profile reads now take turns with this world's calls.
        layer->set_owner_call_mutex(this->p_call_mutex.get());
        this->p_layers.push_back(std::move(layer));
        // The solved structure describes the old stack.
        this->p_reset_solved_state();
        this->p_warm_start_central_pressure = TidalPyConstants::d_NAN;   // a different planet now
    }

    // Why add_layer would refuse `layer`, or an empty string when it would take it: the layer must continue the
    // stack (its inner radius at the previous layer's outer radius), stay inside the world radius, and carry a
    // name no other layer has, since layers are reached by name.
    std::string layer_rejection_reason(const c_BaseLayer& layer) const {
        const double prev_outer = this->p_layers.empty() ? 0.0
                                : this->p_layers.back()->get_radius_outer();
        if (std::abs(layer.get_radius_inner() - prev_outer) > layer_continuity_tol(prev_outer)) {
            return "TidalPy: layer geometry is not continuous. Inner radius does not match the previous layer's "
                   "outer radius (add layers inner-to-outer, innermost starting at radius 0).";
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
    c_BaseLayer* get_layer(std::size_t index) const {
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
    bool   get_geometry_converged() const noexcept { return this->p_geometry_converged; }
    bool   get_thermal_converged() const noexcept { return this->p_thermal_converged; }

    // Rate of change of a layer's temperature [K s-1] from the heat entering, leaving, and generated in it:
    //   M c_p dT/dt = L_in - L_out + H.
    // NaN for a layer with no heat capacity (one that is not a solid-liquid layer).
    double calc_layer_temperature_rate(std::size_t layer_index) const noexcept {
        if (layer_index >= this->p_layer_thermal.size()) { return TidalPyConstants::d_NAN; }
        const c_LayerThermal& thermal = this->p_layer_thermal[layer_index];
        const double mass = this->p_layers[layer_index]->get_mass();
        const double heat_capacity = thermal.heat_capacity;
        if (!(mass > TidalPyConstants::d_EPS) || !(heat_capacity > TidalPyConstants::d_EPS)) {
            return TidalPyConstants::d_NAN;
        }
        return (thermal.heat_flow_in - thermal.heat_flow_out + thermal.heating) / (mass * heat_capacity);
    }

    // Viscoelastic profile queries, post-melt with pre-melt variants. NaN when no layer contains r, the
    // layer is geometry-only, or the EOS has not been solved.
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
    bool get_all_eos_set() const noexcept {
        if (this->p_layers.empty()) { return false; }
        for (const auto& layer : this->p_layers) {
            if (!layer->get_eos_set()) { return false; }
        }
        return true;
    }

    // Whole-planet equation-of-state solve. Integrates gravity, pressure, enclosed mass, and moment of
    // inertia from center to surface with each layer's material EOS model as the local density source, then
    // populates every layer's c_LayerEOSData and mass. All MKS.
    //
    // Assumes spherical symmetry, and that each layer's density comes from its material EOS model: the
    // pressure for the analytic models, the radius for the interpolated one.
    void solve_eos(const c_WorldEOSSolveConfig& cfg) {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        const std::size_t n_layers = this->p_layers.size();
        if (n_layers == 0) {
            throw std::invalid_argument("TidalPy: cannot solve the EOS for a world with no layers");
        }
        if (!this->get_all_eos_set()) {
            throw std::invalid_argument(
                "TidalPy: every layer must have a material EOS model attached (set_eos) "
                "before the world EOS can be solved");
        }
        if (cfg.slices_per_layer < 2) {
            throw std::invalid_argument("TidalPy: slices_per_layer must be at least 2");
        }

        double G_to_use = cfg.G_to_use;
        if (G_to_use < 0.0) {
            G_to_use = c_get_G();
        }

        const std::size_t slices       = cfg.slices_per_layer;
        const std::size_t total_slices = n_layers * slices;

        // A concatenated per-layer ascending radius grid, and the bulk density guess as the volume-weighted
        // EOS density at the surface pressure.
        std::vector<double> full_radius(total_slices);
        std::vector<double> upper_radii(n_layers);
        double total_volume  = 0.0;
        double mass_estimate = 0.0;
        // What a floating layer holds on to: its mass, from its configuration or from the first solve.
        if (cfg.reset_layer_masses) { this->p_reference_mass.clear(); }
        this->p_reference_mass.resize(n_layers, TidalPyConstants::d_NAN);
        for (std::size_t i = 0; i < n_layers; ++i) {
            c_BaseLayer* layer   = this->p_layers[i].get();
            const double r_inner = layer->get_radius_inner();
            const double r_outer = layer->get_radius_outer();
            upper_radii[i]       = r_outer;
            if (!layer->get_is_volume_fixed() && !cfg.reset_layer_masses
                && !std::isfinite(this->p_reference_mass[i])
                && (layer->get_mass() > TidalPyConstants::d_EPS)) {
                this->p_reference_mass[i] = layer->get_mass();
            }
            for (std::size_t s = 0; s < slices; ++s) {
                const double frac = static_cast<double>(s) / static_cast<double>(slices - 1);
                full_radius[i * slices + s] = r_inner + frac * (r_outer - r_inner);
            }
            const double r_mid     = 0.5 * (r_inner + r_outer);
            const double rho_mid   = layer->get_eos()->calc_density(
                cfg.surface_pressure, TidalPyConstants::d_NAN, r_mid);
            const double shell_vol = (4.0 / 3.0) * TidalPyConstants::d_PI
                                   * (r_outer * r_outer * r_outer - r_inner * r_inner * r_inner);
            total_volume  += shell_vol;
            mass_estimate += rho_mid * shell_vol;
        }
        const double planet_bulk_density = (total_volume > TidalPyConstants::d_EPS) ? (mass_estimate / total_volume) : 3500.0;
        const double planet_radius       = upper_radii.back();
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
        for (double& radius : full_radius) { radius /= length_scale; }
        for (double& radius : upper_radii) { radius /= length_scale; }
        const double G_solve = G_to_use / (length_scale * length_scale * length_scale / (mass_scale * second2_scale));
        const double surface_pressure_solve = cfg.surface_pressure / pascal_scale;
        const double bulk_density_solve     = planet_bulk_density / density_scale;

        // Per-layer pre-eval functions and inputs. The solution object is rebuilt once per thermal pass.
        std::vector<PreEvalFunc>    eos_function_vec;
        std::vector<c_EOS_ODEInput> eos_input_vec;
        eos_function_vec.reserve(n_layers);
        eos_input_vec.reserve(n_layers);

        // Everything this solve evaluates its materials with is its own: a copy of each layer's material model,
        // the per-layer inputs the structure ODE reaches it through (the solver keeps a copy of each
        // c_EOS_ODEInput, whose eos_input_ptr points into these), and the heat sources. The solution co-owns
        // it, so a retained or exported solution keeps reading exactly what was solved, whatever later happens
        // to the layers, their models, or another solve.
        auto solve_state = std::make_shared<c_EOSSolveState>();
        solve_state->materials.reserve(n_layers);
        solve_state->inputs.assign(n_layers, c_MaterialEOSInput());

        c_EOS_ODEInput ode_input;
        ode_input.G_to_use      = G_solve;
        ode_input.planet_radius = upper_radii.back();
        ode_input.update_bulk   = false;
        ode_input.update_shear  = false;
        for (std::size_t i = 0; i < n_layers; ++i) {
            solve_state->materials.push_back(c_clone_material_eos(*this->p_layers[i]->get_eos()));
            c_MaterialEOSInput& material_input = solve_state->inputs[i];
            material_input.eos_model_ptr = solve_state->materials[i].get();
            material_input.length_scale  = length_scale;
            material_input.pascal_scale  = pascal_scale;
            material_input.density_scale = density_scale;
            ode_input.eos_input_ptr = reinterpret_cast<char*>(&material_input);
            eos_function_vec.push_back(c_preeval_material_eos);
            eos_input_vec.push_back(ode_input);
        }

        // Thermal layout. Each layer carries its own temperature and its cooling model says how heat moves
        // inside it; a uniform override replaces every layer's value.
        std::vector<c_LayerThermal> layer_thermal;
        c_init_layer_thermal(this->p_layers, layer_thermal, cfg.temperature);
        // Heat sources, prepared once: none of them depends on the solved state. They need the heat flow of a
        // thermal solve to act through, so a solve with temperature switched off leaves them out.
        c_WorldState world_state;
        world_state.time       = cfg.time;
        world_state.layers_ptr = &this->p_layers;
        c_Heating& heating = solve_state->heating;
        heating.update_sources(world_state, length_scale, density_scale);
        const bool heating_active = heating.get_is_active();
        if (heating_active && !cfg.solve_temperature) {
            TIDALPY_LOG_WARN(
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
        const double gravity_scale = length_scale / second2_scale;

        // Pass 0 is isothermal at each layer's own temperature, which is the whole solve for a world with no
        // temperature contrast. Each later pass integrates the profile, then relaxes the boundary layers,
        // interface temperatures, and heat flows against the structure it produced.
        // A layer holding its mass lets its radii float, which moves every layer above it, so the same
        // passes relax the geometry.
        bool geometry_floats = false;
        for (const auto& layer_uptr : this->p_layers) {
            if (!layer_uptr->get_is_volume_fixed()) { geometry_floats = true; }
        }
        std::vector<c_EOSSegment> segment_vec;
        std::shared_ptr<c_EOSSolution> solution;
        // The central pressure the secant iteration starts from, in solve units: the world's last converged solve
        // (kept when a layer change forgets the solved state), then the pass before. A re-solve after a small change
        // then converges in a pass or two. NaN, on a world that has never been solved, starts from a uniform sphere.
        double central_pressure_guess = this->p_warm_start_central_pressure / pascal_scale;
        const std::size_t last_pass =
            (thermal_contrast || geometry_floats) ? cfg.max_thermal_passes : 0;
        std::size_t thermal_passes = 0;
        bool thermal_converged     = !thermal_contrast;
        bool geometry_converged    = !geometry_floats;
        // Floating layers move during the passes; a solve that fails or throws puts them back.
        std::vector<double> radius_inner_before(n_layers);
        std::vector<double> radius_outer_before(n_layers);
        for (std::size_t i = 0; i < n_layers; ++i) {
            radius_inner_before[i] = this->p_layers[i]->get_radius_inner();
            radius_outer_before[i] = this->p_layers[i]->get_radius_outer();
        }
        const double world_radius_before = this->p_radius;
        // The radii each pass solves on. A pass moves floating layers after it solves, so the layers go back onto
        // the last pass's grid at the end: every slice and interface then sits in its own layer.
        std::vector<double> radius_inner_solved = radius_inner_before;
        std::vector<double> radius_outer_solved = radius_outer_before;
        double world_radius_solved = world_radius_before;
        const auto restore_radii = [&]() {
            if (!geometry_floats) { return; }
            for (std::size_t i = 0; i < n_layers; ++i) {
                this->p_layers[i]->set_radii(radius_inner_before[i], radius_outer_before[i]);
            }
            this->p_radius = world_radius_before;
        };
        try {
            for (std::size_t pass = 0; pass <= last_pass; ++pass) {
                const bool integrate_temperature = thermal_contrast && (pass > 0);
                if (pass > 0) {
                    // The layers may have moved, so the grid and the segment bounds follow them.
                    for (std::size_t i = 0; i < n_layers; ++i) {
                        const c_BaseLayer* layer = this->p_layers[i].get();
                        const double r_inner = layer->get_radius_inner() / length_scale;
                        const double r_outer = layer->get_radius_outer() / length_scale;
                        upper_radii[i] = r_outer;
                        for (std::size_t s = 0; s < slices; ++s) {
                            const double frac = static_cast<double>(s) / static_cast<double>(slices - 1);
                            full_radius[i * slices + s] = r_inner + frac * (r_outer - r_inner);
                        }
                    }
                    // The structure derivatives vanish above the planet radius, so a world that grew has to say so.
                    for (std::size_t i = 0; i < n_layers; ++i) {
                        eos_input_vec[i].planet_radius = upper_radii.back();
                    }
                    for (std::size_t i = 0; i < n_layers; ++i) {
                        radius_inner_solved[i] = this->p_layers[i]->get_radius_inner();
                        radius_outer_solved[i] = this->p_layers[i]->get_radius_outer();
                    }
                    world_radius_solved = this->p_radius;
                }
                c_build_thermal_segments(
                    layer_thermal, this->p_layers, integrate_temperature,
                    length_scale, gravity_scale, segment_vec);
                for (std::size_t i = 0; i < n_layers; ++i) {
                    // The viscosity and melt models of the material always see the temperature; its density law
                    // sees it only when the layer asked for a thermal EOS.
                    const auto* physics_layer = dynamic_cast<const c_PhysicsLayer*>(this->p_layers[i].get());
                    const bool thermal_eos = (physics_layer != nullptr) && physics_layer->get_use_thermal_eos();
                    solve_state->inputs[i].temperature           = layer_thermal[i].temperature;
                    solve_state->inputs[i].use_state_temperature = integrate_temperature;
                    solve_state->inputs[i].thermal_density       = thermal_eos;
                }

                solution = std::make_shared<c_EOSSolution>(
                    upper_radii.data(), n_layers, full_radius.data(), total_slices);
                solution->input_keepalive = solve_state;
                c_solve_eos(
                    solution.get(),
                    eos_function_vec,
                    eos_input_vec,
                    bulk_density_solve,
                    surface_pressure_solve,
                    G_solve,
                    cfg.integration_method,
                    cfg.rtol,
                    cfg.atol,
                    cfg.pressure_tol,
                    cfg.max_iters,
                    cfg.verbose,
                    &segment_vec,
                    integrate_temperature,
                    central_pressure_guess
                );

                // Return the solution to SI: the arrays, layer radii, and pressure error are scaled in place, and
                // every later evaluation of the retained integrators (call_si) converts on the way in and out.
                if (cfg.nondimensionalize) {
                    solution->dimensionalize_data(scales_uptr.get(), true);
                }
                if (!solution->success) { break; }
                thermal_passes = pass;
                central_pressure_guess = solution->central_pressure / pascal_scale;

                if (thermal_contrast) {
                    const double thermal_change = c_update_layer_thermal(
                        *solution, this->p_layers, cfg.surface_temperature, integrate_temperature,
                        layer_thermal, &heating);
                    thermal_converged = integrate_temperature && (thermal_change < cfg.thermal_tol);
                }
                if (geometry_floats) {
                    geometry_converged = (this->update_floating_radii(*solution) < cfg.radius_tol);
                }
                if (thermal_converged && geometry_converged) { break; }
            }
        } catch (...) {
            // Nothing half-solved survives a throw: the layers go back where they were and the world is unsolved.
            restore_radii();
            this->p_reset_solved_state();
            throw;
        }

        // Commit. Whatever was solved on top of the previous structure no longer describes the world.
        this->mark_structure_dirty();
        this->p_thermal_passes     = thermal_passes;
        this->p_thermal_converged  = thermal_converged;
        this->p_geometry_converged = geometry_converged;

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

        // Store scalar results.
        this->p_eos_success          = solution->success;
        this->p_eos_message          = solution->message;
        this->p_eos_iterations       = solution->iterations;
        this->p_eos_max_iters_hit    = solution->max_iters_hit;
        this->p_eos_pressure_error   = solution->pressure_error;
        this->p_surface_gravity_eos  = solution->surface_gravity;
        this->p_surface_pressure_eos = solution->surface_pressure;
        this->p_central_pressure     = solution->central_pressure;
        this->p_planet_mass_eos      = solution->mass;
        this->p_planet_moi_eos       = solution->moi;
        this->p_eos_solved           = solution->success && solution->other_vecs_set;

        if (!this->p_eos_solved) {
            // A failed solve leaves the world unsolved rather than holding the previous structure: the layers go
            // back where they were, and no layer profile, thermal state, or radial setup from before survives.
            // The failed solution stays for its diagnostics only.
            restore_radii();
            this->p_layer_thermal.clear();
            this->p_solve_state.reset();
            this->p_love = c_LoveWorkspace();
            for (const auto& layer_uptr : this->p_layers) { layer_uptr->clear_eos_data(); }
            this->p_eos_solution = std::move(solution);
            return;
        }
        this->p_warm_start_central_pressure = this->p_central_pressure;
        this->p_layer_thermal = std::move(layer_thermal);
        this->p_solve_state   = solve_state;

        // Floating layers go back onto the grid the last pass solved on; once converged the move this undoes is
        // below radius_tol.
        if (geometry_floats) {
            for (std::size_t i = 0; i < n_layers; ++i) {
                this->p_layers[i]->set_radii(radius_inner_solved[i], radius_outer_solved[i]);
            }
            this->p_radius = world_radius_solved;
        }
        if (!thermal_converged) {
            TIDALPY_LOG_WARN(
                "TidalPy: world '{}' EOS solve used all {} thermal passes before its interface temperatures and heat "
                "flows settled to thermal_tol = {:.1e}; the temperature profile is from the last pass.",
                this->get_name(), cfg.max_thermal_passes, cfg.thermal_tol);
        }
        if (!geometry_converged) {
            TIDALPY_LOG_WARN(
                "TidalPy: world '{}' EOS solve used all {} passes before its floating layer radii settled to "
                "radius_tol = {:.1e}; the structure is from the last pass.",
                this->get_name(), cfg.max_thermal_passes, cfg.radius_tol);
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

        // A layer whose pressure passes what its material's pressure law represents (Birch-Murnaghan with K0' below
        // 4 turns over at a finite compression) has its density held at the law's largest compression there, and a
        // bulk modulus near zero; the solve still succeeds, so say so.
        for (std::size_t layer_index = 0; layer_index < n_layers; ++layer_index) {
            const double law_max_pressure = solve_state->materials[layer_index]->get_max_pressure();
            const std::size_t base_slice = layer_index * slices;
            if (!std::isfinite(law_max_pressure) || (base_slice >= solution->pressure_array_vec.size())) { continue; }
            const double layer_max_pressure = solution->pressure_array_vec[base_slice];
            if (layer_max_pressure > law_max_pressure) {
                TIDALPY_LOG_WARN(
                    "TidalPy: layer '{}' of world '{}' reaches {:.4e} Pa, past the {:.4e} Pa its material's pressure "
                    "law represents; its density is held at the law's largest compression there. Check the law's "
                    "bulk modulus derivative (a Birch-Murnaghan K0' below 4 turns over).",
                    this->p_layers[layer_index]->get_name(), this->get_name(), layer_max_pressure, law_max_pressure);
            }
        }

        // Populate each layer's structure and viscoelastic profile from its slice of the full arrays.
        {
            const std::size_t total_slices = solution->radius_array_size;
            for (std::size_t layer_index = 0; layer_index < n_layers; ++layer_index) {
                const std::size_t slice_start = layer_index * slices;
                const std::size_t slice_end   = slice_start + slices;
                if (slice_end > total_slices) { break; }
                c_BaseLayer* layer = this->p_layers[layer_index].get();

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
        this->update_radial_segments();
    }

    // The stretches the radial solver integrates, inner to outer, from the last successful EOS solve. Each layer
    // is one stretch unless its partial-melt model makes part of it molten: the post-melt shear modulus sits at the
    // model's liquid_shear floor, or its rigidity mu / (rho g R) is below the config's minimum_solid_rigidity, which
    // takes in the steep weakening just short of the floor. Molten stretches are found on the EOS slices, their
    // edges refined by bisection on the dense profile, and the layer is split there, each edge on the solid side so
    // no solid stretch reaches into the melt. Whether a molten stretch is solved as a liquid is decided per Love
    // solve (ensure_radial_cache), since the layer flags can change without a new EOS solve.
    void update_radial_segments() {
        this->p_radial_segments.clear();
        const c_EOSSolution* solution = this->p_eos_solution.get();
        const std::size_t n_layers = this->p_layers.size();
        if (!this->p_eos_solved || (solution == nullptr) || (n_layers == 0)) { return; }
        const std::size_t slices = solution->radius_array_size / n_layers;
        // The rigidity scale rho g R of the planet: bulk density, surface gravity, radius.
        const double planet_radius = solution->radius;
        const double planet_volume =
            (4.0 / 3.0) * TidalPyConstants::d_PI * planet_radius * planet_radius * planet_radius;
        const double rigidity_scale = (planet_volume > TidalPyConstants::d_EPS)
            ? (this->p_planet_mass_eos / planet_volume) * this->p_surface_gravity_eos * planet_radius : 0.0;
        const double min_rigidity = tidalpy_config_ptr->d_MIN_SOLID_RIGIDITY;
        const double weak_shear = (std::isfinite(min_rigidity) && std::isfinite(rigidity_scale))
            ? min_rigidity * rigidity_scale : 0.0;   // [Pa]

        for (std::size_t layer_i = 0; layer_i < n_layers; ++layer_i) {
            const c_BaseLayer* layer  = this->p_layers[layer_i].get();
            const double radius_inner = layer->get_radius_inner();
            const double radius_outer = layer->get_radius_outer();
            const auto* physics_layer = dynamic_cast<const c_PhysicsLayer*>(layer);
            // The melt model the solve used, which the layer may have replaced since.
            const c_PartialMeltBase* melt_model =
                ((physics_layer != nullptr) && this->p_solve_state && (layer_i < this->p_solve_state->materials.size()))
                ? this->p_solve_state->materials[layer_i]->get_partial_melt_model() : nullptr;
            if ((melt_model == nullptr) || (slices < 2)) {
                this->p_radial_segments.push_back({layer_i, radius_inner, radius_outer, false});
                continue;
            }

            const double molten_shear =
                std::max(melt_model->get_liquid_shear() * (1.0 + d_MOLTEN_SHEAR_RTOL), weak_shear);   // [Pa]
            const auto is_molten = [layer, molten_shear](double radius) {
                return layer->get_shear_modulus(radius) <= molten_shear;
            };
            // An edge is refined until its bracket is as narrow as the tolerance on a layer interface, and a
            // stretch no thicker than that is part of its neighbor.
            const double edge_tolerance = d_LAYER_BOUNDARY_RTOL * radius_outer;
            const std::size_t first_segment = this->p_radial_segments.size();
            const double* slice_radius = solution->radius_array_vec.data() + layer_i * slices;
            double stretch_start  = radius_inner;
            bool   stretch_molten = is_molten(slice_radius[0]);
            for (std::size_t slice_i = 1; slice_i < slices; ++slice_i) {
                const bool molten_here = is_molten(slice_radius[slice_i]);
                if (molten_here == stretch_molten) { continue; }
                double lower = slice_radius[slice_i - 1];
                double upper = slice_radius[slice_i];
                while (upper - lower > edge_tolerance) {
                    const double middle = 0.5 * (lower + upper);
                    if (is_molten(middle) == stretch_molten) { lower = middle; } else { upper = middle; }
                }
                // The bracket end that is solid, so the molten stretch takes the unresolved sliver.
                const double edge = stretch_molten ? upper : lower;
                if (edge - stretch_start > edge_tolerance) {
                    this->p_radial_segments.push_back({layer_i, stretch_start, edge, stretch_molten});
                    stretch_start = edge;
                }
                stretch_molten = molten_here;
            }
            this->p_radial_segments.push_back({layer_i, stretch_start, radius_outer, stretch_molten});

            // A skipped sliver can leave two neighbors in the same state; they are one stretch.
            std::vector<c_RadialSegment> merged;
            for (std::size_t segment_i = first_segment; segment_i < this->p_radial_segments.size(); ++segment_i) {
                const c_RadialSegment& segment = this->p_radial_segments[segment_i];
                if (!merged.empty() && (merged.back().molten == segment.molten)) {
                    merged.back().radius_outer = segment.radius_outer;
                } else {
                    merged.push_back(segment);
                }
            }
            this->p_radial_segments.resize(first_segment);
            this->p_radial_segments.insert(this->p_radial_segments.end(), merged.begin(), merged.end());

            for (const c_RadialSegment& segment : merged) {
                if (segment.molten && physics_layer->get_is_solid()) {
                    TIDALPY_LOG_INFO(
                        "TidalPy: layer '{}' of world '{}' is molten between {:.6e} and {:.6e} m; the radial solver "
                        "treats that stretch as a static liquid.",
                        layer->get_name(), this->get_name(), segment.radius_inner, segment.radius_outer);
                }
            }
        }
    }

    // Empty before a successful EOS solve.
    const std::vector<c_RadialSegment>& get_radial_segments() const noexcept { return this->p_radial_segments; }

    // Whether a radius [m] lies inside a molten stretch of a solid layer, which the radial solver treats as a
    // static liquid.
    bool get_is_molten_at(double radius) const noexcept {
        for (const c_RadialSegment& segment : this->p_radial_segments) {
            if (segment.molten && (radius > segment.radius_inner) && (radius < segment.radius_outer)) { return true; }
        }
        return false;
    }

    // The molten stretches the radial solver treats as static liquids: those of the layers that are solid.
    std::vector<c_RadialSegment> get_molten_regions() const {
        std::vector<c_RadialSegment> regions;
        for (const c_RadialSegment& segment : this->p_radial_segments) {
            if (!segment.molten) { continue; }
            const auto* physics_layer =
                dynamic_cast<const c_PhysicsLayer*>(this->p_layers[segment.world_layer].get());
            if ((physics_layer != nullptr) && physics_layer->get_is_solid()) { regions.push_back(segment); }
        }
        return regions;
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
        this->p_radial_segments.clear();
    }

    // A layer this world owns was moved through a view. The solved profile no longer lines up with the layers, so
    // it is forgotten. A change of material leaves the solved state alone: the solve ran on copies of the
    // materials, and it describes the world until the next solve_eos.
    void update_after_layer_geometry_change() {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        this->p_reset_solved_state();
    }

protected:
    // Forget everything solved: the EOS solution, the layer profiles, the thermal state, and every result built
    // on them. Called when the layers no longer match what was solved (a layer added, a binary load) and when a
    // solve throws, so no reader can mistake an old structure for the current one.
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
    // Solve the EOS and copy its result out before any other thread can start a solve on this world.
    c_WorldEOSReport solve_eos_report(const c_WorldEOSSolveConfig& cfg) {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        this->solve_eos(cfg);
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
        report.max_iters_hit      = this->p_eos_max_iters_hit;
        report.pressure_error     = this->p_eos_pressure_error;
        report.surface_gravity    = this->p_surface_gravity_eos;
        report.surface_pressure   = this->p_surface_pressure_eos;
        report.central_pressure   = this->p_central_pressure;
        report.planet_mass        = this->p_planet_mass_eos;
        report.planet_moi         = this->p_planet_moi_eos;
        report.thermal_passes     = this->p_thermal_passes;
        report.thermal_converged  = this->p_thermal_converged;
        report.geometry_converged = this->p_geometry_converged;

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
        report.layer_temperature_rate.reserve(num_layers);
        report.layer_radius_outer.reserve(num_layers);
        for (std::size_t layer_i = 0; layer_i < num_layers; ++layer_i) {
            report.layer_temperature_rate.push_back(this->calc_layer_temperature_rate(layer_i));
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
    // spin-rate change uses the structure-resolved value rather than the uniform-density one.
    void           set_spin_model(const c_Spin& spin) noexcept { this->p_spin = spin; }
    const c_Spin&  get_spin_model() const noexcept { return this->p_spin; }

    // The EOS-solved value once the EOS has been solved, else the spin model's factor * M R^2 estimate.
    double get_moment_of_inertia() const noexcept {
        if (this->p_eos_solved && std::isfinite(this->p_planet_moi_eos)) {
            return this->p_planet_moi_eos;
        }
        return this->p_spin.calc_moment_of_inertia(this->get_mass(), this->get_radius());
    }

    // M_host * dU/dO / I, from the last calc_tides' dU/dO and the world's moment of inertia.
    double calc_spin_derivative(double host_mass) const {
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

        // The solver's layers: every world layer once, except that a solid layer with molten stretches is split
        // at their edges and each molten stretch is a static liquid (it keeps the layer's compressibility flag).
        // The static formulation reads only density and gravity, so it does not see the melt-weakened bulk
        // modulus there. Geometry-only layers default to static solid. Gathered before the cache check because
        // the flags are user-mutable without an EOS re-solve, so a cache hit is only valid when they still match.
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
            const c_BaseLayer* layer = this->p_layers[world_layer].get();
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

        std::size_t segment_i = 0;
        const std::vector<c_RadialSegment>& segments = this->p_radial_segments;
        for (std::size_t layer_i = 0; layer_i < n_layers; ++layer_i) {
            const c_BaseLayer* layer = this->p_layers[layer_i].get();
            const auto* phys = dynamic_cast<const c_PhysicsLayer*>(layer);
            const bool is_solid          = (phys == nullptr) || phys->get_is_solid();
            const bool is_static         = (phys == nullptr) || phys->get_is_static();
            const bool is_incompressible = (phys != nullptr) && phys->get_is_incompressible();

            std::size_t layer_segment_end = segment_i;
            while ((layer_segment_end < segments.size()) && (segments[layer_segment_end].world_layer == layer_i)) {
                ++layer_segment_end;
            }
            if (!is_solid || (layer_segment_end - segment_i < 2)) {
                // One stretch: the layer's own flags, or a static liquid when a solid layer is molten throughout.
                const bool molten = is_solid && (layer_segment_end > segment_i) && segments[segment_i].molten;
                add_solver_layer(layer_i, (is_solid && !molten) ? 0 : 1, molten || is_static, is_incompressible,
                                 layer->get_radius_inner(), layer->get_radius_outer());
            } else {
                for (std::size_t stretch_i = segment_i; stretch_i < layer_segment_end; ++stretch_i) {
                    const c_RadialSegment& segment = segments[stretch_i];
                    add_solver_layer(layer_i, segment.molten ? 1 : 0, segment.molten || is_static, is_incompressible,
                                     segment.radius_inner, segment.radius_outer);
                }
            }
            segment_i = layer_segment_end;
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
        rt.use_kamata         = cfg.use_kamata;
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
        return rt;
    }

    // The world's own solve: its results are what the Love getters report and release_radial_storage hands out.
    void solve_love_numbers(const c_LoveSolveConfig& cfg) {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
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

    // Composite Simpson intervals per tidal layer for the quasi-homogeneous averages; even, so 129 nodes.
    static constexpr std::size_t homogeneous_quadrature_intervals = 128;

    // One tidal layer as the quasi-homogeneous Love methods see it: its volume-averaged post-melt shear and bulk
    // moduli, its log-volume-averaged post-melt viscosities, and its tidal scale. Frequency independent.
    struct c_HomogeneousLayer {
        const c_PhysicsLayer* layer = nullptr;
        std::size_t layer_index     = 0;
        double shear_modulus   = TidalPyConstants::d_NAN;   // [Pa]
        double bulk_modulus    = TidalPyConstants::d_NAN;   // [Pa]
        double shear_viscosity = TidalPyConstants::d_NAN;   // [Pa s]
        double bulk_viscosity  = TidalPyConstants::d_NAN;   // [Pa s]
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
        solver->set_material_eval(this->make_material_eval(cfg.frequency, workspace.radial_world_layer));

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
        solver->set_material_eval(
            [eos_solution, radius_copy, shear_copy, bulk_copy, in_first_by_layer, in_count_by_layer, world_layer_of](
                    std::size_t solver_layer_index,
                    double radius_si,
                    double* state_out,
                    std::complex<double>& shear_out,
                    std::complex<double>& bulk_out) {
                // The solver counts the stretches of a split layer as layers of their own.
                const std::size_t layer_index = world_layer_of[solver_layer_index];
                // One dense read gives the structure and the density; the supplied profile gives the moduli.
                eos_solution->call_si(layer_index, radius_si, state_out);
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
    // frequency-independent state, then the layer's rheology, the only part that knows the frequency.
    //
    // The callable co-owns the solved EOS and the rheologies, and resolves the layers here, once, rather
    // than at every radius. The solved EOS still reaches each layer's material model through this world, so
    // the world has to outlive the callable; a solution exported to Python holds its world for that reason.
    c_EOSSolution::MaterialEval make_material_eval(
            double frequency,
            const std::vector<std::size_t>& world_layer_of) const {
        const std::size_t n_layers = this->p_layers.size();
        std::vector<std::shared_ptr<const c_RheologyBase>> shear_bylayer(n_layers);
        std::vector<std::shared_ptr<const c_RheologyBase>> bulk_bylayer(n_layers);
        std::vector<char> is_physics_bylayer(n_layers, 0);
        for (std::size_t layer_i = 0; layer_i < n_layers; ++layer_i) {
            const auto* physics_layer = dynamic_cast<const c_PhysicsLayer*>(this->p_layers[layer_i].get());
            if (physics_layer == nullptr) { continue; }
            is_physics_bylayer[layer_i] = 1;
            shear_bylayer[layer_i]      = physics_layer->share_shear_rheology();
            bulk_bylayer[layer_i]       = physics_layer->share_bulk_rheology();
        }
        std::shared_ptr<const c_EOSSolution> eos_solution = this->p_eos_solution;
        return [eos_solution, shear_bylayer, bulk_bylayer, is_physics_bylayer, world_layer_of, frequency](
                std::size_t solver_layer_index,
                double radius_si,
                double* state_out,
                std::complex<double>& shear_out,
                std::complex<double>& bulk_out) {
            // The solver counts the stretches of a split layer as layers of their own.
            const std::size_t layer_index = world_layer_of[solver_layer_index];
            eos_solution->call_si(layer_index, radius_si, state_out);
            // A layer with no material models has no modulus to report.
            if (layer_index >= is_physics_bylayer.size() || !is_physics_bylayer[layer_index]) { return; }
            const double static_shear = state_out[C_EOS_SHEAR_MODULUS_INDEX];
            const double static_bulk  = state_out[C_EOS_BULK_MODULUS_INDEX];
            // Purely real (no dissipation) where no rheology is attached.
            shear_out = shear_bylayer[layer_index]
                ? shear_bylayer[layer_index]->calc_complex_modulus(
                    static_shear, state_out[C_EOS_SHEAR_VISCOSITY_INDEX], frequency)
                : std::complex<double>(static_shear, 0.0);
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
    // A physics layer takes part with the moduli its material gives, so a fluid layer adds its fluid Love number;
    // a geometry-only layer, one that is not tidal, or one with a zero tidal scale takes no part.
    void build_homogeneous_layers(c_HomogeneousLoveCache& cache) const {
        const std::size_t n_intervals = homogeneous_quadrature_intervals;
        const double planet_radius = this->get_radius();
        const double planet_volume =
            (4.0 / 3.0) * TidalPyConstants::d_PI * planet_radius * planet_radius * planet_radius;
        cache.layers.clear();
        for (std::size_t layer_i = 0; layer_i < this->p_layers.size(); ++layer_i) {
            const c_BaseLayer* layer = this->p_layers[layer_i].get();
            const double tidal_scale = layer->calc_tidal_scale(planet_volume);
            if (!(tidal_scale > 0.0)) { continue; }
            const auto* physics = dynamic_cast<const c_PhysicsLayer*>(layer);
            if (physics == nullptr) { continue; }   // geometry-only layer: no material moduli
            const double r_inner = layer->get_radius_inner();
            const double r_outer = layer->get_radius_outer();
            if (!(r_outer > r_inner)) { continue; }

            const double dr = (r_outer - r_inner) / static_cast<double>(n_intervals);
            double weight_sum          = 0.0;
            double shear_sum           = 0.0;
            double bulk_sum            = 0.0;
            double log_shear_visc_sum  = 0.0;
            double log_bulk_visc_sum   = 0.0;
            // One dense evaluation per node gives all four quantities.
            double state[C_EOS_DY_VALUES];
            for (std::size_t i = 0; i <= n_intervals; ++i) {
                const double r = (i == n_intervals) ? r_outer : r_inner + static_cast<double>(i) * dr;
                const double simpson = (i == 0 || i == n_intervals) ? 1.0 : ((i % 2 == 1) ? 4.0 : 2.0);
                const double weight = simpson * r * r;
                physics->get_eos_state(r, state);
                weight_sum         += weight;
                shear_sum          += weight * state[C_EOS_SHEAR_MODULUS_INDEX];     // post-melt
                bulk_sum           += weight * state[C_EOS_BULK_MODULUS_INDEX];
                log_shear_visc_sum += weight * std::log10(state[C_EOS_SHEAR_VISCOSITY_INDEX]);
                log_bulk_visc_sum  += weight * std::log10(state[C_EOS_BULK_VISCOSITY_INDEX]);
            }
            c_HomogeneousLayer averaged;
            averaged.layer           = physics;
            averaged.layer_index     = layer_i;
            averaged.shear_modulus   = shear_sum / weight_sum;
            averaged.bulk_modulus    = bulk_sum / weight_sum;
            averaged.shear_viscosity = std::pow(10.0, log_shear_visc_sum / weight_sum);
            averaged.bulk_viscosity  = std::pow(10.0, log_bulk_visc_sum / weight_sum);
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
                "TidalPy: the homogeneous Love-number methods need at least one tidal layer (is_tidal, with a nonzero "
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
            const c_RheologyBase* rheology = averaged.layer->get_shear_rheology_model();
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
    void set_eos_solver_overrides(const c_EOSSolverOverrides& overrides) noexcept {
        this->p_eos_solver_overrides = overrides;
    }
    void set_radial_solver_overrides(const c_RadialSolverOverrides& overrides) noexcept {
        this->p_radial_solver_overrides = overrides;
    }
    const c_EOSSolverOverrides& get_eos_solver_overrides() const noexcept {
        return this->p_eos_solver_overrides;
    }
    const c_RadialSolverOverrides& get_radial_solver_overrides() const noexcept {
        return this->p_radial_solver_overrides;
    }

    // The [eos_solver] section of the TidalPy configuration with the world's pinned keys on top.
    c_WorldEOSSolveConfig make_eos_solve_config() const {
        c_WorldEOSSolveConfig cfg;
        this->p_eos_solver_overrides.apply(cfg);
        return cfg;
    }

    // The world's configured method and its cpl / ctl parameters from the [tides] config, plus its pinned
    // [radial_solver] keys. The tide paths start from this, so the configured method drives every solve.
    c_LoveSolveConfig make_love_solve_config() const {
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
    // of the layers' complex shear moduli [Pa], and the volume of the layers that took part [m3].
    std::complex<double> get_love_analytic_shear() const noexcept {
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
        if (!this->p_love.analytic_success) { return TidalPyConstants::d_NAN; }
        double volume = 0.0;
        for (const c_LayerLove& part : this->p_love.analytic_layers) { volume += part.volume; }
        return volume;
    }
    // Each tidal layer's part of the last quasi-homogeneous solve; empty otherwise.
    const std::vector<c_LayerLove>& get_love_layer_parts() const noexcept { return this->p_love.analytic_layers; }

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
    // analytic methods. See c_estimate_surface_amplification in RadialSolver_x/boundaries/boundaries_.hpp.
    double get_love_surface_amplification() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_love.get_surface_amplification();
    }
    // Reciprocal condition number of the surface boundary condition system; NaN before a shooting-method solve and
    // for the other methods. See c_estimate_surface_rcond in RadialSolver_x/boundaries/boundaries_.hpp.
    double get_love_surface_rcond() const noexcept {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        return this->p_love.get_surface_rcond();
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

    // Global (1D) tidal dissipation. The tide-model holder, config, results, and the analytic calc_tides path
    // live on c_BaseWorld; this class hides calc_tides to add the rheology path and the per-layer heating
    // distribution, out-of-line in world_tides_.hpp, which carries the heavy global-potential engine.
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
    // efficient way to build a map. num_threads applies to the per-point evaluation after the solves.
    void get_3d_tidal_heating_array(
            const c_TideSolveConfig& state,
            const double* radii,
            const double* colatitudes,
            size_t num_points,
            double* out_heating,
            int num_threads = 1);

    // Instantaneous tidal displacements [m] on the (radius, colatitude, longitude, time) grid; see
    // c_RheologyTide::calc_3d_displacements_grid. Same preconditions as get_3d_tidal_heating.
    void get_3d_displacements_grid(
            const c_TideSolveConfig& state,
            const c_Grid3DAxes& axes,
            double* out_disp,
            int num_threads = 1);

    // Instantaneous stress [Pa] and strain on the same grid; see c_RheologyTide::calc_3d_stress_strain_grid.
    void get_3d_stress_strain_grid(
            const c_TideSolveConfig& state,
            const c_Grid3DAxes& axes,
            double* out_stress,
            double* out_strain,
            int num_threads = 1);

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

    // calc_3d_tides written into caller buffers. The radial solves run on the calling thread and the
    // per-point evaluation on up to cfg.num_threads threads; the result is identical for any thread count.
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
    // orbit-averaged heating density over the layer, scaled so the layers sum to `total_heating` (defined in
    // world_tides_.hpp).
    void calc_layer_tidal_heating_radial(
            const c_TideSolveConfig& state,
            double total_heating,
            std::vector<double>& out,
            const std::vector<c_RetainedRadialSolve>* retained_solves = nullptr);

    // The radial solves calc_tides lends its per-layer heating integral; null outside that call (see
    // c_RetainedRadialSolve). The 3D radial-group solve takes a matching one instead of solving again.
    const std::vector<c_RetainedRadialSolve>* get_retained_radial_solves() const noexcept {
        return this->p_retained_radial_solves;
    }

    // The tidal scale of each layer in the quasi-homogeneous Love methods (c_BaseLayer::calc_tidal_scale).
    double get_layer_tidal_scale(std::size_t index) const {
        const double planet_radius = this->get_radius();
        const double planet_volume =
            (4.0 / 3.0) * TidalPyConstants::d_PI * planet_radius * planet_radius * planet_radius;
        return this->get_layer(index)->calc_tidal_scale(planet_volume);
    }

    // The heating [W] the last calc_tides put in a layer; NaN before one, or when calc_tides was not asked for it.
    double get_layer_tidal_heating(std::size_t index) const noexcept {
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

    // Only solid-liquid layers produce radiogenic heating; the rest contribute zero.
    double calc_internal_heating(double time) const noexcept {
        double total = 0.0;
        for (const auto& layer : this->p_layers) {
            const auto* sl = dynamic_cast<const c_SolidLiquidLayer*>(layer.get());
            if (sl != nullptr) {
                total += sl->calc_radiogenic_heating(time, sl->get_mass());
            }
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

    void write_binary(std::ostream& out) const override {
        this->write_layered_binary(out, static_cast<uint32_t>(BinaryClassID::LayeredWorld));
    }

    void read_binary(std::istream& in, bool force = false) override {
        const c_WorldCallLock call_lock(this->p_call_mutex.get());
        try {
            c_TidalPyBaseClass::read_binary(in, force);
            this->read_world_fields(in);
            if (!in) {
                throw std::runtime_error("TidalPy: failed to read LayeredWorld binary data");
            }
            this->read_tide_section(in, force);
            uint64_t n_layers = 0;
            in.read(reinterpret_cast<char*>(&n_layers), sizeof(uint64_t));
            c_SpinConfig spin_config;
            in.read(reinterpret_cast<char*>(&spin_config.moment_of_inertia_factor), sizeof(double));
            if (!in) {
                throw std::runtime_error("TidalPy: failed to read LayeredWorld binary data");
            }
            try {
                this->p_spin = c_Spin(spin_config);
            } catch (const std::invalid_argument& error) {
                throw std::runtime_error(std::string("TidalPy: corrupt LayeredWorld binary data: ") + error.what());
            }
            this->read_solver_overrides(in);
            check_binary_count(in, n_layers, TIDALPY_BINARY_HEADER_BYTES, "layer");
            // Read every layer before replacing the old ones, so a corrupt record leaves the stack whole.
            std::vector<std::unique_ptr<c_BaseLayer>> loaded_layers;
            loaded_layers.reserve(n_layers);
            for (uint64_t i = 0; i < n_layers; ++i) {
                loaded_layers.push_back(c_layer_from_binary(in, force));
                loaded_layers.back()->set_owner_call_mutex(this->p_call_mutex.get());
            }
            this->p_layers = std::move(loaded_layers);
        } catch (...) {
            // The world fields may already hold the new record's values, so nothing solved describes this world.
            this->p_reset_solved_state();
            throw;
        }
        // Nothing solved describes the loaded layers, and the masses floating layers held belong to the old ones.
        this->p_reference_mass.clear();
        this->p_reset_solved_state();
        this->p_warm_start_central_pressure = TidalPyConstants::d_NAN;
    }

protected:
    // Shared so a subclass reuses the layout with its own BinaryClassID.
    void write_layered_binary(std::ostream& out, uint32_t class_id) const {
        std::ostringstream overrides_stream;
        this->write_solver_overrides(overrides_stream);
        const std::string overrides_bytes = overrides_stream.str();
        const uint64_t payload = this->world_payload_bytes() + tide_section_payload_bytes()
            + sizeof(uint64_t) + sizeof(double) + static_cast<uint64_t>(overrides_bytes.size());
        write_binary_header(out, class_id, payload);
        this->write_world_fields(out);
        this->write_tide_section(out);
        const auto n_layers = static_cast<uint64_t>(this->p_layers.size());
        out.write(reinterpret_cast<const char*>(&n_layers), sizeof(uint64_t));
        const double moi_factor = this->p_spin.get_config().moment_of_inertia_factor;
        out.write(reinterpret_cast<const char*>(&moi_factor), sizeof(double));
        out.write(overrides_bytes.data(), static_cast<std::streamsize>(overrides_bytes.size()));
        if (!out) {
            throw std::runtime_error("TidalPy: failed to write layered-world binary data");
        }
        // Each layer writes its own complete record, its models included.
        for (const auto& layer : this->p_layers) { layer->write_binary(out); }
    }

    // One pinned solver key: a presence flag, then the value as the fixed-width Stored type when set.
    template <typename Stored, typename Value>
    static void write_optional_setting(std::ostream& out, const std::optional<Value>& setting) {
        const uint8_t present = setting.has_value() ? 1 : 0;
        out.write(reinterpret_cast<const char*>(&present), sizeof(uint8_t));
        if (setting.has_value()) {
            const Stored stored = static_cast<Stored>(*setting);
            out.write(reinterpret_cast<const char*>(&stored), sizeof(Stored));
        }
    }

    template <typename Stored, typename Value>
    static void read_optional_setting(std::istream& in, std::optional<Value>& setting) {
        uint8_t present = 0;
        in.read(reinterpret_cast<char*>(&present), sizeof(uint8_t));
        if (present) {
            Stored stored {};
            in.read(reinterpret_cast<char*>(&stored), sizeof(Stored));
            setting = static_cast<Value>(stored);
        } else {
            setting.reset();
        }
    }

    // The [eos_solver] then the [radial_solver] keys, in declaration order: ODE methods as int32_t, counts as
    // uint64_t, switches as uint8_t.
    void write_solver_overrides(std::ostream& out) const {
        const c_EOSSolverOverrides& eos = this->p_eos_solver_overrides;
        write_optional_setting<int32_t>(out, eos.integration_method);
        write_optional_setting<double>(out, eos.rtol);
        write_optional_setting<double>(out, eos.atol);
        write_optional_setting<double>(out, eos.pressure_tol);
        write_optional_setting<uint64_t>(out, eos.max_iters);
        write_optional_setting<uint64_t>(out, eos.slices_per_layer);
        write_optional_setting<uint8_t>(out, eos.nondimensionalize);
        write_optional_setting<uint8_t>(out, eos.solve_temperature);
        const c_RadialSolverOverrides& radial = this->p_radial_solver_overrides;
        write_optional_setting<int32_t>(out, radial.integration_method);
        write_optional_setting<double>(out, radial.rtol);
        write_optional_setting<double>(out, radial.atol);
        write_optional_setting<uint8_t>(out, radial.use_kamata);
        write_optional_setting<double>(out, radial.start_radius_tol);
        write_optional_setting<uint8_t>(out, radial.scale_rtols);
        write_optional_setting<uint64_t>(out, radial.max_num_steps);
        write_optional_setting<uint64_t>(out, radial.expected_size);
        write_optional_setting<uint64_t>(out, radial.max_ram_MB);
        write_optional_setting<uint8_t>(out, radial.nondimensionalize);
    }

    // Reads into locals and commits only when the whole block was read.
    void read_solver_overrides(std::istream& in) {
        c_EOSSolverOverrides eos;
        read_optional_setting<int32_t>(in, eos.integration_method);
        read_optional_setting<double>(in, eos.rtol);
        read_optional_setting<double>(in, eos.atol);
        read_optional_setting<double>(in, eos.pressure_tol);
        read_optional_setting<uint64_t>(in, eos.max_iters);
        read_optional_setting<uint64_t>(in, eos.slices_per_layer);
        read_optional_setting<uint8_t>(in, eos.nondimensionalize);
        read_optional_setting<uint8_t>(in, eos.solve_temperature);
        c_RadialSolverOverrides radial;
        read_optional_setting<int32_t>(in, radial.integration_method);
        read_optional_setting<double>(in, radial.rtol);
        read_optional_setting<double>(in, radial.atol);
        read_optional_setting<uint8_t>(in, radial.use_kamata);
        read_optional_setting<double>(in, radial.start_radius_tol);
        read_optional_setting<uint8_t>(in, radial.scale_rtols);
        read_optional_setting<uint64_t>(in, radial.max_num_steps);
        read_optional_setting<uint64_t>(in, radial.expected_size);
        read_optional_setting<uint64_t>(in, radial.max_ram_MB);
        read_optional_setting<uint8_t>(in, radial.nondimensionalize);
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read the LayeredWorld solver settings binary data");
        }
        this->p_eos_solver_overrides    = eos;
        this->p_radial_solver_overrides = radial;
    }

    // Steps every floating layer toward the mass it holds and carries the layers above it, returning the
    // largest relative radius change, which is what the solve watches to stop.
    //
    // A layer that holds its mass moves its top by the mass it is short of over the slope of the enclosed
    // mass there, dm/dr = 4 pi r^2 rho. Its base has already moved with the layer below, which changes the
    // mass between its faces by that base shell, so both ends enter the step. A layer that holds its volume
    // keeps it, so its outer radius follows from its new base.
    double update_floating_radii(const c_EOSSolution& solution) {
        const std::size_t n_layers = this->p_layers.size();
        double largest_change   = 0.0;
        double radius_inner_new = 0.0;
        for (std::size_t layer_i = 0; layer_i < n_layers; ++layer_i) {
            c_BaseLayer* layer = this->p_layers[layer_i].get();
            const double radius_inner = layer->get_radius_inner();
            const double radius_outer = layer->get_radius_outer();
            double radius_outer_new   = radius_outer;

            if (layer->get_is_volume_fixed()) {
                const double volume_term = radius_outer * radius_outer * radius_outer
                    - radius_inner * radius_inner * radius_inner;
                radius_outer_new = std::cbrt(
                    radius_inner_new * radius_inner_new * radius_inner_new + volume_term);
            } else {
                // Step in enclosed volume, not radius: the mass the layer is short of occupies
                // mass_missing / rho of it. At constant density that lands the radius in one pass however
                // far away it starts, where a step in radius would crawl as the cube root.
                double state[C_EOS_DY_VALUES];
                solution.call_si(layer_i, radius_outer, state);
                const double mass_outer    = state[C_EOS_MASS_INDEX];
                const double density_outer = state[C_EOS_DENSITY_INDEX];
                solution.call_si(layer_i, radius_inner, state);
                const double mass_inner = state[C_EOS_MASS_INDEX];
                // A layer whose reference mass is still unset adopts what it holds inside the boundaries it
                // was given, so it stays put until something else moves.
                if (!std::isfinite(this->p_reference_mass[layer_i])) {
                    this->p_reference_mass[layer_i] = mass_outer - mass_inner;
                }
                const double target_mass = this->p_reference_mass[layer_i];
                const double volume_term = radius_outer * radius_outer * radius_outer
                    - radius_inner * radius_inner * radius_inner;
                double radius_cubed = radius_inner_new * radius_inner_new * radius_inner_new + volume_term;
                if (std::isfinite(target_mass) && (density_outer > TidalPyConstants::d_EPS)) {
                    const double mass_missing = target_mass - (mass_outer - mass_inner);
                    radius_cubed += 3.0 * mass_missing / (4.0 * TidalPyConstants::d_PI * density_outer);
                }
                if (!(radius_cubed > 0.0)) {
                    throw std::runtime_error(
                        "TidalPy: a layer holding its mass shrank past its own base during the EOS solve. "
                        "Check the mass it holds against the density its material EOS model gives.");
                }
                radius_outer_new = std::cbrt(radius_cubed);
            }

            const double scale  = (radius_outer > TidalPyConstants::d_EPS) ? radius_outer : 1.0;
            const double change = std::fabs(radius_outer_new - radius_outer) / scale;
            if (change > largest_change) { largest_change = change; }
            layer->set_radii(radius_inner_new, radius_outer_new);
            radius_inner_new = radius_outer_new;
        }
        // The outermost layer's top is the world radius.
        this->p_radius = radius_inner_new;
        return largest_change;
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
        const c_BaseLayer* layer = this->find_layer_for_radius(radius);
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
            const c_BaseLayer* layer = this->p_layers[layer_i].get();
            if (radius >= layer->get_radius_inner() && radius <= layer->get_radius_outer()) {
                this->p_eos_solution->call_si(layer_i, radius, state_out);
                return;
            }
        }
    }

    // The shear (is_shear) or bulk complex modulus [Pa] at a radius [m] and frequency [rad s-1] from the rheology
    // of the layer that holds the radius; NaN for a geometry-only layer. The caller holds the call lock, which is also
    // the layer's.
    std::complex<double> p_complex_modulus(bool is_shear, double radius, double frequency) const noexcept {
        const auto* physics_layer = dynamic_cast<const c_PhysicsLayer*>(this->find_layer_for_radius(radius));
        if (physics_layer == nullptr) { return std::complex<double>(TidalPyConstants::d_NAN, 0.0); }
        return physics_layer->p_complex_modulus(is_shear, radius, frequency);
    }

    // The layer whose radial span contains radius [m], non-owning. Radii beyond the surface clamp to the
    // outermost layer, radii below the innermost inner radius fall in the innermost one, and the result is
    // null only for a world with no layers. Public because it is a const observer query like
    // get_density(radius); the 3D tidal-heating path reads the layer's solid and incompressible flags.
public:
    const c_BaseLayer* find_layer_for_radius(double radius) const noexcept {
        if (this->p_layers.empty()) { return nullptr; }
        for (const auto& layer : this->p_layers) {
            if (radius <= layer->get_radius_outer()) { return layer.get(); }
        }
        return this->p_layers.back().get();
    }

protected:
    std::vector<std::unique_ptr<c_BaseLayer>> p_layers;

    // Set by solve_eos and not serialized, so a load repopulates them by re-solving.
    bool        p_eos_success          = false;
    bool        p_eos_solved           = false;
    std::string p_eos_message          = "EOS not yet solved.";
    int         p_eos_iterations       = -1;
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
    // The mass each floating layer holds on to [kg]; NaN for a layer that holds its volume instead.
    std::vector<double> p_reference_mass;
    bool p_geometry_converged  = true;
    size_t p_thermal_passes    = 0;
    bool   p_thermal_converged = true;

    // The world's own Love solve (solve_love_numbers): the cached radial solver, rebuilt only when the EOS grid or
    // the layer assumptions change, or the quasi-homogeneous results. Not serialized.
    c_LoveWorkspace p_love;
    // Taken by the calls c_WorldCallLock lists; held through a pointer so the world stays movable.
    std::unique_ptr<std::recursive_mutex> p_call_mutex = std::make_unique<std::recursive_mutex>();
    // Set only inside calc_layer_tidal_heating_radial, under the call lock (get_retained_radial_solves).
    const std::vector<c_RetainedRadialSolve>* p_retained_radial_solves = nullptr;
    // The stretches of the last successful EOS solve (update_radial_segments).
    std::vector<c_RadialSegment> p_radial_segments;

    // Global (1D) tidal dissipation: the model/config/result state lives on c_BaseWorld;
    // c_LayeredWorld adds only the per-layer heating distribution (results not serialized).
    std::vector<double>                  p_layer_tidal_heating;
};

} // namespace tidalpy
