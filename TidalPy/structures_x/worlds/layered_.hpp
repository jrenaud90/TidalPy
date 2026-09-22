#pragma once
/*
 * layered_.hpp: c_LayeredWorld, a world built from an ordered stack of layers (extends c_BaseWorld).
 *
 * Owns its layers as std::unique_ptr<c_BaseLayer> (inner to outer, index 0 = innermost) and provides the
 * whole-planet aggregates (total mass, internal radiogenic heating) and geometry validation. The whole-planet
 * EOS and radial (Love number) solves, which walk every layer, are methods on this class.
 *
 * Binary format (20-byte header + payload):
 *   header: class_id = BinaryClassID::LayeredWorld (201)
 *   payload: [all c_BaseWorld fields] + layer_count (uint64_t)
 *   Each layer's own complete binary record (with its recursively-serialized
 *   physics sub-models) follows, in index order, as separate appended records.
 */

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <istream>
#include <limits>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "base_.hpp"
#include "../layers/factory_.hpp"

#include "constants_.hpp"   // TidalPyConstants::d_PI, tidalpy_config_ptr
#include "../../dynamics_x/spin_.hpp"   // c_Spin (spin-dynamics model attached to the world)
#include "solver_.hpp"      // c_solve_eos, c_EOS_ODEInput, c_EOSSolution, ODEMethod, PreEvalFunc
#include "material_.hpp"    // c_MaterialEOSInput, c_preeval_material_eos
#include "thermal_layout_.hpp"  // c_LayerThermal and the segment layout of a thermal solve

// RadialSolver sub-modules (shooting solver, storage, love numbers) are compiled into this translation unit so
// the shooting CyRK integration runs in the extension that owns the CySolverResult objects; no cross-extension
// call().
#include "../../Utilities_x/math_x/numerics_.hpp"        // c_isclose
#include "../../Utilities_x/dimensions/nondimensional_.hpp"  // c_NonDimensionalScales
#include "../../utilities/arrays/interp_.hpp"        // c_interp_complex
#include "../../RadialSolver_x/rs_constants_.hpp"
#include "../../RadialSolver_x/rs_solution_.hpp"
#include "../../RadialSolver_x/love_.hpp"
#include "../../RadialSolver_x/shooting_.hpp"
#include "../../RadialSolver_x/world_radial_solver_.hpp"  // c_WorldRadialSolver (cached Love-number solver)

// Global (1D) tidal dissipation: the model hierarchy and the stored config and result structs are light (no
// global-potential tables). The calc_tides orchestration that runs the heavy global-potential engine and
// collapse is defined out-of-line in world_tides_.hpp, so the eccentricity and obliquity tables compile into
// only the one extension that includes it.
#include "../../Tides_x/classes/tide_base_.hpp"     // tidalpy::c_TideBase
#include "../../Tides_x/classes/tide_result_.hpp"   // c_TideConfig, c_TideSolveConfig, c_GlobalTideResult
#include "../../Tides_x/love/love_method_.hpp"    // c_LoveMethod, homogeneous-sphere Love numbers
// Relative paths (not bare names) so every extension that includes layered_.hpp resolves
// these without needing Utilities_x/lookups on its include path. Light headers, no tables.
#include "../../Utilities_x/lookups/keys_.hpp"     // c_Key4, c_Key2
#include "../../Utilities_x/lookups/intmap_.hpp"   // c_IntMap (per-mode solver Love-number store)

namespace tidalpy {

// Solver settings from the shared runtime config: every EOS and Love-number solve starts from the [eos_solver]
// and [radial_solver] sections of TidalPy_Configs_x.toml (pushed into the TidalPyConfig singleton by
// update_constants_x), and a caller overrides only the fields it passes. The member initializers of the two
// structs stand in only when the config has not been loaded, which does not happen after a normal import.
// CyRK ODEMethod from the int stored in the shared config; the fallback covers an unloaded config.
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

// Parameters for the whole-planet EOS solve, grouped so c_LayeredWorld::solve_eos takes a single argument.
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
    // Temperature [K] for every layer, overriding the temperature each carries. NaN uses each layer's own.
    double    temperature         = TidalPyConstants::d_NAN;
    // Carry temperature and heat flow through the solve. Without a temperature contrast there is no profile
    // to integrate, so the solve keeps its four structure variables whatever this says.
    bool      solve_temperature   = true;
    // Temperature [K] the outermost layer radiates to. NaN leaves no flow through the surface.
    double    surface_temperature = TidalPyConstants::d_NAN;
    // Time [s] the heat sources are evaluated at, on the clock the radiogenics models share. NaN takes each
    // model's own reference time.
    double    time                = TidalPyConstants::d_NAN;
    // Cap on the passes that relax the boundary layers and interface flows against the structure, and the
    // relative change in the interface temperatures and flows that ends them.
    size_t    max_thermal_passes  = 12;
    double    thermal_tol         = 1.0e-8;
    // Relative change in a floating layer's radius that ends the passes, and whether this solve redefines
    // the mass every floating layer holds (otherwise the first solve of the world's life sets it).
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
    // Surface boundary conditions to solve for, in order: 1 = tidal, 2 = loading, 0 = free.
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
    size_t    expected_size      = 1000;
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

    // Set the boundary conditions from a plain array, for the Cython wrappers.
    void set_bc_models(const int* models_ptr, size_t num_models) {
        if (models_ptr == nullptr || num_models == 0) { return; }
        this->bc_models.assign(models_ptr, models_ptr + num_models);
    }
};

class c_LayeredWorld : public c_BaseWorld {
public:
    // Absolute tolerance on the gap between a layer's inner radius and the previous layer's outer
    // radius, from config_x [numerical].layer_continuity_rtol scaled by that radius.
    static double layer_continuity_tol(double previous_outer_radius) noexcept {
        const double scale = (previous_outer_radius > 1.0) ? previous_outer_radius : 1.0;
        return tidalpy_config_ptr->d_LAYER_CONTINUITY_RTOL * scale;
    }

    // Construction
    c_LayeredWorld() = default;

    explicit c_LayeredWorld(const c_WorldConfig& cfg) : c_BaseWorld(cfg) {}

    ~c_LayeredWorld() override = default;

    // Layer ownership
    // Add a layer, inner to outer. Validates that its inner radius matches the
    // current outermost radius (0 for the first layer) within the relative
    // tolerance. Throws std::invalid_argument on a gap/overlap.
    void add_layer(std::unique_ptr<c_BaseLayer> layer) {
        if (!layer) {
            throw std::invalid_argument("TidalPy: cannot add a null layer to a world");
        }
        const double prev_outer = this->p_layers.empty() ? 0.0
                                : this->p_layers.back()->get_radius_outer();
        const double inner      = layer->get_radius_inner();
        const double tol        = layer_continuity_tol(prev_outer);
        if (std::abs(inner - prev_outer) > tol) {
            throw std::invalid_argument(
                "TidalPy: layer geometry is not continuous. Inner radius does not "
                "match the previous layer's outer radius (add layers inner-to-outer)");
        }
        this->p_layers.push_back(std::move(layer));
    }

    // True if `layer` would continue the stack: its inner radius matches the
    // current outermost radius (0 for the first layer) within tolerance. Lets a
    // caller check continuity before transferring ownership via add_layer.
    bool accepts_layer(const c_BaseLayer& layer) const noexcept {
        const double prev_outer = this->p_layers.empty() ? 0.0
                                : this->p_layers.back()->get_radius_outer();
        const double tol        = layer_continuity_tol(prev_outer);
        return std::abs(layer.get_radius_inner() - prev_outer) <= tol;
    }

    // Non-owning observer pointer to the layer at index (throws out_of_range).
    c_BaseLayer* get_layer(std::size_t index) const {
        if (index >= this->p_layers.size()) {
            throw std::out_of_range("TidalPy: layer index out of range");
        }
        return this->p_layers[index].get();
    }

    std::size_t get_num_layers() const noexcept { return this->p_layers.size(); }

    // Whole-planet EOS profile queries (const, MKS). Once the world EOS solve has populated each layer, these
    // return the radially interpolated density, gravity, and pressure at radius r [m] from the layer that
    // contains r (clamped at the surface). NaN when no layer contains r or the EOS has not been solved.
    double get_density(double radius) const noexcept {
        const c_BaseLayer* layer = this->find_layer_for_radius(radius);
        return (layer != nullptr) ? layer->get_density(radius)
                                   : std::numeric_limits<double>::quiet_NaN();
    }

    double get_gravity(double radius) const noexcept {
        const c_BaseLayer* layer = this->find_layer_for_radius(radius);
        return (layer != nullptr) ? layer->get_gravity(radius)
                                   : std::numeric_limits<double>::quiet_NaN();
    }

    double get_pressure(double radius) const noexcept {
        const c_BaseLayer* layer = this->find_layer_for_radius(radius);
        return (layer != nullptr) ? layer->get_pressure(radius)
                                   : std::numeric_limits<double>::quiet_NaN();
    }

    // Temperature [K] at a radius [m] from the solved profile, and the heat flowing outward through the
    // sphere of that radius [W]. A solve with no temperature contrast reports each layer's own temperature
    // and no flow. NaN before the EOS is solved or outside the world.
    double get_temperature(double radius) const noexcept {
        return this->p_read_eos_state(radius, C_EOS_TEMPERATURE_INDEX);
    }

    double get_heat_flow(double radius) const noexcept {
        return this->p_read_eos_state(radius, C_EOS_HEAT_FLOW_INDEX);
    }

    // Per-layer thermal results of the last solve: temperatures, heat flows, boundary layers, and the
    // Rayleigh and Nusselt numbers of a convecting layer.
    const std::vector<c_LayerThermal>& get_layer_thermal() const noexcept { return this->p_layer_thermal; }
    size_t get_thermal_passes()    const noexcept { return this->p_thermal_passes; }
    bool   get_geometry_converged() const noexcept { return this->p_geometry_converged; }
    bool   get_thermal_converged() const noexcept { return this->p_thermal_converged; }

    // The world's heat sources, as the last thermal solve prepared them.
    const c_Heating& get_heating() const noexcept { return this->p_heating; }

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

    // Viscoelastic profile queries (post-melt, with pre-melt variants). Each delegates to the layer containing
    // radius. NaN when no layer contains r, the layer is geometry-only, or the EOS has not been solved.
    double get_shear_modulus(double radius) const noexcept {
        const c_BaseLayer* layer = this->find_layer_for_radius(radius);
        return (layer != nullptr) ? layer->get_shear_modulus(radius) : TidalPyConstants::d_NAN;
    }
    double get_bulk_modulus(double radius) const noexcept {
        const c_BaseLayer* layer = this->find_layer_for_radius(radius);
        return (layer != nullptr) ? layer->get_bulk_modulus(radius) : TidalPyConstants::d_NAN;
    }
    double get_shear_viscosity(double radius) const noexcept {
        const c_BaseLayer* layer = this->find_layer_for_radius(radius);
        return (layer != nullptr) ? layer->get_shear_viscosity(radius) : TidalPyConstants::d_NAN;
    }
    double get_bulk_viscosity(double radius) const noexcept {
        const c_BaseLayer* layer = this->find_layer_for_radius(radius);
        return (layer != nullptr) ? layer->get_bulk_viscosity(radius) : TidalPyConstants::d_NAN;
    }
    double get_melt_fraction(double radius) const noexcept {
        const c_BaseLayer* layer = this->find_layer_for_radius(radius);
        return (layer != nullptr) ? layer->get_melt_fraction(radius) : TidalPyConstants::d_NAN;
    }

    // Radius-resolved complex moduli [Pa] at a frequency, the only per-frequency step: find the layer and apply
    // its rheology to the stored post-melt static modulus and viscosity. NaN for a geometry-only layer or
    // before the solve.
    std::complex<double> calc_complex_shear_modulus(double radius, double frequency) const noexcept {
        const auto* physics_layer = dynamic_cast<const c_PhysicsLayer*>(this->find_layer_for_radius(radius));
        if (physics_layer == nullptr) { return std::complex<double>(TidalPyConstants::d_NAN, 0.0); }
        return physics_layer->calc_complex_shear_modulus(radius, frequency);
    }
    std::complex<double> calc_complex_bulk_modulus(double radius, double frequency) const noexcept {
        const auto* physics_layer = dynamic_cast<const c_PhysicsLayer*>(this->find_layer_for_radius(radius));
        if (physics_layer == nullptr) { return std::complex<double>(TidalPyConstants::d_NAN, 0.0); }
        return physics_layer->calc_complex_bulk_modulus(radius, frequency);
    }

    // True once the EOS solve has populated the innermost layer's profile.
    bool get_eos_solved() const noexcept {
        return !this->p_layers.empty() && this->p_layers.front()->get_eos_data_populated();
    }

    // True once every layer has a material EOS model attached (a precondition
    // for the world-level EOS solve).
    bool get_all_eos_set() const noexcept {
        if (this->p_layers.empty()) { return false; }
        for (const auto& layer : this->p_layers) {
            if (!layer->get_eos_set()) { return false; }
        }
        return true;
    }

    // Whole-planet equation-of-state solve (populates the layer profiles). Integrates the radial structure
    // (gravity, pressure, enclosed mass, moment of inertia) from center to surface with each layer's attached
    // material EOS model as the local density source, through the Material_x/eos c_solve_eos machinery and the
    // c_preeval_material_eos pre-eval, and populates every layer's c_LayerEOSData and mass on success. All
    // MKS. Throws std::invalid_argument when the world has no layers, a layer lacks an EOS model, or
    // slices_per_layer < 2.
    //
    // Assumptions
    // -----------
    // - Spherical symmetry.
    // - Each layer's density comes from its material EOS model (pressure for the analytic models, radius for
    //   the interpolated model).
    void solve_eos(const c_WorldEOSSolveConfig& cfg) {
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

        // Build a per-layer ascending radius grid (concatenated) and the bulk
        // density guess (volume-weighted EOS density at the surface pressure).
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

        // Unit scales of the solve. The default non-dimensional solve integrates in the length, density, and time
        // units of c_NonDimensionalScales (the radius, the bulk density, and 1/sqrt(pi G rho)), so the tolerances
        // mean the same thing for every planet and the central pressure is of order one. An SI solve keeps every
        // scale at one.
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

        // The per-layer material-EOS inputs live in a world member, not a local. The solver keeps a copy of
        // each c_EOS_ODEInput, and that copy holds a pointer (eos_input_ptr) into this vector through which
        // every later evaluation of the density and moduli reaches the layer's EOS model, so it must outlive
        // solve_eos. The vector lives as long as the world and is re-set on every solve.
        // From here on the solve replaces the structure, so nothing solved on top of the old one may be read.
        this->mark_structure_dirty();
        this->p_eos_material_inputs.assign(n_layers, c_MaterialEOSInput());

        c_EOS_ODEInput ode_input;
        ode_input.G_to_use      = G_solve;
        ode_input.planet_radius = upper_radii.back();
        ode_input.update_bulk   = false;
        ode_input.update_shear  = false;
        for (std::size_t i = 0; i < n_layers; ++i) {
            this->p_eos_material_inputs[i].eos_model_ptr = this->p_layers[i]->get_eos();
            this->p_eos_material_inputs[i].length_scale  = length_scale;
            this->p_eos_material_inputs[i].pascal_scale  = pascal_scale;
            this->p_eos_material_inputs[i].density_scale = density_scale;
            ode_input.eos_input_ptr = reinterpret_cast<char*>(&this->p_eos_material_inputs[i]);
            eos_function_vec.push_back(c_preeval_material_eos);
            eos_input_vec.push_back(ode_input);
        }

        // Thermal layout. Each layer carries its own temperature and its cooling model says how heat moves
        // inside it; a uniform override replaces every layer's value.
        c_init_layer_thermal(this->p_layers, this->p_layer_thermal);
        if (std::isfinite(cfg.temperature)) {
            for (c_LayerThermal& thermal : this->p_layer_thermal) { thermal.temperature = cfg.temperature; }
        }
        // Heat sources, prepared once: none of them depends on the solved state. They need the heat flow of a
        // thermal solve to act through, so a solve with temperature switched off leaves them out.
        c_WorldState world_state;
        world_state.time       = cfg.time;
        world_state.layers_ptr = &this->p_layers;
        this->p_heating.update_sources(world_state, length_scale, density_scale);
        const bool heating_active = this->p_heating.get_is_active();
        if (heating_active && !cfg.solve_temperature) {
            TIDALPY_LOG_WARN(
                "TidalPy: world '{}' has layers with use_heating set, but solve_temperature is off, so this EOS "
                "solve carries no heat flow and the heating is ignored.", this->get_name());
        }
        // A temperature contrast or a heated layer gives the solve a profile to integrate.
        const bool thermal_contrast = cfg.solve_temperature
            && (c_thermal_contrast_present(this->p_layer_thermal, cfg.surface_temperature) || heating_active);
        for (std::size_t i = 0; i < n_layers; ++i) {
            eos_input_vec[i].heating_ptr = thermal_contrast ? &this->p_heating : nullptr;
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
        // The central pressure the secant iteration starts from, in solve units: the world's last converged solve,
        // then the pass before. A re-solve after a small change then converges in a pass or two. NaN, on a world
        // that has never been solved, starts from a uniform sphere.
        double central_pressure_guess = this->p_eos_solved
            ? (this->p_central_pressure / pascal_scale) : TidalPyConstants::d_NAN;
        const std::size_t last_pass =
            (thermal_contrast || geometry_floats) ? cfg.max_thermal_passes : 0;
        this->p_thermal_passes    = 0;
        this->p_thermal_converged = !thermal_contrast;
        bool geometry_converged   = !geometry_floats;
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
            }
            c_build_thermal_segments(
                this->p_layer_thermal, this->p_layers, integrate_temperature,
                length_scale, gravity_scale, segment_vec);
            for (std::size_t i = 0; i < n_layers; ++i) {
                // The viscosity and melt models of the material always see the temperature; its density law
                // sees it only when the layer asked for a thermal EOS.
                const auto* physics_layer = dynamic_cast<const c_PhysicsLayer*>(this->p_layers[i].get());
                const bool thermal_eos = (physics_layer != nullptr) && physics_layer->get_use_thermal_eos();
                this->p_eos_material_inputs[i].temperature           = this->p_layer_thermal[i].temperature;
                this->p_eos_material_inputs[i].use_state_temperature = integrate_temperature;
                this->p_eos_material_inputs[i].thermal_density       = thermal_eos;
            }

            solution = std::make_shared<c_EOSSolution>(
                upper_radii.data(), n_layers, full_radius.data(), total_slices);
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
            this->p_thermal_passes = pass;
            central_pressure_guess = solution->central_pressure / pascal_scale;

            if (thermal_contrast) {
                const double thermal_change = c_update_layer_thermal(
                    *solution, this->p_layers, cfg.surface_temperature, integrate_temperature,
                    this->p_layer_thermal, &this->p_heating);
                this->p_thermal_converged = integrate_temperature && (thermal_change < cfg.thermal_tol);
            }
            if (geometry_floats) {
                geometry_converged = (this->update_floating_radii(*solution) < cfg.radius_tol);
            }
            if (this->p_thermal_converged && geometry_converged) { break; }
        }
        this->p_geometry_converged = geometry_converged;

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

        // Populate each layer's structure + viscoelastic profile from its slice of
        // the full arrays.
        if (this->p_eos_solved) {
            const std::size_t total_slices = solution->radius_array_size;
            for (std::size_t layer_index = 0; layer_index < n_layers; ++layer_index) {
                const std::size_t slice_start = layer_index * slices;
                const std::size_t slice_end   = slice_start + slices;
                if (slice_end > total_slices) { break; }
                c_BaseLayer* layer = this->p_layers[layer_index].get();

                // The layer's mass is the enclosed mass gained across its slices. Adjacent layers share the
                // interface slice, so the layer masses sum exactly to the planet mass.
                layer->set_mass(
                    solution->mass_array_vec[slice_end - 1] - solution->mass_array_vec[slice_start]);

                c_LayerEOSData eos_data;
                eos_data.populate(
                    std::vector<double>(solution->radius_array_vec.begin()   + slice_start, solution->radius_array_vec.begin()   + slice_end),
                    std::vector<double>(solution->density_array_vec.begin()  + slice_start, solution->density_array_vec.begin()  + slice_end),
                    std::vector<double>(solution->gravity_array_vec.begin()  + slice_start, solution->gravity_array_vec.begin()  + slice_end),
                    std::vector<double>(solution->pressure_array_vec.begin() + slice_start, solution->pressure_array_vec.begin() + slice_end));

                // Install the dense-output evaluator for this layer: the solution's SI-radius call, which
                // evaluates the layer's retained CySolverResult (owning its solver) and the layer's EOS function
                // for the density and moduli. The captured shared_ptr co-owns the whole solution, so the dense
                // data and solver stay alive and callable post-solve, and the EOS arguments point at
                // this->p_eos_material_inputs, which also outlives the solve. The lambda is compiled in this
                // (CyRK-owning) extension, so the CySolverResult is only ever called by the CyRK copy that built
                // it. The slice arrays populated above are the fallback for a manual update_eos_data.
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

        // Retain the full solution so callers can read the radial profile arrays.
        this->p_eos_solution = std::move(solution);
    }

    // Everything solved on top of the structure describes the structure it was solved with: the Love numbers,
    // the global tides, the heating handed to each layer, and the cached radial-solver setup (a re-solve changes
    // the structure and moduli even when the grid size does not). solve_eos calls this before it replaces the
    // structure, so none of them can be read against the new one; each comes back with its own next solve.
    void mark_structure_dirty() noexcept {
        this->p_love_solved           = false;
        this->p_love_analytic_success = false;
        this->p_tides_solved          = false;
        this->p_tide_solver_love.clear();
        this->p_layer_tidal_heating.clear();
        for (const auto& layer_uptr : this->p_layers) {
            layer_uptr->set_tidal_heating(TidalPyConstants::d_NAN);
        }
        if (this->p_radial_solver) { this->p_radial_solver->invalidate(); }
    }

    // EOS solve result accessors (valid after solve_eos; NaN/empty otherwise).
    bool               get_eos_success()          const noexcept { return this->p_eos_success; }
    const std::string& get_eos_message()          const noexcept { return this->p_eos_message; }
    int                get_eos_iterations()       const noexcept { return this->p_eos_iterations; }
    // True when the central-pressure iteration stopped at max_iters without meeting pressure_tol; the solution
    // is still populated from the last iteration.
    bool               get_eos_max_iters_hit()    const noexcept { return this->p_eos_max_iters_hit; }
    double             get_eos_pressure_error()   const noexcept { return this->p_eos_pressure_error; }
    double             get_surface_gravity_eos()  const noexcept { return this->p_surface_gravity_eos; }
    double             get_surface_pressure_eos() const noexcept { return this->p_surface_pressure_eos; }
    double             get_central_pressure()     const noexcept { return this->p_central_pressure; }
    double             get_planet_mass_eos()      const noexcept { return this->p_planet_mass_eos; }
    double             get_planet_moi_eos()       const noexcept { return this->p_planet_moi_eos; }

    // Spin dynamics (rates only). The world drives its attached c_Spin model with its own EOS-based moment of
    // inertia, so the spin-rate change uses the structure-resolved value rather than the uniform-density one.
    void           set_spin_model(const c_Spin& spin) noexcept { this->p_spin = spin; }
    const c_Spin&  get_spin_model() const noexcept { return this->p_spin; }

    // Moment of inertia [kg m2]: the EOS-solved value (get_planet_moi_eos) when the EOS has been solved,
    // otherwise the spin model's factor * M R^2 estimate from the world mass and radius.
    double get_moment_of_inertia() const noexcept {
        if (this->p_eos_solved && std::isfinite(this->p_planet_moi_eos)) {
            return this->p_planet_moi_eos;
        }
        return this->p_spin.calc_moment_of_inertia(this->get_mass(), this->get_radius());
    }

    // Tidal spin-rate change [rad s-2] = M_host * dU/dO / I, using the world's stored dU/dO (from the
    // last calc_tides) and its moment of inertia. Requires a completed tidal solve.
    double calc_spin_derivative(double host_mass) const {
        if (!this->get_tides_solved()) {
            throw std::runtime_error(
                "TidalPy: spin derivative needs a tidal solve first: call calc_tides()");
        }
        return this->p_spin.calc_dspin_dt(host_mass, this->get_tidal_dU_dO(), this->get_moment_of_inertia());
    }

    // Synchronous spin rate [rad s-1]: equal to the orbital mean motion.
    double calc_synchronous_spin(double orbital_frequency) const noexcept {
        return this->p_spin.calc_synchronous_spin(orbital_frequency);
    }

    // Non-owning observer pointer to the retained full-planet EOS solution (the
    // source of the radial profile arrays), or nullptr if solve_eos was not run.
    const c_EOSSolution* get_eos_solution() const noexcept { return this->p_eos_solution.get(); }

    // Whole-planet Love-number solve (non-const; throws std::invalid_argument without a prior solve_eos). Uses
    // the world's EOS arrays and the layers' rheology models to build the frequency-dependent complex moduli,
    // then calls the shooting solver directly. The world class is the complete interface;
    // c_RadialSolutionStorage is an internal detail owned by the cached radial solver. Spherical symmetry,
    // all MKS.
    // Ensure the cached radial-solver setup matches the current EOS/config; (re)build it if not. Returns false (and
    // stamps an error on the solver storage) if the per-layer slice partitioning is invalid. Throws on hard
    // precondition failures (no EOS solve, no config, no layers, too few slices).
    bool ensure_radial_cache(const c_LoveSolveConfig& cfg) {
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

        if (!this->p_radial_solver)
            this->p_radial_solver = std::make_unique<::c_WorldRadialSolver>();
        ::c_WorldRadialSolver* solver = this->p_radial_solver.get();

        // Per-layer metadata (geometry-only layers default to static solid). Gathered before
        // the cache check because the flags are user-mutable without an EOS re-solve: a cache
        // hit is only valid when they match what the cached inputs were built from.
        auto layer_types   = std::make_unique<int[]>(n_layers);
        auto is_static_arr = std::make_unique<bool[]>(n_layers);
        auto is_incomp_arr = std::make_unique<bool[]>(n_layers);
        std::vector<double> upper_radii(n_layers);
        for (std::size_t i = 0; i < n_layers; ++i) {
            const auto* phys = dynamic_cast<const c_PhysicsLayer*>(this->p_layers[i].get());
            if (phys != nullptr) {
                layer_types[i]   = phys->get_is_solid() ? 0 : 1;
                is_static_arr[i] = phys->get_is_static();
                is_incomp_arr[i] = phys->get_is_incompressible();
            } else {
                layer_types[i]   = 0;
                is_static_arr[i] = true;
                is_incomp_arr[i] = false;
            }
            upper_radii[i] = this->p_layers[i]->get_radius_outer();
        }

        const std::size_t num_ytypes = cfg.bc_models.empty() ? 1 : cfg.bc_models.size();
        if (solver->cache_matches(n_layers, total_slices, cfg.degree_l, cfg.nondimensionalize, num_ytypes)
            && solver->layer_flags_match(layer_types.get(), is_static_arr.get(), is_incomp_arr.get(), n_layers))
            return true;

        const c_EOSSolution* world_eos = this->p_eos_solution.get();
        const double r_planet = world_eos->radius;
        const double vol      = (4.0 / 3.0) * TidalPyConstants::d_PI * r_planet * r_planet * r_planet;
        const double bulk_rho = (vol > TidalPyConstants::d_EPS) ? this->p_planet_mass_eos / vol : 3500.0;

        return solver->build_cache(
            world_eos->radius_array_vec,
            world_eos->density_array_vec,
            world_eos->gravity_array_vec,
            world_eos->pressure_array_vec,
            world_eos->mass_array_vec,
            world_eos->moi_array_vec,
            upper_radii,
            layer_types.get(),
            is_static_arr.get(),
            is_incomp_arr.get(),
            n_layers,
            r_planet,
            bulk_rho,
            cfg.degree_l,
            cfg.nondimensionalize,
            this->p_eos_solution.get(),
            num_ytypes
        );
    }

    // Build the per-call runtime config from the user-facing solve config.
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
        rt.warnings           = cfg.warnings;
        return rt;
    }

    void solve_love_numbers(const c_LoveSolveConfig& cfg) {
        this->solve_love_numbers(cfg, nullptr);
    }

    // Composite Simpson intervals per tidal layer for the homogeneous volume average (even, 129 nodes).
    static constexpr std::size_t homogeneous_quadrature_intervals = 128;

    // Frequency-independent inputs to the homogeneous volume average for one tidal layer
    struct c_HomogeneousShearLayer {
        const c_PhysicsLayer* layer = nullptr;
        double dr = 0.0;                      // node spacing [m]
        std::vector<double> radius;           // node radii [m]
        std::vector<double> static_modulus;   // post-melt static shear modulus at each node [Pa]
        std::vector<double> viscosity;        // post-melt shear viscosity at each node [Pa s]
    };

    // Reusable state for a run of homogeneous Love solves against one unchanged interior. The node values are
    // read on first use and the averaged modulus is kept per distinct frequency, so a caller solving many
    // (degree, frequency) pairs pays for each frequency once.
    struct c_HomogeneousLoveCache {
        struct c_ShearAtFrequency {
            double frequency;                 // forcing frequency the average was formed at [rad s-1]
            bool use_static;                  // true for the cpl / ctl methods, which average the static modulus
            std::complex<double> shear;       // volume-averaged shear modulus [Pa]
        };
        bool built = false;
        std::vector<c_HomogeneousShearLayer> layers;
        double tidal_volume = 0.0;            // summed volume of the averaged layers [m3]
        std::vector<c_ShearAtFrequency> shear_by_frequency;
    };

    // Same, with reusable state for a caller that solves many frequencies in a row against an unchanged interior.
    // The homogeneous methods use it; the radial-solver methods ignore it.
    void solve_love_numbers(const c_LoveSolveConfig& cfg, c_HomogeneousLoveCache* cache) {
        const c_LoveMethod method = c_love_method_from_int(cfg.love_method);
        this->p_love_method_last = method;
        if (c_love_method_is_homogeneous(method)) {
            this->solve_love_numbers_homogeneous(cfg, method, cache);
            return;
        }
        if (method == c_LoveMethod::LaterallyInhomogeneous) {
            throw std::logic_error(
                "TidalPy: the laterally_inhomogeneous Love-number method is reserved for the 3D Love solver and is "
                "not implemented.");
        }
        if (!this->ensure_radial_cache(cfg)) { this->p_love_solved = false; return; }
        ::c_WorldRadialSolver* solver = this->p_radial_solver.get();

        // Frequency-dependent step: hand the solver a provider that evaluates each layer's material state at the
        // radius the integrator asks for, then solve. Nothing is sampled onto the slice grid, so the Love numbers
        // no longer carry the first-order slice error, and the layer the provider is asked about is the one the
        // solver is integrating: an interface radius belongs to two layers, and a lookup by radius alone would give
        // both copies the lower layer's material. The provider captures this solve's frequency.
        //
        // The provider stays on the storage after the solve, so the solution keeps answering at any radius: the
        // y3 of a dynamic liquid layer is rebuilt from the density and gravity it reads, and an exported solution
        // reports the complex moduli this solve used.
        solver->set_material_eval(this->make_material_eval(cfg.frequency));

        c_LoveSolveRuntimeConfig rt = this->make_runtime_config(cfg);
        solver->solve(rt);
        this->p_love_solved = solver->get_solved();
    }

    // Solve using externally-supplied complex moduli (the standalone array API path) instead of the layer rheology.
    void solve_love_numbers_supplied(
            const c_LoveSolveConfig& cfg,
            const std::complex<double>* shear_in,
            const std::complex<double>* bulk_in,
            const double* radius_in,
            std::size_t n_in) {
        const c_LoveMethod method = c_love_method_from_int(cfg.love_method);
        if (!c_love_method_uses_radial_solver(method)) {
            throw std::invalid_argument(
                "TidalPy: solve_love_numbers_supplied supports only the radial_solver and propagation_matrix "
                "Love-number methods.");
        }
        this->p_love_method_last = method;
        if (!this->ensure_radial_cache(cfg)) { this->p_love_solved = false; return; }
        ::c_WorldRadialSolver* solver = this->p_radial_solver.get();

        // Find, once per layer, the run of supplied points that bounds it. The last point at its base through the
        // first point at its top. The provider below then interpolates inside that run, so a repeated interface
        // radius gives the upper copy to the layer above and the lower copy to the layer below; interpolating
        // across all layers at once would give a lower layer's top the upper layer's value. The bounds match
        // within the relative tolerance the cache uses for interfaces, so copies that differ by rounding still
        // count as the interface.
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
        // The provider owns copies of the supplied profile rather than borrowing the caller's arrays, because it
        // is also what an exported solution reads: it has to outlive this call. The layer runs are captured the
        // same way for the same reason.
        std::vector<double>               radius_copy(radius_in, radius_in + n_in);
        std::vector<std::complex<double>> shear_copy(shear_in, shear_in + n_in);
        std::vector<std::complex<double>> bulk_copy(bulk_in, bulk_in + n_in);
        std::shared_ptr<const c_EOSSolution> eos_solution = this->p_eos_solution;
        solver->set_material_eval(
            [eos_solution, radius_copy, shear_copy, bulk_copy, in_first_by_layer, in_count_by_layer](
                    std::size_t layer_index,
                    double radius_si,
                    double* state_out,
                    std::complex<double>& shear_out,
                    std::complex<double>& bulk_out) {
                // One dense read gives the structure and the density; the supplied profile gives the moduli.
                eos_solution->call_si(layer_index, radius_si, state_out);
                const std::size_t in_first = in_first_by_layer[layer_index];
                const std::size_t in_count = in_count_by_layer[layer_index];
                // The supplied grid is usually uniform within a layer, so the fractional position seeds the search.
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
                // Supplied moduli carry no separate unrelaxed value, so the real part stands in for it, and the
                // profile said nothing about viscosity: it gave the response directly.
                state_out[C_EOS_SHEAR_MODULUS_INDEX]   = shear_out.real();
                state_out[C_EOS_BULK_MODULUS_INDEX]    = bulk_out.real();
                state_out[C_EOS_SHEAR_VISCOSITY_INDEX] = TidalPyConstants::d_NAN;
                state_out[C_EOS_BULK_VISCOSITY_INDEX]  = TidalPyConstants::d_NAN;
            });

        c_LoveSolveRuntimeConfig rt = this->make_runtime_config(cfg);
        rt.redim_eos_arrays = true;
        solver->solve(rt);

        // Report SI EOS scalars on the released storage (the world own EOS solution is dimensional).
        if (solver->get_storage() != nullptr) {
            c_EOSSolution* dst       = solver->get_storage()->get_eos_solution_ptr();
            const c_EOSSolution* src = this->p_eos_solution.get();
            dst->radius           = src->radius;
            dst->mass             = src->mass;
            dst->moi              = src->moi;
            dst->surface_gravity  = src->surface_gravity;
            dst->surface_pressure = src->surface_pressure;
            dst->central_pressure = src->central_pressure;
        }
        this->p_love_solved = solver->get_solved();
    }

    // The per-layer material-state provider the radial solver reads at each integration radius: one dense EOS call
    // for the frequency-independent state, then the layer's rheology, which is the only part that knows the
    // frequency. Fills the evaluation layout of eos_layout_.hpp (SI) and the two complex moduli [Pa].
    //
    // The callable co-owns the solved EOS and the rheologies, and the layers are resolved here, once, rather than at
    // every radius. The solved EOS still reaches each layer's material model through this world, so the world has
    // to outlive the callable; a solution exported to Python holds its world for that reason.
    c_EOSSolution::MaterialEval make_material_eval(double frequency) const {
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
        return [eos_solution, shear_bylayer, bulk_bylayer, is_physics_bylayer, frequency](
                std::size_t layer_index,
                double radius_si,
                double* state_out,
                std::complex<double>& shear_out,
                std::complex<double>& bulk_out) {
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

    // Move the radial-solution storage out of the helper (one-shot export to a RadialSolverSolution).
    std::unique_ptr<::c_RadialSolutionStorage> release_radial_storage() {
        if (!this->p_radial_solver) { return nullptr; }
        return this->p_radial_solver->release_storage();
    }

    // Love-number solve result accessors, valid after solve_love_numbers succeeds and NaN or empty otherwise.
    // The world is the sole interface; c_RadialSolutionStorage is internal.

    // Analytic Love numbers (homogeneous, cpl, ctl): the world is treated as a homogeneous incompressible
    // sphere with the planet's bulk density, EOS surface gravity, and radius, and the volume-averaged shear
    // modulus of the layers flagged is_tidal (composite Simpson rule in radius with the r^2 volume weight,
    // using each layer's radius-resolved moduli and rheology). The homogeneous method averages the complex
    // modulus at the forcing frequency; cpl and ctl average the static (unrelaxed) modulus and then impose the
    // constant phase lag (1 - i/Q) or time lag (1 - i omega dt), taking Q or dt from the solve config or, when
    // unset, from the attached tide model.
    //
    // Assumptions
    // -----------
    // - Incompressible homogeneous-sphere response; layered structure enters only through the volume average.
    // - Gas layers carry no shear modulus and are skipped; liquid layers contribute their (zero) shear modulus.
    // Read the frequency-independent node values of the homogeneous volume average into the cache.
    void build_homogeneous_shear_nodes(c_HomogeneousLoveCache& cache) const {
        const std::size_t n_intervals = homogeneous_quadrature_intervals;
        cache.layers.clear();
        cache.tidal_volume = 0.0;
        for (const auto& layer_ptr : this->p_layers) {
            const c_BaseLayer* layer = layer_ptr.get();
            if (!layer->get_is_tidal()) {
                continue;
            }
            const auto* physics = dynamic_cast<const c_PhysicsLayer*>(layer);
            if (physics == nullptr) {
                continue;   // gas layers: no shear modulus
            }
            const double r_inner = layer->get_radius_inner();
            const double r_outer = layer->get_radius_outer();
            if (!(r_outer > r_inner)) {
                continue;
            }
            c_HomogeneousShearLayer nodes;
            nodes.layer = physics;
            nodes.dr = (r_outer - r_inner) / static_cast<double>(n_intervals);
            nodes.radius.resize(n_intervals + 1);
            nodes.static_modulus.resize(n_intervals + 1);
            nodes.viscosity.resize(n_intervals + 1);
            for (std::size_t i = 0; i <= n_intervals; ++i) {
                const double r = (i == n_intervals) ? r_outer : r_inner + static_cast<double>(i) * nodes.dr;
                nodes.radius[i]         = r;
                nodes.static_modulus[i] = physics->get_shear_modulus(r);    // post-melt
                nodes.viscosity[i]      = physics->get_shear_viscosity(r);  // post-melt
            }
            cache.tidal_volume += layer->get_volume();
            cache.layers.push_back(std::move(nodes));
        }
        cache.built = true;
    }

    // Volume-weighted complex shear modulus over the cached tidal layers at one frequency. The rheology is applied
    // exactly as calc_complex_shear_modulus applies it, in the same summation order, so the result matches a solve
    // that re-reads every node. Returns false and records error -41 on a non-finite modulus.
    bool average_homogeneous_shear(
            const c_HomogeneousLoveCache& cache,
            double frequency,
            bool use_static,
            std::complex<double>& shear_avg) {
        const std::size_t n_intervals = homogeneous_quadrature_intervals;
        std::complex<double> shear_integral(0.0, 0.0);
        for (const c_HomogeneousShearLayer& nodes : cache.layers) {
            const c_RheologyBase* rheology = nodes.layer->get_shear_rheology_model();
            std::complex<double> layer_sum(0.0, 0.0);
            for (std::size_t i = 0; i <= n_intervals; ++i) {
                const double r = nodes.radius[i];
                const double weight = (i == 0 || i == n_intervals) ? 1.0 : ((i % 2 == 1) ? 4.0 : 2.0);
                const std::complex<double> mu = (use_static || rheology == nullptr)
                    ? std::complex<double>(nodes.static_modulus[i], 0.0)
                    : rheology->calc_complex_modulus(nodes.static_modulus[i], nodes.viscosity[i], frequency);
                if (!std::isfinite(mu.real()) || !std::isfinite(mu.imag())) {
                    this->p_love_analytic_error_code = -41;
                    this->p_love_analytic_message =
                        "TidalPy: layer '" + nodes.layer->get_name() + "' returned a non-finite shear modulus at r = "
                        + std::to_string(r) + " m (viscosity or rheology not set?); the homogeneous Love-number "
                        "methods need finite moduli in every tidal layer.";
                    return false;
                }
                layer_sum += weight * mu * (r * r);
            }
            shear_integral += layer_sum * (nodes.dr / 3.0) * (4.0 * TidalPyConstants::d_PI);
        }
        shear_avg = shear_integral / cache.tidal_volume;
        return true;
    }

    void solve_love_numbers_homogeneous(
            const c_LoveSolveConfig& cfg,
            c_LoveMethod method,
            c_HomogeneousLoveCache* cache) {
        this->p_love_solved = false;
        this->p_love_analytic_success = false;
        this->p_love_analytic = c_LoveNumbers();
        this->p_love_analytic_shear = std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        this->p_love_analytic_tidal_volume = TidalPyConstants::d_NAN;
        if (!this->p_eos_solved || !this->p_eos_solution) {
            throw std::invalid_argument("TidalPy: solve_eos() must be called before solve_love_numbers().");
        }
        if (cfg.degree_l < 2) {
            throw std::invalid_argument("TidalPy: the homogeneous Love-number methods need degree_l >= 2.");
        }
        const bool use_static = (method != c_LoveMethod::Homogeneous);

        // Without a caller-owned cache the node values go into a local one that ends with this solve.
        c_HomogeneousLoveCache local_cache;
        c_HomogeneousLoveCache& active = (cache != nullptr) ? *cache : local_cache;
        if (!active.built) {
            this->build_homogeneous_shear_nodes(active);
        }
        if (!(active.tidal_volume > 0.0)) {
            this->p_love_analytic_error_code = -40;
            this->p_love_analytic_message =
                "TidalPy: the homogeneous Love-number methods need at least one tidal layer (is_tidal) with a shear "
                "modulus.";
            return;
        }
        const double tidal_volume = active.tidal_volume;

        // The average depends on the frequency but not on the degree, so it is formed once per distinct frequency.
        std::complex<double> shear_avg(0.0, 0.0);
        bool have_average = false;
        for (const auto& entry : active.shear_by_frequency) {
            if (entry.frequency == cfg.frequency && entry.use_static == use_static) {
                shear_avg = entry.shear;
                have_average = true;
                break;
            }
        }
        if (!have_average) {
            if (!this->average_homogeneous_shear(active, cfg.frequency, use_static, shear_avg)) {
                return;
            }
            active.shear_by_frequency.push_back({cfg.frequency, use_static, shear_avg});
        }

        const double radius = this->get_radius();
        const double density_bulk = this->get_mass() / ((4.0 / 3.0) * TidalPyConstants::d_PI * radius * radius * radius);
        const double gravity = this->p_surface_gravity_eos;
        if (!(gravity > 0.0) || !std::isfinite(gravity)) {
            this->p_love_analytic_error_code = -42;
            this->p_love_analytic_message =
                "TidalPy: the EOS surface gravity is not finite and positive; re-run solve_eos() before the "
                "homogeneous Love-number methods.";
            return;
        }

        c_LoveNumbers love = c_calc_homogeneous_love_numbers(
            shear_avg,
            density_bulk,
            gravity,
            radius,
            cfg.degree_l);
        if (method == c_LoveMethod::HomogeneousCPL) {
            // Precedence: the solve config, then the [tides] config, then the attached tide model.
            double fixed_q = cfg.fixed_q;
            if (!std::isfinite(fixed_q)) {
                fixed_q = this->get_tide_config().love_fixed_q;
            }
            if (!std::isfinite(fixed_q) && this->p_tide) {
                fixed_q = this->p_tide->get_fixed_q(cfg.degree_l);
            }
            if (!(fixed_q > 0.0) || !std::isfinite(fixed_q)) {
                throw std::invalid_argument(
                    "TidalPy: the cpl Love-number method needs a positive fixed_q for degree "
                    + std::to_string(cfg.degree_l) + " (pass fixed_q, set it in the tides config, or attach a tide "
                    "model that carries a fixed Q).");
            }
            love = c_apply_fixed_q(love, fixed_q);
        } else if (method == c_LoveMethod::HomogeneousCTL) {
            double fixed_dt = cfg.fixed_dt;
            if (!std::isfinite(fixed_dt)) {
                fixed_dt = this->get_tide_config().love_fixed_dt;
            }
            if (!std::isfinite(fixed_dt) && this->p_tide) {
                fixed_dt = this->p_tide->get_fixed_dt(cfg.degree_l);
            }
            if (!(fixed_dt >= 0.0) || !std::isfinite(fixed_dt)) {
                throw std::invalid_argument(
                    "TidalPy: the ctl Love-number method needs a non-negative fixed_dt for degree "
                    + std::to_string(cfg.degree_l) + " (pass fixed_dt, set love_fixed_dt_s in the tides config, or "
                    "attach a tide model that carries a fixed time lag).");
            }
            love = c_apply_fixed_dt(love, cfg.frequency, fixed_dt);
        }

        this->p_love_analytic = love;
        this->p_love_analytic_shear = shear_avg;
        this->p_love_analytic_tidal_volume = tidal_volume;
        this->p_love_analytic_success = true;
        this->p_love_analytic_error_code = 0;
        this->p_love_analytic_message = std::string("Homogeneous-sphere Love numbers (") + c_love_method_name(method) + ").";
        this->p_love_solved = true;
    }

    // Love-solve config carrying the world's configured method and its cpl / ctl parameters (from the [tides]
    // config); the tide paths start from this so the configured method drives every Love-number solve.
    c_LoveSolveConfig make_love_solve_config() const {
        c_LoveSolveConfig cfg;
        const c_TideConfig& tide_cfg = this->get_tide_config();
        cfg.love_method = tide_cfg.love_method;
        cfg.fixed_q     = tide_cfg.love_fixed_q;
        cfg.fixed_dt    = tide_cfg.love_fixed_dt;
        return cfg;
    }

    // Same, for paths that need the depth-resolved radial solution (3D stress/strain/heating): the analytic methods
    // have no radial y-functions, so they are rejected with an explanatory error.
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

    bool love_is_analytic() const noexcept { return c_love_method_is_homogeneous(this->p_love_method_last); }
    int  get_love_method_last_int() const noexcept { return static_cast<int>(this->p_love_method_last); }
    // Diagnostics of the last analytic solve: the volume-averaged shear modulus [Pa] and the averaged volume [m3]
    // (NaN after a radial-solver solve).
    std::complex<double> get_love_analytic_shear() const noexcept { return this->p_love_analytic_shear; }
    double get_love_analytic_tidal_volume() const noexcept { return this->p_love_analytic_tidal_volume; }

    // Non-owning pointer to the internal solution storage (owned by the cached
    // radial solver). Null until solve_love_numbers has built the cache.
    const ::c_RadialSolutionStorage* get_love_storage() const noexcept {
        return this->p_radial_solver ? this->p_radial_solver->get_storage() : nullptr;
    }

    bool get_love_solved() const noexcept { return this->p_love_solved; }
    bool get_love_success() const noexcept {
        if (this->love_is_analytic()) return this->p_love_analytic_success;
        const auto* s = this->get_love_storage();
        return (s && this->p_love_solved) ? s->success : false;
    }
    int get_love_error_code() const noexcept {
        if (this->love_is_analytic()) return this->p_love_analytic_error_code;
        const auto* s = this->get_love_storage();
        return s ? s->error_code : -100;
    }
    const std::string& get_love_message() const noexcept {
        static const std::string no_msg = "No love-number solve has been run.";
        if (this->love_is_analytic()) return this->p_love_analytic_message;
        const auto* s = this->get_love_storage();
        return s ? s->message : no_msg;
    }
    std::size_t get_love_num_ytypes() const noexcept {
        if (this->love_is_analytic()) return this->p_love_analytic_success ? 1 : 0;
        const auto* s = this->get_love_storage();
        return s ? s->num_ytypes : 0;
    }
    // Worst-case error amplification of the surface boundary condition solve (shooting method; 0 until a
    // solve has run, and 0 for the analytic methods). See c_estimate_surface_amplification in
    // RadialSolver_x/boundaries/boundaries_.hpp.
    double get_love_surface_amplification() const noexcept {
        if (this->love_is_analytic()) return 0.0;
        const auto* s = this->get_love_storage();
        return s ? s->surface_amplification : 0.0;
    }
    // Primary Love numbers (k, h, l) for the given boundary-condition ytype index (the analytic methods hold a
    // single tidal set at index 0). NaN when no solve describes the current structure: never solved, failed, or
    // followed by a solve_eos.
    std::complex<double> get_love_number_k(std::size_t ytype_idx = 0) const noexcept {
        if (this->love_is_analytic()) {
            return (this->p_love_analytic_success && ytype_idx == 0)
                ? this->p_love_analytic.k : std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        }
        const auto* s = this->get_love_storage();
        if (!s || !this->p_love_solved || ytype_idx >= s->complex_love_vec.size())
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        return s->complex_love_vec[ytype_idx].k;
    }
    std::complex<double> get_love_number_h(std::size_t ytype_idx = 0) const noexcept {
        if (this->love_is_analytic()) {
            return (this->p_love_analytic_success && ytype_idx == 0)
                ? this->p_love_analytic.h : std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        }
        const auto* s = this->get_love_storage();
        if (!s || !this->p_love_solved || ytype_idx >= s->complex_love_vec.size())
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        return s->complex_love_vec[ytype_idx].h;
    }
    std::complex<double> get_love_number_l(std::size_t ytype_idx = 0) const noexcept {
        if (this->love_is_analytic()) {
            return (this->p_love_analytic_success && ytype_idx == 0)
                ? this->p_love_analytic.l : std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        }
        const auto* s = this->get_love_storage();
        if (!s || !this->p_love_solved || ytype_idx >= s->complex_love_vec.size())
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        return s->complex_love_vec[ytype_idx].l;
    }
    // Full radial solution y-value (SI) at the surface for a given ytype and y-index (0..5 -> y1..y6).
    // Returns NaN if not solved or after an analytic solve (no radial functions). Evaluated from the radial
    // solver's dense calling system (shooting) or the gridded solution (matrix) via get_surface_y.
    std::complex<double> get_love_surface_y(
            std::size_t ytype_idx, std::size_t y_idx) const noexcept {
        if (this->love_is_analytic()) return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        const auto* s = this->get_love_storage();
        if (!s || !this->p_love_solved || !s->success || y_idx >= C_MAX_NUM_Y)
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        std::complex<double> surface_y[C_MAX_NUM_Y];
        if (!s->get_surface_y(ytype_idx, surface_y))
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        return surface_y[y_idx];
    }
    // Full radial solution y-value (SI) at an arbitrary radius [m] for a given ytype and y-index. The shooting
    // method evaluates its dense per-layer interpolants at this radius (accurate anywhere, including between EOS
    // grid slices); the matrix method linearly interpolates its constructed grid. Returns NaN if not solved,
    // after an analytic solve, out of range, or below the solver's starting radius.
    std::complex<double> get_radial_solution_y(
            double radius,
            std::size_t ytype_idx,
            std::size_t y_idx) const noexcept {
        if (this->love_is_analytic()) return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        const auto* s = this->get_love_storage();
        if (!s || !this->p_love_solved || !s->success || y_idx >= C_MAX_NUM_Y)
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        std::complex<double> y_at_r[C_MAX_NUM_Y];
        if (!s->get_radial_solution(radius, ytype_idx, y_at_r))
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        return y_at_r[y_idx];
    }

    // Global (1D) tidal dissipation. The tide-model holder, config, results, and the analytic calc_tides path
    // live on c_BaseWorld; c_LayeredWorld hides calc_tides to add the rheology (radial-solver) path and the
    // per-layer heating distribution. It is defined out-of-line in world_tides_.hpp, which carries the heavy
    // global-potential engine.
    void calc_tides(const c_TideSolveConfig& state);

    // On-demand 3D tidal stress, strain, and heating. The tidal potential is built from the world's [tides]
    // truncation config (max degree l, eccentricity and obliquity truncation); there is no potential-model
    // object. The orchestration lives on the rheology tide model (c_RheologyTide) and is defined out-of-line
    // in world_tides_.hpp, which carries the kernel and potential-engine headers.
    //
    // Secular (cycle and orbit-averaged) 3D volumetric heating [W m-3] at (radius, colatitude): the physically
    // time-averaged power density, independent of longitude and time. Requires the rheology tide model and a
    // solved EOS. Its volume integral equals the 1D global heating (get_tidal_heating).
    double get_3d_tidal_heating(
            const c_TideSolveConfig& state,
            double radius,
            double colatitude);

    // Batch form: secular 3D volumetric heating [W m-3] at num_points paired (radii[i], colatitudes[i]),
    // written into out_heating[i]. Same physics and preconditions as the scalar form, but the radial solve is
    // amortized across points (one per unique (l, frequency)), so it is the efficient way to build a map.
    // num_threads applies to the per-point evaluation after the radial solves. Defined in world_tides_.hpp.
    void get_3d_tidal_heating_array(
            const c_TideSolveConfig& state,
            const double* radii,
            const double* colatitudes,
            size_t num_points,
            double* out_heating,
            int num_threads = 1);

    // Instantaneous tidal displacements [m] on the (radius, colatitude, longitude, time) grid; see
    // c_RheologyTide::calc_3d_displacements_grid. Same rheology + solved-EOS preconditions as
    // get_3d_tidal_heating. Defined out-of-line in world_tides_.hpp.
    void get_3d_displacements_grid(
            const c_TideSolveConfig& state,
            const c_Grid3DAxes& axes,
            double* out_disp,
            int num_threads = 1);

    // Instantaneous stress [Pa] and strain on the (radius, colatitude, longitude, time) grid; see
    // c_RheologyTide::calc_3d_stress_strain_grid. Same rheology + solved-EOS preconditions as
    // get_3d_tidal_heating. Defined out-of-line in world_tides_.hpp.
    void get_3d_stress_strain_grid(
            const c_TideSolveConfig& state,
            const c_Grid3DAxes& axes,
            double* out_stress,
            double* out_strain,
            int num_threads = 1);

    // Collapsed (summed/averaged) secular 3D tidal heating (see c_Heating3DCollapseConfig): the radial
    // power profile, colatitude profile, per-layer totals, and/or whole-planet total, per the flags.
    // Same rheology + solved-EOS preconditions as get_3d_tidal_heating. Defined out-of-line in
    // world_tides_.hpp.
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

    // The axes and output shape calc_3d_tides produces for these inputs, with values and layer_totals left empty,
    // so a caller can allocate the buffers calc_3d_tides_into fills. Needs only the layer geometry. Defined
    // out-of-line in world_tides_.hpp.
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

    // calc_3d_tides written into caller buffers: out_values holds as many doubles as the layout's shape, and
    // out_layer_totals n_layers * n_times doubles when all three spatial axes are summed (null otherwise). The
    // radial solves run on the calling thread and the per-point evaluation on up to cfg.num_threads threads; the
    // result is identical for any thread count. Defined out-of-line in world_tides_.hpp.
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

    // Effective per-layer tidal-heating scale for the layer's tidal_scale_method (defined in
    // world_tides_.hpp). Used by calc_tides to distribute the global heating to the layers.
    double effective_tidal_scale(
            const c_BaseLayer* layer, double planet_volume, const c_TideSolveConfig& state) const;

    // Per-layer tidal heating [W] = world heating × the layer's effective tidal scale (0 if the
    // layer is non-tidal). NaN if the index is out of range or tides are unsolved.
    double get_layer_tidal_heating(std::size_t index) const noexcept {
        if (!this->p_tides_solved || index >= this->p_layer_tidal_heating.size()) {
            return TidalPyConstants::d_NAN;
        }
        return this->p_layer_tidal_heating[index];
    }

    // Whole-planet aggregates (const, MKS)
    // Total mass [kg] = sum of layer masses.
    double calc_total_mass() const noexcept {
        double total = 0.0;
        for (const auto& layer : this->p_layers) { total += layer->get_mass(); }
        return total;
    }

    // Total internal radiogenic heating [W] at time. Only SolidLiquidLayers
    // produce radiogenic heating; other layer types contribute zero.
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

    // Returns true if every layer boundary is continuous (inner-to-outer) and the
    // innermost layer starts at radius 0.
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

    // Binary I/O: the world fields and layer count, then each layer's full record.
    void write_binary(std::ostream& out) const override {
        this->write_layered_binary(out, static_cast<uint32_t>(BinaryClassID::LayeredWorld));
    }

    void read_binary(std::istream& in, bool force = false) override {
        c_TidalPyBaseClass::read_binary(in, force);
        this->read_world_fields(in);
        uint64_t n_layers = 0;
        in.read(reinterpret_cast<char*>(&n_layers), sizeof(uint64_t));
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read LayeredWorld binary data");
        }
        this->p_layers.clear();
        this->p_layers.reserve(n_layers);
        for (uint64_t i = 0; i < n_layers; ++i) {
            this->p_layers.push_back(c_layer_from_binary(in, force));
        }
    }

protected:
    // Shared writer so subclasses (e.g. c_GasGiantWorld) reuse the layout with
    // their own BinaryClassID.
    void write_layered_binary(std::ostream& out, uint32_t class_id) const {
        const uint64_t payload = this->world_payload_bytes() + sizeof(uint64_t);
        write_binary_header(out, class_id, payload);
        this->write_world_fields(out);
        const auto n_layers = static_cast<uint64_t>(this->p_layers.size());
        out.write(reinterpret_cast<const char*>(&n_layers), sizeof(uint64_t));
        if (!out) {
            throw std::runtime_error("TidalPy: failed to write layered-world binary data");
        }
        // Each layer writes its own complete record (recursively including its models).
        for (const auto& layer : this->p_layers) { layer->write_binary(out); }
    }

    // Compute and store a layer's frequency-independent viscoelastic state over its radial slice: the pre-melt
    // static moduli (from the layer's static values) and pre-melt viscosities (from the attached viscosity
    // models, NaN if unset), then the post-melt versions (the partial-melt model applied to the shear pair and
    // then the bulk pair; post equals pre without a melt model). No-op for a geometry-only BaseLayer.
    // temperature is the placeholder profile temperature.
    // Step every floating layer toward the mass it holds, and carry the layers above it. Returns the largest
    // relative radius change, which is what the solve watches to stop.
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
                // mass_missing / rho of it, and the shell keeps the rest. At constant density that lands the
                // radius in one pass however far away it starts, where a step in radius would crawl in
                // proportion to the cube root.
                double state[C_EOS_DY_VALUES];
                solution.call_si(layer_i, radius_outer, state);
                const double mass_outer    = state[2];
                const double density_outer = state[4];
                solution.call_si(layer_i, radius_inner, state);
                const double mass_inner = state[2];
                // A layer whose reference mass is still unset adopts what it holds inside the boundaries it was
                // given, so it stays where it is until something else moves.
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

    // One value of the evaluation layout at a radius [m], through the layer that holds it.
    double p_read_eos_state(double radius, std::size_t value_index) const noexcept {
        if (!this->p_eos_solved || !this->p_eos_solution) { return TidalPyConstants::d_NAN; }
        for (std::size_t layer_i = 0; layer_i < this->p_layers.size(); ++layer_i) {
            const c_BaseLayer* layer = this->p_layers[layer_i].get();
            if (radius >= layer->get_radius_inner() && radius <= layer->get_radius_outer()) {
                double state[C_EOS_DY_VALUES];
                this->p_eos_solution->call_si(layer_i, radius, state);
                return state[value_index];
            }
        }
        return TidalPyConstants::d_NAN;
    }

    // Non-owning observer pointer to the layer whose radial span contains radius [m]. Radii beyond the surface
    // clamp to the outermost layer, radii below the innermost inner radius fall in the innermost layer, and
    // the result is nullptr only when the world has no layers. Public because it is a const observer query
    // like get_density(radius); the 3D tidal-heating path reads the layer's solid and incompressible flags.
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

    // EOS solve results (set by solve_eos; not serialized, so repopulate by re-solving).
    bool        p_eos_success          = false;
    bool        p_eos_solved           = false;
    std::string p_eos_message          = "EOS not yet solved.";
    int         p_eos_iterations       = -1;
    bool        p_eos_max_iters_hit    = false;
    double      p_eos_pressure_error   = std::numeric_limits<double>::quiet_NaN();
    double      p_surface_gravity_eos  = std::numeric_limits<double>::quiet_NaN();
    double      p_surface_pressure_eos = std::numeric_limits<double>::quiet_NaN();
    double      p_central_pressure     = std::numeric_limits<double>::quiet_NaN();
    double      p_planet_mass_eos      = std::numeric_limits<double>::quiet_NaN();
    double      p_planet_moi_eos       = std::numeric_limits<double>::quiet_NaN();
    c_Spin      p_spin {};              // spin-dynamics model (uses the world's EOS moment of inertia)
    std::shared_ptr<c_EOSSolution> p_eos_solution;  // retained full-planet solution (co-owned by layer dense evaluators)
    // Per-layer material-EOS inputs referenced by the solver's stored diffeq args;
    // must outlive solve_eos so the dense output's diffeq re-calls stay valid.
    std::vector<c_MaterialEOSInput> p_eos_material_inputs;
    // Thermal description of every layer from the last solve, and how the thermal passes ended.
    std::vector<c_LayerThermal> p_layer_thermal;
    // The world's heat sources. The structure ODE of a thermal solve reads them through a pointer to this member.
    c_Heating p_heating;
    // The mass each floating layer holds on to [kg]; NaN for a layer that holds its volume instead.
    std::vector<double> p_reference_mass;
    bool p_geometry_converged  = true;
    size_t p_thermal_passes    = 0;
    bool   p_thermal_converged = true;

    // Love-number solve: a cached, reusable radial solver (not serialized) holding the frequency-independent
    // setup and the reused solution storage; rebuilt only when the EOS grid or the layer assumptions change.
    bool p_love_solved = false;
    c_LoveMethod  p_love_method_last = c_LoveMethod::RadialSolver;   // method of the most recent Love solve
    // Result store for the analytic (homogeneous / cpl / ctl) methods; the radial-solver methods report from the
    // cached solver's storage instead.
    bool          p_love_analytic_success = false;
    int           p_love_analytic_error_code = -100;
    std::string   p_love_analytic_message = "No love-number solve has been run.";
    c_LoveNumbers p_love_analytic;
    std::complex<double> p_love_analytic_shear = {TidalPyConstants::d_NAN, 0.0};   // volume-averaged shear [Pa]
    double        p_love_analytic_tidal_volume = TidalPyConstants::d_NAN;          // averaged volume [m3]
    std::unique_ptr<::c_WorldRadialSolver> p_radial_solver;

    // Global (1D) tidal dissipation: the model/config/result state lives on c_BaseWorld;
    // c_LayeredWorld adds only the per-layer heating distribution (results not serialized).
    std::vector<double>                  p_layer_tidal_heating;
};

} // namespace tidalpy
