#pragma once
/*
 * layered_.hpp — c_LayeredWorld: a world built from an ordered stack of layers.
 *
 * Inherits c_BaseWorld. Owns its layers as std::unique_ptr<c_BaseLayer> (inner to
 * outer, index 0 = innermost) and provides whole-planet aggregates (total mass,
 * internal radiogenic heating) and geometry validation. The whole-planet EOS and
 * radial solves (which walk all layers) are added as methods on this class in
 * later phases.
 *
 * Binary format (20-byte header + payload):
 *   header: class_id = BinaryClassID::LayeredWorld (201)
 *   payload: [all c_BaseWorld fields] + layer_count (uint64_t)
 *   Each layer's own complete binary record (with its recursively-serialized
 *   physics sub-models) follows, in index order, as separate appended records.
 */

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
#include "solver_.hpp"      // c_solve_eos, c_EOS_ODEInput, c_EOSSolution, ODEMethod, PreEvalFunc
#include "material_.hpp"    // c_MaterialEOSInput, c_preeval_material_eos

// RadialSolver sub-modules: shooting solver, storage, love numbers.
// Compiled into this TU so the shooting CyRK integration runs in the same
// extension that owns the CySolverResult objects — no cross-extension call().
#include "../../utilities/math/numerics_.hpp"        // c_isclose
#include "../../utilities/dimensions/nondimensional_.hpp"  // c_NonDimensionalScales
#include "../../RadialSolver_x/rs_constants_.hpp"
#include "../../RadialSolver_x/rs_solution_.hpp"
#include "../../RadialSolver_x/love_.hpp"
#include "../../RadialSolver_x/shooting_.hpp"

namespace tidalpy {

// -------------------------------------------------------------------------------
// c_WorldEOSSolveConfig — parameters for the whole-planet EOS solve.
// Using a config struct keeps c_LayeredWorld::solve_eos to a single argument.
// -------------------------------------------------------------------------------
struct c_WorldEOSSolveConfig {
    double    surface_pressure    = 0.0;                 // [Pa]
    size_t    slices_per_layer    = 100;                 // radial samples per layer (>= 2)
    double    G_to_use            = -1.0;                // [m^3 kg^-1 s^-2]; < 0 -> TidalPy config G
    ODEMethod integration_method  = ODEMethod::DOP853;
    double    rtol                = 1.0e-6;
    double    atol                = 1.0e-10;
    double    pressure_tol        = 1.0e-3;
    size_t    max_iters           = 100;
    double    temperature         = 0.0;                 // [K]; passed to calc_density (unused yet)
    bool      verbose             = false;
};

// -------------------------------------------------------------------------------
// c_LoveSolveConfig — parameters for the whole-planet Love-number solve.
// -------------------------------------------------------------------------------
struct c_LoveSolveConfig {
    double    frequency_rad_s    = 1.0e-5;            // [rad/s]; tidal forcing frequency
    int       degree_l           = 2;                  // harmonic degree
    bool      solve_tidal        = true;               // include tidal (bc_model=1) boundary condition
    bool      use_kamata         = true;
    bool      nondimensionalize  = true;
    double    starting_radius    = 0.0;                // [m]; 0 → auto
    double    start_radius_tol   = 1.0e-4;
    ODEMethod integration_method = ODEMethod::DOP853;
    double    rtol               = 1.0e-6;
    double    atol               = 1.0e-10;
    bool      scale_rtols        = true;
    size_t    max_num_steps      = 500000;
    size_t    expected_size      = 500;
    size_t    max_ram_MB         = 500;
    double    max_step           = 0.0;
    bool      verbose            = false;
    bool      warnings           = true;
    double    eos_rtol           = 1.0e-6;
    double    eos_atol           = 1.0e-10;
    double    eos_pressure_tol   = 1.0e-3;
    int       eos_max_iters      = 100;
};

// -------------------------------------------------------------------------------
// c_LayeredWorld
// -------------------------------------------------------------------------------
class c_LayeredWorld : public c_BaseWorld {
public:
    // Relative tolerance used when checking layer-boundary continuity.
    static constexpr double d_LAYER_CONTINUITY_RTOL = 1.0e-6;

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------
    c_LayeredWorld() = default;

    explicit c_LayeredWorld(const c_WorldConfig& cfg) : c_BaseWorld(cfg) {}

    ~c_LayeredWorld() override = default;

    // -----------------------------------------------------------------------
    // Layer ownership
    // -----------------------------------------------------------------------
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
        const double tol        = d_LAYER_CONTINUITY_RTOL
                                * (prev_outer > 1.0 ? prev_outer : 1.0);
        if (std::abs(inner - prev_outer) > tol) {
            throw std::invalid_argument(
                "TidalPy: layer geometry is not continuous — inner radius does not "
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
        const double tol        = d_LAYER_CONTINUITY_RTOL
                                * (prev_outer > 1.0 ? prev_outer : 1.0);
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

    // -----------------------------------------------------------------------
    // Whole-planet EOS profile queries (const, MKS)
    //
    // After the world-level EOS solve has populated each layer's EOS data, these
    // return the radially-interpolated density, gravity, and pressure at radius
    // r [m] by delegating to the layer that contains r (clamped at the surface).
    // Return NaN when no layer contains r or the EOS has not been solved.
    // -----------------------------------------------------------------------
    double get_density(double radius_m) const noexcept {
        const c_BaseLayer* layer = this->find_layer_for_radius(radius_m);
        return (layer != nullptr) ? layer->get_density(radius_m)
                                   : std::numeric_limits<double>::quiet_NaN();
    }

    double get_gravity(double radius_m) const noexcept {
        const c_BaseLayer* layer = this->find_layer_for_radius(radius_m);
        return (layer != nullptr) ? layer->get_gravity(radius_m)
                                   : std::numeric_limits<double>::quiet_NaN();
    }

    double get_pressure(double radius_m) const noexcept {
        const c_BaseLayer* layer = this->find_layer_for_radius(radius_m);
        return (layer != nullptr) ? layer->get_pressure(radius_m)
                                   : std::numeric_limits<double>::quiet_NaN();
    }

    // -----------------------------------------------------------------------
    // Viscoelastic profile queries (post-melt by default; pre-melt variants too).
    // Each finds the layer containing radius_m and delegates. NaN when no layer
    // contains r, the layer is geometry-only, or the EOS has not been solved.
    // -----------------------------------------------------------------------
    double get_shear_modulus(double radius_m) const noexcept {
        const c_BaseLayer* layer = this->find_layer_for_radius(radius_m);
        return (layer != nullptr) ? layer->get_shear_modulus(radius_m) : TidalPyConstants::d_NAN;
    }
    double get_bulk_modulus(double radius_m) const noexcept {
        const c_BaseLayer* layer = this->find_layer_for_radius(radius_m);
        return (layer != nullptr) ? layer->get_bulk_modulus(radius_m) : TidalPyConstants::d_NAN;
    }
    double get_shear_viscosity(double radius_m) const noexcept {
        const c_BaseLayer* layer = this->find_layer_for_radius(radius_m);
        return (layer != nullptr) ? layer->get_shear_viscosity(radius_m) : TidalPyConstants::d_NAN;
    }
    double get_bulk_viscosity(double radius_m) const noexcept {
        const c_BaseLayer* layer = this->find_layer_for_radius(radius_m);
        return (layer != nullptr) ? layer->get_bulk_viscosity(radius_m) : TidalPyConstants::d_NAN;
    }
    double get_premelt_shear_modulus(double radius_m) const noexcept {
        const c_BaseLayer* layer = this->find_layer_for_radius(radius_m);
        return (layer != nullptr) ? layer->get_premelt_shear_modulus(radius_m) : TidalPyConstants::d_NAN;
    }
    double get_premelt_bulk_modulus(double radius_m) const noexcept {
        const c_BaseLayer* layer = this->find_layer_for_radius(radius_m);
        return (layer != nullptr) ? layer->get_premelt_bulk_modulus(radius_m) : TidalPyConstants::d_NAN;
    }
    double get_premelt_shear_viscosity(double radius_m) const noexcept {
        const c_BaseLayer* layer = this->find_layer_for_radius(radius_m);
        return (layer != nullptr) ? layer->get_premelt_shear_viscosity(radius_m) : TidalPyConstants::d_NAN;
    }
    double get_premelt_bulk_viscosity(double radius_m) const noexcept {
        const c_BaseLayer* layer = this->find_layer_for_radius(radius_m);
        return (layer != nullptr) ? layer->get_premelt_bulk_viscosity(radius_m) : TidalPyConstants::d_NAN;
    }

    // -----------------------------------------------------------------------
    // Radius-resolved complex moduli [Pa] at frequency_rad_s (the only per-ω step):
    // find the layer, apply its rheology to the stored post-melt static modulus +
    // viscosity. NaN+0i for a geometry-only layer or before the solve.
    // -----------------------------------------------------------------------
    std::complex<double> calc_complex_shear_modulus(double radius_m, double frequency_rad_s) const noexcept {
        const auto* physics_layer = dynamic_cast<const c_PhysicsLayer*>(this->find_layer_for_radius(radius_m));
        if (physics_layer == nullptr) { return std::complex<double>(TidalPyConstants::d_NAN, 0.0); }
        return physics_layer->calc_complex_shear_modulus(radius_m, frequency_rad_s);
    }
    std::complex<double> calc_complex_bulk_modulus(double radius_m, double frequency_rad_s) const noexcept {
        const auto* physics_layer = dynamic_cast<const c_PhysicsLayer*>(this->find_layer_for_radius(radius_m));
        if (physics_layer == nullptr) { return std::complex<double>(TidalPyConstants::d_NAN, 0.0); }
        return physics_layer->calc_complex_bulk_modulus(radius_m, frequency_rad_s);
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

    // -----------------------------------------------------------------------
    // Whole-planet equation-of-state solve (non-const; populates layer profiles)
    //
    // Integrates the planet's radial structure (gravity, pressure, enclosed mass,
    // moment of inertia) from center to surface, using each layer's attached
    // material EOS model as the local density source, and populates every layer's
    // c_LayerEOSData on success. Reuses the Material_x/eos c_solve_eos machinery
    // via the c_preeval_material_eos pre-eval. All quantities MKS.
    //
    // Throws std::invalid_argument if the world has no layers, any layer lacks an
    // EOS model, or slices_per_layer < 2.
    //
    // Assumptions
    // -----------
    // - Spherical symmetry.
    // - Each layer's density comes from its material EOS model (pressure for the
    //   analytic models, radius for the interpolated model).
    // -----------------------------------------------------------------------
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
            G_to_use = (tidalpy_config_ptr != nullptr) ? tidalpy_config_ptr->d_G : 6.674015e-11;
        }

        const std::size_t slices       = cfg.slices_per_layer;
        const std::size_t total_slices = n_layers * slices;

        // Build a per-layer ascending radius grid (concatenated) and the bulk
        // density guess (volume-weighted EOS density at the surface pressure).
        std::vector<double> full_radius(total_slices);
        std::vector<double> upper_radii(n_layers);
        double total_volume  = 0.0;
        double mass_estimate = 0.0;
        for (std::size_t i = 0; i < n_layers; ++i) {
            c_BaseLayer* layer   = this->p_layers[i].get();
            const double r_inner = layer->get_radius_inner();
            const double r_outer = layer->get_radius_outer();
            upper_radii[i]       = r_outer;
            for (std::size_t s = 0; s < slices; ++s) {
                const double frac = static_cast<double>(s) / static_cast<double>(slices - 1);
                full_radius[i * slices + s] = r_inner + frac * (r_outer - r_inner);
            }
            const double r_mid     = 0.5 * (r_inner + r_outer);
            const double rho_mid   = layer->get_eos()->calc_density(cfg.surface_pressure, cfg.temperature, r_mid);
            const double shell_vol = (4.0 / 3.0) * TidalPyConstants::d_PI
                                   * (r_outer * r_outer * r_outer - r_inner * r_inner * r_inner);
            total_volume  += shell_vol;
            mass_estimate += rho_mid * shell_vol;
        }
        const double planet_bulk_density = (total_volume > 0.0) ? (mass_estimate / total_volume) : 3500.0;

        // Build the EOS solution object and the per-layer pre-eval functions/inputs.
        auto solution = std::make_shared<c_EOSSolution>(
            upper_radii.data(), n_layers, full_radius.data(), total_slices);

        std::vector<PreEvalFunc>    eos_function_vec;
        std::vector<c_EOS_ODEInput> eos_input_vec;
        eos_function_vec.reserve(n_layers);
        eos_input_vec.reserve(n_layers);

        // The per-layer material-EOS inputs are stored as a WORLD MEMBER (not a
        // local). The solver copies the c_EOS_ODEInput, but that copy holds a
        // pointer (eos_input_ptr) into this vector, and the dense-output re-calls
        // the diffeq (for the density extra output) post-solve through that pointer
        // — so it must outlive solve_eos. The retained solution co-owns the dense
        // evaluators; this vector lives as long as the world (re-set on every solve).
        this->p_eos_material_inputs.assign(n_layers, c_MaterialEOSInput());

        c_EOS_ODEInput ode_input;
        ode_input.G_to_use      = G_to_use;
        ode_input.planet_radius = upper_radii.back();
        ode_input.final_solve   = false;
        ode_input.update_bulk   = false;
        ode_input.update_shear  = false;
        for (std::size_t i = 0; i < n_layers; ++i) {
            this->p_eos_material_inputs[i].eos_model_ptr = this->p_layers[i]->get_eos();
            this->p_eos_material_inputs[i].temperature_k = cfg.temperature;
            ode_input.eos_input_ptr = reinterpret_cast<char*>(&this->p_eos_material_inputs[i]);
            eos_function_vec.push_back(c_preeval_material_eos);
            eos_input_vec.push_back(ode_input);
        }

        // Solve (CyRK ODE integration with the surface-pressure loop).
        c_solve_eos(
            solution.get(),
            eos_function_vec,
            eos_input_vec,
            planet_bulk_density,
            cfg.surface_pressure,
            G_to_use,
            cfg.integration_method,
            cfg.rtol,
            cfg.atol,
            cfg.pressure_tol,
            cfg.max_iters,
            cfg.verbose
        );

        // Store scalar results.
        this->p_eos_success          = solution->success;
        this->p_eos_message          = solution->message;
        this->p_eos_iterations       = solution->iterations;
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

                c_LayerEOSData eos_data;
                eos_data.populate(
                    std::vector<double>(solution->radius_array_vec.begin()   + slice_start, solution->radius_array_vec.begin()   + slice_end),
                    std::vector<double>(solution->density_array_vec.begin()  + slice_start, solution->density_array_vec.begin()  + slice_end),
                    std::vector<double>(solution->gravity_array_vec.begin()  + slice_start, solution->gravity_array_vec.begin()  + slice_end),
                    std::vector<double>(solution->pressure_array_vec.begin() + slice_start, solution->pressure_array_vec.begin() + slice_end));

                // Install the CyRK dense-output evaluator for this layer. The
                // CySolverResult owns its solver (solver_uptr), and the captured
                // shared_ptr co-owns the whole solution, so the dense data + solver
                // stay alive and callable post-solve (no dangling / leak). The diffeq
                // args it re-evaluates (for the density extra output) point at
                // this->p_eos_material_inputs, which also outlives the solve. The
                // lambda is compiled in this (CyRK-owning) extension, so the CySolverResult
                // is only ever called by the CyRK copy that built it. The slice arrays
                // populated above are the fallback for a manual update_eos_data.
                if (layer_index < solution->cysolver_results_uptr_bylayer_vec.size()) {
                    CySolverResult* layer_solver_result =
                        solution->cysolver_results_uptr_bylayer_vec[layer_index].get();
                    if (layer_solver_result != nullptr) {
                        std::shared_ptr<c_EOSSolution> solution_owner = solution;
                        eos_data.set_dense_eval(
                            [solution_owner, layer_solver_result](double radius_m, double* y_out) {
                                layer_solver_result->call(radius_m, y_out);
                            });
                    }
                }
                layer->update_eos_data(eos_data);

                // Frequency-independent viscoelastic post-pass (PhysicsLayers only).
                this->populate_layer_viscoelastic(
                    layer, *solution, slice_start, slices, cfg.temperature);
            }
        }

        // Retain the full solution so callers can read the radial profile arrays.
        this->p_eos_solution = std::move(solution);
    }

    // EOS solve result accessors (valid after solve_eos; NaN/empty otherwise).
    bool               get_eos_success()          const noexcept { return this->p_eos_success; }
    const std::string& get_eos_message()          const noexcept { return this->p_eos_message; }
    int                get_eos_iterations()       const noexcept { return this->p_eos_iterations; }
    double             get_eos_pressure_error()   const noexcept { return this->p_eos_pressure_error; }
    double             get_surface_gravity_eos()  const noexcept { return this->p_surface_gravity_eos; }
    double             get_surface_pressure_eos() const noexcept { return this->p_surface_pressure_eos; }
    double             get_central_pressure()     const noexcept { return this->p_central_pressure; }
    double             get_planet_mass_eos()      const noexcept { return this->p_planet_mass_eos; }
    double             get_planet_moi_eos()       const noexcept { return this->p_planet_moi_eos; }

    // Non-owning observer pointer to the retained full-planet EOS solution (the
    // source of the radial profile arrays), or nullptr if solve_eos was not run.
    const c_EOSSolution* get_eos_solution() const noexcept { return this->p_eos_solution.get(); }

    // -----------------------------------------------------------------------
    // Whole-planet Love-number solve (non-const; requires a prior solve_eos).
    //
    // Uses the world's already-computed EOS (density, gravity, pressure arrays)
    // and the layers' rheology models to build frequency-dependent complex moduli,
    // then calls the shooting solver directly — no external radial-solver function
    // needed.  The world class is the complete interface; c_RadialSolutionStorage
    // is an internal detail owned by p_love_solution.
    //
    // Throws std::invalid_argument if solve_eos has not been run.
    //
    // Assumptions: spherical symmetry; all quantities MKS.
    // -----------------------------------------------------------------------
    void solve_love_numbers(const c_LoveSolveConfig& cfg) {
        if (!this->p_eos_solved || !this->p_eos_solution) {
            throw std::invalid_argument(
                "TidalPy: solve_eos must succeed before solve_love_numbers");
        }
        if (tidalpy_config_ptr == nullptr) {
            throw std::runtime_error(
                "TidalPy: config not initialised — call initialize_tidalpy_config() first");
        }

        const std::size_t n_layers     = this->p_layers.size();
        const std::size_t total_slices = this->p_eos_solution->radius_array_size;

        if (n_layers == 0)
            throw std::invalid_argument("TidalPy: world has no layers");
        if (total_slices < 5)
            throw std::invalid_argument("TidalPy: EOS solution has too few radial slices (< 5)");

        // ── Per-layer metadata ─────────────────────────────────────────────
        auto layer_types   = std::make_unique<int[]>(n_layers);
        auto is_static_arr = std::make_unique<bool[]>(n_layers);
        auto is_incomp_arr = std::make_unique<bool[]>(n_layers);
        for (std::size_t i = 0; i < n_layers; ++i) {
            const auto* phys = dynamic_cast<const c_PhysicsLayer*>(this->p_layers[i].get());
            if (phys != nullptr) {
                layer_types[i]   = phys->get_is_solid() ? 0 : 1;
                is_static_arr[i] = phys->get_is_static();
                is_incomp_arr[i] = phys->get_is_incompressible();
            } else {
                layer_types[i]   = 0;      // geometry-only layer treated as solid
                is_static_arr[i] = true;
                is_incomp_arr[i] = false;
            }
        }

        // ── Planet-level scalars ───────────────────────────────────────────
        const c_EOSSolution* world_eos = this->p_eos_solution.get();
        const double r_planet    = world_eos->radius;
        const double g_surface   = world_eos->surface_gravity;   // SI value, saved for find_love
        const double vol         = (4.0 / 3.0) * TidalPyConstants::d_PI
                                 * r_planet * r_planet * r_planet;
        const double bulk_rho    = (vol > 0.0) ? this->p_planet_mass_eos / vol : 3500.0;
        const double G_si        = tidalpy_config_ptr->d_G;

        // ── Mutable working copies of all EOS arrays ───────────────────────
        // We non-dimensionalize these in place; the originals in p_eos_solution
        // are never modified.
        std::vector<double> radius_copy  (world_eos->radius_array_vec);
        std::vector<double> density_copy (world_eos->density_array_vec);
        std::vector<double> gravity_copy (world_eos->gravity_array_vec);
        std::vector<double> pressure_copy(world_eos->pressure_array_vec);
        std::vector<double> mass_copy    (world_eos->mass_array_vec);
        std::vector<double> moi_copy     (world_eos->moi_array_vec);

        // Frequency-dependent complex moduli at the forcing frequency.
        std::vector<std::complex<double>> shear_vec(total_slices);
        std::vector<std::complex<double>> bulk_vec(total_slices);
        for (std::size_t i = 0; i < total_slices; ++i) {
            const double r = radius_copy[i];
            shear_vec[i] = this->calc_complex_shear_modulus(r, cfg.frequency_rad_s);
            bulk_vec[i]  = this->calc_complex_bulk_modulus(r, cfg.frequency_rad_s);
        }

        // Upper radii per layer (will be non-dimensionalized below if needed).
        std::vector<double> upper_radii(n_layers);
        for (std::size_t i = 0; i < n_layers; ++i)
            upper_radii[i] = this->p_layers[i]->get_radius_outer();

        // ── Non-dimensionalization ─────────────────────────────────────────
        ::c_NonDimensionalScales non_dim(cfg.frequency_rad_s, r_planet, bulk_rho);

        double G_nd      = G_si;
        double rho_nd    = bulk_rho;
        double freq_nd   = cfg.frequency_rad_s;
        double surf_p_nd = this->p_surface_pressure_eos;
        double start_r   = cfg.starting_radius;

        if (cfg.nondimensionalize) {
            const double L    = non_dim.length_conversion;
            const double M    = non_dim.mass_conversion;
            const double T2   = non_dim.second2_conversion;
            const double T    = non_dim.second_conversion;
            const double Pa   = non_dim.pascal_conversion;
            const double rho  = non_dim.density_conversion;
            const double grav = L / T2;
            const double moi  = M * L * L;

            for (std::size_t i = 0; i < total_slices; ++i) {
                radius_copy[i]   /= L;
                density_copy[i]  /= rho;
                gravity_copy[i]  /= grav;
                pressure_copy[i] /= Pa;
                mass_copy[i]     /= M;
                moi_copy[i]      /= moi;
                shear_vec[i]     /= Pa;
                bulk_vec[i]      /= Pa;
            }
            for (std::size_t i = 0; i < n_layers; ++i)
                upper_radii[i] /= L;

            G_nd      = G_si / (non_dim.length3_conversion / (M * T2));
            rho_nd    = bulk_rho / rho;
            freq_nd   = cfg.frequency_rad_s / (1.0 / T);
            surf_p_nd = this->p_surface_pressure_eos / Pa;
            start_r   = cfg.starting_radius / L;
        }

        // ── Create solution storage with non-dim arrays ────────────────────
        // The storage owns the internal c_EOSSolution; we inject the world's
        // already-computed (and now non-dim) EOS data into it so the shooting
        // ODE derivative functions can call eos_solution->call() and get the
        // correct gravity / density / moduli via simple array interpolation.
        this->p_love_solution = std::make_unique<::c_RadialSolutionStorage>(
            1 /*num_ytypes*/, upper_radii.data(), n_layers,
            radius_copy.data(), total_slices, cfg.degree_l
        );
        ::c_RadialSolutionStorage* storage = this->p_love_solution.get();

        storage->get_eos_solution_ptr()->inject_from_world_eos(
            radius_copy.data(),   gravity_copy.data(),
            pressure_copy.data(), mass_copy.data(),
            moi_copy.data(),      density_copy.data(),
            shear_vec.data(),     bulk_vec.data(),
            total_slices
        );
        // inject_from_world_eos sets p_use_array_interp = true so call() uses
        // linear array interpolation instead of requiring a CyRK dense result.

        // ── Count slices belonging to each layer ───────────────────────────
        std::vector<std::size_t> first_slice_idx(n_layers);
        std::vector<std::size_t> num_slices(n_layers);
        for (std::size_t layer_i = 0; layer_i < n_layers; ++layer_i) {
            first_slice_idx[layer_i] = (layer_i == 0)
                ? 0
                : first_slice_idx[layer_i - 1] + num_slices[layer_i - 1];

            const double layer_r = upper_radii[layer_i];
            std::size_t count = 0, iface = 0;
            for (std::size_t s = first_slice_idx[layer_i]; s < total_slices; ++s) {
                const double r = radius_copy[s];
                if (c_isclose(r, layer_r, 1.0e-9, 0.0)) {
                    if (++iface > 1) break;
                } else if (r > layer_r) {
                    break;
                }
                ++count;
            }
            if (count < 5) {
                storage->error_code = -5;
                storage->message    = "TidalPy: at least 5 slices per layer required";
                storage->success    = false;
                this->p_love_solved = false;
                return;
            }
            num_slices[layer_i] = count;
        }

        // ── Boundary condition model (tidal = 1, free = 0) ────────────────
        int bc_model = cfg.solve_tidal ? 1 : 0;
        const std::size_t num_bc_models = 1;

        // ── Shooting solver ────────────────────────────────────────────────
        ::c_shooting_solver(
            storage,
            freq_nd, rho_nd,
            layer_types.get(), is_static_arr.get(), is_incomp_arr.get(),
            first_slice_idx, num_slices,
            num_bc_models, &bc_model,
            G_nd,
            cfg.degree_l, cfg.use_kamata,
            start_r, cfg.start_radius_tol,
            cfg.integration_method, cfg.rtol, cfg.atol,
            cfg.scale_rtols, cfg.max_num_steps,
            cfg.expected_size, cfg.max_ram_MB, cfg.max_step,
            cfg.verbose
        );

        // ── Redimensionalise ───────────────────────────────────────────────
        if (cfg.nondimensionalize) {
            storage->dimensionalize_data(&non_dim, true);
            // dimensionalize_data scales the *_array_vec members but not the scalar
            // surface_gravity that find_love() needs; restore the SI value here.
            storage->get_eos_solution_ptr()->surface_gravity = g_surface;
        }

        // ── Extract Love numbers ───────────────────────────────────────────
        if (storage->success)
            storage->find_love();

        this->p_love_solved = storage->success;
    }

    // -----------------------------------------------------------------------
    // Love-number solve result accessors
    // All valid after solve_love_numbers succeeds; return NaN / empty otherwise.
    // The world is the sole interface — c_RadialSolutionStorage is internal.
    // -----------------------------------------------------------------------

    bool get_love_solved() const noexcept { return this->p_love_solved; }

    bool get_love_success() const noexcept {
        return this->p_love_solution ? this->p_love_solution->success : false;
    }

    int get_love_error_code() const noexcept {
        return this->p_love_solution ? this->p_love_solution->error_code : -100;
    }

    const std::string& get_love_message() const noexcept {
        static const std::string no_msg = "No love-number solve has been run.";
        return this->p_love_solution ? this->p_love_solution->message : no_msg;
    }

    std::size_t get_love_num_ytypes() const noexcept {
        return this->p_love_solution ? this->p_love_solution->num_ytypes : 0;
    }

    // Primary Love numbers (k, h, l) for the given boundary-condition ytype index.
    std::complex<double> get_love_number_k(std::size_t ytype_idx = 0) const noexcept {
        if (!this->p_love_solution || ytype_idx >= this->p_love_solution->complex_love_vec.size())
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        return this->p_love_solution->complex_love_vec[ytype_idx].k;
    }

    std::complex<double> get_love_number_h(std::size_t ytype_idx = 0) const noexcept {
        if (!this->p_love_solution || ytype_idx >= this->p_love_solution->complex_love_vec.size())
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        return this->p_love_solution->complex_love_vec[ytype_idx].h;
    }

    std::complex<double> get_love_number_l(std::size_t ytype_idx = 0) const noexcept {
        if (!this->p_love_solution || ytype_idx >= this->p_love_solution->complex_love_vec.size())
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        return this->p_love_solution->complex_love_vec[ytype_idx].l;
    }

    // Full radial solution y-value at the surface for a given ytype and y-index.
    // Returns NaN if not solved.  Matches the storage layout: 6 complex y-values
    // (12 doubles) per ytype per slice.
    std::complex<double> get_love_surface_y(
            std::size_t ytype_idx, std::size_t y_idx) const noexcept {
        if (!this->p_love_solution || !this->p_love_solution->success)
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        const std::size_t n_slices  = this->p_love_solution->num_slices;
        const std::size_t n_ytypes  = this->p_love_solution->num_ytypes;
        if (ytype_idx >= n_ytypes || y_idx >= C_MAX_NUM_Y)
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        // Layout: [slice][ytype * C_MAX_NUM_Y_REAL] interleaved real/imag pairs
        const std::size_t stride    = C_MAX_NUM_Y_REAL * n_ytypes;
        const std::size_t base      = (n_slices - 1) * stride
                                    + ytype_idx * C_MAX_NUM_Y_REAL
                                    + y_idx * 2;
        const auto& v = this->p_love_solution->full_solution_vec;
        if (base + 1 >= v.size())
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        return std::complex<double>(v[base], v[base + 1]);
    }

    // -----------------------------------------------------------------------
    // Whole-planet aggregates (const, MKS)
    // -----------------------------------------------------------------------
    // Total mass [kg] = sum of layer masses.
    double calc_total_mass() const noexcept {
        double total = 0.0;
        for (const auto& layer : this->p_layers) { total += layer->get_mass(); }
        return total;
    }

    // Total internal radiogenic heating [W] at time_s. Only SolidLiquidLayers
    // produce radiogenic heating; other layer types contribute zero.
    double calc_internal_heating(double time_s) const noexcept {
        double total = 0.0;
        for (const auto& layer : this->p_layers) {
            const auto* sl = dynamic_cast<const c_SolidLiquidLayer*>(layer.get());
            if (sl != nullptr) {
                total += sl->calc_radiogenic_heating(time_s, sl->get_mass());
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
            const double tol   = d_LAYER_CONTINUITY_RTOL * (prev_outer > 1.0 ? prev_outer : 1.0);
            if (std::abs(inner - prev_outer) > tol) { return false; }
            prev_outer = layer->get_radius_outer();
        }
        return true;
    }

    // -----------------------------------------------------------------------
    // Binary I/O — world fields + layer count, then each layer's full record.
    // -----------------------------------------------------------------------
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

    // Compute and store a layer's frequency-independent viscoelastic state over
    // its radial slice: pre-melt static moduli (from the layer's static values) +
    // pre-melt viscosities (from the attached viscosity models, NaN if unset),
    // then the post-melt versions (the partial-melt model applied to the shear and
    // bulk pairs; post == pre if no melt model). No-op for a geometry-only
    // BaseLayer. temperature_k is the placeholder profile temperature until the
    // thermal pipeline lands.
    void populate_layer_viscoelastic(
            c_BaseLayer* layer,
            const c_EOSSolution& solution,
            std::size_t slice_start,
            std::size_t slice_count,
            double temperature_k) const {
        auto* physics_layer = dynamic_cast<c_PhysicsLayer*>(layer);
        if (physics_layer == nullptr) { return; }

        const double static_shear_pa = physics_layer->get_shear_modulus_static();
        const double static_bulk_pa  = physics_layer->get_bulk_modulus_static();
        c_ViscosityBase*   shear_viscosity_model = physics_layer->get_shear_viscosity_model();
        c_ViscosityBase*   bulk_viscosity_model  = physics_layer->get_bulk_viscosity_model();
        c_PartialMeltBase* partial_melt_model    = physics_layer->get_partial_melt_model();

        std::vector<double> premelt_shear_pa(slice_count, static_shear_pa);
        std::vector<double> premelt_bulk_pa(slice_count, static_bulk_pa);
        std::vector<double> premelt_shear_visc_pas(slice_count);
        std::vector<double> premelt_bulk_visc_pas(slice_count);
        std::vector<double> postmelt_shear_pa(slice_count);
        std::vector<double> postmelt_bulk_pa(slice_count);
        std::vector<double> postmelt_shear_visc_pas(slice_count);
        std::vector<double> postmelt_bulk_visc_pas(slice_count);

        for (std::size_t slice_offset = 0; slice_offset < slice_count; ++slice_offset) {
            const double pressure_pa = solution.pressure_array_vec[slice_start + slice_offset];

            const double premelt_shear_visc = (shear_viscosity_model != nullptr)
                ? shear_viscosity_model->calc_viscosity(temperature_k, pressure_pa)
                : TidalPyConstants::d_NAN;
            const double premelt_bulk_visc = (bulk_viscosity_model != nullptr)
                ? bulk_viscosity_model->calc_viscosity(temperature_k, pressure_pa)
                : TidalPyConstants::d_NAN;
            premelt_shear_visc_pas[slice_offset] = premelt_shear_visc;
            premelt_bulk_visc_pas[slice_offset]  = premelt_bulk_visc;

            if (partial_melt_model != nullptr) {
                // Apply the melt model to the shear pair, then the bulk pair. The
                // liquid viscosity is a placeholder (the pre-melt viscosity) until a
                // dedicated liquid-viscosity model exists.
                c_PartialMeltInputs shear_inputs;
                shear_inputs.temperature_k     = temperature_k;
                shear_inputs.premelt_viscosity = premelt_shear_visc;
                shear_inputs.premelt_shear     = static_shear_pa;
                shear_inputs.liquid_viscosity  = premelt_shear_visc;
                const c_PartialMeltResult shear_result = partial_melt_model->calc_partial_melt(shear_inputs);
                postmelt_shear_pa[slice_offset]       = shear_result.postmelt_shear_modulus;
                postmelt_shear_visc_pas[slice_offset] = shear_result.postmelt_viscosity;

                c_PartialMeltInputs bulk_inputs;
                bulk_inputs.temperature_k     = temperature_k;
                bulk_inputs.premelt_viscosity = premelt_bulk_visc;
                bulk_inputs.premelt_shear     = static_bulk_pa;
                bulk_inputs.liquid_viscosity  = premelt_bulk_visc;
                const c_PartialMeltResult bulk_result = partial_melt_model->calc_partial_melt(bulk_inputs);
                postmelt_bulk_pa[slice_offset]       = bulk_result.postmelt_shear_modulus;
                postmelt_bulk_visc_pas[slice_offset] = bulk_result.postmelt_viscosity;
            } else {
                postmelt_shear_pa[slice_offset]       = static_shear_pa;
                postmelt_bulk_pa[slice_offset]        = static_bulk_pa;
                postmelt_shear_visc_pas[slice_offset] = premelt_shear_visc;
                postmelt_bulk_visc_pas[slice_offset]  = premelt_bulk_visc;
            }
        }

        layer->update_viscoelastic_data(
            premelt_shear_pa, premelt_bulk_pa, premelt_shear_visc_pas, premelt_bulk_visc_pas,
            postmelt_shear_pa, postmelt_bulk_pa, postmelt_shear_visc_pas, postmelt_bulk_visc_pas);
    }

    // Non-owning observer pointer to the layer whose radial span contains
    // radius_m [m]. Radii beyond the surface clamp to the outermost layer; radii
    // below the innermost inner radius fall in the innermost layer. Returns
    // nullptr only when the world has no layers.
    const c_BaseLayer* find_layer_for_radius(double radius_m) const noexcept {
        if (this->p_layers.empty()) { return nullptr; }
        for (const auto& layer : this->p_layers) {
            if (radius_m <= layer->get_radius_outer()) { return layer.get(); }
        }
        return this->p_layers.back().get();
    }

    std::vector<std::unique_ptr<c_BaseLayer>> p_layers;

    // EOS solve results (set by solve_eos; not serialized — repopulate by re-solving).
    bool        p_eos_success          = false;
    bool        p_eos_solved           = false;
    std::string p_eos_message          = "EOS not yet solved.";
    int         p_eos_iterations       = -1;
    double      p_eos_pressure_error   = std::numeric_limits<double>::quiet_NaN();
    double      p_surface_gravity_eos  = std::numeric_limits<double>::quiet_NaN();
    double      p_surface_pressure_eos = std::numeric_limits<double>::quiet_NaN();
    double      p_central_pressure     = std::numeric_limits<double>::quiet_NaN();
    double      p_planet_mass_eos      = std::numeric_limits<double>::quiet_NaN();
    double      p_planet_moi_eos       = std::numeric_limits<double>::quiet_NaN();
    std::shared_ptr<c_EOSSolution> p_eos_solution;  // retained full-planet solution (co-owned by layer dense evaluators)
    // Per-layer material-EOS inputs referenced by the solver's stored diffeq args;
    // must outlive solve_eos so the dense output's diffeq re-calls stay valid.
    std::vector<c_MaterialEOSInput> p_eos_material_inputs;

    // Love-number solve results (set by solve_love_numbers; not serialized).
    bool p_love_solved = false;
    std::unique_ptr<::c_RadialSolutionStorage> p_love_solution;
};

} // namespace tidalpy
