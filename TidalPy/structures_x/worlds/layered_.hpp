#pragma once
/*
 * layered_.hpp — c_LayeredWorld: a world built from an ordered stack of layers.
 *
 * Inherits c_BaseWorld. Owns its layers as std::unique_ptr<c_BaseLayer> (inner to
 * outer, index 0 = innermost) and provides whole-planet aggregates (total mass,
 * internal radiogenic heating) and geometry validation. The whole-planet EOS and
 * radial (Love number) solves, which walk all layers, are methods on this class.
 *
 * Binary format (20-byte header + payload):
 *   header: class_id = BinaryClassID::LayeredWorld (201)
 *   payload: [all c_BaseWorld fields] + layer_count (uint64_t)
 *   Each layer's own complete binary record (with its recursively-serialized
 *   physics sub-models) follows, in index order, as separate appended records.
 */

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

// RadialSolver sub-modules: shooting solver, storage, love numbers.
// Compiled into this TU so the shooting CyRK integration runs in the same
// extension that owns the CySolverResult objects — no cross-extension call().
#include "../../Utilities_x/math_x/numerics_.hpp"        // c_isclose
#include "../../utilities/dimensions/nondimensional_.hpp"  // c_NonDimensionalScales
#include "../../utilities/arrays/interp_.hpp"        // c_interp_complex
#include "../../RadialSolver_x/rs_constants_.hpp"
#include "../../RadialSolver_x/rs_solution_.hpp"
#include "../../RadialSolver_x/love_.hpp"
#include "../../RadialSolver_x/shooting_.hpp"
#include "../../RadialSolver_x/world_radial_solver_.hpp"  // c_WorldRadialSolver (cached Love-number solver)

// Global (1D) tidal dissipation: the model hierarchy + the stored config/result structs are
// LIGHT (no global-potential tables). The calc_tides orchestration that runs the heavy
// global-potential engine + collapse is defined out-of-line in world_tides_.hpp so the
// eccentricity/obliquity tables compile into only the one extension that includes it.
#include "../../Tides_x/classes/tide_base_.hpp"     // tidalpy::c_TideBase
#include "../../Tides_x/classes/tide_result_.hpp"   // c_TideConfig, c_TideSolveConfig, c_GlobalTideResult
// Relative paths (not bare names) so every extension that includes layered_.hpp resolves
// these without needing Utilities_x/lookups on its include path. Light headers, no tables.
#include "../../Utilities_x/lookups/keys_.hpp"     // c_Key4, c_Key2
#include "../../Utilities_x/lookups/intmap_.hpp"   // c_IntMap (per-mode solver Love-number store)

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
    int       bc_model           = 1;                  // surface boundary condition: 1 = tidal, 2 = loading, 0 = free
    bool      use_prop_matrix    = false;              // false = shooting method, true = propagation matrix
    int       core_model         = 0;                  // propagation-matrix core starting condition (0-4)
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
        const double planet_bulk_density = (total_volume > TidalPyConstants::d_EPS) ? (mass_estimate / total_volume) : 3500.0;

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

        // A re-solve changes the structure/moduli arrays even if the grid size is
        // unchanged, so invalidate the cached radial-solver setup.
        if (this->p_radial_solver) this->p_radial_solver->invalidate();
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

    // -----------------------------------------------------------------------
    // Spin dynamics (the c_Spin model attached to the world; rates only)
    //
    // The world holds a spin model and drives it with its own, EOS-based moment of inertia, so the
    // spin-rate change uses the accurate structure-resolved MoI rather than the uniform-density value.
    // -----------------------------------------------------------------------
    void           set_spin_model(const c_Spin& spin) noexcept { this->p_spin = spin; }
    const c_Spin&  get_spin_model() const noexcept { return this->p_spin; }

    // Moment of inertia [kg m2]: the EOS-solved value (get_planet_moi_eos) when the EOS has been solved,
    // otherwise the spin model's uniform-density fallback from the world mass and radius.
    double get_moment_of_inertia() const noexcept {
        if (this->p_eos_solved && std::isfinite(this->p_planet_moi_eos)) {
            return this->p_planet_moi_eos;
        }
        return this->p_spin.calc_moment_of_inertia(this->get_mass(), this->get_radius(), 0.0);
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

        if (solver->cache_matches(n_layers, total_slices, cfg.degree_l, cfg.nondimensionalize)
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
            // Read the solved structure variables (gravity, ...) from the world's dense EOS, not array interpolation.
            world_eos
        );
    }

    // Build the per-call runtime config from the user-facing solve config.
    c_LoveSolveRuntimeConfig make_runtime_config(const c_LoveSolveConfig& cfg) const {
        c_LoveSolveRuntimeConfig rt;
        rt.frequency_rad_s    = cfg.frequency_rad_s;
        rt.bc_model           = cfg.bc_model;
        rt.use_prop_matrix    = cfg.use_prop_matrix;
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
        if (!this->ensure_radial_cache(cfg)) { this->p_love_solved = false; return; }
        ::c_WorldRadialSolver* solver  = this->p_radial_solver.get();
        const std::size_t total_slices = solver->total_slices();

        // Frequency-dependent step: fill the complex moduli from the layer rheology, then solve.
        const std::vector<double>& radius_si = solver->radius_si();
        std::complex<double>* shear_out = solver->shear_scratch_data();
        std::complex<double>* bulk_out  = solver->bulk_scratch_data();
        for (std::size_t i = 0; i < total_slices; ++i) {
            const double r = radius_si[i];
            shear_out[i] = this->calc_complex_shear_modulus(r, cfg.frequency_rad_s);
            bulk_out[i]  = this->calc_complex_bulk_modulus(r, cfg.frequency_rad_s);
        }

        c_LoveSolveRuntimeConfig rt = this->make_runtime_config(cfg);
        solver->solve(rt);
        this->p_love_solved = solver->get_solved();
    }

    // Solve using externally-supplied complex moduli (the standalone array API path) instead of the layer rheology.
    // shear_in / bulk_in are defined at the radii radius_in (length n_in) and are linearly interpolated onto the
    // world EOS radius grid. Runs in export mode: the storage EOS arrays are redimensionalized and its scalar
    // results are copied from the world own EOS so the released solution reports SI values.
    void solve_love_numbers_supplied(
            const c_LoveSolveConfig& cfg,
            const std::complex<double>* shear_in,
            const std::complex<double>* bulk_in,
            const double* radius_in,
            std::size_t n_in) {
        if (!this->ensure_radial_cache(cfg)) { this->p_love_solved = false; return; }
        ::c_WorldRadialSolver* solver  = this->p_radial_solver.get();
        const std::size_t total_slices = solver->total_slices();

        const std::vector<double>& radius_si = solver->radius_si();
        std::complex<double>* shear_out = solver->shear_scratch_data();
        std::complex<double>* bulk_out  = solver->bulk_scratch_data();
        // Interpolate the supplied complex moduli (defined at radius_in) onto the world's EOS grid via the shared
        // utilities interpolator. Both grids ascend, so seeding the search from the aligned fractional
        // position keeps each lookup near O(1).
        for (std::size_t i = 0; i < total_slices; ++i) {
            const double r = radius_si[i];
            const std::size_t seed = (total_slices > 0) ? (i * n_in) / total_slices : 0;
            shear_out[i] = c_interp_complex(r, radius_in, shear_in, n_in, seed);
            bulk_out[i]  = c_interp_complex(r, radius_in, bulk_in, n_in, seed);
        }

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

    // Move the radial-solution storage out of the helper (one-shot export to a RadialSolverSolution).
    std::unique_ptr<::c_RadialSolutionStorage> release_radial_storage() {
        return this->p_radial_solver ? this->p_radial_solver->release_storage() : nullptr;
    }

    // -----------------------------------------------------------------------
    // Love-number solve result accessors
    // All valid after solve_love_numbers succeeds; return NaN / empty otherwise.
    // The world is the sole interface — c_RadialSolutionStorage is internal.
    // -----------------------------------------------------------------------

    // Non-owning pointer to the internal solution storage (owned by the cached
    // radial solver). Null until solve_love_numbers has built the cache.
    const ::c_RadialSolutionStorage* get_love_storage() const noexcept {
        return this->p_radial_solver ? this->p_radial_solver->get_storage() : nullptr;
    }

    bool get_love_solved() const noexcept { return this->p_love_solved; }

    bool get_love_success() const noexcept {
        const auto* s = this->get_love_storage();
        return s ? s->success : false;
    }

    int get_love_error_code() const noexcept {
        const auto* s = this->get_love_storage();
        return s ? s->error_code : -100;
    }

    const std::string& get_love_message() const noexcept {
        static const std::string no_msg = "No love-number solve has been run.";
        const auto* s = this->get_love_storage();
        return s ? s->message : no_msg;
    }

    std::size_t get_love_num_ytypes() const noexcept {
        const auto* s = this->get_love_storage();
        return s ? s->num_ytypes : 0;
    }

    // Worst-case error amplification of the surface boundary condition solve (shooting method; 0 until a
    // solve has run). See c_estimate_surface_amplification in RadialSolver_x/boundaries/boundaries_.hpp.
    double get_love_surface_amplification() const noexcept {
        const auto* s = this->get_love_storage();
        return s ? s->surface_amplification : 0.0;
    }

    // Primary Love numbers (k, h, l) for the given boundary-condition ytype index.
    std::complex<double> get_love_number_k(std::size_t ytype_idx = 0) const noexcept {
        const auto* s = this->get_love_storage();
        if (!s || ytype_idx >= s->complex_love_vec.size())
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        return s->complex_love_vec[ytype_idx].k;
    }

    std::complex<double> get_love_number_h(std::size_t ytype_idx = 0) const noexcept {
        const auto* s = this->get_love_storage();
        if (!s || ytype_idx >= s->complex_love_vec.size())
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        return s->complex_love_vec[ytype_idx].h;
    }

    std::complex<double> get_love_number_l(std::size_t ytype_idx = 0) const noexcept {
        const auto* s = this->get_love_storage();
        if (!s || ytype_idx >= s->complex_love_vec.size())
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        return s->complex_love_vec[ytype_idx].l;
    }

    // Full radial solution y-value (SI) at the surface for a given ytype and y-index (0..5 -> y1..y6).
    // Returns NaN if not solved. Evaluated from the radial solver's dense calling system (shooting) or the
    // gridded solution (matrix) via get_surface_y; no longer reads the gridded buffer directly.
    std::complex<double> get_love_surface_y(
            std::size_t ytype_idx, std::size_t y_idx) const noexcept {
        const auto* s = this->get_love_storage();
        if (!s || !s->success || y_idx >= C_MAX_NUM_Y)
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        std::complex<double> surface_y[C_MAX_NUM_Y];
        if (!s->get_surface_y(ytype_idx, surface_y))
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        return surface_y[y_idx];
    }

    // Full radial solution y-value (SI) at an arbitrary radius [m] for a given ytype and y-index. The shooting
    // method evaluates its dense per-layer interpolants at this radius (accurate anywhere, including between EOS
    // grid slices); the matrix method linearly interpolates its constructed grid. Returns NaN if not solved,
    // out of range, or below the solver's starting radius.
    std::complex<double> get_radial_solution_y(
            double radius_m,
            std::size_t ytype_idx,
            std::size_t y_idx) const noexcept {
        const auto* s = this->get_love_storage();
        if (!s || !s->success || y_idx >= C_MAX_NUM_Y)
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);

        std::complex<double> y_at_r[C_MAX_NUM_Y];
        if (!s->get_radial_solution(radius_m, ytype_idx, y_at_r))
            return std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        return y_at_r[y_idx];
    }

    // -----------------------------------------------------------------------
    // Global (1D) tidal dissipation
    //
    // The tide-model holder, config, results, and the analytic calc_tides path live on
    // c_BaseWorld (common to all world types). c_LayeredWorld extends calc_tides with the
    // rheology (radial-solver) path and the per-layer heating distribution. calc_tides is
    // defined out-of-line in world_tides_.hpp (it needs the heavy global-potential engine).
    // -----------------------------------------------------------------------
    // Run the global tidal solve for the supplied orbital/spin state (defined in
    // world_tides_.hpp). Hides c_BaseWorld::calc_tides to add rheology + layer distribution.
    void calc_tides(const c_TideSolveConfig& state);

    // -----------------------------------------------------------------------
    // On-demand 3D tidal stress/strain/heating
    //
    // The tidal potential is built dynamically from the world's [tides] truncation config (max degree
    // l, eccentricity/obliquity truncation) — there is no potential-model object. The 3D orchestration
    // lives on the rheology tide model (c_RheologyTide); get_3d_tidal_heating delegates to it and is
    // defined out-of-line in world_tides_.hpp (it needs the kernel + potential-engine headers).
    // -----------------------------------------------------------------------
    // On-demand secular (cycle/orbit-averaged) 3D tidal volumetric heating [W m-3] at (radius,
    // colatitude). The physically time-averaged power density (single omega/2, complex amplitudes, no
    // abs); longitude- and time-independent. Requires the rheology tide model and a solved EOS. Its
    // volume integral equals the 1D global heating (get_tidal_heating). Defined out-of-line in
    // world_tides_.hpp.
    double get_3d_tidal_heating(
            const c_TideSolveConfig& state,
            double radius,
            double colatitude);

    // Vectorized batch form: secular 3D volumetric heating [W m-3] at num_points paired
    // (radii[i], colatitudes[i]), written into out_heating[i]. Same physics/preconditions as the scalar
    // get_3d_tidal_heating, but the radial solve is amortized across points (one solve per unique (l,
    // frequency) rather than per point), so it is the efficient way to build a map. Defined out-of-line
    // in world_tides_.hpp.
    void get_3d_tidal_heating_array(
            const c_TideSolveConfig& state,
            const double* radii,
            const double* colatitudes,
            size_t num_points,
            double* out_heating);

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

        const double const_shear_pa = physics_layer->get_shear_modulus_static();
        const double const_bulk_pa  = physics_layer->get_bulk_modulus_static();
        c_ViscosityBase*   shear_viscosity_model = physics_layer->get_shear_viscosity_model();
        c_ViscosityBase*   bulk_viscosity_model  = physics_layer->get_bulk_viscosity_model();
        c_PartialMeltBase* partial_melt_model    = physics_layer->get_partial_melt_model();

        // The static moduli and viscosities are carried as extra outputs of the EOS
        // ODE (filled by each layer's EOS model during the CyRK solve) and stored on
        // the solution. For a layer whose EOS supplies them (the interpolated /PREM
        // model) these radius-varying values are used; a NaN means "not provided", so
        // the layer's constant modulus (or its viscosity model) is used instead.
        const bool have_eos_extras =
            solution.other_vecs_set
            && solution.complex_shear_array_vec.size() == solution.radius_array_size
            && solution.shear_viscosity_array_vec.size() == solution.radius_array_size;

        std::vector<double> premelt_shear_pa(slice_count);
        std::vector<double> premelt_bulk_pa(slice_count);
        std::vector<double> premelt_shear_visc_pas(slice_count);
        std::vector<double> premelt_bulk_visc_pas(slice_count);
        std::vector<double> postmelt_shear_pa(slice_count);
        std::vector<double> postmelt_bulk_pa(slice_count);
        std::vector<double> postmelt_shear_visc_pas(slice_count);
        std::vector<double> postmelt_bulk_visc_pas(slice_count);

        for (std::size_t slice_offset = 0; slice_offset < slice_count; ++slice_offset) {
            const std::size_t global_slice = slice_start + slice_offset;
            const double pressure_pa = solution.pressure_array_vec[global_slice];

            // Static moduli: prefer the EOS extra-output (CyRK solution), else the
            // layer constant.
            double static_shear_pa = const_shear_pa;
            double static_bulk_pa  = const_bulk_pa;
            if (have_eos_extras) {
                const double eos_shear = solution.complex_shear_array_vec[global_slice].real();
                if (std::isfinite(eos_shear)) { static_shear_pa = eos_shear; }
                const double eos_bulk = solution.complex_bulk_array_vec[global_slice].real();
                if (std::isfinite(eos_bulk)) { static_bulk_pa = eos_bulk; }
            }
            premelt_shear_pa[slice_offset] = static_shear_pa;
            premelt_bulk_pa[slice_offset]  = static_bulk_pa;

            // Pre-melt viscosities: viscosity-model value, overridden by an EOS
            // extra-output (CyRK solution) when the EOS provides one.
            double premelt_shear_visc = (shear_viscosity_model != nullptr)
                ? shear_viscosity_model->calc_viscosity(temperature_k, pressure_pa)
                : TidalPyConstants::d_NAN;
            double premelt_bulk_visc = (bulk_viscosity_model != nullptr)
                ? bulk_viscosity_model->calc_viscosity(temperature_k, pressure_pa)
                : TidalPyConstants::d_NAN;
            if (have_eos_extras) {
                const double eos_shear_visc = solution.shear_viscosity_array_vec[global_slice];
                if (std::isfinite(eos_shear_visc)) { premelt_shear_visc = eos_shear_visc; }
                const double eos_bulk_visc = solution.bulk_viscosity_array_vec[global_slice];
                if (std::isfinite(eos_bulk_visc)) { premelt_bulk_visc = eos_bulk_visc; }
            }
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
    // nullptr only when the world has no layers. Public (a const observer query, like
    // get_density(radius)); the 3D tidal-heating path reads the layer's solid/incompressible flags.
public:
    const c_BaseLayer* find_layer_for_radius(double radius_m) const noexcept {
        if (this->p_layers.empty()) { return nullptr; }
        for (const auto& layer : this->p_layers) {
            if (radius_m <= layer->get_radius_outer()) { return layer.get(); }
        }
        return this->p_layers.back().get();
    }

protected:
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
    c_Spin      p_spin {};              // spin-dynamics model (uses the world's EOS moment of inertia)
    std::shared_ptr<c_EOSSolution> p_eos_solution;  // retained full-planet solution (co-owned by layer dense evaluators)
    // Per-layer material-EOS inputs referenced by the solver's stored diffeq args;
    // must outlive solve_eos so the dense output's diffeq re-calls stay valid.
    std::vector<c_MaterialEOSInput> p_eos_material_inputs;

    // Love-number solve: cached, reusable radial solver (set by solve_love_numbers;
    // not serialized). Holds the frequency-independent setup + the reused solution
    // storage; rebuilt only when the EOS grid/assumptions change.
    bool p_love_solved = false;
    std::unique_ptr<::c_WorldRadialSolver> p_radial_solver;

    // Global (1D) tidal dissipation: the model/config/result state lives on c_BaseWorld;
    // c_LayeredWorld adds only the per-layer heating distribution (results not serialized).
    std::vector<double>                  p_layer_tidal_heating;
};

} // namespace tidalpy
