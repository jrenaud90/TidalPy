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
        auto solution = std::make_unique<c_EOSSolution>(
            upper_radii.data(), n_layers, full_radius.data(), total_slices);

        std::vector<PreEvalFunc>        eos_function_vec;
        std::vector<c_EOS_ODEInput>     eos_input_vec;
        std::vector<c_MaterialEOSInput> material_inputs(n_layers);
        eos_function_vec.reserve(n_layers);
        eos_input_vec.reserve(n_layers);

        c_EOS_ODEInput ode_input;
        ode_input.G_to_use      = G_to_use;
        ode_input.planet_radius = upper_radii.back();
        ode_input.final_solve   = false;
        ode_input.update_bulk   = false;
        ode_input.update_shear  = false;
        for (std::size_t i = 0; i < n_layers; ++i) {
            material_inputs[i].eos_model_ptr = this->p_layers[i]->get_eos();
            material_inputs[i].temperature_k = cfg.temperature;
            ode_input.eos_input_ptr          = reinterpret_cast<char*>(&material_inputs[i]);
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

        // Populate each layer's EOS profile from its slice of the full arrays.
        if (this->p_eos_solved) {
            const std::size_t n = solution->radius_array_size;
            for (std::size_t i = 0; i < n_layers; ++i) {
                const std::size_t start = i * slices;
                const std::size_t end   = start + slices;
                if (end > n) { break; }
                c_LayerEOSData eos_data;
                eos_data.populate(
                    std::vector<double>(solution->radius_array_vec.begin()   + start, solution->radius_array_vec.begin()   + end),
                    std::vector<double>(solution->density_array_vec.begin()  + start, solution->density_array_vec.begin()  + end),
                    std::vector<double>(solution->gravity_array_vec.begin()  + start, solution->gravity_array_vec.begin()  + end),
                    std::vector<double>(solution->pressure_array_vec.begin() + start, solution->pressure_array_vec.begin() + end));
                this->p_layers[i]->update_eos_data(eos_data);
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
    std::unique_ptr<c_EOSSolution> p_eos_solution;  // retained full-planet solution
};

} // namespace tidalpy
