#pragma once
/*
 * heating_.hpp: the heat generated inside a layered world, as the thermal structure solve reads it.
 *
 * c_Heating belongs to a world and sums that world's heat sources into a volumetric heating h(r) [W m-3]. The
 * structure ODE integrates dL/dr = 4 pi r^2 h, so the heat flow L(r) and the conductive temperature profile
 * follow the heating at every radius. Only layers with `use_heating` set are heated; the rest report zero.
 *
 * A source is prepared once per solve (solve_radial_heat) and then evaluated inside the right-hand side, so
 * its evaluation allocates nothing and casts nothing. A source whose heating depends on the solved interior
 * says so through is_state_dependent, which is what tells a solve it has to iterate on the heating.
 *
 * Sources
 * -------
 * radiogenic   each layer's radiogenics model gives a specific rate [W kg-1] at the solve's time, and the
 *              heating is that rate times the local density. It does not depend on the state, so it is exact
 *              on every pass, also for a layer whose mass is an output of the solve.
 *
 * All quantities are MKS unless a method says it works in the units of the solve.
 */

#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

#include "constants_.hpp"
#include "ode_.hpp"                        // c_EOSHeatingBase
#include "../layers/base_.hpp"
#include "../layers/physics_.hpp"
#include "../layers/solidliquid_.hpp"
#include "../../radiogenics_x/radiogenics_base_.hpp"

namespace tidalpy {

// What the heat sources are told about the world before a solve.
struct c_WorldState {
    // Time [s] on the clock the radiogenics models share. NaN asks each model for its own reference time.
    double time = TidalPyConstants::d_NAN;
    // The world's layers, inner to outer (non-owning).
    const std::vector<std::unique_ptr<c_BaseLayer>>* layers_ptr = nullptr;
};

// -------------------------------------------------------------------------------
// c_HeatSourceBase: one physical source of internal heat.
// -------------------------------------------------------------------------------
class c_HeatSourceBase {
public:
    virtual ~c_HeatSourceBase() = default;

    // True when the heating depends on the solved interior, so a solve has to iterate on it.
    virtual bool is_state_dependent() const noexcept = 0;

    // Prepare the source for a solve of the world in this state. It takes the world's state, not one
    // layer's, because a source can span the planet.
    virtual void solve_radial_heat(const c_WorldState& state) = 0;

    // Volumetric heating [W m-3] at a radius [m] of one layer, where the local density is `density` [kg m-3].
    virtual double calc_heating(std::size_t layer_index, double radius, double density) const noexcept = 0;
};

// -------------------------------------------------------------------------------
// c_RadiogenicHeatSource: the layers' radiogenics models. Each model's decay law is linear in the mass it is
// handed, so its specific rate is its heating of one kilogram.
// -------------------------------------------------------------------------------
class c_RadiogenicHeatSource : public c_HeatSourceBase {
public:
    bool is_state_dependent() const noexcept override { return false; }

    void solve_radial_heat(const c_WorldState& state) override {
        this->p_specific_heating_bylayer.clear();
        if (state.layers_ptr == nullptr) { return; }
        const std::size_t n_layers = state.layers_ptr->size();
        this->p_specific_heating_bylayer.assign(n_layers, 0.0);
        for (std::size_t layer_i = 0; layer_i < n_layers; ++layer_i) {
            const auto* solidliquid_layer =
                dynamic_cast<const c_SolidLiquidLayer*>((*state.layers_ptr)[layer_i].get());
            if (solidliquid_layer == nullptr) { continue; }
            const c_RadiogenicsBase* radiogenics_model = solidliquid_layer->get_radiogenics_model();
            if (radiogenics_model == nullptr) { continue; }
            const double time = std::isfinite(state.time) ? state.time : radiogenics_model->get_ref_time();
            const double specific_heating = radiogenics_model->calc_heating(time, 1.0);
            if (std::isfinite(specific_heating)) {
                this->p_specific_heating_bylayer[layer_i] = specific_heating;
            }
        }
    }

    double calc_heating(std::size_t layer_index, double /*radius*/, double density) const noexcept override {
        if (layer_index >= this->p_specific_heating_bylayer.size()) { return 0.0; }
        return this->p_specific_heating_bylayer[layer_index] * density;
    }

    // Specific heating rate [W kg-1] of a layer at the time of the last solve_radial_heat.
    double get_specific_heating(std::size_t layer_index) const noexcept {
        return (layer_index < this->p_specific_heating_bylayer.size())
            ? this->p_specific_heating_bylayer[layer_index] : 0.0;
    }

protected:
    std::vector<double> p_specific_heating_bylayer;   // [W kg-1]
};

// -------------------------------------------------------------------------------
// c_Heating: a world's heat sources, summed.
// -------------------------------------------------------------------------------
class c_Heating : public c_EOSHeatingBase {
public:
    c_Heating() { this->p_sources.push_back(&this->p_radiogenic_source); }
    ~c_Heating() override = default;

    // The source list points into this object, so it is not copied.
    c_Heating(const c_Heating&) = delete;
    c_Heating& operator=(const c_Heating&) = delete;

    // Prepare every source for a solve of the world in `state`, and record which layers are heated and the
    // length and density units the solve runs in (one for an SI solve).
    void update_sources(const c_WorldState& state, double length_scale, double density_scale) {
        this->p_length_scale  = length_scale;
        this->p_density_scale = density_scale;
        this->p_use_heating_bylayer.clear();
        if (state.layers_ptr != nullptr) {
            for (const auto& layer_uptr : *state.layers_ptr) {
                const auto* physics_layer = dynamic_cast<const c_PhysicsLayer*>(layer_uptr.get());
                this->p_use_heating_bylayer.push_back(
                    ((physics_layer != nullptr) && physics_layer->get_use_heating()) ? 1 : 0);
            }
        }
        for (c_HeatSourceBase* source_ptr : this->p_sources) { source_ptr->solve_radial_heat(state); }
    }

    // True when any layer is heated, so a solve has a heat flow to integrate.
    bool get_is_active() const noexcept {
        for (const char use_heating : this->p_use_heating_bylayer) {
            if (use_heating) { return true; }
        }
        return false;
    }

    // True when a source depends on the solved interior.
    bool get_is_state_dependent() const noexcept {
        for (const c_HeatSourceBase* source_ptr : this->p_sources) {
            if (source_ptr->is_state_dependent()) { return true; }
        }
        return false;
    }

    // Volumetric heating [W m-3] at a radius [m] of one layer, where the local density is `density` [kg m-3]:
    // every source summed, and zero in a layer that is not heated.
    double calc_heating(std::size_t layer_index, double radius, double density) const noexcept {
        if ((layer_index >= this->p_use_heating_bylayer.size()) || !this->p_use_heating_bylayer[layer_index]) {
            return 0.0;
        }
        double heating = 0.0;
        for (const c_HeatSourceBase* source_ptr : this->p_sources) {
            heating += source_ptr->calc_heating(layer_index, radius, density);
        }
        return heating;
    }

    // dL/dr in the units the solve runs in (see c_EOSHeatingBase). The heat flow itself stays in Watts.
    double calc_heat_flow_gradient(std::size_t layer_index, double radius, double density) const noexcept override {
        const double radius_si = radius * this->p_length_scale;
        return 4.0 * TidalPyConstants::d_PI * radius_si * radius_si * this->p_length_scale
            * this->calc_heating(layer_index, radius_si, density * this->p_density_scale);
    }

    const c_RadiogenicHeatSource& get_radiogenic_source() const noexcept { return this->p_radiogenic_source; }

protected:
    c_RadiogenicHeatSource p_radiogenic_source;
    // Non-owning pointers to the sources above, in the order they are summed.
    std::vector<c_HeatSourceBase*> p_sources;
    std::vector<char> p_use_heating_bylayer;
    double p_length_scale  = 1.0;
    double p_density_scale = 1.0;
};

}  // namespace tidalpy
