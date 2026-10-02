#pragma once
/*
 * heating_.hpp: the heat generated inside a world, as the thermal structure solve reads it.
 *
 * c_Heating sums a world's heat sources into a volumetric heating h(r) [W m-3]. The structure ODE integrates
 * dL/dr = 4 pi r^2 h, so the heat flow L(r) and the conductive temperature profile follow the heating at every radius.
 * Only layers with `use_heating` set are heated; the rest report zero. Each EOS solve builds its own c_Heating from
 * the world's state (c_WorldState), so a solve and the solution it leaves behind keep the sources it was solved with.
 *
 * A source is prepared once per solve (solve_radial_heat) and then evaluated inside the right-hand side, so its
 * evaluation allocates nothing and casts nothing. A source that spreads a power [W] over a layer by mass needs the
 * layer's mass, which is an output of the solve; the solve hands every source the layers as its last pass left them
 * (update_layers), so on a converged thermal solve the layer receives that power to the solve's thermal tolerance.
 *
 * Sources
 * -------
 * radiogenic   each layer's radiogenics model gives a specific rate [W kg-1] at the solve's time, and the heating is
 *              that rate times the local density. It does not depend on the state.
 * tidal        the per-layer heating [W] of the world's last calc_tides (c_TidalHeatingRecord). Where calc_tides
 *              integrated the radial solution over the layer, the heating follows that radial profile; otherwise it
 *              is spread by mass. It depends on the solved interior through the tides, so an evolution alternates
 *              calc_tides and solve_eos (each solve takes the tides of the structure before it).
 * prescribed   a power [W], spread by mass, or a specific rate [W kg-1], per layer, set by the caller
 *              (c_PrescribedHeating).
 *
 * All quantities are MKS unless a method says it works in the units of the solve.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "constants_.hpp"
#include "ode_.hpp"                        // c_EOSHeatingBase
#include "../layers/layer_.hpp"
#include "../../Radiogenics/radiogenics_base_.hpp"

namespace tidalpy {

// The heat sources, in the order c_Heating sums them.
enum class c_HeatSourceKind : uint8_t {
    Radiogenic = 0,
    Tidal = 1,
    Prescribed = 2
};
inline constexpr std::size_t C_NUM_HEAT_SOURCES = 3;

// A layer as the sources that spread a power over it see it during a solve: its radii [m] and mass [kg] as the last
// pass left them (its own radii and current mass before the first).
struct c_HeatedLayer {
    double radius_inner = 0.0;
    double radius_outer = 0.0;
    double mass = TidalPyConstants::d_NAN;
};

// The last calc_tides heating of each layer, which the tidal source spreads. A layer with a radial profile has the
// heating density [W m-3] (the shell power dP/dr over 4 pi r^2) at nodes x in (0, 1) across the layer
// (x = (r - r_inner) / (r_outer - r_inner)); one without is spread by mass.
struct c_TidalHeatingRecord {
    std::vector<double> layer_power;                     // [W]; NaN or empty for no tidal heating
    std::vector<std::vector<double>> profile_fraction;   // node positions x, ascending
    std::vector<std::vector<double>> profile_density;    // heating density at the nodes [W m-3]

    bool empty() const noexcept { return this->layer_power.empty(); }
};

// What the caller prescribes for one layer: a power [W] spread by mass, or a specific rate [W kg-1]. NaN leaves it
// unset.
struct c_PrescribedLayerHeating {
    double power = TidalPyConstants::d_NAN;
    double specific_rate = TidalPyConstants::d_NAN;

    bool is_set() const noexcept { return std::isfinite(this->power) || std::isfinite(this->specific_rate); }
};

// What the heat sources are told about the world before a solve.
struct c_WorldState {
    // Time [s] on the clock the radiogenics models share. NaN asks each model for its own reference time.
    double time = TidalPyConstants::d_NAN;
    // The world's layers, inner to outer (non-owning).
    const std::vector<std::unique_ptr<c_Layer>>* layers_ptr = nullptr;
    // The world's tidal heating record and prescribed heating, per layer (non-owning; null or empty for none).
    const c_TidalHeatingRecord* tidal_ptr = nullptr;
    const std::vector<c_PrescribedLayerHeating>* prescribed_ptr = nullptr;
};

// c_HeatSourceBase: one physical source of internal heat.
class c_HeatSourceBase {
public:
    virtual ~c_HeatSourceBase() = default;

    // True when the heating depends on the solved interior, so a solve has to iterate on it.
    virtual bool is_state_dependent() const noexcept = 0;

    // Prepare the source for a solve of the world in this state. It takes the world's state, not one
    // layer's, because a source can span the planet.
    virtual void solve_radial_heat(const c_WorldState& state) = 0;

    // The layers' radii and masses as the last pass left them, for a source that spreads a power over a layer.
    virtual void update_layers(const std::vector<c_HeatedLayer>& layers) { (void)layers; }

    // Volumetric heating [W m-3] at a radius [m] of one layer, where the local density is `density` [kg m-3].
    virtual double calc_heating(std::size_t layer_index, double radius, double density) const noexcept = 0;

    // The heat [W] the source generates in a layer of mass `mass` [kg]: the volume integral of calc_heating over
    // the layer, from the source's own definition rather than a quadrature.
    virtual double calc_layer_power(std::size_t layer_index, double mass) const noexcept = 0;
};

// c_RadiogenicHeatSource: the layers' radiogenics models. Each model's decay law is linear in the mass it is
// handed, so its specific rate is its heating of one kilogram.
class c_RadiogenicHeatSource : public c_HeatSourceBase {
public:
    bool is_state_dependent() const noexcept override { return false; }

    void solve_radial_heat(const c_WorldState& state) override {
        this->p_specific_heating_bylayer.clear();
        if (state.layers_ptr == nullptr) { return; }
        const std::size_t n_layers = state.layers_ptr->size();
        this->p_specific_heating_bylayer.assign(n_layers, 0.0);
        for (std::size_t layer_i = 0; layer_i < n_layers; ++layer_i) {
            const c_RadiogenicsBase* radiogenics_model = (*state.layers_ptr)[layer_i]->get_radiogenics_model();
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

    double calc_layer_power(std::size_t layer_index, double mass) const noexcept override {
        return this->get_specific_heating(layer_index) * mass;
    }

    // Specific heating rate [W kg-1] of a layer at the time of the last solve_radial_heat.
    double get_specific_heating(std::size_t layer_index) const noexcept {
        return (layer_index < this->p_specific_heating_bylayer.size())
            ? this->p_specific_heating_bylayer[layer_index] : 0.0;
    }

protected:
    std::vector<double> p_specific_heating_bylayer;   // [W kg-1]
};

// c_PowerHeatSource: the shared part of the sources that hold a power [W] or a specific rate [W kg-1] per layer. A
// power is spread by mass, h = (P / M) rho, with M the layer's mass from update_layers; a specific rate s gives
// h = s rho directly.
class c_PowerHeatSource : public c_HeatSourceBase {
public:
    void update_layers(const std::vector<c_HeatedLayer>& layers) override {
        const std::size_t n_layers = std::min(layers.size(), this->p_power_bylayer.size());
        for (std::size_t layer_i = 0; layer_i < n_layers; ++layer_i) {
            const double power = this->p_power_bylayer[layer_i];
            const double mass = layers[layer_i].mass;
            if (!std::isfinite(power)) { continue; }
            this->p_specific_heating_bylayer[layer_i] =
                (std::isfinite(mass) && (mass > TidalPyConstants::d_EPS)) ? power / mass : 0.0;
        }
    }

    double calc_heating(std::size_t layer_index, double /*radius*/, double density) const noexcept override {
        if (layer_index >= this->p_specific_heating_bylayer.size()) { return 0.0; }
        return this->p_specific_heating_bylayer[layer_index] * density;
    }

    double calc_layer_power(std::size_t layer_index, double mass) const noexcept override {
        if (layer_index >= this->p_specific_heating_bylayer.size()) { return 0.0; }
        return this->p_specific_heating_bylayer[layer_index] * mass;
    }

protected:
    // The layers' powers [W] (NaN for none) and the specific rates [W kg-1] they give (or that were set directly).
    void p_reset(std::size_t n_layers) {
        this->p_power_bylayer.assign(n_layers, TidalPyConstants::d_NAN);
        this->p_specific_heating_bylayer.assign(n_layers, 0.0);
    }

    std::vector<double> p_power_bylayer;   // [W]
    std::vector<double> p_specific_heating_bylayer;   // [W kg-1]
};

// The integral of 4 pi r^2 h(r) [W] over [radius_lower, radius_upper] [m] for a heating density h [W m-3] linear in r
// between h_lower and h_upper.
inline double c_linear_shell_power(double radius_lower, double radius_upper, double h_lower, double h_upper) noexcept {
    const double width = radius_upper - radius_lower;
    if (!(width > 0.0)) { return 0.0; }
    const double slope = (h_upper - h_lower) / width;
    const double cube_difference = radius_upper * radius_upper * radius_upper
        - radius_lower * radius_lower * radius_lower;
    const double fourth_difference = radius_upper * radius_upper * radius_upper * radius_upper
        - radius_lower * radius_lower * radius_lower * radius_lower;
    return 4.0 * TidalPyConstants::d_PI
        * ((h_lower - slope * radius_lower) * cube_difference / 3.0 + slope * fourth_difference / 4.0);
}

// c_TidalHeatSource: the world's last calc_tides heating (c_TidalHeatingRecord). A layer with a radial profile takes
// h(r) = s g(x), with g the piecewise-linear interpolant of its heating densities over x in [0, 1] (flat past the end
// nodes; the density is finite at the center of a layer that reaches it) and s the factor that makes the integral of
// 4 pi r^2 h over the layer, on the radii it has in the current pass, the layer's heating. So the layer receives
// exactly its tidal heating wherever the solve puts it. A layer without a profile is spread by mass.
class c_TidalHeatSource : public c_PowerHeatSource {
public:
    bool is_state_dependent() const noexcept override { return true; }

    void solve_radial_heat(const c_WorldState& state) override {
        const std::size_t n_layers = (state.layers_ptr == nullptr) ? 0 : state.layers_ptr->size();
        this->p_reset(n_layers);
        this->p_fraction_bylayer.assign(n_layers, std::vector<double>());
        this->p_density_bylayer.assign(n_layers, std::vector<double>());
        this->p_inner_bylayer.assign(n_layers, 0.0);
        this->p_span_bylayer.assign(n_layers, 0.0);
        this->p_scale_bylayer.assign(n_layers, 0.0);
        this->p_layer_power_bylayer.assign(n_layers, 0.0);
        if ((state.tidal_ptr == nullptr) || state.tidal_ptr->empty()) { return; }
        const c_TidalHeatingRecord& record = *state.tidal_ptr;
        for (std::size_t layer_i = 0; layer_i < std::min(n_layers, record.layer_power.size()); ++layer_i) {
            const double power = record.layer_power[layer_i];
            if (!std::isfinite(power)) { continue; }
            const bool has_profile = (layer_i < record.profile_fraction.size())
                && (layer_i < record.profile_density.size())
                && !record.profile_fraction[layer_i].empty()
                && (record.profile_fraction[layer_i].size() == record.profile_density[layer_i].size());
            if (!has_profile) {
                this->p_power_bylayer[layer_i] = power;
                continue;
            }
            // The profile with its ends; update_layers scales it to the layer's power on the pass's radii.
            const std::vector<double>& nodes  = record.profile_fraction[layer_i];
            const std::vector<double>& values = record.profile_density[layer_i];
            std::vector<double> fraction = {0.0};
            std::vector<double> density = {values.front()};
            fraction.insert(fraction.end(), nodes.begin(), nodes.end());
            density.insert(density.end(), values.begin(), values.end());
            fraction.push_back(1.0);
            density.push_back(values.back());
            this->p_fraction_bylayer[layer_i]    = std::move(fraction);
            this->p_density_bylayer[layer_i]     = std::move(density);
            this->p_layer_power_bylayer[layer_i] = power;
        }
    }

    void update_layers(const std::vector<c_HeatedLayer>& layers) override {
        c_PowerHeatSource::update_layers(layers);
        for (std::size_t layer_i = 0; layer_i < std::min(layers.size(), this->p_span_bylayer.size()); ++layer_i) {
            const double radius_inner = layers[layer_i].radius_inner;
            const double span = layers[layer_i].radius_outer - radius_inner;
            this->p_inner_bylayer[layer_i] = radius_inner;
            this->p_span_bylayer[layer_i]  = span;
            const std::vector<double>& fraction = this->p_fraction_bylayer[layer_i];
            const std::vector<double>& density  = this->p_density_bylayer[layer_i];
            if (fraction.empty() || !(span > 0.0)) { continue; }
            double integral = 0.0;
            for (std::size_t node_i = 1; node_i < fraction.size(); ++node_i) {
                integral += c_linear_shell_power(
                    radius_inner + fraction[node_i - 1] * span, radius_inner + fraction[node_i] * span,
                    density[node_i - 1], density[node_i]);
            }
            this->p_scale_bylayer[layer_i] = (std::isfinite(integral) && (std::fabs(integral) > 0.0))
                ? this->p_layer_power_bylayer[layer_i] / integral : 0.0;
        }
    }

    double calc_heating(std::size_t layer_index, double radius, double density) const noexcept override {
        if ((layer_index >= this->p_density_bylayer.size()) || this->p_density_bylayer[layer_index].empty()) {
            return c_PowerHeatSource::calc_heating(layer_index, radius, density);
        }
        const double span = this->p_span_bylayer[layer_index];
        if (!(span > 0.0)) { return 0.0; }
        const std::vector<double>& fraction = this->p_fraction_bylayer[layer_index];
        const std::vector<double>& profile  = this->p_density_bylayer[layer_index];
        const double x = std::clamp((radius - this->p_inner_bylayer[layer_index]) / span, 0.0, 1.0);
        const std::size_t upper = std::min<std::size_t>(
            std::max<std::size_t>(std::upper_bound(fraction.begin(), fraction.end(), x) - fraction.begin(), 1),
            fraction.size() - 1);
        const double width = fraction[upper] - fraction[upper - 1];
        const double weight = (width > 0.0) ? (x - fraction[upper - 1]) / width : 0.0;
        return this->p_scale_bylayer[layer_index]
            * (profile[upper - 1] + weight * (profile[upper] - profile[upper - 1]));
    }

    double calc_layer_power(std::size_t layer_index, double mass) const noexcept override {
        if ((layer_index < this->p_density_bylayer.size()) && !this->p_density_bylayer[layer_index].empty()) {
            return (this->p_scale_bylayer[layer_index] != 0.0) ? this->p_layer_power_bylayer[layer_index] : 0.0;
        }
        return c_PowerHeatSource::calc_layer_power(layer_index, mass);
    }

protected:
    std::vector<std::vector<double>> p_fraction_bylayer;
    std::vector<std::vector<double>> p_density_bylayer;   // heating density at the nodes [W m-3], unscaled
    std::vector<double> p_inner_bylayer;                  // [m]
    std::vector<double> p_span_bylayer;                   // [m]
    std::vector<double> p_scale_bylayer;                  // [dimensionless]
    std::vector<double> p_layer_power_bylayer;            // a profiled layer's heating [W]
};

// c_PrescribedHeatSource: the caller's per-layer power or specific rate (c_PrescribedLayerHeating).
class c_PrescribedHeatSource : public c_PowerHeatSource {
public:
    bool is_state_dependent() const noexcept override { return false; }

    void solve_radial_heat(const c_WorldState& state) override {
        const std::size_t n_layers = (state.layers_ptr == nullptr) ? 0 : state.layers_ptr->size();
        this->p_reset(n_layers);
        if (state.prescribed_ptr == nullptr) { return; }
        const std::vector<c_PrescribedLayerHeating>& prescribed = *state.prescribed_ptr;
        for (std::size_t layer_i = 0; layer_i < std::min(n_layers, prescribed.size()); ++layer_i) {
            if (std::isfinite(prescribed[layer_i].specific_rate)) {
                this->p_specific_heating_bylayer[layer_i] = prescribed[layer_i].specific_rate;
            } else if (std::isfinite(prescribed[layer_i].power)) {
                this->p_power_bylayer[layer_i] = prescribed[layer_i].power;
            }
        }
    }
};

// c_Heating: a world's heat sources, summed.
class c_Heating : public c_EOSHeatingBase {
public:
    c_Heating() {
        this->p_sources = {&this->p_radiogenic_source, &this->p_tidal_source, &this->p_prescribed_source};
    }
    ~c_Heating() override = default;

    // The source list points into this object, so it is not copied.
    c_Heating(const c_Heating&) = delete;
    c_Heating& operator=(const c_Heating&) = delete;

    // Prepare every source for a solve of the world in `state`, and record which layers are heated and the
    // length and density units the solve runs in (one for an SI solve). The layers start at their own radii and
    // current masses.
    void update_sources(const c_WorldState& state, double length_scale, double density_scale) {
        this->p_length_scale = length_scale;
        this->p_density_scale = density_scale;
        this->p_use_heating_bylayer.clear();
        std::vector<c_HeatedLayer> layers;
        if (state.layers_ptr != nullptr) {
            for (const auto& layer_uptr : *state.layers_ptr) {
                this->p_use_heating_bylayer.push_back(layer_uptr->get_use_heating() ? 1 : 0);
                c_HeatedLayer heated;
                heated.radius_inner = layer_uptr->get_radius_inner();
                heated.radius_outer = layer_uptr->get_radius_outer();
                heated.mass         = layer_uptr->get_mass();
                layers.push_back(heated);
            }
        }
        for (c_HeatSourceBase* source_ptr : this->p_sources) {
            source_ptr->solve_radial_heat(state);
            source_ptr->update_layers(layers);
        }
    }

    // The layers' radii and masses as the last pass of a solve left them.
    void update_layers(const std::vector<c_HeatedLayer>& layers) {
        for (c_HeatSourceBase* source_ptr : this->p_sources) { source_ptr->update_layers(layers); }
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
        if (!this->p_heats(layer_index)) { return 0.0; }
        double heating = 0.0;
        for (const c_HeatSourceBase* source_ptr : this->p_sources) {
            heating += source_ptr->calc_heating(layer_index, radius, density);
        }
        return heating;
    }

    // The same from one source alone.
    double calc_source_heating(
            c_HeatSourceKind kind,
            std::size_t layer_index,
            double radius,
            double density) const noexcept {
        if (!this->p_heats(layer_index)) { return 0.0; }
        return this->p_sources[static_cast<std::size_t>(kind)]->calc_heating(layer_index, radius, density);
    }

    // The heat [W] one source generates in a layer of mass `mass` [kg]; zero in a layer that is not heated.
    double calc_layer_power(c_HeatSourceKind kind, std::size_t layer_index, double mass) const noexcept {
        if (!this->p_heats(layer_index)) { return 0.0; }
        const double power = this->p_sources[static_cast<std::size_t>(kind)]->calc_layer_power(layer_index, mass);
        return std::isfinite(power) ? power : 0.0;
    }

    // dL/dr in the units the solve runs in (see c_EOSHeatingBase). The heat flow itself stays in Watts.
    double calc_heat_flow_gradient(std::size_t layer_index, double radius, double density) const noexcept override {
        const double radius_si = radius * this->p_length_scale;
        return 4.0 * TidalPyConstants::d_PI * radius_si * radius_si * this->p_length_scale
            * this->calc_heating(layer_index, radius_si, density * this->p_density_scale);
    }

    const c_RadiogenicHeatSource& get_radiogenic_source() const noexcept { return this->p_radiogenic_source; }

protected:
    bool p_heats(std::size_t layer_index) const noexcept {
        return (layer_index < this->p_use_heating_bylayer.size()) && this->p_use_heating_bylayer[layer_index];
    }

    c_RadiogenicHeatSource p_radiogenic_source;
    c_TidalHeatSource p_tidal_source;
    c_PrescribedHeatSource p_prescribed_source;
    // Non-owning pointers to the sources above, in c_HeatSourceKind order.
    std::vector<c_HeatSourceBase*> p_sources;
    std::vector<char> p_use_heating_bylayer;
    double p_length_scale = 1.0;
    double p_density_scale = 1.0;
};

}  // namespace tidalpy
