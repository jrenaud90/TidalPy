#pragma once
/*
 * thermal_layout_.hpp: the temperature and heat-flow layout of a whole-planet EOS solve.
 *
 * Each layer carries its own temperature, and its attached cooling model says how heat moves inside it. This
 * header turns those into the radial segments the structure solve integrates (c_EOSSegment) and into the heat
 * flow across every interface. All quantities are MKS except where a coefficient is marked as carrying the unit
 * conversion of the solve.
 *
 * The thermal network is a chain of resistances. A layer's temperature is an input, so the heat flow between two
 * layers follows from their temperatures and the resistance of the material between them:
 *
 *     R = (1 / (4 pi k)) (1/r_inner - 1/r_outer)     [K W-1], a conducting spherical shell
 *     L = (T_lower - T_node) / R_lower               [W], the same through both sides of an interface
 *
 * The node temperature between two layers is the series-resistance temperature of their two facing resistances.
 * A layer with no cooling model is isothermal and perfectly conducting: it pins the node to its own temperature.
 * The center carries no flow (a regular solution), and the surface node is the world's surface temperature.
 * Where neither side of an interface has a resistance (an isothermal layer under another, or under the surface),
 * nothing holds a contrast across it, so the lower layer stores nothing: it passes on the heat that enters it and
 * the heat it generates, and the node keeps its temperature.
 *
 * A layer with no temperature of its own (one whose temperature is not a positive number, such as the 0 K
 * default) takes no part: its interfaces carry no flow, so it is neither a heat sink nor a
 * source for its neighbors, and each neighbor keeps its own temperature at the shared interface.
 *
 * Each layer's cooling model builds that layer's profile (c_CoolingBase::build_profile, Cooling/cooling_.hpp): which
 * stretches conduct and which are adiabatic, their resistances, the temperature at the base of a convecting interior,
 * and its Rayleigh and Nusselt numbers. It reads the solved structure and the layer's material through a probe
 * (c_LayerProbe below), so it chooses where to evaluate them. In short:
 *   off / none   one isothermal segment (so does a layer with no cooling model).
 *   conduction   two conducting halves, meeting at the layer's mid-radius where its temperature applies.
 *   convection   a conducting boundary layer at the base and the top, from the model's Nusselt scaling, around an
 *                adiabatic interior. The layer's temperature applies at the top of that interior, under the upper
 *                boundary layer: the upper-mantle temperature of parameterized convection (Stevenson et al. 1983;
 *                Schubert et al. 2001), where the model also takes its viscosity. The adiabat warms downward from
 *                it, so the base of the interior sits at T exp(integral of alpha g / c_p dr). A layer whose base
 *                carries no heat (the innermost layer, or one above a layer outside the network) has no boundary
 *                layer there: its interior reaches down to its base.
 *
 * The heat flow steps between the base and the top of a convecting interior: the difference is the heat the
 * lumped interior stores or releases, which is what makes its temperature evolve.
 *
 * Heat generated inside a conducting stretch (c_Heating) changes both of its ends. With H(r) the heat generated
 * between the base of the stretch and r, the flow leaving its top is the flow entering plus H(top), and the
 * temperature drop across it is
 *
 *     T_base - T_top = L_base R + drop,     drop = integral of H(r) / (4 pi r^2 k) dr
 *
 * Both follow from the heating and the solved density alone, so they are found by quadrature against the last
 * structure and the network stays a chain of resistances between effective temperatures.
 *
 * A material that melts at one temperature carries no latent heat in its heat capacity (there is no range to spread
 * it over), so a boundary between the solid and liquid zones of a layer carries it instead: as the layer's temperature
 * changes the boundary moves and melts or freezes mass, which adds a latent-heat capacity to the layer's temperature
 * rate (a Stefan condition, c_zone_boundary_latent_capacity).
 *
 * References
 * ----------
 * - Turcotte and Schubert (2002), Geodynamics: boundary-layer convection and the Nusselt scaling.
 * - Solomatov (1995): the stagnant-lid regime of strongly temperature-dependent viscosity.
 * - Stevenson, Spohn, and Schubert (1983); Schubert, Turcotte, and Olson (2001): the upper-mantle temperature of
 *   parameterized convection, at the top of the adiabatic interior.
 */

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "constants_.hpp"
#include "ode_.hpp"             // c_EOSSegment, c_TemperatureKind
#include "eos_solution_.hpp"    // c_EOSSolution
#include "../layers/layer_.hpp"
#include "../../Cooling/cooling_base_.hpp"
#include "../../Utilities/math/quadrature_.hpp"   // c_gauss_legendre_nodes
#include "heating_.hpp"                               // c_Heating

namespace tidalpy {

// Gauss-Legendre nodes for the integrals over one stretch of a layer: the heating, which is smooth geometry times a
// radiogenic heating that follows the density, and the gravity along an adiabat. This is far past what a layer's
// profile needs.
inline constexpr int d_HEATING_QUADRATURE_NODES = 16;

// The stretch quadrature rule on [-1, 1], built once (thread-safe static initialization).
struct c_StretchQuadrature {
    std::vector<double> nodes;
    std::vector<double> weights;
    c_StretchQuadrature() { c_gauss_legendre_nodes(d_HEATING_QUADRATURE_NODES, this->nodes, this->weights); }
};

inline const c_StretchQuadrature& c_stretch_quadrature() {
    static const c_StretchQuadrature rule;
    return rule;
}

// c_LayerThermal: the thermal description of one layer during a solve. Nothing here is stored on the layer; the
// solve builds it, iterates on it, and reports what the caller asked for.
struct c_LayerThermal {
    c_TemperatureKind kind = c_TemperatureKind::Isothermal;

    // The layer's radii [m] in the current pass: its own, or where the last pass put it when a layer below holds its
    // mass (and so moves every layer above it).
    double radius_inner = 0.0;
    double radius_outer = 0.0;

    // False for a layer with no temperature of its own: its interfaces are insulating (see the header comment).
    bool in_network = true;
    // True when the cooling model gave no usable boundary-layer thickness (a NaN viscosity, say), so the
    // boundary layers fell back to the largest share of the layer they may take.
    bool boundary_fallback = false;

    double temperature      = 0.0;  // [K] the layer's own temperature, an input of the solve
    double top_temperature  = 0.0;  // [K] at the top of its interior: the layer's own temperature
    double base_temperature = 0.0;  // [K] at the base of its interior (a convecting interior is warmer at its base)
    double node_temperature = 0.0;  // [K] at the interface above this layer

    double heat_flow_in  = 0.0;   // [W] entering the base
    double heat_flow_out = 0.0;   // [W] leaving the top

    double boundary_thickness = 0.0;   // [m] conducting boundary layer of a convecting layer
    double resistance_bottom  = 0.0;   // [K W-1] base to the layer's own temperature (zero with no base boundary layer)
    double resistance_top     = 0.0;   // [K W-1] the layer's own temperature to its top

    double rayleigh_number = 0.0;
    double nusselt_number  = 1.0;
    // True when a convecting interior was liquid at its reference point and took the liquid scaling.
    bool   magma_ocean     = false;
    // Where a convecting layer's model evaluated its viscosity: the top of its interior, at the layer's temperature.
    double reference_pressure      = TidalPyConstants::d_NAN;   // [Pa]
    double reference_viscosity     = TidalPyConstants::d_NAN;   // [Pa s]
    double reference_melt_fraction = TidalPyConstants::d_NAN;   // [m3 m-3]

    // Heat generated inside the layer [W], and the part of it generated inside the conducting stretch below and
    // above the layer's own temperature, with the temperature drop [K] each part adds across its stretch.
    double heating             = 0.0;
    double heating_bottom      = 0.0;
    double heating_top         = 0.0;
    double heating_drop_bottom = 0.0;
    double heating_drop_top    = 0.0;

    // The layer's material at its own temperature: at zero pressure before a solve, then at its mid-radius pressure.
    double conductivity      = TidalPyConstants::d_NAN;   // k     [W m-1 K-1]
    double thermal_expansion = TidalPyConstants::d_NAN;   // alpha [K-1]
    double heat_capacity     = TidalPyConstants::d_NAN;   // c_p   [J kg-1 K-1], latent heat included

    // The latent heat the boundaries between the layer's solid and liquid zones absorb per kelvin of the layer's
    // temperature [J K-1], where its material melts at one temperature (c_zone_boundary_latent_capacity).
    double latent_capacity = 0.0;
};

// True when no heat crosses a layer's base: the innermost layer (the center carries no flow) or one above a layer
// outside the network. A thermal boundary layer forms only where heat crosses a boundary, so a convecting layer
// has none there.
inline bool c_base_is_insulated(const std::vector<c_LayerThermal>& thermal_vec, std::size_t layer_index) noexcept {
    return (layer_index == 0) || !thermal_vec[layer_index - 1].in_network;
}

// The temperature kind a layer's cooling model gives its layer. A layer without one is isothermal.
inline c_TemperatureKind c_layer_temperature_kind(const c_Layer* layer) noexcept {
    const c_CoolingBase* cooling_model = layer->get_cooling_model();
    return (cooling_model == nullptr) ? c_TemperatureKind::Isothermal : cooling_model->get_temperature_kind();
}

// The layer's material at a point and its own lumped temperature [K], for the thermal network.
inline void c_layer_thermal_state(
        const c_Layer* layer,
        double pressure,
        double temperature,
        double radius,
        c_MaterialState& state) noexcept {
    c_ThermoPoint point;
    point.pressure    = pressure;
    point.temperature = temperature;
    point.radius      = radius;
    layer->calc_state(point, state);
}

// Build the per-layer thermal description from the layers themselves: the temperature each carries, the kind its
// cooling model asks for, and its thermal material properties at zero pressure. A finite temperature_override
// replaces every layer's own temperature. The resistances and flows are filled in later, against a solved
// structure.
inline void c_init_layer_thermal(
        const std::vector<std::unique_ptr<c_Layer>>& layers,
        std::vector<c_LayerThermal>& out,
        double temperature_override = TidalPyConstants::d_NAN) {
    const std::size_t n_layers = layers.size();
    out.assign(n_layers, c_LayerThermal());
    for (std::size_t layer_i = 0; layer_i < n_layers; ++layer_i) {
        const c_Layer* layer = layers[layer_i].get();
        c_LayerThermal& thermal  = out[layer_i];

        thermal.radius_inner = layer->get_radius_inner();
        thermal.radius_outer = layer->get_radius_outer();
        thermal.temperature = std::isfinite(temperature_override) ? temperature_override : layer->get_temperature();
        thermal.in_network  = std::isfinite(thermal.temperature) && (thermal.temperature > 0.0);
        thermal.top_temperature  = thermal.temperature;
        thermal.base_temperature = thermal.temperature;
        thermal.node_temperature = thermal.temperature;
        thermal.kind             = thermal.in_network ? c_layer_temperature_kind(layer) : c_TemperatureKind::Isothermal;

        // The thermal properties belong to the material, so any layer with one has them.
        if (layer->get_material_set()) {
            c_MaterialState state;
            c_layer_thermal_state(layer, 0.0, thermal.temperature, layer->get_radius_mid(), state);
            thermal.conductivity      = state.thermal_conductivity;
            thermal.thermal_expansion = state.thermal_expansion;
            thermal.heat_capacity     = state.heat_capacity;
        }
        if (!(thermal.conductivity > TidalPyConstants::d_EPS)) {
            // Without a conductivity there is no gradient to integrate.
            thermal.kind = c_TemperatureKind::Isothermal;
        }
    }
}

// True when the world carries a temperature contrast, so the solve has a profile to integrate. Everything at one
// temperature is isothermal whatever the cooling models say, and the solve keeps its four structure variables.
inline bool c_thermal_contrast_present(
        const std::vector<c_LayerThermal>& thermal_vec,
        double surface_temperature) noexcept {
    // Only the layers in the network count: a layer with no temperature exchanges no heat, so it sets up no
    // contrast however far its placeholder temperature is from its neighbors'.
    const c_LayerThermal* first = nullptr;
    for (const c_LayerThermal& thermal : thermal_vec) {
        if (!thermal.in_network) { continue; }
        if (first == nullptr) { first = &thermal; continue; }
        if (!c_isclose(thermal.temperature, first->temperature, 1.0e-12, 1.0e-12)) { return true; }
    }
    if (first == nullptr) { return false; }
    const c_LayerThermal& top = thermal_vec.back();
    if (top.in_network && std::isfinite(surface_temperature)
        && !c_isclose(surface_temperature, top.temperature, 1.0e-12, 1.0e-12)) {
        return true;
    }
    return false;
}

// Heat generated between two radii of a layer [W], and the temperature drop [K] it adds across that stretch when
// the stretch conducts (zero for a conductivity that is not positive). With the order of integration swapped,
//     drop = integral of 4 pi s^2 h(s) R(s, r_top) ds,     R(s, r_top) = (1/s - 1/r_top) / (4 pi k),
// so both come from one pass over the same nodes. The density is read from the solved structure.
inline void c_stretch_heating(
        const c_EOSSolution& solution,
        const c_Heating& heating,
        std::size_t layer_index,
        double radius_lower,
        double radius_upper,
        double conductivity,
        double& heat_out,
        double& drop_out) {
    heat_out = 0.0;
    drop_out = 0.0;
    if (!(radius_upper > radius_lower)) { return; }

    const std::vector<double>& nodes   = c_stretch_quadrature().nodes;
    const std::vector<double>& weights = c_stretch_quadrature().weights;
    const double half_width = 0.5 * (radius_upper - radius_lower);
    const double midpoint   = 0.5 * (radius_upper + radius_lower);
    const bool conducts     = (conductivity > TidalPyConstants::d_EPS);
    double state[C_EOS_DY_VALUES];
    for (std::size_t node_i = 0; node_i < nodes.size(); ++node_i) {
        const double radius = midpoint + half_width * nodes[node_i];
        solution.call_si(layer_index, radius, state);
        const double shell_heat = 4.0 * TidalPyConstants::d_PI * radius * radius
            * heating.calc_heating(layer_index, radius, state[C_EOS_DENSITY_INDEX]) * half_width * weights[node_i];
        if (!std::isfinite(shell_heat)) { continue; }
        heat_out += shell_heat;
        if (conducts) {
            drop_out += shell_heat * (1.0 / radius - 1.0 / radius_upper)
                / (4.0 * TidalPyConstants::d_PI * conductivity);
        }
    }
}

// The shear modulus [Pa] at or below which a material behaves as a liquid: the radial solver's minimum_solid_rigidity
// times the rigidity scale rho g R = 3 G M^2 / (4 pi R^4) of a planet of mass M [kg] and radius R [m] (bulk density,
// surface gravity, radius). Zero when the config or the planet gives no scale.
inline double c_liquid_shear_threshold(double mass, double radius, double G) noexcept {
    const double rigidity_scale = (radius > TidalPyConstants::d_EPS)
        ? 3.0 * G * mass * mass / (4.0 * TidalPyConstants::d_PI * radius * radius * radius * radius) : 0.0;
    const double min_rigidity = (tidalpy_config_ptr != nullptr)
        ? tidalpy_config_ptr->d_MIN_SOLID_RIGIDITY : TidalPyConstants::d_NAN;
    const double threshold = min_rigidity * rigidity_scale;
    return (std::isfinite(threshold) && (threshold > 0.0)) ? threshold : 0.0;
}

// The threshold a solve took its zones with (c_EOSSolution::liquid_shear_threshold): the same one, from the world's
// stated mass and radius, so the Love solve and the thermal network agree with the zones the integration found.
inline double c_liquid_shear_threshold(const c_EOSSolution& solution) noexcept {
    return solution.liquid_shear_threshold;
}

// The thermal network's view of one layer for its cooling model: the solved structure at a radius, and the layer's
// material at a point (c_LayerThermalProbe, Cooling/cooling_base_.hpp). The material is liquid where the Love solve
// takes it as liquid: everywhere in a liquid layer, and where a layer that can change state is fully molten or has a
// post-melt shear modulus at or below liquid_shear (c_liquid_shear_threshold).
class c_LayerProbe final : public c_LayerThermalProbe {
public:
    c_LayerProbe(
            const c_EOSSolution& solution,
            const c_Layer& layer,
            std::size_t layer_index,
            double temperature,
            double liquid_shear) noexcept
        : p_solution(solution), p_layer(layer), p_layer_index(layer_index), p_temperature(temperature),
          p_liquid_shear(liquid_shear) {}

    void calc_structure(double radius, double& gravity, double& pressure) const override {
        double structure[C_EOS_Y_VALUES];
        this->p_solution.call_y_si(this->p_layer_index, radius, structure);
        gravity  = structure[C_EOS_GRAVITY_INDEX];
        pressure = structure[C_EOS_PRESSURE_INDEX];
    }

    void calc_transport_state(
            double pressure,
            double temperature,
            double radius,
            c_TransportState& out) const override {
        c_MaterialState state;
        c_layer_thermal_state(&this->p_layer, pressure, temperature, radius, state);
        out.density              = state.density;
        out.thermal_conductivity = state.thermal_conductivity;
        out.heat_capacity        = state.heat_capacity;
        out.thermal_expansion    = state.thermal_expansion;
        out.shear_viscosity      = state.shear_viscosity;
        out.melt_fraction        = state.melt_fraction;
        out.is_liquid            = this->p_layer.get_is_liquid()
            || (this->p_layer.get_can_change_state()
                && ((state.phase == c_MaterialPhase::Liquid) || (state.shear_modulus <= this->p_liquid_shear)));
    }

    // The integral of alpha g / c_p from radius_lower to radius_upper, using the solved gravity and pressure and the
    // layer's material at its own temperature (whose expansivity can fall with compression), as the integrated
    // adiabat does. Zero for an empty stretch or a layer without a material, which leaves the interior isothermal.
    double calc_adiabat_exponent(double radius_lower, double radius_upper) const override {
        if (!(radius_upper > radius_lower) || !this->p_layer.get_material_set()) { return 0.0; }
        const std::vector<double>& nodes   = c_stretch_quadrature().nodes;
        const std::vector<double>& weights = c_stretch_quadrature().weights;
        const double half_width = 0.5 * (radius_upper - radius_lower);
        const double midpoint   = 0.5 * (radius_upper + radius_lower);
        double integral = 0.0;
        double gravity  = 0.0;
        double pressure = 0.0;
        c_TransportState material;
        for (std::size_t node_i = 0; node_i < nodes.size(); ++node_i) {
            const double radius = midpoint + half_width * nodes[node_i];
            this->calc_structure(radius, gravity, pressure);
            this->calc_transport_state(pressure, this->p_temperature, radius, material);
            if (!(material.heat_capacity > TidalPyConstants::d_EPS)) { continue; }
            const double term = gravity * material.thermal_expansion / material.heat_capacity;
            if (std::isfinite(term)) { integral += term * half_width * weights[node_i]; }
        }
        return integral;
    }

private:
    const c_EOSSolution& p_solution;
    const c_Layer&       p_layer;
    std::size_t          p_layer_index;
    double               p_temperature;
    double               p_liquid_shear;
};

// Update the boundary layers, resistances, interface temperatures, and heat flows against a solved structure.
// Returns the largest relative change in the interface temperatures and flows, which is what the solve watches
// to decide it has converged.
inline double c_update_layer_thermal(
        const c_EOSSolution& solution,
        const std::vector<std::unique_ptr<c_Layer>>& layers,
        double surface_temperature,
        std::vector<c_LayerThermal>& thermal_vec,
        const c_Heating* heating_ptr = nullptr) {
    const std::size_t n_layers = layers.size();
    double largest_change = 0.0;
    const bool heated = (heating_ptr != nullptr) && heating_ptr->get_is_active();
    const double liquid_shear = c_liquid_shear_threshold(solution);

    // Each layer's cooling model builds its profile against this pass's structure, its neighbors, and its own last
    // pass (the base of its adiabat and its boundary layers lag by one).
    for (std::size_t layer_i = 0; layer_i < n_layers; ++layer_i) {
        const c_Layer* layer = layers[layer_i].get();
        c_LayerThermal& thermal  = thermal_vec[layer_i];
        thermal.top_temperature  = thermal.temperature;
        const c_CoolingBase* cooling_model = layer->get_cooling_model();
        if ((thermal.kind == c_TemperatureKind::Isothermal) || (cooling_model == nullptr)) {
            thermal.boundary_thickness = 0.0;
            thermal.resistance_bottom  = 0.0;
            thermal.resistance_top     = 0.0;
            thermal.boundary_fallback  = false;
            thermal.base_temperature   = thermal.temperature;
            continue;
        }

        c_LayerThermalContext context;
        context.radius_inner       = thermal.radius_inner;
        context.radius_outer       = thermal.radius_outer;
        context.temperature        = thermal.temperature;
        context.base_temperature   = thermal.base_temperature;
        context.boundary_thickness = thermal.boundary_thickness;
        context.insulated_base     = c_base_is_insulated(thermal_vec, layer_i);
        // A neighbor outside the network exchanges no heat, so there is no drop across that boundary layer.
        context.inner_temperature = ((layer_i > 0) && thermal_vec[layer_i - 1].in_network)
            ? thermal_vec[layer_i - 1].top_temperature
            : thermal.base_temperature;
        context.outer_temperature = thermal.temperature;
        if (layer_i + 1 < n_layers) {
            if (thermal_vec[layer_i + 1].in_network) {
                context.outer_temperature = thermal_vec[layer_i + 1].base_temperature;
            }
        } else if (std::isfinite(surface_temperature)) {
            context.outer_temperature = surface_temperature;
        }

        const c_LayerProbe probe(solution, *layer, layer_i, thermal.temperature, liquid_shear);
        c_LayerThermalProfile profile;
        cooling_model->build_profile(context, probe, profile);

        thermal.boundary_thickness      = profile.boundary_thickness;
        thermal.resistance_bottom       = profile.resistance_bottom;
        thermal.resistance_top          = profile.resistance_top;
        thermal.base_temperature        = profile.base_temperature;
        thermal.boundary_fallback       = profile.boundary_fallback;
        thermal.rayleigh_number         = profile.rayleigh_number;
        thermal.nusselt_number          = profile.nusselt_number;
        thermal.magma_ocean             = profile.magma_ocean;
        thermal.reference_pressure      = profile.reference_pressure;
        thermal.reference_viscosity     = profile.reference_viscosity;
        thermal.reference_melt_fraction = profile.reference_melt_fraction;
        // The material the profile used, where it gave one.
        if (std::isfinite(profile.conductivity))      { thermal.conductivity      = profile.conductivity; }
        if (std::isfinite(profile.heat_capacity))     { thermal.heat_capacity     = profile.heat_capacity; }
        if (std::isfinite(profile.thermal_expansion)) { thermal.thermal_expansion = profile.thermal_expansion; }
    }

    // Heat generated in each layer, and in the conducting stretches on either side of its own temperature.
    for (std::size_t layer_i = 0; layer_i < n_layers; ++layer_i) {
        c_LayerThermal& thermal = thermal_vec[layer_i];
        thermal.heating             = 0.0;
        thermal.heating_bottom      = 0.0;
        thermal.heating_top         = 0.0;
        thermal.heating_drop_bottom = 0.0;
        thermal.heating_drop_top    = 0.0;
        if (!heated) { continue; }

        const double radius_inner = thermal.radius_inner;
        const double radius_outer = thermal.radius_outer;
        double unused_drop = 0.0;
        c_stretch_heating(
            solution, *heating_ptr, layer_i, radius_inner, radius_outer, 0.0, thermal.heating, unused_drop);
        if (thermal.kind == c_TemperatureKind::Isothermal) { continue; }

        // Where the layer's own temperature applies: the mid-radius of a conducting layer, the two ends of the
        // interior of a convecting one (whose interior starts at its base when the base carries no heat).
        const bool convects = (thermal.kind == c_TemperatureKind::Adiabatic);
        const double boundary_bottom = c_base_is_insulated(thermal_vec, layer_i) ? 0.0 : thermal.boundary_thickness;
        const double bottom_end = convects
            ? (radius_inner + boundary_bottom) : 0.5 * (radius_inner + radius_outer);
        const double top_start = convects
            ? (radius_outer - thermal.boundary_thickness) : 0.5 * (radius_inner + radius_outer);
        c_stretch_heating(
            solution, *heating_ptr, layer_i, radius_inner, bottom_end, thermal.conductivity,
            thermal.heating_bottom, thermal.heating_drop_bottom);
        c_stretch_heating(
            solution, *heating_ptr, layer_i, top_start, radius_outer, thermal.conductivity,
            thermal.heating_top, thermal.heating_drop_top);
    }

    // Interface nodes, from the center outward. The flow through an interface is the same on both sides, so the
    // node sits where the two facing resistances balance; a zero resistance pins it to that layer. Heating
    // inside a stretch shifts the temperature its far end sees, which keeps the balance a resistance chain:
    //     lower, top stretch:     (T - drop + H R) - T_node = L_node R
    //     upper, bottom stretch:  T_node - (T + drop)       = L_node R
    for (std::size_t layer_i = 0; layer_i < n_layers; ++layer_i) {
        c_LayerThermal& lower = thermal_vec[layer_i];
        const bool at_surface = (layer_i + 1 == n_layers);

        const double resistance_lower = lower.resistance_top;
        const double temperature_lower =
            lower.top_temperature - lower.heating_drop_top + lower.heating_top * lower.resistance_top;
        double resistance_upper  = 0.0;
        double temperature_upper = 0.0;
        if (at_surface) {
            if (!std::isfinite(surface_temperature)) {
                // No surface temperature: nothing drives a flow out of the world.
                lower.node_temperature = temperature_lower;
                lower.heat_flow_out    = 0.0;
                continue;
            }
            temperature_upper = surface_temperature;
        } else {
            resistance_upper  = thermal_vec[layer_i + 1].resistance_bottom;
            temperature_upper =
                thermal_vec[layer_i + 1].base_temperature + thermal_vec[layer_i + 1].heating_drop_bottom;
        }

        double node = 0.0;
        double flow = 0.0;
        // Any positive resistance conducts. A resistance in K/W has no natural scale to compare with: the
        // millimetre boundary layers of a molten, vigorously convecting layer are about 1e-18 K/W, and must still
        // hold that layer's own temperature. c_shell_resistance gives exactly zero for no shell at all.
        const bool lower_conducts = (resistance_lower > 0.0);
        const bool upper_conducts = (resistance_upper > 0.0);
        const bool upper_in_network = at_surface || thermal_vec[layer_i + 1].in_network;
        if (!lower.in_network || !upper_in_network) {
            // A layer with no temperature closes the interface: no flow, and the side in the network keeps its
            // own temperature there.
            node = lower.in_network ? temperature_lower : temperature_upper;
            flow = 0.0;
        } else if (lower_conducts && upper_conducts) {
            const double conductance_lower = 1.0 / resistance_lower;
            const double conductance_upper = 1.0 / resistance_upper;
            node = (temperature_lower * conductance_lower + temperature_upper * conductance_upper)
                 / (conductance_lower + conductance_upper);
            flow = (temperature_lower - node) * conductance_lower;
        } else if (lower_conducts) {
            // The layer above is perfectly conducting, so it holds the interface at its own temperature.
            node = temperature_upper;
            flow = (temperature_lower - node) / resistance_lower;
        } else if (upper_conducts) {
            node = temperature_lower;
            flow = (node - temperature_upper) / resistance_upper;
        } else {
            // Neither side has a resistance, so nothing holds a contrast across the interface: the lower layer
            // stores nothing and passes on the heat entering it plus the heat it generates, which is the flow its
            // own profile reaches at its top. The node keeps the lower layer's temperature, where that profile ends.
            node = temperature_lower;
            flow = lower.heat_flow_in + lower.heating;
        }

        const double reference = std::fabs(node) + std::fabs(flow) + 1.0;
        const double change = (std::fabs(node - lower.node_temperature) + std::fabs(flow - lower.heat_flow_out))
                            / reference;
        if (change > largest_change) { largest_change = change; }

        lower.node_temperature = node;
        lower.heat_flow_out    = flow;
        if (!at_surface) { thermal_vec[layer_i + 1].heat_flow_in = flow; }
    }
    if (!thermal_vec.empty()) { thermal_vec.front().heat_flow_in = 0.0; }
    return largest_change;
}

// Turn the per-layer thermal description into the radial segments the solve integrates, with the radii in the units
// the solve runs in, on each layer's radii in the current pass. The gradient within a segment comes from the material
// at each point of the integration.
inline void c_build_thermal_segments(
        const std::vector<c_LayerThermal>& thermal_vec,
        bool integrate_temperature,
        double length_scale,
        std::vector<c_EOSSegment>& out) {
    const std::size_t n_layers = thermal_vec.size();
    out.clear();
    out.reserve(3 * n_layers);

    for (std::size_t layer_i = 0; layer_i < n_layers; ++layer_i) {
        const c_LayerThermal& thermal = thermal_vec[layer_i];
        const double radius_inner = thermal.radius_inner;
        const double radius_outer = thermal.radius_outer;

        c_EOSSegment segment;
        segment.layer_index      = layer_i;
        segment.upper_radius     = radius_outer / length_scale;
        segment.start_heat_flow  = thermal.heat_flow_in;

        const bool isothermal = (!integrate_temperature) || (thermal.kind == c_TemperatureKind::Isothermal);
        if (isothermal) {
            // One segment holding this layer's own temperature. It breaks the profile at its base, which is what
            // a perfectly conducting layer does to its neighbors.
            segment.temperature_kind  = c_TemperatureKind::Isothermal;
            segment.start_temperature = thermal.temperature;
            out.push_back(segment);
            continue;
        }

        // The innermost layer has no node below to continue from, and nor does one above a layer outside the
        // network (whose profile holds its own placeholder temperature), so its base starts where the stretch has
        // to start to reach the base of the layer's interior (its own temperature, for a conducting layer).
        const bool starts_fresh = c_base_is_insulated(thermal_vec, layer_i);
        const double start_temperature = starts_fresh
            ? (thermal.base_temperature + thermal.heat_flow_in * thermal.resistance_bottom
               + thermal.heating_drop_bottom)
            : TidalPyConstants::d_NAN;

        // A convecting layer whose base carries no heat has no boundary layer there.
        const double boundary = thermal.boundary_thickness;
        const double boundary_bottom = starts_fresh ? 0.0 : boundary;
        const bool has_interior = (thermal.kind == c_TemperatureKind::Adiabatic)
            && (boundary > 0.0)
            && (radius_outer - boundary > radius_inner + boundary_bottom);

        // The base of a convecting interior (where its adiabat starts) or the mid-radius of a conducting layer
        // splits the layer, so the two stretches around it carry their own heat flow.
        const double split_lower = has_interior
            ? (radius_inner + boundary_bottom) : 0.5 * (radius_inner + radius_outer);

        const bool has_lower_stretch = !(has_interior && starts_fresh);
        if (has_lower_stretch) {
            c_EOSSegment lower = segment;
            lower.temperature_kind  = c_TemperatureKind::Conductive;
            lower.start_temperature = start_temperature;
            lower.start_heat_flow   = thermal.heat_flow_in;
            lower.upper_radius      = split_lower / length_scale;
            out.push_back(lower);
        }

        if (has_interior) {
            c_EOSSegment interior = segment;
            interior.temperature_kind  = c_TemperatureKind::Adiabatic;
            // With no lower stretch the interior is the first segment of the layer and starts at its temperature.
            interior.start_temperature = has_lower_stretch ? TidalPyConstants::d_NAN : start_temperature;
            interior.start_heat_flow   = thermal.heat_flow_in + thermal.heating_bottom;
            interior.upper_radius      = (radius_outer - boundary) / length_scale;
            out.push_back(interior);
        }

        c_EOSSegment upper = segment;
        upper.temperature_kind  = c_TemperatureKind::Conductive;
        upper.start_temperature = TidalPyConstants::d_NAN;
        // The flow leaving the top is the flow entering the stretch plus the heat generated inside it.
        upper.start_heat_flow   = thermal.heat_flow_out - thermal.heating_top;
        upper.upper_radius      = radius_outer / length_scale;
        out.push_back(upper);
    }
}

// How far the temperature [K] at a radius [m] of a layer moves per kelvin of the layer's own temperature, with its
// neighbors' interface temperatures held: one in an isothermal layer, T(r) / T along an adiabatic interior (which
// scales with the temperature at its top), and in a conducting stretch the fraction of the way from the interface at
// its far end to the end the layer holds, scaled by how that end moves. A stretch whose far end is an insulated base
// moves with the layer.
inline double c_temperature_sensitivity(
        const std::vector<c_LayerThermal>& thermal_vec,
        std::size_t layer_index,
        double radius,
        double temperature) noexcept {
    const c_LayerThermal& thermal = thermal_vec[layer_index];
    const double layer_temperature = thermal.temperature;
    if ((thermal.kind == c_TemperatureKind::Isothermal) || !(layer_temperature > 0.0)) { return 1.0; }
    const bool insulated_base = c_base_is_insulated(thermal_vec, layer_index);
    const double node_below = insulated_base ? TidalPyConstants::d_NAN
        : thermal_vec[layer_index - 1].node_temperature;
    const double node_above = thermal.node_temperature;

    // The fraction of the way from a stretch's far end to the end the layer holds, times how that end moves.
    const auto stretch = [temperature](double far_end, double held_end, double held_scale) {
        const double span = held_end - far_end;
        if (!std::isfinite(far_end) || !(std::fabs(span) > TidalPyConstants::d_EPS)) { return held_scale; }
        return held_scale * (temperature - far_end) / span;
    };
    const double interior_scale = temperature / layer_temperature;
    double sensitivity = interior_scale;
    if (thermal.kind == c_TemperatureKind::Conductive) {
        const double radius_mid = 0.5 * (thermal.radius_inner + thermal.radius_outer);
        sensitivity = (radius >= radius_mid)
            ? stretch(node_above, layer_temperature, 1.0)
            : (insulated_base ? interior_scale : stretch(node_below, layer_temperature, 1.0));
    } else {
        const double boundary_bottom = insulated_base ? 0.0 : thermal.boundary_thickness;
        if (radius > thermal.radius_outer - thermal.boundary_thickness) {
            sensitivity = stretch(node_above, layer_temperature, 1.0);
        } else if (radius < thermal.radius_inner + boundary_bottom) {
            sensitivity = stretch(node_below, thermal.base_temperature, thermal.base_temperature / layer_temperature);
        }
    }
    return std::isfinite(sensitivity) ? std::max(sensitivity, 0.0) : 1.0;
}

// The latent heat [J K-1] a boundary between a solid and a liquid zone of one layer absorbs per kelvin of the layer's
// temperature, where the layer's material melts at one temperature T_m(P) (a step; a material with a melting range
// spreads its latent heat over its heat capacity instead). The boundary sits where G(r) = T(r) - T_m(P(r)) = 0, and
// moves by dr_b/dT = -(dG/dT) / (dG/dr), so the liquid mass changes by 4 pi r_b^2 rho |dr_b/dT| per kelvin:
//     C_latent = L 4 pi r_b^2 rho_solid S / |dG/dr|,
// with S the temperature sensitivity at the boundary (c_temperature_sensitivity) and dG/dr a central difference
// across the boundary. Zero where the material melts over a range or the boundary does not move.
//
// Assumptions: the neighbors' interface temperatures hold while the layer's own temperature changes, and the melt
// that forms takes the solid's density at the boundary.
struct c_ZoneBoundary {
    std::size_t layer_index = 0;
    double radius = 0.0;        // [m]
    double half_width = 0.0;    // of the central difference across it, inside both zones [m]
    bool solid_below = true;    // the solid zone is the one beneath it
};

inline double c_zone_boundary_latent_capacity(
        const c_EOSSolution& solution,
        const c_Material& material,
        const c_MaterialSwitches& switches,
        const std::vector<c_LayerThermal>& thermal_vec,
        const c_ZoneBoundary& boundary) {
    const std::size_t layer_index = boundary.layer_index;
    const double boundary_radius  = boundary.radius;
    const double half_width       = boundary.half_width;
    const bool solid_below        = boundary.solid_below;
    const double latent_heat = material.get_latent_heat();
    if (!(latent_heat > 0.0) || !(half_width > 0.0)) { return 0.0; }
    double state[C_EOS_DY_VALUES];
    solution.call_si(layer_index, boundary_radius, state);
    double solidus  = TidalPyConstants::d_NAN;
    double liquidus = TidalPyConstants::d_NAN;
    material.calc_melting_range(state[C_EOS_PRESSURE_INDEX], switches, solidus, liquidus);
    if (!std::isfinite(solidus) || (liquidus - solidus > TidalPyConstants::d_EPS)) { return 0.0; }

    // G on each side of the boundary.
    const auto margin = [&](double radius) {
        double side[C_EOS_DY_VALUES];
        solution.call_si(layer_index, radius, side);
        double side_solidus  = TidalPyConstants::d_NAN;
        double side_liquidus = TidalPyConstants::d_NAN;
        material.calc_melting_range(side[C_EOS_PRESSURE_INDEX], switches, side_solidus, side_liquidus);
        return side[C_EOS_TEMPERATURE_INDEX] - side_solidus;
    };
    const double margin_slope = (margin(boundary_radius + half_width) - margin(boundary_radius - half_width))
        / (2.0 * half_width);   // [K m-1]
    if (!(std::fabs(margin_slope) > TidalPyConstants::d_EPS)) { return 0.0; }

    // The solid's density beside the boundary, where the material is unambiguously solid.
    double solid_state[C_EOS_DY_VALUES];
    const double solid_radius = solid_below ? boundary_radius - half_width : boundary_radius + half_width;
    solution.call_si(layer_index, solid_radius, solid_state);
    const double sensitivity = c_temperature_sensitivity(
        thermal_vec, layer_index, boundary_radius, state[C_EOS_TEMPERATURE_INDEX]);
    const double capacity = latent_heat * 4.0 * TidalPyConstants::d_PI * boundary_radius * boundary_radius
        * solid_state[C_EOS_DENSITY_INDEX] * sensitivity / std::fabs(margin_slope);
    return std::isfinite(capacity) ? capacity : 0.0;
}

}  // namespace tidalpy
