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
 * What each cooling model makes of a layer:
 *   off / none   one isothermal segment.
 *   conduction   two conducting halves, meeting at the layer's mid-radius where its temperature applies.
 *   convection   a conducting boundary layer at the base and the top, from the model's Nusselt scaling, around an
 *                adiabatic interior. The layer's temperature applies at the top of that interior, under the upper
 *                boundary layer: the upper-mantle temperature of parameterized convection (Stevenson et al. 1983;
 *                Schubert et al. 2001). The adiabat warms downward from it, so the base of the interior sits at
 *                T exp(integral of alpha g / c_p dr). A layer whose base carries no heat (the innermost layer, or
 *                one above a layer outside the network) has no boundary layer there: its interior reaches down to
 *                its base.
 *
 * The cooling model's boundary-layer thickness, D / Nu, is the conducting thickness that carries its flux,
 * Nu k dT / D, across the whole drop dT. Two boundary layers each carry that flux across their own share of the
 * drop, so each is half as thick; a layer with only the upper one keeps the whole thickness. At Nu = 1 (a
 * sub-critical layer) the two would be the halves of a conducting layer; d_MAX_BOUNDARY_FRACTION keeps each a
 * little thinner so that an interior remains.
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

// Largest share of a layer's thickness one conducting boundary layer may take, so a convecting layer keeps an
// interior to be adiabatic in.
inline constexpr double d_MAX_BOUNDARY_FRACTION = 0.4;

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
};

// Thermal resistance [K W-1] of a conducting spherical shell. Zero for a degenerate shell or conductivity.
inline double c_shell_resistance(double radius_inner, double radius_outer, double conductivity) noexcept {
    if (!(conductivity > TidalPyConstants::d_EPS) || !(radius_outer > radius_inner)) { return 0.0; }
    const double inner = (radius_inner > TidalPyConstants::d_EPS) ? (1.0 / radius_inner) : (1.0 / radius_outer);
    return (inner - 1.0 / radius_outer) / (4.0 * TidalPyConstants::d_PI * conductivity);
}

// True when no heat crosses a layer's base: the innermost layer (the center carries no flow) or one above a layer
// outside the network. A thermal boundary layer forms only where heat crosses a boundary, so a convecting layer
// has none there.
inline bool c_base_is_insulated(const std::vector<c_LayerThermal>& thermal_vec, std::size_t layer_index) noexcept {
    return (layer_index == 0) || !thermal_vec[layer_index - 1].in_network;
}

// The temperature kind a layer's cooling model asks for. A layer without one is isothermal.
inline c_TemperatureKind c_layer_temperature_kind(const c_Layer* layer) noexcept {
    const c_CoolingBase* cooling_model = layer->get_cooling_model();
    if (cooling_model == nullptr) { return c_TemperatureKind::Isothermal; }
    switch (cooling_model->get_model_type()) {
        case c_CoolingModel::Conduction: return c_TemperatureKind::Conductive;
        case c_CoolingModel::Convection: return c_TemperatureKind::Adiabatic;
        case c_CoolingModel::Off:        return c_TemperatureKind::Isothermal;
    }
    return c_TemperatureKind::Isothermal;
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

// Exponent of the adiabatic warming across a convecting interior, the integral of alpha g / c_p from its base to its
// top, using the solved gravity and pressure and the layer's material at its own temperature (whose expansivity can
// fall with compression), as the integrated adiabat does. With dT/dr = -alpha g T / c_p the base sits at
// T_top exp(exponent). Zero for an empty stretch or a layer without a material, which leaves the interior
// isothermal.
inline double c_adiabat_exponent(
        const c_EOSSolution& solution,
        const c_Layer* layer,
        std::size_t layer_index,
        double radius_lower,
        double radius_upper,
        const c_LayerThermal& thermal) {
    if (!(radius_upper > radius_lower) || !layer->get_material_set()) {
        return 0.0;
    }
    const std::vector<double>& nodes   = c_stretch_quadrature().nodes;
    const std::vector<double>& weights = c_stretch_quadrature().weights;
    const double half_width = 0.5 * (radius_upper - radius_lower);
    const double midpoint   = 0.5 * (radius_upper + radius_lower);
    double integral = 0.0;
    double structure[C_EOS_Y_VALUES];
    c_MaterialState material;
    for (std::size_t node_i = 0; node_i < nodes.size(); ++node_i) {
        const double radius = midpoint + half_width * nodes[node_i];
        solution.call_y_si(layer_index, radius, structure);
        c_layer_thermal_state(layer, structure[C_EOS_PRESSURE_INDEX], thermal.temperature, radius, material);
        if (!(material.heat_capacity > TidalPyConstants::d_EPS)) { continue; }
        const double term = structure[C_EOS_GRAVITY_INDEX] * material.thermal_expansion / material.heat_capacity;
        if (std::isfinite(term)) { integral += term * half_width * weights[node_i]; }
    }
    return integral;
}

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

    for (std::size_t layer_i = 0; layer_i < n_layers; ++layer_i) {
        const c_Layer* layer = layers[layer_i].get();
        c_LayerThermal& thermal  = thermal_vec[layer_i];
        const double radius_inner = layer->get_radius_inner();
        const double radius_outer = layer->get_radius_outer();
        const double thickness    = radius_outer - radius_inner;
        const double radius_mid   = 0.5 * (radius_inner + radius_outer);

        thermal.boundary_thickness = 0.0;
        thermal.resistance_bottom  = 0.0;
        thermal.resistance_top     = 0.0;
        thermal.boundary_fallback  = false;
        thermal.top_temperature    = thermal.temperature;
        if (thermal.kind == c_TemperatureKind::Isothermal) {
            thermal.base_temperature = thermal.temperature;
            continue;
        }

        if (thermal.kind == c_TemperatureKind::Conductive) {
            thermal.resistance_bottom = c_shell_resistance(radius_inner, radius_mid, thermal.conductivity);
            thermal.resistance_top    = c_shell_resistance(radius_mid, radius_outer, thermal.conductivity);
            thermal.base_temperature  = thermal.temperature;
            continue;
        }

        // Convecting layer: the cooling model sizes both boundary layers from the temperature drop across the
        // layer, the local state, and the viscosity at the layer's own temperature. The drop is the sum of the two
        // boundary layers' drops: from the top of the layer below (its own temperature) to the base of this layer's
        // adiabat, and from this layer's temperature (the top of that adiabat) to the base of the layer above or the
        // surface. A base that carries no heat (the center, or a layer outside the network) has no boundary layer,
        // so that layer has only the upper drop. The adiabat's base comes from the last pass, so it lags by one.
        double structure[C_EOS_Y_VALUES];
        solution.call_y_si(layer_i, radius_mid, structure);
        const double gravity  = structure[C_EOS_GRAVITY_INDEX];
        const double pressure = structure[C_EOS_PRESSURE_INDEX];

        const c_CoolingBase* cooling_model = layer->get_cooling_model();
        // The material at the layer's own (lumped) temperature, which is not the local temperature of the
        // profile at this radius, so it is asked of the material rather than read from the solution.
        c_MaterialState material;
        c_layer_thermal_state(layer, pressure, thermal.temperature, radius_mid, material);
        if (std::isfinite(material.thermal_conductivity)) { thermal.conductivity = material.thermal_conductivity; }
        if (std::isfinite(material.heat_capacity))        { thermal.heat_capacity = material.heat_capacity; }
        if (std::isfinite(material.thermal_expansion))    { thermal.thermal_expansion = material.thermal_expansion; }

        // A neighbor outside the network exchanges no heat, so there is no drop across that boundary layer.
        const double inner_temperature = ((layer_i > 0) && thermal_vec[layer_i - 1].in_network)
            ? thermal_vec[layer_i - 1].top_temperature
            : thermal.base_temperature;
        double outer_temperature = thermal.temperature;
        if (layer_i + 1 < n_layers) {
            if (thermal_vec[layer_i + 1].in_network) { outer_temperature = thermal_vec[layer_i + 1].base_temperature; }
        } else if (std::isfinite(surface_temperature)) {
            outer_temperature = surface_temperature;
        }

        c_CoolingInputs cooling_inputs;
        cooling_inputs.delta_temp = std::fabs(inner_temperature - thermal.base_temperature)
                                  + std::fabs(thermal.temperature - outer_temperature);
        cooling_inputs.thickness  = thickness;
        cooling_inputs.gravity    = gravity;
        cooling_inputs.density    = std::isfinite(material.density) ? material.density : layer->get_density_bulk();
        cooling_inputs.viscosity  = material.shear_viscosity;
        cooling_inputs.thermal_conductivity = thermal.conductivity;
        // The diffusivity uses the density the material has here, not a separate reference density.
        cooling_inputs.thermal_diffusivity  =
            thermal.conductivity / (cooling_inputs.density * thermal.heat_capacity);
        // The expansivity at the layer's own state, which can fall with compression.
        cooling_inputs.thermal_expansion    = thermal.thermal_expansion;
        const c_CoolingResult cooling_result = cooling_model->calc_cooling(cooling_inputs);

        thermal.rayleigh_number = cooling_result.rayleigh_number;
        thermal.nusselt_number  = cooling_result.nusselt_number;
        // The model's thickness carries its flux across the whole drop; two boundary layers split the drop at that
        // flux, so each takes half of it (see the header comment).
        const bool insulated_base = c_base_is_insulated(thermal_vec, layer_i);
        double boundary = insulated_base ? cooling_result.blt : 0.5 * cooling_result.blt;
        if (!(boundary > 0.0) || !std::isfinite(boundary)) {
            // The solve reports it; the layer's viscosity is the usual cause.
            thermal.boundary_fallback = true;
            boundary = d_MAX_BOUNDARY_FRACTION * thickness;
        }
        if (boundary > d_MAX_BOUNDARY_FRACTION * thickness) { boundary = d_MAX_BOUNDARY_FRACTION * thickness; }
        thermal.boundary_thickness = boundary;

        thermal.resistance_bottom = insulated_base ? 0.0 : c_shell_resistance(
            radius_inner, radius_inner + boundary, thermal.conductivity);
        thermal.resistance_top = c_shell_resistance(
            radius_outer - boundary, radius_outer, thermal.conductivity);

        // The layer's temperature holds at the top of the interior, and the adiabat warms downward from it to the
        // interior's base, along the gravity of the solved structure.
        const double interior_base = radius_inner + (insulated_base ? 0.0 : boundary);
        thermal.base_temperature = thermal.temperature * std::exp(c_adiabat_exponent(
            solution, layer, layer_i, interior_base, radius_outer - boundary, thermal));
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

        const double radius_inner = layers[layer_i]->get_radius_inner();
        const double radius_outer = layers[layer_i]->get_radius_outer();
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
// the solve runs in. The gradient within a segment comes from the material at each point of the integration.
inline void c_build_thermal_segments(
        const std::vector<c_LayerThermal>& thermal_vec,
        const std::vector<std::unique_ptr<c_Layer>>& layers,
        bool integrate_temperature,
        double length_scale,
        std::vector<c_EOSSegment>& out) {
    const std::size_t n_layers = layers.size();
    out.clear();
    out.reserve(3 * n_layers);

    for (std::size_t layer_i = 0; layer_i < n_layers; ++layer_i) {
        const c_LayerThermal& thermal = thermal_vec[layer_i];
        const double radius_inner = layers[layer_i]->get_radius_inner();
        const double radius_outer = layers[layer_i]->get_radius_outer();

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

}  // namespace tidalpy
