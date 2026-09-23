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
 *
 * A layer with no temperature of its own (a geometry-only layer, or one whose temperature is not a positive
 * number, such as the 0 K default) takes no part: its interfaces carry no flow, so it is neither a heat sink nor a
 * source for its neighbors, and each neighbor keeps its own temperature at the shared interface.
 *
 * What each cooling model makes of a layer:
 *   off / none   one isothermal segment.
 *   conduction   two conducting halves, meeting at the layer's mid-radius where its temperature applies.
 *   convection   a conducting boundary layer at the base and the top, from the model's Nusselt scaling, around an
 *                adiabatic interior. The layer's temperature applies at the base of that interior.
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
 */

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "constants_.hpp"
#include "ode_.hpp"             // c_EOSSegment, c_TemperatureKind
#include "eos_solution_.hpp"    // c_EOSSolution
#include "../layers/base_.hpp"
#include "../layers/physics_.hpp"
#include "../layers/solidliquid_.hpp"
#include "../../cooling_x/cooling_base_.hpp"
#include "../../Utilities_x/math_x/quadrature_.hpp"   // c_gauss_legendre_nodes
#include "heating_.hpp"                               // c_Heating

namespace tidalpy {

// Largest share of a layer's thickness one conducting boundary layer may take, so a convecting layer keeps an
// interior to be adiabatic in.
inline constexpr double d_MAX_BOUNDARY_FRACTION = 0.4;

// Gauss-Legendre nodes for the heating integrals over one stretch of a layer. The integrand is the heating times
// smooth geometry, and a radiogenic heating follows the density, so this is far past what a layer's profile needs.
inline constexpr int d_HEATING_QUADRATURE_NODES = 16;

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
    double top_temperature  = 0.0;  // [K] at the top of its interior (an adiabat arrives colder than its base)
    double node_temperature = 0.0;  // [K] at the interface above this layer

    double heat_flow_in  = 0.0;   // [W] entering the base
    double heat_flow_out = 0.0;   // [W] leaving the top

    double boundary_thickness = 0.0;   // [m] conducting boundary layer of a convecting layer
    double resistance_bottom  = 0.0;   // [K W-1] base to the layer's own temperature
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

    // Material properties at the layer's mid-radius.
    double conductivity      = TidalPyConstants::d_NAN;   // k     [W m-1 K-1]
    double thermal_expansion = TidalPyConstants::d_NAN;   // alpha [K-1]
    double heat_capacity     = TidalPyConstants::d_NAN;   // c_p   [J kg-1 K-1]
};

// Thermal resistance [K W-1] of a conducting spherical shell. Zero for a degenerate shell or conductivity.
inline double c_shell_resistance(double radius_inner, double radius_outer, double conductivity) noexcept {
    if (!(conductivity > TidalPyConstants::d_EPS) || !(radius_outer > radius_inner)) { return 0.0; }
    const double inner = (radius_inner > TidalPyConstants::d_EPS) ? (1.0 / radius_inner) : (1.0 / radius_outer);
    return (inner - 1.0 / radius_outer) / (4.0 * TidalPyConstants::d_PI * conductivity);
}

// The temperature kind a layer's cooling model asks for. A layer that cannot hold one is isothermal.
inline c_TemperatureKind c_layer_temperature_kind(const c_BaseLayer* layer) noexcept {
    const auto* solidliquid_layer = dynamic_cast<const c_SolidLiquidLayer*>(layer);
    if (solidliquid_layer == nullptr) { return c_TemperatureKind::Isothermal; }
    const c_CoolingBase* cooling_model = solidliquid_layer->get_cooling_model();
    if (cooling_model == nullptr) { return c_TemperatureKind::Isothermal; }
    switch (cooling_model->get_model_type()) {
        case c_CoolingModel::Conduction: return c_TemperatureKind::Conductive;
        case c_CoolingModel::Convection: return c_TemperatureKind::Adiabatic;
        case c_CoolingModel::Off:        return c_TemperatureKind::Isothermal;
    }
    return c_TemperatureKind::Isothermal;
}

// Build the per-layer thermal description from the layers themselves: the temperature each carries, the kind its
// cooling model asks for, and its thermal material properties. A finite temperature_override replaces every
// layer's own temperature (a geometry-only layer then has one too). The resistances and flows are filled in later,
// against a solved structure.
inline void c_init_layer_thermal(
        const std::vector<std::unique_ptr<c_BaseLayer>>& layers,
        std::vector<c_LayerThermal>& out,
        double temperature_override = TidalPyConstants::d_NAN) {
    const std::size_t n_layers = layers.size();
    out.assign(n_layers, c_LayerThermal());
    for (std::size_t layer_i = 0; layer_i < n_layers; ++layer_i) {
        const c_BaseLayer* layer = layers[layer_i].get();
        c_LayerThermal& thermal  = out[layer_i];

        const auto* physics_layer = dynamic_cast<const c_PhysicsLayer*>(layer);
        if (std::isfinite(temperature_override)) {
            thermal.temperature = temperature_override;
        } else {
            thermal.temperature = (physics_layer != nullptr) ? physics_layer->get_temperature() : 0.0;
        }
        thermal.in_network = ((physics_layer != nullptr) || std::isfinite(temperature_override))
            && std::isfinite(thermal.temperature) && (thermal.temperature > 0.0);
        thermal.top_temperature = thermal.temperature;
        thermal.node_temperature = thermal.temperature;
        thermal.kind            = thermal.in_network ? c_layer_temperature_kind(layer) : c_TemperatureKind::Isothermal;

        // The thermal constants belong to the material, so any layer with an EOS model has them.
        if (const c_MaterialEOSBase* eos_model = layer->get_eos()) {
            thermal.conductivity      = eos_model->get_thermal_conductivity();
            thermal.thermal_expansion = eos_model->get_thermal_expansion();
            thermal.heat_capacity     = eos_model->get_heat_capacity();
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

    std::vector<double> nodes;
    std::vector<double> weights;
    c_gauss_legendre_nodes(d_HEATING_QUADRATURE_NODES, nodes, weights);
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

// Update the boundary layers, resistances, interface temperatures, and heat flows against a solved structure.
// Returns the largest relative change in the interface temperatures and flows, which is what the solve watches
// to decide it has converged.
inline double c_update_layer_thermal(
        const c_EOSSolution& solution,
        const std::vector<std::unique_ptr<c_BaseLayer>>& layers,
        double surface_temperature,
        bool solution_carries_temperature,
        std::vector<c_LayerThermal>& thermal_vec,
        const c_Heating* heating_ptr = nullptr) {
    const std::size_t n_layers = layers.size();
    double largest_change = 0.0;
    const bool heated = (heating_ptr != nullptr) && heating_ptr->get_is_active();

    for (std::size_t layer_i = 0; layer_i < n_layers; ++layer_i) {
        const c_BaseLayer* layer = layers[layer_i].get();
        c_LayerThermal& thermal  = thermal_vec[layer_i];
        const double radius_inner = layer->get_radius_inner();
        const double radius_outer = layer->get_radius_outer();
        const double thickness    = radius_outer - radius_inner;
        const double radius_mid   = 0.5 * (radius_inner + radius_outer);

        thermal.boundary_thickness = 0.0;
        thermal.resistance_bottom  = 0.0;
        thermal.resistance_top     = 0.0;
        thermal.boundary_fallback  = false;
        if (thermal.kind == c_TemperatureKind::Isothermal) { continue; }

        if (thermal.kind == c_TemperatureKind::Conductive) {
            thermal.resistance_bottom = c_shell_resistance(radius_inner, radius_mid, thermal.conductivity);
            thermal.resistance_top    = c_shell_resistance(radius_mid, radius_outer, thermal.conductivity);
            thermal.top_temperature   = thermal.temperature;
            continue;
        }

        // Convecting layer: the cooling model sizes both boundary layers from the temperature drop across the
        // layer, the local state, and the viscosity at the layer's own temperature. The drop is the sum of the two
        // boundary layers' drops: from the top of the layer below (the end of its adiabat, when it convects) to this
        // layer's temperature, and from that temperature to the layer above or the surface. The center carries no
        // flow, so the innermost layer has only the upper drop.
        double structure[C_EOS_Y_VALUES];
        solution.call_y_si(layer_i, radius_mid, structure);
        const double gravity  = structure[0];
        const double pressure = structure[1];

        const auto* solidliquid_layer = dynamic_cast<const c_SolidLiquidLayer*>(layer);
        const c_CoolingBase* cooling_model = solidliquid_layer->get_cooling_model();
        // The material at the layer's own (lumped) temperature, which is not the local temperature of the
        // profile at this radius, so it is asked of the EOS model rather than read from the solution.
        c_MaterialState material;
        if (const c_MaterialEOSBase* eos_model = layer->get_eos()) {
            eos_model->calc_material_state(
                pressure, thermal.temperature, solidliquid_layer->get_use_thermal_eos(), radius_mid, material);
        }

        // A neighbor outside the network exchanges no heat, so there is no drop across that boundary layer.
        const double inner_temperature = ((layer_i > 0) && thermal_vec[layer_i - 1].in_network)
            ? thermal_vec[layer_i - 1].top_temperature
            : thermal.temperature;
        double outer_temperature = thermal.temperature;
        if (layer_i + 1 < n_layers) {
            if (thermal_vec[layer_i + 1].in_network) { outer_temperature = thermal_vec[layer_i + 1].temperature; }
        } else if (std::isfinite(surface_temperature)) {
            outer_temperature = surface_temperature;
        }

        c_CoolingInputs cooling_inputs;
        cooling_inputs.delta_temp = std::fabs(inner_temperature - thermal.temperature)
                                  + std::fabs(thermal.temperature - outer_temperature);
        cooling_inputs.thickness  = thickness;
        cooling_inputs.gravity    = gravity;
        cooling_inputs.density    = std::isfinite(material.density) ? material.density : layer->get_density_bulk();
        cooling_inputs.viscosity  = material.shear_viscosity;
        cooling_inputs.thermal_conductivity = thermal.conductivity;
        // The diffusivity uses the density the material has here, not a separate reference density.
        cooling_inputs.thermal_diffusivity  =
            thermal.conductivity / (cooling_inputs.density * thermal.heat_capacity);
        cooling_inputs.thermal_expansion    = thermal.thermal_expansion;
        const c_CoolingResult cooling_result = cooling_model->calc_cooling(cooling_inputs);

        thermal.rayleigh_number = cooling_result.rayleigh_number;
        thermal.nusselt_number  = cooling_result.nusselt_number;
        double boundary = cooling_result.blt;
        if (!(boundary > 0.0) || !std::isfinite(boundary)) {
            // The solve reports it; the layer's viscosity is the usual cause.
            thermal.boundary_fallback = true;
            boundary = d_MAX_BOUNDARY_FRACTION * thickness;
        }
        if (boundary > d_MAX_BOUNDARY_FRACTION * thickness) { boundary = d_MAX_BOUNDARY_FRACTION * thickness; }
        thermal.boundary_thickness = boundary;

        thermal.resistance_bottom = c_shell_resistance(
            radius_inner, radius_inner + boundary, thermal.conductivity);
        thermal.resistance_top = c_shell_resistance(
            radius_outer - boundary, radius_outer, thermal.conductivity);

        // The adiabat arrives at the top of the interior colder than its base; the solved profile says by how
        // much, and the first pass (with no profile yet) starts from no drop at all.
        if (solution_carries_temperature) {
            double state[C_EOS_DY_VALUES];
            solution.call_si(layer_i, radius_outer - boundary, state);
            const double solved_top = state[C_EOS_TEMPERATURE_INDEX];
            if (std::isfinite(solved_top)) { thermal.top_temperature = solved_top; }
        }
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
        // interior of a convecting one.
        const bool convects = (thermal.kind == c_TemperatureKind::Adiabatic);
        const double bottom_end = convects
            ? (radius_inner + thermal.boundary_thickness) : 0.5 * (radius_inner + radius_outer);
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
            temperature_upper = thermal_vec[layer_i + 1].temperature + thermal_vec[layer_i + 1].heating_drop_bottom;
        }

        double node = 0.0;
        double flow = 0.0;
        const bool lower_conducts = (resistance_lower > TidalPyConstants::d_EPS);
        const bool upper_conducts = (resistance_upper > TidalPyConstants::d_EPS);
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
            // Neither side resolves a gradient: the interface carries no modeled flow.
            node = 0.5 * (temperature_lower + temperature_upper);
            flow = 0.0;
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

// Turn the per-layer thermal description into the radial segments the solve integrates. The radii and the two
// gradient coefficients are converted into the units the solve runs in.
inline void c_build_thermal_segments(
        const std::vector<c_LayerThermal>& thermal_vec,
        const std::vector<std::unique_ptr<c_BaseLayer>>& layers,
        bool integrate_temperature,
        double length_scale,
        double gravity_scale,
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
        segment.conduction_coeff = (thermal.conductivity > TidalPyConstants::d_EPS)
            ? 1.0 / (4.0 * TidalPyConstants::d_PI * thermal.conductivity * length_scale) : 0.0;
        segment.adiabat_coeff = (thermal.heat_capacity > TidalPyConstants::d_EPS)
            ? thermal.thermal_expansion * gravity_scale * length_scale / thermal.heat_capacity : 0.0;

        const bool isothermal = (!integrate_temperature) || (thermal.kind == c_TemperatureKind::Isothermal);
        if (isothermal) {
            // One segment holding this layer's own temperature. It breaks the profile at its base, which is what
            // a perfectly conducting layer does to its neighbors.
            segment.temperature_kind  = c_TemperatureKind::Isothermal;
            segment.start_temperature = thermal.temperature;
            out.push_back(segment);
            continue;
        }

        const double boundary = thermal.boundary_thickness;
        const bool has_interior = (thermal.kind == c_TemperatureKind::Adiabatic)
            && (boundary > 0.0)
            && (radius_outer - boundary > radius_inner + boundary);

        // The base of the interior (or the mid-radius of a conducting layer) is where the layer's own
        // temperature applies, so the two stretches around it carry their own heat flow.
        const double split_lower = has_interior ? (radius_inner + boundary) : 0.5 * (radius_inner + radius_outer);

        // The innermost layer has no node below to continue from, and nor does one above a layer outside the
        // network (whose profile holds its own placeholder temperature), so its base starts where the stretch has
        // to start to reach the layer's own temperature at its top.
        const bool starts_fresh = (layer_i == 0) || !thermal_vec[layer_i - 1].in_network;
        c_EOSSegment lower = segment;
        lower.temperature_kind  = c_TemperatureKind::Conductive;
        lower.start_temperature = starts_fresh
            ? (thermal.temperature + thermal.heat_flow_in * thermal.resistance_bottom + thermal.heating_drop_bottom)
            : TidalPyConstants::d_NAN;
        lower.start_heat_flow   = thermal.heat_flow_in;
        lower.upper_radius      = split_lower / length_scale;
        out.push_back(lower);

        if (has_interior) {
            c_EOSSegment interior = segment;
            interior.temperature_kind  = c_TemperatureKind::Adiabatic;
            interior.start_temperature = TidalPyConstants::d_NAN;
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
