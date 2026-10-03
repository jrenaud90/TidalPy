#pragma once
/* TidalPy's cooling (heat-transport) models. All quantities MKS.
 *
 * Each model declares its parameters in one table (c_SpecModel, spec_model_.hpp), which gives its constructor,
 * validation, config entries, and binary record, and builds its layer's thermal profile (build_profile, see
 * cooling_base_.hpp).
 *
 * References
 * ----------
 * - Turcotte and Schubert (2002), Geodynamics: Rayleigh and Nusselt convection scaling.
 * - Solomatov (1995); Schubert, Turcotte, and Olson (2001): boundary-layer theory.
 * - Stevenson, Spohn, and Schubert (1983): the upper-mantle temperature of parameterized convection, at the top of
 *   the adiabatic interior.
 * - Solomatov (2000), Fluid dynamics of a terrestrial magma ocean; Lebrun et al. (2013): Nu = 0.089 Ra^(1/3) for a
 *   liquid (soft-turbulence) interior.
 */

#include <algorithm>
#include <cmath>
#include <istream>
#include <memory>
#include <string>
#include <vector>

#include "binary_.hpp"
#include "constants_.hpp"
#include "cooling_base_.hpp"
#include "registry_.hpp"
#include "spec_model_.hpp"

namespace tidalpy {

// No heat transport inside the layer: zero flux, and no boundary layer, so its thickness is NaN (not applicable).
TIDALPY_FORCE_INLINE c_CoolingResult cool_off(const c_CoolingInputs& /*in*/) noexcept {
    c_CoolingResult result;
    result.cooling_flux    = 0.0;
    result.blt             = TidalPyConstants::d_NAN;
    result.rayleigh_number = 0.0;
    result.nusselt_number  = 1.0;
    return result;
}

// Conduction across the whole layer: flux = k * delta_temp / thickness.
TIDALPY_FORCE_INLINE c_CoolingResult cool_conduction(const c_CoolingInputs& in) noexcept {
    c_CoolingResult result;
    result.blt             = in.thickness;
    result.cooling_flux    = in.thermal_conductivity * in.delta_temp / c_guard_denominator(in.thickness);
    result.rayleigh_number = 0.0;
    result.nusselt_number  = 1.0;
    return result;
}

// Parameterized convection via the Rayleigh number.
//
//   Ra = expansion * density * gravity * delta_temp * thickness^3 / (viscosity * diffusivity)
//   Nu = max(alpha * (Ra / Ra_crit)^beta, Nu_min)
//   boundary layer = thickness / Nu
//   flux = k * delta_temp / boundary_layer = Nu * k * delta_temp / thickness
//
// delta_temp is the whole drop across the layer (the sum of its boundary layers' drops), and Nu is the flux over
// that of conduction across the whole layer (Turcotte and Schubert 2002), so Nu = 1 is conduction. The boundary
// layer is the conducting thickness that carries the flux across the whole drop: a layer with a boundary layer at
// its base and its top splits the drop between them at that flux, so each is half this thick. Nu_min is the
// [numerical] minimum_nusselt setting; at its default of 1 a sub-critical or rigid layer conducts. Degenerate
// inputs (no temperature contrast, or a vanishingly thin layer) collapse to Ra = 0 and Nu = Nu_min. Each test is
// made once so that Ra, Nu, and the boundary layer stay consistent with one another. A liquid interior (a magma
// ocean) takes the liquid scaling, Nu = alpha_liquid Ra^beta_liquid (Ra_crit = 1).
TIDALPY_FORCE_INLINE c_CoolingResult cool_convection(
        const c_CoolingInputs& in,
        double alpha,
        double beta,
        double critical_rayleigh) noexcept {
    const double eps = TidalPyConstants::d_EPS;
    // An unwired config gives NaN limits, so the result shows the missing initialization.
    const bool config_wired    = (tidalpy_config_ptr != nullptr);
    const double min_thickness = config_wired ? tidalpy_config_ptr->d_MIN_THICKNESS : TidalPyConstants::d_NAN;
    const double min_nusselt   = config_wired ? tidalpy_config_ptr->d_MIN_NUSSELT : TidalPyConstants::d_NAN;
    c_CoolingResult result;

    const double rate_heat_loss   = in.thermal_diffusivity / c_guard_denominator(in.thickness);
    const double parcel_rise_rate = in.thermal_expansion * in.density * in.gravity
                                  * in.delta_temp * in.thickness * in.thickness
                                  / c_guard_denominator(in.viscosity);

    const bool no_contrast = !(in.delta_temp > eps);
    const bool too_thin    = !(in.thickness > min_thickness);

    double rayleigh = parcel_rise_rate / c_guard_denominator(rate_heat_loss);
    if (no_contrast || too_thin) { rayleigh = 0.0; }

    double nusselt = alpha * std::pow(rayleigh / c_guard_denominator(critical_rayleigh), beta);
    // A NaN from another input (the viscosity, say) falls through so it reaches the caller.
    if (no_contrast || too_thin || (nusselt <= min_nusselt)) { nusselt = min_nusselt; }

    // With no contrast the boundary layer still follows the floored Nusselt number: the thermal network uses it
    // as the resistance between this layer and its neighbors, whose temperatures can differ from this one's.
    double blt = in.thickness / c_guard_denominator(nusselt);
    if (too_thin) { blt = in.thickness; }

    result.cooling_flux    = in.thermal_conductivity * in.delta_temp / c_guard_denominator(blt);
    result.blt             = blt;
    result.rayleigh_number = rayleigh;
    result.nusselt_number  = nusselt;
    return result;
}

// No heat transport inside the layer (alias "none"): the layer holds one temperature.
class c_OffCooling final : public c_SpecModel<c_OffCooling, c_CoolingBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::OffCooling;

    static const std::vector<c_ParamSpec<c_OffCooling>>& parameter_specs() {
        static const std::vector<c_ParamSpec<c_OffCooling>> specs = {};
        return specs;
    }

    c_OffCooling() : c_OffCooling(c_ParamMap{}) {}
    explicit c_OffCooling(const c_ParamMap& params) : c_SpecModel("off") { this->p_initialize(params); }

    c_TemperatureKind get_temperature_kind() const noexcept override { return c_TemperatureKind::Isothermal; }

    c_CoolingResult calc_cooling(const c_CoolingInputs& inputs) const override { return cool_off(inputs); }

    void build_profile(
            const c_LayerThermalContext& context,
            const c_LayerThermalProbe& /*probe*/,
            c_LayerThermalProfile& out) const override {
        out = c_LayerThermalProfile();
        out.base_temperature = context.temperature;
    }

    void calc_cooling_vectorize(
            const std::vector<double>& delta_temp,
            const std::vector<double>& viscosity,
            const c_CoolingInputs& base_inputs,
            std::size_t num_points,
            double* out_cooling_flux,
            double* out_blt,
            double* out_rayleigh,
            double* out_nusselt) const override {
        p_vectorize_kernel(
            [](const c_CoolingInputs& inputs) { return cool_off(inputs); },
            delta_temp, viscosity, base_inputs, num_points, out_cooling_flux, out_blt, out_rayleigh, out_nusselt);
    }
};

// Conduction through the whole layer (alias "conductive"): two conducting halves meeting at the mid-radius, where
// the layer's temperature applies.
class c_ConductiveCooling final : public c_SpecModel<c_ConductiveCooling, c_CoolingBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::ConductiveCooling;

    static const std::vector<c_ParamSpec<c_ConductiveCooling>>& parameter_specs() {
        static const std::vector<c_ParamSpec<c_ConductiveCooling>> specs = {};
        return specs;
    }

    c_ConductiveCooling() : c_ConductiveCooling(c_ParamMap{}) {}
    explicit c_ConductiveCooling(const c_ParamMap& params) : c_SpecModel("conduction") { this->p_initialize(params); }

    c_TemperatureKind get_temperature_kind() const noexcept override { return c_TemperatureKind::Conductive; }

    c_CoolingResult calc_cooling(const c_CoolingInputs& inputs) const override { return cool_conduction(inputs); }

    void build_profile(
            const c_LayerThermalContext& context,
            const c_LayerThermalProbe& probe,
            c_LayerThermalProfile& out) const override {
        out = c_LayerThermalProfile();
        double gravity  = 0.0;
        double pressure = 0.0;
        p_mid_layer_state(context, probe, gravity, pressure, out);
        const double radius_mid = 0.5 * (context.radius_inner + context.radius_outer);
        out.resistance_bottom = probe.calc_shell_resistance(
            context.radius_inner, radius_mid, context.inner_node_temperature, context.temperature);
        out.resistance_top = probe.calc_shell_resistance(
            radius_mid, context.radius_outer, context.temperature, context.outer_node_temperature);
        out.base_temperature  = context.temperature;
    }

    void calc_cooling_vectorize(
            const std::vector<double>& delta_temp,
            const std::vector<double>& viscosity,
            const c_CoolingInputs& base_inputs,
            std::size_t num_points,
            double* out_cooling_flux,
            double* out_blt,
            double* out_rayleigh,
            double* out_nusselt) const override {
        p_vectorize_kernel(
            [](const c_CoolingInputs& inputs) { return cool_conduction(inputs); },
            delta_temp, viscosity, base_inputs, num_points, out_cooling_flux, out_blt, out_rayleigh, out_nusselt);
    }
};

// The upper boundary layer of a convecting layer and the reference point its viscosity is taken at (just under that
// layer) depend on each other. On a solved structure each is a cheap function of the other, so the model iterates the
// pair to self-consistency within one profile: at most this many times, until the boundary thickness changes by less
// than this relative amount. The iteration starts from the last pass's boundary layer (the layer's top before there is
// one), so it settles on the same boundary layer the passes would, without a structure integration per step; near the
// melt onset, where the viscosity follows the pressure through the melt fraction, the passes alone took tens.
inline constexpr int d_REFERENCE_ITERATIONS = 200;
inline constexpr double d_REFERENCE_TOLERANCE = 1.0e-12;

// Parameterized boundary-layer convection (alias "convective"): a conducting boundary layer at the base and the top,
// from the Nusselt scaling, around an adiabatic interior whose top is at the layer's temperature (the upper-mantle
// temperature of parameterized convection). The viscosity of the Rayleigh number is the material's at that top,
// the layer's temperature and the pressure there; a material liquid there, as the Love solve takes it
// (c_LayerThermalProbe::calc_transport_state), is a magma ocean and takes the liquid scaling. Gravity, density, and
// the thermal constants are the mid-layer values, as they are for a conducting layer.
class c_ConvectiveCooling final : public c_SpecModel<c_ConvectiveCooling, c_CoolingBase> {
public:
    static constexpr BinaryClassID C_CLASS_ID = BinaryClassID::ConvectiveCooling;

    static const std::vector<c_ParamSpec<c_ConvectiveCooling>>& parameter_specs() {
        using Self = c_ConvectiveCooling;
        static const std::vector<c_ParamSpec<Self>> specs = {
            {"convection_alpha", "convection_alpha", &Self::p_convection_alpha, 1.0, c_ParamBounds::Positive,
             "Nusselt prefactor of a solid interior: Nu = alpha (Ra / Ra_crit)^beta."},
            {"convection_beta", "convection_beta", &Self::p_convection_beta, 1.0 / 3.0, c_ParamBounds::Positive,
             "Nusselt exponent of a solid interior."},
            {"critical_rayleigh", "critical_rayleigh", &Self::p_critical_rayleigh, 1100.0, c_ParamBounds::Positive,
             "Critical Rayleigh number of a solid interior."},
            {"liquid_convection_alpha", "liquid_convection_alpha", &Self::p_liquid_convection_alpha, 0.089,
             c_ParamBounds::Positive,
             "Nusselt prefactor of a liquid interior (a magma ocean): Nu = alpha Ra^beta (Solomatov 2000)."},
            {"liquid_convection_beta", "liquid_convection_beta", &Self::p_liquid_convection_beta, 1.0 / 3.0,
             c_ParamBounds::Positive, "Nusselt exponent of a liquid interior."},
        };
        return specs;
    }

    c_ConvectiveCooling() : c_ConvectiveCooling(c_ParamMap{}) {}
    explicit c_ConvectiveCooling(const c_ParamMap& params) : c_SpecModel("convection") {
        this->p_initialize(params);
    }

    c_TemperatureKind get_temperature_kind() const noexcept override { return c_TemperatureKind::Adiabatic; }

    c_CoolingResult calc_cooling(const c_CoolingInputs& inputs) const override {
        return inputs.liquid
            ? cool_convection(inputs, this->p_liquid_convection_alpha, this->p_liquid_convection_beta, 1.0)
            : cool_convection(inputs, this->p_convection_alpha, this->p_convection_beta, this->p_critical_rayleigh);
    }

    void build_profile(
            const c_LayerThermalContext& context,
            const c_LayerThermalProbe& probe,
            c_LayerThermalProfile& out) const override {
        out = c_LayerThermalProfile();
        double gravity      = 0.0;
        double pressure_mid = 0.0;
        const c_TransportState mid = p_mid_layer_state(context, probe, gravity, pressure_mid, out);
        const double thickness = context.radius_outer - context.radius_inner;

        // The drop is the sum of the two boundary layers' drops: from the top of the layer below to the base of this
        // layer's adiabat, and from this layer's temperature to the base of the layer above or the surface. The
        // adiabat's base comes from the last pass, so it lags by one. Only an unstable drop (hotter below) drives
        // convection, so a boundary layer the wrong way round (a core colder than the base of the mantle above it, a
        // layer heated from above) adds none rather than counting as if it were heated from below.
        c_CoolingInputs inputs;
        inputs.delta_temp = std::max(context.inner_temperature - context.base_temperature, 0.0)
                          + std::max(context.temperature - context.outer_temperature, 0.0);
        inputs.thickness            = thickness;
        inputs.gravity              = gravity;
        inputs.density              = mid.density;
        inputs.thermal_conductivity = out.conductivity;
        // The diffusivity uses the density the material has here, not a separate reference density, and the sensible
        // heat capacity: a melting range's latent heat buffers the temperature but does not slow heat diffusion.
        inputs.thermal_diffusivity  = out.conductivity / (mid.density * mid.sensible_heat_capacity);
        // The expansivity at the layer's own state, which can fall with compression.
        inputs.thermal_expansion    = out.thermal_expansion;

        // The reference point is the top of the interior, under the upper boundary layer, where the layer's own
        // temperature applies; the boundary layer follows from the viscosity there (d_REFERENCE_ITERATIONS). The flux
        // law's thickness carries its flux across the whole drop; two boundary layers split the drop at that flux, so
        // each takes half of it. A base that carries no heat has no boundary layer.
        double boundary = context.boundary_thickness;
        for (int iteration = 0; iteration < d_REFERENCE_ITERATIONS; ++iteration) {
            const double trial = boundary;
            double reference_radius = context.radius_outer - trial;
            if (reference_radius < context.radius_inner) { reference_radius = context.radius_inner; }
            double reference_gravity  = 0.0;
            double reference_pressure = 0.0;
            probe.calc_structure(reference_radius, reference_gravity, reference_pressure);
            c_TransportState reference;
            probe.calc_transport_state(reference_pressure, context.temperature, reference_radius, reference);
            out.reference_pressure      = reference_pressure;
            out.reference_viscosity     = reference.shear_viscosity;
            out.reference_melt_fraction = reference.melt_fraction;
            out.magma_ocean             = reference.is_liquid;

            inputs.viscosity = reference.shear_viscosity;
            inputs.liquid    = reference.is_liquid;
            const c_CoolingResult result = this->calc_cooling(inputs);
            out.rayleigh_number = result.rayleigh_number;
            out.nusselt_number  = result.nusselt_number;

            boundary = context.insulated_base ? result.blt : 0.5 * result.blt;
            out.boundary_fallback = !(boundary > 0.0) || !std::isfinite(boundary);
            // The solve reports a fallback; the layer's viscosity is the usual cause.
            if (out.boundary_fallback) { boundary = d_MAX_BOUNDARY_FRACTION * thickness; }
            if (boundary > d_MAX_BOUNDARY_FRACTION * thickness) { boundary = d_MAX_BOUNDARY_FRACTION * thickness; }
            if (std::fabs(boundary - trial) <= d_REFERENCE_TOLERANCE * boundary) { break; }
        }
        out.boundary_thickness = boundary;

        // The layer's temperature holds at the top of the interior, and the adiabat warms downward from it to the
        // interior's base, along the gravity of the solved structure.
        const double interior_base = context.radius_inner + (context.insulated_base ? 0.0 : boundary);
        out.base_temperature = probe.calc_adiabat_base_temperature(
            interior_base, context.radius_outer - boundary, context.temperature);
        out.resistance_bottom = context.insulated_base ? 0.0 : probe.calc_shell_resistance(
            context.radius_inner, interior_base, context.inner_node_temperature, out.base_temperature);
        out.resistance_top = probe.calc_shell_resistance(
            context.radius_outer - boundary, context.radius_outer, context.temperature, context.outer_node_temperature);
    }

    void calc_cooling_vectorize(
            const std::vector<double>& delta_temp,
            const std::vector<double>& viscosity,
            const c_CoolingInputs& base_inputs,
            std::size_t num_points,
            double* out_cooling_flux,
            double* out_blt,
            double* out_rayleigh,
            double* out_nusselt) const override {
        const bool liquid = base_inputs.liquid;
        const double alpha = liquid ? this->p_liquid_convection_alpha : this->p_convection_alpha;
        const double beta  = liquid ? this->p_liquid_convection_beta : this->p_convection_beta;
        const double critical_rayleigh = liquid ? 1.0 : this->p_critical_rayleigh;
        p_vectorize_kernel(
            [alpha, beta, critical_rayleigh](const c_CoolingInputs& inputs) {
                return cool_convection(inputs, alpha, beta, critical_rayleigh);
            },
            delta_temp, viscosity, base_inputs, num_points, out_cooling_flux, out_blt, out_rayleigh, out_nusselt);
    }

protected:
    double p_convection_alpha        = 0.0;
    double p_convection_beta         = 0.0;
    double p_critical_rayleigh       = 0.0;
    double p_liquid_convection_alpha = 0.0;
    double p_liquid_convection_beta  = 0.0;
};

inline const c_ModelRegistry<c_CoolingBase>& c_cooling_registry() {
    static const c_ModelRegistry<c_CoolingBase> registry = {
        {{"off", "none"},
         BinaryClassID::OffCooling,        &c_make_entry<c_CoolingBase, c_OffCooling>},
        {{"convection", "convective"},
         BinaryClassID::ConvectiveCooling, &c_make_entry<c_CoolingBase, c_ConvectiveCooling>},
        {{"conduction", "conductive"},
         BinaryClassID::ConductiveCooling, &c_make_entry<c_CoolingBase, c_ConductiveCooling>},
    };
    return registry;
}

// The family's entry points, each one line over the generic registry functions.
inline std::unique_ptr<c_CoolingBase> c_find_cooling(const std::string& model_name, const c_ParamMap& params) {
    return c_make_model(c_cooling_registry(), model_name, params);
}

inline std::unique_ptr<c_CoolingBase> c_cooling_from_binary(std::istream& in, bool force = false) {
    return c_model_from_binary(c_cooling_registry(), in, force);
}

inline std::string c_cooling_canonical_name(const std::string& model_name) {
    return c_canonical_model_name(c_cooling_registry(), model_name);
}

inline std::vector<std::string> c_cooling_model_names() {
    return c_model_names(c_cooling_registry());
}

}  // namespace tidalpy
