#pragma once
/* Abstract base for TidalPy cooling (heat-transport) models, and the interface through which a model builds its
 * layer's thermal profile. Concrete models live in cooling_.hpp. All quantities MKS.
 *
 * A cooling model says how heat moves through its layer. During a thermal EOS solve the thermal network
 * (Structures/worlds/thermal_layout_.hpp) asks each layer's model for its profile (build_profile): which stretches
 * conduct and which are adiabatic, the resistance of each conducting stretch, the temperature at the base of a
 * convecting interior, and the Rayleigh and Nusselt numbers. The model reads the solved structure and its layer's
 * material through a c_LayerThermalProbe, so it chooses for itself where to evaluate its material: the convection
 * model takes its viscosity at the top of its adiabatic interior, where the layer's temperature applies. The network
 * then joins the layers' profiles into a chain of resistances.
 *
 * References
 * ----------
 * - Turcotte and Schubert (2002), Geodynamics: boundary-layer convection and the Nusselt scaling.
 * - Solomatov (1995): the stagnant-lid regime of strongly temperature-dependent viscosity.
 */

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "broadcast_.hpp"
#include "constants_.hpp"
#include "physics_base_.hpp"
#include "thermo_point_.hpp"   // c_TemperatureKind

namespace tidalpy {

// Largest share of a layer's thickness one conducting boundary layer may take, so a convecting layer keeps an
// interior to be adiabatic in.
inline constexpr double d_MAX_BOUNDARY_FRACTION = 0.4;

// The uniform conductivity [W m-1 K-1] that gives a spherical shell a resistance [K W-1]: a uniform shell's resistance
// is (1 / (4 pi k)) (1/r_inner - 1/r_outer). `fallback` for a shell from the center or one without a resistance.
inline double c_shell_conductivity(
        double radius_inner,
        double radius_outer,
        double resistance,
        double fallback) noexcept
{
    if (!(radius_inner > TidalPyConstants::d_EPS) ||
        !(radius_outer > radius_inner) ||
        !(resistance > 0.0) ||
        !std::isfinite(resistance))
    {
        return fallback;
    }
    return (1.0 / radius_inner - 1.0 / radius_outer) / (4.0 * TidalPyConstants::d_PI * resistance);
}

// Physical state passed to a cooling model's flux law.
struct c_CoolingInputs {
    double delta_temp           = 0.0;    // temperature drop across the layer, all its boundary layers [K]
    double thickness            = 0.0;    // layer (or sub-layer) thickness [m]
    double gravity              = 0.0;    // gravitational acceleration [m/s^2]
    double density              = 0.0;    // bulk density [kg/m^3]
    double viscosity            = 0.0;    // dynamic viscosity [Pa·s]
    double thermal_conductivity = 0.0;    // thermal conductivity [W/m/K]
    double thermal_diffusivity  = 0.0;    // thermal diffusivity [m^2/s]
    double thermal_expansion    = 0.0;    // thermal expansivity [1/K]
    bool   liquid               = false;  // the convecting interior is liquid (a magma ocean)
};

struct c_CoolingResult {
    double cooling_flux    = 0.0;  // heat flux leaving the layer [W/m^2]
    double blt             = 0.0;  // conducting thickness carrying the flux across the whole drop [m]
    double rayleigh_number = 0.0;  // Rayleigh number [dimensionless]
    double nusselt_number  = 1.0;  // Nusselt number [dimensionless]
};

// What a cooling model reads of its layer's material at a point.
struct c_TransportState {
    double density              = TidalPyConstants::d_NAN;   // [kg m-3]
    double thermal_conductivity = TidalPyConstants::d_NAN;   // [W m-1 K-1]
    double heat_capacity        = TidalPyConstants::d_NAN;   // effective, latent heat included [J kg-1 K-1]
    // [J kg-1 K-1] without the latent heat of a melting range: what heat diffusion sees
    double sensible_heat_capacity = TidalPyConstants::d_NAN;
    double thermal_expansion    = TidalPyConstants::d_NAN;   // [K-1]
    // [K-1] the latent heat's share of an adiabat's expansivity in a pressure-dependent melting range (c_MaterialState)
    double latent_expansion = 0.0;
    double shear_viscosity      = TidalPyConstants::d_NAN;   // post-melt [Pa s]
    double melt_fraction        = 0.0;                       // [m3 m-3]
    bool   is_liquid            = false;                     // liquid to the layer (see calc_transport_state)
};

// A layer and its surroundings during one pass of a thermal solve.
struct c_LayerThermalContext {
    double radius_inner       = 0.0;    // [m]
    double radius_outer       = 0.0;    // [m]
    double temperature        = 0.0;    // the layer's own temperature [K]
    // The temperature across the lower boundary layer from the top of the layer below, and across the upper one to
    // the base of the layer above or the surface [K]. A neighbor that exchanges no heat gives the value on this side
    // of that boundary layer (the base of the interior below, the layer's own temperature above), so no drop.
    double inner_temperature  = 0.0;
    double outer_temperature  = 0.0;
    // From the last pass: the interfaces below and above the layer [K] (each the layer's own temperature before any
    // pass), the far ends of its conducting stretches.
    double inner_node_temperature = 0.0;
    double outer_node_temperature = 0.0;
    // From the last pass: the base of this layer's interior [K] (the layer's temperature before any pass).
    double base_temperature   = 0.0;
    // No heat crosses the base (the innermost layer, or one above a layer outside the network), so no boundary layer
    // forms there.
    bool   insulated_base     = false;
};

// Read access to the solved structure and the layer's material, given to a cooling model by the thermal network.
class c_LayerThermalProbe {
public:
    virtual ~c_LayerThermalProbe() = default;

    // Gravity [m s-2] and pressure [Pa] of the solved structure at a radius of the layer [m].
    virtual void calc_structure(double radius, double& gravity, double& pressure) const = 0;

    // The layer's material at a pressure [Pa], temperature [K], and radius [m], with the layer's switches. It is liquid
    // where the Love solve takes it as liquid: everywhere in a liquid layer (a liquid-only material, or a forced
    // liquid state), and, in a layer that can change state, where it is fully molten or its post-melt shear modulus is
    // at or below the radial solver's solid threshold (minimum_solid_rigidity times the planet's rho g R).
    virtual void calc_transport_state(
            double pressure,
            double temperature,
            double radius,
            c_TransportState& out) const = 0;

    // The conduction resistance [K W-1] of the layer between two radii [m] whose ends are at two temperatures [K], as
    // steady conduction through the material gives it: (1/r_inner - 1/r_outer) / (4 pi k_mean), with k_mean the
    // material's conductivity averaged over the temperature between the ends (the Kirchhoff transform, exact for a
    // conductivity that follows the temperature). Zero for a shell from the center (no heat crosses it there), an
    // empty shell, or a material that does not conduct.
    virtual double calc_shell_resistance(
            double radius_inner,
            double radius_outer,
            double temperature_inner,
            double temperature_outer) const = 0;

    // The temperature [K] at the base (radius_lower [m]) of an adiabatic stretch whose top (radius_upper [m]) is at
    // top_temperature [K]: dT/dr = -(alpha + alpha_L) g T / c_p over the solved structure, with alpha_L the material's
    // latent expansion. top_temperature for an empty stretch.
    virtual double calc_adiabat_base_temperature(
            double radius_lower,
            double radius_upper,
            double top_temperature) const = 0;
};

// What a cooling model makes of its layer for one pass of the thermal network.
// The model's kind (get_temperature_kind) says how to read it. Isothermal: one stretch at the layer's temperature.
// Conductive: two conducting halves meeting at the mid-radius, where the layer's temperature applies. Adiabatic:
// conducting boundary layers around an adiabatic interior whose top is at the layer's temperature.
struct c_LayerThermalProfile {
    double boundary_thickness = 0.0;   // [m] each conducting boundary layer of a convecting layer
    double resistance_bottom  = 0.0;   // [K W-1] base to the layer's own temperature (zero with no base stretch)
    double resistance_top     = 0.0;   // [K W-1] the layer's own temperature to its top
    double base_temperature   = 0.0;   // [K] the base of the interior (a convecting interior is warmer at its base)

    double rayleigh_number = 0.0;
    double nusselt_number  = 1.0;
    // True when the flux law gave no usable boundary-layer thickness (a NaN viscosity, say), so each boundary layer
    // took the largest share of the layer it may.
    bool   boundary_fallback = false;
    // True when the convecting interior was liquid at its reference point and took the liquid scaling.
    bool   magma_ocean = false;

    // The material the profile used for its conduction and its adiabat: at the mid-radius pressure and the layer's own
    // temperature. NaN where the model read none (an isothermal layer).
    double conductivity      = TidalPyConstants::d_NAN;   // [W m-1 K-1]
    double heat_capacity     = TidalPyConstants::d_NAN;   // [J kg-1 K-1], latent heat included
    double thermal_expansion = TidalPyConstants::d_NAN;   // [K-1]

    // Where a convecting layer's viscosity was evaluated: the top of its interior, at the layer's own temperature.
    double reference_pressure      = TidalPyConstants::d_NAN;   // [Pa]
    double reference_viscosity     = TidalPyConstants::d_NAN;   // [Pa s]
    double reference_melt_fraction = TidalPyConstants::d_NAN;   // [m3 m-3]
};

class c_CoolingBase : public c_PhysicsBase {
public:
    static constexpr const char* C_FAMILY_NAME = "cooling";

    c_CoolingBase() = default;

    explicit c_CoolingBase(const std::string& model_name) : c_PhysicsBase(model_name) {}

    ~c_CoolingBase() override = default;

    // The kind of profile the model gives its layer, known before any structure is solved.
    virtual c_TemperatureKind get_temperature_kind() const noexcept = 0;

    // The model's flux law for one state. Assumes steady-state boundary-layer theory.
    virtual c_CoolingResult calc_cooling(const c_CoolingInputs& inputs) const = 0;

    // The layer's thermal profile for one pass of the network (see the header comment).
    virtual void build_profile(
            const c_LayerThermalContext& context,
            const c_LayerThermalProbe& probe,
            c_LayerThermalProfile& out) const = 0;

    // Element-wise over the temperature drop and the viscosity at an otherwise fixed state (base_inputs). Each holds
    // one value per point or a single value used at every point. The four outputs are caller-owned buffers of
    // num_points values, the broadcast length of the two inputs. Each model implements it with p_vectorize_kernel, so
    // its kernel inlines into the loop rather than costing a virtual call per point.
    virtual void calc_cooling_vectorize(
            const std::vector<double>& delta_temp,
            const std::vector<double>& viscosity,
            const c_CoolingInputs& base_inputs,
            std::size_t num_points,
            double* out_cooling_flux,
            double* out_blt,
            double* out_rayleigh,
            double* out_nusselt) const = 0;

protected:
    // The shared loop over a model's kernel, a callable taking c_CoolingInputs and returning a c_CoolingResult. A model
    // passes a lambda holding a copy of its parameters: the compiler inlines it and, unlike the model's own members,
    // the copy cannot alias the output buffers, so the parameters stay in registers across the loop.
    template <class Kernel>
    static void p_vectorize_kernel(
            const Kernel& kernel,
            const std::vector<double>& delta_temp,
            const std::vector<double>& viscosity,
            const c_CoolingInputs& base_inputs,
            std::size_t num_points,
            double* out_cooling_flux,
            double* out_blt,
            double* out_rayleigh,
            double* out_nusselt)
    {
        if (c_broadcast_length({delta_temp.size(), viscosity.size()}, "calc_cooling_vectorize") != num_points)
        {
            throw std::invalid_argument("TidalPy::calc_cooling_vectorize: the output buffers do not match the inputs");
        }
        if (num_points == 0)
        {
            return;
        }
        // Which inputs vary is fixed at compile time: a runtime broadcast stride in the index costs a sweep about an
        // eighth of its time.
        const bool delta_temp_varies = c_broadcast_stride(delta_temp.size()) == 1;
        const bool viscosity_varies  = c_broadcast_stride(viscosity.size()) == 1;
        if (delta_temp_varies && viscosity_varies)
        {
            p_sweep<true, true>(kernel, delta_temp, viscosity, base_inputs, num_points,
                                out_cooling_flux, out_blt, out_rayleigh, out_nusselt);
        }
        else if (delta_temp_varies)
        {
            p_sweep<true, false>(kernel, delta_temp, viscosity, base_inputs, num_points,
                                 out_cooling_flux, out_blt, out_rayleigh, out_nusselt);
        }
        else if (viscosity_varies)
        {
            p_sweep<false, true>(kernel, delta_temp, viscosity, base_inputs, num_points,
                                 out_cooling_flux, out_blt, out_rayleigh, out_nusselt);
        }
        else
        {
            p_sweep<false, false>(kernel, delta_temp, viscosity, base_inputs, num_points,
                                  out_cooling_flux, out_blt, out_rayleigh, out_nusselt);
        }
    }

    // The material properties a model's profile reports, read at the layer's mid-radius pressure and own temperature;
    // a property the material does not give is left NaN. Returns the transport state for the caller's further use.
    static c_TransportState p_mid_layer_state(
            const c_LayerThermalContext& context,
            const c_LayerThermalProbe& probe,
            double& gravity_out,
            double& pressure_out,
            c_LayerThermalProfile& out) {
        const double radius_mid = 0.5 * (context.radius_inner + context.radius_outer);
        probe.calc_structure(radius_mid, gravity_out, pressure_out);
        c_TransportState state;
        probe.calc_transport_state(pressure_out, context.temperature, radius_mid, state);
        out.conductivity      = state.thermal_conductivity;
        out.heat_capacity     = state.heat_capacity;
        out.thermal_expansion = state.thermal_expansion;
        return state;
    }

private:
    template <bool DELTA_TEMP_VARIES, bool VISCOSITY_VARIES, class Kernel>
    static void p_sweep(
            const Kernel& kernel,
            const std::vector<double>& delta_temp,
            const std::vector<double>& viscosity,
            const c_CoolingInputs& base_inputs,
            std::size_t num_points,
            double* out_cooling_flux,
            double* out_blt,
            double* out_rayleigh,
            double* out_nusselt)
    {
        const double* const delta_temp_data = delta_temp.data();
        const double* const viscosity_data  = viscosity.data();
        c_CoolingInputs inputs = base_inputs;
        inputs.delta_temp = delta_temp_data[0];
        inputs.viscosity  = viscosity_data[0];

        for (std::size_t i = 0; i < num_points; ++i)
        {
            if constexpr (DELTA_TEMP_VARIES) { inputs.delta_temp = delta_temp_data[i]; }
            if constexpr (VISCOSITY_VARIES)  { inputs.viscosity  = viscosity_data[i]; }
            const c_CoolingResult result = kernel(inputs);
            out_cooling_flux[i] = result.cooling_flux;
            out_blt[i]          = result.blt;
            out_rayleigh[i]     = result.rayleigh_number;
            out_nusselt[i]      = result.nusselt_number;
        }
    }
};

} // namespace tidalpy
