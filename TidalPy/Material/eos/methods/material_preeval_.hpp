#pragma once
/* The EOS function of a layer during the whole-planet structure solve: it asks the layer's material (c_Material) for
 * its properties at the current radius, pressure, and temperature, with the layer's switches. While the solve iterates
 * it needs the density alone (and the thermal properties when it integrates temperature); the full state follows only
 * when the input's full_state flag is set, which is every dense evaluation of the finished solution. The material's
 * rigidity margin is the event that splits a layer that can change state into solid and liquid zones.
 */

#include <complex>

#include "../ode_.hpp"
#include "../../material_.hpp"
#include "../../../constants_.hpp"


/// The material works in SI, so a non-dimensional solve scales the radius and pressure up before the call and the
/// density and moduli down afterwards; the viscosities and thermal properties come back in SI either way.
struct c_MaterialPreevalInput
{
    const tidalpy::c_Material* material_ptr = nullptr;
    tidalpy::c_MaterialSwitches switches;
    // The layer's temperature [K]; a solve that integrates temperature takes the local state value instead.
    double temperature           = 0.0;
    bool   use_state_temperature = false;
    double length_scale  = 1.0;
    double pascal_scale  = 1.0;
    double density_scale = 1.0;
};


/// The material point of a structure-ODE state, in SI: the state temperature in a solve that integrates it.
inline tidalpy::c_ThermoPoint c_preeval_point(
        const c_MaterialPreevalInput& eos_data,
        double radius,
        const double* radial_solutions) noexcept
{
    tidalpy::c_ThermoPoint point;
    point.radius      = radius * eos_data.length_scale;
    point.pressure    = radial_solutions[1] * eos_data.pascal_scale;
    point.temperature = eos_data.use_state_temperature ? radial_solutions[4] : eos_data.temperature;
    return point;
}


/// CyRK EventFunc: the rigidity margin ln(mu / mu_min) of the layer's material at the state, with mu its post-melt
/// shear modulus and mu_min the input's minimum_shear_modulus. Positive while the material is solid enough for the
/// radial solver's solid equations, negative (-inf with no shear modulus) where it is a liquid to them. A threshold
/// of zero makes the material a liquid only where its shear modulus vanishes (c_Material::calc_rigidity_margin).
inline double c_preeval_rigidity_margin(double radius, double* radial_solutions, char* input_args) noexcept
{
    const c_EOS_ODEInput* ode_args = reinterpret_cast<const c_EOS_ODEInput*>(input_args);
    const c_MaterialPreevalInput* eos_data = reinterpret_cast<const c_MaterialPreevalInput*>(ode_args->eos_input_ptr);
    const tidalpy::c_ThermoPoint point = c_preeval_point(*eos_data, radius, radial_solutions);
    return eos_data->material_ptr->calc_rigidity_margin(point, eos_data->switches, ode_args->minimum_shear_modulus);
}


/// CyRK PreEvalFunc signature.
inline void c_preeval_material(
        char* preeval_output,
        double radius,
        double* radial_solutions,
        char* preeval_input
        ) noexcept
{
    c_EOS_ODEInput*         ode_args = reinterpret_cast<c_EOS_ODEInput*>(preeval_input);
    c_MaterialPreevalInput* eos_data = reinterpret_cast<c_MaterialPreevalInput*>(ode_args->eos_input_ptr);
    c_EOSOutput* output = reinterpret_cast<c_EOSOutput*>(preeval_output);
    const tidalpy::c_Material& material = *eos_data->material_ptr;

    // The state temperature sits at index 4 of a thermal solve; the flag is only set for one.
    const tidalpy::c_ThermoPoint point = c_preeval_point(*eos_data, radius, radial_solutions);

    output->shear_modulus   = std::complex<double>(TidalPyConstants::d_NAN, 0.0);
    output->bulk_modulus    = std::complex<double>(TidalPyConstants::d_NAN, 0.0);
    output->shear_viscosity = TidalPyConstants::d_NAN;
    output->bulk_viscosity  = TidalPyConstants::d_NAN;
    output->melt_fraction   = TidalPyConstants::d_NAN;

    if (!(ode_args->full_state || ode_args->thermal_state))
    {
        output->density = material.calc_density(point, eos_data->switches) / eos_data->density_scale;
        return;
    }

    tidalpy::c_MaterialState state;
    if (ode_args->full_state) { material.calc_state(point, eos_data->switches, state); }
    else                      { material.calc_thermal(point, eos_data->switches, state); }
    output->density            = state.density / eos_data->density_scale;
    output->melt_fraction      = state.melt_fraction;
    output->thermal_expansion  = state.thermal_expansion;
    output->heat_capacity      = state.heat_capacity;
    output->latent_expansion   = state.latent_expansion;
    output->thermal_conductivity = state.thermal_conductivity;
    if (!ode_args->full_state) { return; }

    // A tidal deformation is adiabatic, so the radial solver and the getters see the adiabatic bulk modulus.
    output->shear_modulus   = std::complex<double>(state.shear_modulus / eos_data->pascal_scale, 0.0);
    output->bulk_modulus    = std::complex<double>(state.adiabatic_bulk_modulus / eos_data->pascal_scale, 0.0);
    output->shear_viscosity = state.shear_viscosity;
    output->bulk_viscosity  = state.bulk_viscosity;

    // N^2 = -g (rho' + rho^2 g / K_S) / rho, with rho' = d rho / dr by the chain rule through the material's own
    // density along the structure ODE's dP/dr = -rho g and dT/dr (zero where the temperature is not integrated, which
    // is uniform in each segment). rho^2 g / K_S = -(rho / K_S) dP/dr joins the pressure slope before the sum, so a
    // material whose density follows its bulk modulus (d rho / dP = rho / K_S, no other slope) gives exactly 0 rather
    // than the roundoff of two large terms, which would set a dynamic liquid's response at very long periods.
    const double pressure_gradient_si = -output->density * radial_solutions[0] * eos_data->pascal_scale;
    const double temperature_gradient = (eos_data->use_state_temperature && (radius > TidalPyConstants::d_EPS_10))
        ? c_eos_temperature_gradient(radius, radial_solutions, *output, *ode_args) : 0.0;
    const double adiabatic_pressure_slope =
        (state.adiabatic_bulk_modulus > 0.0) ? state.density / state.adiabatic_bulk_modulus : 0.0;
    const double gradient_off_pressure = state.density_temperature_slope * temperature_gradient
        + state.density_radius_slope * eos_data->length_scale;
    const double gradient_off_adiabat =
        ((state.density_pressure_slope - adiabatic_pressure_slope) * pressure_gradient_si + gradient_off_pressure)
        / eos_data->density_scale;
    output->buoyancy_frequency_squared =
        (output->density > 0.0) ? -radial_solutions[0] * gradient_off_adiabat / output->density : 0.0;
    output->density_gradient =
        (state.density_pressure_slope * pressure_gradient_si + gradient_off_pressure) / eos_data->density_scale;
}
