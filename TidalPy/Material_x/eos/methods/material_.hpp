#pragma once
/*
 * material_.hpp: EOS pre-evaluation that asks a layer's material EOS model (c_MaterialEOSBase) for the density at
 * the current radius and pressure (y[1]) during the whole-planet structure solve. The static moduli and
 * viscosities are requested only when the solve's update_shear / update_bulk flags are set (the final pass).
 */

#include <complex>

#include "../ode_.hpp"             // c_EOS_ODEInput, c_EOSOutput, PreEvalFunc signature
#include "../material_eos_.hpp"    // tidalpy::c_MaterialEOSBase
#include "../../../constants_.hpp" // TidalPyConstants::d_NAN


/// Input for the material-EOS pre-evaluation: a non-owning model pointer, the (unused) temperature, and the unit
/// scales of the solve. The models work in SI, so a non-dimensional solve scales the radius and pressure up before
/// the call and the density and moduli down afterwards; the viscosities are returned in SI either way.
struct c_MaterialEOSInput
{
    tidalpy::c_MaterialEOSBase* eos_model_ptr = nullptr;
    // Temperature [K] handed to the EOS model, or NaN for its athermal density. A thermal solve overrides it
    // with the local state temperature when the layer asked for a thermal EOS.
    double temperature   = 0.0;
    bool   use_state_temperature = false;
    double length_scale  = 1.0;
    double pascal_scale  = 1.0;
    double density_scale = 1.0;
};


/// Material-EOS pre-evaluation (CyRK PreEvalFunc signature).
inline void c_preeval_material_eos(
        char* preeval_output,
        double radius,
        double* radial_solutions,
        char* preeval_input
        ) noexcept
{
    c_EOS_ODEInput*     ode_args = reinterpret_cast<c_EOS_ODEInput*>(preeval_input);
    c_MaterialEOSInput* eos_data = reinterpret_cast<c_MaterialEOSInput*>(ode_args->eos_input_ptr);
    c_EOSOutput* output = reinterpret_cast<c_EOSOutput*>(preeval_output);

    const double radius_si   = radius * eos_data->length_scale;
    const double pressure_si = radial_solutions[1] * eos_data->pascal_scale;
    // The state temperature sits at index 4 of a thermal solve; the flag is only set for one.
    const double temperature = eos_data->use_state_temperature
        ? radial_solutions[4] : eos_data->temperature;

    tidalpy::c_MaterialEOSBase* eos_model = eos_data->eos_model_ptr;

    output->density = eos_model->calc_density(
        pressure_si, temperature, radius_si) / eos_data->density_scale;

    if (ode_args->update_shear)
    {
        output->shear_modulus = std::complex<double>(
            eos_model->calc_static_shear_modulus(radius_si) / eos_data->pascal_scale, 0.0);
        output->shear_viscosity = eos_model->calc_shear_viscosity(radius_si);
    }
    else
    {
        output->shear_modulus   = std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        output->shear_viscosity = TidalPyConstants::d_NAN;
    }
    if (ode_args->update_bulk)
    {
        output->bulk_modulus = std::complex<double>(
            eos_model->calc_static_bulk_modulus(radius_si) / eos_data->pascal_scale, 0.0);
        output->bulk_viscosity = eos_model->calc_bulk_viscosity(radius_si);
    }
    else
    {
        output->bulk_modulus   = std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        output->bulk_viscosity = TidalPyConstants::d_NAN;
    }
}
