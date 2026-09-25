#pragma once
/* EOS pre-evaluation: asks a layer's material for its properties at the current radius, pressure, and
 * temperature during the whole-planet structure solve. While the solve iterates it needs the density alone;
 * the full frequency-independent state follows only when the update_shear / update_bulk flags are set, which
 * is every dense evaluation of the finished solution.
 */

#include <complex>

#include "../ode_.hpp"
#include "../material_eos_.hpp"
#include "../../../constants_.hpp"


/// The models work in SI, so a non-dimensional solve scales the radius and pressure up before the call and
/// the density and moduli down afterwards; the viscosities come back in SI either way.
struct c_MaterialEOSInput
{
    tidalpy::c_MaterialEOSBase* eos_model_ptr = nullptr;
    // A solve that integrates temperature takes the local state value instead.
    double temperature   = 0.0;
    bool   use_state_temperature = false;
    // Whether the density law sees the temperature; the viscosity and partial-melt models always do.
    bool   thermal_density = false;
    double length_scale  = 1.0;
    double pascal_scale  = 1.0;
    double density_scale = 1.0;
};


/// CyRK PreEvalFunc signature.
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

    if (!(ode_args->update_shear || ode_args->update_bulk))
    {
        // The structure iteration needs the density alone, so skip the viscosity and melt models.
        const double density_temperature = eos_data->thermal_density ? temperature : TidalPyConstants::d_NAN;
        output->density = eos_model->calc_density(
            pressure_si, density_temperature, radius_si) / eos_data->density_scale;
        output->shear_modulus   = std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        output->bulk_modulus    = std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        output->shear_viscosity = TidalPyConstants::d_NAN;
        output->bulk_viscosity  = TidalPyConstants::d_NAN;
        output->melt_fraction   = TidalPyConstants::d_NAN;
        return;
    }

    tidalpy::c_MaterialState state;
    eos_model->calc_material_state(pressure_si, temperature, eos_data->thermal_density, radius_si, state);
    output->density         = state.density / eos_data->density_scale;
    output->shear_modulus   = std::complex<double>(state.shear_modulus / eos_data->pascal_scale, 0.0);
    output->bulk_modulus    = std::complex<double>(state.bulk_modulus / eos_data->pascal_scale, 0.0);
    output->shear_viscosity = state.shear_viscosity;
    output->bulk_viscosity  = state.bulk_viscosity;
    output->melt_fraction   = state.melt_fraction;
}
