#pragma once
/*
 * material_.hpp — EOS pre-evaluation that dispatches to a layer's material EOS
 * model (c_MaterialEOSBase) during the whole-planet radial structure solve.
 *
 * This is the class-based replacement for the array-interpolation pre-eval in
 * methods/interpolate_.hpp. Instead of reading pre-tabulated density arrays, it
 * asks each layer's attached material EOS model for the density at the current
 * radius and pressure. The radial solver carries pressure as the y[1] state
 * variable, so an analytic density(pressure) model (Birch-Murnaghan, Vinet) is
 * evaluated inline; the interpolated and constant models ignore pressure.
 *
 * All quantities are MKS. Besides density, the model is also asked for the static
 * shear modulus, static bulk modulus, and shear/bulk viscosity at this radius; the
 * interpolated model returns radius-varying values (so they become extra outputs of
 * the EOS ODE and are carried by the CyRK solution), while the analytic/constant
 * models return NaN for these (the world solve then falls back to the layer's
 * constant). The moduli are gated by the solve's update_shear / update_bulk flags
 * (only populated on the final solve pass).
 */

#include <complex>

#include "../ode_.hpp"             // c_EOS_ODEInput, c_EOSOutput, PreEvalFunc signature
#include "../material_eos_.hpp"    // tidalpy::c_MaterialEOSBase
#include "../../../constants_.hpp" // TidalPyConstants::d_NAN


/// Input data for the material-EOS pre-evaluation function.
///
/// Holds a non-owning observer pointer to the layer's material EOS model and the
/// (currently unused) temperature passed through to calc_density.
struct c_MaterialEOSInput
{
    tidalpy::c_MaterialEOSBase* eos_model_ptr = nullptr;
    double temperature_k = 0.0;
};


// The function signature matches CyRK's PreEvalFunc:
//   void(char* preeval_output, double radius, double* radial_solutions, char* preeval_input)
inline void c_preeval_material_eos(
        // Values that will be updated by the function
        char* preeval_output,
        // Input that is used by the pre-eval
        double radius,
        double* radial_solutions,
        char* preeval_input
        ) noexcept
{
    // Cast input args to the proper structures.
    c_EOS_ODEInput*     ode_args = reinterpret_cast<c_EOS_ODEInput*>(preeval_input);
    c_MaterialEOSInput* eos_data = reinterpret_cast<c_MaterialEOSInput*>(ode_args->eos_input_ptr);

    // Cast output to the proper structure.
    c_EOSOutput* output = reinterpret_cast<c_EOSOutput*>(preeval_output);

    // y[1] carries the local pressure [Pa] as integrated by the radial solver.
    const double pressure_pa = radial_solutions[1];

    tidalpy::c_MaterialEOSBase* eos_model = eos_data->eos_model_ptr;

    // Density from the layer's material EOS model.
    output->density = eos_model->calc_density(
        pressure_pa, eos_data->temperature_k, radius);

    // Static moduli from the model (radius-varying for the interpolated model; NaN
    // for analytic/constant models). Only populated on the final solve pass.
    if (ode_args->update_shear)
    {
        output->shear_modulus = std::complex<double>(
            eos_model->calc_static_shear_modulus(radius), 0.0);
        output->shear_viscosity = eos_model->calc_shear_viscosity(radius);
    }
    else
    {
        output->shear_modulus   = std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        output->shear_viscosity = TidalPyConstants::d_NAN;
    }
    if (ode_args->update_bulk)
    {
        output->bulk_modulus = std::complex<double>(
            eos_model->calc_static_bulk_modulus(radius), 0.0);
        output->bulk_viscosity = eos_model->calc_bulk_viscosity(radius);
    }
    else
    {
        output->bulk_modulus   = std::complex<double>(TidalPyConstants::d_NAN, 0.0);
        output->bulk_viscosity = TidalPyConstants::d_NAN;
    }
}
