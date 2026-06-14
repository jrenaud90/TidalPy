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
 * All quantities are MKS. The material EOS models do not (yet) provide complex
 * moduli, so the shear/bulk modulus outputs are set to NaN — the world-level
 * solve consumes only the density (plus the integrated gravity / pressure /
 * mass / moment of inertia) when populating each layer's EOS profile.
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

    // Density from the layer's material EOS model.
    output->density = eos_data->eos_model_ptr->calc_density(
        pressure_pa, eos_data->temperature_k, radius);

    // Material EOS models do not provide complex moduli; mark them NaN.
    output->bulk_modulus  = std::complex<double>(TidalPyConstants::d_NAN, 0.0);
    output->shear_modulus = std::complex<double>(TidalPyConstants::d_NAN, 0.0);
}
