#pragma once

#include <complex>

#include "c_common.hpp"   // CyRK: PreEvalFunc typedef
#include "constants_.hpp" // TidalPy: TidalPyConstants


// EOS ODE layout: y = [0] gravity, [1] pressure, [2] mass, [3] moment of inertia; extras tracked on the final
// solve = [4] density, [5, 6] static shear modulus (real, imag), [7, 8] static bulk modulus (real, imag),
// [9] shear viscosity, [10] bulk viscosity. The extras come from the layer's EOS model (NaN for analytic models).
static const size_t C_EOS_Y_VALUES     = 4;
static const size_t C_EOS_EXTRA_VALUES = 7;
static const size_t C_EOS_DY_VALUES    = 11;

static const double C_FOUR_PI = 4.0 * TidalPyConstants::d_PI;


/// EOS evaluation output at a radius.
struct c_EOSOutput
{
    double density                        = TidalPyConstants::d_NAN;
    std::complex<double> bulk_modulus     = {TidalPyConstants::d_NAN, 0.0};
    std::complex<double> shear_modulus    = {TidalPyConstants::d_NAN, 0.0};
    double shear_viscosity                = TidalPyConstants::d_NAN;
    double bulk_viscosity                 = TidalPyConstants::d_NAN;
};


/// Input parameters for the EOS ODE solver.
struct c_EOS_ODEInput
{
    double G_to_use       = 0.0;
    double planet_radius  = 0.0;
    char*  eos_input_ptr  = nullptr;
    bool   final_solve    = false;
    bool   update_bulk    = false;
    bool   update_shear   = false;
};


/// Hydrostatic structure ODE (CyRK DiffeqFuncType signature) for a self-gravitating spherically symmetric body.
inline void c_eos_diffeq(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* input_args,
        PreEvalFunc eos_function) noexcept
{
    c_EOS_ODEInput* eos_input_ptr = reinterpret_cast<c_EOS_ODEInput*>(input_args);

    const double r2 = radius * radius;
    const double grav_coeff = C_FOUR_PI * eos_input_ptr->G_to_use;

    c_EOSOutput eos_output;
    eos_function(
        reinterpret_cast<char*>(&eos_output),
        radius,
        y_ptr,
        input_args
    );

    const double rho = eos_output.density;

    // The gravity equation has a 1/r singularity: all derivatives are zero at the origin and beyond the planet.
    if ((radius < TidalPyConstants::d_EPS_10) || (radius > eos_input_ptr->planet_radius))
    {
        dy_ptr[0] = 0.0;
        dy_ptr[1] = 0.0;
        dy_ptr[2] = 0.0;
        dy_ptr[3] = 0.0;
    }
    else
    {
        dy_ptr[0] = grav_coeff * rho - 2.0 * y_ptr[0] * (1.0 / radius);
        dy_ptr[1] = -rho * y_ptr[0];
        dy_ptr[2] = C_FOUR_PI * rho * r2;
        dy_ptr[3] = (2.0 / 3.0) * dy_ptr[2] * r2;
    }

    // The extras are stored on the final solve only and are not integrated.
    if (eos_input_ptr->final_solve)
    {
        dy_ptr[4] = eos_output.density;

        dy_ptr[5] = eos_output.shear_modulus.real();
        dy_ptr[6] = eos_output.shear_modulus.imag();

        dy_ptr[7] = eos_output.bulk_modulus.real();
        dy_ptr[8] = eos_output.bulk_modulus.imag();

        dy_ptr[9]  = eos_output.shear_viscosity;
        dy_ptr[10] = eos_output.bulk_viscosity;
    }
}
