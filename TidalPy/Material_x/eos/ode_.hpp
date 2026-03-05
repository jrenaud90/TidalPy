#pragma once

#include <complex>

#include "c_common.hpp"   // CyRK: PreEvalFunc typedef
#include "constants_.hpp" // TidalPy: TidalPyConstants


// Number of dependent y values in the EOS ODE
static const size_t C_EOS_Y_VALUES     = 4;
// Number of extra parameters tracked during final solve
static const size_t C_EOS_EXTRA_VALUES = 5;
// Total number of dy values (y + extra)
static const size_t C_EOS_DY_VALUES    = 9;
// EOS ODE dependent y-values:
// 0 : gravity
// 1 : pressure
// 2 : total mass
// 3 : moment of inertia
// Additionally there are 5 extra parameters tracked during the final solve:
// 4 : density
// 5 : Real part of complex shear modulus
// 6 : Imag part of complex shear modulus
// 7 : Real part of complex bulk modulus
// 8 : Imag part of complex bulk modulus

static const double C_FOUR_PI = 4.0 * TidalPyConstants::d_PI;


/// Output structure from EOS evaluation at a given radius.
struct c_EOSOutput
{
    double density                        = 0.0;
    std::complex<double> bulk_modulus     = {0.0, 0.0};
    std::complex<double> shear_modulus    = {0.0, 0.0};
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


/// Solve for EOS components as a function of radius.
///
/// Matches CyRK DiffeqFuncType signature:
///   void(double* dy_ptr, double radius, double* y_ptr, char* input_args, PreEvalFunc eos_function)
///
/// References
/// ----------
/// Standard interior structure equations for self-gravitating bodies.
inline void c_eos_diffeq(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* input_args,
        PreEvalFunc eos_function) noexcept
{
    // Cast input args to correct structure type
    c_EOS_ODEInput* eos_input_ptr = reinterpret_cast<c_EOS_ODEInput*>(input_args);

    // Other constants
    const double r2 = radius * radius;
    const double grav_coeff = C_FOUR_PI * eos_input_ptr->G_to_use;

    // Update viscoelastic parameters using the user-provided equation of state
    c_EOSOutput eos_output;
    eos_function(
        reinterpret_cast<char*>(&eos_output),
        radius,
        y_ptr,
        input_args
    );

    const double rho = eos_output.density;

    // Solve for the dependent variables
    // Gravity is proportionate to 1 / r so there is a singularity at r=0.
    // Set all derivatives equal to zero near the origin or beyond the planet.
    if ((radius < TidalPyConstants::d_EPS_10) || (radius > eos_input_ptr->planet_radius))
    {
        // Acceleration due to Gravity
        dy_ptr[0] = 0.0;

        // Pressure
        dy_ptr[1] = 0.0;

        // Total mass
        dy_ptr[2] = 0.0;

        // Moment of inertia
        dy_ptr[3] = 0.0;
    }
    else
    {
        // Acceleration due to Gravity
        dy_ptr[0] = grav_coeff * rho - 2.0 * y_ptr[0] * (1.0 / radius);

        // Pressure
        dy_ptr[1] = -rho * y_ptr[0];

        // Mass and MOI assume spherical symmetry
        // Total mass
        dy_ptr[2] = C_FOUR_PI * rho * r2;

        // Moment of inertia (r^2 multiplied by dm which was found above)
        dy_ptr[3] = (2.0 / 3.0) * dy_ptr[2] * r2;
    }

    // Store other parameters
    // TODO: Track the static shear and bulk as well as the bulk and shear viscosity as additional outputs.
    if (eos_input_ptr->final_solve)
    {
        // There are 4 dependent y values and then 5 additional parameters that are saved but
        // not used during integration but which the user may want for reference.
        dy_ptr[4] = eos_output.density;

        dy_ptr[5] = eos_output.shear_modulus.real();
        dy_ptr[6] = eos_output.shear_modulus.imag();

        dy_ptr[7] = eos_output.bulk_modulus.real();
        dy_ptr[8] = eos_output.bulk_modulus.imag();
    }

    // Done
}
