#pragma once

#include <complex>

#include "c_common.hpp"    // CyRK: PreEvalFunc
#include "constants_.hpp"
#include "eos_layout_.hpp"



static const double C_FOUR_PI = 4.0 * TidalPyConstants::d_PI;


struct c_EOSOutput
{
    double density                        = TidalPyConstants::d_NAN;
    std::complex<double> bulk_modulus     = {TidalPyConstants::d_NAN, 0.0};
    std::complex<double> shear_modulus    = {TidalPyConstants::d_NAN, 0.0};
    double shear_viscosity                = TidalPyConstants::d_NAN;
    double bulk_viscosity                 = TidalPyConstants::d_NAN;
    double melt_fraction                  = 0.0;
};

struct c_EOSMaterialState
{
    double gravity = TidalPyConstants::d_NAN;
    double density = TidalPyConstants::d_NAN;
    std::complex<double> shear_modulus = {TidalPyConstants::d_NAN, TidalPyConstants::d_NAN};
    std::complex<double> bulk_modulus  = {TidalPyConstants::d_NAN, TidalPyConstants::d_NAN};
};


/// How a radial segment sets its temperature gradient.
enum class c_TemperatureKind : uint8_t
{
    Isothermal = 0,   // dT/dr = 0
    Conductive = 1,   // dT/dr = -conduction_coeff * L / r^2 (Fourier's law in a spherical shell)
    Adiabatic  = 2,   // dT/dr = -adiabat_coeff * g * T
};

/// A stretch of a layer over which the temperature gradient keeps one form. A layer is one segment unless
/// its temperature profile has a kink, which is where the adaptive stepper would otherwise lose its order.
struct c_EOSSegment
{
    double            upper_radius      = 0.0;                            // segment top [solve units]
    size_t            layer_index       = 0;                              // the layer this segment belongs to
    c_TemperatureKind temperature_kind  = c_TemperatureKind::Isothermal;
    double            start_temperature = TidalPyConstants::d_NAN;        // [K]; NaN continues from below
    double            start_heat_flow   = 0.0;                            // [W] entering the segment's base
    double            conduction_coeff  = 0.0;                            // 1 / (4 pi k length_scale)
    double            adiabat_coeff     = 0.0;                            // alpha g_scale length_scale / c_p
};

/// Heat generated inside the planet, as the thermal structure ODE reads it. Abstract so this header stays
/// free of the layer and world classes; the world implements it from its heat sources.
class c_EOSHeatingBase
{
public:
    virtual ~c_EOSHeatingBase() = default;

    /// dL/dr = 4 pi r^2 h, in the units the solve runs in.
    virtual double calc_heat_flow_gradient(size_t layer_index, double radius, double density) const noexcept = 0;
};

/// The temperature fields are set per segment by the solver.
struct c_EOS_ODEInput
{
    double G_to_use       = 0.0;
    double planet_radius  = 0.0;
    char*  eos_input_ptr  = nullptr;
    bool   update_bulk    = false;
    bool   update_shear   = false;
    c_TemperatureKind temperature_kind = c_TemperatureKind::Isothermal;
    double conduction_coeff = 0.0;
    double adiabat_coeff    = 0.0;
    // Heat sources of a thermal solve, non-owning and null for none.
    const c_EOSHeatingBase* heating_ptr = nullptr;
    size_t layer_index = 0;
};


/// The four structure derivatives of a self-gravitating spherically symmetric body in hydrostatic
/// equilibrium. Returns the local density, which the thermal ODE needs for its heat sources.
inline double c_eos_structure_derivatives(
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

    // The gravity equation has a removable 1/r singularity at the origin. There g = (4/3) pi G rho r, so its
    // derivative is the limit (4/3) pi G rho (a zero there collapses the integrator's first step and costs it tens
    // of steps to recover), while pressure, mass, and moment of inertia start flat. Everything is still past the
    // surface.
    if (radius > eos_input_ptr->planet_radius)
    {
        dy_ptr[0] = 0.0;
        dy_ptr[1] = 0.0;
        dy_ptr[2] = 0.0;
        dy_ptr[3] = 0.0;
    }
    else if (radius < TidalPyConstants::d_EPS_10)
    {
        dy_ptr[0] = grav_coeff * rho / 3.0;
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
    return rho;
}


/// CyRK DiffeqFuncType signature.
inline void c_eos_diffeq(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* input_args,
        PreEvalFunc eos_function) noexcept
{
    c_eos_structure_derivatives(dy_ptr, radius, y_ptr, input_args, eos_function);
}


/// As c_eos_diffeq, with the density evaluated at the local temperature when the layer's EOS is thermal,
/// plus two extra states:
///   dT/dr = 0, -conduction_coeff L / r^2, or -adiabat_coeff g T, by the segment's temperature kind, and
///   dL/dr = 4 pi r^2 h, the heat the world's sources generate at this radius.
/// The temperature is in Kelvin and the heat flow in Watts whatever units the rest of the solve runs in;
/// the coefficients carry the conversion.
inline void c_eos_diffeq_thermal(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* input_args,
        PreEvalFunc eos_function) noexcept
{
    c_EOS_ODEInput* eos_input_ptr = reinterpret_cast<c_EOS_ODEInput*>(input_args);

    const double rho = c_eos_structure_derivatives(dy_ptr, radius, y_ptr, input_args, eos_function);

    if ((radius < TidalPyConstants::d_EPS_10) || (radius > eos_input_ptr->planet_radius))
    {
        dy_ptr[4] = 0.0;
        dy_ptr[5] = 0.0;
        return;
    }

    switch (eos_input_ptr->temperature_kind)
    {
        case c_TemperatureKind::Conductive:
            dy_ptr[4] = -eos_input_ptr->conduction_coeff * y_ptr[5] / (radius * radius);
            break;
        case c_TemperatureKind::Adiabatic:
            dy_ptr[4] = -eos_input_ptr->adiabat_coeff * y_ptr[0] * y_ptr[4];
            break;
        default:
            dy_ptr[4] = 0.0;
            break;
    }
    dy_ptr[5] = (eos_input_ptr->heating_ptr != nullptr)
        ? eos_input_ptr->heating_ptr->calc_heat_flow_gradient(eos_input_ptr->layer_index, radius, rho) : 0.0;
}
