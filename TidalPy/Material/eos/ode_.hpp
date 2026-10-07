#pragma once

#include <complex>

#include "c_common.hpp"    // CyRK: PreEvalFunc
#include "constants_.hpp"
#include "eos_layout_.hpp"
#include "../../Utilities/classes/thermo_point_.hpp"   // c_TemperatureKind

// The EOS integration headers live in the global namespace; the temperature kind is the cooling models' (tidalpy).
using tidalpy::c_TemperatureKind;



static const double C_FOUR_PI = 4.0 * TidalPyConstants::d_PI;


/// What a layer's EOS function reports at a point, in the units of the solve except where marked SI. The structure
/// integration needs only the density; a thermal integration adds the three thermal properties, and a dense
/// evaluation of the finished solution asks for everything.
struct c_EOSOutput
{
    double density                        = TidalPyConstants::d_NAN;
    std::complex<double> bulk_modulus     = {TidalPyConstants::d_NAN, 0.0};   // adiabatic
    std::complex<double> shear_modulus    = {TidalPyConstants::d_NAN, 0.0};
    double shear_viscosity                = TidalPyConstants::d_NAN;   // SI
    double bulk_viscosity                 = TidalPyConstants::d_NAN;   // SI
    double melt_fraction                  = 0.0;
    double thermal_expansion              = TidalPyConstants::d_NAN;   // SI [1/K]
    double heat_capacity                  = TidalPyConstants::d_NAN;   // SI [J kg-1 K-1], latent heat included
    double thermal_conductivity           = TidalPyConstants::d_NAN;   // SI [W m-1 K-1]
    // SI [1/K]: the latent heat's share of an adiabat's expansivity inside a pressure-dependent melting range.
    double latent_expansion = 0.0;
};

struct c_EOSMaterialState
{
    double gravity = TidalPyConstants::d_NAN;
    double density = TidalPyConstants::d_NAN;
    std::complex<double> shear_modulus = {TidalPyConstants::d_NAN, TidalPyConstants::d_NAN};
    std::complex<double> bulk_modulus  = {TidalPyConstants::d_NAN, TidalPyConstants::d_NAN};
};


/// A stretch of a layer over which the temperature gradient keeps one form. A layer is one segment unless
/// its temperature profile has a kink, which is where the adaptive stepper would otherwise lose its order.
///
/// The top is a radius on the geometry the segments were laid out on (c_EOSLayerBounds). A layer that holds its mass
/// ends each segment where its enclosed mass reaches a fraction of that mass instead, so its segments move with it.
struct c_EOSSegment
{
    double            upper_radius        = 0.0;                          // segment top [solve units]
    size_t            layer_index         = 0;                            // the layer this segment belongs to
    double            upper_mass_fraction = 1.0;                          // of a mass-held layer's mass, at the top
    c_TemperatureKind temperature_kind    = c_TemperatureKind::Isothermal;
    double            start_temperature   = TidalPyConstants::d_NAN;      // [K]; NaN continues from below
    double            start_heat_flow     = 0.0;                          // [W] entering the segment's base
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
    // Ask the EOS function for the density, moduli, viscosities, and melt fraction (a dense evaluation of the finished
    // solution) rather than the density alone.
    bool   full_state     = false;
    // Ask it for the thermal properties too (a solve that integrates temperature).
    bool   thermal_state  = false;
    c_TemperatureKind temperature_kind = c_TemperatureKind::Isothermal;
    // The solve's length [m] and gravity [m s-2] per unit, which turn the SI thermal gradients into its units.
    double length_scale  = 1.0;
    double gravity_scale = 1.0;
    // Heat sources of a thermal solve, non-owning and null for none.
    const c_EOSHeatingBase* heating_ptr = nullptr;
    size_t layer_index = 0;
    // The two events that can end a piece of the integration (c_solve_eos): the enclosed mass [solve units] at which a
    // layer holding its mass, or one of its segments, ends, and the shear modulus [Pa] at or below which the
    // material counts as a liquid (the rigidity margin, read by the layer's state event).
    double target_mass           = TidalPyConstants::d_NAN;
    double minimum_shear_modulus = 0.0;
};


/// CyRK EventFunc: zero where the enclosed mass reaches the input's target mass.
inline double c_eos_mass_event(double radius, double* y_ptr, char* input_args) noexcept
{
    (void)radius;
    const c_EOS_ODEInput* eos_input_ptr = reinterpret_cast<const c_EOS_ODEInput*>(input_args);
    return y_ptr[2] - eos_input_ptr->target_mass;
}


/// The four structure derivatives of a self-gravitating spherically symmetric body in hydrostatic
/// equilibrium. Fills `eos_output` with what the EOS function reports at the point, which the thermal ODE reads.
inline void c_eos_structure_derivatives(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* input_args,
        PreEvalFunc eos_function,
        c_EOSOutput& eos_output) noexcept
{
    c_EOS_ODEInput* eos_input_ptr = reinterpret_cast<c_EOS_ODEInput*>(input_args);

    const double r2 = radius * radius;
    const double grav_coeff = C_FOUR_PI * eos_input_ptr->G_to_use;

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
}


/// CyRK DiffeqFuncType signature.
inline void c_eos_diffeq(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* input_args,
        PreEvalFunc eos_function) noexcept
{
    c_EOSOutput eos_output;
    c_eos_structure_derivatives(dy_ptr, radius, y_ptr, input_args, eos_function, eos_output);
}


/// As c_eos_diffeq, with the material evaluated at the local temperature, plus two extra states:
///   dT/dr = 0, -L / (4 pi r^2 k), or -(alpha + alpha_L) g T / c_p, by the segment's temperature kind, with the
///   material's conductivity, expansivity, latent expansion alpha_L, and heat capacity at the point, and
///   dL/dr = 4 pi r^2 h, the heat the world's sources generate at this radius.
/// The temperature is in Kelvin and the heat flow in Watts whatever units the rest of the solve runs in; the input's
/// length and gravity scales carry the conversion.
inline void c_eos_diffeq_thermal(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* input_args,
        PreEvalFunc eos_function) noexcept
{
    c_EOS_ODEInput* eos_input_ptr = reinterpret_cast<c_EOS_ODEInput*>(input_args);

    c_EOSOutput eos_output;
    c_eos_structure_derivatives(dy_ptr, radius, y_ptr, input_args, eos_function, eos_output);

    if ((radius < TidalPyConstants::d_EPS_10) || (radius > eos_input_ptr->planet_radius))
    {
        dy_ptr[4] = 0.0;
        dy_ptr[5] = 0.0;
        return;
    }

    switch (eos_input_ptr->temperature_kind)
    {
        case c_TemperatureKind::Conductive:
        {
            const double conductivity = eos_output.thermal_conductivity;
            dy_ptr[4] = (conductivity > TidalPyConstants::d_EPS)
                ? -y_ptr[5] / (C_FOUR_PI * radius * radius * conductivity * eos_input_ptr->length_scale) : 0.0;
            break;
        }
        case c_TemperatureKind::Adiabatic:
        {
            const double heat_capacity = eos_output.heat_capacity;
            const double expansion = eos_output.thermal_expansion + eos_output.latent_expansion;
            dy_ptr[4] = (heat_capacity > TidalPyConstants::d_EPS)
                ? -expansion * y_ptr[0] * eos_input_ptr->gravity_scale
                    * eos_input_ptr->length_scale * y_ptr[4] / heat_capacity
                : 0.0;
            break;
        }
        default:
            dy_ptr[4] = 0.0;
            break;
    }
    dy_ptr[5] = (eos_input_ptr->heating_ptr != nullptr)
        ? eos_input_ptr->heating_ptr->calc_heat_flow_gradient(eos_input_ptr->layer_index, radius, eos_output.density)
        : 0.0;
}
