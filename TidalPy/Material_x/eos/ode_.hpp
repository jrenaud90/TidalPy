#pragma once

#include <complex>

#include "c_common.hpp"    // CyRK: PreEvalFunc typedef
#include "constants_.hpp"  // TidalPy: TidalPyConstants
#include "eos_layout_.hpp" // TidalPy: the state and evaluation layouts



static const double C_FOUR_PI = 4.0 * TidalPyConstants::d_PI;


/// EOS evaluation output at a radius.
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

/// One radial segment of the structure solve: a stretch of a layer over which the temperature gradient keeps
/// one form. A layer is one segment unless its temperature profile has a kink, which is where the adaptive
/// stepper would otherwise lose its order. The coefficients are already in the units the solve runs in.
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

/// Heat generated inside the planet, as the thermal structure ODE reads it. The world implements this from its
/// heat sources; it is abstract here so this header stays free of the layer and world classes.
class c_EOSHeatingBase
{
public:
    virtual ~c_EOSHeatingBase() = default;

    /// dL/dr = 4 pi r^2 h at a radius of one layer, where the local density is `density`. The radius, the
    /// density, and the length in the returned Watts per length are all in the units the solve runs in.
    virtual double calc_heat_flow_gradient(size_t layer_index, double radius, double density) const noexcept = 0;
};

/// Input parameters for the EOS ODE solver. The temperature fields are set per segment by the solver.
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
    // Heat sources of a thermal solve (non-owning; null for none) and the layer this input belongs to.
    const c_EOSHeatingBase* heating_ptr = nullptr;
    size_t layer_index = 0;
};


/// The four structure derivatives of a self-gravitating spherically symmetric body in hydrostatic equilibrium.
/// Returns the local density, which the thermal ODE needs for its heat sources.
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
    return rho;
}


/// Hydrostatic structure ODE (CyRK DiffeqFuncType signature).
inline void c_eos_diffeq(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* input_args,
        PreEvalFunc eos_function) noexcept
{
    c_eos_structure_derivatives(dy_ptr, radius, y_ptr, input_args, eos_function);
}


/// Hydrostatic structure ODE carrying temperature and heat flow (CyRK DiffeqFuncType signature).
///
/// The four structure derivatives are those of c_eos_diffeq, with the density evaluated at the local
/// temperature when the layer's EOS is thermal. The two extra states are
///   dT/dr = 0, -conduction_coeff L / r^2, or -adiabat_coeff g T, by the segment's temperature kind, and
///   dL/dr = 4 pi r^2 h, the heat generated at this radius by the world's heat sources (zero without any).
/// The temperature is in Kelvin whatever units the rest of the solve runs in; the coefficients carry the
/// conversion. The heat flow is in Watts.
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
