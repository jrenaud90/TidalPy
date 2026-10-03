#pragma once
/* c_ThermoPoint: where a material property is evaluated, the argument every property law (equation of state, shear
 * modulus, viscosity) takes, so no law can be handed its arguments in the wrong order. All MKS.
 *
 * c_TemperatureKind: how a stretch of a layer sets its temperature gradient, shared by the cooling models (which say
 * what each stretch of their layer is) and the structure ODE (which integrates it).
 */

#include <cstdint>

#include "constants_.hpp"

namespace tidalpy {

enum class c_TemperatureKind : uint8_t {
    Isothermal = 0,   // dT/dr = 0
    Conductive = 1,   // dT/dr = -L / (4 pi r^2 k), Fourier's law in a spherical shell
    Adiabatic  = 2,   // dT/dr = -alpha g T / c_p
};

struct c_ThermoPoint {
    double pressure    = 0.0;                        // [Pa]
    double temperature = TidalPyConstants::d_NAN;    // [K]; NaN for a law evaluated without a temperature
    double radius      = TidalPyConstants::d_NAN;    // [m]; read only by laws tabulated in radius
};

} // namespace tidalpy
