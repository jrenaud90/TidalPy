#pragma once
/* c_ThermoPoint: where a material property is evaluated, the argument every property law (equation of state, shear
 * modulus, viscosity) takes, so no law can be handed its arguments in the wrong order. All MKS.
 */

#include "constants_.hpp"

namespace tidalpy {

struct c_ThermoPoint {
    double pressure    = 0.0;                        // [Pa]
    double temperature = TidalPyConstants::d_NAN;    // [K]; NaN for a law evaluated without a temperature
    double radius      = TidalPyConstants::d_NAN;    // [m]; read only by laws tabulated in radius
};

} // namespace tidalpy
