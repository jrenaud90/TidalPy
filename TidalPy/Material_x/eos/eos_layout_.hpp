#pragma once
/* The state and evaluation layouts of the whole-planet EOS solve.
 *
 * Held apart from ode_.hpp, which pulls in CyRK, so a consumer that only needs the sizes and indices (a
 * layer's dense-output buffer, say) can have them without that dependency.
 *
 * Structure ODE state: [0] gravity, [1] pressure, [2] mass, [3] moment of inertia, and, when the solve
 * carries temperature, [4] temperature [K] and [5] heat flow through the sphere of this radius [W].
 *
 * Evaluation layout, what c_EOSSolution::call hands back at a radius: the four structure variables, then
 * [4] density, [5] static shear modulus, [6] static bulk modulus, [7] shear viscosity, [8] bulk viscosity,
 * [9] temperature, [10] heat flow, [11] melt fraction.
 *
 * Every value here is frequency independent and real; the moduli are the material's unrelaxed values after
 * its partial-melt model.
 */

#include <cstddef>

static const std::size_t C_EOS_Y_VALUES         = 4;
static const std::size_t C_EOS_THERMAL_Y_VALUES = 6;
static const std::size_t C_EOS_EXTRA_VALUES     = 5;
static const std::size_t C_EOS_DY_VALUES        = 12;
static const std::size_t C_EOS_DENSITY_INDEX         = 4;
static const std::size_t C_EOS_SHEAR_MODULUS_INDEX   = 5;
static const std::size_t C_EOS_BULK_MODULUS_INDEX    = 6;
static const std::size_t C_EOS_SHEAR_VISCOSITY_INDEX = 7;
static const std::size_t C_EOS_BULK_VISCOSITY_INDEX  = 8;
static const std::size_t C_EOS_TEMPERATURE_INDEX     = 9;
static const std::size_t C_EOS_HEAT_FLOW_INDEX       = 10;
static const std::size_t C_EOS_MELT_FRACTION_INDEX   = 11;
