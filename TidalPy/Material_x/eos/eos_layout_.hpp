#pragma once
/*
 * eos_layout_.hpp: the state and evaluation layouts of the whole-planet EOS solve.
 *
 * Held apart from ode_.hpp, which pulls in CyRK, so a consumer that only needs the sizes and indices (a layer's
 * dense-output buffer, for instance) can have them without that dependency.
 *
 * Structure ODE state: [0] gravity, [1] pressure, [2] mass, [3] moment of inertia, and, when the solve carries
 * temperature, [4] temperature [K] and [5] heat flow through the sphere of this radius [W].
 *
 * Evaluation layout (what c_EOSSolution hands back at a radius): the four structure variables, then [4] density,
 * [5, 6] static shear modulus (real, imag), [7, 8] static bulk modulus (real, imag), [9] shear viscosity,
 * [10] bulk viscosity, [11] temperature, [12] heat flow.
 */

#include <cstddef>

static const std::size_t C_EOS_Y_VALUES         = 4;
static const std::size_t C_EOS_THERMAL_Y_VALUES = 6;
static const std::size_t C_EOS_EXTRA_VALUES     = 7;
static const std::size_t C_EOS_DY_VALUES        = 13;
// Evaluation-layout indices of the two thermal outputs.
static const std::size_t C_EOS_TEMPERATURE_INDEX = 11;
static const std::size_t C_EOS_HEAT_FLOW_INDEX   = 12;
