// reversed_.hpp - Top-to-bottom interface boundary conditions (for collapse phase)
// Ported from TidalPy/RadialSolver/interfaces/reversed.pyx
//
// References
// ----------
// S74  : Saito (1974)
// TS72 : Takeuchi & Saito (1972)
#pragma once

#include <cmath>
#include <complex>
#include <limits>

#include "../../constants_.hpp"


inline void c_top_to_bottom_interface_bc(
        std::complex<double>* constant_vector_ptr,
        std::complex<double>* layer_above_constant_vector_ptr,
        std::complex<double>* uppermost_y_per_solution_ptr,
        double gravity_upper,
        double layer_above_lower_gravity,
        double density_upper,
        double layer_above_lower_density,
        int layer_type,
        int layer_above_type,
        bool layer_is_static,
        bool layer_above_is_static,
        bool layer_is_incomp,
        bool layer_above_is_incomp,
        size_t num_sols,
        size_t max_num_y) noexcept
{
    const double nan_val = std::numeric_limits<double>::quiet_NaN();
    const std::complex<double> cmplx_NAN(nan_val, nan_val);
    const double g_const = 4.0 * TidalPyConstants::d_PI;  // Note: G_to_use is folded into gravity values

    // Interfaces are defined at the bottom of the layer in question. However, this function is calculating
    // the transition at the top of each layer as it works its way down.
    // So, for interface values, we actually need the ones of the layer above us.
    const double interface_gravity = 0.5 * (gravity_upper + layer_above_lower_gravity);
    double liquid_density_at_interface = nan_val;

    const bool layer_is_solid = (layer_type == 0);
    const bool layer_above_is_solid = (layer_above_type == 0);

    if (!layer_is_solid) {
        if (layer_is_static) {
            liquid_density_at_interface = density_upper;
        } else if (!layer_above_is_solid && layer_above_is_static) {
            liquid_density_at_interface = layer_above_lower_density;
        } else {
            liquid_density_at_interface = density_upper;
        }
    } else if (!layer_above_is_solid) {
        liquid_density_at_interface = layer_above_lower_density;
    }

    // Solve for constant vector
    std::complex<double> y4_frac_1  = cmplx_NAN;
    std::complex<double> y4_frac_2  = cmplx_NAN;
    std::complex<double> gamma_1    = cmplx_NAN;
    std::complex<double> gamma_2    = cmplx_NAN;
    std::complex<double> lambda_1   = cmplx_NAN;
    std::complex<double> lambda_2   = cmplx_NAN;
    std::complex<double> lower_s1y1 = cmplx_NAN;
    std::complex<double> lower_s1y2 = cmplx_NAN;
    std::complex<double> lower_s1y5 = cmplx_NAN;
    std::complex<double> lower_s1y6 = cmplx_NAN;
    std::complex<double> lower_s2y1 = cmplx_NAN;
    std::complex<double> lower_s2y2 = cmplx_NAN;
    std::complex<double> lower_s2y5 = cmplx_NAN;
    std::complex<double> lower_s2y6 = cmplx_NAN;

    if (layer_is_solid) {
        // Solid Layer
        if (layer_above_is_solid) {
            // Both layers are solid. Constants are the same.
            for (size_t solution_i = 0; solution_i < num_sols; ++solution_i) {
                constant_vector_ptr[solution_i] = layer_above_constant_vector_ptr[solution_i];
            }
        } else {
            // Create helper functions
            y4_frac_1 = (
                    -uppermost_y_per_solution_ptr[0 * max_num_y + 3] /
                    uppermost_y_per_solution_ptr[2 * max_num_y + 3]
                );
            y4_frac_2 = (
                    -uppermost_y_per_solution_ptr[1 * max_num_y + 3] /
                    uppermost_y_per_solution_ptr[2 * max_num_y + 3]
                );

            if (layer_above_is_static) {
                // Need to find 3 solid constants from 1 liquid constant
                // S74, Page 131
                constant_vector_ptr[0] = layer_above_constant_vector_ptr[0];
                // Derived by JPR based on Eq 21 (2nd line) of S74
                gamma_1 =
                    (
                        uppermost_y_per_solution_ptr[0 * max_num_y + 1] +
                        y4_frac_1 * uppermost_y_per_solution_ptr[2 * max_num_y + 1]
                    ) -
                    (
                        liquid_density_at_interface *
                        (
                            interface_gravity *
                            (
                                uppermost_y_per_solution_ptr[0 * max_num_y + 0] +
                                y4_frac_1 * uppermost_y_per_solution_ptr[2 * max_num_y + 0]
                            ) -
                            (
                                uppermost_y_per_solution_ptr[0 * max_num_y + 4] +
                                y4_frac_1 * uppermost_y_per_solution_ptr[2 * max_num_y + 4]
                            )
                        )
                    );
                gamma_2 =
                    (
                        uppermost_y_per_solution_ptr[1 * max_num_y + 1] +
                        y4_frac_2 * uppermost_y_per_solution_ptr[2 * max_num_y + 1]
                    ) -
                    (
                        liquid_density_at_interface *
                        (
                            interface_gravity *
                            (
                                uppermost_y_per_solution_ptr[1 * max_num_y + 0] +
                                y4_frac_2 * uppermost_y_per_solution_ptr[2 * max_num_y + 0]
                            ) -
                            (
                                uppermost_y_per_solution_ptr[1 * max_num_y + 4] +
                                y4_frac_2 * uppermost_y_per_solution_ptr[2 * max_num_y + 4]
                            )
                        )
                    );

                constant_vector_ptr[1] = (-gamma_1 / gamma_2) * constant_vector_ptr[0];
                // TS72, Eq. 142 (utilizes y_4 = 0)
                constant_vector_ptr[2] = y4_frac_1 * constant_vector_ptr[0] + y4_frac_2 * constant_vector_ptr[1];

            } else {
                // Need to find 3 solid constants from 2 liquid constants
                // TS72, Eq. 144
                constant_vector_ptr[0] = layer_above_constant_vector_ptr[0];
                constant_vector_ptr[1] = layer_above_constant_vector_ptr[1];
                // TS72, Eq. 142 (utilizes y_4 = 0)
                constant_vector_ptr[2] = y4_frac_1 * constant_vector_ptr[0] + y4_frac_2 * constant_vector_ptr[1];
            }
        }
    } else {
        if (layer_is_static) {
            if (!layer_above_is_solid) {
                // Liquid layer above
                if (layer_above_is_static) {
                    constant_vector_ptr[0] = layer_above_constant_vector_ptr[0];
                } else {
                    constant_vector_ptr[0] = layer_above_constant_vector_ptr[0];
                }
            } else {
                // Solid layer above
                constant_vector_ptr[0] = layer_above_constant_vector_ptr[0];
            }
        } else {
            if (!layer_above_is_solid) {
                // Liquid layer above
                if (layer_above_is_static) {
                    // Need to find 2 liquid (dynamic) constants from 1 liquid (static) constant
                    constant_vector_ptr[0] = layer_above_constant_vector_ptr[0];
                    // Pull out ys
                    lower_s1y1 = uppermost_y_per_solution_ptr[0 * max_num_y + 0];
                    lower_s1y2 = uppermost_y_per_solution_ptr[0 * max_num_y + 1];
                    lower_s1y5 = uppermost_y_per_solution_ptr[0 * max_num_y + 2];
                    lower_s1y6 = uppermost_y_per_solution_ptr[0 * max_num_y + 3];
                    lower_s2y1 = uppermost_y_per_solution_ptr[1 * max_num_y + 0];
                    lower_s2y2 = uppermost_y_per_solution_ptr[1 * max_num_y + 1];
                    lower_s2y5 = uppermost_y_per_solution_ptr[1 * max_num_y + 2];
                    lower_s2y6 = uppermost_y_per_solution_ptr[1 * max_num_y + 3];
                    // lambda_j = (y_2j - rho * ( g * y_1j - y_5j))
                    lambda_1 = lower_s1y2 - liquid_density_at_interface *
                            (interface_gravity * lower_s1y1 - lower_s1y5);
                    lambda_2 = lower_s2y2 - liquid_density_at_interface *
                            (interface_gravity * lower_s2y1 - lower_s2y5);
                    constant_vector_ptr[1] = (-lambda_1 / lambda_2) * constant_vector_ptr[0];
                } else {
                    // Both layers are dynamic liquids. Constants are the same.
                    for (size_t solution_i = 0; solution_i < num_sols; ++solution_i) {
                        constant_vector_ptr[solution_i] = layer_above_constant_vector_ptr[solution_i];
                    }
                }
            } else {
                // Solid layer above
                // TS72 Eqs. 148-149
                constant_vector_ptr[0] = layer_above_constant_vector_ptr[0];
                constant_vector_ptr[1] = layer_above_constant_vector_ptr[1];
            }
        }
    }
}
