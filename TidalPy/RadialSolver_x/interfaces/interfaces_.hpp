// interfaces_.hpp - Interface conditions between layers
// Ported from TidalPy/RadialSolver/interfaces/interfaces.pyx
//
// References
// ----------
// S74  : Saito (1974; J. Phy. Earth; DOI: 10.4294/jpe1952.22.123)
// TS72 : Takeuchi & Saito (1972), Seismic surface waves, Methods Comput. Phys., 11, 217-295.
#pragma once

#include <cmath>
#include <complex>
#include <limits>

#include "../../constants_.hpp"


inline void c_solve_upper_y_at_interface(
        std::complex<double>* lower_layer_y_ptr,
        std::complex<double>* upper_layer_y_ptr,
        size_t num_sols_lower,
        size_t num_sols_upper,
        size_t max_num_y,
        int lower_layer_type,
        bool lower_is_static,
        bool lower_is_incompressible,
        int upper_layer_type,
        bool upper_is_static,
        bool upper_is_incompressible,
        double interface_gravity,
        double liquid_density,
        double G_to_use) noexcept
{
    const double nan_val = std::numeric_limits<double>::quiet_NaN();
    const std::complex<double> cmplx_NAN(nan_val, nan_val);
    const std::complex<double> cmplx_zero(0.0, 0.0);
    const std::complex<double> cmplx_one(1.0, 0.0);

    const bool upper_solid = (upper_layer_type == 0);
    const bool lower_solid = (lower_layer_type == 0);

    const bool solid_solid   = lower_solid && upper_solid;
    const bool solid_liquid  = lower_solid && !upper_solid;
    const bool liquid_solid  = !lower_solid && upper_solid;
    const bool liquid_liquid = !lower_solid && !upper_solid;

    const bool static_static   = lower_is_static && upper_is_static;
    const bool static_dynamic  = lower_is_static && !upper_is_static;
    const bool dynamic_static  = !lower_is_static && upper_is_static;
    const bool dynamic_dynamic = !lower_is_static && !upper_is_static;

    // TODO: Compressibility is not currently taken into account. This needs to be checked asap!

    // Other constants that may be needed
    std::complex<double> lambda_1 = cmplx_NAN;
    std::complex<double> lambda_2 = cmplx_NAN;
    std::complex<double> coeff_1  = cmplx_NAN;
    std::complex<double> coeff_2  = cmplx_NAN;
    std::complex<double> coeff_3  = cmplx_NAN;
    std::complex<double> coeff_4  = cmplx_NAN;
    std::complex<double> coeff_5  = cmplx_NAN;
    std::complex<double> coeff_6  = cmplx_NAN;
    std::complex<double> frac_1   = cmplx_NAN;
    std::complex<double> frac_2   = cmplx_NAN;
    std::complex<double> const_1  = cmplx_NAN;

    const double g_const = 4.0 * TidalPyConstants::d_PI * G_to_use;

    // Initialize upper y to nan
    for (size_t yi_upper = 0; yi_upper < 18; ++yi_upper)
    {
        upper_layer_y_ptr[yi_upper] = cmplx_NAN;
    }

    if (solid_solid)
    {
        // Does not matter if the layers are static or dynamic, solid-solid exchange perfectly.
        for (size_t yi_lower = 0; yi_lower < max_num_y; ++yi_lower)
        {
            size_t yi_upper = yi_lower;
            for (size_t soli_lower = 0; soli_lower < num_sols_lower; ++soli_lower)
            {
                size_t soli_upper = soli_lower;
                upper_layer_y_ptr[soli_upper * max_num_y + yi_upper] =
                    lower_layer_y_ptr[soli_lower * max_num_y + yi_lower];
            }
        }
    } else if (liquid_liquid)
    {
        if (static_static || dynamic_dynamic)
        {
            for (size_t yi_lower = 0; yi_lower < max_num_y; ++yi_lower)
            {
                size_t yi_upper = yi_lower;
                for (size_t soli_lower = 0; soli_lower < num_sols_lower; ++soli_lower)
                {
                    size_t soli_upper = soli_lower;
                    upper_layer_y_ptr[soli_upper * max_num_y + yi_upper] =
                        lower_layer_y_ptr[soli_lower * max_num_y + yi_lower];
                }
            }
        } else if (static_dynamic)
        {
            // Solution 1
            upper_layer_y_ptr[0] = cmplx_zero;
            upper_layer_y_ptr[1] = -liquid_density * lower_layer_y_ptr[0];
            upper_layer_y_ptr[2] = lower_layer_y_ptr[0];
            upper_layer_y_ptr[3] =
                lower_layer_y_ptr[1] + (g_const * liquid_density / interface_gravity) * lower_layer_y_ptr[0];

            // Solution 2
            upper_layer_y_ptr[1 * max_num_y + 0] = cmplx_one;
            upper_layer_y_ptr[1 * max_num_y + 1] =
                liquid_density * interface_gravity * upper_layer_y_ptr[1 * max_num_y + 0];
            upper_layer_y_ptr[1 * max_num_y + 2] = cmplx_zero;
            upper_layer_y_ptr[1 * max_num_y + 3] =
                -g_const * liquid_density * upper_layer_y_ptr[1 * max_num_y + 0];
        } else if (dynamic_static)
        {
            lambda_1 =
                lower_layer_y_ptr[1] -
                liquid_density * (interface_gravity * lower_layer_y_ptr[0] - lower_layer_y_ptr[2]);
            lambda_2 =
                lower_layer_y_ptr[1 * max_num_y + 1] -
                liquid_density * (interface_gravity * lower_layer_y_ptr[1 * max_num_y + 0] -
                                  lower_layer_y_ptr[1 * max_num_y + 2]);

            coeff_1 = cmplx_one;
            coeff_2 = -(lambda_1 / lambda_2) * coeff_1;

            coeff_3 = lower_layer_y_ptr[3] + g_const * lower_layer_y_ptr[1] / interface_gravity;
            coeff_4 =
                lower_layer_y_ptr[1 * max_num_y + 3] +
                g_const * lower_layer_y_ptr[1 * max_num_y + 1] / interface_gravity;

            upper_layer_y_ptr[0] = coeff_1 * lower_layer_y_ptr[2] + coeff_2 * lower_layer_y_ptr[1 * max_num_y + 2];
            upper_layer_y_ptr[1] = coeff_1 * coeff_3 + coeff_2 * coeff_4;
        }
    } else if (liquid_solid)
    {
        if (dynamic_dynamic || dynamic_static)
        {
            // See Eqs. 148-149 in TS72
            for (size_t soli_upper = 0; soli_upper < 3; ++soli_upper)
            {
                if (soli_upper == 0 || soli_upper == 1) {
                    upper_layer_y_ptr[soli_upper * max_num_y + 0] = lower_layer_y_ptr[soli_upper * max_num_y + 0];
                    upper_layer_y_ptr[soli_upper * max_num_y + 1] = lower_layer_y_ptr[soli_upper * max_num_y + 1];
                    upper_layer_y_ptr[soli_upper * max_num_y + 4] = lower_layer_y_ptr[soli_upper * max_num_y + 2];
                    upper_layer_y_ptr[soli_upper * max_num_y + 5] = lower_layer_y_ptr[soli_upper * max_num_y + 3];

                    upper_layer_y_ptr[soli_upper * max_num_y + 2] = cmplx_zero;
                    upper_layer_y_ptr[soli_upper * max_num_y + 3] = cmplx_zero;
                } else
                {
                    upper_layer_y_ptr[soli_upper * max_num_y + 0] = cmplx_zero;
                    upper_layer_y_ptr[soli_upper * max_num_y + 1] = cmplx_zero;
                    upper_layer_y_ptr[soli_upper * max_num_y + 2] = cmplx_one;
                    upper_layer_y_ptr[soli_upper * max_num_y + 3] = cmplx_zero;
                    upper_layer_y_ptr[soli_upper * max_num_y + 4] = cmplx_zero;
                    upper_layer_y_ptr[soli_upper * max_num_y + 5] = cmplx_zero;
                }
            }
        } else if (static_dynamic || static_static)
        {
            // Eqs. 20 in S74
            upper_layer_y_ptr[0] = cmplx_zero;
            upper_layer_y_ptr[1] = -liquid_density * lower_layer_y_ptr[0];
            upper_layer_y_ptr[2] = cmplx_zero;
            upper_layer_y_ptr[3] = cmplx_zero;
            upper_layer_y_ptr[4] = lower_layer_y_ptr[0];
            upper_layer_y_ptr[5] =
                lower_layer_y_ptr[1] + (g_const * liquid_density / interface_gravity) * lower_layer_y_ptr[0];

            upper_layer_y_ptr[1 * max_num_y + 0] = cmplx_one;
            upper_layer_y_ptr[1 * max_num_y + 1] =
                liquid_density * interface_gravity * upper_layer_y_ptr[1 * max_num_y + 0];
            upper_layer_y_ptr[1 * max_num_y + 2] = cmplx_zero;
            upper_layer_y_ptr[1 * max_num_y + 3] = cmplx_zero;
            upper_layer_y_ptr[1 * max_num_y + 4] = cmplx_zero;
            upper_layer_y_ptr[1 * max_num_y + 5] =
                -g_const * liquid_density * upper_layer_y_ptr[1 * max_num_y + 0];

            upper_layer_y_ptr[2 * max_num_y + 0] = cmplx_zero;
            upper_layer_y_ptr[2 * max_num_y + 1] = cmplx_zero;
            upper_layer_y_ptr[2 * max_num_y + 2] = cmplx_one;
            upper_layer_y_ptr[2 * max_num_y + 3] = cmplx_zero;
            upper_layer_y_ptr[2 * max_num_y + 4] = cmplx_zero;
            upper_layer_y_ptr[2 * max_num_y + 5] = cmplx_zero;
        }
    } else if (solid_liquid)
    {
        if (dynamic_dynamic || static_dynamic)
        {
            // Eqs. 140-144 in TS72
            for (size_t soli_upper = 0; soli_upper < 2; ++soli_upper)
            {
                coeff_1 = lower_layer_y_ptr[soli_upper * max_num_y + 3] / lower_layer_y_ptr[2 * max_num_y + 3];

                upper_layer_y_ptr[soli_upper * max_num_y + 0] =
                    lower_layer_y_ptr[soli_upper * max_num_y + 0] - coeff_1 * lower_layer_y_ptr[2 * max_num_y + 0];
                upper_layer_y_ptr[soli_upper * max_num_y + 1] =
                    lower_layer_y_ptr[soli_upper * max_num_y + 1] - coeff_1 * lower_layer_y_ptr[2 * max_num_y + 1];
                upper_layer_y_ptr[soli_upper * max_num_y + 2] =
                    lower_layer_y_ptr[soli_upper * max_num_y + 4] - coeff_1 * lower_layer_y_ptr[2 * max_num_y + 4];
                upper_layer_y_ptr[soli_upper * max_num_y + 3] =
                    lower_layer_y_ptr[soli_upper * max_num_y + 5] - coeff_1 * lower_layer_y_ptr[2 * max_num_y + 5];
            }
        } else if (dynamic_static || static_static)
        {
            // Eq. 21 in S74
            frac_1 = -lower_layer_y_ptr[0 * max_num_y + 3] / lower_layer_y_ptr[2 * max_num_y + 3];
            frac_2 = -lower_layer_y_ptr[1 * max_num_y + 3] / lower_layer_y_ptr[2 * max_num_y + 3];

            lambda_1 =
                lower_layer_y_ptr[1] + frac_1 * lower_layer_y_ptr[2 * max_num_y + 1] -
                liquid_density * (
                        interface_gravity * (
                            lower_layer_y_ptr[0] + frac_1 * lower_layer_y_ptr[2 * max_num_y + 0]) -
                        (lower_layer_y_ptr[4] + frac_1 * lower_layer_y_ptr[2 * max_num_y + 4])
                );
            lambda_2 =
                lower_layer_y_ptr[1 * max_num_y + 1] + frac_2 * lower_layer_y_ptr[2 * max_num_y + 1] -
                liquid_density * (
                        interface_gravity * (
                            lower_layer_y_ptr[1 * max_num_y + 0] + frac_2 * lower_layer_y_ptr[2 * max_num_y + 0]) -
                        (lower_layer_y_ptr[1 * max_num_y + 4] + frac_2 * lower_layer_y_ptr[2 * max_num_y + 4])
                );

            coeff_1 = std::complex<double>(1.0, 0.0);
            coeff_2 = -(lambda_1 / lambda_2) * coeff_1;
            coeff_3 = frac_1 * coeff_1 + frac_2 * coeff_2;

            const_1 = (g_const / interface_gravity);

            coeff_4 = lower_layer_y_ptr[0 * max_num_y + 5] + const_1 * lower_layer_y_ptr[0 * max_num_y + 1];
            coeff_5 = lower_layer_y_ptr[1 * max_num_y + 5] + const_1 * lower_layer_y_ptr[1 * max_num_y + 1];
            coeff_6 = lower_layer_y_ptr[2 * max_num_y + 5] + const_1 * lower_layer_y_ptr[2 * max_num_y + 1];

            upper_layer_y_ptr[0] =
                coeff_1 * lower_layer_y_ptr[0 * max_num_y + 4] +
                coeff_2 * lower_layer_y_ptr[1 * max_num_y + 4] +
                coeff_3 * lower_layer_y_ptr[2 * max_num_y + 4];
            upper_layer_y_ptr[1] =
                coeff_1 * coeff_4 + coeff_2 * coeff_5 + coeff_3 * coeff_6;
        }
    }
}
