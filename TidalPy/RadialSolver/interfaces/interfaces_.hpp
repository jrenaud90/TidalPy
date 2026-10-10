// interfaces_.hpp: interface conditions between layers.
//
// References
// ----------
// S74  : Saito (1974; J. Phy. Earth; DOI: 10.4294/jpe1952.22.123)
// TS72 : Takeuchi & Saito (1972), Seismic surface waves, Methods Comput. Phys., 11, 217-295.
#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

#include "../../constants_.hpp"
#include "../rs_constants_.hpp"   // C_MAX_NUM_SOL
#include "../layer_kind_.hpp"     // c_layer_num_solutions


// One side of a layer interface: its layer's type and assumptions and its material at the interface radius.
struct c_InterfaceSide {
    int    layer_type;   // 0 = solid, 1 = liquid
    bool   is_static;
    double gravity;
    double density;
};

// The material an interface condition reads.
struct c_InterfaceValues {
    double gravity;
    double liquid_density;
};

// Gravity at an interface is the mean of the two sides. The liquid density is the liquid side of a solid-liquid
// pair and the static side of a static-dynamic liquid pair; NaN where neither interface condition reads it (a
// solid-solid pair, and liquid pairs of the same kind, which pass their solutions through). The upward shooting
// integration and the downward collapse both take their interface values from here.
inline c_InterfaceValues c_interface_values(const c_InterfaceSide& lower, const c_InterfaceSide& upper) noexcept
{
    const bool lower_solid = (lower.layer_type == 0);
    const bool upper_solid = (upper.layer_type == 0);

    double liquid_density = TidalPyConstants::d_NAN;
    if (lower_solid && !upper_solid)
    {
        liquid_density = upper.density;
    }
    else if (!lower_solid && upper_solid)
    {
        liquid_density = lower.density;
    }
    else if (!lower_solid && !upper_solid)
    {
        if (upper.is_static && !lower.is_static)
        {
            liquid_density = upper.density;
        }
        else if (lower.is_static && !upper.is_static)
        {
            liquid_density = lower.density;
        }
    }
    return c_InterfaceValues{0.5 * (upper.gravity + lower.gravity), liquid_density};
}


// The interface's linear map between the two layers' solutions, [upper_solution * C_MAX_NUM_SOL + lower_solution]:
// each upper solution that continues the lower ones is the combination sum_i T[j][i] y_lower_i of them, and a row of
// zeros is a solution that starts at the interface (y3 of a solid over a liquid, say). The same map runs the collapse
// downward: the lower layer's constants are c_lower_i = sum_j c_upper_j T[j][i] (c_collapse_through_interface).
inline constexpr std::size_t C_INTERFACE_TRANSFER_SIZE = C_MAX_NUM_SOL * C_MAX_NUM_SOL;

// Starting y of the upper layer's solutions from the top-of-lower-layer y, by layer type and static/dynamic
// pairing (TS72 Eqs. 140-149; S74 Eqs. 20-21). Unused entries are NaN. transfer_ptr, when given, receives the map
// above (C_INTERFACE_TRANSFER_SIZE entries).
inline void c_solve_upper_y_at_interface(
        std::complex<double>* lower_layer_y_ptr,
        std::complex<double>* upper_layer_y_ptr,
        size_t num_sols_lower,
        size_t num_sols_upper,
        size_t max_num_y,
        int lower_layer_type,
        bool lower_is_static,
        int upper_layer_type,
        bool upper_is_static,
        double interface_gravity,
        double liquid_density,
        double G_to_use,
        std::complex<double>* transfer_ptr = nullptr) noexcept
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

    // Compressibility does not enter: the conditions are the continuity of the ys across the interface (and a liquid's
    // vanishing shear stress), which hold for compressible and incompressible layers alike.

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

    // The caller's buffer holds num_sols_upper rows of max_num_y.
    for (size_t yi_upper = 0; yi_upper < (num_sols_upper * max_num_y); ++yi_upper)
    {
        upper_layer_y_ptr[yi_upper] = cmplx_NAN;
    }
    // Every entry of the map not set below is zero; a caller that wants none writes into a local.
    std::array<std::complex<double>, C_INTERFACE_TRANSFER_SIZE> transfer_local{};
    std::complex<double>* transfer = (transfer_ptr != nullptr) ? transfer_ptr : transfer_local.data();
    for (size_t entry_i = 0; entry_i < C_INTERFACE_TRANSFER_SIZE; ++entry_i)
    {
        transfer[entry_i] = cmplx_zero;
    }
    const auto set_transfer = [transfer](size_t upper_solution, size_t lower_solution, std::complex<double> value)
    {
        transfer[upper_solution * C_MAX_NUM_SOL + lower_solution] = value;
    };

    // Solid-solid passes every y through, static or dynamic, and so does a liquid pair of the same kind.
    if (solid_solid || (liquid_liquid && (static_static || dynamic_dynamic)))
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
        for (size_t soli = 0; soli < num_sols_lower; ++soli) { set_transfer(soli, soli, cmplx_one); }
    } else if (liquid_liquid)
    {
        if (static_dynamic)
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
            // Solution 1 continues the static liquid's; solution 2 starts here.
            set_transfer(0, 0, cmplx_one);
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
            set_transfer(0, 0, coeff_1);
            set_transfer(0, 1, coeff_2);
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
            // Solutions 1 and 2 continue the liquid's; solution 3 (y3) starts here.
            set_transfer(0, 0, cmplx_one);
            set_transfer(1, 1, cmplx_one);
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
            // Solution 1 continues the static liquid's; solutions 2 and 3 start here.
            set_transfer(0, 0, cmplx_one);
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
                set_transfer(soli_upper, soli_upper, cmplx_one);
                set_transfer(soli_upper, 2, -coeff_1);
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
            set_transfer(0, 0, coeff_1);
            set_transfer(0, 1, coeff_2);
            set_transfer(0, 2, coeff_3);
        }
    }
}


// The lower layer's constants from the upper layer's, through the interface's map (c_solve_upper_y_at_interface's
// transfer): c_lower_i = sum_j c_upper_j T[j][i] over the upper layer's solutions. This is the shooting method's
// collapse at an interface, the same algebra as the upward pass by construction.
inline void c_collapse_through_interface(
        std::complex<double>* lower_constants_ptr,
        const std::complex<double>* upper_constants_ptr,
        const std::complex<double>* transfer_ptr,
        size_t num_sols_lower,
        size_t num_sols_upper) noexcept
{
    for (size_t lower_i = 0; lower_i < num_sols_lower; ++lower_i)
    {
        std::complex<double> constant(0.0, 0.0);
        for (size_t upper_j = 0; upper_j < num_sols_upper; ++upper_j)
        {
            constant += upper_constants_ptr[upper_j] * transfer_ptr[upper_j * C_MAX_NUM_SOL + lower_i];
        }
        lower_constants_ptr[lower_i] = constant;
    }
}


// A layer's constants from the constants of the layer above, given the layer's own solutions at its top (TS72 Eqs.
// 142-149; S74 Eqs. 20-21): the interface map is rebuilt from those solutions (c_solve_upper_y_at_interface) and
// applied downward (c_collapse_through_interface). The shooting method records the map on its upward pass instead;
// this is the standalone form the Python wrapper exposes. num_sols is this layer's, max_num_y the row stride.
inline void c_top_to_bottom_interface_bc(
        std::complex<double>* constant_vector_ptr,
        const std::complex<double>* layer_above_constant_vector_ptr,
        std::complex<double>* uppermost_y_per_solution_ptr,
        double gravity_upper,
        double layer_above_lower_gravity,
        double density_upper,
        double layer_above_lower_density,
        int layer_type,
        int layer_above_type,
        bool layer_is_static,
        bool layer_above_is_static,
        size_t num_sols,
        size_t max_num_y)
{
    // This layer is the lower side of the interface and the layer above the upper side.
    const c_InterfaceValues interface_values = c_interface_values(
        c_InterfaceSide{layer_type, layer_is_static, gravity_upper, density_upper},
        c_InterfaceSide{layer_above_type, layer_above_is_static, layer_above_lower_gravity, layer_above_lower_density});
    const size_t num_sols_above = c_layer_num_solutions(layer_above_type, layer_above_is_static);
    std::vector<std::complex<double>> upper_y(C_MAX_NUM_SOL * max_num_y);
    std::array<std::complex<double>, C_INTERFACE_TRANSFER_SIZE> transfer{};
    // The map does not depend on G (only the upper ys it is built with do), so any G serves.
    const double any_G = 1.0;
    c_solve_upper_y_at_interface(
        uppermost_y_per_solution_ptr, upper_y.data(), num_sols, num_sols_above, max_num_y, layer_type, layer_is_static,
        layer_above_type, layer_above_is_static, interface_values.gravity, interface_values.liquid_density, any_G,
        transfer.data());
    c_collapse_through_interface(
        constant_vector_ptr, layer_above_constant_vector_ptr, transfer.data(), num_sols, num_sols_above);
    for (size_t solution_i = num_sols; solution_i < C_MAX_NUM_SOL; ++solution_i)
    {
        constant_vector_ptr[solution_i] = std::complex<double>(
            std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
    }
}
