#pragma once

#include <complex>
#include <cmath>

#include "c_common.hpp"        // CyRK: DiffeqFuncType, PreEvalFunc
#include "eos_solution_.hpp"   // TidalPy: c_EOSSolution
#include "../../constants_.hpp"


/// Arguments passed to each radial solver ODE function via the char* args_ptr.
struct c_RadialSolverArgs
{
    double degree_l;
    double lp1;          // l + 1
    double lm1;          // l - 1
    double llp1;         // l * (l + 1)
    double G;            // Gravitational constant
    double grav_coeff;   // 4 * pi * G
    double frequency;
    size_t layer_index;
    c_EOSSolution* eos_solution_ptr;
};


// ============================================================================
//  Helpers: EOS lookup and y-vector packing
// ============================================================================

/// Gravity, density, and the complex moduli at a radius, in one EOS evaluation. The moduli are complex because
/// they are the material's response at this solve's forcing frequency; see c_EOSSolution::call_material.
static inline void c_read_eos(
        c_RadialSolverArgs* rs_args_ptr,
        double radius,
        double& gravity,
        double& density,
        std::complex<double>& shear_modulus,
        std::complex<double>& bulk_modulus
        ) noexcept
{
    c_EOSMaterialState material_state;
    rs_args_ptr->eos_solution_ptr->call_material(rs_args_ptr->layer_index, radius, material_state);

    gravity       = material_state.gravity;
    density       = material_state.density;
    shear_modulus = material_state.shear_modulus;
    bulk_modulus  = material_state.bulk_modulus;
}


/// Read the six solid-layer ys from CyRK's real-pair array.
static inline void c_read_y6(
        double* y_ptr,
        std::complex<double>& y1,
        std::complex<double>& y2,
        std::complex<double>& y3,
        std::complex<double>& y4,
        std::complex<double>& y5,
        std::complex<double>& y6
        ) noexcept
{
    y1 = std::complex<double>(y_ptr[0],  y_ptr[1]);
    y2 = std::complex<double>(y_ptr[2],  y_ptr[3]);
    y3 = std::complex<double>(y_ptr[4],  y_ptr[5]);
    y4 = std::complex<double>(y_ptr[6],  y_ptr[7]);
    y5 = std::complex<double>(y_ptr[8],  y_ptr[9]);
    y6 = std::complex<double>(y_ptr[10], y_ptr[11]);
}

/// Read the four dynamic-liquid ys (y1, y2, y5, y6).
static inline void c_read_y4_liquid(
        double* y_ptr,
        std::complex<double>& y1,
        std::complex<double>& y2,
        std::complex<double>& y5,
        std::complex<double>& y6
        ) noexcept
{
    y1 = std::complex<double>(y_ptr[0], y_ptr[1]);
    y2 = std::complex<double>(y_ptr[2], y_ptr[3]);
    y5 = std::complex<double>(y_ptr[4], y_ptr[5]);
    y6 = std::complex<double>(y_ptr[6], y_ptr[7]);
}

/// Read the two static-liquid ys (y5, y7).
static inline void c_read_y4_static_liquid(
        double* y_ptr,
        std::complex<double>& y5,
        std::complex<double>& y7
        ) noexcept
{
    y5 = std::complex<double>(y_ptr[0], y_ptr[1]);
    y7 = std::complex<double>(y_ptr[2], y_ptr[3]);
}

/// Write six complex dy values to the real-pair output array.
static inline void c_write_dy6(
        double* dy_ptr,
        const std::complex<double>& dy1,
        const std::complex<double>& dy2,
        const std::complex<double>& dy3,
        const std::complex<double>& dy4,
        const std::complex<double>& dy5,
        const std::complex<double>& dy6
        ) noexcept
{
    dy_ptr[0]  = dy1.real(); dy_ptr[1]  = dy1.imag();
    dy_ptr[2]  = dy2.real(); dy_ptr[3]  = dy2.imag();
    dy_ptr[4]  = dy3.real(); dy_ptr[5]  = dy3.imag();
    dy_ptr[6]  = dy4.real(); dy_ptr[7]  = dy4.imag();
    dy_ptr[8]  = dy5.real(); dy_ptr[9]  = dy5.imag();
    dy_ptr[10] = dy6.real(); dy_ptr[11] = dy6.imag();
}

/// Write four complex dy values (dynamic liquid).
static inline void c_write_dy4(
        double* dy_ptr,
        const std::complex<double>& dy1,
        const std::complex<double>& dy2,
        const std::complex<double>& dy5,
        const std::complex<double>& dy6
        ) noexcept
{
    dy_ptr[0] = dy1.real(); dy_ptr[1] = dy1.imag();
    dy_ptr[2] = dy2.real(); dy_ptr[3] = dy2.imag();
    dy_ptr[4] = dy5.real(); dy_ptr[5] = dy5.imag();
    dy_ptr[6] = dy6.real(); dy_ptr[7] = dy6.imag();
}

/// Write two complex dy values (static liquid).
static inline void c_write_dy2(
        double* dy_ptr,
        const std::complex<double>& dy5,
        const std::complex<double>& dy7
        ) noexcept
{
    dy_ptr[0] = dy5.real(); dy_ptr[1] = dy5.imag();
    dy_ptr[2] = dy7.real(); dy_ptr[3] = dy7.imag();
}

// ============================================================================
//  Solid Compressible (static and dynamic)
// ============================================================================

/// Radial derivative equations for a solid, compressible layer. The static form drops the inertial term, which
/// appears only in dy2 and dy4.
///
/// References: TS72 Eq. 82, KMN15 Eqs. 4--9, B15 Eqs. 13--18
template <bool dynamic>
inline void c_solid_compressible_body(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr
        ) noexcept
{
    c_RadialSolverArgs* rs_args_ptr = reinterpret_cast<c_RadialSolverArgs*>(args_ptr);

    double gravity, density;
    std::complex<double> shear_modulus, bulk_modulus;
    c_read_eos(rs_args_ptr, radius, gravity, density, shear_modulus, bulk_modulus);

    std::complex<double> y1, y2, y3, y4, y5, y6;
    c_read_y6(y_ptr, y1, y2, y3, y4, y5, y6);

    const std::complex<double> lame = bulk_modulus - (2.0 / 3.0) * shear_modulus;

    const double r_inverse       = 1.0 / radius;
    const double density_gravity = density * gravity;
    const double grav_term       = rs_args_ptr->grav_coeff * density;

    const std::complex<double> lame_2mu         = lame + 2.0 * shear_modulus;
    const std::complex<double> lame_2mu_inverse = 1.0 / lame_2mu;
    const std::complex<double> two_shear_r_inv  = 2.0 * shear_modulus * r_inverse;
    const std::complex<double> y1_y3_term       = 2.0 * y1 - rs_args_ptr->llp1 * y3;

    const std::complex<double> dy1 =
        lame_2mu_inverse * (
            y1_y3_term * -lame * r_inverse +
            y2
        );

    std::complex<double> dy2;
    std::complex<double> dy4;
    if constexpr (dynamic)
    {
        const double dynamic_term = -rs_args_ptr->frequency * rs_args_ptr->frequency * density * radius;

        dy2 =
            r_inverse * (
                y1 * (dynamic_term - 2.0 * density_gravity) +
                y2 * -2.0 +
                y4 * rs_args_ptr->llp1 +
                y5 * density * rs_args_ptr->lp1 +
                y6 * -density * radius +
                dy1 * 2.0 * lame +
                y1_y3_term * (2.0 * (lame + shear_modulus) * r_inverse - density_gravity)
            );

        dy4 =
            r_inverse * (
                y1 * (density_gravity + two_shear_r_inv) +
                y3 * (dynamic_term - two_shear_r_inv) +
                y4 * -3.0 +
                y5 * -density +
                dy1 * -lame +
                y1_y3_term * -lame_2mu * r_inverse
            );
    }
    else
    {
        dy2 =
            r_inverse * (
                y1 * -2.0 * density_gravity +
                y2 * -2.0 +
                y4 * rs_args_ptr->llp1 +
                y5 * density * rs_args_ptr->lp1 +
                y6 * -density * radius +
                dy1 * 2.0 * lame +
                y1_y3_term * (2.0 * (lame + shear_modulus) * r_inverse - density_gravity)
            );

        dy4 =
            r_inverse * (
                y1 * (density_gravity + two_shear_r_inv) +
                y3 * -two_shear_r_inv +
                y4 * -3.0 +
                y5 * -density +
                dy1 * -lame +
                y1_y3_term * -lame_2mu * r_inverse
            );
    }

    const std::complex<double> dy3 =
        y1 * -r_inverse +
        y3 * r_inverse +
        y4 * (1.0 / shear_modulus);

    const std::complex<double> dy5 =
        y1 * grav_term +
        y5 * -rs_args_ptr->lp1 * r_inverse +
        y6;

    const std::complex<double> dy6 =
        r_inverse * (
            y1 * grav_term * rs_args_ptr->lm1 +
            y6 * rs_args_ptr->lm1 +
            y1_y3_term * grav_term
        );

    c_write_dy6(dy_ptr, dy1, dy2, dy3, dy4, dy5, dy6);
}

/// Radial derivative equations for a solid, dynamic, compressible layer.
inline void c_solid_dynamic_compressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused
        ) noexcept
{
    c_solid_compressible_body<true>(dy_ptr, radius, y_ptr, args_ptr);
}

/// Radial derivative equations for a solid, static, compressible layer.
inline void c_solid_static_compressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused
        ) noexcept
{
    c_solid_compressible_body<false>(dy_ptr, radius, y_ptr, args_ptr);
}


// ============================================================================
//  Solid Incompressible (static and dynamic)
// ============================================================================

/// Radial derivative equations for a solid, incompressible layer. The static form drops the inertial term, which
/// appears only in dy2 and dy4.
///
/// References: TS72 Eq. 82, KMN15 Eqs. 4--9, B15 Eqs. 13--18
template <bool dynamic>
inline void c_solid_incompressible_body(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr
        ) noexcept
{
    c_RadialSolverArgs* rs_args_ptr = reinterpret_cast<c_RadialSolverArgs*>(args_ptr);

    double gravity, density;
    std::complex<double> shear_modulus, bulk_modulus;
    c_read_eos(rs_args_ptr, radius, gravity, density, shear_modulus, bulk_modulus);

    std::complex<double> y1, y2, y3, y4, y5, y6;
    c_read_y6(y_ptr, y1, y2, y3, y4, y5, y6);

    const double r_inverse       = 1.0 / radius;
    const double density_gravity = density * gravity;
    const double grav_term       = rs_args_ptr->grav_coeff * density;
    const std::complex<double> two_shear_r_inv = 2.0 * shear_modulus * r_inverse;
    const std::complex<double> y1_y3_term      = 2.0 * y1 - rs_args_ptr->llp1 * y3;

    const std::complex<double> dy1 =
        y1_y3_term * -1.0 * r_inverse;

    std::complex<double> dy2;
    std::complex<double> dy4;
    if constexpr (dynamic)
    {
        const double dynamic_term = -rs_args_ptr->frequency * rs_args_ptr->frequency * density * radius;

        dy2 =
            r_inverse * (
                y1 * (dynamic_term + 12.0 * shear_modulus * r_inverse - 4.0 * density_gravity) +
                y3 * rs_args_ptr->llp1 * (density_gravity - 6.0 * shear_modulus * r_inverse) +
                y4 * rs_args_ptr->llp1 +
                y5 * density * rs_args_ptr->lp1 +
                y6 * -density * radius
            );

        dy4 =
            r_inverse * (
                y1 * (density_gravity - 3.0 * two_shear_r_inv) +
                y2 * -1.0 +
                y3 * (dynamic_term + two_shear_r_inv * (2.0 * rs_args_ptr->llp1 - 1.0)) +
                y4 * -3.0 +
                y5 * -density
            );
    }
    else
    {
        dy2 =
            r_inverse * (
                y1 * (12.0 * shear_modulus * r_inverse - 4.0 * density_gravity) +
                y3 * rs_args_ptr->llp1 * (density_gravity - 6.0 * shear_modulus * r_inverse) +
                y4 * rs_args_ptr->llp1 +
                y5 * density * rs_args_ptr->lp1 +
                y6 * -density * radius
            );

        dy4 =
            r_inverse * (
                y1 * (density_gravity - 3.0 * two_shear_r_inv) +
                y2 * -1.0 +
                y3 * (two_shear_r_inv * (2.0 * rs_args_ptr->llp1 - 1.0)) +
                y4 * -3.0 +
                y5 * -density
            );
    }

    const std::complex<double> dy3 =
        y1 * -r_inverse +
        y3 * r_inverse +
        y4 * (1.0 / shear_modulus);

    const std::complex<double> dy5 =
        y1 * grav_term +
        y5 * -rs_args_ptr->lp1 * r_inverse +
        y6;

    const std::complex<double> dy6 =
        r_inverse * (
            y1 * grav_term * rs_args_ptr->lm1 +
            y6 * rs_args_ptr->lm1 +
            y1_y3_term * grav_term
        );

    c_write_dy6(dy_ptr, dy1, dy2, dy3, dy4, dy5, dy6);
}

/// Radial derivative equations for a solid, dynamic, incompressible layer.
inline void c_solid_dynamic_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused
        ) noexcept
{
    c_solid_incompressible_body<true>(dy_ptr, radius, y_ptr, args_ptr);
}

/// Radial derivative equations for a solid, static, incompressible layer.
inline void c_solid_static_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused
        ) noexcept
{
    c_solid_incompressible_body<false>(dy_ptr, radius, y_ptr, args_ptr);
}


// ============================================================================
//  Liquid Dynamic Compressible
// ============================================================================

/// Radial derivative equations for a liquid, dynamic, compressible layer.
/// y4 = 0 always. y3 is computed analytically from y1, y2, y5.
/// Active y-values: y1, y2, y5, y6 (stored as 4 values, 8 doubles).
///
/// References: TS72 Eq. 87, KMN15 Eqs. 11--14
inline void c_liquid_dynamic_compressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused
        ) noexcept
{
    c_RadialSolverArgs* rs_args_ptr = reinterpret_cast<c_RadialSolverArgs*>(args_ptr);

    double gravity, density;
    std::complex<double> shear_modulus, bulk_modulus;
    c_read_eos(rs_args_ptr, radius, gravity, density, shear_modulus, bulk_modulus);

    std::complex<double> y1, y2, y5, y6;
    c_read_y4_liquid(y_ptr, y1, y2, y5, y6);

    const double r_inverse         = 1.0 / radius;
    const double density_gravity   = density * gravity;
    const double f2                = rs_args_ptr->frequency * rs_args_ptr->frequency;
    const double dynamic_term_no_r = -f2 * density;
    const double grav_term         = rs_args_ptr->grav_coeff * density;

    // Liquid: shear = 0 so lame = bulk modulus.
    const std::complex<double> lame_inverse = 1.0 / bulk_modulus;

    // y1_y3_term formed without an explicit y3 for numerical stability.
    const std::complex<double> coeff_r = rs_args_ptr->llp1 / (f2 * radius);
    const std::complex<double> y1_y3_term =
        y1 * (2.0 - gravity * coeff_r) +
        y2 * coeff_r / density +
        y5 * coeff_r;

    const std::complex<double> dy1 =
        y2 * lame_inverse -
        y1_y3_term * r_inverse;

    // The solid form carries a 2 (lame + mu) / r coefficient on y1_y3_term in TS72; it vanishes for mu = 0.
    const std::complex<double> dy2 =
        y1 * (dynamic_term_no_r - 2.0 * density_gravity * r_inverse) +
        y5 * density * rs_args_ptr->lp1 * r_inverse -
        y6 * density -
        y1_y3_term * density_gravity * r_inverse;

    const std::complex<double> dy5 =
        y1 * grav_term -
        y5 * rs_args_ptr->lp1 * r_inverse +
        y6;

    const std::complex<double> dy6 =
        r_inverse * (
            rs_args_ptr->lm1 * (y1 * grav_term + y6) +
            y1_y3_term * grav_term
        );

    c_write_dy4(dy_ptr, dy1, dy2, dy5, dy6);
}


// ============================================================================
//  Liquid Dynamic Incompressible
// ============================================================================

/// Radial derivative equations for a liquid, dynamic, incompressible layer.
/// y4 = 0 always. y3 computed analytically. div(u) = 0.
/// Active y-values: y1, y2, y5, y6 (stored as 4 values, 8 doubles).
///
/// References: TS72 Eq. 87
inline void c_liquid_dynamic_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused
        ) noexcept
{
    c_RadialSolverArgs* rs_args_ptr = reinterpret_cast<c_RadialSolverArgs*>(args_ptr);

    double gravity, density;
    std::complex<double> shear_modulus, bulk_modulus;
    c_read_eos(rs_args_ptr, radius, gravity, density, shear_modulus, bulk_modulus);

    std::complex<double> y1, y2, y5, y6;
    c_read_y4_liquid(y_ptr, y1, y2, y5, y6);

    const double r_inverse       = 1.0 / radius;
    const double density_gravity = density * gravity;
    const double dynamic_term    = -rs_args_ptr->frequency * rs_args_ptr->frequency * density * radius;
    const double grav_term       = rs_args_ptr->grav_coeff * density;

    const std::complex<double> y3 =
        (1.0 / dynamic_term) * (y2 + density * y5 - density_gravity * y1);
    const std::complex<double> y1_y3_term = 2.0 * y1 - rs_args_ptr->llp1 * y3;

    const std::complex<double> dy1 =
        y1_y3_term * -r_inverse;

    const std::complex<double> dy2 =
        r_inverse * (
            y1 * (dynamic_term - 2.0 * density_gravity) +
            y5 * density * rs_args_ptr->lp1 +
            y6 * -density * radius +
            y1_y3_term * -density_gravity  // the solid form's lame term vanishes for mu = 0
        );

    const std::complex<double> dy5 =
        y1 * grav_term +
        y5 * -rs_args_ptr->lp1 * r_inverse +
        y6;

    const std::complex<double> dy6 =
        r_inverse * (
            y1 * grav_term * rs_args_ptr->lm1 +
            y6 * rs_args_ptr->lm1 +
            y1_y3_term * grav_term
        );

    c_write_dy4(dy_ptr, dy1, dy2, dy5, dy6);
}


// ============================================================================
//  Liquid Static Incompressible
// ============================================================================

/// Radial derivative equations for a liquid, static, incompressible layer.
/// Only y5 and y7 are defined. y7 = y6 + (4*pi*G/g)*y2.
/// Active y-values: y5, y7 (stored as 2 values, 4 doubles).
///
/// References: S74 Eq. 18
inline void c_liquid_static_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused
        ) noexcept
{
    c_RadialSolverArgs* rs_args_ptr = reinterpret_cast<c_RadialSolverArgs*>(args_ptr);

    double gravity, density;
    std::complex<double> shear_modulus, bulk_modulus;
    c_read_eos(rs_args_ptr, radius, gravity, density, shear_modulus, bulk_modulus);

    std::complex<double> y5, y7;
    c_read_y4_static_liquid(y_ptr, y5, y7);

    const double r_inverse = 1.0 / radius;
    const double grav_term = rs_args_ptr->grav_coeff * density / gravity;

    const std::complex<double> dy5 =
        y5 * (grav_term - rs_args_ptr->lp1 * r_inverse) +
        y7;

    const std::complex<double> dy7 =
        y5 * 2.0 * rs_args_ptr->lm1 * r_inverse * grav_term +
        y7 * (rs_args_ptr->lm1 * r_inverse - grav_term);

    c_write_dy2(dy_ptr, dy5, dy7);
}


// ============================================================================
//  Dispatch Functions
// ============================================================================

/// Return the correct ODE function for the given layer type and assumptions.
inline DiffeqFuncType c_find_layer_diffeq(
        int layer_type,
        int layer_is_static,
        int layer_is_incomp
        ) noexcept
{
    if (layer_type == 0)
    {
        if (layer_is_static == 1)
        {
            if (layer_is_incomp == 1)
            {
                return c_solid_static_incompressible;
            }
            else
            {
                return c_solid_static_compressible;
            }
        }
        else
        {
            if (layer_is_incomp == 1)
            {
                return c_solid_dynamic_incompressible;
            }
            else
            {
                return c_solid_dynamic_compressible;
            }
        }
    }
    else
    {
        if (layer_is_static == 1)
        {
            // TODO: A compressible static liquid uses the incompressible function; confirm this is correct.
            return c_liquid_static_incompressible;
        }
        else
        {
            if (layer_is_incomp == 1)
            {
                return c_liquid_dynamic_incompressible;
            }
            else
            {
                return c_liquid_dynamic_compressible;
            }
        }
    }
}


/// Number of independent shooting solutions: 3 solid, 2 dynamic liquid, 1 static liquid.
inline size_t c_find_num_shooting_solutions(
        int layer_type,
        int layer_is_static,
        int layer_is_incomp
        ) noexcept
{
    if (layer_type == 0)
    {
        return 3;
    }
    else
    {
        if (layer_is_static == 1)
        {
            return 1;
        }
        else
        {
            return 2;
        }
    }
}
