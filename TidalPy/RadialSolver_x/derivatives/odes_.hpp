#pragma once

#include <complex>
#include <cmath>

#include "c_common.hpp"        // CyRK: DiffeqFuncType, PreEvalFunc
#include "eos_solution_.hpp"   // TidalPy: c_EOSSolution
#include "../constants_.hpp"   // RadialSolver_x: C_MAX_NUM_Y, etc.


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
//  Helper: read EOS outputs at a given radius
// ============================================================================

/// Read EOS array called from a CyRK Dense output and extract gravity,
/// density, shear, and bulk modulus.
///
/// The EOS stores 4 doubles and two complex:
///   0: Gravity
///   1: Pressure
///   2: Mass
///   3: MOI
///   4: Density
///   5a: Complex Shear (real)
///   5b: Complex Shear (imag)
///   6a: Complex Bulk (real)
///   6b: Complex Bulk (imag)
static inline void c_read_eos(
        c_RadialSolverArgs* rs_args_ptr,
        double radius,
        double& gravity,
        double& density,
        std::complex<double>& shear_modulus,
        std::complex<double>& bulk_modulus
        ) noexcept
{
    double eos_array[9];
    rs_args_ptr->eos_solution_ptr->call(rs_args_ptr->layer_index, radius, &eos_array[0]);

    gravity       = eos_array[0];
    density       = eos_array[4];
    shear_modulus = std::complex<double>(eos_array[5], eos_array[6]);
    bulk_modulus  = std::complex<double>(eos_array[7], eos_array[8]);
}


/// Read y-values from the double-pair array (CyRK format) into 6 complex values.
/// Read y-values for solid (6 active: y1, y2, y3, y4, y5, y6); applies to static or dynamic; compressible or incomp.
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

/// Read y-values for dynamic liquid (4 active: y1, y2, y5, y6). Applies to compressible and incomp.
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

/// Read y-values for static liquid (2 active: y5, y7). Applies to incomp only.
static inline void c_read_y4_static_liquid(
        double* y_ptr,
        std::complex<double>& y5,
        std::complex<double>& y7
        ) noexcept
{
    y5 = std::complex<double>(y_ptr[0], y_ptr[1]);
    y7 = std::complex<double>(y_ptr[2], y_ptr[3]);
}

/// Write 6 complex dy values back to the double-pair output array.
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

/// Write 4 complex dy values back to the double-pair output array (dynamic liquid).
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

/// Write 2 complex dy values back to the double-pair output array (static liquid).
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
//  Solid Dynamic Compressible
// ============================================================================

/// Radial derivative equations for a solid, dynamic, compressible layer.
///
/// References: TS72 Eq. 82, KMN15 Eqs. 4--9, B15 Eqs. 13--18
inline void c_solid_dynamic_compressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused
        ) noexcept
{
    c_RadialSolverArgs* rs_args_ptr = reinterpret_cast<c_RadialSolverArgs*>(args_ptr);

    // EOS lookup
    double gravity, density;
    std::complex<double> shear_modulus, bulk_modulus;
    c_read_eos(rs_args_ptr, radius, gravity, density, shear_modulus, bulk_modulus);

    // y-values
    std::complex<double> y1, y2, y3, y4, y5, y6;
    c_read_y6(y_ptr, y1, y2, y3, y4, y5, y6);

    // Lame parameter
    const std::complex<double> lame = bulk_modulus - (2.0 / 3.0) * shear_modulus;

    // Precomputed terms
    const double r_inverse       = 1.0 / radius;
    const double density_gravity = density * gravity;
    const double dynamic_term    = -rs_args_ptr->frequency * rs_args_ptr->frequency * density * radius;
    const double grav_term       = rs_args_ptr->grav_coeff * density;

    const std::complex<double> lame_2mu         = lame + 2.0 * shear_modulus;
    const std::complex<double> lame_2mu_inverse = 1.0 / lame_2mu;
    const std::complex<double> two_shear_r_inv  = 2.0 * shear_modulus * r_inverse;
    const std::complex<double> y1_y3_term       = 2.0 * y1 - rs_args_ptr->llp1 * y3;

    // Derivatives (TS72 Eq. 82; KMN15 Eqs. 4-9; B15 Eqs. 13-18)
    const std::complex<double> dy1 =
        lame_2mu_inverse * (
            y1_y3_term * -lame * r_inverse +
            y2
        );

    const std::complex<double> dy2 =
        r_inverse * (
            y1 * (dynamic_term - 2.0 * density_gravity) +
            y2 * -2.0 +
            y4 * rs_args_ptr->llp1 +
            y5 * density * rs_args_ptr->lp1 +
            y6 * -density * radius +
            dy1 * 2.0 * lame +
            y1_y3_term * (2.0 * (lame + shear_modulus) * r_inverse - density_gravity)
        );

    const std::complex<double> dy3 =
        y1 * -r_inverse +
        y3 * r_inverse +
        y4 * (1.0 / shear_modulus);

    const std::complex<double> dy4 =
        r_inverse * (
            y1 * (density_gravity + two_shear_r_inv) +
            y3 * (dynamic_term - two_shear_r_inv) +
            y4 * -3.0 +
            y5 * -density +
            dy1 * -lame +
            y1_y3_term * -lame_2mu * r_inverse
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

    c_write_dy6(dy_ptr, dy1, dy2, dy3, dy4, dy5, dy6);
}


// ============================================================================
//  Solid Dynamic Incompressible
// ============================================================================

/// Radial derivative equations for a solid, dynamic, incompressible layer.
///
/// References: TS72 Eq. 82, KMN15 Eqs. 4--9, B15 Eqs. 13--18
inline void c_solid_dynamic_incompressible(
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

    std::complex<double> y1, y2, y3, y4, y5, y6;
    c_read_y6(y_ptr, y1, y2, y3, y4, y5, y6);

    const double r_inverse               = 1.0 / radius;
    const double density_gravity         = density * gravity;
    const double dynamic_term            = -rs_args_ptr->frequency * rs_args_ptr->frequency * density * radius;
    const double grav_term               = rs_args_ptr->grav_coeff * density;
    const std::complex<double> two_shear_r_inv = 2.0 * shear_modulus * r_inverse;
    const std::complex<double> y1_y3_term      = 2.0 * y1 - rs_args_ptr->llp1 * y3;

    const std::complex<double> dy1 =
        y1_y3_term * -1.0 * r_inverse;

    const std::complex<double> dy2 =
        r_inverse * (
            y1 * (dynamic_term + 12.0 * shear_modulus * r_inverse - 4.0 * density_gravity) +
            y3 * rs_args_ptr->llp1 * (density_gravity - 6.0 * shear_modulus * r_inverse) +
            y4 * rs_args_ptr->llp1 +
            y5 * density * rs_args_ptr->lp1 +
            y6 * -density * radius
        );

    const std::complex<double> dy3 =
        y1 * -r_inverse +
        y3 * r_inverse +
        y4 * (1.0 / shear_modulus);

    const std::complex<double> dy4 =
        r_inverse * (
            y1 * (density_gravity - 3.0 * two_shear_r_inv) +
            y2 * -1.0 +
            y3 * (dynamic_term + two_shear_r_inv * (2.0 * rs_args_ptr->llp1 - 1.0)) +
            y4 * -3.0 +
            y5 * -density
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

    c_write_dy6(dy_ptr, dy1, dy2, dy3, dy4, dy5, dy6);
}


// ============================================================================
//  Solid Static Compressible
// ============================================================================

/// Radial derivative equations for a solid, static, compressible layer.
/// The static case sets all frequency dependence to zero.
///
/// References: TS72 Eq. 82, KMN15 Eqs. 4--9, B15 Eqs. 13--18
inline void c_solid_static_compressible(
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

    const std::complex<double> dy2 =
        r_inverse * (
            y1 * -2.0 * density_gravity +
            y2 * -2.0 +
            y4 * rs_args_ptr->llp1 +
            y5 * density * rs_args_ptr->lp1 +
            y6 * -density * radius +
            dy1 * 2.0 * lame +
            y1_y3_term * (2.0 * (lame + shear_modulus) * r_inverse - density_gravity)
        );

    const std::complex<double> dy3 =
        y1 * -r_inverse +
        y3 * r_inverse +
        y4 * (1.0 / shear_modulus);

    const std::complex<double> dy4 =
        r_inverse * (
            y1 * (density_gravity + two_shear_r_inv) +
            y3 * -two_shear_r_inv +
            y4 * -3.0 +
            y5 * -density +
            dy1 * -lame +
            y1_y3_term * -lame_2mu * r_inverse
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

    c_write_dy6(dy_ptr, dy1, dy2, dy3, dy4, dy5, dy6);
}


// ============================================================================
//  Solid Static Incompressible
// ============================================================================

/// Radial derivative equations for a solid, static, incompressible layer.
///
/// References: TS72 Eq. 82, KMN15 Eqs. 4--9, B15 Eqs. 13--18
inline void c_solid_static_incompressible(
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

    std::complex<double> y1, y2, y3, y4, y5, y6;
    c_read_y6(y_ptr, y1, y2, y3, y4, y5, y6);

    const double r_inverse       = 1.0 / radius;
    const double density_gravity = density * gravity;
    const double grav_term       = rs_args_ptr->grav_coeff * density;
    const std::complex<double> two_shear_r_inv = 2.0 * shear_modulus * r_inverse;
    const std::complex<double> y1_y3_term      = 2.0 * y1 - rs_args_ptr->llp1 * y3;

    const std::complex<double> dy1 =
        -1.0 * y1_y3_term * r_inverse;

    const std::complex<double> dy2 =
        r_inverse * (
            y1 * (12.0 * shear_modulus * r_inverse - 4.0 * density_gravity) +
            y3 * rs_args_ptr->llp1 * (density_gravity - 6.0 * shear_modulus * r_inverse) +
            y4 * rs_args_ptr->llp1 +
            y5 * density * rs_args_ptr->lp1 +
            y6 * -density * radius
        );

    const std::complex<double> dy3 =
        y1 * -r_inverse +
        y3 * r_inverse +
        y4 * (1.0 / shear_modulus);

    const std::complex<double> dy4 =
        r_inverse * (
            y1 * (density_gravity - 3.0 * two_shear_r_inv) +
            y2 * -1.0 +
            y3 * (two_shear_r_inv * (2.0 * rs_args_ptr->llp1 - 1.0)) +
            y4 * -3.0 +
            y5 * -density
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

    c_write_dy6(dy_ptr, dy1, dy2, dy3, dy4, dy5, dy6);
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

    // 4 active y-values for dynamic liquid
    std::complex<double> y1, y2, y5, y6;
    c_read_y4_liquid(y_ptr, y1, y2, y5, y6);

    const double r_inverse         = 1.0 / radius;
    const double density_gravity   = density * gravity;
    const double f2                = rs_args_ptr->frequency * rs_args_ptr->frequency;
    const double dynamic_term_no_r = -f2 * density;
    const double dynamic_term      = dynamic_term_no_r * radius;
    const double grav_term         = rs_args_ptr->grav_coeff * density;

    // For liquid layers, shear = 0 so lame = bulk modulus
    const std::complex<double> lame_inverse = 1.0 / bulk_modulus;

    // Compute y1_y3_term analytically (avoiding explicit y3 computation for numerical stability)
    const std::complex<double> coeff_r = rs_args_ptr->llp1 / (f2 * radius);
    const std::complex<double> y1_y3_term =
        y1 * (2.0 - gravity * coeff_r) +
        y2 * coeff_r / density +
        y5 * coeff_r;

    // Derivatives (TS72 Eq. 87)
    const std::complex<double> dy1 =
        y2 * lame_inverse -
        y1_y3_term * r_inverse;

    // FIX: In TS72 the solid version has a [2*(lame+mu)*r_inv] coefficient for y1_y3_term;
    //      for liquid (mu=0) the first term vanishes. The original code has a TODO about this.
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

    // y3 for incompressible dynamic liquid
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
            // FIX: Same TODO as compressible liquid re: missing lame term from solid version
            y1_y3_term * -density_gravity
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
/// Only y5 and y7 are defined. y7 = y6 + (4*pi*G*rho/g)*y5.
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

    // Only y5 and y7 for static liquid
    std::complex<double> y5, y7;
    c_read_y4_static_liquid(y_ptr, y5, y7);

    const double r_inverse = 1.0 / radius;
    const double grav_term = rs_args_ptr->grav_coeff * density / gravity;

    // S74 Eq. 18
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
        // Solid
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
        // Liquid
        if (layer_is_static == 1)
        {
            // FIX: Compressible static liquid uses same function as incompressible. Check if this is correct.
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


/// Return the number of independent shooting solutions for the given layer type.
///
/// Solid layers: 3 solutions
/// Dynamic liquid layers: 2 solutions (compressible or incompressible)
/// Static liquid layers: 1 solution
inline size_t c_find_num_shooting_solutions(
        int layer_type,
        int layer_is_static,
        int layer_is_incomp
        ) noexcept
{
    if (layer_type == 0)
    {
        // Solid: always 3 independent solutions
        return 3;
    }
    else
    {
        // Liquid
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
