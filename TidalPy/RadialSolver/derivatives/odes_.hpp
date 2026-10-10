#pragma once
/* The radial ODEs of the shooting method. A layer's independent solutions are integrated as one system: each
 * right-hand-side call reads the material once (c_read_eos) and applies the layer kind's equations to every solution
 * in turn, stored one after another, two reals per complex y (layer_kind_.hpp gives the counts).
 *
 * A dynamic liquid carries the pressure-like variable P = y2 - rho g y1 + rho y5 in place of y2 (its slot), so its ys
 * are (y1, P, y5, y6). P is minus the Eulerian pressure perturbation plus rho times the potential perturbation, and
 * y3 = -P / (rho omega^2 r). TS72's y2 form couples y1 and y2 through l(l+1) / (omega^2 r), which pairs eigenvalues
 * about +-sqrt(l(l+1) rho g^2 / K) / (omega r) and shrinks an explicit integrator's steps with the forcing period. In P
 * the large coupling enters through P alone, and the remaining rates are the buoyancy ones,
 * sqrt(-l(l+1) N^2) / (omega r) with N^2 = -g (rho' / rho + rho g / K). The change of variables is exact given
 * rho' = d rho / dr of the same density the equations read (from the EOS layout's N^2) and the structure's
 * dg/dr = 4 pi G rho - 2 g / r. P is stored in the unit c_pressure_unit, omega^2 / (pi G rho_bar) of a stress at long
 * periods, where it is of order rho r y3. The shooting method converts y2 to P where a liquid layer starts and back at
 * its top, and the dense solution converts back too (c_y2_from_pressure_variable), so every y the solver reports is
 * TS72's.
 *
 * References: TS72 (Takeuchi and Saito 1972), KMN15 (Kamata et al. 2015), B15 (Beuthe 2015), S74 (Saito 1974)
 */

#include <array>
#include <complex>
#include <cmath>

#include "c_common.hpp"        // CyRK: DiffeqFuncType, PreEvalFunc
#include "eos_solution_.hpp"   // TidalPy: c_EOSSolution
#include "../../constants_.hpp"
#include "../layer_kind_.hpp"   // c_LayerKindLayout, c_layer_num_solutions
#include "solid_dy1_.hpp"      // c_solid_dy1_compressible, c_solid_dy1_incompressible


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
    double pressure_unit;   // a dynamic liquid's stored P per stress (c_pressure_unit)
    size_t layer_index;
    c_EOSSolution* eos_solution_ptr;
    // The layer's solutions and the complex ys each carries, each stored y's weight (c_layer_y_weights),
    // and the normalized Gram determinant at which the current integration segment ends so the shooting method can
    // re-orthonormalize them (c_solution_independence_event).
    size_t num_solutions;
    size_t num_ys;
    std::array<double, 6> y_weights;
    double independence_floor;
};


/// The material the radial ODEs read at one radius, once per right-hand-side call.
struct c_RadialMaterial
{
    double gravity;
    double density;
    // The EOS's buoyancy frequency squared, the static bulk modulus K_S it is measured against, and d rho / dr; read
    // by the dynamic liquids.
    double buoyancy_frequency_squared;
    double static_bulk_modulus;
    double density_gradient;
    std::complex<double> shear_modulus;
    std::complex<double> bulk_modulus;
};


// ============================================================================
//  Helpers: EOS lookup and y-vector packing
// ============================================================================

/// Gravity, density, N^2, the static bulk modulus, rho', and the complex moduli at a radius, in one EOS evaluation. The
/// complex moduli are the material's response at this solve's forcing frequency; see c_EOSSolution::call_material.
static inline c_RadialMaterial c_read_eos(const c_RadialSolverArgs* rs_args_ptr, double radius) noexcept
{
    c_EOSMaterialState material_state;
    rs_args_ptr->eos_solution_ptr->call_material(rs_args_ptr->layer_index, radius, material_state);
    return c_RadialMaterial{
        material_state.gravity, material_state.density, material_state.buoyancy_frequency_squared,
        material_state.static_bulk_modulus, material_state.density_gradient, material_state.shear_modulus,
        material_state.bulk_modulus};
}


/// Read the six solid-layer ys from CyRK's real-pair array.
static inline void c_read_y6(
        const double* y_ptr,
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

/// Read the four dynamic-liquid ys (y1, P, y5, y6).
static inline void c_read_y4_liquid(
        const double* y_ptr,
        std::complex<double>& y1,
        std::complex<double>& pressure_variable,
        std::complex<double>& y5,
        std::complex<double>& y6
        ) noexcept
{
    y1                = std::complex<double>(y_ptr[0], y_ptr[1]);
    pressure_variable = std::complex<double>(y_ptr[2], y_ptr[3]);
    y5                = std::complex<double>(y_ptr[4], y_ptr[5]);
    y6                = std::complex<double>(y_ptr[6], y_ptr[7]);
}

/// Read the two static-liquid ys (y5, y7).
static inline void c_read_y4_static_liquid(
        const double* y_ptr,
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
//  Per-solution equations
// ============================================================================

/// Radial derivative equations for one solution of a solid, compressible layer. The static form drops the inertial
/// term, which appears only in dy2 and dy4.
///
/// References: TS72 Eq. 82, KMN15 Eqs. 4--9, B15 Eqs. 13--18
template <bool dynamic>
inline void c_solid_compressible_kernel(
        const c_RadialSolverArgs* rs_args_ptr,
        const c_RadialMaterial& material,
        double radius,
        const double* y_ptr,
        double* dy_ptr
        ) noexcept
{
    const double gravity = material.gravity;
    const double density = material.density;
    const std::complex<double>& shear_modulus = material.shear_modulus;
    const std::complex<double>& bulk_modulus  = material.bulk_modulus;

    std::complex<double> y1, y2, y3, y4, y5, y6;
    c_read_y6(y_ptr, y1, y2, y3, y4, y5, y6);

    const std::complex<double> lame = bulk_modulus - (2.0 / 3.0) * shear_modulus;

    const double r_inverse       = 1.0 / radius;
    const double density_gravity = density * gravity;
    const double grav_term       = rs_args_ptr->grav_coeff * density;

    const std::complex<double> lame_2mu         = lame + 2.0 * shear_modulus;
    const std::complex<double> two_shear_r_inv  = 2.0 * shear_modulus * r_inverse;
    const std::complex<double> y1_y3_term       = 2.0 * y1 - rs_args_ptr->llp1 * y3;

    const std::complex<double> dy1 = c_solid_dy1_compressible(y1_y3_term, y2, lame, lame_2mu, r_inverse);

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


/// Radial derivative equations for one solution of a solid, incompressible layer. The static form drops the inertial
/// term, which appears only in dy2 and dy4.
///
/// References: TS72 Eq. 82, KMN15 Eqs. 4--9, B15 Eqs. 13--18
template <bool dynamic>
inline void c_solid_incompressible_kernel(
        const c_RadialSolverArgs* rs_args_ptr,
        const c_RadialMaterial& material,
        double radius,
        const double* y_ptr,
        double* dy_ptr
        ) noexcept
{
    const double gravity = material.gravity;
    const double density = material.density;
    const std::complex<double>& shear_modulus = material.shear_modulus;

    std::complex<double> y1, y2, y3, y4, y5, y6;
    c_read_y6(y_ptr, y1, y2, y3, y4, y5, y6);

    const double r_inverse       = 1.0 / radius;
    const double density_gravity = density * gravity;
    const double grav_term       = rs_args_ptr->grav_coeff * density;
    const std::complex<double> two_shear_r_inv = 2.0 * shear_modulus * r_inverse;
    const std::complex<double> y1_y3_term      = 2.0 * y1 - rs_args_ptr->llp1 * y3;

    const std::complex<double> dy1 = c_solid_dy1_incompressible(y1_y3_term, r_inverse);

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


/// A dynamic liquid's stratification rho' + rho^2 g / K = -rho N^2 / g in its own bulk modulus K, from the EOS's N_S^2
/// against its static K_S, with rho' = -rho N_S^2 / g - rho^2 g / K_S: -rho N_S^2 / g + rho^2 g (1 - K / K_S) / K. A
/// neutral liquid's is then exactly 0 rather than the roundoff of two large terms, which would set its response at very
/// long periods. K differs from K_S where a rheology makes it complex or the moduli are supplied, and an EOS with no
/// bulk modulus (K_S = 0) has no K_S term. An incompressible liquid (1 / K = 0) reads rho' itself, which a constant
/// density gives as exactly 0.
template <bool incompressible>
inline std::complex<double> c_liquid_stratification(const c_RadialMaterial& material) noexcept
{
    if constexpr (incompressible)
    {
        return material.density_gradient;
    }
    else
    {
        const bool   has_static_bulk       = material.static_bulk_modulus > 0.0;
        const double static_stratification =
            -material.density * material.buoyancy_frequency_squared / material.gravity;
        const double weight_term           = material.density * material.density * material.gravity;
        const std::complex<double> off_static_bulk =
            has_static_bulk ? 1.0 - material.bulk_modulus / material.static_bulk_modulus : std::complex<double>(1.0);
        return static_stratification + weight_term * off_static_bulk / material.bulk_modulus;
    }
}

/// Radial derivative equations for one solution of a dynamic liquid layer in (y1, P, y5, y6), P = y2 - rho g y1 +
/// rho y5, stored in the pressure unit (c_pressure_unit). y4 = 0, and y3 = -P / (rho omega^2 r) enters through
/// 2 y1 - l(l+1) y3 = 2 y1 + l(l+1) P / (rho omega^2 r). From TS72's y2 form (TS72 Eq. 87, KMN15 Eqs. 11--14),
/// dP/dr = dy2/dr - d(rho g)/dr y1 - rho g dy1/dr + rho' y5 + rho dy5/dr, in which dg/dr = 4 pi G rho - 2 g / r cancels
/// every term in 1 / r and in G. With y2 = P + rho g y1 - rho y5 this leaves
///     dP/dr = -omega^2 rho y1 + (rho' + rho^2 g / K) (y5 - g y1) - (rho g / K) P,
/// with the stratification rho' + rho^2 g / K from c_liquid_stratification. An incompressible liquid is the limit
/// 1 / K = 0 (div u = 0).
template <bool incompressible>
inline void c_liquid_dynamic_kernel(
        const c_RadialSolverArgs* rs_args_ptr,
        const c_RadialMaterial& material,
        double radius,
        const double* y_ptr,
        double* dy_ptr
        ) noexcept
{
    const double gravity = material.gravity;
    const double density = material.density;

    std::complex<double> y1, pressure_variable, y5, y6;
    c_read_y4_liquid(y_ptr, y1, pressure_variable, y5, y6);

    const double r_inverse = 1.0 / radius;
    const double f2        = rs_args_ptr->frequency * rs_args_ptr->frequency;
    const double grav_term = rs_args_ptr->grav_coeff * density;
    // P is stored in the pressure unit, and so are its rates: omega^2 and the stratification enter divided by it.
    const double pressure_unit = rs_args_ptr->pressure_unit;
    const double inertia       = f2 / pressure_unit;

    // 2 y1 - l(l+1) y3, with y3 = -P / (rho omega^2 r).
    const double coeff_r = rs_args_ptr->llp1 / (inertia * radius);
    const std::complex<double> y1_y3_term = 2.0 * y1 + (coeff_r / density) * pressure_variable;

    const std::complex<double> dy5 =
        y1 * grav_term -
        y5 * rs_args_ptr->lp1 * r_inverse +
        y6;

    const std::complex<double> dy6 =
        r_inverse * (
            rs_args_ptr->lm1 * (y1 * grav_term + y6) +
            y1_y3_term * grav_term
        );

    // The buoyancy term -(rho N^2 / g) (y5 - g y1) and the inertia -omega^2 rho y1.
    const std::complex<double> stratification = c_liquid_stratification<incompressible>(material);
    std::complex<double> dy1;
    std::complex<double> dpressure;
    if constexpr (incompressible)
    {
        dy1       = -y1_y3_term * r_inverse;
        dpressure = -inertia * density * y1 + (stratification / pressure_unit) * (y5 - gravity * y1);
    }
    else
    {
        // The liquid's shear modulus is zero, so the Lame parameter is K and the volume change is y2 / K.
        const std::complex<double> inverse_bulk = 1.0 / material.bulk_modulus;
        const std::complex<double> y2 =
            c_y2_from_pressure_variable(y1, pressure_variable, y5, gravity, density, pressure_unit);
        dy1       = y2 * inverse_bulk - y1_y3_term * r_inverse;
        dpressure = -inertia * density * y1 + (stratification / pressure_unit) * (y5 - gravity * y1)
            - density * gravity * inverse_bulk * pressure_variable;
    }

    c_write_dy4(dy_ptr, dy1, dpressure, dy5, dy6);
}


/// Radial derivative equations for the one solution of a static liquid layer. Only y5 and y7 are defined.
/// y7 = y6 + (4*pi*G/g)*y2. Compressibility does not enter.
///
/// References: S74 Eq. 18
inline void c_liquid_static_kernel(
        const c_RadialSolverArgs* rs_args_ptr,
        const c_RadialMaterial& material,
        double radius,
        const double* y_ptr,
        double* dy_ptr
        ) noexcept
{
    std::complex<double> y5, y7;
    c_read_y4_static_liquid(y_ptr, y5, y7);

    const double r_inverse = 1.0 / radius;
    const double grav_term = rs_args_ptr->grav_coeff * material.density / material.gravity;

    const std::complex<double> dy5 =
        y5 * (grav_term - rs_args_ptr->lp1 * r_inverse) +
        y7;

    const std::complex<double> dy7 =
        y5 * 2.0 * rs_args_ptr->lm1 * r_inverse * grav_term +
        y7 * (rs_args_ptr->lm1 * r_inverse - grav_term);

    c_write_dy2(dy_ptr, dy5, dy7);
}


// ============================================================================
//  Layer systems: every solution of a layer, one material read per call
// ============================================================================

using c_RadialKernel = void (*)(const c_RadialSolverArgs*, const c_RadialMaterial&, double, const double*, double*);

/// A layer kind's solutions, layout.num_solutions of layout.num_ys complex ys each, stored in turn, from one material
/// read.
template <const c_LayerKindLayout& layout, c_RadialKernel kernel>
inline void c_layer_system(double* dy_ptr, double radius, double* y_ptr, char* args_ptr) noexcept
{
    constexpr size_t num_solutions = layout.num_solutions;
    constexpr size_t num_ys        = layout.num_ys;
    const c_RadialSolverArgs* rs_args_ptr = reinterpret_cast<const c_RadialSolverArgs*>(args_ptr);
    const c_RadialMaterial material = c_read_eos(rs_args_ptr, radius);
    for (size_t solution_i = 0; solution_i < num_solutions; ++solution_i)
    {
        kernel(rs_args_ptr, material, radius, y_ptr + 2 * num_ys * solution_i, dy_ptr + 2 * num_ys * solution_i);
    }
}

/// Solid, dynamic, compressible layer: three solutions of y1 to y6.
inline void c_solid_dynamic_compressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused) noexcept
{
    c_layer_system<C_SOLID_LAYOUT, c_solid_compressible_kernel<true>>(dy_ptr, radius, y_ptr, args_ptr);
}

/// Solid, static, compressible layer: three solutions of y1 to y6.
inline void c_solid_static_compressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused) noexcept
{
    c_layer_system<C_SOLID_LAYOUT, c_solid_compressible_kernel<false>>(dy_ptr, radius, y_ptr, args_ptr);
}

/// Solid, dynamic, incompressible layer: three solutions of y1 to y6.
inline void c_solid_dynamic_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused) noexcept
{
    c_layer_system<C_SOLID_LAYOUT, c_solid_incompressible_kernel<true>>(dy_ptr, radius, y_ptr, args_ptr);
}

/// Solid, static, incompressible layer: three solutions of y1 to y6.
inline void c_solid_static_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused) noexcept
{
    c_layer_system<C_SOLID_LAYOUT, c_solid_incompressible_kernel<false>>(dy_ptr, radius, y_ptr, args_ptr);
}

/// Liquid, dynamic, compressible layer: two solutions of (y1, P, y5, y6).
inline void c_liquid_dynamic_compressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused) noexcept
{
    c_layer_system<C_DYNAMIC_LIQUID_LAYOUT, c_liquid_dynamic_kernel<false>>(dy_ptr, radius, y_ptr, args_ptr);
}

/// Liquid, dynamic, incompressible layer: two solutions of (y1, P, y5, y6).
inline void c_liquid_dynamic_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused) noexcept
{
    c_layer_system<C_DYNAMIC_LIQUID_LAYOUT, c_liquid_dynamic_kernel<true>>(dy_ptr, radius, y_ptr, args_ptr);
}

/// Liquid, static layer, compressible or not: one solution of (y5, y7).
inline void c_liquid_static_incompressible(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* args_ptr,
        PreEvalFunc unused) noexcept
{
    c_layer_system<C_STATIC_LIQUID_LAYOUT, c_liquid_static_kernel>(dy_ptr, radius, y_ptr, args_ptr);
}


// ============================================================================
//  Dispatch Functions
// ============================================================================

/// Return the ODE system of a layer kind (every solution of the layer together).
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
            // A static liquid's equations (S74 Eq. 18) do not read the bulk modulus.
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
    return c_layer_num_solutions(layer_type, layer_is_static == 1);
}
