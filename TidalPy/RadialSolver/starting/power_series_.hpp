// power_series_.hpp: starting conditions from a power series about the center, after Martens (2016) and Smylie (2013).
//
// Near the center the layer is taken as homogeneous (constant density and moduli, g = gamma r with
// gamma = 4 pi G rho / 3), the same assumption as the Takeuchi and Saito (1972) and Kamata et al. (2015) closed forms.
// Scaling each y_i as y_i = r^(l - 2 + s_i + nu0) z_i turns every layer's equations (derivatives/odes_.hpp) into
// r dz/dr = (B0 + B2 r^2) z. A regular solution is z = sum_k c_k r^(2k), with (nu0 I - B0) c_0 = 0 and
// ((nu0 + 2k) I - B0) c_k = B2 c_(k-1); terms are added until the next one falls below machine precision, so the
// series is summed to convergence rather than cut at Martens' three terms. Solids have three regular solutions
// (nu0 = 0, 0, 2) and dynamic liquids two (nu0 = 0, 0); the dynamic incompressible liquid's series ends after its
// first term (B2 = 0). Static liquids need no series: its first term is Saito's (1974) solution, used for every
// method (driver_.hpp).
//
// The solid leading vectors come from Martens' free-constant solutions (A1,1, A6,1, A4,0; thesis Eqs. 4.96-4.104),
// converted to TidalPy's y6 (TS72) with y6_TS = y6_Martens + (l + 1) y5 / r; y1 to y5 are the same in both. Rotation
// is left out (Martens' Omega = 0). The free constants are arbitrary, and two choices differ from Martens':
//   - The first solution is TS72's polynomial solution (l A1,1 + (l gamma - w^2 - 3 gamma) A6,1), exact as one term.
//   - In a compressible solid with gravity, the other two are the series of TS72's two wavenumber solutions: each
//     starts as (alpha^2 f - (l + 1) beta^2) times the A6,1 vector, with the free part of its resonant k = 1 step
//     set to TS72's z4 = mu (1 - ((l - 1) h + 2 (f + 1)) / (2l + 3)). Each then carries one wavenumber, so the fast
//     wave of a weak solid does not leave two solutions growing along it together.
// The k = 1 step of the nu0 = 0 solid solutions is resonant with the nu0 = 2 solution (A4,0); otherwise its free part
// is set to zero, as Martens sets A4,2 = 0.
//
// References
// ----------
// Martens (2016) PhD thesis, Caltech, Sec. 4.2.8; Martens et al. (2019) LoadDef; Smylie (2013) Earth Dynamics;
// Crossley (1975) GJRAS 41, 153; TS72: Takeuchi & Saito (1972).
#pragma once

#include <cmath>
#include <complex>
#include <cstddef>

#include <Eigen/Dense>

#include "../../constants_.hpp"
#include "common_.hpp"   // c_solid_wave_numbers

// The series stops when a term falls below machine precision relative to the sum. It fails rather than return an
// inaccurate start when that needs more than C_POWER_SERIES_MAX_TERMS terms, when any component's terms grew so far
// beyond its sum that cancellation took more than half its digits, or when the solution grew from its first terms by
// more than 1 / sqrt(machine epsilon): the start then lies deep in a layer's exponential regime (|k r| above about 20;
// a weak solid, or a dynamic liquid at long periods), where the closed-form starts, which carry that growth
// analytically, apply.
constexpr size_t C_POWER_SERIES_MAX_TERMS = 100;

template <int N>
using c_SeriesMatrix = Eigen::Matrix<std::complex<double>, N, N>;
template <int N>
using c_SeriesVector = Eigen::Matrix<std::complex<double>, N, 1>;


// Sum one regular solution's series at radius: z = sum_k d_k with d_0 = lead and
// ((nu0 + 2k) I - B0) d_k = r^2 B2 d_(k-1). resonant_ptr, when given, is the nu0 = 2 solution's leading vector
// (z4 = 1): the singular nu = 2 step of a nu0 = 0 solid solution sets its component along it so that the step's z4
// coefficient is resonant_z4 (in the units of z4, before the r^2 of the term). scale holds each z component's typical
// size (a modulus for the stresses, gamma for the potential terms): the sum runs in z / scale, where the matrices are
// free of the moduli's units, so the linear solves and the convergence test hold in SI as well as non-dimensional
// units. False when the series refused (see C_POWER_SERIES_MAX_TERMS).
template <int N>
inline bool c_power_series_sum(
        const c_SeriesMatrix<N>& b0,
        const c_SeriesMatrix<N>& b2,
        const c_SeriesVector<N>& lead,
        int nu0,
        const c_SeriesVector<N>* resonant_ptr,
        std::complex<double> resonant_z4,
        const double (&scale)[N],
        double radius,
        c_SeriesVector<N>& z_out) noexcept
{
    const double r2 = radius * radius;
    const c_SeriesMatrix<N> identity = c_SeriesMatrix<N>::Identity();

    // S^-1 B S for S = diag(scale), and S^-1 on the vectors.
    c_SeriesMatrix<N> b0_scaled;
    c_SeriesMatrix<N> b2_scaled;
    c_SeriesVector<N> lead_scaled;
    c_SeriesVector<N> resonant_scaled;
    for (int row = 0; row < N; ++row)
    {
        for (int col = 0; col < N; ++col)
        {
            b0_scaled(row, col) = b0(row, col) * (scale[col] / scale[row]);
            b2_scaled(row, col) = b2(row, col) * (scale[col] / scale[row]);
        }
        lead_scaled(row) = lead(row) / scale[row];
        if (resonant_ptr != nullptr) { resonant_scaled(row) = (*resonant_ptr)(row) / scale[row]; }
    }

    c_SeriesVector<N> term = lead_scaled;
    c_SeriesVector<N> sum = lead_scaled;
    // The largest term of each component, which the cancellation test measures that component's sum against: a
    // component small next to the others (the stresses of a very weak solid) can lose its digits unseen in a norm.
    Eigen::Matrix<double, N, 1> largest_terms = term.cwiseAbs();
    // The size of the first terms, which the growth test measures the sum against.
    double first_terms = term.cwiseAbs().maxCoeff();
    for (size_t k = 1; k < C_POWER_SERIES_MAX_TERMS; ++k)
    {
        const double nu = static_cast<double>(nu0 + 2 * static_cast<int>(k));
        const c_SeriesMatrix<N> matrix = nu * identity - b0_scaled;
        const c_SeriesVector<N> rhs = r2 * (b2_scaled * term);
        if ((resonant_ptr != nullptr) && (nu == 2.0))
        {
            term = matrix.completeOrthogonalDecomposition().solve(rhs);
            const std::complex<double> target_z4 = resonant_z4 * r2 / scale[3];
            term += ((target_z4 - term(3)) / resonant_scaled(3)) * resonant_scaled;
        }
        else
        {
            term = matrix.partialPivLu().solve(rhs);
        }
        sum += term;
        const double term_size = term.cwiseAbs().maxCoeff();
        largest_terms = largest_terms.cwiseMax(term.cwiseAbs());
        if (k == 1) { first_terms = std::fmax(first_terms, term_size); }
        const double sum_size = sum.cwiseAbs().maxCoeff();
        if (!std::isfinite(sum_size))
        {
            return false;
        }
        if (term_size <= TidalPyConstants::d_EPS * sum_size)
        {
            for (int row = 0; row < N; ++row) { z_out(row) = sum(row) * scale[row]; }
            const double half_digits = std::sqrt(TidalPyConstants::d_EPS);
            // Refuse growth past 1 / sqrt(eps), or more than half the digits of any component lost to cancellation
            // (a component whose every term is zero is zero by construction).
            if (sum_size * half_digits > first_terms) { return false; }
            for (int row = 0; row < N; ++row)
            {
                if (largest_terms(row) * half_digits > std::abs(sum(row))) { return false; }
            }
            return true;
        }
    }
    return false;
}


// The size of the potential terms z5, z6 (and, times density, of a liquid's z2): gamma + omega^2 [s-2], or 1 for a
// static layer without gravity.
inline double c_power_series_gravity_scale(double gamma, double frequency2) noexcept
{
    const double gravity_scale = gamma + frequency2;
    return (gravity_scale > 0.0) ? gravity_scale : 1.0;
}


// Write each solution's y_i = r^(l - 2 + s_i + nu0) z_i into starting_conditions_ptr (num_ys per solution).
template <int N>
inline void c_power_series_write(
        const c_SeriesVector<N>& z,
        const double (&powers)[N],
        int nu0,
        double radius,
        int degree_l,
        size_t solution_i,
        size_t num_ys,
        std::complex<double>* starting_conditions_ptr) noexcept
{
    const double base_power = static_cast<double>(degree_l) - 2.0 + static_cast<double>(nu0);
    for (int y_i = 0; y_i < N; ++y_i)
    {
        starting_conditions_ptr[solution_i * num_ys + y_i] = z(y_i) * std::pow(radius, base_power + powers[y_i]);
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Solid Layers
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Starting conditions at the bottom of a solid layer: three independent solutions of y1..y6. frequency = 0 gives the
// static form (the inertial terms appear only in B2). Units: density [kg m-3], moduli [Pa], frequency [rad s-1].
// False when a series did not converge.
inline bool c_power_series_solid(
        const double frequency,
        const double radius,
        const double density,
        const std::complex<double>& bulk_modulus,
        const std::complex<double>& shear_modulus,
        const bool is_incompressible,
        const int degree_l,
        const double G_to_use,
        const size_t num_ys,
        std::complex<double>* starting_conditions_ptr) noexcept
{
    const double l     = static_cast<double>(degree_l);
    const double lp1   = l + 1.0;
    const double llp1  = l * lp1;
    const double gamma = 4.0 * TidalPyConstants::d_PI * G_to_use * density / 3.0;
    const double w2    = frequency * frequency;
    const std::complex<double>& mu = shear_modulus;

    c_SeriesMatrix<6> b0 = c_SeriesMatrix<6>::Zero();
    c_SeriesVector<6> lead_c;
    if (is_incompressible)
    {
        b0(0, 0) = -lp1;
        b0(0, 2) = llp1;
        b0(1, 0) = 12.0 * mu;
        b0(1, 1) = 2.0 - l;
        b0(1, 2) = -6.0 * llp1 * mu;
        b0(3, 0) = -6.0 * mu;
        b0(3, 1) = -1.0;
        b0(3, 2) = 2.0 * mu * (2.0 * llp1 - 1.0);

        // The nu0 = 2 null vector of B0, normalized to z4 = 1.
        lead_c << lp1 / (2.0 * mu * (l + 2.0)),
                  (l * l - l - 3.0) * lp1 / (l * (l + 2.0)),
                  (l + 3.0) / (2.0 * mu * l * (l + 2.0)),
                  1.0,
                  0.0,
                  -3.0 * gamma * lp1 / (2.0 * mu * (l + 2.0));
    }
    else
    {
        const std::complex<double> lame           = bulk_modulus - (2.0 / 3.0) * mu;
        const std::complex<double> lame_2mu       = lame + 2.0 * mu;
        const std::complex<double> three_lame_2mu = 3.0 * lame + 2.0 * mu;
        b0(0, 0) = -(l * lame + 2.0 * l * mu + lame - 2.0 * mu) / lame_2mu;
        b0(0, 1) = 1.0 / lame_2mu;
        b0(0, 2) = llp1 * lame / lame_2mu;
        b0(1, 0) = 4.0 * mu * three_lame_2mu / lame_2mu;
        b0(1, 1) = -(l * lame + 2.0 * l * mu - 2.0 * lame) / lame_2mu;
        b0(1, 2) = -2.0 * llp1 * mu * three_lame_2mu / lame_2mu;
        b0(3, 0) = -2.0 * mu * three_lame_2mu / lame_2mu;
        b0(3, 1) = -lame / lame_2mu;
        b0(3, 2) = 2.0 * mu * (2.0 * llp1 * (lame + mu) - lame - 2.0 * mu) / lame_2mu;

        // Martens' A4,0 solution (thesis Eqs. 4.103-4.106) with A4,0 = 1.
        const std::complex<double> p1  = 2.0 * l * (l * (l + 2.0) * lame + (l * (l + 2.0) - 1.0) * mu);
        const std::complex<double> p2  = l * (l + 5.0) + l * (l + 3.0) * lame / mu;
        const std::complex<double> q1  = (llp1 + l * (l + 3.0)) * lame + 2.0 * llp1 * mu;
        const std::complex<double> q2  = 2.0 * lp1 + (l + 3.0) * lame / mu;
        const std::complex<double> a11 = 1.0 / mu - l * p2 / p1;
        const std::complex<double> a31 = p2 / p1;
        const std::complex<double> a52 = (3.0 * gamma / (2.0 * (2.0 * l + 3.0))) * ((l + 3.0) * a11 - llp1 * a31);
        const std::complex<double> a61 = (l + 2.0) * a52 - 3.0 * gamma * a11;
        lead_c << a11, q2 - q1 * p2 / p1, a31, 1.0, a52, a61 + lp1 * a52;
    }
    b0(1, 3) = llp1;
    b0(2, 0) = -1.0;
    b0(2, 2) = 2.0 - l;
    b0(2, 3) = 1.0 / mu;
    b0(3, 3) = -lp1;
    b0(4, 0) = 3.0 * gamma;
    b0(4, 4) = -(2.0 * l + 1.0);
    b0(4, 5) = 1.0;
    b0(5, 0) = 3.0 * gamma * lp1;
    b0(5, 2) = -3.0 * gamma * llp1;

    c_SeriesMatrix<6> b2 = c_SeriesMatrix<6>::Zero();
    b2(1, 0) = -density * (4.0 * gamma + w2);
    b2(1, 2) = gamma * density * llp1;
    b2(1, 4) = density * lp1;
    b2(1, 5) = -density;
    b2(3, 0) = gamma * density;
    b2(3, 2) = -w2 * density;
    b2(3, 4) = -density;

    // Martens' A6,1 solution (thesis Eq. 4.101, times l), converted, and TS72's polynomial solution
    // (y1 = l r^(l - 1)), which is l A1,1 + (l gamma - w^2 - 3 gamma) A6,1 of her solutions (Eqs. 4.96 and 4.101).
    // B2 annihilates the polynomial solution, so it carries none of the growth the other solutions share when
    // |k^2 r^2| is large.
    c_SeriesVector<6> lead_b;
    lead_b << 0.0, 0.0, 0.0, 0.0, 1.0, 2.0 * l + 1.0;
    const double lead_a_y5 = l * gamma - w2;
    c_SeriesVector<6> lead_a;
    lead_a << l, 2.0 * mu * (l - 1.0) * l, 1.0, 2.0 * mu * (l - 1.0), lead_a_y5,
              (2.0 * l + 1.0) * lead_a_y5 - 3.0 * l * gamma;

    static constexpr double powers[6] = {1.0, 0.0, 1.0, 0.0, 2.0, 1.0};
    const double modulus_scale = std::abs(mu);
    const double gravity_scale = c_power_series_gravity_scale(gamma, w2);
    const double scale[6] = {1.0, modulus_scale, 1.0, modulus_scale, gravity_scale, gravity_scale};
    c_SeriesVector<6> z;
    // The polynomial solution is exact as its first term. Summing its series would only turn the roundoff left by
    // B2 lead_a = 0 into growth along the other solutions.
    c_power_series_write<6>(lead_a, powers, 0, radius, degree_l, 0, num_ys, starting_conditions_ptr);

    // A wavenumber solution's series is summed alongside the roundoff of the others: its error grows as
    // eps exp(2 |k| r) (once with the faster wave, once through the slower solution's own cancellation; measured
    // against extended-precision Takeuchi solutions). Past |k| r = ln(1 / sqrt(eps)) / 2 that takes more than half the
    // digits, so the series refuses there.
    const double wave_bound = 0.25 * std::log(1.0 / TidalPyConstants::d_EPS);
    const double wave_bound2 = wave_bound * wave_bound;
    const double r2 = radius * radius;

    if (is_incompressible)
    {
        // One wavenumber, the shear wave's, carried by the A4,0 solution alone. The A6,1 solution is the exact
        // pressure-potential solution (0, -rho r^l, 0, 0, r^l, (2l + 1) r^(l - 1)), written without a series for
        // the reason the polynomial solution is.
        if (std::abs(w2 / (mu / density)) * r2 > wave_bound2) { return false; }
        z = lead_b;
        z(1) = -density * r2;
        c_power_series_write<6>(z, powers, 0, radius, degree_l, 1, num_ys, starting_conditions_ptr);
        if (!c_power_series_sum<6>(b0, b2, lead_c, 2, nullptr, 0.0, scale, radius, z)) { return false; }
        c_power_series_write<6>(z, powers, 2, radius, degree_l, 2, num_ys, starting_conditions_ptr);
        return true;
    }
    if (!(gamma > 0.0))
    {
        // Without gravity f and h are undefined: Martens' A6,1 and A4,0 solutions.
        if (!c_power_series_sum<6>(b0, b2, lead_b, 0, &lead_c, 0.0, scale, radius, z)) { return false; }
        c_power_series_write<6>(z, powers, 0, radius, degree_l, 1, num_ys, starting_conditions_ptr);
        if (!c_power_series_sum<6>(b0, b2, lead_c, 2, nullptr, 0.0, scale, radius, z)) { return false; }
        c_power_series_write<6>(z, powers, 2, radius, degree_l, 2, num_ys, starting_conditions_ptr);
        return true;
    }

    // The series of TS72's two wavenumber solutions (f and h from c_solid_wave_numbers, as the closed forms use).
    const std::complex<double> lame   = bulk_modulus - (2.0 / 3.0) * mu;
    const std::complex<double> alpha2 = (lame + 2.0 * mu) / density;
    const std::complex<double> beta2  = mu / density;
    const c_SolidWaveNumbers waves = c_solid_wave_numbers(w2, gamma, alpha2, beta2, degree_l);
    if (std::fmax(std::abs(waves.k2_pos), std::abs(waves.k2_neg)) * r2 > wave_bound2) { return false; }
    const std::complex<double> f_values[2] = {waves.f_neg, waves.f_pos};
    const std::complex<double> h_values[2] = {waves.h_neg, waves.h_pos};
    for (size_t wave_i = 0; wave_i < 2; ++wave_i)
    {
        const std::complex<double>& f = f_values[wave_i];
        const std::complex<double>& h = h_values[wave_i];
        const c_SeriesVector<6> lead_wave = (alpha2 * f - lp1 * beta2) * lead_b;
        const std::complex<double> wave_z4 = mu * (1.0 - ((l - 1.0) * h + 2.0 * (f + 1.0)) / (2.0 * l + 3.0));
        if (!c_power_series_sum<6>(b0, b2, lead_wave, 0, &lead_c, wave_z4, scale, radius, z)) { return false; }
        c_power_series_write<6>(z, powers, 0, radius, degree_l, wave_i + 1, num_ys, starting_conditions_ptr);
    }
    return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Liquid Layers
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Starting conditions at the bottom of a dynamic liquid layer: two independent solutions of (y1, y2, y5, y6).
// Units: density [kg m-3], bulk modulus [Pa], frequency [rad s-1]. False when a series did not converge.
inline bool c_power_series_liquid_dynamic(
        const double frequency,
        const double radius,
        const double density,
        const std::complex<double>& bulk_modulus,
        const bool is_incompressible,
        const int degree_l,
        const double G_to_use,
        const size_t num_ys,
        std::complex<double>* starting_conditions_ptr) noexcept
{
    const double l     = static_cast<double>(degree_l);
    const double lp1   = l + 1.0;
    const double llp1  = l * lp1;
    const double gamma = 4.0 * TidalPyConstants::d_PI * G_to_use * density / 3.0;
    const double w2    = frequency * frequency;
    const double lgw   = gamma * l - w2;

    c_SeriesMatrix<4> b0 = c_SeriesMatrix<4>::Zero();
    b0(0, 0) = lp1 * lgw / w2;
    b0(0, 1) = -llp1 / (w2 * density);
    b0(0, 2) = -llp1 / w2;
    b0(1, 0) = density * (gamma * gamma * llp1 - 4.0 * gamma * w2 - w2 * w2) / w2;
    b0(1, 1) = -l * (gamma * lp1 + w2) / w2;
    b0(1, 2) = -density * lp1 * lgw / w2;
    b0(1, 3) = -density;
    b0(2, 0) = 3.0 * gamma;
    b0(2, 2) = -(2.0 * l + 1.0);
    b0(2, 3) = 1.0;
    b0(3, 0) = -3.0 * gamma * lp1 * lgw / w2;
    b0(3, 1) = 3.0 * gamma * llp1 / (w2 * density);
    b0(3, 2) = 3.0 * gamma * llp1 / w2;

    // The incompressible liquid's series ends after its first term.
    c_SeriesMatrix<4> b2 = c_SeriesMatrix<4>::Zero();
    if (!is_incompressible) { b2(0, 1) = 1.0 / bulk_modulus; }

    // As in a solid: one solution with z1 = 0 at the center, and TS72's polynomial solution (y1 = l r^(l - 1),
    // y2 = 0), whose series ends after one term.
    c_SeriesVector<4> lead_b;
    lead_b << 0.0, -density, 1.0, 2.0 * l + 1.0;
    const double lead_a_y5 = l * gamma - w2;
    c_SeriesVector<4> lead_a;
    lead_a << l, 0.0, lead_a_y5, (2.0 * l + 1.0) * lead_a_y5 - 3.0 * l * gamma;

    static constexpr double powers[4] = {1.0, 2.0, 2.0, 1.0};
    const double gravity_scale = c_power_series_gravity_scale(gamma, w2);
    const double scale[4] = {1.0, density * gravity_scale, gravity_scale, gravity_scale};
    c_SeriesVector<4> z;
    // The polynomial solution is exact as its first term (see c_power_series_solid).
    c_power_series_write<4>(lead_a, powers, 0, radius, degree_l, 0, num_ys, starting_conditions_ptr);
    if (!c_power_series_sum<4>(b0, b2, lead_b, 0, nullptr, 0.0, scale, radius, z)) { return false; }
    c_power_series_write<4>(z, powers, 0, radius, degree_l, 1, num_ys, starting_conditions_ptr);
    return true;
}
