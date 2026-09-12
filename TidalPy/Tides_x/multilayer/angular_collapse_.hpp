// angular_collapse_.hpp - Analytic colatitude collapse of the secular 3D tidal heating.
//
// The numerical collapse (world_tides_.hpp) integrates the secular heating density over colatitude on a
// Gauss-Legendre theta grid. This header does that integral analytically, using the precomputed angular Gram
// table (angular_gram_.hpp) for two waves of the same degree and a Gauss-Legendre cross-degree Gram matrix
// (exact: the integrands are polynomials in cos theta) for two waves of different degrees, so no theta grid or
// per-point Legendre evaluation is needed on the collapse itself.
//
// Per coherent wave the reduced (unit-amplitude, phasor-free) strain/stress components are
//     eps~_k = sum_j Ce[k][j] f_j(theta),    sigma~_k = sum_j Cs[k][j] f_j(theta)
// with complex radial coefficients Ce, Cs and the real six-function angular basis of angular_gram_.hpp:
//   f1 = P, f2 = dP/dtheta, f3 = d2P/dtheta2, f4 = P/sin,
//   f5 = -m^2 P/sin^2 + cot dP  (phi-phi operator), f6 = (dP - cot P)/sin  (theta-phi operator).
// The i*mu factors of the r-phi / theta-phi strains (mu = azimuthal_sign * m, the signed azimuthal wavenumber of
// the wave) are folded into the coefficients, keeping the basis real.
//
// The secular density at a frequency is (|omega|/2) Im(sigma_c : conj(eps_c)) of the total amplitude at that
// frequency, so with several waves a, b sharing (|omega|, mu) the colatitude integral of the longitude-mean density
// is the sum over ordered pairs of
//     int_0^pi Im( c_a conj(c_b) sum_k w_k sigma~_k^(a) conj(eps~_k^(b)) ) sin(theta) dtheta
//         = sum_k w_k sum_ij Im( c_a conj(c_b) Cs_a[k][i] conj(Ce_b[k][j]) ) G_ij(l_a, l_b, m)
// with G_ij(l_a, l_b, m) = int_0^pi f_i^{l_a} f_j^{l_b} sin(theta) dtheta. The single-wave integral is the
// a == b term with c_a conj(c_b) = |c_a|^2.
#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

#include "strain_radial_.hpp"
#include "angular_gram_.hpp"
#include "legendre_driver_.hpp"   // tidalpy::c_legendre (P_lm and its theta derivatives)
#include "constants_.hpp"         // TidalPyConstants::d_PI


namespace tidalpy {
namespace tides {

// Gauss-Legendre nodes/weights on [-1, 1] (Newton-Raphson on the Legendre polynomial). Used for the colatitude
// integral via the substitution x = cos(theta): int_0^pi f(theta) sin(theta) dtheta = sum_i w[i] f(acos(x[i])).
inline void c_gauss_legendre_nodes(
        int num_nodes,
        std::vector<double>& nodes,
        std::vector<double>& weights) {
    nodes.resize(num_nodes);
    weights.resize(num_nodes);
    const double pi = TidalPyConstants::d_PI;
    const int half = (num_nodes + 1) / 2;
    for (int i = 0; i < half; ++i) {
        double x = std::cos(pi * (static_cast<double>(i) + 0.75) / (static_cast<double>(num_nodes) + 0.5));
        double dp = 1.0;
        for (int iter = 0; iter < 100; ++iter) {
            double p0 = 1.0;   // P_0
            double p1 = x;     // P_1
            for (int k = 2; k <= num_nodes; ++k) {
                const double p2 = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / static_cast<double>(k);
                p0 = p1;
                p1 = p2;
            }
            dp = static_cast<double>(num_nodes) * (x * p1 - p0) / (x * x - 1.0);
            const double dx = -p1 / dp;
            x += dx;
            if (std::abs(dx) <= 1.0e-15) { break; }
        }
        nodes[i]                = -x;
        nodes[num_nodes - 1 - i] = x;
        const double w = 2.0 / ((1.0 - x * x) * dp * dp);
        weights[i]                = w;
        weights[num_nodes - 1 - i] = w;
    }
}

// The six angular basis functions of one (l, m) at a colatitude. For m = 0 the f4 and f6 functions are
// unbounded at the poles but enter the strain only through the factor i*m = 0, so they are returned as zero
// (matching the precomputed table).
inline void c_angular_basis_at(int degree_l, int order_m, double colatitude, double basis[6]) {
    const c_LegendreValue legendre = c_legendre(degree_l, order_m, colatitude);
    const double sin_t = std::sin(colatitude);
    const double cot_t = std::cos(colatitude) / sin_t;
    basis[0] = legendre.p;
    basis[1] = legendre.dp_dtheta;
    basis[2] = legendre.d2p_dtheta2;
    if (order_m == 0) {
        basis[3] = 0.0;
        basis[4] = cot_t * legendre.dp_dtheta;
        basis[5] = 0.0;
    } else {
        const double m2 = static_cast<double>(order_m) * static_cast<double>(order_m);
        basis[3] = legendre.p / sin_t;
        basis[4] = -m2 * legendre.p / (sin_t * sin_t) + cot_t * legendre.dp_dtheta;
        basis[5] = (legendre.dp_dtheta - cot_t * legendre.p) / sin_t;
    }
}

// Gauss-Legendre order for the cross-degree Gram integrals. The integrands f_i^{l_a} f_j^{l_b} sin(theta) are
// polynomials in cos(theta) of degree at most l_a + l_b + 2 <= 22 for l <= 10, so 32 nodes (exact to degree 63)
// integrate them exactly; the margin keeps the near-pole cancellations in f5/f6 well away from the nodes.
inline constexpr int C_ANGULAR_GRAM_CROSS_NODES = 32;

// The angular Gram matrix G_ij(l_a, l_b, m) = int_0^pi f_i^{l_a}(theta) f_j^{l_b}(theta) sin(theta) dtheta between
// two degrees at the same order m. Same degree -> the precomputed table; different degrees -> Gauss-Legendre
// quadrature in x = cos(theta) (exact, see above). Not symmetric for l_a != l_b: G_ij(l_a, l_b) = G_ji(l_b, l_a).
// Returns false if a degree or the order is out of the tabulated range (l = 2..10, m = 0..min(l_a, l_b)).
inline bool c_angular_gram_pair(int degree_a, int degree_b, int order_m, double gram[6][6]) {
    if (degree_a == degree_b) {
        return c_angular_gram(degree_a, order_m, gram);
    }
    if (degree_a < C_ANGULAR_GRAM_MIN_L || degree_a > C_ANGULAR_GRAM_MAX_L
        || degree_b < C_ANGULAR_GRAM_MIN_L || degree_b > C_ANGULAR_GRAM_MAX_L
        || order_m < 0 || order_m > degree_a || order_m > degree_b) {
        return false;
    }
    std::vector<double> nodes;
    std::vector<double> weights;
    c_gauss_legendre_nodes(C_ANGULAR_GRAM_CROSS_NODES, nodes, weights);
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) { gram[i][j] = 0.0; }
    }
    double basis_a[6];
    double basis_b[6];
    for (int n = 0; n < C_ANGULAR_GRAM_CROSS_NODES; ++n) {
        const double colatitude = std::acos(nodes[n]);
        c_angular_basis_at(degree_a, order_m, colatitude, basis_a);
        c_angular_basis_at(degree_b, order_m, colatitude, basis_b);
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                gram[i][j] += weights[n] * basis_a[i] * basis_b[j];
            }
        }
    }
    return true;
}

// Reduced strain (Ce) and stress (Cs) coefficient matrices [component k][basis j] of a unit-amplitude wave with
// signed azimuthal wavenumber mu, from its radial coefficients. Only the nonzero entries are set.
inline void c_reduced_strain_stress_coeffs(
        const c_StrainRadialCoeffs& radial,
        double mu,
        std::complex<double> strain_coeff[6][6],
        std::complex<double> stress_coeff[6][6]) {
    const std::complex<double> imag_unit(0.0, 1.0);
    const std::complex<double> c_zero(0.0, 0.0);
    for (int k = 0; k < 6; ++k) {
        for (int j = 0; j < 6; ++j) { strain_coeff[k][j] = c_zero; }
    }
    strain_coeff[0][0] = radial.dy1_dr;                            // eps_rr   = dy1_dr * f1
    strain_coeff[1][0] = radial.y1_over_r;                         // eps_thth = y1/r * f1
    strain_coeff[1][2] = radial.y3_over_r;                         //           + y3/r * f3
    strain_coeff[2][0] = radial.y1_over_r;                         // eps_phph = y1/r * f1
    strain_coeff[2][4] = radial.y3_over_r;                         //           + y3/r * f5
    strain_coeff[3][1] = radial.y4_over_2mu;                       // eps_rth  = y4/2mu * f2
    strain_coeff[4][3] = imag_unit * mu * radial.y4_over_2mu;      // eps_rph  = i mu y4/2mu * f4
    strain_coeff[5][5] = 2.0 * imag_unit * mu * radial.y3_over_2r; // eps_thph = 2 i mu y3/2r * f6

    // Trace coefficients T_j = Ce[0][j] + Ce[1][j] + Ce[2][j].
    std::complex<double> trace[6];
    for (int j = 0; j < 6; ++j) {
        trace[j] = strain_coeff[0][j] + strain_coeff[1][j] + strain_coeff[2][j];
    }
    // Stress: sigma = 2 mu eps + lame tr(eps) delta -> Cs[k][j] = 2 mu Ce[k][j] + (k<3) lame T[j].
    const std::complex<double> two_mu = 2.0 * radial.shear;
    for (int k = 0; k < 6; ++k) {
        for (int j = 0; j < 6; ++j) {
            stress_coeff[k][j] = two_mu * strain_coeff[k][j] + ((k < 3) ? radial.lame * trace[j] : c_zero);
        }
    }
}

// Colatitude integral of the secular heating between two coherent waves a and b that share the same order m,
// azimuthal sign, and |omega| (possibly different degrees):
//     sum_k w_k sum_ij Im( c_pair * Cs_a[k][i] conj(Ce_b[k][j]) ) G_ij
// with c_pair = amplitude_a * conj(amplitude_b) and gram the (cross-)degree Gram matrix G_ij(l_a, l_b, m). The
// caller multiplies by 0.5*|omega| (and r^2, 2*pi) for the radial power contribution and sums over the ordered
// pairs of its group. a == b with c_pair = |amplitude|^2 is the single-wave integral.
inline double c_theta_integrated_heating_pair(
        const c_StrainRadialCoeffs& radial_a,
        const c_StrainRadialCoeffs& radial_b,
        int order_m,
        int azimuthal_sign,
        std::complex<double> c_pair,
        const double gram[6][6]) noexcept
{
    const double mu = static_cast<double>(azimuthal_sign) * static_cast<double>(order_m);
    std::complex<double> strain_a[6][6];
    std::complex<double> stress_a[6][6];
    std::complex<double> strain_b[6][6];
    std::complex<double> stress_b[6][6];
    c_reduced_strain_stress_coeffs(radial_a, mu, strain_a, stress_a);
    c_reduced_strain_stress_coeffs(radial_b, mu, strain_b, stress_b);

    double total = 0.0;
    for (int k = 0; k < 6; ++k) {
        const double weight = (k < 3) ? 1.0 : 2.0;
        std::complex<double> accum(0.0, 0.0);
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                accum += stress_a[k][i] * std::conj(strain_b[k][j]) * gram[i][j];
            }
        }
        total += weight * std::imag(c_pair * accum);
    }
    return total;
}

// Colatitude integral of one (l, m) wave's secular heating angular factor on its own (unit amplitude, positive
// azimuthal sign): sum_k w_k sum_ij Im(Cs[k][i] conj(Ce[k][j])) G_ij(l, m). The caller multiplies by
// 0.5*|omega|*|amplitude|^2 (and r^2, 2*pi). NaN if (l, m) is out of table range. Kept for single-wave use; the
// coherent collapse goes through c_theta_integrated_heating_pair.
inline double c_theta_integrated_heating(
        const c_StrainRadialCoeffs& radial,
        int degree_l,
        int order_m) noexcept
{
    double gram[6][6];
    if (!c_angular_gram(degree_l, order_m, gram))
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return c_theta_integrated_heating_pair(radial, radial, order_m, 1, std::complex<double>(1.0, 0.0), gram);
}

}  // namespace tides
}  // namespace tidalpy
