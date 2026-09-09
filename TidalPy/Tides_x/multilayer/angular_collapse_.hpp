// angular_collapse_.hpp - Analytic colatitude collapse of the secular 3D tidal heating.
//
// The numerical collapse (world_tides_.hpp) integrates the secular heating density over colatitude on a
// Gauss-Legendre theta grid. This header does that integral ANALYTICALLY for one (l, m) mode, using the
// precomputed angular Gram table (angular_gram_.hpp), so no theta grid or Legendre evaluation is needed.
//
// Per mode the secular density factorizes as
//     h_mode(r, theta) / amplitude^2 = sum_k w_k Im( sigma_tilde_k(r, theta) conj(eps_tilde_k(r, theta)) )
// where each reduced (phasor-free) strain/stress component is (complex radial coefficient) x (real
// angular function of theta). The six angular functions form the bounded basis of angular_gram_.hpp:
//   f1 = P, f2 = dP/dtheta, f3 = d2P/dtheta2, f4 = P/sin,
//   f5 = -m^2 P/sin^2 + cot dP  (phi-phi operator), f6 = (dP - cot P)/sin  (theta-phi operator).
// Writing eps_tilde_k = sum_j Ce[k][j] f_j and sigma_tilde_k = sum_j Cs[k][j] f_j (complex radial
// coefficients), the colatitude integral is
//     int_0^pi (h_mode / amplitude^2) sin(theta) dtheta
//         = sum_k w_k sum_ij Im( Cs[k][i] conj(Ce[k][j]) ) G_ij(l, m)
// with G_ij(l, m) = int_0^pi f_i f_j sin(theta) dtheta the Gram table. The i*m factors of the rphi/thphi
// strains are folded into the (complex) radial coefficients, keeping the basis real.
#pragma once

#include <complex>
#include <limits>

#include "strain_radial_.hpp"
#include "angular_gram_.hpp"


namespace tidalpy {
namespace tides {

// Colatitude integral of one (l, m) mode's secular heating angular factor (see the header comment).
// Returns sum_k w_k sum_ij Im(Cs[k][i] conj(Ce[k][j])) G_ij; the caller multiplies by 0.5*|omega|*
// amplitude^2 (and r^2, 2*pi) to get the radial power contribution. NaN if (l, m) is out of table range.
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

    const double order_m_d = static_cast<double>(order_m);
    const std::complex<double> imag_unit(0.0, 1.0);
    const std::complex<double> c_zero(0.0, 0.0);

    // Reduced strain coefficients eps_tilde_k = sum_j Ce[k][j] f_j (only the nonzero entries are set).
    std::complex<double> strain_coeff[6][6];
    for (int k = 0; k < 6; ++k)
    {
        for (int j = 0; j < 6; ++j) { strain_coeff[k][j] = c_zero; }
    }
    strain_coeff[0][0] = radial.dy1_dr;                          // eps_rr   = dy1_dr * f1
    strain_coeff[1][0] = radial.y1_over_r;                       // eps_thth = y1/r * f1
    strain_coeff[1][2] = radial.y3_over_r;                       //           + y3/r * f3
    strain_coeff[2][0] = radial.y1_over_r;                       // eps_phph = y1/r * f1
    strain_coeff[2][4] = radial.y3_over_r;                       //           + y3/r * f5
    strain_coeff[3][1] = radial.y4_over_2mu;                     // eps_rth  = y4/2mu * f2
    strain_coeff[4][3] = imag_unit * order_m_d * radial.y4_over_2mu;        // eps_rph  = i m y4/2mu * f4
    strain_coeff[5][5] = 2.0 * imag_unit * order_m_d * radial.y3_over_2r;   // eps_thph = 2 i m y3/2r * f6

    // Trace coefficients T_j = Ce[0][j] + Ce[1][j] + Ce[2][j].
    std::complex<double> trace[6];
    for (int j = 0; j < 6; ++j)
    {
        trace[j] = strain_coeff[0][j] + strain_coeff[1][j] + strain_coeff[2][j];
    }

    // Stress: sigma = 2 mu eps + lame tr(eps) delta -> Cs[k][j] = 2 mu Ce[k][j] + (k<3) lame T[j].
    std::complex<double> stress_coeff[6][6];
    const std::complex<double> two_mu = 2.0 * radial.shear;
    for (int k = 0; k < 6; ++k)
    {
        for (int j = 0; j < 6; ++j)
        {
            stress_coeff[k][j] = two_mu * strain_coeff[k][j]
                               + ((k < 3) ? radial.lame * trace[j] : c_zero);
        }
    }

    double total = 0.0;
    for (int k = 0; k < 6; ++k)
    {
        const double weight = (k < 3) ? 1.0 : 2.0;
        double accum = 0.0;
        for (int i = 0; i < 6; ++i)
        {
            for (int j = 0; j < 6; ++j)
            {
                accum += std::imag(stress_coeff[k][i] * std::conj(strain_coeff[k][j])) * gram[i][j];
            }
        }
        total += weight * accum;
    }
    return total;
}

}  // namespace tides
}  // namespace tidalpy
