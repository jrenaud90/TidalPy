// takeuchi_.hpp - Starting conditions based on Takeuchi & Saito (1972)
// Ported from TidalPy/RadialSolver/starting/takeuchi.pyx
//
// References
// ----------
// TS72: Takeuchi & Saito (1972) Eqs. 95-102
#pragma once

#include <cmath>
#include <complex>

#include "../../constants_.hpp"
#include "common_.hpp"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Solid Layers
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Calculate the starting guess at the bottom of a solid layer using the dynamic assumption.
// TS72 Eqs. 95-102.
// Three independent solutions (sn1, sn2, sn3).
inline void c_takeuchi_solid_dynamic_compressible(
        const double frequency,
        const double radius,
        const double density,
        const std::complex<double>& bulk_modulus,
        const std::complex<double>& shear_modulus,
        const int degree_l,
        const double G_to_use,
        const size_t num_ys,
        std::complex<double>* starting_conditions_ptr) noexcept
{
    // Convert compressibility parameters
    const std::complex<double> lame = bulk_modulus - (2.0 / 3.0) * shear_modulus;

    // Constants
    const double gamma           = 4.0 * TidalPyConstants::d_PI * G_to_use * density / 3.0;
    const double dynamic_term    = frequency * frequency;
    const std::complex<double> alpha2 = (lame + 2.0 * shear_modulus) / density;
    const std::complex<double> beta2  = shear_modulus / density;

    // Optimizations
    const double r_inverse    = 1.0 / radius;
    const double r2           = radius * radius;
    const double degree_l_dbl = static_cast<double>(degree_l);
    const double lp1          = degree_l_dbl + 1.0;
    const double lm1          = degree_l_dbl - 1.0;
    const double dlp1         = 2.0 * degree_l_dbl + 1.0;
    const double dlp3         = 2.0 * degree_l_dbl + 3.0;
    const double llp1         = degree_l_dbl * lp1;

    // Helper functions - See Eq. 99 of TS72
    const std::complex<double> k2_quad_pos  = (dynamic_term / beta2) + ((dynamic_term + 4.0 * gamma) / alpha2);
    const std::complex<double> k2_quad_neg  = (dynamic_term / beta2) - ((dynamic_term + 4.0 * gamma) / alpha2);
    const std::complex<double> k2_quad      = (k2_quad_neg * k2_quad_neg) + ((4.0 * llp1 * gamma * gamma) / (alpha2 * beta2));
    const std::complex<double> k2_quad_sqrt = std::sqrt(k2_quad);

    // QUESTION: (Issue #43) TS74 has these flipped compared to KMN15. Going with TS74 for this func.
    // See the -/+ order in TS72 EQ. 99
    const size_t neg_index = 0;
    const size_t pos_index = 1;
    const std::complex<double> k2_neg = (1.0 / 2.0) * (k2_quad_pos - k2_quad_sqrt);
    const std::complex<double> k2_pos = (1.0 / 2.0) * (k2_quad_pos + k2_quad_sqrt);

    const std::complex<double> f_pos = (beta2 * k2_pos - dynamic_term) / gamma;
    const std::complex<double> f_neg = (beta2 * k2_neg - dynamic_term) / gamma;

    const std::complex<double> h_pos = f_pos - lp1;
    const std::complex<double> h_neg = f_neg - lp1;

    // Calculate Takeuchi and Saito functions
    const std::complex<double> z2_pos = k2_pos * r2;
    const std::complex<double> z2_neg = k2_neg * r2;

    std::complex<double> phi_pos, phi_lp1_pos, psi_pos;
    std::complex<double> phi_neg, phi_lp1_neg, psi_neg;
    c_takeuchi_phi_psi(z2_pos, degree_l, &phi_pos, &phi_lp1_pos, &psi_pos);
    c_takeuchi_phi_psi(z2_neg, degree_l, &phi_neg, &phi_lp1_neg, &psi_neg);

    // See Eq. 102 in TS72

    // y1, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 0] =
        ((-std::pow(radius, lp1)) / dlp3) * (0.5 * degree_l_dbl * h_pos * psi_pos + f_pos * phi_lp1_pos);
    starting_conditions_ptr[neg_index * num_ys + 0] =
        ((-std::pow(radius, lp1)) / dlp3) * (0.5 * degree_l_dbl * h_neg * psi_neg + f_neg * phi_lp1_neg);
    starting_conditions_ptr[2 * num_ys + 0] = degree_l_dbl * std::pow(radius, lm1);

    // y2, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 1] =
        -(lame + 2.0 * shear_modulus) * std::pow(radius, degree_l_dbl) * f_pos * phi_pos +
        (shear_modulus * std::pow(radius, degree_l_dbl) / dlp3) * (
            -degree_l_dbl * lm1 * h_pos * psi_pos + 2.0 * (2.0 * f_pos + llp1) * phi_lp1_pos
            );
    starting_conditions_ptr[neg_index * num_ys + 1] =
        -(lame + 2.0 * shear_modulus) * std::pow(radius, degree_l_dbl) * f_neg * phi_neg +
        (shear_modulus * std::pow(radius, degree_l_dbl) / dlp3) * (
            -degree_l_dbl * lm1 * h_neg * psi_neg + 2.0 * (2.0 * f_neg + llp1) * phi_lp1_neg
            );
    starting_conditions_ptr[2 * num_ys + 1] =
        2.0 * shear_modulus * degree_l_dbl * lm1 * std::pow(radius, degree_l_dbl - 2.0);

    // y3, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 2] =
        (-std::pow(radius, lp1) / dlp3) * (0.5 * h_pos * psi_pos - phi_lp1_pos);
    starting_conditions_ptr[neg_index * num_ys + 2] =
        (-std::pow(radius, lp1) / dlp3) * (0.5 * h_neg * psi_neg - phi_lp1_neg);
    starting_conditions_ptr[2 * num_ys + 2] =
        std::pow(radius, lm1);

    // y4, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 3] =
        shear_modulus * std::pow(radius, degree_l_dbl) * (
            phi_pos - (1.0 / dlp3) * (lm1 * h_pos * psi_pos + 2.0 * (f_pos + 1.0) * phi_lp1_pos)
            );
    starting_conditions_ptr[neg_index * num_ys + 3] =
        shear_modulus * std::pow(radius, degree_l_dbl) * (
            phi_neg - (1.0 / dlp3) * (lm1 * h_neg * psi_neg + 2.0 * (f_neg + 1.0) * phi_lp1_neg)
            );
    starting_conditions_ptr[2 * num_ys + 3] =
        2.0 * shear_modulus * lm1 * std::pow(radius, degree_l_dbl - 2.0);

    // y5, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 4] =
        std::pow(radius, degree_l_dbl + 2.0) * (
            (alpha2 * f_pos - lp1 * beta2) / r2 - (3.0 * gamma * f_pos / (2.0 * dlp3)) * psi_pos
            );
    starting_conditions_ptr[neg_index * num_ys + 4] =
        std::pow(radius, degree_l_dbl + 2.0) * (
            (alpha2 * f_neg - lp1 * beta2) / r2 - (3.0 * gamma * f_neg / (2.0 * dlp3)) * psi_neg
            );
    starting_conditions_ptr[2 * num_ys + 4] =
        (degree_l_dbl * gamma - dynamic_term) * std::pow(radius, degree_l_dbl);

    // y6, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 5] =
        dlp1 * r_inverse * starting_conditions_ptr[pos_index * num_ys + 4] +
        (3.0 * degree_l_dbl * gamma * h_pos * std::pow(radius, lp1) / (2.0 * dlp3)) * psi_pos;
    starting_conditions_ptr[neg_index * num_ys + 5] =
        dlp1 * r_inverse * starting_conditions_ptr[neg_index * num_ys + 4] +
        (3.0 * degree_l_dbl * gamma * h_neg * std::pow(radius, lp1) / (2.0 * dlp3)) * psi_neg;
    starting_conditions_ptr[2 * num_ys + 5] =
        dlp1 * r_inverse * starting_conditions_ptr[2 * num_ys + 4] -
        3.0 * degree_l_dbl * gamma * std::pow(radius, lm1);
}


// Calculate the starting guess at the bottom of a solid layer using the static assumption.
// TS72 Eqs. 95-102 (with w=0).
// Three independent solutions (sn1, sn2, sn3).
inline void c_takeuchi_solid_static_compressible(
        const double radius,
        const double density,
        const std::complex<double>& bulk_modulus,
        const std::complex<double>& shear_modulus,
        const int degree_l,
        const double G_to_use,
        const size_t num_ys,
        std::complex<double>* starting_conditions_ptr) noexcept
{
    // Convert compressibility parameters
    const std::complex<double> lame = bulk_modulus - (2.0 / 3.0) * shear_modulus;

    // Constants
    const double gamma          = 4.0 * TidalPyConstants::d_PI * G_to_use * density / 3.0;
    const std::complex<double> alpha2 = (lame + 2.0 * shear_modulus) / density;
    const std::complex<double> beta2  = shear_modulus / density;

    // Optimizations
    const double r_inverse    = 1.0 / radius;
    const double r2           = radius * radius;
    const double degree_l_dbl = static_cast<double>(degree_l);
    const double lp1          = degree_l_dbl + 1.0;
    const double lm1          = degree_l_dbl - 1.0;
    const double dlp1         = 2.0 * degree_l_dbl + 1.0;
    const double dlp3         = 2.0 * degree_l_dbl + 3.0;
    const double llp1         = degree_l_dbl * lp1;

    // Helper functions - See Eq. 99 of TS72
    const std::complex<double> k2_quad_pos  = (4.0 * gamma / alpha2);
    const std::complex<double> k2_quad      = k2_quad_pos * k2_quad_pos + ((4.0 * llp1 * gamma * gamma) / (alpha2 * beta2));
    const std::complex<double> k2_quad_sqrt = std::sqrt(k2_quad);

    // QUESTION: (Issue #43) TS74 has these flipped compared to KMN15. Going with TS74 for this func.
    const size_t neg_index = 0;
    const size_t pos_index = 1;
    const std::complex<double> k2_neg = (1.0 / 2.0) * (k2_quad_pos - k2_quad_sqrt);
    const std::complex<double> k2_pos = (1.0 / 2.0) * (k2_quad_pos + k2_quad_sqrt);

    const std::complex<double> f_pos = (beta2 * k2_pos) / gamma;
    const std::complex<double> f_neg = (beta2 * k2_neg) / gamma;

    const std::complex<double> h_pos = f_pos - lp1;
    const std::complex<double> h_neg = f_neg - lp1;

    // Calculate Takeuchi and Saito functions
    const std::complex<double> z2_pos = k2_pos * r2;
    const std::complex<double> z2_neg = k2_neg * r2;

    std::complex<double> phi_pos, phi_lp1_pos, psi_pos;
    std::complex<double> phi_neg, phi_lp1_neg, psi_neg;
    c_takeuchi_phi_psi(z2_pos, degree_l, &phi_pos, &phi_lp1_pos, &psi_pos);
    c_takeuchi_phi_psi(z2_neg, degree_l, &phi_neg, &phi_lp1_neg, &psi_neg);

    // See Eq. 102 in TS72
    // y1, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 0] =
        ((-std::pow(radius, lp1)) / dlp3) * (0.5 * degree_l_dbl * h_pos * psi_pos + f_pos * phi_lp1_pos);
    starting_conditions_ptr[neg_index * num_ys + 0] =
        ((-std::pow(radius, lp1)) / dlp3) * (0.5 * degree_l_dbl * h_neg * psi_neg + f_neg * phi_lp1_neg);
    starting_conditions_ptr[2 * num_ys + 0] =
        degree_l_dbl * std::pow(radius, lm1);

    // y2, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 1] =
        -(lame + 2.0 * shear_modulus) * std::pow(radius, degree_l_dbl) * f_pos * phi_pos +
        (shear_modulus * std::pow(radius, degree_l_dbl) / dlp3) * (
            -degree_l_dbl * lm1 * h_pos * psi_pos + 2.0 * (2.0 * f_pos + llp1) * phi_lp1_pos
            );
    starting_conditions_ptr[neg_index * num_ys + 1] =
        -(lame + 2.0 * shear_modulus) * std::pow(radius, degree_l_dbl) * f_neg * phi_neg +
        (shear_modulus * std::pow(radius, degree_l_dbl) / dlp3) * (
            -degree_l_dbl * lm1 * h_neg * psi_neg + 2.0 * (2.0 * f_neg + llp1) * phi_lp1_neg
            );
    starting_conditions_ptr[2 * num_ys + 1] =
        2.0 * shear_modulus * degree_l_dbl * lm1 * std::pow(radius, degree_l_dbl - 2.0);

    // y3, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 2] =
        (-std::pow(radius, lp1) / dlp3) * (0.5 * h_pos * psi_pos - phi_lp1_pos);
    starting_conditions_ptr[neg_index * num_ys + 2] =
        (-std::pow(radius, lp1) / dlp3) * (0.5 * h_neg * psi_neg - phi_lp1_neg);
    starting_conditions_ptr[2 * num_ys + 2] =
        std::pow(radius, lm1);

    // y4, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 3] =
        shear_modulus * std::pow(radius, degree_l_dbl) *
        (phi_pos - (1.0 / dlp3) * (lm1 * h_pos * psi_pos + 2.0 * (f_pos + 1.0) * phi_lp1_pos));
    starting_conditions_ptr[neg_index * num_ys + 3] =
        shear_modulus * std::pow(radius, degree_l_dbl) *
        (phi_neg - (1.0 / dlp3) * (lm1 * h_neg * psi_neg + 2.0 * (f_neg + 1.0) * phi_lp1_neg));
    starting_conditions_ptr[2 * num_ys + 3] =
        2.0 * shear_modulus * lm1 * std::pow(radius, degree_l_dbl - 2.0);

    // y5, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 4] =
        std::pow(radius, degree_l_dbl + 2.0) *
        ((alpha2 * f_pos - lp1 * beta2) / r2 - (3.0 * gamma * f_pos / (2.0 * dlp3)) * psi_pos);
    starting_conditions_ptr[neg_index * num_ys + 4] =
        std::pow(radius, degree_l_dbl + 2.0) *
        ((alpha2 * f_neg - lp1 * beta2) / r2 - (3.0 * gamma * f_neg / (2.0 * dlp3)) * psi_neg);
    starting_conditions_ptr[2 * num_ys + 4] =
        degree_l_dbl * gamma * std::pow(radius, degree_l_dbl);

    // y6, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 5] =
        dlp1 * r_inverse * starting_conditions_ptr[pos_index * num_ys + 4] +
        (3.0 * degree_l_dbl * gamma * h_pos * std::pow(radius, lp1) / (2.0 * dlp3)) * psi_pos;
    starting_conditions_ptr[neg_index * num_ys + 5] =
        dlp1 * r_inverse * starting_conditions_ptr[neg_index * num_ys + 4] +
        (3.0 * degree_l_dbl * gamma * h_neg * std::pow(radius, lp1) / (2.0 * dlp3)) * psi_neg;
    starting_conditions_ptr[2 * num_ys + 5] =
        dlp1 * r_inverse * starting_conditions_ptr[2 * num_ys + 4] -
        3.0 * degree_l_dbl * gamma * std::pow(radius, lm1);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Liquid Layers
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Calculate the starting guess at the bottom of a liquid layer using the dynamic assumption.
// TS72 Eqs. 95-102 (mu=0).
// Two independent solutions (sn1, sn2).
inline void c_takeuchi_liquid_dynamic_compressible(
        const double frequency,
        const double radius,
        const double density,
        const std::complex<double>& bulk_modulus,
        const int degree_l,
        const double G_to_use,
        const size_t num_ys,
        std::complex<double>* starting_conditions_ptr) noexcept
{
    // For liquid, shear modulus = 0 so lame = bulk_modulus
    const std::complex<double> lame = bulk_modulus;

    // Constants
    const double dynamic_term   = frequency * frequency;
    const double gamma          = 4.0 * TidalPyConstants::d_PI * G_to_use * density / 3.0;
    const std::complex<double> alpha2 = lame / density;

    // Optimizations
    const double r_inverse    = 1.0 / radius;
    const double r2           = radius * radius;
    const double degree_l_dbl = static_cast<double>(degree_l);
    const double lp1          = degree_l_dbl + 1.0;
    const double lm1          = degree_l_dbl - 1.0;
    const double dlp1         = 2.0 * degree_l_dbl + 1.0;
    const double dlp3         = 2.0 * degree_l_dbl + 3.0;
    const double llp1         = degree_l_dbl * lp1;

    // k2, h, and f no longer depend on k2. See Eq. 101 of TS72
    const double f  = -dynamic_term / gamma;
    const double h  = f - lp1;
    const std::complex<double> k2 = (1.0 / alpha2) * (dynamic_term + 4.0 * gamma - llp1 * gamma * gamma / dynamic_term);

    // Calculate Takeuchi and Saito functions
    const std::complex<double> z = k2 * r2;
    std::complex<double> phi, phi_lp1, psi;
    c_takeuchi_phi_psi(z, degree_l, &phi, &phi_lp1, &psi);

    // Found by setting mu=0 in Eq. 102 of TS72
    // y1, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 0] =
        ((-std::pow(radius, lp1)) / dlp3) * (0.5 * degree_l_dbl * h * psi + f * phi_lp1);
    starting_conditions_ptr[1 * num_ys + 0] =
        degree_l_dbl * std::pow(radius, lm1);

    // y2, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 1] =
        -lame * std::pow(radius, degree_l_dbl) * f * phi;
    starting_conditions_ptr[1 * num_ys + 1] =
        std::complex<double>(0.0, 0.0);

    // y5, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 2] =
        std::pow(radius, degree_l_dbl + 2.0) * ((alpha2 * f / r2) - (3.0 * gamma * f / (2.0 * dlp3)) * psi);
    starting_conditions_ptr[1 * num_ys + 2] =
        (degree_l_dbl * gamma - dynamic_term) * std::pow(radius, degree_l_dbl);

    // y6, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 3] =
        dlp1 * r_inverse * starting_conditions_ptr[0 * num_ys + 2] + (3.0 * degree_l_dbl * gamma * h * std::pow(radius, lp1) / (2.0 * dlp3)) * psi;
    starting_conditions_ptr[1 * num_ys + 3] =
        dlp1 * r_inverse * starting_conditions_ptr[1 * num_ys + 2] - 3.0 * degree_l_dbl * gamma * std::pow(radius, lm1);
}
