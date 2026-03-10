// kamata_.hpp - Starting conditions based on Kamata et al. (2015)
// Ported from TidalPy/RadialSolver/starting/kamata.pyx
//
// References
// ----------
// KMN15: Kamata, Matsuyama, & Nimmo (2015; JGR:P) Eqs. B1-B37
#pragma once

#include <cmath>
#include <complex>

#include "../../constants_.hpp"
#include "common_.hpp"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Solid Layers
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Calculate the starting guess at the bottom of a solid layer using the dynamic assumption.
// KMN15 Eqs. B1-B16.
// Three independent solutions (sn1, sn2, sn3).
inline void c_kamata_solid_dynamic_compressible(
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
    const double gamma          = 4.0 * TidalPyConstants::d_PI * G_to_use * density / 3.0;
    const double dynamic_term   = frequency * frequency;
    const std::complex<double> alpha2 = (lame + 2.0 * shear_modulus) / density;
    const std::complex<double> beta2  = shear_modulus / density;

    // Optimizations
    const double r_inverse    = 1.0 / radius;
    const double r2_inverse   = r_inverse * r_inverse;
    const double r2           = radius * radius;
    const double degree_l_dbl = static_cast<double>(degree_l);
    const double lp1          = degree_l_dbl + 1.0;
    const double dlp1         = 2.0 * degree_l_dbl + 1.0;
    const double llp1         = degree_l_dbl * lp1;

    // Helper functions
    const std::complex<double> k2_quad_pos =(dynamic_term / beta2) + ((dynamic_term + 4.0 * gamma) / alpha2);
    const std::complex<double> k2_quad_neg = (dynamic_term / beta2) - ((dynamic_term + 4.0 * gamma) / alpha2);
    const std::complex<double> k2_quad     = (k2_quad_neg * k2_quad_neg) +
        ((4.0 * degree_l_dbl * (degree_l_dbl + 1.0) * (gamma * gamma)) / (alpha2 * beta2));

    // QUESTION: (Issue #43) KMN15 has these flipped compared to TS72. Going with KMN15 for this func.
    const size_t neg_index = 1;
    const size_t pos_index = 0;
    const std::complex<double> k2_quad_sqrt = std::sqrt(k2_quad);
    const std::complex<double> k2_pos = (1.0 / 2.0) * (k2_quad_pos + k2_quad_sqrt);
    const std::complex<double> k2_neg = (1.0 / 2.0) * (k2_quad_pos - k2_quad_sqrt);

    const std::complex<double> f_k2_pos = (beta2 * k2_pos - dynamic_term) / gamma;
    const std::complex<double> f_k2_neg = (beta2 * k2_neg - dynamic_term) / gamma;

    const std::complex<double> h_k2_pos = f_k2_pos - lp1;
    const std::complex<double> h_k2_neg = f_k2_neg - lp1;

    const std::complex<double> z_k2_pos = c_z_calc(k2_pos * r2, degree_l);
    const std::complex<double> z_k2_neg = c_z_calc(k2_neg * r2, degree_l);

    // See Eqs. B1-B12 of KMN15

    // y1, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 0] =
        -f_k2_pos * z_k2_pos * r_inverse;
    starting_conditions_ptr[neg_index * num_ys + 0] =
        -f_k2_neg * z_k2_neg * r_inverse;
    starting_conditions_ptr[2 * num_ys + 0] =
        degree_l_dbl * r_inverse;

    // y2, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 1] =
        -density * f_k2_pos * alpha2 * k2_pos + (2.0 * shear_modulus * r2_inverse) * (2.0 * f_k2_pos + llp1) * z_k2_pos;
    starting_conditions_ptr[neg_index * num_ys + 1] =
        -density * f_k2_neg * alpha2 * k2_neg + (2.0 * shear_modulus * r2_inverse) * (2.0 * f_k2_neg + llp1) * z_k2_neg;
    starting_conditions_ptr[2 * num_ys + 1] =
        2.0 * shear_modulus * degree_l_dbl * (degree_l_dbl - 1.0) * r2_inverse;

    // y3, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 2] =
        z_k2_pos * r_inverse;
    starting_conditions_ptr[neg_index * num_ys + 2] =
        z_k2_neg * r_inverse;
    starting_conditions_ptr[2 * num_ys + 2] =
        r_inverse;

    // y4, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 3] =
        shear_modulus * k2_pos - (2.0 * shear_modulus * r2_inverse) * (f_k2_pos + 1.0) * z_k2_pos;
    starting_conditions_ptr[neg_index * num_ys + 3] =
        shear_modulus * k2_neg - (2.0 * shear_modulus * r2_inverse) * (f_k2_neg + 1.0) * z_k2_neg;
    starting_conditions_ptr[2 * num_ys + 3] =
        2.0 * shear_modulus * (degree_l_dbl - 1.0) * r2_inverse;

    // y5, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 4] =
        3.0 * gamma * f_k2_pos - h_k2_pos * (degree_l_dbl * gamma - dynamic_term);
    starting_conditions_ptr[neg_index * num_ys + 4] =
        3.0 * gamma * f_k2_neg - h_k2_neg * (degree_l_dbl * gamma - dynamic_term);
    starting_conditions_ptr[2 * num_ys + 4] =
        degree_l_dbl * gamma - dynamic_term;

    // y6, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 5] =
        dlp1 * starting_conditions_ptr[pos_index * num_ys + 4] * r_inverse;
    starting_conditions_ptr[neg_index * num_ys + 5] =
        dlp1 * starting_conditions_ptr[neg_index * num_ys + 4] * r_inverse;
    starting_conditions_ptr[2 * num_ys + 5] =
        dlp1 * starting_conditions_ptr[2 * num_ys + 4] * r_inverse - (3.0 * degree_l_dbl * gamma * r_inverse);
}


// Calculate the starting guess at the bottom of a solid layer using the static assumption.
// KMN15 Eqs. B1-B16 (with w=0).
// Three independent solutions (sn1, sn2, sn3).
inline void c_kamata_solid_static_compressible(
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
    const double r2_inverse   = r_inverse * r_inverse;
    const double r2           = radius * radius;
    const double degree_l_dbl = static_cast<double>(degree_l);
    const double lp1          = degree_l_dbl + 1.0;
    const double dlp1         = 2.0 * degree_l_dbl + 1.0;
    const double llp1         = degree_l_dbl * lp1;

    // Helper functions
    const std::complex<double> k2_quad_pos = 4.0 * gamma / alpha2;
    const std::complex<double> k2_quad_neg = -k2_quad_pos;
    const std::complex<double> k2_quad = k2_quad_neg * k2_quad_neg + ((4.0 * degree_l_dbl * lp1 * gamma * gamma) / (alpha2 * beta2));

    // QUESTION: (Issue #43) KMN15 has these flipped compared to TS72. Going with KMN15 for this func.
    const size_t neg_index = 1;
    const size_t pos_index = 0;
    const std::complex<double> k2_quad_sqrt = std::sqrt(k2_quad);
    const std::complex<double> k2_pos = (1.0 / 2.0) * (k2_quad_pos + k2_quad_sqrt);
    const std::complex<double> k2_neg = (1.0 / 2.0) * (k2_quad_pos - k2_quad_sqrt);

    const std::complex<double> f_k2_pos = beta2 * k2_pos / gamma;
    const std::complex<double> f_k2_neg = beta2 * k2_neg / gamma;

    const std::complex<double> h_k2_pos = f_k2_pos - lp1;
    const std::complex<double> h_k2_neg = f_k2_neg - lp1;

    const std::complex<double> z_k2_pos = c_z_calc(k2_pos * r2, degree_l);
    const std::complex<double> z_k2_neg = c_z_calc(k2_neg * r2, degree_l);

    // See Eqs. B1-B12 of KMN15

    // y1, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 0] =
        -f_k2_pos * z_k2_pos * r_inverse;
    starting_conditions_ptr[neg_index * num_ys + 0] =
        -f_k2_neg * z_k2_neg * r_inverse;
    starting_conditions_ptr[2 * num_ys + 0] =
        degree_l_dbl * r_inverse;

    // y2, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 1] =
        -density * f_k2_pos * alpha2 * k2_pos + (2.0 * shear_modulus * r2_inverse) * (2.0 * f_k2_pos + llp1) * z_k2_pos;
    starting_conditions_ptr[neg_index * num_ys + 1] =
        -density * f_k2_neg * alpha2 * k2_neg + (2.0 * shear_modulus * r2_inverse) * (2.0 * f_k2_neg + llp1) * z_k2_neg;
    starting_conditions_ptr[2 * num_ys + 1] =
        2.0 * shear_modulus * degree_l_dbl * (degree_l_dbl - 1.0) * r2_inverse;

    // y3, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 2] =
        z_k2_pos * r_inverse;
    starting_conditions_ptr[neg_index * num_ys + 2] =
        z_k2_neg * r_inverse;
    starting_conditions_ptr[2 * num_ys + 2] =
        r_inverse;

    // y4, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 3] =
        shear_modulus * k2_pos - (2.0 * shear_modulus * r2_inverse) * (f_k2_pos + 1.0) * z_k2_pos;
    starting_conditions_ptr[neg_index * num_ys + 3] =
        shear_modulus * k2_neg - (2.0 * shear_modulus * r2_inverse) * (f_k2_neg + 1.0) * z_k2_neg;
    starting_conditions_ptr[2 * num_ys + 3] =
        2.0 * shear_modulus * (degree_l_dbl - 1.0) * r2_inverse;

    // y5, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 4] =
        3.0 * gamma * f_k2_pos - h_k2_pos * (degree_l_dbl * gamma);
    starting_conditions_ptr[neg_index * num_ys + 4] =
        3.0 * gamma * f_k2_neg - h_k2_neg * (degree_l_dbl * gamma);
    starting_conditions_ptr[2 * num_ys + 4] =
        degree_l_dbl * gamma;

    // y6, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 5] =
        dlp1 * starting_conditions_ptr[pos_index * num_ys + 4] * r_inverse;
    starting_conditions_ptr[neg_index * num_ys + 5] =
        dlp1 * starting_conditions_ptr[neg_index * num_ys + 4] * r_inverse;
    starting_conditions_ptr[2 * num_ys + 5] =
        dlp1 * starting_conditions_ptr[2 * num_ys + 4] * r_inverse - (3.0 * degree_l_dbl * gamma * r_inverse);
}


// Calculate the starting guess at the bottom of a solid layer using the dynamic and incompressible assumption.
// KMN15 Eqs. B17-B28.
// Three independent solutions (sn1, sn2, sn3).
inline void c_kamata_solid_dynamic_incompressible(
        const double frequency,
        const double radius,
        const double density,
        const std::complex<double>& shear_modulus,
        const int degree_l,
        const double G_to_use,
        const size_t num_ys,
        std::complex<double>* starting_conditions_ptr) noexcept
{
    // Constants
    const double gamma         = 4.0 * TidalPyConstants::d_PI * G_to_use * density / 3.0;
    const double dynamic_term  = frequency * frequency;
    const std::complex<double> beta2 = shear_modulus / density;

    // Optimizations
    const double r_inverse    = 1.0 / radius;
    const double r2_inverse   = r_inverse * r_inverse;
    const double r2           = radius * radius;
    const double degree_l_dbl = static_cast<double>(degree_l);
    const double lp1          = degree_l_dbl + 1.0;
    const double lm1          = degree_l_dbl - 1.0;
    const double dlp1         = 2.0 * degree_l_dbl + 1.0;
    const double llp1         = degree_l_dbl * lp1;

    // QUESTION: (Issue #43) KMN15 has these flipped compared to TS72. Going with KMN15 for this func.
    const size_t neg_index = 1;
    const size_t pos_index = 0;
    const std::complex<double> k2_pos   = dynamic_term / beta2;
    const std::complex<double> f_k2_neg = -dynamic_term / gamma;
    const std::complex<double> h_k2_neg = f_k2_neg - lp1;
    const std::complex<double> z_k2_pos = c_z_calc(k2_pos * r2, degree_l);

    // See Eqs. B17-B28 of KMN15

    // y1, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 0] =
        std::complex<double>(0.0, 0.0);
    starting_conditions_ptr[neg_index * num_ys + 0] =
        std::complex<double>(0.0, 0.0);
    starting_conditions_ptr[2 * num_ys + 0] =
        degree_l_dbl * r_inverse;

    // y2, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 1] =
        llp1 * (-density * gamma + 2.0 * shear_modulus * z_k2_pos * r2_inverse);
    starting_conditions_ptr[neg_index * num_ys + 1] =
        density * ((dynamic_term / gamma) * (dynamic_term + 4.0 * gamma) - llp1 * gamma);
    starting_conditions_ptr[2 * num_ys + 1] =
        2.0 * shear_modulus * degree_l_dbl * lm1 * r2_inverse;

    // y3, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 2] =
        z_k2_pos * r_inverse;
    starting_conditions_ptr[neg_index * num_ys + 2] =
        std::complex<double>(0.0, 0.0);
    starting_conditions_ptr[2 * num_ys + 2] =
        r_inverse;

    // y4, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 3] =
        shear_modulus * (dynamic_term / beta2 - 2.0 * r2_inverse * z_k2_pos);
    starting_conditions_ptr[neg_index * num_ys + 3] =
        std::complex<double>(0.0, 0.0);
    starting_conditions_ptr[2 * num_ys + 3] =
        2.0 * shear_modulus * lm1 * r2_inverse;

    // y5, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 4] =
        lp1 * (degree_l_dbl * gamma - dynamic_term);
    starting_conditions_ptr[neg_index * num_ys + 4] =
        (h_k2_neg - 3.0) * dynamic_term - h_k2_neg * degree_l_dbl * gamma;
    starting_conditions_ptr[2 * num_ys + 4] =
        degree_l_dbl * gamma - dynamic_term;

    // y6, solutions 1--3
    starting_conditions_ptr[pos_index * num_ys + 5] =
        dlp1 * starting_conditions_ptr[pos_index * num_ys + 4] * r_inverse;
    starting_conditions_ptr[neg_index * num_ys + 5] =
        dlp1 * starting_conditions_ptr[neg_index * num_ys + 4] * r_inverse;
    starting_conditions_ptr[2 * num_ys + 5] =
        dlp1 * starting_conditions_ptr[2 * num_ys + 4] * r_inverse - (3.0 * degree_l_dbl * gamma * r_inverse);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Liquid Layers
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Calculate the starting guess at the bottom of a liquid layer using the dynamic assumption.
// KMN15 Eqs. B29-B37.
// Two independent solutions (sn1, sn2).
inline void c_kamata_liquid_dynamic_compressible(
        const double frequency,
        const double radius,
        const double density,
        const std::complex<double>& bulk_modulus,
        const int degree_l,
        const double G_to_use,
        const size_t num_ys,
        std::complex<double>* starting_conditions_ptr) noexcept
{
    // For liquid layer, shear modulus is zero so 1st Lame parameter = bulk modulus
    const std::complex<double> lame = bulk_modulus;

    // Optimizations
    const double dynamic_term = frequency * frequency;
    const double r_inverse    = 1.0 / radius;
    const double r2           = radius * radius;
    const double degree_l_dbl = static_cast<double>(degree_l);

    // Helper functions
    const double gamma          = (4.0 * TidalPyConstants::d_PI * G_to_use * density / 3.0);
    const double f              = -dynamic_term / gamma;
    const double h              = f - (degree_l_dbl + 1.0);
    const std::complex<double> alpha2 = lame / density;
    const std::complex<double> k2     = (1.0 / alpha2) * (dynamic_term + 4.0 * gamma -
        (degree_l_dbl * (degree_l_dbl + 1.0) * gamma * gamma / dynamic_term));

    // See Eqs. B33--B36 in KMN15
    // y1, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 0] =
        -f * r_inverse * c_z_calc(k2 * r2, degree_l);
    starting_conditions_ptr[1 * num_ys + 0] =
        degree_l_dbl * r_inverse;

    // y2, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 1] =
        -density * (f * (dynamic_term + 4.0 * gamma) + degree_l_dbl * (degree_l_dbl + 1.0) * gamma);
    starting_conditions_ptr[1 * num_ys + 1] =
        std::complex<double>(0.0, 0.0);

    // y5, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 2] =
        3.0 * gamma * f - h * (degree_l_dbl * gamma - dynamic_term);
    starting_conditions_ptr[1 * num_ys + 2] =
        degree_l_dbl * gamma - dynamic_term;

    // y6, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 3] =
        (2.0 * degree_l_dbl + 1.0) * starting_conditions_ptr[0 * num_ys + 2] * r_inverse;
    starting_conditions_ptr[1 * num_ys + 3] =
        ((2.0 * degree_l_dbl + 1.0) * starting_conditions_ptr[1 * num_ys + 2] * r_inverse) -
        ((3.0 * degree_l_dbl * gamma) * r_inverse);
}


// Calculate the starting guess at the bottom of a liquid layer using the dynamic and incompressible assumption.
// KMN15 Eqs. B29-B37 (incompressible limit).
// Two independent solutions (sn1, sn2).
inline void c_kamata_liquid_dynamic_incompressible(
        const double frequency,
        const double radius,
        const double density,
        const int degree_l,
        const double G_to_use,
        const size_t num_ys,
        std::complex<double>* starting_conditions_ptr) noexcept
{
    // Optimizations
    const double dynamic_term = frequency * frequency;
    const double r_inverse    = 1.0 / radius;
    const double degree_l_dbl = static_cast<double>(degree_l);
    const double lp1          = degree_l_dbl + 1.0;
    const double llp1         = degree_l_dbl * lp1;
    const double dlp1         = 2.0 * degree_l_dbl + 1.0;

    // Helper functions
    const double gamma = (4.0 * TidalPyConstants::d_PI * G_to_use * density / 3.0);
    const double f     = -dynamic_term / gamma;
    const double h     = f - lp1;

    // See Eqs. B33--B36 in KMN15

    // y1, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 0] =
        std::complex<double>(0.0, 0.0);
    starting_conditions_ptr[1 * num_ys + 0] =
        degree_l_dbl * r_inverse;

    // y2, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 1] =
        -density * (f * (dynamic_term + 4.0 * gamma) + llp1 * gamma);
    starting_conditions_ptr[1 * num_ys + 1] =
        std::complex<double>(0.0, 0.0);

    // y5, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 2] =
        3.0 * gamma * f - h * (degree_l_dbl * gamma - dynamic_term);
    starting_conditions_ptr[1 * num_ys + 2] =
        degree_l_dbl * gamma - dynamic_term;

    // y6, solutions 1--2
    starting_conditions_ptr[0 * num_ys + 3] =
        dlp1 * starting_conditions_ptr[0 * num_ys + 2] * r_inverse;
    starting_conditions_ptr[1 * num_ys + 3] =
        (dlp1 * starting_conditions_ptr[1 * num_ys + 2] * r_inverse) - ((3.0 * degree_l_dbl * gamma) * r_inverse);
}
