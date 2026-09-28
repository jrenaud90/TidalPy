// multilayer_bind_.hpp - Flat-pointer wrappers around the on-demand 3D kernel.
//
// The kernel/potential headers use std::complex and small structs; these wrappers expose them with
// plain double pointers so the Cython layer can bind them directly. Complex values are passed and
// returned as adjacent (real, imag) double pairs.
#pragma once

#include <complex>
#include <cstddef>

#include "kernel_.hpp"
#include "angular_gram_.hpp"


namespace tidalpy {
namespace tides {

// Fill the symmetric 6x6 angular Gram matrix for (l, m) into gram36 (row-major). Returns 1 on success,
// 0 if (l, m) is out of the tabulated range (l = 2..10, m = 0..l). Used to test the table from Python.
inline int c_angular_gram_flat(int degree_l, int order_m, double* gram36) noexcept
{
    double gram[6][6];
    if (!c_angular_gram(degree_l, order_m, gram)) { return 0; }
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < 6; ++j) { gram36[i * 6 + j] = gram[i][j]; }
    }
    return 1;
}

// Cycle-averaged volumetric heating [W m-3] at a forcing frequency [rad s-1] from the 6 complex stress and 6 complex
// strain amplitudes at that frequency, each passed as 12 doubles (re, im per component): (|omega|/2) times the
// magnitude of the weighted Im(stress conj strain), the same factor the world path applies to its summed amplitudes.
inline double c_volumetric_heating_flat(const double* stress12, const double* strain12, double frequency) noexcept
{
    c_Tensor6 stress, strain;
    for (std::size_t k = 0; k < 6; ++k)
    {
        stress.c[k] = std::complex<double>(stress12[2 * k], stress12[2 * k + 1]);
        strain.c[k] = std::complex<double>(strain12[2 * k], strain12[2 * k + 1]);
    }
    return 0.5 * std::abs(frequency) * c_volumetric_heating(stress, strain);
}

// Complex potential point from 12 doubles: (real, imaginary) of U, dU/dtheta, dU/dphi, d2U/dtheta2, d2U/dphi2, and
// d2U/dtheta_dphi, in the phasor convention of the 3D potential engine (U(t) = Re[U_c e^{i omega t}]).
inline c_PotentialPointC c_potential_point_from_flat(const double* potential12) noexcept
{
    return c_PotentialPointC(
        std::complex<double>(potential12[0], potential12[1]),
        std::complex<double>(potential12[2], potential12[3]),
        std::complex<double>(potential12[4], potential12[5]),
        std::complex<double>(potential12[6], potential12[7]),
        std::complex<double>(potential12[8], potential12[9]),
        std::complex<double>(potential12[10], potential12[11]));
}

// Strain/stress/heating at one point. y_ri = 12 doubles (y1re,y1im,...,y6re,y6im; only y1..y4 used).
// potential12 = one mode's complex potential row as 12 doubles (see c_potential_point_from_flat).
// strain12/stress12 = 12 doubles each (6 complex amplitudes). heating1 = 1 double, this mode's cycle-averaged
// volumetric heating [W m-3] at the forcing frequency [rad s-1] (see c_volumetric_heating_flat).
inline void c_strain_stress_heating(
        const double* y_ri,
        double shear_re,
        double shear_im,
        double bulk_re,
        double bulk_im,
        double radius,
        double degree_l,
        double frequency,
        int is_solid,
        int is_incomp,
        const double* potential12,
        double colatitude,
        double* strain12,
        double* stress12,
        double* heating1) noexcept
{
    const std::complex<double> y1(y_ri[0], y_ri[1]);
    const std::complex<double> y2(y_ri[2], y_ri[3]);
    const std::complex<double> y3(y_ri[4], y_ri[5]);
    const std::complex<double> y4(y_ri[6], y_ri[7]);
    const std::complex<double> shear(shear_re, shear_im);
    const std::complex<double> bulk(bulk_re, bulk_im);

    const c_StrainRadialCoeffs R = c_compute_strain_radial_coeffs(
        y1,
        y2,
        y3,
        y4,
        shear,
        bulk,
        radius,
        degree_l,
        is_solid != 0,
        is_incomp != 0);

    const c_PotentialPointC P = c_potential_point_from_flat(potential12);

    c_Tensor6 strain, stress;
    c_compute_strain_stress(R, P, colatitude, strain, stress);
    for (std::size_t k = 0; k < 6; ++k)
    {
        strain12[2 * k] = strain.c[k].real(); strain12[2 * k + 1] = strain.c[k].imag();
        stress12[2 * k] = stress.c[k].real(); stress12[2 * k + 1] = stress.c[k].imag();
    }
    *heating1 = 0.5 * std::abs(frequency) * c_volumetric_heating(stress, strain);
}

// Displacements at one point. y_ri = 12 doubles (y1re, y1im, ..., y6re, y6im; only y1 and y3 used), potential12 =
// one mode's complex potential row as 12 doubles, disp6 = 6 doubles (3 complex amplitudes: radial, polar,
// azimuthal) [m].
inline void c_displacements_flat(
        const double* y_ri,
        const double* potential12,
        double colatitude,
        double* disp6) noexcept
{
    const std::complex<double> y1(y_ri[0], y_ri[1]);
    const std::complex<double> y3(y_ri[4], y_ri[5]);
    const c_PotentialPointC P = c_potential_point_from_flat(potential12);
    c_Vector3 u;
    c_compute_displacements(y1, y3, P, colatitude, u);
    for (std::size_t k = 0; k < 3; ++k)
    {
        disp6[2 * k] = u.c[k].real(); disp6[2 * k + 1] = u.c[k].imag();
    }
}

}  // namespace tides
}  // namespace tidalpy
