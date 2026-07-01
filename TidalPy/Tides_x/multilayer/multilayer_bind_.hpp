// multilayer_bind_.hpp - Flat-pointer wrappers around the on-demand 3D kernel.
//
// The kernel/potential headers use std::complex and small structs; these wrappers expose them with
// plain double pointers so the Cython layer can bind them directly. Complex values are passed and
// returned as adjacent (real, imag) double pairs.
#pragma once

#include <complex>
#include <cstddef>

#include "kernel_.hpp"


namespace tidalpy {
namespace tides {

// Volumetric tidal heating [W m-3] from the 6 complex stress and 6 complex strain components, each passed
// as 12 doubles (re, im per component). Used to heat once from the mode-summed stress/strain tensors.
inline double c_volumetric_heating_flat(const double* stress12, const double* strain12) noexcept
{
    c_Tensor6 stress, strain;
    for (std::size_t k = 0; k < 6; ++k)
    {
        stress.c[k] = std::complex<double>(stress12[2 * k], stress12[2 * k + 1]);
        strain.c[k] = std::complex<double>(strain12[2 * k], strain12[2 * k + 1]);
    }
    return c_volumetric_heating(stress, strain);
}

// Strain/stress/heating at one point. y_ri = 12 doubles (y1re,y1im,...,y6re,y6im; only y1..y4 used).
// pot6 = the 6 potential values above. strain12/stress12 = 12 doubles each (6 complex). heating1 = 1 double.
inline void c_strain_stress_heating(
        const double* y_ri,
        double shear_re, double shear_im, double bulk_re, double bulk_im,
        double radius, double degree_l, int is_solid, int is_incomp,
        const double* pot6, double colatitude,
        double* strain12, double* stress12, double* heating1) noexcept
{
    const std::complex<double> y1(y_ri[0], y_ri[1]);
    const std::complex<double> y2(y_ri[2], y_ri[3]);
    const std::complex<double> y3(y_ri[4], y_ri[5]);
    const std::complex<double> y4(y_ri[6], y_ri[7]);
    const std::complex<double> shear(shear_re, shear_im);
    const std::complex<double> bulk(bulk_re, bulk_im);

    const c_StrainRadialCoeffs R = c_compute_strain_radial_coeffs(
        y1, y2, y3, y4, shear, bulk, radius, degree_l, is_solid != 0, is_incomp != 0
    );

    c_PotentialPoint P{pot6[0], pot6[1], pot6[2], pot6[3], pot6[4], pot6[5]};

    c_Tensor6 strain, stress;
    c_compute_strain_stress(R, P, colatitude, strain, stress);
    for (std::size_t k = 0; k < 6; ++k)
    {
        strain12[2 * k] = strain.c[k].real(); strain12[2 * k + 1] = strain.c[k].imag();
        stress12[2 * k] = stress.c[k].real(); stress12[2 * k + 1] = stress.c[k].imag();
    }
    *heating1 = c_volumetric_heating(stress, strain);
}

}  // namespace tides
}  // namespace tidalpy
