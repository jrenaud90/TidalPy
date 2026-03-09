// saito_.hpp - Starting conditions based on Saito (1974)
// Ported from TidalPy/RadialSolver/starting/saito.pyx
//
// References
// ----------
// S74: Saito (1974) Eq. 19
#pragma once

#include <cmath>
#include <complex>


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Liquid Layers
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Calculate the initial guess at the bottom of a liquid layer using the static assumption.
// S74 Eq. 19.
// One independent solution (sn1).
//
// Using the static assumption in a liquid layer results in one independent solution for the radial derivative.
// Allows for a general tidal harmonic l, for static tides (w = 0).
// Compressibility and all dissipation dependence is lost due to no dependence on bulk or shear moduli.
inline void c_saito_liquid_static_incompressible(
        double radius,
        int degree_l,
        size_t num_ys,
        std::complex<double>* starting_conditions_ptr) noexcept
{
    double degree_l_dbl = static_cast<double>(degree_l);

    // See Eq. 19 in Saito 1974
    // y5 solution 0
    starting_conditions_ptr[0] = std::pow(radius, degree_l_dbl);

    // y7 solution 0
    starting_conditions_ptr[1] = 2.0 * (degree_l_dbl - 1.0) * std::pow(radius, degree_l_dbl - 1.0);
}
