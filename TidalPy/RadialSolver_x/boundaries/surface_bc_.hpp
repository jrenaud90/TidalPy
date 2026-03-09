// surface_bc_.hpp - Surface boundary conditions
// Ported from TidalPy/RadialSolver/boundaries/surface_bc.pyx
//
// References
// ----------
// Beuthe (2015)
// Saito (1974)
#pragma once

#include <cstddef>
#include <limits>


// Populates a boundary condition array at a planet's surface.
//
// Parameters
// ----------
// boundary_conditions_ptr : double*
//     Output array, must have space for 15 elements (5 max models * 3 BCs per model).
// bc_model_ptr : int*
//     Array of BC model types: 0=free surface, 1=tidal potential, 2=loading potential.
// num_bcs : size_t
//     Number of BC models. Must be 1-5.
// radius_to_use : double
//     Planet radius [m].
// bulk_density_to_use : double
//     Planet bulk density [kg m-3].
// degree_l_dbl : double
//     Tidal harmonic degree.
//
// Returns
// -------
// int : 0=success, -1=num_bcs>5, -2=num_bcs<=0, -3=unknown model.
inline int c_get_surface_bc(
        double* boundary_conditions_ptr,
        int* bc_model_ptr,
        size_t num_bcs,
        double radius_to_use,
        double bulk_density_to_use,
        double degree_l_dbl) noexcept
{
    const double nan_val = std::numeric_limits<double>::quiet_NaN();

    if (num_bcs > 5) [[unlikely]]
    {
        return -1;
    }
    if (num_bcs <= 0) [[unlikely]]
    {
        return -2;
    }

    // Initialize all boundary conditions to NaN
    // 15 = 5 (max_num_solutions) * 3 (number of surface conditions)
    for (size_t i = 0; i < 15; ++i) {
        boundary_conditions_ptr[i] = nan_val;
    }

    for (size_t j = 0; j < num_bcs; ++j) {
        if (bc_model_ptr[j] == 0) {
            // Free Surface
            boundary_conditions_ptr[j * 3 + 0] = 0.0;
            boundary_conditions_ptr[j * 3 + 1] = 0.0;
            boundary_conditions_ptr[j * 3 + 2] = 0.0;
        } else if (bc_model_ptr[j] == 1) {
            // Tidal Potential
            boundary_conditions_ptr[j * 3 + 0] = 0.0;
            boundary_conditions_ptr[j * 3 + 1] = 0.0;
            boundary_conditions_ptr[j * 3 + 2] = (2.0 * degree_l_dbl + 1.0) / radius_to_use;
        } else if (bc_model_ptr[j] == 2) {
            // Loading Potential
            // See Eq. 6 in Beuthe (2015) and Eq. 9 of Saito (1974)
            boundary_conditions_ptr[j * 3 + 0] = (-1.0 / 3.0) * (2.0 * degree_l_dbl + 1.0) * bulk_density_to_use;
            boundary_conditions_ptr[j * 3 + 1] = 0.0;
            boundary_conditions_ptr[j * 3 + 2] = (2.0 * degree_l_dbl + 1.0) / radius_to_use;
        } else {
            return -3;
        }
    }
    return 0;
}
