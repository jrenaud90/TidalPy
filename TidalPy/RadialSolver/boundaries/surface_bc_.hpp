// surface_bc_.hpp - Surface boundary conditions for the radial solver.
//
// References
// ----------
// Beuthe (2015)
// Saito (1974)
#pragma once

#include <cstddef>
#include <limits>
#include <string>

#include "../rs_constants_.hpp"   // C_MAX_NUM_YTYPES, C_MAX_SURFACE_CONDITIONS


// Populates a boundary condition array at a planet's surface.
//
// Parameters
// ----------
// boundary_conditions_ptr : double*
//     Output array, must have space for C_MAX_SURFACE_CONDITIONS elements (3 per model, at most C_MAX_NUM_YTYPES).
// bc_model_ptr : int*
//     Array of BC model types: 0=free surface, 1=tidal potential, 2=loading potential.
// num_bcs : size_t
//     Number of BC models, 1 to C_MAX_NUM_YTYPES.
// radius_to_use : double
//     Planet radius [m].
// bulk_density_to_use : double
//     Planet bulk density [kg m-3].
// degree_l_dbl : double
//     Tidal harmonic degree.
//
// Returns
// -------
// int : 0=success, -1=too many models, -2=num_bcs<=0, -3=unknown model (c_surface_bc_error_message says which).
inline int c_get_surface_bc(
        double* boundary_conditions_ptr,
        const int* bc_model_ptr,
        size_t num_bcs,
        double radius_to_use,
        double bulk_density_to_use,
        double degree_l_dbl) noexcept
{
    const double nan_val = std::numeric_limits<double>::quiet_NaN();

    if (num_bcs > C_MAX_NUM_YTYPES) [[unlikely]]
    {
        return -1;
    }
    if (num_bcs <= 0) [[unlikely]]
    {
        return -2;
    }

    for (size_t i = 0; i < C_MAX_SURFACE_CONDITIONS; ++i) {
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


// Why c_get_surface_bc refused a model list (its non-zero code), for a solve's message.
inline std::string c_surface_bc_error_message(int code)
{
    return "Invalid surface boundary conditions (code " + std::to_string(code) + "): between 1 and " +
        std::to_string(C_MAX_NUM_YTYPES) + " models are allowed, each free (0), tidal (1), or loading (2).\n";
}
