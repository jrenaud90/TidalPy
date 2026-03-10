// driver_.hpp - Dispatcher for starting conditions
// Ported from TidalPy/RadialSolver/starting/driver.pyx
#pragma once

#include <cmath>
#include <complex>
#include <string>

#include "kamata_.hpp"
#include "takeuchi_.hpp"
#include "saito_.hpp"


// Dispatch to the correct starting condition function based on layer type, static/dynamic, and compressibility.
//
// Parameters
// ----------
// success_ptr : bool*, output
//     Set to true on success, false on failure.
// message : std::string&, output
//     Error message on failure.
// layer_type : int
//     0 = solid, 1 = liquid.
// is_static : int
//     1 = static, 0 = dynamic.
// is_incompressible : int
//     1 = incompressible, 0 = compressible.
// use_kamata : bool
//     true = use Kamata (2015) starting conditions, false = use Takeuchi & Saito (1972).
// frequency : double
//     Forcing frequency [rad s-1]. Only used for dynamic cases.
// radius : double
//     Radius where the radial functions are calculated [m].
// density : double
//     Density at radius [kg m-3].
// bulk_modulus : complex
//     Bulk modulus at radius [Pa].
// shear_modulus : complex
//     Shear modulus at radius [Pa].
// degree_l : int
//     Tidal harmonic order.
// G_to_use : double
//     Gravitational constant.
// num_ys : size_t
//     Number of y-values per solution.
// starting_conditions_ptr : complex*, output
//     Output array for starting conditions.
// run_y_checks : bool
//     If true, validate num_ys matches expected value for the chosen method.
inline void c_find_starting_conditions(
        bool* success_ptr,
        std::string& message,
        const int layer_type,
        const bool is_static,
        const bool is_incompressible,
        const bool use_kamata,
        const double frequency,
        const double radius,
        const double density,
        const std::complex<double> bulk_modulus,
        const std::complex<double> shear_modulus,
        const int degree_l,
        const double G_to_use,
        const size_t num_ys,
        std::complex<double>* starting_conditions_ptr,
        const bool run_y_checks = true) noexcept
{
    size_t num_ys_for_assumption;

    // Assume success and adjust if not
    *success_ptr = true;

    // For static liquid layers, no matter the other assumptions, we use Saito's method.
    if ((layer_type != 0) && is_static)
    {
        // Liquid Static Layer
        if (run_y_checks)
        {
            num_ys_for_assumption = 2;
            if (num_ys_for_assumption != num_ys)
            {
                *success_ptr = false;
                message = "RadialSolver::Shooting::FindStartingConditions: Incorrect number of ys for given the starting condition assumptions.";
            }
        }
        if (*success_ptr)
        {
            c_saito_liquid_static_incompressible(
                radius, degree_l, num_ys, starting_conditions_ptr
                );
        }

    // Work through the Kamata models
    } else if (use_kamata)
    {
        // Kamata solid layer
        if (layer_type == 0) {
            // Solid layer
            if (is_static && is_incompressible)
            {
                *success_ptr = false;
                message = "RadialSolver::Shooting::FindStartingConditions: Incompressibility is not implemented for Kamata starting conditions for static-solid layers.\nRecommend using dynamic-incompressible instead.";
            } else if (is_static && (!is_incompressible))
            {
                if (run_y_checks)
                {
                    num_ys_for_assumption = 6;
                    if (num_ys_for_assumption != num_ys)
                    {
                        *success_ptr = false;
                        message = "RadialSolver::Shooting::FindStartingConditions: Incorrect number of ys for given the starting condition assumptions.";
                    }
                }
                if (*success_ptr)
                {
                    c_kamata_solid_static_compressible(
                        radius, density, bulk_modulus, shear_modulus, degree_l, G_to_use, num_ys, starting_conditions_ptr
                    );
                }
            } else if ((!is_static) && is_incompressible)
            {
                if (run_y_checks)
                {
                    num_ys_for_assumption = 6;
                    if (num_ys_for_assumption != num_ys)
                    {
                        *success_ptr = false;
                        message = "RadialSolver::Shooting::FindStartingConditions: Incorrect number of ys for given the starting condition assumptions.";
                    }
                }
                if (*success_ptr)
                {
                    c_kamata_solid_dynamic_incompressible(
                        frequency, radius, density, shear_modulus, degree_l, G_to_use, num_ys, starting_conditions_ptr
                        );
                }
            } else
            {
                if (run_y_checks)
                {
                    num_ys_for_assumption = 6;
                    if (num_ys_for_assumption != num_ys)
                    {
                        *success_ptr = false;
                        message = "RadialSolver::Shooting::FindStartingConditions: Incorrect number of ys for given the starting condition assumptions.";
                    }
                }
                if (*success_ptr)
                {
                    c_kamata_solid_dynamic_compressible(
                        frequency, radius, density, bulk_modulus, shear_modulus, degree_l, G_to_use, num_ys,
                        starting_conditions_ptr
                        );
                }
            }
        } else
        {
            // Kamata liquid layer
            if (is_static)
            {
                // Covered by Saito method above
            } else if ((!is_static) && is_incompressible)
            {
                if (run_y_checks)
                {
                    num_ys_for_assumption = 4;
                    if (num_ys_for_assumption != num_ys)
                    {
                        *success_ptr = false;
                        message = "RadialSolver::Shooting::FindStartingConditions: Incorrect number of ys for given the starting condition assumptions.";
                    }
                }
                if (*success_ptr)
                {
                    c_kamata_liquid_dynamic_incompressible(
                        frequency, radius, density, degree_l, G_to_use, num_ys, starting_conditions_ptr
                        );
                }
            } else
            {
                if (run_y_checks)
                {
                    num_ys_for_assumption = 4;
                    if (num_ys_for_assumption != num_ys)
                    {
                        *success_ptr = false;
                        message = "RadialSolver::Shooting::FindStartingConditions: Incorrect number of ys for given the starting condition assumptions.";
                    }
                }
                if (*success_ptr)
                {
                    c_kamata_liquid_dynamic_compressible(
                        frequency, radius, density, bulk_modulus, degree_l, G_to_use, num_ys, starting_conditions_ptr
                        );
                }
            }
        }

    // Work through the Takeuchi models
    } else
    {
        if (is_incompressible)
        {
            *success_ptr = false;
            message = "RadialSolver::Shooting::FindStartingConditions: Incompressibility is not implemented for most of the Takeuchi starting conditions. \nRecommend using Kamata (set use_kamata=True) instead.";
        } else {
            if (layer_type == 0)
            {
                // Solid layer
                if (is_static) {
                    if (run_y_checks) {
                        num_ys_for_assumption = 6;
                        if (num_ys_for_assumption != num_ys)
                        {
                            *success_ptr = false;
                            message = "RadialSolver::Shooting::FindStartingConditions: Incorrect number of ys for given the starting condition assumptions.";
                        }
                    }
                    if (*success_ptr)
                    {
                        c_takeuchi_solid_static_compressible(
                            radius, density, bulk_modulus, shear_modulus, degree_l, G_to_use, num_ys, starting_conditions_ptr
                        );
                    }
                } else
                {
                    if (run_y_checks)
                    {
                        num_ys_for_assumption = 6;
                        if (num_ys_for_assumption != num_ys)
                        {
                            *success_ptr = false;
                            message = "RadialSolver::Shooting::FindStartingConditions: Incorrect number of ys for given the starting condition assumptions.";
                        }
                    }
                    if (*success_ptr)
                    {
                        c_takeuchi_solid_dynamic_compressible(
                            frequency, radius, density, bulk_modulus, shear_modulus,
                            degree_l, G_to_use, num_ys, starting_conditions_ptr
                            );
                    }
                }
            } else
            {
                // Liquid layer
                if (is_static)
                {
                    // Handled by Saito above
                } else
                {
                    if (run_y_checks)
                    {
                        num_ys_for_assumption = 4;
                        if (num_ys_for_assumption != num_ys)
                        {
                            *success_ptr = false;
                            message = "RadialSolver::Shooting::FindStartingConditions: Incorrect number of ys for given the starting condition assumptions.";
                        }
                    }
                    if (*success_ptr)
                    {
                        c_takeuchi_liquid_dynamic_compressible(
                            frequency, radius, density, bulk_modulus, degree_l, G_to_use, num_ys, starting_conditions_ptr
                            );
                    }
                }
            }
        }
    }
}
