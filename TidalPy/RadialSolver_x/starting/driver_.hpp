// driver_.hpp: dispatcher for the shooting method's starting conditions.
#pragma once

#include <cmath>
#include <complex>
#include <string>

#include "kamata_.hpp"
#include "takeuchi_.hpp"
#include "saito_.hpp"


// Fill starting_conditions_ptr (num_ys per solution) for the layer type and assumptions at radius [m], from
// Kamata et al. (2015) when use_kamata, else Takeuchi and Saito (1972); static liquids always use Saito (1974).
// Sets *success_ptr and message on an unsupported combination or, when run_y_checks, a wrong num_ys.
// layer_type: 0 = solid, 1 = liquid. Units: density [kg m-3], moduli [Pa], frequency [rad s-1].
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
    *success_ptr = true;

    const bool is_liquid = (layer_type != 0);

    // Saito (1974) covers every static liquid, so only the other layers can ask for a combination with no
    // starting conditions.
    if (!(is_liquid && is_static))
    {
        if (use_kamata && !is_liquid && is_static && is_incompressible)
        {
            *success_ptr = false;
            message = "RadialSolver::Shooting::FindStartingConditions: Incompressibility is not implemented for "
                "Kamata starting conditions for static-solid layers.\nRecommend using dynamic-incompressible "
                "instead.";
            return;
        }
        if (!use_kamata && is_incompressible)
        {
            *success_ptr = false;
            message = "RadialSolver::Shooting::FindStartingConditions: Incompressibility is not implemented for "
                "most of the Takeuchi starting conditions. \nRecommend using Kamata (set use_kamata=True) "
                "instead.";
            return;
        }
    }

    if (run_y_checks)
    {
        // Two ys for a static liquid, four for a dynamic one, six for a solid.
        const size_t num_ys_for_assumption = is_liquid ? (is_static ? 2 : 4) : 6;
        if (num_ys_for_assumption != num_ys)
        {
            *success_ptr = false;
            message = "RadialSolver::Shooting::FindStartingConditions: Incorrect number of ys for given the starting "
                "condition assumptions.";
            return;
        }
    }

    if (is_liquid && is_static)
    {
        c_saito_liquid_static_incompressible(
            radius, degree_l, num_ys, starting_conditions_ptr
            );
    }
    else if (use_kamata)
    {
        if (!is_liquid)
        {
            if (is_static)
            {
                c_kamata_solid_static_compressible(
                    radius,
                    density,
                    bulk_modulus,
                    shear_modulus,
                    degree_l,
                    G_to_use,
                    num_ys,
                    starting_conditions_ptr);
            }
            else if (is_incompressible)
            {
                c_kamata_solid_dynamic_incompressible(
                    frequency,
                    radius,
                    density,
                    shear_modulus,
                    degree_l,
                    G_to_use,
                    num_ys,
                    starting_conditions_ptr);
            }
            else
            {
                c_kamata_solid_dynamic_compressible(
                    frequency,
                    radius,
                    density,
                    bulk_modulus,
                    shear_modulus,
                    degree_l,
                    G_to_use,
                    num_ys,
                    starting_conditions_ptr);
            }
        }
        else if (is_incompressible)
        {
            c_kamata_liquid_dynamic_incompressible(
                frequency,
                radius,
                density,
                degree_l,
                G_to_use,
                num_ys,
                starting_conditions_ptr);
        }
        else
        {
            c_kamata_liquid_dynamic_compressible(
                frequency,
                radius,
                density,
                bulk_modulus,
                degree_l,
                G_to_use,
                num_ys,
                starting_conditions_ptr);
        }
    }
    else if (!is_liquid)
    {
        if (is_static)
        {
            c_takeuchi_solid_static_compressible(
                radius,
                density,
                bulk_modulus,
                shear_modulus,
                degree_l,
                G_to_use,
                num_ys,
                starting_conditions_ptr);
        }
        else
        {
            c_takeuchi_solid_dynamic_compressible(
                frequency,
                radius,
                density,
                bulk_modulus,
                shear_modulus,
                degree_l,
                G_to_use,
                num_ys,
                starting_conditions_ptr);
        }
    }
    else
    {
        c_takeuchi_liquid_dynamic_compressible(
            frequency,
            radius,
            density,
            bulk_modulus,
            degree_l,
            G_to_use,
            num_ys,
            starting_conditions_ptr);
    }
}
