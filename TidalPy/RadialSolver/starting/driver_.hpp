// driver_.hpp: dispatcher for the shooting method's starting conditions.
#pragma once

#include <cmath>
#include <complex>
#include <string>

#include "kamata_.hpp"
#include "takeuchi_.hpp"
#include "saito_.hpp"
#include "power_series_.hpp"
#include "unity_.hpp"
#include "starting_method_.hpp"
#include "../layer_kind_.hpp"


// Fill starting_conditions_ptr (num_ys per solution) for the layer type and assumptions at radius [m], by
// starting_method (c_StartingMethod as an int): Takeuchi and Saito (1972), Kamata et al. (2015), the Martens (2016)
// power series, or unit vectors. Every method covers every layer type; static liquids use Saito (1974) for every
// method but unity.
// Sets *success_ptr and message for an unknown method, a power series that refused to start, or, when run_y_checks,
// a wrong num_ys.
// layer_type: 0 = solid, 1 = liquid. Units: density [kg m-3], moduli [Pa], frequency [rad s-1].
inline void c_find_starting_conditions(
        bool* success_ptr,
        std::string& message,
        const int layer_type,
        const bool is_static,
        const bool is_incompressible,
        const int starting_method,
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
    const c_StartingMethod method = static_cast<c_StartingMethod>(starting_method);
    if ((starting_method < static_cast<int>(c_StartingMethod::Takeuchi)) ||
        (starting_method > static_cast<int>(c_StartingMethod::Unity)))
    {
        *success_ptr = false;
        message = "RadialSolver::Shooting::FindStartingConditions: Unknown starting method (" +
            std::to_string(starting_method) + ").";
        return;
    }

    if (run_y_checks)
    {
        if (c_layer_num_ys(layer_type, is_static) != num_ys)
        {
            *success_ptr = false;
            message = "RadialSolver::Shooting::FindStartingConditions: Incorrect number of ys for given the starting "
                "condition assumptions.";
            return;
        }
    }

    if (method == c_StartingMethod::Unity)
    {
        c_unity_starting_conditions(layer_type, is_static, num_ys, starting_conditions_ptr);
    }
    else if (is_liquid && is_static)
    {
        c_saito_liquid_static_incompressible(
            radius, degree_l, num_ys, starting_conditions_ptr
            );
    }
    else if (method == c_StartingMethod::PowerSeries)
    {
        bool converged;
        if (!is_liquid)
        {
            converged = c_power_series_solid(
                is_static ? 0.0 : frequency,
                radius,
                density,
                bulk_modulus,
                shear_modulus,
                is_incompressible,
                degree_l,
                G_to_use,
                num_ys,
                starting_conditions_ptr);
        }
        else
        {
            converged = c_power_series_liquid_dynamic(
                frequency,
                radius,
                density,
                bulk_modulus,
                is_incompressible,
                degree_l,
                G_to_use,
                num_ys,
                starting_conditions_ptr);
        }
        if (!converged)
        {
            *success_ptr = false;
            message = "RadialSolver::Shooting::FindStartingConditions: The power series starting conditions refused "
                "to start: the series did not converge at the starting radius, or the solutions grow too steeply "
                "there (a weak solid, or a dynamic liquid at long periods).\n"
                "Recommend starting_method='takeuchi' or 'kamata', or a smaller starting radius.";
            return;
        }
    }
    else if (method == c_StartingMethod::Kamata)
    {
        if (!is_liquid)
        {
            if (is_static && is_incompressible)
            {
                c_kamata_solid_static_incompressible(
                    radius,
                    density,
                    shear_modulus,
                    degree_l,
                    G_to_use,
                    num_ys,
                    starting_conditions_ptr);
            }
            else if (is_static)
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
        if (is_static && is_incompressible)
        {
            c_takeuchi_solid_static_incompressible(
                radius,
                density,
                shear_modulus,
                degree_l,
                G_to_use,
                num_ys,
                starting_conditions_ptr);
        }
        else if (is_static)
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
        else if (is_incompressible)
        {
            c_takeuchi_solid_dynamic_incompressible(
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
    else if (is_incompressible)
    {
        c_takeuchi_liquid_dynamic_incompressible(
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
