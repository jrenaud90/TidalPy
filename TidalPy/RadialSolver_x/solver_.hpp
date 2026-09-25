// solver_.hpp: input checks for the standalone array-based radial solver (solver.pyx).
//
// c_validate_and_prep_radial_inputs validates the supplied arrays and maps the layer-type and boundary-condition
// names to integers; c_parse_ode_method maps an integration method name to CyRK's enum. The solve itself runs on a
// world built from the profile, the same path a built world takes.
#pragma once

#include <cctype>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// CyRK imports
#include "cysolver.hpp"   // ODEMethod

// TidalPy imports
#include "../constants_.hpp"
#include "../Utilities_x/math_x/numerics_.hpp"

// RadialSolver imports
#include "rs_constants_.hpp"

#include "../Utilities_x/classes_x/model_names_.hpp"   // c_to_lower

// CyRK's enum for an integration method name, in any case. `method_label` names the argument in the error an
// unknown name raises.
inline ODEMethod c_parse_ode_method(const std::string& method_name, const std::string& method_label)
{
    const std::string method_lower = tidalpy::c_to_lower(method_name);
    if (method_lower == "rk45")   return ODEMethod::RK45;
    if (method_lower == "rk23")   return ODEMethod::RK23;
    if (method_lower == "dop853") return ODEMethod::DOP853;
    if (method_lower == "bdf")    return ODEMethod::BDF;
    if (method_lower == "lsoda")  return ODEMethod::LSODA;
    if (method_lower == "radau")  return ODEMethod::RADAU;
    throw std::invalid_argument(
        "Unsupported " + method_label + " provided: " + method_lower +
        ". Supported: rk23, rk45, dop853, bdf, lsoda, radau.");
}

inline void c_validate_and_prep_radial_inputs(
    size_t total_slices,
    const double* radius_array,
    double frequency,
    size_t num_layers,
    const std::vector<std::string>& layer_types,
    const bool* is_static_bylayer,
    const bool* is_incompressible_bylayer,
    const double* upper_radius_bylayer_array,
    bool use_prop_matrix,
    double starting_radius,
    const std::vector<std::string>& solve_for,
    const std::vector<std::string>& eos_method_bylayer,
    bool warnings,
    int* layer_types_out_ptr,
    int* bc_models_out_ptr,
    size_t& num_bc_models_out
)
{
    if (num_layers == 0)
        throw std::invalid_argument("At least one layer is required.");
    if (total_slices < 2)
        throw std::invalid_argument("`radius_array` needs at least two values.");
    if (layer_types.size() != num_layers)
        throw std::invalid_argument("Number of `layer_types` must match `num_layers`.");
    if (!eos_method_bylayer.empty() && (eos_method_bylayer.size() != num_layers))
        throw std::invalid_argument("`eos_method_bylayer` must give one method per layer.");
    // The solver builds the planet out to the top layer's upper radius, so a longer profile would be cut silently.
    if (!c_isclose(upper_radius_bylayer_array[num_layers - 1], radius_array[total_slices - 1], 1.0e-9, 0.0))
        throw std::invalid_argument(
            "The top layer's upper radius must equal the last value of `radius_array` (the planet radius).");

    double last_layer_r = 0.0;
    for (size_t layer_i = 0; layer_i < num_layers; ++layer_i)
    {
        if (upper_radius_bylayer_array[layer_i] <= last_layer_r)
        {
            throw std::invalid_argument("`upper_radius_bylayer_array` must be in ascending order.");
        }
        last_layer_r = upper_radius_bylayer_array[layer_i];
    }

    if (std::abs(frequency) < tidalpy_config_ptr->d_MIN_FREQUENCY)
        throw std::invalid_argument("Forcing frequency is too small (are you sure you are in rad s-1?).");
    else if (std::abs(frequency) > tidalpy_config_ptr->d_MAX_FREQUENCY)
        throw std::invalid_argument("Forcing frequency is too large (are you sure you are in rad s-1?).");

    if (use_prop_matrix)
    {
        if (num_layers > 1)
            throw std::logic_error("Currently, TidalPy's propagation matrix technique only works for 1-layer worlds.");
        if (tidalpy::c_to_lower(layer_types[0]) != "solid")
            throw std::invalid_argument("The Propagation matrix technique only works for solid layers.");
        if (not is_static_bylayer[0])
            throw std::invalid_argument("The Propagation matrix technique does not allow for dynamic layers.");
        if (not is_incompressible_bylayer[0])
            throw std::invalid_argument("The Propagation matrix technique does not allow for compressible layers.");
    }

    // The same cap the solver's automatic choice uses, so a caller is never refused a starting radius
    // the solver would have picked itself.
    const double max_start_radius_frac = tidalpy_config_ptr->d_MAX_START_RADIUS_FRAC;
    if ((starting_radius != 0.0) &&
        (starting_radius > max_start_radius_frac * radius_array[total_slices - 1]))
    {
        // Default stream precision prints 0.90 as "90", not "90.000000".
        std::ostringstream message;
        message << "Starting radius is above " << (100.0 * max_start_radius_frac)
                << "% of the planet radius (config_x [numerical].max_start_radius_fraction)."
                << " Try a lower radius.";
        throw std::invalid_argument(message.str());
    }

    if (radius_array[0] != 0.0)
    {
        throw std::invalid_argument("Radius array must start at zero.");
    }

    double last_layer_radius = 0.0;
    for (size_t layer_i = 0; layer_i < num_layers; ++layer_i)
    {
        bool top_layer      = (layer_i == num_layers - 1);
        double layer_radius = upper_radius_bylayer_array[layer_i];

        size_t slice_check       = 0;
        size_t layer_check       = 0;
        double last_radius_check = 0.0;

        for (size_t slice_i = 0; slice_i < total_slices; ++slice_i)
        {
            double radius_check = radius_array[slice_i];
            if (radius_check < 0.0)
                throw std::invalid_argument("A negative radius value was found in `radius_array`.");
            if (radius_check < last_radius_check)
                throw std::invalid_argument("Radius array must be in ascending order.");
            if (c_isclose(radius_check, layer_radius, 1.0e-9, 0.0))
                layer_check += 1;
            if (last_layer_radius <= radius_check && radius_check <= layer_radius)
                slice_check += 1;
            
            last_radius_check = radius_check;
        }

        last_layer_radius = layer_radius;

        if (slice_check < C_RS_MIN_SLICES_PER_LAYER)
            throw std::invalid_argument(
                "A minimum of " + std::to_string(C_RS_MIN_SLICES_PER_LAYER) +
                " sub-slices (including top and bottom) are required for each layer.");

        if (top_layer)
        {
            if (layer_check != 1)
                throw std::invalid_argument("Radius of layer " + std::to_string(layer_i) + " found " + std::to_string(layer_check) + " times. Expected 1 time (non-interface layer).");
        }
        else
        {
            if (layer_check != 2)
                throw std::invalid_argument("Radius of layer " + std::to_string(layer_i) + " found " + std::to_string(layer_check) + " times. Expected 2 times (interface layer).");
        }
    }

    bool dynamic_liquid = false;
    for (size_t layer_i = 0; layer_i < num_layers; ++layer_i)
    {
        const std::string l_type = tidalpy::c_to_lower(layer_types[layer_i]);

        if (l_type == "solid")
        {
            layer_types_out_ptr[layer_i] = 0;
        }
        else if (l_type == "liquid")
        {
            layer_types_out_ptr[layer_i] = 1;
            if (not is_static_bylayer[layer_i])
            {
                dynamic_liquid = true;
            }
        }
        else
        {
            layer_types_out_ptr[layer_i] = -1;
            throw std::invalid_argument("Layer type " + layer_types[layer_i] + " is not supported.");
        }
    }

    if (dynamic_liquid && std::abs(frequency) < 2.5e-5 && warnings)
    {
        printf("WARNING: Dynamic liquid layer detected in RadialSolver for a small frequency. Results may be unstable. Extra care is advised!\n");
    }

    // Every layer interpolates its supplied profile, the only EOS method this API offers.
    for (const std::string& eos_method : eos_method_bylayer)
    {
        if (tidalpy::c_to_lower(eos_method) != "interpolate")
        {
            throw std::logic_error("Unknown EOS method provided: " + eos_method);
        }
    }

    num_bc_models_out = solve_for.size();
    if (num_bc_models_out == 0)
    {
        num_bc_models_out = 1;
        bc_models_out_ptr[0] = 1;
    }
    else
    {
        for (size_t i = 0; i < num_bc_models_out; ++i)
        {
            std::string solve_tmp = tidalpy::c_to_lower(solve_for[i]);
            if (solve_tmp == "free")         bc_models_out_ptr[i] = 0;
            else if (solve_tmp == "tidal")   bc_models_out_ptr[i] = 1;
            else if (solve_tmp == "loading") bc_models_out_ptr[i] = 2;
            else throw std::invalid_argument("Unsupported value provided for `solve_for`: " + solve_tmp);
        }
    }
}