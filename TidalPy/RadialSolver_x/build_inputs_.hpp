#pragma once
/*
 * build_inputs_.hpp: Native input builders for the `_x` radial solver.
 *
 * Both builders assemble the array-based inputs that `radial_solver` expects (radius, density,
 * complex bulk/shear moduli, per-layer upper radii, planet bulk density) from a layer description:
 *
 *   c_build_rs_input_homogeneous_layers : every layer has constant properties; the radial grid is
 *                                         generated from the layer thickness fractions and slice counts.
 *   c_build_rs_input_from_data          : properties are supplied on a user radius grid; the grid is
 *                                         repaired so that r = 0 is present, every interface radius
 *                                         appears twice (top of the lower layer, base of the upper
 *                                         layer), and every layer's upper radius is present.
 *
 * The complex moduli come from the `rheology_x` models: the caller passes one `c_RheologyBase*` per
 * layer for shear and for bulk (a single model applied to every layer is resolved before reaching
 * this header, so the C++ side always sees per-layer arrays).
 *
 * Invalid input throws std::invalid_argument (Cython maps it to ValueError). Grid repairs are logged
 * as warnings through the TidalPy spdlog logger when `warnings` is set.
 */

#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "constants_.hpp"
#include "logger_.hpp"
#include "numerics_.hpp"
#include "rheology_base_.hpp"

namespace tidalpy {

// Minimum number of radial slices a layer must have for the shooting method's starting conditions.
inline constexpr std::size_t C_RS_MIN_SLICES_PER_LAYER = 5;

// ---------------------------------------------------------------------------------------------------------------------
// c_RadialSolverInputs
//
// Output container for both builders. All vectors are sized to the final number of slices (or
// layers); the planet bulk density is the mass-weighted mean density of the assembled structure.
// ---------------------------------------------------------------------------------------------------------------------
struct c_RadialSolverInputs {
    std::vector<double> radius;                                   // [m]      one entry per slice
    std::vector<double> density;                              // [kg/m3]  one entry per slice
    std::vector<std::complex<double>> complex_bulk_modulus;      // [Pa]     one entry per slice
    std::vector<std::complex<double>> complex_shear_modulus;     // [Pa]     one entry per slice
    std::vector<double> upper_radius_bylayer;                     // [m]      one entry per layer
    std::vector<std::size_t> slices_bylayer;                        //          one entry per layer
    double forcing_frequency = TidalPyConstants::d_NAN;       // [rad/s]
    double planet_bulk_density = TidalPyConstants::d_NAN;     // [kg/m3]
};

namespace detail {

inline void c_check_layer_vector_size(
        const std::vector<double>& values,
        std::size_t num_layers,
        const char* name)
{
    if (values.size() != num_layers) {
        throw std::invalid_argument(
            std::string("TidalPy: build_rs_input — `") + name +
            "` must have one entry per layer (" + std::to_string(num_layers) + "), found " +
            std::to_string(values.size()) + ".");
    }
}

inline void c_check_rheology_vector(
        const std::vector<const c_RheologyBase*>& models,
        std::size_t num_layers,
        const char* name)
{
    if (models.size() != num_layers) {
        throw std::invalid_argument(
            std::string("TidalPy: build_rs_input — `") + name +
            "` must have one rheology model per layer (" + std::to_string(num_layers) + "), found " +
            std::to_string(models.size()) + ".");
    }
    for (std::size_t layer_i = 0; layer_i < num_layers; ++layer_i) {
        if (models[layer_i] == nullptr) {
            throw std::invalid_argument(
                std::string("TidalPy: build_rs_input — `") + name + "` entry " +
                std::to_string(layer_i) + " is a null rheology model.");
        }
    }
}

// Evaluate one layer's complex moduli in place: slices [first, first + count) of `out` are filled
// from the static modulus and viscosity vectors using the layer's rheology model.
inline void c_fill_layer_complex_modulus(
        const c_RheologyBase& model,
        const std::vector<double>& static_modulus,
        const std::vector<double>& viscosity,
        double frequency,
        std::size_t first,
        std::size_t count,
        std::vector<std::complex<double>>& out)
{
    for (std::size_t slice_i = first; slice_i < first + count; ++slice_i) {
        out[slice_i] = model.calc_complex_modulus(
            static_modulus[slice_i], viscosity[slice_i], frequency);
    }
}

} // namespace detail

// ---------------------------------------------------------------------------------------------------------------------
// Layer-size conversions (homogeneous builder helpers)
// ---------------------------------------------------------------------------------------------------------------------

// Radius fractions are the cumulative layer upper radii divided by the planet radius (strictly
// increasing, last entry 1). Converts them to per-layer thickness fractions.
inline void c_thickness_from_radius_fractions(
        const std::vector<double>& radius_fraction_bylayer,
        std::vector<double>& out_thickness_fraction_bylayer)
{
    const std::size_t num_layers = radius_fraction_bylayer.size();
    out_thickness_fraction_bylayer.resize(num_layers);
    double last_fraction = 0.0;
    for (std::size_t layer_i = 0; layer_i < num_layers; ++layer_i) {
        const double fraction = radius_fraction_bylayer[layer_i];
        if (fraction <= last_fraction) {
            throw std::invalid_argument(
                "TidalPy: build_rs_input — `radius_fraction_tuple` entries must increase from the "
                "planet center to the surface with no repeated values.");
        }
        out_thickness_fraction_bylayer[layer_i] = fraction - last_fraction;
        last_fraction = fraction;
    }
    if (!c_isclose(last_fraction, 1.0, 1.0e-9, 0.0)) {
        throw std::invalid_argument(
            "TidalPy: build_rs_input — the last entry of `radius_fraction_tuple` must equal 1 (the "
            "planet surface); found " + std::to_string(last_fraction) + ".");
    }
}

// Volume fractions are each layer's share of the total planet volume (sum to 1). Converts them to
// per-layer thickness fractions.
inline void c_thickness_from_volume_fractions(
        double planet_radius,
        const std::vector<double>& volume_fraction_bylayer,
        std::vector<double>& out_thickness_fraction_bylayer)
{
    const std::size_t num_layers = volume_fraction_bylayer.size();
    out_thickness_fraction_bylayer.resize(num_layers);
    const double planet_radius3 = planet_radius * planet_radius * planet_radius;
    double radius_below = 0.0;
    for (std::size_t layer_i = 0; layer_i < num_layers; ++layer_i) {
        const double volume_fraction = volume_fraction_bylayer[layer_i];
        if (volume_fraction <= 0.0) {
            throw std::invalid_argument(
                "TidalPy: build_rs_input — `volume_fraction_tuple` entries must be positive.");
        }
        const double layer_radius = std::cbrt(volume_fraction * planet_radius3 +
                                              radius_below * radius_below * radius_below);
        out_thickness_fraction_bylayer[layer_i] = (layer_radius - radius_below) / planet_radius;
        radius_below = layer_radius;
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// c_build_rs_input_homogeneous_layers
//
// Each layer has constant density, static moduli, and viscosities. The radial grid of layer `i`
// is `slices_bylayer[i]` evenly spaced points from the layer base to its top (both inclusive), so
// interface radii appear twice in the assembled radius array as the solver requires.
// ---------------------------------------------------------------------------------------------------------------------
inline void c_build_rs_input_homogeneous_layers(
        double planet_radius,
        double forcing_frequency,
        const std::vector<double>& density_bylayer,
        const std::vector<double>& static_bulk_modulus_bylayer,
        const std::vector<double>& static_shear_modulus_bylayer,
        const std::vector<double>& bulk_viscosity_bylayer,
        const std::vector<double>& shear_viscosity_bylayer,
        const std::vector<double>& thickness_fraction_bylayer,
        const std::vector<std::size_t>& slices_bylayer,
        const std::vector<const c_RheologyBase*>& shear_rheology_bylayer,
        const std::vector<const c_RheologyBase*>& bulk_rheology_bylayer,
        c_RadialSolverInputs& out)
{
    const std::size_t num_layers = density_bylayer.size();
    if (num_layers == 0) {
        throw std::invalid_argument("TidalPy: build_rs_input — at least one layer is required.");
    }
    if (!(planet_radius > 0.0) || !std::isfinite(planet_radius)) {
        throw std::invalid_argument("TidalPy: build_rs_input — `planet_radius` must be positive and finite.");
    }
    detail::c_check_layer_vector_size(static_bulk_modulus_bylayer, num_layers, "static_bulk_modulus_tuple");
    detail::c_check_layer_vector_size(static_shear_modulus_bylayer, num_layers, "static_shear_modulus_tuple");
    detail::c_check_layer_vector_size(bulk_viscosity_bylayer, num_layers, "bulk_viscosity_tuple");
    detail::c_check_layer_vector_size(shear_viscosity_bylayer, num_layers, "shear_viscosity_tuple");
    detail::c_check_layer_vector_size(thickness_fraction_bylayer, num_layers, "thickness_fraction_tuple");
    detail::c_check_rheology_vector(shear_rheology_bylayer, num_layers, "shear_rheology_model_tuple");
    detail::c_check_rheology_vector(bulk_rheology_bylayer, num_layers, "bulk_rheology_model_tuple");
    if (slices_bylayer.size() != num_layers) {
        throw std::invalid_argument(
            "TidalPy: build_rs_input — `slices_tuple` must have one entry per layer (" +
            std::to_string(num_layers) + "), found " + std::to_string(slices_bylayer.size()) + ".");
    }

    // Validate the layer sizes and count the slices.
    double total_thickness_fraction = 0.0;
    std::size_t total_slices = 0;
    for (std::size_t layer_i = 0; layer_i < num_layers; ++layer_i) {
        const double thickness_fraction = thickness_fraction_bylayer[layer_i];
        if (!(thickness_fraction > 0.0)) {
            throw std::invalid_argument(
                "TidalPy: build_rs_input — layer " + std::to_string(layer_i) +
                " has a negative or zero thickness fraction.");
        }
        total_thickness_fraction += thickness_fraction;
        const std::size_t layer_slices = slices_bylayer[layer_i];
        if (layer_slices < C_RS_MIN_SLICES_PER_LAYER) {
            throw std::invalid_argument(
                "TidalPy: build_rs_input — layer " + std::to_string(layer_i) + " has " +
                std::to_string(layer_slices) + " slices when at least " +
                std::to_string(C_RS_MIN_SLICES_PER_LAYER) + " are required.");
        }
        total_slices += layer_slices;
    }
    if (!c_isclose(total_thickness_fraction, 1.0, 1.0e-9, 0.0)) {
        throw std::invalid_argument(
            "TidalPy: build_rs_input — layer thickness fractions must sum to 1 (found " +
            std::to_string(total_thickness_fraction) + ").");
    }

    // Allocate the outputs plus the per-slice static properties the rheologies consume.
    out.radius.resize(total_slices);
    out.density.resize(total_slices);
    out.complex_bulk_modulus.resize(total_slices);
    out.complex_shear_modulus.resize(total_slices);
    out.upper_radius_bylayer.resize(num_layers);
    out.slices_bylayer = slices_bylayer;
    out.forcing_frequency = forcing_frequency;
    std::vector<double> static_shear(total_slices);
    std::vector<double> static_bulk(total_slices);
    std::vector<double> shear_viscosity(total_slices);
    std::vector<double> bulk_viscosity(total_slices);

    // Populate layer by layer.
    const double planet_radius3 = planet_radius * planet_radius * planet_radius;
    double planet_bulk_density = 0.0;
    double last_layer_radius = 0.0;
    double last_layer_radius3 = 0.0;
    std::size_t first_slice_in_layer = 0;
    for (std::size_t layer_i = 0; layer_i < num_layers; ++layer_i) {
        const std::size_t layer_slices = slices_bylayer[layer_i];
        const double layer_density = density_bylayer[layer_i];
        const double layer_thickness = thickness_fraction_bylayer[layer_i] * planet_radius;
        const double layer_radius = last_layer_radius + layer_thickness;
        // layer_slices >= 5 was checked above, so no divide by zero here. The "- 1" makes the grid
        // inclusive of the layer top.
        const double dr = layer_thickness / static_cast<double>(layer_slices - 1);

        out.upper_radius_bylayer[layer_i] = layer_radius;
        const double layer_radius3 = layer_radius * layer_radius * layer_radius;
        planet_bulk_density += layer_density * (layer_radius3 - last_layer_radius3);
        last_layer_radius3 = layer_radius3;

        for (std::size_t slice_i = 0; slice_i < layer_slices; ++slice_i) {
            const std::size_t full_index = first_slice_in_layer + slice_i;
            if (slice_i == 0) {
                out.radius[full_index] = last_layer_radius;
            } else if (slice_i == layer_slices - 1) {
                out.radius[full_index] = last_layer_radius + layer_thickness;
            } else {
                out.radius[full_index] = last_layer_radius + static_cast<double>(slice_i) * dr;
            }
            out.density[full_index] = layer_density;
            static_shear[full_index] = static_shear_modulus_bylayer[layer_i];
            static_bulk[full_index] = static_bulk_modulus_bylayer[layer_i];
            shear_viscosity[full_index] = shear_viscosity_bylayer[layer_i];
            bulk_viscosity[full_index] = bulk_viscosity_bylayer[layer_i];
        }

        detail::c_fill_layer_complex_modulus(
            *shear_rheology_bylayer[layer_i], static_shear, shear_viscosity,
            forcing_frequency, first_slice_in_layer, layer_slices, out.complex_shear_modulus);
        detail::c_fill_layer_complex_modulus(
            *bulk_rheology_bylayer[layer_i], static_bulk, bulk_viscosity,
            forcing_frequency, first_slice_in_layer, layer_slices, out.complex_bulk_modulus);

        first_slice_in_layer += layer_slices;
        last_layer_radius = layer_radius;
    }

    out.planet_bulk_density = planet_bulk_density / planet_radius3;
}

// ---------------------------------------------------------------------------------------------------------------------
// c_build_rs_input_from_data
//
// Properties are given on an ascending user radius grid; `layer_upper_radius_bylayer` marks the
// layer tops (the last entry must be the planet radius, the last grid point). The grid is copied
// and repaired: a slice is inserted at r = 0 if missing, at each layer base when the previous
// layer's top is not repeated, and at each layer top when it is absent. Inserted slices copy the
// neighboring slice's properties: a base copies the layer's first provided slice, a top copies the
// slice below it. An interface radius listed only once is taken as the top of the lower layer (with
// the properties it carries), matching the classic builder. The planet bulk density is the mass of
// the piecewise-constant shells (each shell takes the density of its upper slice) divided by the
// planet volume.
// ---------------------------------------------------------------------------------------------------------------------
inline void c_build_rs_input_from_data(
        double forcing_frequency,
        const std::vector<double>& radius,
        const std::vector<double>& density,
        const std::vector<double>& static_bulk_modulus,
        const std::vector<double>& static_shear_modulus,
        const std::vector<double>& bulk_viscosity,
        const std::vector<double>& shear_viscosity,
        const std::vector<double>& layer_upper_radius_bylayer,
        const std::vector<const c_RheologyBase*>& shear_rheology_bylayer,
        const std::vector<const c_RheologyBase*>& bulk_rheology_bylayer,
        bool warnings,
        c_RadialSolverInputs& out)
{
    const std::size_t num_slices_input = radius.size();
    const std::size_t num_layers = layer_upper_radius_bylayer.size();
    if (num_slices_input == 0) {
        throw std::invalid_argument("TidalPy: build_rs_input — `radius_array` is empty.");
    }
    if (num_layers == 0) {
        throw std::invalid_argument("TidalPy: build_rs_input — at least one layer is required.");
    }
    const double planet_radius = radius[num_slices_input - 1];
    if (!(planet_radius > 0.0) || !std::isfinite(planet_radius)) {
        throw std::invalid_argument(
            "TidalPy: build_rs_input — the last entry of `radius_array` (the planet radius) must be "
            "positive and finite.");
    }
    if (!c_isclose(layer_upper_radius_bylayer[num_layers - 1], planet_radius, 1.0e-9, 0.0)) {
        throw std::invalid_argument(
            "TidalPy: build_rs_input — the upper radius of the last layer must equal the planet radius "
            "(the last entry of `radius_array`). Expected " + std::to_string(planet_radius) +
            ", found " + std::to_string(layer_upper_radius_bylayer[num_layers - 1]) + ".");
    }
    detail::c_check_rheology_vector(shear_rheology_bylayer, num_layers, "shear_rheology_model_tuple");
    detail::c_check_rheology_vector(bulk_rheology_bylayer, num_layers, "bulk_rheology_model_tuple");

    auto check_slice_vector = [num_slices_input](const std::vector<double>& values, const char* name) {
        if (values.size() != num_slices_input) {
            throw std::invalid_argument(
                std::string("TidalPy: build_rs_input — `") + name +
                "` must have the same length as `radius_array` (" + std::to_string(num_slices_input) +
                "), found " + std::to_string(values.size()) + ".");
        }
    };
    check_slice_vector(density, "density_array");
    check_slice_vector(static_bulk_modulus, "static_bulk_modulus_array");
    check_slice_vector(static_shear_modulus, "static_shear_modulus_array");
    check_slice_vector(bulk_viscosity, "bulk_viscosity_array");
    check_slice_vector(shear_viscosity, "shear_viscosity_array");
    for (std::size_t slice_i = 1; slice_i < num_slices_input; ++slice_i) {
        if (radius[slice_i - 1] > radius[slice_i]) {
            throw std::invalid_argument("TidalPy: build_rs_input — `radius_array` must be in ascending order.");
        }
    }
    for (std::size_t layer_i = 1; layer_i < num_layers; ++layer_i) {
        if (layer_upper_radius_bylayer[layer_i] <= layer_upper_radius_bylayer[layer_i - 1]) {
            throw std::invalid_argument(
                "TidalPy: build_rs_input — `layer_upper_radius_tuple` must increase from the planet "
                "center to the surface.");
        }
    }

    // Worst case: one slice inserted at the base and one at the top of every layer.
    const std::size_t max_slices = num_slices_input + 2 * num_layers;
    std::vector<double> radius_use; radius_use.reserve(max_slices);
    std::vector<double> density_use; density_use.reserve(max_slices);
    std::vector<double> shear_use; shear_use.reserve(max_slices);
    std::vector<double> bulk_use; bulk_use.reserve(max_slices);
    std::vector<double> shear_visc_use; shear_visc_use.reserve(max_slices);
    std::vector<double> bulk_visc_use; bulk_visc_use.reserve(max_slices);
    out.slices_bylayer.assign(num_layers, 0);
    out.upper_radius_bylayer = layer_upper_radius_bylayer;
    out.forcing_frequency = forcing_frequency;

    const double pi_4_3 = (4.0 / 3.0) * TidalPyConstants::d_PI;
    double growing_mass = 0.0;
    double last_r3 = 0.0;

    auto push_slice = [&](double radius, std::size_t source_index) {
        radius_use.push_back(radius);
        density_use.push_back(density[source_index]);
        shear_use.push_back(static_shear_modulus[source_index]);
        bulk_use.push_back(static_bulk_modulus[source_index]);
        shear_visc_use.push_back(shear_viscosity[source_index]);
        bulk_visc_use.push_back(bulk_viscosity[source_index]);
        const double r3 = radius * radius * radius;
        growing_mass += pi_4_3 * (r3 - last_r3) * density[source_index];
        last_r3 = r3;
    };

    std::size_t slice_i_input = 0;
    double last_layer_upper_radius = 0.0;
    for (std::size_t layer_i = 0; layer_i < num_layers; ++layer_i) {
        const double layer_upper_radius = layer_upper_radius_bylayer[layer_i];
        const std::size_t first_output_slice = radius_use.size();
        if (slice_i_input >= num_slices_input) {
            throw std::invalid_argument(
                "TidalPy: build_rs_input — `radius_array` ran out of slices before layer " +
                std::to_string(layer_i) + " (check `layer_upper_radius_tuple`).");
        }

        // Base of the layer: the previous layer's top (or r = 0) must appear again.
        if (!c_isclose(radius[slice_i_input], last_layer_upper_radius, 1.0e-9, 0.0)) {
            if (warnings) {
                if (layer_i == 0) {
                    TIDALPY_LOG_WARN(
                        "build_rs_input_from_data: radius array must start at zero; inserting r = 0 at "
                        "slice 0 with the properties of the first provided slice.");
                } else {
                    TIDALPY_LOG_WARN(
                        "build_rs_input_from_data: layer {} does not start at the previous layer's upper "
                        "radius; interface radii must appear twice. Inserting a slice at r = {} m with the "
                        "properties of this layer's first provided slice.",
                        layer_i, last_layer_upper_radius);
                }
            }
            push_slice(last_layer_upper_radius, slice_i_input);
        }

        // Copy the provided slices that belong to this layer.
        bool top_present = false;
        while (slice_i_input < num_slices_input) {
            const double slice_radius = radius[slice_i_input];
            if (c_isclose(slice_radius, layer_upper_radius, 1.0e-9, 0.0)) {
                top_present = true;
            } else if (slice_radius > layer_upper_radius) {
                break;  // First slice of the next layer.
            }
            push_slice(slice_radius, slice_i_input);
            ++slice_i_input;
            if (top_present) {
                break;
            }
        }

        // Top of the layer: insert it when the grid skipped it.
        if (!top_present) {
            if (warnings) {
                TIDALPY_LOG_WARN(
                    "build_rs_input_from_data: layer {} does not have its upper radius ({} m) in the radius "
                    "array; inserting it with the properties of the slice below.",
                    layer_i, layer_upper_radius);
            }
            push_slice(layer_upper_radius, slice_i_input - 1);
        }

        const std::size_t layer_slices = radius_use.size() - first_output_slice;
        if (layer_slices < C_RS_MIN_SLICES_PER_LAYER) {
            throw std::invalid_argument(
                "TidalPy: build_rs_input — layer " + std::to_string(layer_i) + " has " +
                std::to_string(layer_slices) + " slices when at least " +
                std::to_string(C_RS_MIN_SLICES_PER_LAYER) + " are required.");
        }
        out.slices_bylayer[layer_i] = layer_slices;
        last_layer_upper_radius = layer_upper_radius;
    }

    const std::size_t total_slices = radius_use.size();
    out.radius = std::move(radius_use);
    out.density = std::move(density_use);
    out.complex_bulk_modulus.resize(total_slices);
    out.complex_shear_modulus.resize(total_slices);
    out.planet_bulk_density = growing_mass / (pi_4_3 * planet_radius * planet_radius * planet_radius);

    std::size_t first_slice_in_layer = 0;
    for (std::size_t layer_i = 0; layer_i < num_layers; ++layer_i) {
        const std::size_t layer_slices = out.slices_bylayer[layer_i];
        detail::c_fill_layer_complex_modulus(
            *shear_rheology_bylayer[layer_i], shear_use, shear_visc_use, forcing_frequency,
            first_slice_in_layer, layer_slices, out.complex_shear_modulus);
        detail::c_fill_layer_complex_modulus(
            *bulk_rheology_bylayer[layer_i], bulk_use, bulk_visc_use, forcing_frequency,
            first_slice_in_layer, layer_slices, out.complex_bulk_modulus);
        first_slice_in_layer += layer_slices;
    }
}

} // namespace tidalpy
