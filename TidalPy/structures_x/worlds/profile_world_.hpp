#pragma once
/*
 * profile_world_.hpp: build a layered world from a radial profile, entirely in C++.
 *
 * The standalone radial solver is handed a planet as contiguous arrays (radius, density, and the two static
 * moduli) plus its layer boundaries and per-layer solver assumptions. This turns those into a c_LayeredWorld
 * whose layers each carry an interpolated material EOS over their own slice of the profile, which is what the
 * world-attached EOS and Love solves already consume.
 *
 * It is the C++ twin of `build_world_from_layered_profile` in structures_x/configs/world_builder.py, which
 * reaches the same place through a configuration dict. That route turns every slice into a Python float,
 * validates a schema, and merges per-material defaults on every call; this one copies the slices straight
 * into the EOS tables. The two must agree layer for layer, so the rules they share are written the same way
 * here: the slice partition comes from c_partition_radius_by_layer (the one place that rule lives), a layer
 * takes the last radius of its own slice as its outer radius rather than the declared boundary, is_tidal
 * follows is_solid, and the material name is left empty so no per-material defaults are pulled in.
 *
 * The world built here anchors the lifetime of a solve and is never handed back to Python, so it takes none
 * of the world-level extras construct_world attaches (tide model, albedo and emissivity defaults, the
 * retained source configuration). A caller that wants those builds a world the ordinary way.
 */

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "layered_.hpp"                                   // c_LayeredWorld, c_WorldConfig
#include "constants_.hpp"                                 // TidalPyConstants::d_PI
#include "../layers/solidliquid_.hpp"                     // c_SolidLiquidLayer, c_SolidLiquidConfig
#include "../../Material_x/eos/material_eos_.hpp"         // c_InterpolatedEOS, c_MaterialEOSConfig
#include "../../Utilities_x/arrays/layer_partition_.hpp"  // c_partition_radius_by_layer

namespace tidalpy {

/// Smallest number of profile slices an interpolated layer can be built from.
inline constexpr std::size_t d_PROFILE_MIN_SLICES_PER_LAYER = 2;

/// Build a layered world whose layers interpolate their own slice of a radial profile.
///
/// Parameters
/// ----------
/// radius_ptr, density_ptr, shear_modulus_ptr, bulk_modulus_ptr : const double*
///     The profile [m, kg m-3, Pa, Pa], num_slices long, ascending in radius with each interface radius
///     repeated. The moduli are the static (unrelaxed) ones; a viscoelastic response is supplied separately
///     to the Love solve.
/// num_slices : size_t
///     Length of the profile arrays.
/// upper_radius_bylayer_ptr : const double*
///     Upper radius of each layer [m], inner to outer, num_layers long.
/// layer_type_ptr : const int*
///     Per-layer type code, 0 for solid (the encoding the radial-solver input check produces).
/// is_static_ptr, is_incompressible_ptr : const bool*
///     Per-layer radial-solver assumptions.
/// num_layers : size_t
///     Number of layers.
/// planet_bulk_density : double
///     Bulk density [kg m-3], which fixes the world's mass and so its non-dimensional scales.
/// name : const std::string&
///     Name for the constructed world.
///
/// Returns
/// -------
/// shared_ptr<c_LayeredWorld>
///     The world, with every layer added inner to outer and its EOS attached.
///
/// Throws
/// ------
/// std::invalid_argument
///     If a pointer is null, there are no layers or slices, or a layer holds fewer than two slices.
inline std::shared_ptr<c_LayeredWorld> c_build_world_from_layered_profile(
        const double* radius_ptr,
        const double* density_ptr,
        const double* shear_modulus_ptr,
        const double* bulk_modulus_ptr,
        const std::size_t num_slices,
        const double* upper_radius_bylayer_ptr,
        const int* layer_type_ptr,
        const bool* is_static_ptr,
        const bool* is_incompressible_ptr,
        const std::size_t num_layers,
        const double planet_bulk_density,
        const std::string& name)
{
    if ((radius_ptr == nullptr) || (density_ptr == nullptr) || (shear_modulus_ptr == nullptr)
        || (bulk_modulus_ptr == nullptr) || (upper_radius_bylayer_ptr == nullptr)
        || (layer_type_ptr == nullptr) || (is_static_ptr == nullptr) || (is_incompressible_ptr == nullptr))
    {
        throw std::invalid_argument("TidalPy: a null array was passed to the layered-profile world builder.");
    }
    if ((num_slices == 0) || (num_layers == 0))
    {
        throw std::invalid_argument(
            "TidalPy: the layered-profile world builder needs at least one layer and one slice.");
    }

    // The one place the interface-duplication rule lives, so this and the solves cannot disagree about which
    // copy of a boundary radius belongs to which layer.
    std::vector<std::size_t> first_slice_bylayer;
    std::vector<std::size_t> num_slices_bylayer;
    c_partition_radius_by_layer(
        radius_ptr, num_slices, upper_radius_bylayer_ptr, num_layers,
        first_slice_bylayer, num_slices_bylayer);

    const double planet_radius = radius_ptr[num_slices - 1];
    const double planet_mass   = planet_bulk_density * (4.0 / 3.0) * TidalPyConstants::d_PI
                               * planet_radius * planet_radius * planet_radius;

    c_WorldConfig world_cfg;
    world_cfg.name           = name;
    world_cfg.world_type_str = "layered";
    world_cfg.radius         = planet_radius;
    world_cfg.mass           = planet_mass;

    std::shared_ptr<c_LayeredWorld> world = std::make_shared<c_LayeredWorld>(world_cfg);

    double radius_inner = 0.0;
    for (std::size_t layer_i = 0; layer_i < num_layers; ++layer_i)
    {
        const std::size_t first = first_slice_bylayer[layer_i];
        const std::size_t count = num_slices_bylayer[layer_i];
        if (count < d_PROFILE_MIN_SLICES_PER_LAYER)
        {
            throw std::invalid_argument(
                std::string("TidalPy: layer ") + std::to_string(layer_i)
                + " of the supplied profile holds fewer than two points; a layer needs at least two to "
                  "interpolate across.");
        }
        const std::size_t stop = first + count;

        // The material is the layer's slice of the profile. The viscosity tables stay empty: the profile
        // carries none, and the Love solve is handed complex moduli instead.
        c_MaterialEOSConfig eos_cfg;
        eos_cfg.radius.assign(radius_ptr + first, radius_ptr + stop);
        eos_cfg.density.assign(density_ptr + first, density_ptr + stop);
        eos_cfg.shear_modulus.assign(shear_modulus_ptr + first, shear_modulus_ptr + stop);
        eos_cfg.bulk_modulus.assign(bulk_modulus_ptr + first, bulk_modulus_ptr + stop);

        const bool is_solid = (layer_type_ptr[layer_i] == 0);

        c_SolidLiquidConfig layer_cfg;
        layer_cfg.name         = std::string("layer_") + std::to_string(layer_i);
        layer_cfg.layer_index  = static_cast<int>(layer_i);
        layer_cfg.radius_inner = radius_inner;
        // The slice's own top, not the declared boundary: the two agree when the boundary falls on the grid,
        // and where they do not the material and the layer must still end at the same radius.
        layer_cfg.radius_outer = radius_ptr[stop - 1];
        // The mass has no meaningful value yet; every successful EOS solve overwrites it with the solved one.
        layer_cfg.mass = 0.0;
        // An empty material name pulls in no per-material defaults: the profile is the material, so the layer
        // gets no viscosity law, no partial-melt model, and no rheology it did not ask for.
        layer_cfg.material_name     = std::string();
        layer_cfg.is_tidal          = is_solid;
        layer_cfg.is_solid          = is_solid;
        layer_cfg.is_static         = is_static_ptr[layer_i];
        layer_cfg.is_incompressible = is_incompressible_ptr[layer_i];

        std::unique_ptr<c_SolidLiquidLayer> layer = std::make_unique<c_SolidLiquidLayer>(layer_cfg);
        layer->set_eos(std::make_unique<c_InterpolatedEOS>(eos_cfg));
        // Checks continuity against the layer below before taking ownership.
        world->add_layer(std::move(layer));

        radius_inner = layer_cfg.radius_outer;
    }

    return world;
}

}  // namespace tidalpy
