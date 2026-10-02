#pragma once
/*
 * profile_world_.hpp: build a world from a radial profile, entirely in C++.
 *
 * The standalone radial solver is handed a planet as contiguous arrays (radius, density, and the two static
 * moduli) plus its layer boundaries and per-layer solver assumptions. This turns those into a c_BaseWorld
 * whose layers are each made of a material that interpolates their own slice of the profile (a tabulated equation of
 * state and shear modulus), which is what the world-attached EOS and Love solves already consume.
 *
 * `build_world_from_layered_profile` in Structures/configs/world_builder.py calls this routine, so a profile becomes
 * layers in one place: the slice partition comes from c_partition_radius_by_layer (the one place that rule lives), a
 * layer takes the last radius of its own slice as its outer radius rather than the declared boundary, and use_tides
 * follows whether the layer is solid. A world file's data_file (_interpolated_layer_config) also turns a profile
 * into radius-tabulated layers, but through a configuration dict, and it makes a liquid layer's material
 * liquid-only, where this gives every layer a solid phase holding the profile's shear modulus (zero in a liquid) and
 * sets a liquid layer's state to liquid.
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

#include "base_.hpp"                                    // c_BaseWorld, c_WorldConfig
#include "constants_.hpp"                               // TidalPyConstants::d_PI
#include "../layers/layer_.hpp"                         // c_Layer, c_LayerConfig, c_Material
#include "../../Utilities/arrays/layer_partition_.hpp"  // c_partition_radius_by_layer

namespace tidalpy {

/// Smallest number of profile slices an interpolated layer can be built from.
inline constexpr std::size_t d_PROFILE_MIN_SLICES_PER_LAYER = 2;

/// Build a world whose layers interpolate their own slice of a radial profile.
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
/// shared_ptr<c_BaseWorld>
///     The world, with every layer added inner to outer and its material attached.
///
/// Throws
/// ------
/// std::invalid_argument
///     If a pointer is null, there are no layers or slices, or a layer holds fewer than two slices.
inline std::shared_ptr<c_BaseWorld> c_build_world_from_layered_profile(
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

    std::shared_ptr<c_BaseWorld> world = std::make_shared<c_BaseWorld>(world_cfg);

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

        // The material is the layer's slice of the profile. It has no viscosity law: the profile carries none,
        // and the Love solve is handed complex moduli instead.
        const std::vector<double> layer_radius(radius_ptr + first, radius_ptr + stop);
        c_PhaseComponents phase_components;
        phase_components.eos = std::make_shared<const c_InterpolatedEOS>(c_ParamMap{
            {"radius_m", layer_radius},
            {"density_kg_m3", std::vector<double>(density_ptr + first, density_ptr + stop)},
            {"bulk_modulus_pa", std::vector<double>(bulk_modulus_ptr + first, bulk_modulus_ptr + stop)}});
        phase_components.shear_modulus = std::make_shared<const c_InterpolatedShearModulus>(c_ParamMap{
            {"radius_m", layer_radius},
            {"shear_modulus_pa", std::vector<double>(shear_modulus_ptr + first, shear_modulus_ptr + stop)}});
        c_MaterialComponents material_components;
        material_components.solid = std::make_shared<const c_Phase>(c_ParamMap{}, phase_components);

        const bool is_solid = (layer_type_ptr[layer_i] == 0);

        c_LayerConfig layer_cfg;
        layer_cfg.name         = std::string("layer_") + std::to_string(layer_i);
        layer_cfg.layer_index  = static_cast<int>(layer_i);
        layer_cfg.radius_inner = radius_inner;
        // The slice's own top, not the declared boundary: the two agree when the boundary falls on the grid,
        // and where they do not the material and the layer must still end at the same radius.
        layer_cfg.radius_outer = radius_ptr[stop - 1];
        // The mass has no meaningful value yet; every successful EOS solve overwrites it with the solved one.
        layer_cfg.mass = 0.0;
        layer_cfg.use_tides         = is_solid;
        layer_cfg.state             = is_solid ? c_LayerState::Solid : c_LayerState::Liquid;
        layer_cfg.is_static         = is_static_ptr[layer_i];
        layer_cfg.is_incompressible = is_incompressible_ptr[layer_i];

        std::unique_ptr<c_Layer> layer = std::make_unique<c_Layer>(layer_cfg);
        layer->set_material(std::make_shared<const c_Material>(c_ParamMap{}, material_components));
        // Checks continuity against the layer below before taking ownership.
        world->add_layer(std::move(layer));

        radius_inner = layer_cfg.radius_outer;
    }

    return world;
}

}  // namespace tidalpy
