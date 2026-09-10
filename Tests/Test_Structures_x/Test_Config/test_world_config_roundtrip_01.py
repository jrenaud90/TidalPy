"""Round trip of a live world's configuration dict through the TOML world builder.

``get_config_dict()`` must emit exactly what the builder reads (the schema 0.2.0 world table), so a world can be
rebuilt from its own dict, and a directly constructed world's ``save_to_toml`` fallback must write a buildable file.
"""
import math

import pytest

from TidalPy.structures_x import build_world
from TidalPy.structures_x.configs import build_system
from TidalPy.structures_x.configs.world_builder import construct_world
from TidalPy.structures_x.configs.toml_loader import SCHEMA_VERSION, WORLD_TYPES, validate_world_config
from TidalPy.structures_x.worlds.base import BUILDER_WORLD_TYPES, BaseWorld
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.structures_x.layers.solidliquid import SolidLiquidLayer
from TidalPy.Material_x.eos.material_eos import BirchMurnaghanEOS, ConstantDensityEOS
from TidalPy.rheology_x.rheology import Andrade, Elastic
from TidalPy.viscosity_x import make_viscosity
from TidalPy.partial_melt_x import make_partial_melt
from TidalPy.cooling_x import make_cooling
from TidalPy.radiogenics_x.radiogenics import IsotopeRadiogenics
from TidalPy.Tides_x.classes import make_tide


def _nan_equal(left, right):
    """Dict equality that treats NaN as equal to NaN (unset static viscosities are NaN)."""
    if isinstance(left, dict) and isinstance(right, dict):
        return left.keys() == right.keys() and all(_nan_equal(left[key], right[key]) for key in left)
    if isinstance(left, (list, tuple)) and isinstance(right, (list, tuple)):
        return len(left) == len(right) and all(_nan_equal(a, b) for a, b in zip(left, right))
    if isinstance(left, float) and isinstance(right, float):
        if math.isnan(left) and math.isnan(right):
            return True
        return math.isclose(left, right, rel_tol=1e-12, abs_tol=0.0)
    return left == right


def test_builder_world_types_match_loader():
    assert BUILDER_WORLD_TYPES == WORLD_TYPES


@pytest.mark.parametrize("world_name", ["earth_simple", "jupiter_simple", "sol"])
def test_bundled_world_rebuilds_from_config_dict(world_name):
    world = build_world(world_name)
    cfg = world.get_config_dict()
    assert cfg["schema_version"] == SCHEMA_VERSION
    assert "world_type" not in cfg
    assert "num_layers" not in cfg
    validate_world_config(cfg)

    rebuilt = construct_world(cfg)
    assert type(rebuilt) is type(world)
    assert rebuilt.name == world.name
    assert rebuilt.radius == pytest.approx(world.radius)
    assert rebuilt.mass == pytest.approx(world.mass)
    assert rebuilt.get_tide_config() == world.get_tide_config()
    assert _nan_equal(rebuilt.get_config_dict(), cfg)
    if isinstance(world, LayeredWorld):
        assert list(cfg["layers"]) == [layer.name for layer in world]
        world.solve_eos(verbose=False)
        rebuilt.solve_eos(verbose=False)
        radius = 0.5 * world.radius
        assert rebuilt.get_density(radius) == pytest.approx(world.get_density(radius), rel=1e-9)
    else:
        assert rebuilt.effective_temperature == pytest.approx(world.effective_temperature)
        assert rebuilt.luminosity == pytest.approx(world.luminosity)


@pytest.mark.parametrize("world_name", ["earth_simple", "sol"])
def test_directly_built_world_writes_a_buildable_file(world_name, tmp_path):
    world = build_world(world_name)
    world.source_config = None   # force the get_config_dict fallback of save_to_toml
    path = tmp_path / f"{world_name}_live.toml"
    world.save_to_toml(str(path))
    rebuilt = build_world(str(path))
    assert rebuilt.name == world.name
    assert _nan_equal(rebuilt.get_config_dict(), world.get_config_dict())


def test_hand_built_world_rebuilds_from_config_dict():
    radius = 6.0e6
    mass = (4.0 / 3.0) * math.pi * radius ** 3 * 4000.0
    world = LayeredWorld("handmade", radius, mass)
    core = PhysicsLayer(
        "core", 0, 0.0, 0.5 * radius, 0.4 * mass,
        shear_modulus_static_pa=1.0e11, bulk_modulus_static_pa=3.0e11,
        shear_viscosity_static_pas=1.0e22, bulk_viscosity_static_pas=1.0e30)
    core.set_eos(ConstantDensityEOS(reference_density_kg_m3=8000.0))
    core.set_shear_rheology(Elastic())
    mantle = SolidLiquidLayer(
        "mantle", 1, 0.5 * radius, radius, 0.6 * mass,
        shear_modulus_static_pa=6.0e10, bulk_modulus_static_pa=1.3e11,
        shear_viscosity_static_pas=1.0e21, bulk_viscosity_static_pas=1.0e30)
    mantle.set_eos(BirchMurnaghanEOS(
        reference_density_kg_m3=3300.0, reference_bulk_modulus_pa=1.3e11, bulk_modulus_derivative=4.0))
    mantle.set_shear_rheology(Andrade(0.3, 1.0))
    mantle.set_bulk_rheology(Elastic())
    mantle.set_shear_viscosity(make_viscosity("reference", {"reference_viscosity": 1.0e21}))
    mantle.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity": 1.0e30}))
    mantle.set_partial_melt(make_partial_melt("henning"))
    mantle.set_cooling(make_cooling("convective"))
    mantle.set_radiogenics(IsotopeRadiogenics.from_dataset("modern_day_chondritic"))
    world.add_layer(core)
    world.add_layer(mantle)
    tide = make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [50.0]})
    tide_name = tide.model_name
    world.set_tide_model(tide)
    world.set_tide_config(min_degree_l=2, max_degree_l=3, eccentricity_truncation=4, obliquity_truncation=2)

    cfg = world.get_config_dict()
    # The stored world_type label ("terrestrial" by default) is a builder type, so it is kept as is.
    assert cfg["type"] == world.world_type
    assert cfg["type"] in WORLD_TYPES
    assert list(cfg["layers"]) == ["core", "mantle"]
    assert cfg["layers"]["core"]["class"] == "physics"
    assert cfg["layers"]["mantle"]["class"] == "solidliquid"
    assert cfg["layers"]["mantle"]["radiogenics"]["isotope_names"][:2] == ["U238", "U235"]
    assert cfg["tides"]["global_tidal_model"] == tide_name
    assert cfg["tides"]["eccentricity_trunc_lvl"] == 4
    assert cfg["tides"]["obliquity_trunc_lvl"] == 2
    validate_world_config(cfg)

    rebuilt = construct_world(cfg)
    assert _nan_equal(rebuilt.get_config_dict(), cfg)
    assert rebuilt.calc_internal_heating(0.0) == pytest.approx(world.calc_internal_heating(0.0), rel=1e-12)
    world.solve_eos(verbose=False)
    rebuilt.solve_eos(verbose=False)
    query_radius = 0.75 * radius
    assert rebuilt.get_density(query_radius) == pytest.approx(world.get_density(query_radius), rel=1e-9)


def test_bare_base_world_fallback_save_is_rejected(tmp_path):
    """A world with no layers cannot be rebuilt, so the validated fallback refuses to write it."""
    world = BaseWorld("bare", 1.0e6, 1.0e20)
    with pytest.raises(ValueError):
        world.save_to_toml(str(tmp_path / "bare.toml"))
    # The verbatim base-class dump is still available for inspection.
    world.save_config(str(tmp_path / "bare_verbatim.toml"))


def test_duplicate_layer_names_are_rejected():
    world = LayeredWorld("dup", 2.0e6, 1.0e22)
    world.add_layer(PhysicsLayer("shell", 0, 0.0, 1.0e6, 5.0e21))
    with pytest.raises(ValueError):
        world.add_layer(PhysicsLayer("shell", 1, 1.0e6, 2.0e6, 5.0e21))
        world.get_config_dict()


def test_system_with_directly_built_worlds_round_trips(tmp_path):
    """System save falls back to each world's get_config_dict, which must rebuild through build_system."""
    system = build_system("sol_system")
    system.source_config = None
    for world in system:
        world.source_config = None
    path = tmp_path / "live_system.toml"
    system.save_to_toml(str(path))
    rebuilt = build_system(str(path))
    assert [world.name for world in rebuilt] == [world.name for world in system]
