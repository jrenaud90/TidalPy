"""A live world's ``get_config_dict`` is exactly what the builder reads, and its ``save_to_toml`` fallback builds."""
import math

import pytest

from TidalPy.Structures import build_world
from TidalPy.Structures.configs import build_system
from TidalPy.Structures.configs.world_builder import construct_world
from TidalPy.Structures.configs.toml_loader import SCHEMA_VERSION, WORLD_TYPES, validate_world_config
from TidalPy.Structures.worlds.base import BUILDER_WORLD_TYPES, BaseWorld
from TidalPy.Structures.worlds.terrestrial import TerrestrialWorld
from TidalPy.Structures.layers import Layer
from TidalPy.Material import Material, Phase
from TidalPy.Rheology.rheology import Andrade, Elastic
from TidalPy.Viscosity import make_viscosity
from TidalPy.Cooling import make_cooling
from TidalPy.Radiogenics.radiogenics import IsotopeRadiogenics
from TidalPy.Tides.classes import make_tide


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
    if len(world) > 0:
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
    # Forces the get_config_dict fallback of save_to_toml.
    world.source_config = None
    path = tmp_path / f"{world_name}_live.toml"
    world.save_to_toml(str(path))
    rebuilt = build_world(str(path))
    assert rebuilt.name == world.name
    assert _nan_equal(rebuilt.get_config_dict(), world.get_config_dict())


def test_hand_built_world_rebuilds_from_config_dict():
    radius = 6.0e6
    mass = (4.0 / 3.0) * math.pi * radius ** 3 * 4000.0
    world = TerrestrialWorld("handmade", radius, mass)
    core_material = Material(solid=Phase(
        eos={"model": "constant", "reference_density_kg_m3": 8000.0, "bulk_modulus_pa": 3.0e11},
        shear_modulus={"model": "constant", "shear_modulus_pa": 1.0e11},
        shear_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e22},
        bulk_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e30},
    ))
    core = Layer(
        "core",
        0,
        0.0,
        0.5 * radius,
        0.4 * mass,
        material=core_material,
        shear_rheology=Elastic(),
    )
    # A melting mantle: a liquid phase, melting curves, and Henning weakening.
    mantle_material = Material(
        solid=Phase(
            eos={"model": "birch_murnaghan", "reference_density_kg_m3": 3300.0, "reference_bulk_modulus_pa": 1.3e11,
                 "bulk_modulus_derivative": 4.0},
            shear_modulus={"model": "constant", "shear_modulus_pa": 6.0e10},
            shear_viscosity=make_viscosity("reference", {"reference_viscosity_pas": 1.0e21}),
            bulk_viscosity=make_viscosity("constant", {"reference_viscosity_pas": 1.0e30}),
        ),
        liquid=Phase(
            eos={"model": "constant", "reference_density_kg_m3": 2750.0, "bulk_modulus_pa": 2.0e10},
            shear_viscosity={"model": "constant", "reference_viscosity_pas": 0.2},
        ),
        solidus={"model": "constant", "temperature_k": 1600.0},
        liquidus={"model": "constant", "temperature_k": 2000.0},
        weakening="henning",
    )
    mantle = Layer(
        "mantle",
        1,
        0.5 * radius,
        radius,
        0.6 * mass,
        material=mantle_material,
        use_melting=True,
        shear_rheology=Andrade(0.3, 1.0),
        bulk_rheology=Elastic(),
    )
    mantle.cooling = make_cooling("convective")
    mantle.radiogenics = IsotopeRadiogenics(isotopes="modern_day_chondritic")
    world.add_layer(core)
    world.add_layer(mantle)
    tide = make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [50.0]})
    tide_name = tide.model_name
    world.set_tide_model(tide)
    world.set_tide_config(min_degree_l=2, max_degree_l=3, eccentricity_truncation=4, obliquity_truncation=2)

    cfg = world.get_config_dict()
    assert cfg["type"] == world.world_type
    assert cfg["type"] in WORLD_TYPES
    assert list(cfg["layers"]) == ["core", "mantle"]
    # One layer class: no class or material type is written, and each layer carries its material table.
    assert all("class" not in layer and "type" not in layer for layer in cfg["layers"].values())
    assert cfg["layers"]["core"]["material"]["solid"]["eos"]["reference_density_kg_m3"] == 8000.0
    assert cfg["layers"]["mantle"]["use_melting"] is True
    assert cfg["layers"]["mantle"]["material"]["melting"]["weakening"]["model"] == "henning"
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


# An iron core and a silicate mantle, each with the values the earlier schema's iron and mantle_rock defaults gave.
_IRON = {"solid": {
    "thermal_conductivity_w_mk": 7.95,
    "heat_capacity_j_kgk": 840.0,
    "eos": {"model": "constant", "reference_density_kg_m3": 8000.0, "bulk_modulus_pa": 1.6e11,
            "thermal_expansion_1_k": 1.2e-5},
    "shear_modulus": {"model": "constant", "shear_modulus_pa": 5.25e10},
    "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e20},
    "bulk_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e22}}}
_MANTLE_ROCK = {"solid": {
    "thermal_conductivity_w_mk": 3.75,
    "heat_capacity_j_kgk": 1200.0,
    "eos": {"model": "constant", "reference_density_kg_m3": 3500.0, "bulk_modulus_pa": 2.0e11,
            "thermal_expansion_1_k": 5.2e-5},
    "shear_modulus": {"model": "constant", "shear_modulus_pa": 6.0e10},
    "shear_viscosity": {"model": "reference", "reference_viscosity_pas": 1.0e22, "reference_temperature_k": 1000.0,
                        "molar_activation_energy_j_mol": 3.0e5, "molar_activation_volume_m3_mol": 0.0},
    "bulk_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e22}}}


def _liquid_core_config(**core_flags):
    """A two-layer world whose core carries the given radial-solver flags."""
    core = {
        "layer_index": 0,
        "radius_outer_m": 3.0e6,
        "use_tides": False,
        "material": _IRON,
        "shear_rheology": {"model": "maxwell"},
        "cooling": {"model": "off"},
        "radiogenics": {"model": "off"},
    }
    core.update(core_flags)
    return {
        "schema_version": SCHEMA_VERSION,
        "name": "liquid_core",
        "type": "terrestrial",
        "radius_m": 6.0e6,
        "mass_kg": 5.0e24,
        "layers": {
            "core": core,
            "mantle": {
                "layer_index": 1,
                "radius_fraction": 1.0,
                "material": _MANTLE_ROCK,
                "shear_rheology": {"model": "andrade", "alpha": 0.3, "zeta": 1.0},
                "cooling": {"model": "convection", "convection_alpha": 1.0, "convection_beta": 1.0 / 3.0,
                            "critical_rayleigh": 1100.0},
                "radiogenics": {"model": "isotope", "isotopes": "modern_day_chondritic"},
            },
        },
    }


def test_layer_assumption_flags_round_trip_through_the_builder(tmp_path):
    world = construct_world(_liquid_core_config(state="liquid", is_incompressible=True))
    assert (world.core.is_liquid, world.core.is_static, world.core.is_incompressible) == (True, True, True)
    assert (world.mantle.is_liquid, world.mantle.is_static, world.mantle.is_incompressible) == (False, True, False)

    cfg = world.get_config_dict()
    assert (cfg["layers"]["core"]["state"], cfg["layers"]["core"]["is_incompressible"]) == ("liquid", True)
    validate_world_config(cfg)
    assert _nan_equal(construct_world(cfg).get_config_dict(), cfg)

    # Forces the get_config_dict fallback of save_to_toml.
    world.source_config = None
    path = tmp_path / "liquid_core.toml"
    world.save_to_toml(str(path))
    reloaded = build_world(str(path))
    assert (reloaded.core.is_liquid, reloaded.core.is_incompressible) == (True, True)


def test_liquid_layer_key_matches_setting_the_flag_on_the_built_layer():
    """A liquid layer declared in the config gives the same Love number as flagging the built layer."""
    frequency = 2.0 * math.pi / 86400.0
    declared = construct_world(_liquid_core_config(state="liquid"))
    flagged = construct_world(_liquid_core_config())
    flagged.core.state = "liquid"
    solid = construct_world(_liquid_core_config())
    for world in (declared, flagged, solid):
        world.solve_eos(verbose=False)
        world.solve_love_numbers(frequency)

    assert declared.love_success and flagged.love_success and solid.love_success
    assert declared.love_number_k == flagged.love_number_k
    assert declared.love_number_k.real > 1.1 * solid.love_number_k.real


def test_bare_base_world_fallback_save_is_rejected(tmp_path):
    """A world with no layers cannot be rebuilt, so the fallback refuses it; the verbatim dump still writes."""
    world = BaseWorld("bare", 1.0e6, 1.0e20)
    with pytest.raises(ValueError):
        world.save_to_toml(str(tmp_path / "bare.toml"))
    world.save_config(str(tmp_path / "bare_verbatim.toml"))


def test_duplicate_layer_names_are_rejected():
    world = BaseWorld("dup", 2.0e6, 1.0e22)
    world.add_layer(Layer(
        "shell",
        0,
        0.0,
        1.0e6,
        5.0e21,
    ))
    with pytest.raises(ValueError):
        world.add_layer(Layer(
            "shell",
            1,
            1.0e6,
            2.0e6,
            5.0e21,
        ))
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
