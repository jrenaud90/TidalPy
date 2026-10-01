"""The world builder: family dispatch, model wiring, the default material, geometry specs, bundled worlds, and
saving."""
import copy
import math
import warnings

import numpy as np
import pytest
import toml

import TidalPy
from TidalPy.constants import G
from TidalPy.Structures import build_world, construct_world, available_worlds, save_world_to_toml
from TidalPy.Structures.configs import world_builder, config_kind, worldpack
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Structures.worlds.terrestrial import TerrestrialWorld
from TidalPy.Structures.worlds.gasgiant import GasGiantWorld
from TidalPy.Structures.worlds.stellar import StarWorld


# =====================================================================================================================
# Helpers
# =====================================================================================================================
def _constant_material(density):
    """A solid of constant density with the moduli the earlier schema's default silicate block gave [kg m-3]."""
    return {"solid": {
        "eos": {"model": "constant", "reference_density_kg_m3": density, "bulk_modulus_pa": 2.0e11},
        "shear_modulus": {"model": "constant", "shear_modulus_pa": 6.0e10}}}


def _terrestrial_dict():
    """A two-layer terrestrial world; the mantle lists every model so nothing depends on the configured default."""
    return {
        "schema_version": "0.2.0",
        "name": "TestEarth",
        "type": "terrestrial",
        "radius_m": 6.0e6,
        "mass_kg": 5.0e24,
        "spin_frequency_rad_s": 7.0e-5,
        "layers": {
            "core": {"layer_index": 0, "radius_outer_m": 3.0e6, "use_tides": False,
                     "material": _constant_material(9000.0)},
            "mantle": {"layer_index": 1, "radius_outer_m": 6.0e6, "mass_kg": 3.0e24, "use_tides": True,
                       "material": {"solid": {
                           "thermal_conductivity_w_mk": 3.75,
                           "heat_capacity_j_kgk": 1200.0,
                           "eos": {"model": "constant", "reference_density_kg_m3": 4000.0, "bulk_modulus_pa": 2.0e11,
                                   "thermal_expansion_1_k": 5.2e-5},
                           "shear_modulus": {"model": "constant", "shear_modulus_pa": 8.0e10},
                           "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e21},
                           "bulk_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e22}}},
                       "shear_rheology": {"model": "maxwell"}, "bulk_rheology": {"model": "elastic"},
                       "cooling": {"model": "convection"},
                       "radiogenics": {"model": "fixed", "fixed_heat_production_w_kg": 5.0e-12}},
        },
    }


def _terrestrial_edited(path, value):
    config = _terrestrial_dict()
    target = config
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    return config


def _single_layer_world(layer_overrides):
    """A one-layer terrestrial world; ``layer_overrides`` is merged into the layer."""
    layer = {
        "layer_index": 0,
        "radius_outer_m": 6.0e6,
    }
    layer.update(layer_overrides)
    return {
        "schema_version": "0.2.0",
        "name": "OneLayer",
        "type": "terrestrial",
        "radius_m": 6.0e6,
        "mass_kg": 5.0e24,
        "layers": {"only": layer},
    }


def _star_dict(**extra):
    return {"schema_version": "0.2.0", "name": "S", "type": "star", "radius_m": 6.957e8, "mass_kg": 1.988e30, **extra}


# =====================================================================================================================
# Dispatch: each world family builds the right underlying class
# =====================================================================================================================
_GASGIANT = {
    "name": "G", "type": "gasgiant", "radius_m": 7.0e7, "mass_kg": 1.9e27,
    "layers": {"env": {"radius_outer_m": 7.0e7,
                       "material": {"liquid": {"eos": {"model": "constant", "reference_density_kg_m3": 1300.0}}}}},
}


@pytest.mark.parametrize("config, expected_class, num_layers", [
    pytest.param(_terrestrial_dict(), TerrestrialWorld, 2, id="terrestrial"),
    pytest.param(_terrestrial_edited(("type",), "layered"), BaseWorld, 2, id="layered-is-the-base-class"),
    pytest.param(_GASGIANT, GasGiantWorld, 1, id="gasgiant"),
])
def test_construct_layered_family_class(config, expected_class, num_layers):
    world = construct_world(config)
    assert type(world) is expected_class
    assert world.num_layers == num_layers


def test_construct_star_class():
    config = {"name": "S", "type": "star", "radius_m": 7.0e8, "mass_kg": 2.0e30, "effective_temperature_k": 5772.0}
    world = construct_world(config)
    assert isinstance(world, StarWorld)
    assert math.isclose(world.effective_temperature, 5772.0)
    assert world.luminosity > 0.0


def test_baseworld_build_dispatches_to_subclass():
    star = BaseWorld.build(
        {"name": "S", "type": "star", "radius_m": 7.0e8, "mass_kg": 2.0e30, "effective_temperature_k": 5772.0})
    assert isinstance(star, StarWorld)


def test_build_world_returns_cython_world():
    """build_world returns the Cython world itself, not a wrapper, and keeps the normalized config."""
    world = build_world(_terrestrial_dict())
    assert isinstance(world, TerrestrialWorld)
    assert world.name == "TestEarth"
    assert world.world_type == "terrestrial"
    assert world.num_layers == 2
    assert world.config["name"] == "TestEarth"
    assert world.solve_eos(G_to_use=G, verbose=False)["success"]


def test_layers_sorted_by_index_regardless_of_declaration_order():
    config = _terrestrial_dict()
    config["layers"] = {
        "mantle": config["layers"]["mantle"],
        "core": config["layers"]["core"],
    }
    world = construct_world(config)
    # A wrong order would leave the geometry discontinuous and the build or solve would fail.
    assert world.num_layers == 2
    assert world.solve_eos(G_to_use=G, verbose=False)["success"]


# =====================================================================================================================
# Model wiring
# =====================================================================================================================
def test_eos_wired_and_solves():
    world = construct_world(_terrestrial_dict())
    assert world.all_materials_set is True
    result = world.solve_eos(G_to_use=G, verbose=False)
    assert result["success"] is True
    assert math.isclose(world.get_density(1.0e6), 9000.0, rel_tol=0.05)
    assert math.isclose(world.get_density(4.5e6), 4000.0, rel_tol=0.05)


def test_radiogenics_wired_produces_heating():
    assert construct_world(_terrestrial_dict()).calc_internal_heating(0.0) > 0.0


def test_rheology_and_viscosity_wired_give_complex_modulus():
    """Maxwell rheology with a finite viscosity gives a dissipative complex shear modulus in the mantle."""
    world = construct_world(_terrestrial_dict())
    world.solve_eos(G_to_use=G, verbose=False)
    mu = world.calc_complex_shear_modulus(4.5e6, 1.0e-5)
    assert np.isfinite(mu.real) and np.isfinite(mu.imag)
    assert mu.imag != 0.0


@pytest.fixture
def private_config(monkeypatch):
    """A private copy of the configuration that a test may edit."""
    private = copy.deepcopy(TidalPy.config)
    monkeypatch.setattr(TidalPy, "config", private)
    return private


def test_a_layer_without_a_material_needs_the_configured_default(private_config):
    """A layer that names no material takes ``[layers] material``; with none configured the build refuses."""
    config = _terrestrial_dict()
    del config["layers"]["core"]["material"]
    del private_config["layers"]["material"]
    with pytest.raises(ValueError, match="Layer 'core' names no material"):
        construct_world(config)


def test_a_layer_without_a_material_takes_the_configured_default():
    """A layer that names no material takes ``[layers] material`` of the configuration, simple_rock."""
    config = _terrestrial_dict()
    del config["layers"]["core"]["material"]
    world = construct_world(config)
    assert world.all_materials_set is True
    world.solve_eos(G_to_use=G, verbose=False)
    assert math.isclose(world.get_density(1.0e6), 3300.0, rel_tol=1e-6)


@pytest.mark.parametrize("config, match", [
    pytest.param(
        _terrestrial_edited(("layers", "mantle", "shear_rheology"), {"model": "not_a_rheology"}),
        None,
        id="unknown-model-name"),
    pytest.param(
        _terrestrial_edited(
            ("layers", "mantle", "material", "solid", "shear_viscosity"),
            {"model": "constant", "reference_viscosty_pas": 1.0e21}),
        r"\[layers\.mantle\.material\].*'reference_viscosty_pas'",
        id="misspelled-material-key-names-its-table"),
    pytest.param(
        _terrestrial_edited(("layers", "mantle", "shear_rheology"), {"model": "andrade", "alpah": 0.3}),
        r"\[layers\.mantle\.shear_rheology\].*'alpah'",
        id="misspelled-model-key-names-its-table"),
    pytest.param(
        _terrestrial_edited(("layers", "mantle", "material"), "no_such_material"),
        r"\[layers\.mantle\.material\].*no MatPack material named 'no_such_material'",
        id="unknown-matpack-name"),
    pytest.param(
        _terrestrial_edited(("layers", "mantle", "material"), {"model": "constant", "reference_density_kg_m3": 3300.0}),
        r"\[layers\.mantle\.material\] .*'solid\.eos'",
        id="flat-material-table-names-the-phase"),
    pytest.param(
        {"name": "X", "type": "terrestrial", "radius_m": 1.0, "mass_kg": 1.0}, None, id="no-layers"),
])
def test_construct_world_rejects_a_bad_config(config, match):
    with pytest.raises(ValueError, match=match):
        construct_world(config)


# =====================================================================================================================
# The layer's material: a MatPack name, a preset with overrides, a full table, or the configured default
# =====================================================================================================================
@pytest.mark.parametrize("layer_overrides, density", [
    pytest.param({}, 3300.0, id="configured-default-material"),
    pytest.param({"material": "simple_ice"}, 920.0, id="matpack-name"),
    pytest.param(
        {"material": {"preset": "simple_rock", "solid": {"eos": {"reference_density_kg_m3": 5200.0}}}},
        5200.0,
        id="user-value-wins"),
    pytest.param({"material": _constant_material(1000.0)}, 1000.0, id="full-table"),
])
def test_single_layer_eos_follows_the_layer_material(layer_overrides, density):
    world = construct_world(_single_layer_world(layer_overrides))
    assert world.all_materials_set is True
    assert world.solve_eos(G_to_use=G, verbose=False)["success"]
    assert math.isclose(world.get_density(3.0e6), density, rel_tol=0.05)


def test_the_configured_default_material_is_editable(private_config):
    private_config["layers"]["material"] = "simple_ice"
    world = construct_world(_single_layer_world({}))
    assert world.only.get_config_dict()["material"]["solid"]["eos"]["reference_density_kg_m3"] == 920.0
    # A table works as the default too.
    private_config["layers"]["material"] = {
        "preset": "simple_ice", "solid": {"eos": {"reference_density_kg_m3": 930.0}}}
    world = construct_world(_single_layer_world({}))
    assert world.only.get_config_dict()["material"]["solid"]["eos"]["reference_density_kg_m3"] == 930.0


def test_an_error_in_the_configured_default_material_names_the_configuration(private_config):
    """A misspelled default is reported against TidalPy_Configs.toml, not against the layer, which names none."""
    private_config["layers"]["material"] = "simple_rok"
    with pytest.raises(ValueError, match=r"TidalPy_Configs\.toml \[layers\] material .*layer 'only'.*simple_rok"):
        construct_world(_single_layer_world({}))


def test_an_isotope_table_naming_no_dataset_takes_the_configured_one():
    """A radiogenics table that names no dataset or arrays takes [radiogenics] isotopes, keeping its other keys."""
    world = construct_world(_single_layer_world({"radiogenics": {"model": "isotope", "ref_time_s": 1.0e17}}))
    held = world.only.get_config_dict()["radiogenics"]
    assert held["isotope_names"] == ["U238", "U235", "Th232", "K40"]
    assert held["ref_time_s"] == 1.0e17


# =====================================================================================================================
# Layer geometry: derived inner radius and the three outer-radius specifiers
# =====================================================================================================================
def _layer_outer_radii(world):
    return [layer["radius_outer_m"] for layer in world.get_config_dict()["layers"].values()]


def _two_layer_world(name, inner_spec):
    """A 6000 km world whose inner layer takes ``inner_spec``; the outer layer fills the rest."""
    inner = {"layer_index": 0}
    inner.update(inner_spec)
    return {
        "name": name, "type": "terrestrial", "radius_m": 6.0e6, "mass_kg": 5.0e24,
        "layers": {
            "inner": inner,
            "outer": {"layer_index": 1, "radius_fraction": 1.0},
        }}


@pytest.mark.parametrize("name, inner_spec, expected, rel_tol", [
    ("O", {"radius_outer_m": 4.0e6}, 4.0e6, 1e-12),
    ("F", {"radius_fraction": 0.5}, 3.0e6, 1e-12),
    # Innermost layer: outer = R f^(1/3), so f = 0.125 gives R / 2.
    ("V", {"volume_fraction": 0.125}, 3.0e6, 1e-9),
])
def test_outer_radius_spec(name, inner_spec, expected, rel_tol):
    outers = _layer_outer_radii(construct_world(_two_layer_world(name, inner_spec)))
    assert math.isclose(outers[0], expected, rel_tol=rel_tol)
    assert math.isclose(outers[1], 6.0e6, rel_tol=1e-12)


def test_inner_radius_derived_from_previous_layer():
    config = {
        "name": "D", "type": "terrestrial", "radius_m": 6.0e6, "mass_kg": 5.0e24,
        "layers": {
            "core": {"layer_index": 0, "radius_outer_m": 2.0e6, "material": _constant_material(8000.0)},
            "mid": {"layer_index": 1, "radius_fraction": 0.75, "material": _constant_material(5000.0)},
            "shell": {"layer_index": 2, "volume_fraction": 0.125, "material": _constant_material(3000.0)},
            "crust": {"layer_index": 3, "radius_fraction": 1.0, "material": _constant_material(2800.0)},
        }}
    world = construct_world(config)
    outers = _layer_outer_radii(world)
    expected_shell = (4.5e6 ** 3 + 0.125 * 6.0e6 ** 3) ** (1.0 / 3.0)
    assert math.isclose(outers[0], 2.0e6, rel_tol=1e-12)
    assert math.isclose(outers[1], 4.5e6, rel_tol=1e-12)
    assert math.isclose(outers[2], expected_shell, rel_tol=1e-9)
    assert math.isclose(outers[3], 6.0e6, rel_tol=1e-12)
    assert world.solve_eos(G_to_use=G, verbose=False)["success"]


# =====================================================================================================================
# Truncation keys and their aliases
# =====================================================================================================================
def _terrestrial_with_tides(tides_toml_keys):
    config = {"schema_version": "0.2.0", "name": "Tides", "type": "terrestrial",
              "radius_m": 6.371e6, "mass_kg": 5.972e24,
              "layers": {"mantle": {"radius_fraction": 1.0}},
              "tides": dict(tides_toml_keys)}
    return build_world(config)


@pytest.mark.parametrize("key, value, stored_key, expected", [
    # The alias used to be masked by the config default.
    ("eccentricity_trunc_lvl", 6, "eccentricity_trunc_lvl", 6),
    ("eccentricity_truncation", 6, "eccentricity_trunc_lvl", 6),
    ("obliquity_trunc_lvl", 2, "obliquity_trunc_lvl", 2),
    ("obliquity_truncation", 2, "obliquity_trunc_lvl", 2),
    ("obliquity_trunc_lvl", "off", "obliquity_trunc_lvl", 0),
    ("obliquity_trunc_lvl", "gen", "obliquity_trunc_lvl", "gen"),
])
def test_truncation_spelling_takes_effect(key, value, stored_key, expected):
    tides = _terrestrial_with_tides({key: value}).get_config_dict()["tides"]
    assert tides[stored_key] == expected


def test_both_truncation_spellings_raises():
    with pytest.raises(ValueError, match="both 'eccentricity_trunc_lvl' and its alias"):
        _terrestrial_with_tides({"eccentricity_trunc_lvl": 4, "eccentricity_truncation": 6})


# =====================================================================================================================
# World-level defaults tier ([worlds] in the configuration)
# =====================================================================================================================
def _bare_terrestrial():
    """A world config that names no world-level property, so the defaults tier decides."""
    return {"schema_version": "0.2.0", "name": "Bare", "type": "terrestrial",
            "radius_m": 6.0e6, "mass_kg": 5.0e24,
            "layers": {"mantle": {"radius_fraction": 1.0}}}


@pytest.fixture
def worlds_defaults():
    """Override a key in the configuration's [worlds] block and restore it afterwards."""
    block = TidalPy.config["worlds"]
    originals = {}

    def set_default(key, value, world_type=None):
        table = block if world_type is None else block[world_type]
        originals.setdefault((world_type, key), table.get(key))
        table[key] = value

    yield set_default
    for (world_type, key), value in originals.items():
        table = block if world_type is None else block[world_type]
        if value is None:
            table.pop(key, None)
        else:
            table[key] = value


def test_world_defaults_come_from_the_config():
    world = build_world(_bare_terrestrial())
    assert world.albedo == pytest.approx(0.3)
    assert world.emissivity == pytest.approx(1.0)


@pytest.mark.parametrize("user_albedo, expected", [
    pytest.param(None, 0.55, id="config-default"),
    pytest.param(0.11, 0.11, id="user-value-wins"),
])
def test_world_defaults_are_configurable(worlds_defaults, user_albedo, expected):
    worlds_defaults("albedo", 0.55)
    config = _bare_terrestrial()
    if user_albedo is not None:
        config["albedo"] = user_albedo
    assert build_world(config).albedo == pytest.approx(expected)


def test_star_only_defaults_do_not_reach_other_worlds(worlds_defaults):
    """[worlds.star] keys apply to stars; a terrestrial world never sees them."""
    assert build_world(_star_dict()).effective_temperature == pytest.approx(5772.0)
    worlds_defaults("effective_temperature_k", 4000.0, world_type="star")
    # A terrestrial world has no effective temperature to set, so it must still build.
    assert build_world(_bare_terrestrial()).albedo == pytest.approx(0.3)
    assert build_world(_star_dict()).effective_temperature == pytest.approx(4000.0)


# =====================================================================================================================
# Bundled worlds
# =====================================================================================================================
@pytest.mark.parametrize("name, listed", [
    ("earth_simple", True),
    ("jupiter_simple", True),
    ("sol", True),
    # Systems share the world pack directory but are not buildable worlds.
    ("sol_system", False),
])
def test_available_worlds_lists_bundled_worlds_only(name, listed):
    assert (name in available_worlds()) is listed


@pytest.mark.parametrize("config, kind", [
    ({"type": "star", "name": "Sol"}, "world"),
    ({"name": "Sol System", "worlds": {"sun": {"world": "sol"}}}, "system"),
])
def test_config_kind_tells_the_two_apart(config, kind):
    assert config_kind(config) == kind


def test_building_a_system_as_a_world_names_build_system():
    with pytest.raises(ValueError, match="system configuration.*build_system"):
        build_world("sol_system")


@pytest.mark.parametrize("world_name", ["earth_simple", "jupiter_simple", "sol"])
def test_bundled_worlds_build(world_name):
    world = build_world(world_name)
    assert isinstance(world, BaseWorld)
    assert world.name


def test_bundled_earth_simple_solves_with_its_material_tables():
    """earth_simple gives each layer a full material table."""
    world = build_world("earth_simple")
    assert world.all_materials_set is True
    assert world.solve_eos(G_to_use=G, verbose=False)["success"]


def test_unknown_bundled_name_raises():
    with pytest.raises(FileNotFoundError):
        build_world("not_a_bundled_world")


# =====================================================================================================================
# WorldPack install and data-directory-first resolution
# =====================================================================================================================
def test_install_worldpack_copies_packaged_worlds(tmp_path, monkeypatch):
    monkeypatch.setattr(worldpack, "get_worlds_dir", lambda: str(tmp_path))
    worldpack.install_worldpack()
    for name in ("earth_simple", "jupiter_simple", "sol"):
        assert (tmp_path / f"{name}.toml").is_file()


def test_install_worldpack_does_not_clobber_user_edits(tmp_path, monkeypatch):
    monkeypatch.setattr(worldpack, "get_worlds_dir", lambda: str(tmp_path))
    user_copy = tmp_path / "earth_simple.toml"
    user_copy.write_text('schema_version = "0.2.0"\nname = "edited"\n', encoding="utf-8")
    worldpack.install_worldpack()
    assert 'name = "edited"' in user_copy.read_text(encoding="utf-8")
    worldpack.install_worldpack(force=True)
    assert 'name = "edited"' not in user_copy.read_text(encoding="utf-8")


def test_build_world_prefers_data_dir_copy(tmp_path, monkeypatch):
    """A bare-name build resolves the user data-directory TOML over the packaged one."""
    monkeypatch.setattr(worldpack, "get_worlds_dir", lambda: str(tmp_path))
    edited = {
        "schema_version": "0.2.0", "name": "Edited-Earth", "type": "star",
        "radius_m": 1.0e6, "mass_kg": 1.0e24, "effective_temperature_k": 4000.0,
    }
    with open(tmp_path / "earth_simple.toml", "w") as handle:
        toml.dump(edited, handle)
    world = build_world("earth_simple")
    assert world.name == "Edited-Earth"
    assert world.world_type == "star"


def test_schema_minor_mismatch_warns_unless_forced():
    config = _terrestrial_edited(("schema_version",), "0.1.0")
    with pytest.warns(UserWarning):
        build_world(config)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        build_world(config, force=True)


def test_schema_major_mismatch_raises_unless_forced():
    config = _terrestrial_edited(("schema_version",), "1.0.0")
    with pytest.raises(ValueError, match="major versions differ"):
        build_world(config)
    assert build_world(config, force=True).num_layers == 2


# =====================================================================================================================
# Round trip: build, save, reload
# =====================================================================================================================
def test_save_and_reload_round_trip(tmp_path):
    world = build_world(_terrestrial_dict())
    path = tmp_path / "roundtrip.toml"
    world.save_to_toml(str(path))
    assert path.is_file()

    reloaded = build_world(str(path))
    assert reloaded.name == world.name
    assert reloaded.world_type == world.world_type
    assert reloaded.num_layers == world.num_layers
    assert reloaded.config["schema_version"] == "0.2.0"

    world.solve_eos(G_to_use=G, verbose=False)
    reloaded.solve_eos(G_to_use=G, verbose=False)
    assert math.isclose(world.surface_gravity_eos, reloaded.surface_gravity_eos, rel_tol=1.0e-9)


def test_save_world_to_toml_requires_toml_extension(tmp_path):
    with pytest.raises(ValueError):
        save_world_to_toml(_terrestrial_dict(), str(tmp_path / "world.txt"))


def test_save_world_to_toml_no_overwrite(tmp_path):
    config = _terrestrial_dict()
    path = str(tmp_path / "world.toml")
    save_world_to_toml(config, path)
    with pytest.raises(FileExistsError):
        save_world_to_toml(config, path, overwrite=False)
