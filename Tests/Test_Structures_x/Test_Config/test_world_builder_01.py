"""
Tests for the structures_x world builder
(``TidalPy.structures_x.configs.world_builder`` and the convenience exports on
``TidalPy.structures_x``).

Covers building each world family from a dict, the layer ``class`` / material
``type`` split, the three-tier default-resolution chain (user -> per-material
``_x`` config -> constructor), bundled-name resolution + WorldPack_x install,
physics-model wiring (verified behaviorally), error handling, ``build_world``
returning the Cython world directly, and a build -> save -> reload round-trip.

Requires the Cython extensions to be compiled first::

    uv pip install -v <repo_root>
"""

import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.structures_x import (
    build_world, construct_world, available_worlds, save_world_to_toml)
from TidalPy.structures_x.configs import world_builder, config_kind
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.worlds.gasgiant import GasGiantWorld
from TidalPy.structures_x.worlds.stellar import StarWorld


# =====================================================================================================================
# Helpers
# =====================================================================================================================
def _terrestrial_dict():
    """A two-layer terrestrial world with fully explicit mantle physics models.

    No material ``type`` is set, so every value comes from the dict itself (tier 1), then the
    ``[layers.default]`` block of the ``_x`` config (tier 2), then a constructor default (tier 3).
    The mantle lists every model it uses so the behavioral assertions do not depend on the
    default block (tested separately).
    """
    return {
        "schema_version": "0.2.0",
        "name": "TestEarth",
        "type": "terrestrial",
        "radius_m": 6.0e6,
        "mass_kg": 5.0e24,
        "spin_frequency_rad_s": 7.0e-5,
        "layers": {
            "core": {"class": "physics", "layer_index": 0, "radius_outer_m": 3.0e6, "is_tidal": False,
                     "material": {"model": "constant", "reference_density_kg_m3": 9000.0}},
            "mantle": {"class": "solidliquid", "layer_index": 1, "radius_outer_m": 6.0e6, "mass_kg": 3.0e24,
                       "is_tidal": True,
                       "material": {"model": "constant", "reference_density_kg_m3": 4000.0,
                                    "shear_modulus_static_pa": 8.0e10, "bulk_modulus_static_pa": 2.0e11,
                                    "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e21},
                                    "bulk_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e22},
                                    "partial_melt": {"model": "off"}}, "shear_rheology": {"model": "maxwell"},
                       "bulk_rheology": {"model": "elastic"}, "cooling": {"model": "convection"},
                       "radiogenics": {"model": "fixed", "fixed_heat_production_w_kg": 5.0e-12}},
        },
    }


def _single_layer_world(layer_overrides):
    """A one-layer terrestrial world; ``layer_overrides`` is merged into the layer."""
    layer = {
        "class": "solidliquid",
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


# =====================================================================================================================
# Dispatch: each world family builds the right underlying class
# =====================================================================================================================
def test_construct_terrestrial_class():
    world = construct_world(_terrestrial_dict())
    assert isinstance(world, LayeredWorld)
    assert world.num_layers == 2


def test_construct_gasgiant_class():
    config = {
        "name": "G", "type": "gasgiant", "radius_m": 7.0e7, "mass_kg": 1.9e27,
        "layers": {"env": {"class": "gas", "radius_outer_m": 7.0e7,
                           "material": {"model": "constant", "reference_density_kg_m3": 1300.0}}},
    }
    world = construct_world(config)
    assert isinstance(world, GasGiantWorld)
    assert world.num_layers == 1


def test_construct_star_class():
    config = {"name": "S", "type": "star", "radius_m": 7.0e8, "mass_kg": 2.0e30,
              "effective_temperature_k": 5772.0}
    world = construct_world(config)
    assert isinstance(world, StarWorld)
    assert math.isclose(world.effective_temperature, 5772.0)
    assert world.luminosity > 0.0


def test_layered_alias_maps_to_layered_world():
    config = _terrestrial_dict()
    config["type"] = "layered"
    world = construct_world(config)
    assert isinstance(world, LayeredWorld)


# =====================================================================================================================
# Layer ordering
# =====================================================================================================================
def test_layers_sorted_by_index_regardless_of_declaration_order():
    config = _terrestrial_dict()
    # Reverse the declaration order; layer_index must still drive the stack.
    config["layers"] = {
        "mantle": config["layers"]["mantle"],
        "core": config["layers"]["core"],
    }
    world = construct_world(config)
    # If ordering were wrong the geometry would be discontinuous and the build
    # would have raised; a successful 2-layer build confirms correct ordering.
    assert world.num_layers == 2
    result = world.solve_eos(G_to_use=G, verbose=False)
    assert result["success"]


# =====================================================================================================================
# Model wiring (behavioral, explicit configs)
# =====================================================================================================================
def test_eos_wired_and_solves():
    world = construct_world(_terrestrial_dict())
    assert world.all_eos_set is True
    result = world.solve_eos(G_to_use=G, verbose=False)
    assert result["success"] is True
    # Core (9000) is denser than mantle (4000).
    assert math.isclose(world.get_density(1.0e6), 9000.0, rel_tol=0.05)
    assert math.isclose(world.get_density(4.5e6), 4000.0, rel_tol=0.05)


def test_radiogenics_wired_produces_heating():
    world = construct_world(_terrestrial_dict())
    # The mantle's fixed radiogenics model contributes internal heating.
    assert world.calc_internal_heating(0.0) > 0.0


def test_unrecognized_model_key_names_the_table():
    """A misspelled key in a layer model table raises ValueError naming the TOML table and the key."""
    config = _terrestrial_dict()
    config["layers"]["mantle"]["material"]["partial_melt"] = {"model": "henning", "solidus": 1500.0}
    with pytest.raises(ValueError, match=r"\[layers\.mantle\.material\].*partial_melt.*'solidus'"):
        construct_world(config)


def test_rheology_and_viscosity_wired_give_complex_modulus():
    world = construct_world(_terrestrial_dict())
    world.solve_eos(G_to_use=G, verbose=False)
    # Maxwell shear rheology + finite viscosity -> complex shear modulus with a
    # non-zero imaginary (dissipative) part inside the mantle.
    mu = world.calc_complex_shear_modulus(4.5e6, 1.0e-5)
    assert np.isfinite(mu.real) and np.isfinite(mu.imag)
    assert mu.imag != 0.0


def test_missing_eos_on_a_layer_blocks_solve():
    """A layer without an EOS model blocks the solve; ``type = "none"`` keeps the builder from supplying one."""
    config = _terrestrial_dict()
    del config["layers"]["core"]["material"]
    config["layers"]["core"]["type"] = "none"
    world = construct_world(config)
    assert world.all_eos_set is False
    with pytest.raises(ValueError):
        world.solve_eos(verbose=False)


def test_typeless_layer_takes_the_default_eos():
    """A layer that names no material type takes its EOS from the ``[layers.default]`` block."""
    config = _terrestrial_dict()
    del config["layers"]["core"]["material"]
    world = construct_world(config)
    assert world.all_eos_set is True
    world.solve_eos(G_to_use=G, verbose=False)
    # The default block is a copy of mantle_rock (3500 kg/m^3).
    assert math.isclose(world.get_density(1.0e6), 3500.0, rel_tol=1e-6)


# =====================================================================================================================
# Per-material defaults (tier 2: the _x config)
# =====================================================================================================================
def test_material_defaults_filtered_to_class():
    # solidliquid mantle_rock keeps cooling/radiogenics; physics drops them.
    sl = world_builder._material_type_defaults("mantle_rock", "solidliquid")
    assert "cooling" in sl and "radiogenics" in sl
    assert "shear_rheology" in sl
    ph = world_builder._material_type_defaults("mantle_rock", "physics")
    assert "cooling" not in ph and "radiogenics" not in ph
    # The mantle_rock shear rheology default is Andrade with a zeta parameter.
    assert sl["shear_rheology"]["model"] == "andrade"
    assert "zeta" in sl["shear_rheology"]


def test_no_material_type_uses_the_default_block():
    """None selects the ``[layers.default]`` block (a copy of mantle_rock); ``"none"`` selects nothing."""
    from_none = world_builder._material_type_defaults(None, "solidliquid")
    assert from_none == world_builder._material_type_defaults("default", "solidliquid")
    assert from_none == world_builder._material_type_defaults("mantle_rock", "solidliquid")
    assert from_none["material"]["model"] == "constant"
    assert world_builder._material_type_defaults("none", "solidliquid") == {}


def test_material_type_supplies_eos_without_explicit_config():
    # A mantle_rock layer with no explicit eos still gets one from the _x config.
    world = construct_world(_single_layer_world({"type": "mantle_rock"}))
    assert world.all_eos_set is True
    result = world.solve_eos(G_to_use=G, verbose=False)
    assert result["success"]
    # mantle_rock default density is 3500 kg/m^3.
    assert math.isclose(world.get_density(3.0e6), 3500.0, rel_tol=0.05)


def test_user_value_overrides_material_default():
    # An explicit eos density (tier 1) beats the mantle_rock default (tier 2).
    world = construct_world(_single_layer_world({
        "type": "mantle_rock",
        "material": {"model": "constant", "reference_density_kg_m3": 5200.0},
    }))
    world.solve_eos(G_to_use=G, verbose=False)
    assert math.isclose(world.get_density(3.0e6), 5200.0, rel_tol=0.05)


def test_incompatible_class_and_material_type_is_safe():
    # A physics layer with an ice material type: ice's cooling/radiogenics defaults
    # are silently dropped (physics cannot hold them) but its eos still applies.
    world = construct_world(_single_layer_world({"class": "physics", "type": "ice"}))
    assert world.all_eos_set is True
    world.solve_eos(G_to_use=G, verbose=False)
    # ice default density is 1000 kg/m^3.
    assert math.isclose(world.get_density(3.0e6), 1000.0, rel_tol=0.05)


# =====================================================================================================================
# Layer geometry: derived inner radius + the three outer-radius specifiers
# =====================================================================================================================
def _layer_outer_radii(world):
    """Outer radii [m] of a built world's layers, in order (from the C++ geometry)."""
    return [layer["radius_outer_m"] for layer in world.get_config_dict()["layers"].values()]


def _two_layer_world(name, inner_spec):
    """A 6000 km world whose inner layer takes ``inner_spec``; the outer layer fills the rest, as it has to."""
    inner = {"class": "solidliquid", "layer_index": 0}
    inner.update(inner_spec)
    return {
        "name": name, "type": "terrestrial", "radius_m": 6.0e6, "mass_kg": 5.0e24,
        "layers": {
            "inner": inner,
            "outer": {"class": "solidliquid", "layer_index": 1, "radius_fraction": 1.0},
        }}


def test_radius_outer_m_spec():
    world = construct_world(_two_layer_world("O", {"radius_outer_m": 4.0e6}))
    assert math.isclose(_layer_outer_radii(world)[0], 4.0e6, rel_tol=1e-12)
    assert math.isclose(_layer_outer_radii(world)[1], 6.0e6, rel_tol=1e-12)


def test_radius_fraction_spec():
    # outer = radius_fraction * world radius (6e6) = 3e6.
    world = construct_world(_two_layer_world("F", {"radius_fraction": 0.5}))
    assert math.isclose(_layer_outer_radii(world)[0], 3.0e6, rel_tol=1e-12)


def test_volume_fraction_spec():
    # Innermost layer (inner=0): outer = (volume_fraction * R^3)^(1/3) = R * f^(1/3).
    # f = 0.125 -> outer = 0.5 * 6e6 = 3e6.
    world = construct_world(_two_layer_world("V", {"volume_fraction": 0.125}))
    assert math.isclose(_layer_outer_radii(world)[0], 3.0e6, rel_tol=1e-9)


def test_inner_radius_derived_from_previous_layer():
    # Inner radii are never specified; each layer's inner = previous outer (0 first).
    config = {
        "name": "D", "type": "terrestrial", "radius_m": 6.0e6, "mass_kg": 5.0e24,
        "layers": {
            "core": {"class": "physics", "layer_index": 0, "radius_outer_m": 2.0e6,
                     "material": {"model": "constant", "reference_density_kg_m3": 8000.0}},
            "mid": {"class": "physics", "layer_index": 1, "radius_fraction": 0.75,
                    "material": {"model": "constant", "reference_density_kg_m3": 5000.0}},
            "shell": {"class": "physics", "layer_index": 2, "volume_fraction": 0.125,
                      "material": {"model": "constant", "reference_density_kg_m3": 3000.0}},
            "crust": {"class": "physics", "layer_index": 3, "radius_fraction": 1.0,
                      "material": {"model": "constant", "reference_density_kg_m3": 2800.0}},
        }}
    world = construct_world(config)
    outers = _layer_outer_radii(world)
    # core: 2e6; mid: 0.75*6e6 = 4.5e6; shell: (4.5e6^3 + 0.125*6e6^3)^(1/3); crust: the surface.
    expected_shell = (4.5e6 ** 3 + 0.125 * 6.0e6 ** 3) ** (1.0 / 3.0)
    assert math.isclose(outers[0], 2.0e6, rel_tol=1e-12)
    assert math.isclose(outers[1], 4.5e6, rel_tol=1e-12)
    assert math.isclose(outers[2], expected_shell, rel_tol=1e-9)
    assert math.isclose(outers[3], 6.0e6, rel_tol=1e-12)
    # Geometry is continuous, so the structure solve succeeds.
    assert world.solve_eos(G_to_use=G, verbose=False)["success"]


# =====================================================================================================================
# Error handling
# =====================================================================================================================
def test_unknown_model_name_raises():
    config = _terrestrial_dict()
    config["layers"]["mantle"]["shear_rheology"] = {"model": "not_a_rheology"}
    with pytest.raises(ValueError):
        construct_world(config)


def test_construct_world_validates():
    with pytest.raises(ValueError):
        construct_world({"name": "X", "type": "terrestrial", "radius_m": 1.0,
                         "mass_kg": 1.0})  # no layers


# =====================================================================================================================
# build_world returns the Cython world directly + bundled worlds
# =====================================================================================================================
def test_build_world_returns_cython_world():
    world = build_world(_terrestrial_dict())
    # build_world returns the underlying Cython world (a BaseWorld subclass), not a wrapper.
    assert isinstance(world, LayeredWorld)
    assert world.name == "TestEarth"
    assert world.world_type == "terrestrial"
    assert world.num_layers == 2
    # The normalized build config is retained on the world.
    assert world.config["name"] == "TestEarth"
    result = world.solve_eos(G_to_use=G, verbose=False)
    assert result["success"]


def test_baseworld_build_dispatches_to_subclass():
    # BaseWorld.build is the build logic; it returns the type-appropriate subclass.
    from TidalPy.structures_x.worlds.base import BaseWorld
    star = BaseWorld.build({"name": "S", "type": "star", "radius_m": 7.0e8,
                            "mass_kg": 2.0e30, "effective_temperature_k": 5772.0})
    assert isinstance(star, StarWorld)


def _terrestrial_with_tides(tides_toml_keys):
    """Build a one-layer terrestrial world whose ``[tides]`` table holds the given keys."""
    config = {"schema_version": "0.2.0", "name": "Tides", "type": "terrestrial",
              "radius_m": 6.371e6, "mass_kg": 5.972e24,
              "layers": {"mantle": {"class": "solidliquid", "type": "mantle_rock",
                                    "radius_fraction": 1.0}},
              "tides": dict(tides_toml_keys)}
    return build_world(config)


@pytest.mark.parametrize("key", ["eccentricity_trunc_lvl", "eccentricity_truncation"])
def test_eccentricity_truncation_alias_takes_effect(key):
    """Both spellings reach the world. The alias used to be masked by the config_x default."""
    tides = _terrestrial_with_tides({key: 6}).get_config_dict()["tides"]
    assert tides["eccentricity_trunc_lvl"] == 6


@pytest.mark.parametrize("key", ["obliquity_trunc_lvl", "obliquity_truncation"])
def test_obliquity_truncation_alias_takes_effect(key):
    tides = _terrestrial_with_tides({key: 2}).get_config_dict()["tides"]
    assert tides["obliquity_trunc_lvl"] == 2


def test_obliquity_truncation_words_resolve():
    """``"off"`` is 0 and ``"gen"`` is the exact, untruncated form."""
    assert _terrestrial_with_tides({"obliquity_trunc_lvl": "off"}).get_config_dict()[
        "tides"]["obliquity_trunc_lvl"] == 0
    assert _terrestrial_with_tides({"obliquity_trunc_lvl": "gen"}).get_config_dict()[
        "tides"]["obliquity_trunc_lvl"] == 10


def test_both_truncation_spellings_raises():
    """Two spellings of one truncation leave the intended level ambiguous."""
    with pytest.raises(ValueError, match="both 'eccentricity_trunc_lvl' and its alias"):
        _terrestrial_with_tides({"eccentricity_trunc_lvl": 4, "eccentricity_truncation": 6})


# =====================================================================================================
# World-level defaults tier ([worlds] in the _x config)
# =====================================================================================================
def _bare_terrestrial():
    """A world config that names no world-level property, so the defaults tier decides."""
    return {"schema_version": "0.2.0", "name": "Bare", "type": "terrestrial",
            "radius_m": 6.0e6, "mass_kg": 5.0e24,
            "layers": {"mantle": {"class": "solidliquid", "type": "mantle_rock",
                                  "radius_fraction": 1.0}}}


@pytest.fixture
def worlds_defaults():
    """Override a key in the _x config's [worlds] block and restore it afterwards."""
    import TidalPy
    block = TidalPy.config_x["worlds"]
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
    """A world naming no albedo or emissivity picks them up from [worlds]."""
    world = build_world(_bare_terrestrial())
    assert world.albedo == pytest.approx(0.3)
    assert world.emissivity == pytest.approx(1.0)


def test_world_defaults_are_configurable(worlds_defaults):
    worlds_defaults("albedo", 0.55)
    assert build_world(_bare_terrestrial()).albedo == pytest.approx(0.55)


def test_user_world_wins_over_the_defaults(worlds_defaults):
    """Tier 1 still beats tier 2."""
    worlds_defaults("albedo", 0.55)
    config = _bare_terrestrial()
    config["albedo"] = 0.11
    assert build_world(config).albedo == pytest.approx(0.11)


def test_star_only_defaults_do_not_reach_other_worlds(worlds_defaults):
    """[worlds.star] keys apply to stars; a terrestrial world never sees them."""
    star = build_world({"schema_version": "0.2.0", "name": "S", "type": "star",
                        "radius_m": 6.957e8, "mass_kg": 1.988e30})
    assert star.effective_temperature == pytest.approx(5772.0)
    # A terrestrial world has no effective temperature to set, so building it must still work.
    worlds_defaults("effective_temperature_k", 4000.0, world_type="star")
    assert build_world(_bare_terrestrial()).albedo == pytest.approx(0.3)
    assert build_world({"schema_version": "0.2.0", "name": "S", "type": "star",
                        "radius_m": 6.957e8, "mass_kg": 1.988e30}
                       ).effective_temperature == pytest.approx(4000.0)


def test_available_worlds_lists_bundled():
    worlds = available_worlds()
    assert "earth_simple" in worlds
    assert "jupiter_simple" in worlds
    assert "sol" in worlds


def test_available_worlds_excludes_systems():
    """Systems share the world pack directory but are not buildable worlds."""
    assert "sol_system" not in available_worlds()


def test_config_kind_tells_the_two_apart():
    assert config_kind({"type": "star", "name": "Sol"}) == "world"
    assert config_kind({"name": "Sol System", "worlds": {"sun": {"world": "sol"}}}) == "system"


def test_building_a_system_as_a_world_names_build_system():
    with pytest.raises(ValueError, match="system configuration.*build_system"):
        build_world("sol_system")


@pytest.mark.parametrize("world_name", ["earth_simple", "jupiter_simple", "sol"])
def test_bundled_worlds_build(world_name):
    world = build_world(world_name)
    from TidalPy.structures_x.worlds.base import BaseWorld
    assert isinstance(world, BaseWorld)
    assert world.name


def test_bundled_earth_simple_solves_via_material_defaults():
    # earth_simple omits explicit models; both layers must still get an EOS from
    # their material types and the world must solve.
    world = build_world("earth_simple")
    assert world.all_eos_set is True
    result = world.solve_eos(G_to_use=G, verbose=False)
    assert result["success"]


def test_unknown_bundled_name_raises():
    with pytest.raises(FileNotFoundError):
        build_world("not_a_bundled_world")


# =====================================================================================================================
# WorldPack_x install + data-directory-first resolution
# =====================================================================================================================
def test_install_worldpack_x_copies_packaged_worlds(tmp_path, monkeypatch):
    from TidalPy.structures_x.configs import worldpack
    monkeypatch.setattr(worldpack, "get_worlds_x_dir", lambda: str(tmp_path))
    worldpack.install_worldpack_x()
    for name in ("earth_simple", "jupiter_simple", "sol"):
        assert (tmp_path / f"{name}.toml").is_file()


def test_install_worldpack_x_does_not_clobber_user_edits(tmp_path, monkeypatch):
    from TidalPy.structures_x.configs import worldpack
    monkeypatch.setattr(worldpack, "get_worlds_x_dir", lambda: str(tmp_path))
    user_copy = tmp_path / "earth_simple.toml"
    user_copy.write_text('schema_version = "0.2.0"\nname = "edited"\n', encoding="utf-8")
    worldpack.install_worldpack_x()
    assert 'name = "edited"' in user_copy.read_text(encoding="utf-8")
    # force=True re-copies the packaged version, discarding the edit.
    worldpack.install_worldpack_x(force=True)
    assert 'name = "edited"' not in user_copy.read_text(encoding="utf-8")


def test_build_world_prefers_data_dir_copy(tmp_path, monkeypatch):
    """A bare-name build resolves the user data-directory TOML over the packaged one."""
    import toml
    from TidalPy.structures_x.configs import worldpack
    monkeypatch.setattr(worldpack, "get_worlds_x_dir", lambda: str(tmp_path))
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
    config = _terrestrial_dict()
    config["schema_version"] = "0.1.0"   # minor difference -> warn, still builds
    with pytest.warns(UserWarning):
        build_world(config)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        build_world(config, force=True)


def test_schema_major_mismatch_raises_unless_forced():
    config = _terrestrial_dict()
    config["schema_version"] = "1.0.0"   # major difference -> refuse to build
    with pytest.raises(ValueError, match="major versions differ"):
        build_world(config)
    # force bypasses the refusal.
    world = build_world(config, force=True)
    assert world.num_layers == 2


# =====================================================================================================================
# Round-trip: build -> save -> reload
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
    assert math.isclose(world.surface_gravity_eos, reloaded.surface_gravity_eos,
                        rel_tol=1.0e-9)


def test_save_world_to_toml_requires_toml_extension(tmp_path):
    config = _terrestrial_dict()
    with pytest.raises(ValueError):
        save_world_to_toml(config, str(tmp_path / "world.txt"))


def test_save_world_to_toml_no_overwrite(tmp_path):
    config = _terrestrial_dict()
    path = str(tmp_path / "world.toml")
    save_world_to_toml(config, path)
    with pytest.raises(FileExistsError):
        save_world_to_toml(config, path, overwrite=False)
