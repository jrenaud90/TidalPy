"""MatPack (TidalPy.Material.matpack): every bundled material builds and behaves, presets and overrides resolve, and the
data-directory copies are installed, preferred, and checked like the WorldPack's."""

import os
import warnings

import pytest
import toml

import TidalPy
from TidalPy.schema import SCHEMA_VERSION
from TidalPy.Material import (
    CATEGORIES,
    Material,
    available_materials,
    load_material,
    make_material,
    material_config,
    material_info,
)
from TidalPy.Material import matpack

# Listed from the package, not through available_materials: collection runs before the session fixture moves the
# data directory, and listing through the pack would install into the user's.
_ALL = sorted(entry[:-len(".toml")] for entry in os.listdir(matpack.PACKAGED_MATPACK_DIR) if entry.endswith(".toml"))


@pytest.fixture()
def data_dir(tmp_path, monkeypatch):
    """A private materials data directory holding fresh copies of the packaged files."""
    monkeypatch.setattr(matpack, "get_materials_dir", lambda: str(tmp_path))
    monkeypatch.setattr(matpack.MAT_PACK, "p_warned_stale_copies", set())
    matpack.install_matpack()
    return tmp_path


def _write_material(directory, name, table):
    with open(os.path.join(str(directory), name + ".toml"), "w", encoding="utf-8", newline="\n") as file:
        toml.dump({"schema_version": SCHEMA_VERSION, **table}, file)


# =====================================================================================================================
# The catalog
# =====================================================================================================================
def test_catalog_has_every_category_and_the_simplified_materials():
    assert len(_ALL) >= 29
    for category in CATEGORIES:
        assert available_materials(category)
    assert set(available_materials("simplified")) == {
        "simple_rock", "simple_ice", "simple_iron_core", "simple_liquid_iron", "simple_water", "simple_gas"}
    assert sorted(sum((available_materials(category) for category in CATEGORIES), [])) == _ALL


@pytest.mark.parametrize("name", _ALL)
def test_every_material_builds_and_round_trips(name, tmp_path):
    info = material_info(name)
    assert info["description"] and info["category"] in CATEGORIES
    material = load_material(name)
    config = material.get_config_dict()
    assert make_material(config).get_config_dict() == config
    path = str(tmp_path / "material.tpyb")
    material.save_binary(path)
    loaded = Material()
    loaded.load_binary(path)
    assert loaded.get_config_dict() == config


@pytest.mark.parametrize("name", _ALL)
def test_every_material_has_a_physical_state(name):
    """At 100 MPa (a polytrope has no density at zero pressure) and 250 K, under every switch."""
    material = load_material(name)
    for switches in ({}, dict(use_thermal_expansion=True, use_melting=True, use_pressure_melting=True,
                              use_melt_density=True)):
        state = material.calc_state(1.0e8, 250.0, **switches)
        assert state["density"] > 0.0
        assert state["thermal_conductivity"] > 0.0 and state["heat_capacity"] > 0.0
        assert state["bulk_modulus"] >= 0.0
        assert state["shear_modulus"] >= 0.0
        if state["phase"] == "solid":
            assert state["shear_modulus"] > 0.0
            assert state["shear_viscosity"] > 0.0


@pytest.mark.parametrize("name", _ALL)
def test_melting_curves_are_positive_and_do_not_cross(name):
    material = load_material(name)
    if not material.can_melt:
        pytest.skip("cannot melt")
    for pressure in (0.0, 1.0e8, 1.0e9, 1.0e10, 5.0e10):
        solidus, liquidus = material.calc_melting_range(pressure)
        assert solidus > 0.0, (pressure, solidus)
        assert liquidus >= solidus - 1.0e-9, (pressure, solidus, liquidus)


def test_liquid_materials_are_liquid():
    liquid_only = [name for name in _ALL if load_material(name).is_liquid_only]
    assert {"water", "brine", "liquid_iron", "simple_water", "simple_gas", "h2_he_molecular"} <= set(liquid_only)
    for name in liquid_only:
        state = load_material(name).calc_state(1.0e8, 300.0, use_melting=True)
        assert state["phase"] == "liquid" and state["melt_fraction"] == 1.0 and state["shear_modulus"] == 0.0


def test_values_match_their_sources():
    peridotite = load_material("peridotite")
    assert peridotite.calc_melting_range(0.0) == (1661.2, 1982.1)
    # The two Monteux et al. (2016) branches join at 20 GPa to within a few kelvin.
    solidus_low = 1661.2 * (1.0 + 2.0e10 / 1.336e9) ** (1.0 / 7.437)
    assert peridotite.calc_melting_range(2.0e10)[0] == pytest.approx(solidus_low, abs=3.0)
    ice = load_material("ice_ih")
    assert ice.calc_melting_range(0.0)[0] == 273.16
    assert ice.calc_melting_range(1.0e8)[0] < 273.16
    # Past the ice Ih, ice III, and liquid triple point (IAPWS R14-08: 251.165 K at 208.566 MPa) the curve holds.
    assert ice.calc_melting_range(2.08566e8)[0] == pytest.approx(251.165, abs=0.2)
    assert ice.calc_melting_range(1.0e9)[0] == ice.calc_melting_range(2.08566e8)[0]
    water = load_material("water").calc_state(0.0, 273.15)
    assert water["density"] == pytest.approx(999.84)
    # K_S / K_T = 1 + alpha gamma T reproduces Anderson and Ahrens (1994).
    iron = load_material("liquid_iron").calc_state(0.0, 1811.0, use_thermal_expansion=True)
    assert iron["adiabatic_bulk_modulus"] == pytest.approx(1.097e11, rel=2.0e-3)
    assert iron["density"] == pytest.approx(7019.0)


# A point inside each solid material's range, below its solidus: (pressure [Pa], temperature [K]).
_IN_RANGE = {
    "simple_rock": (3.0e9, 1600.0), "simple_ice": (1.0e7, 250.0), "simple_iron_core": (2.0e10, 2000.0),
    "peridotite": (3.0e9, 1500.0), "lower_mantle": (6.0e10, 2500.0), "olivine": (3.0e9, 1600.0),
    "basalt": (5.0e8, 1100.0), "felsic_crust": (3.0e8, 800.0), "chondrite": (1.0e9, 1200.0),
    "serpentinite": (1.0e9, 700.0), "iron": (1.0e11, 3000.0), "iron_sulfide": (5.0e9, 1200.0),
    "ice_ih": (5.0e7, 250.0), "ice_iii": (2.5e8, 245.0), "ice_v": (5.0e8, 255.0), "ice_vi": (1.5e9, 300.0),
    "ice_vii": (1.0e10, 550.0), "ammonia_water": (5.0e7, 170.0), "methane_clathrate": (1.0e7, 260.0),
    "nitrogen_ice": (1.0e5, 50.0)}


@pytest.mark.parametrize("name", sorted(_IN_RANGE))
def test_solids_are_solid_inside_their_range(name):
    pressure, temperature = _IN_RANGE[name]
    state = load_material(name).calc_state(
        pressure, temperature, use_thermal_expansion=True, use_melting=True, use_pressure_melting=True)
    assert state["phase"] == "solid"
    assert state["shear_modulus"] > 1.0e8
    assert state["shear_viscosity"] > 1.0e9


@pytest.mark.parametrize("name, pressures", [
    ("ice_vi", (6.4e8, 1.0e9, 2.0e9, 2.2e9)),
    ("ice_vii", (2.3e9, 5.0e9, 1.0e10, 2.0e10)),
])
def test_high_pressure_ice_viscosity_along_its_melting_curve(name, pressures):
    """Melting-point viscosities of 1e13 to 1e17 Pa s are used for the high-pressure ices."""
    material = load_material(name)
    for pressure in pressures:
        melting_temperature = material.calc_melting_range(pressure)[0]
        viscosity = material.calc_state(pressure, melting_temperature - 1.0)["shear_viscosity"]
        assert 5.0e12 < viscosity < 1.0e17, (pressure, viscosity)


def test_peridotite_branches_meet_and_reach_the_published_cmb_solidus():
    """The two branches of each Monteux et al. (2016) fit join at 20 GPa, and the solidus is about 4150 K at the
    core-mantle boundary (135 GPa; Andrault et al. 2011), below the liquidus."""
    peridotite = load_material("peridotite")
    below = peridotite.calc_melting_range(20.0e9 * (1.0 - 1.0e-9))
    above = peridotite.calc_melting_range(20.0e9 * (1.0 + 1.0e-9))
    assert below == pytest.approx(above, rel=2.0e-3)
    solidus, liquidus = peridotite.calc_melting_range(135.0e9)
    assert solidus == pytest.approx(4150.0, rel=0.01)
    assert liquidus > solidus


def test_phase_presets_share_one_phase():
    assert load_material("ice_ih").liquid.get_config_dict() == load_material("water").liquid.get_config_dict()
    assert load_material("iron").liquid.get_config_dict() == load_material("liquid_iron").liquid.get_config_dict()


# =====================================================================================================================
# Presets and overrides
# =====================================================================================================================
def test_overrides_change_only_what_they_name():
    base = material_config("peridotite")
    changed = material_config("peridotite", {"latent_heat_j_kg": 3.0e5,
                                             "solid": {"eos": {"reference_density_kg_m3": 3300.0}}})
    assert changed["latent_heat_j_kg"] == 3.0e5
    assert changed["solid"]["eos"]["reference_density_kg_m3"] == 3300.0
    changed["solid"]["eos"]["reference_density_kg_m3"] = base["solid"]["eos"]["reference_density_kg_m3"]
    changed["latent_heat_j_kg"] = base["latent_heat_j_kg"]
    assert changed == base


def test_an_alias_of_the_same_model_merges():
    config = material_config("peridotite", {"solid": {"eos": {"model": "bm", "bulk_modulus_derivative": 4.0}}})
    assert config["solid"]["eos"]["reference_bulk_modulus_pa"] == 1.25e11
    assert config["solid"]["eos"]["bulk_modulus_derivative"] == 4.0


def test_another_model_replaces_the_table():
    material = load_material("peridotite", solid={"shear_viscosity": {"model": "constant",
                                                                     "reference_viscosity_pas": 1.0e20}})
    assert material.solid.shear_viscosity.get_config_dict() == {"model": "constant", "reference_viscosity_pas": 1.0e20}


def test_none_removes_a_slot():
    material = load_material("peridotite", liquid=None, melting=None, latent_heat_j_kg=0.0)
    assert not material.can_melt and material.liquid is None
    with pytest.raises(ValueError, match="MatPack material 'peridotite'.*not both a solid"):
        load_material("peridotite", liquid=None)


def test_a_preset_table_and_a_phase_preset():
    salty = load_material({
        "preset": "ice_ih",
        "liquid": {"preset": "brine"},
        "melting": {"solidus": {"model": "constant", "temperature_k": 251.9},
                    "liquidus": {"model": "constant", "temperature_k": 271.2}}})
    assert salty.liquid.get_config_dict() == load_material("brine").liquid.get_config_dict()
    state = salty.calc_state(1.0e7, 260.0, use_melting=True)
    assert state["phase"] == "partial"
    assert state["melt_fraction"] == pytest.approx((260.0 - 251.9) / (271.2 - 251.9))
    with pytest.raises(ValueError, match="has no 'solid' phase"):
        material_config({"preset": "peridotite", "solid": {"preset": "water"}})


def test_a_material_keeps_a_phase():
    with pytest.raises(ValueError, match="neither a 'solid' nor a 'liquid'"):
        load_material("water", liquid=None)
    with pytest.raises(ValueError, match="neither a 'solid' nor a 'liquid'"):
        load_material("water").replace(liquid=None)


def test_a_preset_only_names_a_material_a_phase_or_a_melting_table():
    with pytest.raises(ValueError, match="'solidus' table .* names a 'preset'"):
        material_config({"preset": "peridotite", "melting": {"solidus": {"preset": "basalt"}}})
    with pytest.raises(ValueError, match="'eos' table .* names a 'preset'"):
        material_config({"preset": "peridotite", "solid": {"eos": {"preset": "basalt"}}})


def test_a_melting_table_names_a_preset():
    """A melting table may start from another material's, with its own keys merged over it."""
    config = material_config({"preset": "simple_rock", "liquid": {"preset": "peridotite"},
                              "melting": {"preset": "peridotite", "weakening": {"model": "spohn"}}})
    peridotite = material_config("peridotite")
    assert config["melting"]["solidus"] == peridotite["melting"]["solidus"]
    assert config["melting"]["liquidus"] == peridotite["melting"]["liquidus"]
    assert config["melting"]["weakening"]["model"] == "spohn"
    with pytest.raises(ValueError, match="has no 'melting' table"):
        material_config({"preset": "simple_rock", "melting": {"preset": "simple_rock"}})


@pytest.mark.parametrize("override", [{"latent_heat": 1.0}, {"latent_heat_j_kg": 1.0}])
def test_an_override_wins_under_either_spelling(override):
    """An override given by argument name replaces the preset's value given by config key, and the other way round."""
    assert load_material("peridotite", **override).latent_heat == 1.0
    assert material_config({"preset": "peridotite", **override})["latent_heat_j_kg"] == 1.0
    nested = material_config({"preset": "peridotite", "solid": {"eos": {"reference_density": 1000.0}}})
    assert nested["solid"]["eos"]["reference_density_kg_m3"] == 1000.0
    assert "reference_density" not in nested["solid"]["eos"]


@pytest.mark.parametrize("name", ["../MatPack/water", "MatPack/water", "..", ""])
def test_a_name_holds_no_path(name):
    with pytest.raises(ValueError, match="holds no path|no MatPack material"):
        load_material(name)


def test_metadata_is_not_part_of_the_material():
    assert not set(matpack.METADATA_KEYS) & set(material_config("ice_ih"))


@pytest.mark.parametrize("source, overrides, error, message", [
    ("peridotit", None, ValueError, "did you mean 'peridotite'"),
    (42, None, TypeError, "a name or a table"),
    ("peridotite", {"preset": "basalt"}, ValueError, "cannot name a 'preset'"),
])
def test_bad_sources(source, overrides, error, message):
    with pytest.raises(error, match=message):
        material_config(source, overrides)


def test_a_rejected_override_names_the_material():
    with pytest.raises(ValueError, match="MatPack material 'peridotite'.*did you mean 'reference_density_kg_m3'"):
        load_material("peridotite", solid={"eos": {"reference_densty_kg_m3": 3300.0}})


def test_a_preset_cycle_is_an_error(data_dir):
    _write_material(data_dir, "loop_a", {"preset": "loop_b"})
    _write_material(data_dir, "loop_b", {"preset": "loop_a"})
    with pytest.raises(ValueError, match="loop_a -> loop_b -> loop_a"):
        load_material("loop_a")


def test_a_user_material_can_start_from_a_preset(data_dir):
    _write_material(data_dir, "io_mantle", {
        "description": "Io's mantle.", "category": "rocky", "preset": "peridotite",
        "solid": {"shear_rheology": {"model": "andrade", "alpha": 0.2}}})
    assert "io_mantle" in available_materials("rocky")
    assert load_material("io_mantle").solid.shear_rheology.alpha == 0.2


# =====================================================================================================================
# The data directory
# =====================================================================================================================
def test_installed_copies_are_preferred(data_dir):
    assert os.path.isfile(os.path.join(str(data_dir), "peridotite.toml"))
    assert material_info("peridotite")["path"] == os.path.join(str(data_dir), "peridotite.toml")


def test_an_edited_copy_is_used_and_warned_about_once(data_dir):
    copy_path = os.path.join(str(data_dir), "simple_rock.toml")
    with open(copy_path, encoding="utf-8") as file:
        text = file.read()
    with open(copy_path, "w", encoding="utf-8", newline="\n") as file:
        file.write(text.replace("reference_density_kg_m3 = 3300.0", "reference_density_kg_m3 = 3100.0"))
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        assert load_material("simple_rock").calc_state(0.0)["density"] == 3100.0
        load_material("simple_rock")
    stale = [entry for entry in record if "differs from the one packaged" in str(entry.message)]
    assert len(stale) == 1
    assert "stale_matpack_copy" in str(stale[0].message)
    assert "install_matpack(force=True)" in str(stale[0].message)


def test_without_a_data_directory_the_package_is_read(monkeypatch):
    monkeypatch.setattr(matpack, "get_materials_dir", lambda: None)
    assert matpack.install_matpack() is None
    assert os.path.dirname(material_info("water")["path"]) == matpack.PACKAGED_MATPACK_DIR
    assert "water" in available_materials()


def test_files_are_lf_and_parse(data_dir):
    for name in _ALL:
        path = os.path.join(matpack.PACKAGED_MATPACK_DIR, name + ".toml")
        with open(path, "rb") as file:
            assert b"\r" not in file.read()
        assert toml.load(path)["schema_version"] == SCHEMA_VERSION


def test_names_match_without_regard_to_case(data_dir):
    _write_material(data_dir, "Io_Mantle", {"description": "Io.", "category": "rocky", "preset": "simple_rock"})
    assert "io_mantle" in available_materials()
    assert load_material("io_mantle").get_config_dict() == load_material("Io_Mantle").get_config_dict()
    assert load_material("WATER").is_liquid_only


def _switched_off(monkeypatch, name):
    config = dict(TidalPy.config)
    config["warnings"] = {name: False}
    monkeypatch.setattr(TidalPy, "config", config)


def test_the_stale_copy_warning_can_be_switched_off_and_force_restores(data_dir, monkeypatch):
    copy_path = os.path.join(str(data_dir), "simple_ice.toml")
    with open(copy_path, "a", encoding="utf-8", newline="\n") as file:
        file.write("\n# A local edit.\n")
    _switched_off(monkeypatch, "stale_matpack_copy")
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        load_material("simple_ice")
    assert not [entry for entry in record if "differs from the one packaged" in str(entry.message)]
    matpack.install_matpack(force=True)
    packaged_path = os.path.join(matpack.PACKAGED_MATPACK_DIR, "simple_ice.toml")
    with open(copy_path, "rb") as file, open(packaged_path, "rb") as packaged:
        assert file.read() == packaged.read()


def test_a_broken_user_file_is_skipped_in_listings(data_dir):
    with open(os.path.join(str(data_dir), "broken.toml"), "w", encoding="utf-8", newline="\n") as file:
        file.write("this is = = not toml\n")
    with pytest.warns(UserWarning, match="'broken' could not be read"):
        icy = available_materials("icy")
    assert "ice_ih" in icy and "broken" not in icy
    assert "broken" in available_materials()
    with pytest.raises(ValueError, match="could not parse"):
        load_material("broken")
