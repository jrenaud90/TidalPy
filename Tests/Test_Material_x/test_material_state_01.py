"""Tests for the material state: ``MaterialEOSBase.calc_material_state`` and the parameters and models that feed it.

The material (a layer's EOS model) owns everything frequency-independent: the density law, the static shear law,
the static bulk modulus and viscosities, and the optional viscosity and partial-melt models.
``calc_material_state`` is the one place that maps a pressure and temperature onto all of it, and it is what the
whole-planet EOS solve evaluates as it integrates, so each stage is checked against the model it delegates to.
Nothing here knows a forcing frequency: complex moduli are the rheology's.

Requires the Cython extensions to be compiled first::

    uv pip install -v <repo_root>
"""
import math

import pytest

from TidalPy.Material_x.eos import make_material_eos
from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
from TidalPy.partial_melt_x import make_partial_melt
from TidalPy.structures_x.configs.world_builder import construct_world
from TidalPy.structures_x.layers.gas import GasLayer
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.structures_x.layers.solidliquid import SolidLiquidLayer
from TidalPy.viscosity_x import make_viscosity

_SHEAR = 6.0e10          # [Pa]
_BULK = 1.2e11           # [Pa]
_PRESSURE = 2.0e10       # [Pa]
_TEMPERATURE = 1500.0    # [K]

_VISCOSITY_CONFIG = {
    "reference_viscosity_pas": 1.0e19,
    "reference_temperature_k": 1400.0,
    "molar_activation_energy_j_mol": 3.0e5,
    "molar_activation_volume_m3_mol": 1.0e-6}


def _make_material(**kwargs):
    kwargs.setdefault("shear_modulus_static", _SHEAR)
    kwargs.setdefault("bulk_modulus_static", _BULK)
    return ConstantDensityEOS(reference_density=3300.0, **kwargs)


def test_default_material_parameters():
    material = ConstantDensityEOS()
    assert material.shear_modulus_static == 0.0
    assert material.bulk_modulus_static == 0.0
    assert math.isnan(material.shear_viscosity_static)
    assert math.isnan(material.bulk_viscosity_static)
    assert material.shear_modulus_pressure_derivative == 0.0
    assert material.shear_modulus_temperature_derivative == 0.0
    assert material.shear_modulus_reference_temperature == pytest.approx(300.0)
    assert (material.shear_viscosity_set, material.bulk_viscosity_set, material.partial_melt_set) == (
        False, False, False)


def test_thermal_constants_belong_to_the_material():
    material = ConstantDensityEOS()
    assert material.thermal_conductivity == pytest.approx(4.0)
    assert material.heat_capacity == pytest.approx(1200.0)
    assert material.thermal_expansion == 0.0

    material = ConstantDensityEOS(thermal_conductivity=3.0, heat_capacity=1000.0, thermal_expansion=2.0e-5)
    assert material.thermal_conductivity == pytest.approx(3.0)
    assert material.heat_capacity == pytest.approx(1000.0)
    # kappa = k / (rho c_p), at whatever density it is asked about.
    assert material.calc_thermal_diffusivity(4000.0) == pytest.approx(3.0 / (4000.0 * 1000.0))
    assert math.isnan(material.calc_thermal_diffusivity(0.0))

    config = material.get_config_dict()
    assert config["thermal_conductivity_w_mk"] == pytest.approx(3.0)
    assert config["heat_capacity_j_kgk"] == pytest.approx(1000.0)
    assert config["thermal_expansion_1_k"] == pytest.approx(2.0e-5)
    rebuilt = make_material_eos(config.pop("model"), config)
    assert rebuilt.thermal_conductivity == pytest.approx(3.0)
    assert rebuilt.heat_capacity == pytest.approx(1000.0)


def test_one_expansivity_serves_the_density_law_only_when_asked():
    """The material has a single alpha. thermal_density (a layer's use_thermal_eos) decides if the density sees T."""
    material = ConstantDensityEOS(reference_density=4000.0, thermal_expansion=3.0e-5)
    athermal = material.calc_material_state(0.0, 2000.0, thermal_density=False)
    thermal = material.calc_material_state(0.0, 2000.0, thermal_density=True)
    assert athermal["density"] == pytest.approx(4000.0)
    assert thermal["density"] == pytest.approx(4000.0 * math.exp(-3.0e-5 * (2000.0 - 300.0)))


def test_unknown_material_keyword_is_rejected():
    with pytest.raises(TypeError, match="shear_modulus"):
        ConstantDensityEOS(shear_modulus=1.0e10)


def test_bare_material_returns_its_constants():
    """No models attached: the static constants, no melt, and the density of the law."""
    state = _make_material().calc_material_state(_PRESSURE, _TEMPERATURE)
    assert set(state) == {
        "density", "melt_fraction", "shear_modulus", "bulk_modulus", "shear_viscosity", "bulk_viscosity"}
    assert state["density"] == pytest.approx(3300.0)
    assert state["melt_fraction"] == 0.0
    assert state["shear_modulus"] == pytest.approx(_SHEAR)
    assert state["bulk_modulus"] == pytest.approx(_BULK)
    assert math.isnan(state["shear_viscosity"])


def test_static_viscosity_is_the_fallback_without_a_viscosity_model():
    material = _make_material(shear_viscosity_static=1.0e21, bulk_viscosity_static=2.0e21)
    state = material.calc_material_state(_PRESSURE, _TEMPERATURE)
    assert state["shear_viscosity"] == pytest.approx(1.0e21)
    assert state["bulk_viscosity"] == pytest.approx(2.0e21)


def test_static_constants_are_settable():
    material = _make_material()
    material.shear_modulus_static = 4.0e10
    material.shear_viscosity_static = 3.0e20
    state = material.calc_material_state(0.0, 300.0)
    assert state["shear_modulus"] == pytest.approx(4.0e10)
    assert state["shear_viscosity"] == pytest.approx(3.0e20)


@pytest.mark.parametrize("pressure,temperature", [(0.0, 300.0), (2.0e10, 1500.0), (1.0e11, 3000.0)])
def test_linear_shear_law(pressure, temperature):
    material = _make_material(
        shear_modulus_pressure_derivative=1.5,
        shear_modulus_temperature_derivative=-1.0e7,
        shear_modulus_reference_temperature=300.0)
    expected = _SHEAR + 1.5 * pressure - 1.0e7 * (temperature - 300.0)
    assert material.calc_material_state(pressure, temperature)["shear_modulus"] == pytest.approx(
        expected, rel=1e-13)


def test_shear_law_is_floored_at_the_minimum_modulus():
    material = _make_material(shear_modulus_temperature_derivative=-1.0e9)
    shear = material.calc_material_state(0.0, 5000.0)["shear_modulus"]
    assert 0.0 < shear < 1.0e-3 * _SHEAR


def test_viscosity_model_is_evaluated_at_the_temperature_and_pressure():
    material = _make_material()
    material.set_shear_viscosity(make_viscosity("reference", _VISCOSITY_CONFIG))
    assert material.shear_viscosity_set is True
    expected = make_viscosity("reference", _VISCOSITY_CONFIG).calc_viscosity(_TEMPERATURE, _PRESSURE)
    state = material.calc_material_state(_PRESSURE, _TEMPERATURE)
    assert state["shear_viscosity"] == pytest.approx(expected, rel=1e-13)
    # The cold limit is rigid, so a layer left at 0 K still solves as an elastic body.
    assert math.isinf(material.calc_material_state(_PRESSURE, 0.0)["shear_viscosity"])


def test_partial_melt_weakens_the_shear_pair_and_reports_the_melt_fraction():
    melt_config = {"solidus_k": 1600.0, "liquidus_k": 2000.0}
    material = _make_material()
    material.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e20}))
    material.set_partial_melt(make_partial_melt("henning", melt_config))

    solid = material.calc_material_state(_PRESSURE, 1500.0)
    assert solid["melt_fraction"] == 0.0
    assert solid["shear_modulus"] == pytest.approx(_SHEAR)

    molten = material.calc_material_state(_PRESSURE, 1900.0)
    reference = make_partial_melt("henning", melt_config).calc_partial_melt(1900.0, 1.0e20, _SHEAR, 1.0e20)
    assert molten["melt_fraction"] == pytest.approx(0.75)
    assert molten["shear_viscosity"] == pytest.approx(reference[1], rel=1e-13)
    assert molten["shear_modulus"] == pytest.approx(reference[2], rel=1e-13)
    assert molten["shear_modulus"] < _SHEAR


def test_pressure_law_bulk_modulus_takes_precedence_over_the_constant():
    eos_config = {
        "reference_density_kg_m3": 3500.0,
        "reference_bulk_modulus_pa": 1.3e11,
        "bulk_modulus_derivative": 4.5,
        "thermal_expansion_1_k": 3.0e-5,
        "bulk_modulus_static_pa": _BULK}
    material = make_material_eos("bm", eos_config)

    # thermal_density=False is what a layer without use_thermal_eos asks for: the density law is athermal.
    athermal = material.calc_material_state(_PRESSURE, _TEMPERATURE, thermal_density=False)
    assert athermal["density"] == pytest.approx(material.calc_density(_PRESSURE), rel=1e-13)
    assert athermal["bulk_modulus"] == pytest.approx(material.calc_bulk_modulus(_PRESSURE), rel=1e-13)
    assert athermal["bulk_modulus"] != pytest.approx(_BULK)

    thermal = material.calc_material_state(_PRESSURE, _TEMPERATURE)
    assert thermal["density"] == pytest.approx(material.calc_density(_PRESSURE, _TEMPERATURE), rel=1e-13)
    assert thermal["bulk_modulus"] == pytest.approx(
        material.calc_bulk_modulus(_PRESSURE, _TEMPERATURE), rel=1e-13)
    assert thermal["density"] < athermal["density"]


def test_interpolated_tables_take_precedence():
    material = make_material_eos("interpolate", {
        "radius_m": [0.0, 1.0e6],
        "density_kg_m3": [5000.0, 4000.0],
        "shear_modulus_pa": [8.0e10, 4.0e10],
        "shear_viscosity_pas": [1.0e22, 1.0e20],
        "shear_modulus_static_pa": 1.0e9,
        "bulk_modulus_static_pa": _BULK})
    state = material.calc_material_state(_PRESSURE, _TEMPERATURE, radius=5.0e5)
    assert state["density"] == pytest.approx(4500.0)
    assert state["shear_modulus"] == pytest.approx(6.0e10)      # the table, not the 1e9 constant
    assert state["shear_viscosity"] == pytest.approx(5.05e21)
    assert state["bulk_modulus"] == pytest.approx(_BULK)        # no bulk table: the constant
    assert material.get_tabulated_shear_modulus(5.0e5) == pytest.approx(6.0e10)
    assert math.isnan(material.get_tabulated_bulk_modulus(5.0e5))


def test_factory_builds_the_nested_models():
    material = make_material_eos("constant", {
        "reference_density_kg_m3": 3300.0,
        "shear_modulus_static_pa": _SHEAR,
        "shear_viscosity_static_pas": 1.0e21,
        "bulk_viscosity_static_pas": 2.0e21,   # set so the dict comparison below holds no NaN
        "shear_viscosity": dict(model="reference", **_VISCOSITY_CONFIG),
        "partial_melt": {"model": "henning", "solidus_k": 1600.0, "liquidus_k": 2000.0}})
    assert (material.shear_viscosity_set, material.bulk_viscosity_set, material.partial_melt_set) == (
        True, False, True)
    config = material.get_config_dict()
    assert config["shear_modulus_static_pa"] == pytest.approx(_SHEAR)
    assert config["shear_viscosity"]["model"] == "reference"
    assert config["partial_melt"]["solidus_k"] == pytest.approx(1600.0)
    assert "bulk_viscosity" not in config
    # The config dict rebuilds the same material.
    rebuilt = make_material_eos(config.pop("model"), config)
    assert rebuilt.calc_material_state(_PRESSURE, 1900.0) == material.calc_material_state(_PRESSURE, 1900.0)


def test_nested_model_table_needs_a_model_key():
    with pytest.raises(ValueError, match="model"):
        make_material_eos("constant", {"shear_viscosity": {"reference_viscosity_pas": 1.0e20}})


# =====================================================================================================================
# The layer hands viscosity and partial-melt models to its material
# =====================================================================================================================
def test_layer_helpers_put_the_models_on_the_material():
    layer = PhysicsLayer("mantle", 0, 0.0, 1.0e6, 1.0e22)
    with pytest.raises(ValueError, match="EOS"):
        layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e20}))
    layer.set_eos(_make_material())
    assert layer.shear_viscosity_set is False
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e20}))
    layer.set_partial_melt(make_partial_melt("henning"))
    assert (layer.shear_viscosity_set, layer.bulk_viscosity_set, layer.partial_melt_set) == (True, False, True)
    material = layer.get_config_dict()["material"]
    assert material["shear_viscosity"]["reference_viscosity_pas"] == pytest.approx(1.0e20)
    assert material["partial_melt"]["model"] == "henning"
    assert layer.shear_modulus_static == pytest.approx(_SHEAR)


@pytest.mark.parametrize("layer_class", [PhysicsLayer, SolidLiquidLayer, GasLayer])
def test_material_survives_config_and_binary_roundtrips(layer_class, tmp_path):
    layer = layer_class("shell", 0, 0.0, 1.0e6, 1.0e22, temperature=1700.0, use_thermal_eos=True)
    layer.set_eos(_make_material(
        shear_viscosity_static=5.0e20,
        bulk_viscosity_static=6.0e20,   # set so the dict comparison below holds no NaN
        shear_modulus_pressure_derivative=1.4,
        shear_modulus_temperature_derivative=-8.0e6,
        shear_modulus_reference_temperature=1600.0,
        thermal_conductivity=3.5, heat_capacity=900.0))
    layer.set_shear_viscosity(make_viscosity("reference", _VISCOSITY_CONFIG))
    layer.set_partial_melt(make_partial_melt("henning", {"solidus_k": 1600.0, "liquidus_k": 2000.0}))

    cfg = layer.get_config_dict()
    assert cfg["temperature_k"] == pytest.approx(1700.0)
    assert cfg["use_thermal_eos"] is True
    material = cfg["material"]
    assert material["shear_modulus_static_pa"] == pytest.approx(_SHEAR)
    assert material["shear_modulus_pressure_derivative"] == pytest.approx(1.4)
    assert material["shear_modulus_temperature_derivative_pa_k"] == pytest.approx(-8.0e6)
    assert material["shear_modulus_reference_temperature_k"] == pytest.approx(1600.0)
    for moved in ("shear_modulus_static_pa", "shear_viscosity", "partial_melt", "eos"):
        assert moved not in cfg

    path = tmp_path / "layer.tpyb"
    layer.save_binary(str(path))
    loaded = layer_class("placeholder", 0, 0.0, 1.0, 1.0)
    loaded.load_binary(str(path))
    assert loaded.temperature == pytest.approx(1700.0)
    assert loaded.use_thermal_eos is True
    assert loaded.get_config_dict()["material"] == material
    assert (loaded.shear_viscosity_set, loaded.partial_melt_set) == (True, True)
    # The thermal constants ride in the same record.
    assert material["thermal_conductivity_w_mk"] == pytest.approx(3.5)
    assert material["heat_capacity_j_kgk"] == pytest.approx(900.0)


def test_material_table_builds_through_the_world_builder():
    config = {
        "schema_version": "0.2.0",
        "name": "material_keys",
        "type": "terrestrial",
        "radius_m": 3.0e6,
        "mass_kg": 5.0e23,
        "layers": {
            "mantle": {
                "class": "solidliquid",
                "type": "mantle_rock",
                "radius_fraction": 1.0,
                "temperature_k": 1650.0,
                "use_thermal_eos": True,
                "material": {
                    "shear_modulus_pressure_derivative": 1.4,
                    "shear_modulus_temperature_derivative_pa_k": -8.0e6,
                    "shear_modulus_reference_temperature_k": 1600.0,
                    # One key of a nested default table is overridden; the rest of the table stands.
                    "shear_viscosity": {"reference_viscosity_pas": 3.0e21},
                },
            },
        },
    }
    world = construct_world(config)
    assert world.mantle.temperature == pytest.approx(1650.0)
    assert world.mantle.use_thermal_eos is True
    material = world.mantle.get_config_dict()["material"]
    assert material["shear_modulus_temperature_derivative_pa_k"] == pytest.approx(-8.0e6)
    assert material["shear_viscosity"]["reference_viscosity_pas"] == pytest.approx(3.0e21)
    assert material["shear_viscosity"]["model"] == "reference"          # from the mantle_rock defaults
    assert material["partial_melt"]["model"] == "henning"               # likewise
    rebuilt = construct_world(world.get_config_dict())
    assert rebuilt.mantle.get_config_dict()["material"] == material


@pytest.mark.parametrize("key, where", [
    ("shear_modulus_static_pa", "material"), ("eos", "material"), ("shear_viscosity", "material.shear_viscosity")])
def test_keys_left_on_the_layer_say_where_they_moved(key, where):
    layer = {"class": "physics", "radius_fraction": 1.0}
    layer[key] = 6.0e10 if key.endswith("_pa") else {"model": "constant"}
    config = {"schema_version": "0.2.0", "name": "moved", "type": "terrestrial", "radius_m": 1.0e6,
              "mass_kg": 1.0e22, "layers": {"mantle": layer}}
    with pytest.raises(ValueError, match=r"layers\.mantle\." + where.replace(".", r"\.")):
        construct_world(config)
