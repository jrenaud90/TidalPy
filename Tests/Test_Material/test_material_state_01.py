"""Material parameters and ``calc_material_state``: density, shear and bulk laws, viscosity and partial-melt models."""
import math

import pytest

from TidalPy.Material.eos import make_material_eos
from TidalPy.Material.eos.material_eos import ConstantDensityEOS
from TidalPy.PartialMelt import make_partial_melt
from TidalPy.Structures.configs.world_builder import construct_world
from TidalPy.Structures.layers.gas import GasLayer
from TidalPy.Structures.layers.base import BaseLayer
from TidalPy.Structures.layers.solidliquid import SolidLiquidLayer
from TidalPy.Viscosity import make_viscosity

_SHEAR = 6.0e10          # [Pa]
_BULK = 1.2e11           # [Pa]
_PRESSURE = 2.0e10       # [Pa]
_TEMPERATURE = 1500.0    # [K]
_MELT_CONFIG = {"solidus_k": 1600.0, "liquidus_k": 2000.0}

_VISCOSITY_CONFIG = {
    "reference_viscosity_pas": 1.0e19,
    "reference_temperature_k": 1400.0,
    "molar_activation_energy_j_mol": 3.0e5,
    "molar_activation_volume_m3_mol": 1.0e-6}


def _make_material(**kwargs):
    kwargs.setdefault("shear_modulus_static", _SHEAR)
    kwargs.setdefault("bulk_modulus_static", _BULK)
    return ConstantDensityEOS(reference_density=3300.0, **kwargs)


def _melting_material():
    """A material with a constant 1e20 Pa s viscosity and a Henning partial-melt model."""
    material = _make_material()
    material.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e20}))
    material.set_partial_melt(make_partial_melt("henning", _MELT_CONFIG))
    return material


def _layer(
        layer_class,
        name,
        radius_outer=1.0e6,
        mass=1.0e22,
        **kwargs):
    return layer_class(
        name,
        0,
        0.0,
        radius_outer,
        mass,
        **kwargs)


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
    """Conductivity, heat capacity, and expansivity live on the material, set its diffusivity, and round trip."""
    material = ConstantDensityEOS()
    assert material.thermal_conductivity == pytest.approx(4.0)
    assert material.heat_capacity == pytest.approx(1200.0)
    assert material.thermal_expansion == 0.0

    material = ConstantDensityEOS(thermal_conductivity=3.0, heat_capacity=1000.0, thermal_expansion=2.0e-5)
    assert material.thermal_conductivity == pytest.approx(3.0)
    assert material.heat_capacity == pytest.approx(1000.0)
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
    """The single alpha affects the density only when thermal_density (a layer's use_thermal_eos) is set."""
    material = ConstantDensityEOS(reference_density=4000.0, thermal_expansion=3.0e-5)
    athermal = material.calc_material_state(0.0, 2000.0, thermal_density=False)
    thermal = material.calc_material_state(0.0, 2000.0, thermal_density=True)
    assert athermal["density"] == pytest.approx(4000.0)
    assert thermal["density"] == pytest.approx(4000.0 * math.exp(-3.0e-5 * (2000.0 - 300.0)))


def test_unknown_material_keyword_is_rejected():
    with pytest.raises(TypeError, match="shear_modulus"):
        ConstantDensityEOS(shear_modulus=1.0e10)


def test_bare_material_returns_its_constants():
    """With no models attached the state is the static constants, no melt, and the law's density."""
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
    material = _melting_material()

    solid = material.calc_material_state(_PRESSURE, 1500.0)
    assert solid["melt_fraction"] == 0.0
    assert solid["shear_modulus"] == pytest.approx(_SHEAR)

    molten = material.calc_material_state(_PRESSURE, 1900.0)
    reference = make_partial_melt("henning", _MELT_CONFIG).calc_partial_melt(1900.0, 1.0e20, _SHEAR)
    assert molten["melt_fraction"] == pytest.approx(0.75)
    assert molten["shear_viscosity"] == pytest.approx(reference[1], rel=1e-13)
    assert molten["shear_modulus"] == pytest.approx(reference[2], rel=1e-13)
    assert molten["shear_modulus"] < _SHEAR
    # Past breakdown the viscosity is the model's liquid viscosity.
    assert molten["shear_viscosity"] == pytest.approx(0.2)

    partial = material.calc_material_state(_PRESSURE, 1700.0)
    assert partial["shear_viscosity"] == pytest.approx(1.0e20 * math.exp(-13.5 * 0.25), rel=1e-12)


def test_partial_melt_leaves_the_bulk_modulus_unless_switched_on():
    material = _melting_material()
    bulk_premelt = material.calc_material_state(_PRESSURE, 1500.0)["bulk_modulus"]
    assert material.calc_material_state(_PRESSURE, 1900.0)["bulk_modulus"] == bulk_premelt

    # A constant melt bulk modulus: at 20 GPa a K' of 5 would make the melt as stiff as the solid.
    weakening_config = dict(_MELT_CONFIG, bulk_melt_weakening=True, liquid_bulk_modulus_pa=2.0e10,
                            liquid_bulk_modulus_derivative=0.0)
    material.set_partial_melt(make_partial_melt("henning", weakening_config))
    state = material.calc_material_state(_PRESSURE, 1700.0)
    reference = make_partial_melt("henning", weakening_config)
    expected = reference.calc_bulk_modulus_melt(1700.0, _PRESSURE, bulk_premelt, state["shear_modulus"])
    assert state["bulk_modulus"] == pytest.approx(expected, rel=1e-13)
    assert 2.0e10 < state["bulk_modulus"] < bulk_premelt


def test_partial_melt_leaves_the_density_and_bulk_viscosity_unless_switched_on():
    material = _melting_material()
    material.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e22}))
    cold = material.calc_material_state(_PRESSURE, 1500.0)
    hot = material.calc_material_state(_PRESSURE, 1800.0)
    assert hot["density"] == cold["density"] == pytest.approx(3300.0)
    assert hot["bulk_viscosity"] == cold["bulk_viscosity"] == pytest.approx(1.0e22)

    switched = dict(_MELT_CONFIG, density_melt_mixing=True, liquid_density_kg_m3=2750.0,
                    bulk_viscosity_melt_weakening=True)
    material.set_partial_melt(make_partial_melt("henning", switched))
    # At 1 GPa; by 20 GPa this constant-density solid would be lighter than the compressed melt.
    state = material.calc_material_state(1.0e9, 1800.0)
    reference = make_partial_melt("henning", switched)
    assert state["melt_fraction"] == pytest.approx(0.5)
    assert state["density"] == pytest.approx(reference.calc_mixture_density(1800.0, 1.0e9, 3300.0), rel=1e-14)
    assert state["density"] < 3300.0
    assert state["bulk_viscosity"] == pytest.approx(
        reference.calc_bulk_viscosity_melt(1800.0, 1.0e22, state["shear_viscosity"]), rel=1e-14)
    assert state["bulk_viscosity"] < 1.0e22


def test_an_attached_model_wrapper_is_an_empty_shell():
    """Attaching a partial-melt model moves it into the material; the wrapper then raises instead of reading freed
    memory. A viscosity model is copied in, so its wrapper stays usable."""
    material = _make_material()
    melt = make_partial_melt("henning")
    viscosity = make_viscosity("constant", {"reference_viscosity_pas": 1.0e20})
    material.set_partial_melt(melt)
    material.set_shear_viscosity(viscosity)
    with pytest.raises(RuntimeError):
        melt.calc_bulk_modulus_melt(1700.0, 0.0, 1.0e11, 1.0e10)
    assert viscosity.calc_viscosity(1000.0, 0.0) == pytest.approx(1.0e20)


def test_non_finite_temperature_skips_the_melt_model():
    material = _make_material()
    material.set_partial_melt(make_partial_melt("henning", _MELT_CONFIG))
    state = material.calc_material_state(_PRESSURE, math.nan)
    assert math.isnan(state["melt_fraction"])
    assert state["shear_modulus"] == pytest.approx(_SHEAR)


def test_pressure_law_bulk_modulus_takes_precedence_over_the_constant():
    eos_config = {
        "reference_density_kg_m3": 3500.0,
        "reference_bulk_modulus_pa": 1.3e11,
        "bulk_modulus_derivative": 4.5,
        "thermal_expansion_1_k": 3.0e-5,
        "bulk_modulus_static_pa": _BULK}
    material = make_material_eos("bm", eos_config)

    # thermal_density=False is what a layer without use_thermal_eos asks for.
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
    """Tabulated values win over the static constants; untabulated ones fall back to the constants."""
    material = make_material_eos("interpolate", {
        "radius_m": [0.0, 1.0e6],
        "density_kg_m3": [5000.0, 4000.0],
        "shear_modulus_pa": [8.0e10, 4.0e10],
        "shear_viscosity_pas": [1.0e22, 1.0e20],
        "shear_modulus_static_pa": 1.0e9,
        "bulk_modulus_static_pa": _BULK})
    state = material.calc_material_state(_PRESSURE, _TEMPERATURE, radius=5.0e5)
    assert state["density"] == pytest.approx(4500.0)
    assert state["shear_modulus"] == pytest.approx(6.0e10)
    assert state["shear_viscosity"] == pytest.approx(5.05e21)
    assert state["bulk_modulus"] == pytest.approx(_BULK)
    assert material.get_tabulated_shear_modulus(5.0e5) == pytest.approx(6.0e10)
    assert math.isnan(material.get_tabulated_bulk_modulus(5.0e5))


def test_factory_builds_the_nested_models():
    material = make_material_eos("constant", {
        "reference_density_kg_m3": 3300.0,
        "shear_modulus_static_pa": _SHEAR,
        "shear_viscosity_static_pas": 1.0e21,
        "bulk_viscosity_static_pas": 2.0e21,   # set so the dict comparison below holds no NaN
        "shear_viscosity": dict(model="reference", **_VISCOSITY_CONFIG),
        "partial_melt": dict(model="henning", **_MELT_CONFIG)})
    assert (material.shear_viscosity_set, material.bulk_viscosity_set, material.partial_melt_set) == (
        True, False, True)
    config = material.get_config_dict()
    assert config["shear_modulus_static_pa"] == pytest.approx(_SHEAR)
    assert config["shear_viscosity"]["model"] == "reference"
    assert config["partial_melt"]["solidus_k"] == pytest.approx(1600.0)
    assert "bulk_viscosity" not in config
    rebuilt = make_material_eos(config.pop("model"), config)
    assert rebuilt.calc_material_state(_PRESSURE, 1900.0) == material.calc_material_state(_PRESSURE, 1900.0)


def test_nested_model_table_needs_a_model_key():
    with pytest.raises(ValueError, match="model"):
        make_material_eos("constant", {"shear_viscosity": {"reference_viscosity_pas": 1.0e20}})


# =====================================================================================================================
# The layer hands viscosity and partial-melt models to its material
# =====================================================================================================================
def test_layer_helpers_put_the_models_on_the_material():
    layer = _layer(BaseLayer, "mantle")
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


@pytest.mark.parametrize("layer_class", [BaseLayer, SolidLiquidLayer, GasLayer])
def test_material_survives_config_and_binary_roundtrips(layer_class, tmp_path):
    layer = _layer(layer_class, "shell", temperature=1700.0, use_thermal_eos=True)
    layer.set_eos(_make_material(
        shear_viscosity_static=5.0e20,
        bulk_viscosity_static=6.0e20,   # set so the dict comparison below holds no NaN
        shear_modulus_pressure_derivative=1.4,
        shear_modulus_temperature_derivative=-8.0e6,
        shear_modulus_reference_temperature=1600.0,
        thermal_conductivity=3.5,
        heat_capacity=900.0))
    layer.set_shear_viscosity(make_viscosity("reference", _VISCOSITY_CONFIG))
    layer.set_partial_melt(make_partial_melt("henning", _MELT_CONFIG))

    cfg = layer.get_config_dict()
    assert cfg["temperature_k"] == pytest.approx(1700.0)
    assert cfg["use_thermal_eos"] is True
    material = cfg["material"]
    assert material["shear_modulus_static_pa"] == pytest.approx(_SHEAR)
    assert material["shear_modulus_pressure_derivative"] == pytest.approx(1.4)
    assert material["shear_modulus_temperature_derivative_pa_k"] == pytest.approx(-8.0e6)
    assert material["shear_modulus_reference_temperature_k"] == pytest.approx(1600.0)
    assert material["thermal_conductivity_w_mk"] == pytest.approx(3.5)
    assert material["heat_capacity_j_kgk"] == pytest.approx(900.0)
    for moved in ("shear_modulus_static_pa", "shear_viscosity", "partial_melt", "eos"):
        assert moved not in cfg

    path = tmp_path / "layer.tpyb"
    layer.save_binary(str(path))
    loaded = _layer(layer_class, "placeholder", radius_outer=1.0, mass=1.0)
    loaded.load_binary(str(path))
    assert loaded.temperature == pytest.approx(1700.0)
    assert loaded.use_thermal_eos is True
    assert loaded.get_config_dict()["material"] == material
    assert (loaded.shear_viscosity_set, loaded.partial_melt_set) == (True, True)


def test_material_table_builds_through_the_world_builder():
    """A layer's material table builds, merges with the layer type's defaults, and round trips."""
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
                    # Overrides one key of the default table; the rest of the table stands.
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
    # These two come from the mantle_rock defaults.
    assert material["shear_viscosity"]["model"] == "reference"
    assert material["partial_melt"]["model"] == "henning"
    rebuilt = construct_world(world.get_config_dict())
    assert rebuilt.mantle.get_config_dict()["material"] == material


@pytest.mark.parametrize("key, where", [
    ("shear_modulus_static_pa", "material"), ("eos", "material"), ("shear_viscosity", "material.shear_viscosity")])
def test_keys_left_on_the_layer_say_where_they_moved(key, where):
    layer = {"class": "base", "radius_fraction": 1.0}
    layer[key] = 6.0e10 if key.endswith("_pa") else {"model": "constant"}
    config = {"schema_version": "0.2.0", "name": "moved", "type": "terrestrial", "radius_m": 1.0e6,
              "mass_kg": 1.0e22, "layers": {"mantle": layer}}
    with pytest.raises(ValueError, match=r"layers\.mantle\." + where.replace(".", r"\.")):
        construct_world(config)
