"""Tests for the layer material state: ``PhysicsLayer.calc_material_state`` and the keys that feed it.

``calc_material_state`` is the one place that maps a point's pressure and temperature onto a layer's attached
models, so each stage is checked against the model it delegates to: the EOS density and bulk modulus, the linear
static shear-modulus law, the viscosity models, the partial-melt model, and the rheologies.

Requires the Cython extensions to be compiled first::

    uv pip install -v <repo_root>
"""
import math

import pytest

from TidalPy.Material_x.eos import make_material_eos
from TidalPy.partial_melt_x import make_partial_melt
from TidalPy.rheology_x.rheology import Maxwell, maxwell
from TidalPy.structures_x.configs.world_builder import construct_world
from TidalPy.structures_x.layers.gas import GasLayer
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.structures_x.layers.solidliquid import SolidLiquidLayer
from TidalPy.viscosity_x import make_viscosity

_SHEAR = 6.0e10          # [Pa]
_BULK = 1.2e11           # [Pa]
_PRESSURE = 2.0e10       # [Pa]
_TEMPERATURE = 1500.0    # [K]
_FREQUENCY = 1.0e-5      # [rad/s]

_VISCOSITY_CONFIG = {
    "reference_viscosity_pas": 1.0e19,
    "reference_temperature_k": 1400.0,
    "molar_activation_energy_j_mol": 3.0e5,
    "molar_activation_volume_m3_mol": 1.0e-6}


def _make_layer(**kwargs):
    return PhysicsLayer(
        "mantle", 0, 0.0, 1.0e6, 1.0e22,
        shear_modulus_static=_SHEAR,
        bulk_modulus_static=_BULK,
        **kwargs)


def test_default_material_parameters():
    layer = _make_layer()
    assert layer.temperature == 0.0
    assert layer.use_thermal_eos is False
    assert layer.shear_modulus_pressure_derivative == 0.0
    assert layer.shear_modulus_temperature_derivative == 0.0
    assert layer.shear_modulus_reference_temperature == pytest.approx(300.0)


def test_bare_layer_returns_its_constants():
    """No models attached: the static constants, no melt, no density, and real complex moduli."""
    state = _make_layer().calc_material_state(_PRESSURE, _TEMPERATURE, _FREQUENCY)
    assert math.isnan(state["density"])
    assert state["melt_fraction"] == 0.0
    assert state["shear_modulus"] == state["premelt_shear_modulus"] == pytest.approx(_SHEAR)
    assert state["bulk_modulus"] == state["premelt_bulk_modulus"] == pytest.approx(_BULK)
    assert math.isnan(state["shear_viscosity"])
    assert state["complex_shear_modulus"] == pytest.approx(complex(_SHEAR, 0.0))
    assert state["complex_bulk_modulus"] == pytest.approx(complex(_BULK, 0.0))


def test_static_viscosity_is_the_fallback_without_a_viscosity_model():
    layer = _make_layer(shear_viscosity_static=1.0e21, bulk_viscosity_static=2.0e21)
    state = layer.calc_material_state(_PRESSURE, _TEMPERATURE)
    assert state["shear_viscosity"] == pytest.approx(1.0e21)
    assert state["bulk_viscosity"] == pytest.approx(2.0e21)


@pytest.mark.parametrize("pressure,temperature", [(0.0, 300.0), (2.0e10, 1500.0), (1.0e11, 3000.0)])
def test_linear_shear_law(pressure, temperature):
    layer = _make_layer(
        shear_modulus_pressure_derivative=1.5,
        shear_modulus_temperature_derivative=-1.0e7,
        shear_modulus_reference_temperature=300.0)
    expected = _SHEAR + 1.5 * pressure - 1.0e7 * (temperature - 300.0)
    assert layer.calc_material_state(pressure, temperature)["shear_modulus"] == pytest.approx(expected, rel=1e-13)


def test_shear_law_is_floored_at_the_minimum_modulus():
    layer = _make_layer(shear_modulus_temperature_derivative=-1.0e9)
    shear = layer.calc_material_state(0.0, 5000.0)["shear_modulus"]
    assert 0.0 < shear < 1.0e-3 * _SHEAR


def test_layer_temperature_is_the_default_temperature():
    layer = _make_layer(shear_modulus_temperature_derivative=-1.0e7, temperature=_TEMPERATURE)
    at_layer_temperature = layer.calc_material_state(0.0)["shear_modulus"]
    assert at_layer_temperature == layer.calc_material_state(0.0, _TEMPERATURE)["shear_modulus"]
    layer.temperature = 800.0
    assert layer.calc_material_state(0.0)["shear_modulus"] == layer.calc_material_state(0.0, 800.0)["shear_modulus"]
    assert layer.calc_material_state(0.0)["shear_modulus"] > at_layer_temperature


def test_viscosity_model_is_evaluated_at_the_temperature_and_pressure():
    layer = _make_layer()
    layer.set_shear_viscosity(make_viscosity("reference", _VISCOSITY_CONFIG))
    expected = make_viscosity("reference", _VISCOSITY_CONFIG).calc_viscosity(_TEMPERATURE, _PRESSURE)
    state = layer.calc_material_state(_PRESSURE, _TEMPERATURE)
    assert state["shear_viscosity"] == pytest.approx(expected, rel=1e-13)
    # The cold limit is rigid, so a layer left at 0 K still solves as an elastic body.
    assert math.isinf(layer.calc_material_state(_PRESSURE, 0.0)["shear_viscosity"])


def test_rheology_gives_the_complex_modulus_at_the_frequency():
    layer = _make_layer()
    layer.set_shear_viscosity(make_viscosity("reference", _VISCOSITY_CONFIG))
    layer.set_shear_rheology(Maxwell())
    state = layer.calc_material_state(_PRESSURE, _TEMPERATURE, _FREQUENCY)
    expected = maxwell(state["shear_modulus"], state["shear_viscosity"], _FREQUENCY)
    assert state["complex_shear_modulus"] == pytest.approx(expected, rel=1e-13)
    assert state["complex_shear_modulus"].imag > 0.0
    # Without a frequency the rheology is skipped.
    static = layer.calc_material_state(_PRESSURE, _TEMPERATURE)["complex_shear_modulus"]
    assert static == pytest.approx(complex(_SHEAR, 0.0))


def test_partial_melt_weakens_the_shear_pair_and_reports_the_melt_fraction():
    melt_config = {"solidus_k": 1600.0, "liquidus_k": 2000.0}
    layer = _make_layer()
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e20}))
    layer.set_partial_melt(make_partial_melt("henning", melt_config))

    solid = layer.calc_material_state(_PRESSURE, 1500.0)
    assert solid["melt_fraction"] == 0.0
    assert solid["shear_modulus"] == pytest.approx(solid["premelt_shear_modulus"])

    molten = layer.calc_material_state(_PRESSURE, 1900.0)
    reference = make_partial_melt("henning", melt_config).calc_partial_melt(1900.0, 1.0e20, _SHEAR, 1.0e20)
    assert molten["melt_fraction"] == pytest.approx(0.75)
    assert molten["premelt_shear_modulus"] == pytest.approx(_SHEAR)
    assert molten["shear_viscosity"] == pytest.approx(reference[1], rel=1e-13)
    assert molten["shear_modulus"] == pytest.approx(reference[2], rel=1e-13)
    assert molten["shear_modulus"] < _SHEAR


def test_eos_bulk_modulus_takes_precedence_over_the_layer_constant():
    eos_config = {
        "reference_density_kg_m3": 3500.0,
        "reference_bulk_modulus_pa": 1.3e11,
        "bulk_modulus_derivative": 4.5,
        "thermal_expansion_1_k": 3.0e-5}
    reference = make_material_eos("bm", eos_config)
    layer = _make_layer()
    layer.set_eos(make_material_eos("bm", eos_config))

    # Athermal by default: the EOS does not see the temperature.
    state = layer.calc_material_state(_PRESSURE, _TEMPERATURE)
    assert state["density"] == pytest.approx(reference.calc_density(_PRESSURE), rel=1e-13)
    assert state["bulk_modulus"] == pytest.approx(reference.calc_bulk_modulus(_PRESSURE), rel=1e-13)
    assert state["bulk_modulus"] != pytest.approx(_BULK)

    layer.use_thermal_eos = True
    thermal = layer.calc_material_state(_PRESSURE, _TEMPERATURE)
    assert thermal["density"] == pytest.approx(reference.calc_density(_PRESSURE, _TEMPERATURE), rel=1e-13)
    assert thermal["bulk_modulus"] == pytest.approx(
        reference.calc_bulk_modulus(_PRESSURE, _TEMPERATURE), rel=1e-13)
    assert thermal["density"] < state["density"]


def test_constant_density_eos_leaves_the_layer_bulk_modulus():
    layer = _make_layer()
    layer.set_eos(make_material_eos("constant", {"reference_density_kg_m3": 3300.0}))
    state = layer.calc_material_state(_PRESSURE, _TEMPERATURE)
    assert state["density"] == pytest.approx(3300.0)
    assert state["bulk_modulus"] == pytest.approx(_BULK)


def test_interpolated_eos_tables_take_precedence():
    layer = _make_layer()
    layer.set_eos(make_material_eos("interpolate", {
        "radius_m": [0.0, 1.0e6],
        "density_kg_m3": [5000.0, 4000.0],
        "shear_modulus_pa": [8.0e10, 4.0e10],
        "shear_viscosity_pas": [1.0e22, 1.0e20]}))
    state = layer.calc_material_state(_PRESSURE, _TEMPERATURE, radius=5.0e5)
    assert state["density"] == pytest.approx(4500.0)
    assert state["shear_modulus"] == pytest.approx(6.0e10)
    assert state["shear_viscosity"] == pytest.approx(5.05e21)
    assert state["bulk_modulus"] == pytest.approx(_BULK)   # no bulk table: the layer constant


_MATERIAL_KWARGS = dict(
    temperature=1700.0,
    shear_modulus_pressure_derivative=1.4,
    shear_modulus_temperature_derivative=-8.0e6,
    shear_modulus_reference_temperature=1600.0,
    use_thermal_eos=True)


@pytest.mark.parametrize("layer_class", [PhysicsLayer, SolidLiquidLayer, GasLayer])
def test_material_parameters_survive_config_and_binary_roundtrips(layer_class, tmp_path):
    layer = layer_class("shell", 0, 0.0, 1.0e6, 1.0e22, shear_modulus_static=_SHEAR, **_MATERIAL_KWARGS)
    cfg = layer.get_config_dict()
    assert cfg["temperature_k"] == pytest.approx(1700.0)
    assert cfg["shear_modulus_pressure_derivative"] == pytest.approx(1.4)
    assert cfg["shear_modulus_temperature_derivative_pa_k"] == pytest.approx(-8.0e6)
    assert cfg["shear_modulus_reference_temperature_k"] == pytest.approx(1600.0)
    assert cfg["use_thermal_eos"] is True

    path = tmp_path / "layer.tpyb"
    layer.save_binary(str(path))
    loaded = layer_class("placeholder", 0, 0.0, 1.0, 1.0)
    loaded.load_binary(str(path))
    assert loaded.temperature == pytest.approx(1700.0)
    assert loaded.use_thermal_eos is True
    assert loaded.shear_modulus_pressure_derivative == pytest.approx(1.4)
    assert loaded.shear_modulus_temperature_derivative == pytest.approx(-8.0e6)
    assert loaded.shear_modulus_reference_temperature == pytest.approx(1600.0)


def test_material_keys_build_through_the_world_builder():
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
                "shear_modulus_pressure_derivative": 1.4,
                "shear_modulus_temperature_derivative_pa_k": -8.0e6,
                "shear_modulus_reference_temperature_k": 1600.0,
                "use_thermal_eos": True,
            },
        },
    }
    world = construct_world(config)
    assert world.mantle.temperature == pytest.approx(1650.0)
    assert world.mantle.use_thermal_eos is True
    assert world.mantle.shear_modulus_temperature_derivative == pytest.approx(-8.0e6)
    rebuilt = construct_world(world.get_config_dict())
    assert rebuilt.mantle.temperature == pytest.approx(1650.0)
    assert rebuilt.mantle.shear_modulus_reference_temperature == pytest.approx(1600.0)
