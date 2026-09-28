"""GasLayer: construction, inherited geometry and EOS, config dict, TOML save, and binary round trip."""
import math

import pytest

from TidalPy.Material.eos.material_eos import ConstantDensityEOS
from TidalPy.Rheology.rheology import Maxwell
from TidalPy.Structures.layers.base import BaseLayer
from TidalPy.Structures.layers.gas import GasLayer
from TidalPy.Structures.layers.physics import PhysicsLayer
from TidalPy.Utilities.classes.classes import StructureBase, TidalPyBaseClass

_R_INNER = 0.0
_R_OUTER = 7.0e7                 # [m]
_MASS = 1.0e27                   # [kg]
_MOLECULAR_WEIGHT = 2.0e-3       # [kg/mol] hydrogen
_ADIABATIC_INDEX = 1.4
_REFERENCE_TEMPERATURE = 300.0   # [K]
_REFERENCE_DENSITY = 1.2         # [kg/m^3]

_MATERIAL_KEYS = ("shear_modulus_static", "bulk_modulus_static", "shear_viscosity_static", "bulk_viscosity_static")

_ALL_KEYS = (
    "name", "layer_index", "radius_inner_m", "radius_outer_m", "mass_kg",
    "material_name", "is_tidal", "tidal_scale",
    "love_number_k_re", "love_number_k_im",
    "love_number_h_re", "love_number_h_im",
    "love_number_l_re", "love_number_l_im",
    "mean_molecular_weight_kg_mol", "adiabatic_index",
    "reference_temperature_k",
)


def _make_layer(**kwargs):
    defaults = dict(
        name="atmosphere",
        layer_index=0,
        radius_inner=_R_INNER,
        radius_outer=_R_OUTER,
        mass=_MASS,
        material_name="hydrogen",
        is_tidal=False,
        tidal_scale=1.0,
        mean_molecular_weight=_MOLECULAR_WEIGHT,
        adiabatic_index=_ADIABATIC_INDEX,
        reference_temperature=_REFERENCE_TEMPERATURE,
        reference_density=_REFERENCE_DENSITY,
    )
    defaults.update(kwargs)
    # The static constants belong to the material, so they go to the layer's EOS model.
    material = {key: defaults.pop(key) for key in _MATERIAL_KEYS if key in defaults}
    layer = GasLayer(**defaults)
    layer.set_eos(ConstantDensityEOS(**material))
    return layer


def _roundtrip(layer, tmp_path, placeholder):
    path = str(tmp_path / "layer.tpyb")
    layer.save_binary(path)
    placeholder.load_binary(path)
    return placeholder


def _read_toml(path):
    try:
        import tomllib
    except ImportError:
        tomllib = pytest.importorskip("tomli")
    with open(path, "rb") as toml_file:
        return tomllib.load(toml_file)


def test_construction():
    """GasLayer stores its constructor values."""
    layer = _make_layer()
    assert layer.name == "atmosphere"
    assert layer.layer_index == 0
    assert layer.mean_molecular_weight == pytest.approx(_MOLECULAR_WEIGHT)
    assert layer.adiabatic_index == pytest.approx(_ADIABATIC_INDEX)
    assert layer.reference_temperature == pytest.approx(_REFERENCE_TEMPERATURE)
    assert layer.reference_density == pytest.approx(_REFERENCE_DENSITY)


def test_defaults():
    """Optional gas parameters and solver flags take their documented defaults."""
    layer = GasLayer("test", 0, 0.0, 1e7, 1e27)
    assert layer.mean_molecular_weight == pytest.approx(2.0e-3)
    assert layer.adiabatic_index == pytest.approx(1.4)
    assert layer.reference_temperature == pytest.approx(300.0)
    assert layer.reference_density == pytest.approx(1.0)
    # A gas carries no shear stress, so the radial solver treats it as a static liquid by default.
    assert (layer.is_solid, layer.is_static, layer.is_incompressible) == (False, True, False)
    assert layer.get_config_dict()["is_solid"] is False


def test_assumption_flags_from_constructor():
    layer = _make_layer(is_solid=True, is_static=False, is_incompressible=True)
    assert (layer.is_solid, layer.is_static, layer.is_incompressible) == (True, False, True)


def test_inherits_geometry():
    """Geometry from BaseLayer resolves on a GasLayer."""
    layer = _make_layer(radius_inner=1e7, radius_outer=7e7)
    assert layer.thickness == pytest.approx(7e7 - 1e7)
    assert layer.volume == pytest.approx((4.0 / 3.0) * math.pi * (7e7 ** 3 - 1e7 ** 3), rel=1e-9)
    assert layer.radius == pytest.approx(7e7)


def test_inherits_eos():
    """The EOS profile from BaseLayer works on a GasLayer."""
    layer = _make_layer()
    assert layer.eos_data_populated is False
    assert math.isnan(layer.get_density(_R_OUTER))
    layer.update_eos_data([0.0, _R_OUTER], [5.0, 1.0], [0.0, 25.0], [1e5, 0.0])
    assert layer.eos_data_populated is True
    assert layer.get_density(0.0) == pytest.approx(5.0)


def test_get_config_dict_has_all_keys():
    config = _make_layer().get_config_dict()
    for key in _ALL_KEYS:
        assert key in config, f"Missing key: {key}"


def test_get_config_dict_values():
    """get_config_dict values match the constructor arguments."""
    config = _make_layer().get_config_dict()
    assert config["name"] == "atmosphere"
    assert config["mean_molecular_weight_kg_mol"] == pytest.approx(_MOLECULAR_WEIGHT)
    assert config["adiabatic_index"] == pytest.approx(_ADIABATIC_INDEX)
    assert config["reference_temperature_k"] == pytest.approx(_REFERENCE_TEMPERATURE)
    # The density is the material's, not a layer key.
    assert "reference_density_kg_m3" not in config


def test_get_config_dict_class_name():
    assert _make_layer().get_config_dict()["class"] == "gas"


def test_save_config(tmp_path):
    """save_config writes the gas keys to TOML."""
    path = str(tmp_path / "layer.toml")
    _make_layer().save_config(path)
    data = _read_toml(path)
    assert data["name"] == "atmosphere"
    assert data["mean_molecular_weight_kg_mol"] == pytest.approx(_MOLECULAR_WEIGHT)
    assert data["adiabatic_index"] == pytest.approx(_ADIABATIC_INDEX)


def test_binary_roundtrip(tmp_path):
    """save_binary then load_binary restores every stored field."""
    loaded = _roundtrip(_make_layer(), tmp_path, GasLayer("placeholder", 0, 0.0, 1.0, 1.0))
    assert loaded.name == "atmosphere"
    assert loaded.layer_index == 0
    assert loaded.radius_outer == pytest.approx(_R_OUTER)
    assert loaded.mass == pytest.approx(_MASS)
    assert loaded.mean_molecular_weight == pytest.approx(_MOLECULAR_WEIGHT)
    assert loaded.adiabatic_index == pytest.approx(_ADIABATIC_INDEX)
    assert loaded.reference_temperature == pytest.approx(_REFERENCE_TEMPERATURE)
    assert loaded.reference_density == pytest.approx(_REFERENCE_DENSITY)


def test_binary_roundtrip_derived_fields(tmp_path):
    """Derived geometry is recomputed after load_binary."""
    layer = _make_layer(radius_inner=1e7, radius_outer=7e7)
    loaded = _roundtrip(layer, tmp_path, GasLayer("placeholder", 0, 0.0, 1.0, 1.0))
    assert loaded.thickness == pytest.approx(7e7 - 1e7)


def test_binary_load_file_not_found():
    with pytest.raises(FileNotFoundError):
        _make_layer().load_binary("/nonexistent/path/xyz.tpyb")


def test_attach_rheology_sets_flag():
    """GasLayer inherits set_shear_rheology from PhysicsLayer."""
    layer = _make_layer(shear_modulus_static=1.0e9, shear_viscosity_static=1.0e18)
    assert layer.shear_rheology_set is False
    layer.set_shear_rheology(Maxwell())
    assert layer.shear_rheology_set is True


def test_binary_roundtrip_with_rheology(tmp_path):
    """An attached rheology survives a binary round trip."""
    frequency = 1.0e-5
    layer = _make_layer(shear_modulus_static=1.0e9, shear_viscosity_static=1.0e18)
    layer.set_shear_rheology(Maxwell())
    modulus_before = layer.calc_complex_shear_modulus(frequency)
    loaded = _roundtrip(layer, tmp_path, _make_layer(name="placeholder"))
    assert loaded.shear_rheology_set is True
    modulus_after = loaded.calc_complex_shear_modulus(frequency)
    assert modulus_after.real == pytest.approx(modulus_before.real, rel=1e-12)
    assert modulus_after.imag == pytest.approx(modulus_before.imag, rel=1e-12)


@pytest.mark.parametrize(
    "parent",
    [PhysicsLayer, BaseLayer, StructureBase, TidalPyBaseClass],
    ids=lambda cls: cls.__name__)
def test_is_instance_of_parents(parent):
    assert isinstance(_make_layer(), parent)
