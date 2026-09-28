"""PhysicsLayer: construction, inherited geometry, complex moduli, rheology, config dict, and binary round trip."""
import math

import pytest

from TidalPy.Material.eos.material_eos import ConstantDensityEOS
from TidalPy.Rheology import rheology
from TidalPy.Structures.layers.base import BaseLayer
from TidalPy.Structures.layers.physics import PhysicsLayer
from TidalPy.Utilities.classes.classes import StructureBase, TidalPyBaseClass

_R_INNER = 3.485e6    # [m]
_R_OUTER = 6.371e6    # [m]
_MASS = 4.043e24      # [kg]
_SHEAR = 1.67e11      # [Pa]
_BULK = 3.57e11       # [Pa]
_VISCOSITY = 1.0e21   # [Pa s]
_FREQUENCY = 1.0e-5   # [rad s-1]

_CONFIG_KEYS = (
    "name", "layer_index", "radius_inner_m", "radius_outer_m",
    "mass_kg", "material_name", "is_tidal", "tidal_scale",
    "is_solid", "is_static", "is_incompressible",
    "love_number_k_re", "love_number_k_im",
    "love_number_h_re", "love_number_h_im",
    "love_number_l_re", "love_number_l_im",
)


def _material(shear=_SHEAR, bulk=_BULK, shear_viscosity=_VISCOSITY, bulk_viscosity=_VISCOSITY):
    return ConstantDensityEOS(
        reference_density=4500.0,
        shear_modulus_static=shear,
        bulk_modulus_static=bulk,
        shear_viscosity_static=shear_viscosity,
        bulk_viscosity_static=bulk_viscosity)


def _make_mantle():
    layer = PhysicsLayer(
        name="mantle",
        layer_index=1,
        radius_inner=_R_INNER,
        radius_outer=_R_OUTER,
        mass=_MASS,
        material_name="perovskite",
        is_tidal=True,
        tidal_scale=1.0,
        love_number_k=0.0 + 0.0j,
        love_number_h=0.0 + 0.0j,
        love_number_l=0.0 + 0.0j,
    )
    layer.set_eos(_material())
    return layer


def _roundtrip(layer, tmp_path):
    path = str(tmp_path / "layer.tpyb")
    layer.save_binary(path)
    loaded = PhysicsLayer("placeholder", 99, 0.0, 1.0, 1.0)
    loaded.load_binary(path)
    return loaded


def _read_toml(path):
    try:
        import tomllib
    except ImportError:
        tomllib = pytest.importorskip("tomli")
    with open(path, "rb") as toml_file:
        return tomllib.load(toml_file)


def test_construction():
    """PhysicsLayer stores its values; the static constants are read from its material."""
    layer = PhysicsLayer("core", 0, 0.0, 3.485e6, 1.932e24)
    layer.set_eos(_material(5e10, 2e11, 1e20, 1e22))
    assert layer.name == "core"
    assert layer.layer_index == 0
    assert layer.radius_inner == pytest.approx(0.0)
    assert layer.radius_outer == pytest.approx(3.485e6)
    assert layer.mass == pytest.approx(1.932e24)
    assert layer.shear_modulus_static == pytest.approx(5e10)
    assert layer.bulk_modulus_static == pytest.approx(2e11)
    assert layer.shear_viscosity_static == pytest.approx(1e20)
    assert layer.bulk_viscosity_static == pytest.approx(1e22)
    assert layer.love_number_k == pytest.approx(0.0 + 0.0j)
    assert layer.love_number_h == pytest.approx(0.0 + 0.0j)
    assert layer.love_number_l == pytest.approx(0.0 + 0.0j)


def test_defaults():
    """Without a material the constants are NaN; a default material's moduli are 0 and viscosities NaN."""
    layer = PhysicsLayer("test", 0, 0.0, 1e6, 1e20)
    assert math.isnan(layer.shear_modulus_static)
    layer.set_eos(ConstantDensityEOS())
    assert layer.shear_modulus_static == pytest.approx(0.0)
    assert layer.bulk_modulus_static == pytest.approx(0.0)
    assert math.isnan(layer.shear_viscosity_static)
    assert math.isnan(layer.bulk_viscosity_static)
    assert layer.love_number_k == pytest.approx(0.0 + 0.0j)


def test_unset_static_viscosity_fails_loudly():
    """Without a static viscosity a viscous rheology's modulus is NaN; an elastic one is unaffected."""
    layer = PhysicsLayer("test", 0, 0.0, 1e6, 1e20)
    layer.set_eos(ConstantDensityEOS(shear_modulus_static=_SHEAR))
    layer.set_shear_rheology(rheology.Elastic())
    assert layer.calc_complex_shear_modulus(_FREQUENCY) == pytest.approx(_SHEAR + 0.0j)
    layer.set_shear_rheology(rheology.Maxwell())
    modulus = layer.calc_complex_shear_modulus(_FREQUENCY)
    assert math.isnan(modulus.real) and math.isnan(modulus.imag)
    assert layer.love_number_h == pytest.approx(0.0 + 0.0j)
    assert layer.love_number_l == pytest.approx(0.0 + 0.0j)


def test_inherits_base_fields():
    """Stored fields from BaseLayer resolve on a PhysicsLayer."""
    layer = _make_mantle()
    assert layer.name == "mantle"
    assert layer.layer_index == 1
    assert layer.material_name == "perovskite"
    assert layer.is_tidal is True
    assert layer.tidal_scale == pytest.approx(1.0)


def test_inherits_geometry():
    """Derived geometry from BaseLayer resolves on a PhysicsLayer."""
    layer = _make_mantle()
    assert layer.radius == pytest.approx(_R_OUTER)
    assert layer.mass == pytest.approx(_MASS)
    assert layer.thickness == pytest.approx(_R_OUTER - _R_INNER)
    assert layer.volume == pytest.approx((4.0 / 3.0) * math.pi * (_R_OUTER ** 3 - _R_INNER ** 3), rel=1e-9)
    assert layer.surface_area_outer == pytest.approx(4.0 * math.pi * _R_OUTER ** 2, rel=1e-9)
    assert layer.surface_area_inner == pytest.approx(4.0 * math.pi * _R_INNER ** 2, rel=1e-9)


def test_inherits_eos():
    """The EOS profile from BaseLayer works on a PhysicsLayer."""
    layer = _make_mantle()
    assert layer.eos_data_populated is False
    assert math.isnan(layer.get_density(_R_INNER))
    layer.update_eos_data([_R_INNER, _R_OUTER], [5000.0, 4000.0], [10.0, 9.81], [1e11, 0.0])
    assert layer.eos_data_populated is True
    assert layer.get_density(_R_INNER) == pytest.approx(5000.0)


@pytest.mark.parametrize("frequency", [1e-10, 1e-5, 1.0, 1e5])
@pytest.mark.parametrize("method, expected", [
    ("calc_complex_shear_modulus", _SHEAR),
    ("calc_complex_bulk_modulus", _BULK),
])
def test_complex_modulus_without_rheology_is_static(method, expected, frequency):
    """Without a rheology the complex modulus is the real static modulus at every frequency."""
    modulus = getattr(_make_mantle(), method)(frequency)
    assert modulus.real == pytest.approx(expected)
    assert modulus.imag == pytest.approx(0.0)


def test_complex_modulus_zero_modulus():
    layer = PhysicsLayer("test", 0, 0.0, 1e6, 1e20)
    layer.set_eos(ConstantDensityEOS(shear_modulus_static=0.0, bulk_modulus_static=0.0))
    assert layer.calc_complex_shear_modulus(_FREQUENCY).real == pytest.approx(0.0)
    assert layer.calc_complex_bulk_modulus(_FREQUENCY).real == pytest.approx(0.0)


def test_rheology_not_set_initially():
    layer = _make_mantle()
    assert layer.shear_rheology_set is False
    assert layer.bulk_rheology_set is False


def test_get_config_dict_has_all_keys():
    config = _make_mantle().get_config_dict()
    for key in _CONFIG_KEYS:
        assert key in config, f"Missing key: {key}"


def test_get_config_dict_values():
    """Layer values sit at the top level and the static constants in the material table."""
    config = _make_mantle().get_config_dict()
    material = config["material"]
    assert config["name"] == "mantle"
    assert config["layer_index"] == 1
    assert config["radius_inner_m"] == pytest.approx(_R_INNER)
    assert config["radius_outer_m"] == pytest.approx(_R_OUTER)
    assert config["mass_kg"] == pytest.approx(_MASS)
    assert material["shear_modulus_static_pa"] == pytest.approx(_SHEAR)
    assert material["bulk_modulus_static_pa"] == pytest.approx(_BULK)
    assert material["shear_viscosity_static_pas"] == pytest.approx(_VISCOSITY)
    assert material["bulk_viscosity_static_pas"] == pytest.approx(_VISCOSITY)
    assert config["love_number_k_re"] == pytest.approx(0.0)
    assert config["love_number_k_im"] == pytest.approx(0.0)
    assert config["love_number_h_re"] == pytest.approx(0.0)
    assert config["love_number_l_re"] == pytest.approx(0.0)
    # A default physics layer is a compressible static solid.
    assert config["is_solid"] is True
    assert config["is_static"] is True
    assert config["is_incompressible"] is False


def test_assumption_flags_from_constructor():
    """Solver flags given to the constructor are reported by the properties and get_config_dict."""
    layer = PhysicsLayer("ocean", 0, 0.0, 1e6, 1e20, is_solid=False, is_static=False, is_incompressible=True)
    assert (layer.is_solid, layer.is_static, layer.is_incompressible) == (False, False, True)
    config = layer.get_config_dict()
    assert (config["is_solid"], config["is_static"], config["is_incompressible"]) == (False, False, True)


def test_get_config_dict_class_and_model_tables():
    """Attached models appear as sub-tables keyed by model name; viscosity and melt sit in the material table."""
    from TidalPy.PartialMelt import make_partial_melt
    from TidalPy.Viscosity import make_viscosity
    layer = _make_mantle()
    config = layer.get_config_dict()
    assert config["class"] == "physics"
    for key in ("shear_rheology", "bulk_rheology"):
        assert key not in config
    for key in ("shear_viscosity", "bulk_viscosity", "partial_melt"):
        assert key not in config["material"]
    elastic_name = rheology.Elastic().model_name
    layer.set_shear_rheology(rheology.Andrade(0.4, 1.5))
    layer.set_bulk_rheology(rheology.Elastic())
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e20}))
    layer.set_partial_melt(make_partial_melt("henning"))
    config = layer.get_config_dict()
    assert config["shear_rheology"]["model"] == "andrade"
    assert config["shear_rheology"]["alpha"] == pytest.approx(0.4)
    assert config["bulk_rheology"] == {"model": elastic_name}
    assert config["material"]["shear_viscosity"]["reference_viscosity_pas"] == pytest.approx(1.0e20)
    assert "bulk_viscosity" not in config["material"]
    assert config["material"]["partial_melt"]["model"] == "henning"
    assert "shear_viscosity" not in config


def test_save_config(tmp_path):
    """save_config writes the layer and material keys to TOML."""
    path = str(tmp_path / "layer.toml")
    _make_mantle().save_config(path)
    data = _read_toml(path)
    material = data["material"]
    assert data["name"] == "mantle"
    assert material["shear_modulus_static_pa"] == pytest.approx(_SHEAR)
    assert material["bulk_modulus_static_pa"] == pytest.approx(_BULK)
    assert material["shear_viscosity_static_pas"] == pytest.approx(_VISCOSITY)
    assert material["bulk_viscosity_static_pas"] == pytest.approx(_VISCOSITY)
    assert data["love_number_k_re"] == pytest.approx(0.0)
    assert data["love_number_k_im"] == pytest.approx(0.0)
    assert data["love_number_h_re"] == pytest.approx(0.0)
    assert data["love_number_l_re"] == pytest.approx(0.0)


def test_binary_roundtrip(tmp_path):
    """save_binary then load_binary restores the geometry and physics fields."""
    loaded = _roundtrip(_make_mantle(), tmp_path)
    assert loaded.name == "mantle"
    assert loaded.layer_index == 1
    assert loaded.radius_inner == pytest.approx(_R_INNER)
    assert loaded.radius_outer == pytest.approx(_R_OUTER)
    assert loaded.mass == pytest.approx(_MASS)
    assert loaded.material_name == "perovskite"
    assert loaded.shear_modulus_static == pytest.approx(_SHEAR)
    assert loaded.bulk_modulus_static == pytest.approx(_BULK)
    assert loaded.shear_viscosity_static == pytest.approx(_VISCOSITY)
    assert loaded.bulk_viscosity_static == pytest.approx(_VISCOSITY)
    assert loaded.love_number_k == pytest.approx(0.0 + 0.0j)
    assert loaded.love_number_h == pytest.approx(0.0 + 0.0j)
    assert loaded.love_number_l == pytest.approx(0.0 + 0.0j)


def test_binary_roundtrip_derived_fields(tmp_path):
    """Derived geometry is recomputed after load_binary."""
    loaded = _roundtrip(_make_mantle(), tmp_path)
    expected_volume = (4.0 / 3.0) * math.pi * (_R_OUTER ** 3 - _R_INNER ** 3)
    assert loaded.thickness == pytest.approx(_R_OUTER - _R_INNER, rel=1e-9)
    assert loaded.volume == pytest.approx(expected_volume, rel=1e-9)


def test_binary_roundtrip_complex_modulus(tmp_path):
    """After load_binary the complex modulus uses the restored static modulus."""
    modulus = _roundtrip(_make_mantle(), tmp_path).calc_complex_shear_modulus(_FREQUENCY)
    assert modulus.real == pytest.approx(_SHEAR)
    assert modulus.imag == pytest.approx(0.0)


def test_binary_load_file_not_found():
    with pytest.raises(FileNotFoundError):
        _make_mantle().load_binary("/nonexistent/path/xyz.tpyb")


def test_attach_shear_rheology_sets_flag():
    layer = _make_mantle()
    assert layer.shear_rheology_set is False
    layer.set_shear_rheology(rheology.Maxwell())
    assert layer.shear_rheology_set is True
    assert layer.bulk_rheology_set is False


def test_attach_shear_rheology_changes_complex_modulus():
    """With Maxwell attached the complex modulus matches the model function."""
    layer = _make_mantle()
    layer.set_shear_rheology(rheology.Maxwell())
    modulus = layer.calc_complex_shear_modulus(_FREQUENCY)
    expected = rheology.maxwell(_SHEAR, _VISCOSITY, _FREQUENCY)
    assert modulus.real == pytest.approx(expected.real, rel=1e-9)
    assert modulus.imag == pytest.approx(expected.imag, rel=1e-9)
    assert modulus.imag != 0.0


def test_attach_rheology_consumes_wrapper():
    """A rheology model moves into the layer, so it cannot be attached twice."""
    layer = _make_mantle()
    model = rheology.Maxwell()
    layer.set_shear_rheology(model)
    with pytest.raises(ValueError):
        layer.set_bulk_rheology(model)


@pytest.mark.parametrize("shear_model, bulk_model", [
    pytest.param(None, None, id="none"),
    pytest.param(rheology.Maxwell, None, id="shear_only"),
    pytest.param(rheology.Maxwell, lambda: rheology.Andrade(alpha=0.25, zeta=2.0), id="shear_and_bulk"),
])
def test_binary_roundtrip_with_rheology(tmp_path, shear_model, bulk_model):
    """Attached rheology models, and their absence, survive a binary round trip."""
    layer = _make_mantle()
    if shear_model is not None:
        layer.set_shear_rheology(shear_model())
    if bulk_model is not None:
        layer.set_bulk_rheology(bulk_model())
    shear_before = layer.calc_complex_shear_modulus(_FREQUENCY)
    bulk_before = layer.calc_complex_bulk_modulus(_FREQUENCY)

    loaded = _roundtrip(layer, tmp_path)
    assert loaded.shear_rheology_set is (shear_model is not None)
    assert loaded.bulk_rheology_set is (bulk_model is not None)
    shear_after = loaded.calc_complex_shear_modulus(_FREQUENCY)
    bulk_after = loaded.calc_complex_bulk_modulus(_FREQUENCY)
    assert shear_after.real == pytest.approx(shear_before.real, rel=1e-12)
    assert shear_after.imag == pytest.approx(shear_before.imag, rel=1e-12)
    assert bulk_after.real == pytest.approx(bulk_before.real, rel=1e-12)
    assert bulk_after.imag == pytest.approx(bulk_before.imag, rel=1e-12)


@pytest.mark.parametrize("parent", [BaseLayer, StructureBase, TidalPyBaseClass], ids=lambda cls: cls.__name__)
def test_is_instance_of_parents(parent):
    assert isinstance(_make_mantle(), parent)
