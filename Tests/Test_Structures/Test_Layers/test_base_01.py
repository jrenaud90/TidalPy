"""BaseLayer: construction, geometry, EOS profile, config dict, TOML save, and binary round trip."""
import math

import numpy as np
import pytest

from TidalPy.Structures.layers.base import BaseLayer
from TidalPy.Utilities.classes.classes import StructureBase, TidalPyBaseClass

_R_INNER = 3.485e6  # [m] about the core-mantle boundary
_R_OUTER = 6.371e6  # [m]
_MASS = 4.043e24    # [kg]
_VOLUME = (4.0 / 3.0) * math.pi * (_R_OUTER ** 3 - _R_INNER ** 3)

_MANTLE_KWARGS = dict(
    name="mantle",
    layer_index=1,
    radius_inner=_R_INNER,
    radius_outer=_R_OUTER,
    mass=_MASS,
    material_name="perovskite",
    is_tidal=True,
    tidal_scale=1.0,
)
_MANTLE_EXPECTED = dict(_MANTLE_KWARGS)
_MANTLE_CONFIG = dict(
    name="mantle",
    layer_index=1,
    radius_inner_m=_R_INNER,
    radius_outer_m=_R_OUTER,
    mass_kg=_MASS,
    material_name="perovskite",
    is_tidal=True,
    tidal_scale=1.0,
    is_solid=True,
    is_static=True,
    is_incompressible=False,
    temperature_k=0.0,
    use_thermal_eos=False,
    use_heating=False,
)
_PROFILE = dict(gravity=[10.0, 9.8], pressure=[1e11, 0.0])


def _make_mantle():
    return BaseLayer(**_MANTLE_KWARGS)


def _assert_matches(actual, expected, label):
    if expected is None or isinstance(expected, bool):
        assert actual is expected, label
    elif isinstance(expected, (float, complex)):
        assert actual == pytest.approx(expected), label
    else:
        assert actual == expected, label


def _roundtrip(layer, tmp_path):
    path = str(tmp_path / "layer.tpyb")
    layer.save_binary(path)
    loaded = BaseLayer("placeholder", 99, 0.0, 1.0, 1.0)
    loaded.load_binary(path)
    return loaded


def _read_toml(path):
    try:
        import tomllib
    except ImportError:
        tomllib = pytest.importorskip("tomli")
    with open(path, "rb") as toml_file:
        return tomllib.load(toml_file)


@pytest.mark.parametrize("kwargs, expected", [
    pytest.param(
        dict(name="core", layer_index=0, radius_inner=0.0, radius_outer=3.485e6, mass=1.932e24),
        dict(name="core", layer_index=0, radius_inner=0.0, radius_outer=3.485e6, mass=1.932e24, material_name="",
             is_tidal=True, tidal_scale=None),  # an unset tidal scale means the layer's volume fraction
        id="defaults"),
    pytest.param(_MANTLE_KWARGS, _MANTLE_EXPECTED, id="all_kwargs"),
    pytest.param(
        dict(name="gas", layer_index=2, radius_inner=6e6, radius_outer=7e7, mass=1e26, is_tidal=False,
             tidal_scale=0.5),
        dict(is_tidal=False, tidal_scale=0.5),
        id="not_tidal"),
])
def test_construction(kwargs, expected):
    """BaseLayer stores its constructor values and defaults."""
    layer = BaseLayer(**kwargs)
    for key, value in expected.items():
        _assert_matches(getattr(layer, key), value, key)


@pytest.mark.parametrize("attribute, expected, rel", [
    ("density_bulk", _MASS / _VOLUME, 1e-12),
    ("radius", _R_OUTER, None),
    ("thickness", _R_OUTER - _R_INNER, None),
    ("volume", _VOLUME, 1e-9),
    ("surface_area_outer", 4.0 * math.pi * _R_OUTER ** 2, 1e-9),
    ("surface_area_inner", 4.0 * math.pi * _R_INNER ** 2, 1e-9),
])
def test_derived_geometry(attribute, expected, rel):
    """Derived geometry properties follow from the radii and mass."""
    assert getattr(_make_mantle(), attribute) == pytest.approx(expected, rel=rel)


def test_zero_inner_radius():
    """A full sphere has no inner surface and a thickness equal to its radius."""
    layer = BaseLayer("core", 0, 0.0, 3.485e6, 1.932e24)
    assert layer.surface_area_inner == pytest.approx(0.0)
    assert layer.thickness == pytest.approx(layer.radius_outer)


def test_schema_version():
    assert _make_mantle().get_schema_version_str() == "0.2.0"


@pytest.mark.parametrize("getter", ["get_density", "get_gravity", "get_pressure"])
def test_eos_getters_are_nan_before_update(getter):
    """A fresh layer has no EOS profile and its profile getters return NaN."""
    layer = _make_mantle()
    assert layer.eos_data_populated is False
    assert math.isnan(getattr(layer, getter)(_R_INNER))


@pytest.mark.parametrize("density, radius, expected, rel", [
    pytest.param([5000.0, 4000.0], _R_INNER, 5000.0, None, id="inner"),
    pytest.param([5000.0, 4000.0], _R_OUTER, 4000.0, None, id="outer"),
    pytest.param([5000.0, 4000.0], _R_INNER - 1.0, 5000.0, None, id="below_clamped"),
    pytest.param([5000.0, 4000.0], _R_OUTER + 1.0, 4000.0, None, id="above_clamped"),
    pytest.param([5000.0, 3000.0], 0.5 * (_R_INNER + _R_OUTER), 4000.0, 1e-12, id="midpoint"),
])
def test_eos_linear_interpolation(density, radius, expected, rel):
    """update_eos_data populates the profile, which interpolates linearly and clamps outside its range."""
    layer = _make_mantle()
    layer.update_eos_data([_R_INNER, _R_OUTER], density, _PROFILE["gravity"], _PROFILE["pressure"])
    assert layer.eos_data_populated is True
    assert layer.get_density(radius) == pytest.approx(expected, rel=rel)


@pytest.mark.parametrize("n_points", [3, 10, 100])
def test_eos_multi_point_profiles(n_points):
    """A profile with many samples reproduces its end values."""
    layer = _make_mantle()
    radius = list(np.linspace(_R_INNER, _R_OUTER, n_points))
    fraction = [(value - _R_INNER) / (_R_OUTER - _R_INNER) for value in radius]
    density = [5000.0 - 1000.0 * value for value in fraction]
    pressure = [1e11 * (1.0 - value) for value in fraction]
    layer.update_eos_data(radius, density, [10.0] * n_points, pressure)
    assert layer.eos_data_populated is True
    assert layer.get_density(radius[0]) == pytest.approx(density[0], rel=1e-9)
    assert layer.get_density(radius[-1]) == pytest.approx(density[-1], rel=1e-9)


@pytest.mark.parametrize("method, expected", [
    ("calc_surface_area", 4.0 * math.pi),
    ("calc_volume_sphere", 4.0 * math.pi / 3.0),
])
def test_inherited_structure_calcs(method, expected):
    """StructureBase geometry helpers work on a BaseLayer (evaluated at unit radius)."""
    assert getattr(_make_mantle(), method)(1.0) == pytest.approx(expected)


def test_get_config_dict():
    """get_config_dict holds every constructor value."""
    config = _make_mantle().get_config_dict()
    for key, value in _MANTLE_CONFIG.items():
        assert key in config, f"Missing key: {key}"
        _assert_matches(config[key], value, key)


def test_get_config_dict_class_and_eos_table():
    """The dict names its builder class and carries a material table only once an EOS is attached."""
    from TidalPy.Material.eos.material_eos import ConstantDensityEOS
    layer = _make_mantle()
    assert layer.get_config_dict()["class"] == "base"
    assert "material" not in layer.get_config_dict()
    eos = ConstantDensityEOS(reference_density=4400.0)
    expected_model = eos.model_name
    layer.set_eos(eos)
    eos_config = layer.get_config_dict()["material"]
    assert eos_config["model"] == expected_model
    assert eos_config["reference_density_kg_m3"] == pytest.approx(4400.0)


def test_save_config_writes_toml(tmp_path):
    """save_config writes the geometry keys to TOML."""
    path = str(tmp_path / "layer.toml")
    _make_mantle().save_config(path)
    data = _read_toml(path)
    for key in ("name", "radius_inner_m", "radius_outer_m", "mass_kg"):
        _assert_matches(data[key], _MANTLE_CONFIG[key], key)


def test_binary_roundtrip(tmp_path):
    """save_binary then load_binary restores every stored field."""
    loaded = _roundtrip(_make_mantle(), tmp_path)
    for key, value in _MANTLE_EXPECTED.items():
        _assert_matches(getattr(loaded, key), value, key)


def test_binary_roundtrip_derived_fields(tmp_path):
    """Derived geometry is recomputed after load_binary."""
    loaded = _roundtrip(_make_mantle(), tmp_path)
    assert loaded.thickness == pytest.approx(_R_OUTER - _R_INNER, rel=1e-9)
    assert loaded.volume == pytest.approx(_VOLUME, rel=1e-9)


def test_binary_roundtrip_eos_profile_not_preserved(tmp_path):
    """The EOS profile is derived data, so it is not serialized."""
    layer = _make_mantle()
    layer.update_eos_data([_R_INNER, _R_OUTER], [5000.0, 4000.0], _PROFILE["gravity"], _PROFILE["pressure"])
    assert layer.eos_data_populated is True
    assert _roundtrip(layer, tmp_path).eos_data_populated is False


def test_binary_roundtrip_preserves_eos_model(tmp_path):
    """An attached material EOS model is saved with the layer and rebuilt on load."""
    from TidalPy.Material.eos import BirchMurnaghanEOS
    layer = _make_mantle()
    layer.set_eos(BirchMurnaghanEOS(3300.0, 1.2e11, 4.2))
    expected = layer.get_config_dict()["material"]
    path = str(tmp_path / "layer.tpyb")
    layer.save_binary(path)
    loaded = BaseLayer("placeholder", 0, 0.0, 1.0, 1.0)
    assert loaded.eos_set is False
    loaded.load_binary(path)
    assert loaded.eos_set is True
    assert loaded.get_config_dict()["material"] == expected


def test_binary_load_file_not_found():
    with pytest.raises(FileNotFoundError):
        _make_mantle().load_binary("/nonexistent/path/xyz.tpyb")


@pytest.mark.parametrize("parent", [StructureBase, TidalPyBaseClass], ids=lambda cls: cls.__name__)
def test_is_instance_of_parents(parent):
    assert isinstance(_make_mantle(), parent)
