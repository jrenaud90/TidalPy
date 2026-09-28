"""SolidLiquidLayer: construction, inherited physics, thermal calcs, sub-models, config dict, and binary round trip."""
import math

import pytest

from TidalPy.Cooling.cooling import ConvectiveCooling
from TidalPy.Material.eos.material_eos import ConstantDensityEOS
from TidalPy.Radiogenics.radiogenics import FixedRadiogenics, IsotopeRadiogenics
from TidalPy.Rheology.rheology import Maxwell
from TidalPy.Structures.layers.base import BaseLayer
from TidalPy.Structures.layers.physics import PhysicsLayer
from TidalPy.Structures.layers.solidliquid import SolidLiquidLayer
from TidalPy.Utilities.classes.classes import StructureBase, TidalPyBaseClass

_R_INNER = 3.485e6           # [m]
_R_OUTER = 6.371e6           # [m]
_MASS = 4.043e24             # [kg]
_SHEAR = 1.67e11             # [Pa]
_SHEAR_VISCOSITY = 1.0e21    # [Pa s]
_BULK_VISCOSITY = 2.0e21     # [Pa s]
_CONDUCTIVITY = 4.5          # [W m-1 K-1]
_EXPANSION = 2.0e-5          # [K-1]
_HEAT_CAPACITY = 1200.0      # [J kg-1 K-1]
_DENSITY = 4000.0            # [kg m-3]

# Static and thermal constants all belong to the material, so they go to the layer's EOS model.
_MATERIAL_KEYS = ("shear_modulus_static", "bulk_modulus_static", "shear_viscosity_static",
                  "bulk_viscosity_static", "thermal_conductivity", "thermal_expansion", "heat_capacity",
                  "reference_density")

_ALL_KEYS = (
    "name", "layer_index", "radius_inner_m", "radius_outer_m", "mass_kg",
    "material_name", "is_tidal", "tidal_scale",
    "love_number_k_re", "love_number_k_im",
    "love_number_h_re", "love_number_h_im",
    "love_number_l_re", "love_number_l_im",
)
_MATERIAL_TABLE_KEYS = ("thermal_conductivity_w_mk", "thermal_expansion_1_k", "heat_capacity_j_kgk")


def _make_layer(**kwargs):
    defaults = dict(
        name="mantle",
        layer_index=1,
        radius_inner=_R_INNER,
        radius_outer=_R_OUTER,
        mass=_MASS,
        material_name="perovskite",
        is_tidal=True,
        tidal_scale=1.0,
        shear_modulus_static=_SHEAR,
        bulk_modulus_static=3.57e11,
        shear_viscosity_static=_SHEAR_VISCOSITY,
        bulk_viscosity_static=_BULK_VISCOSITY,
        thermal_conductivity=_CONDUCTIVITY,
        thermal_expansion=_EXPANSION,
        heat_capacity=_HEAT_CAPACITY,
        reference_density=_DENSITY,
    )
    defaults.update(kwargs)
    material = {key: defaults.pop(key) for key in _MATERIAL_KEYS if key in defaults}
    layer = SolidLiquidLayer(**defaults)
    layer.set_eos(ConstantDensityEOS(**material))
    return layer


def _bare_layer():
    return SolidLiquidLayer("placeholder", 0, 0.0, 1.0, 1.0)


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
    """SolidLiquidLayer stores its values; the thermal constants are read from its material."""
    layer = _make_layer()
    assert layer.name == "mantle"
    assert layer.layer_index == 1
    assert layer.shear_modulus_static == pytest.approx(_SHEAR)
    assert layer.shear_viscosity_static == pytest.approx(_SHEAR_VISCOSITY)
    assert layer.bulk_viscosity_static == pytest.approx(_BULK_VISCOSITY)
    assert layer.thermal_conductivity == pytest.approx(_CONDUCTIVITY)
    assert layer.thermal_expansion == pytest.approx(_EXPANSION)
    assert layer.heat_capacity == pytest.approx(_HEAT_CAPACITY)


def test_defaults():
    """Without a material the constants are NaN; a default material gives k = 4, alpha = 0, c_p = 1200."""
    layer = SolidLiquidLayer("test", 0, 0.0, 1e6, 1e20)
    assert math.isnan(layer.thermal_conductivity)
    assert math.isnan(layer.shear_viscosity_static)
    layer.set_eos(ConstantDensityEOS())
    assert layer.thermal_conductivity == pytest.approx(4.0)
    assert layer.thermal_expansion == 0.0  # exactly zero keeps the density law athermal
    assert layer.heat_capacity == pytest.approx(1200.0)
    assert math.isnan(layer.shear_viscosity_static)
    assert math.isnan(layer.bulk_viscosity_static)
    assert (layer.is_solid, layer.is_static, layer.is_incompressible) == (True, True, False)


def test_assumption_flags_from_constructor():
    """Solver flags given to the constructor are reported by the properties and get_config_dict."""
    layer = _make_layer(is_solid=False, is_static=False, is_incompressible=True)
    assert (layer.is_solid, layer.is_static, layer.is_incompressible) == (False, False, True)
    config = layer.get_config_dict()
    assert (config["is_solid"], config["is_static"], config["is_incompressible"]) == (False, False, True)


def test_inherits_geometry():
    """Geometry from BaseLayer resolves on a SolidLiquidLayer."""
    layer = _make_layer()
    assert layer.thickness == pytest.approx(_R_OUTER - _R_INNER)
    assert layer.volume == pytest.approx((4.0 / 3.0) * math.pi * (_R_OUTER ** 3 - _R_INNER ** 3), rel=1e-9)
    assert layer.radius == pytest.approx(_R_OUTER)
    assert layer.radius_inner == pytest.approx(_R_INNER)
    assert layer.mass == pytest.approx(_MASS)


def test_inherits_eos():
    """The EOS profile from BaseLayer works on a SolidLiquidLayer."""
    layer = _make_layer()
    assert layer.eos_data_populated is False
    assert math.isnan(layer.get_density(_R_INNER))
    layer.update_eos_data([_R_INNER, _R_OUTER], [5000.0, 4000.0], [10.0, 9.81], [1e11, 0.0])
    assert layer.eos_data_populated is True
    assert layer.get_density(_R_INNER) == pytest.approx(5000.0)


def test_inherits_love_numbers():
    """Love numbers from PhysicsLayer are stored and bundled."""
    love_k, love_h, love_l = 0.3 - 0.01j, 0.6 - 0.02j, 0.1 - 0.005j
    layer = SolidLiquidLayer(
        "test",
        0,
        0.0,
        1e6,
        1e20,
        love_number_k=love_k,
        love_number_h=love_h,
        love_number_l=love_l)
    assert layer.love_number_k == pytest.approx(love_k)
    assert layer.love_number_h == pytest.approx(love_h)
    assert layer.love_number_l == pytest.approx(love_l)
    love_numbers = layer.love_numbers
    assert love_numbers.k == pytest.approx(love_k)
    assert love_numbers.h == pytest.approx(love_h)
    assert love_numbers.l == pytest.approx(love_l)


@pytest.mark.parametrize("temperature", [1000.0, 4000.0])
def test_thermal_conductivity_constant(temperature):
    assert _make_layer().calc_thermal_conductivity(temperature) == pytest.approx(_CONDUCTIVITY)


def test_thermal_diffusivity():
    """Diffusivity is k / (rho c_p) with the layer's bulk density."""
    layer = _make_layer()
    expected = _CONDUCTIVITY / (layer.density_bulk * _HEAT_CAPACITY)
    assert layer.calc_thermal_diffusivity(2000.0) == pytest.approx(expected)


def test_adiabatic_gradient_without_eos_profile_is_zero():
    assert _make_layer().calc_adiabatic_temperature_gradient(2000.0) == pytest.approx(0.0)


@pytest.mark.parametrize("temperature_base, temperature_top", [(3000.0, 1000.0), (2000.0, 2000.0)])
def test_heat_flux_conductive(temperature_base, temperature_top):
    """Conductive flux is k dT / thickness (zero when isothermal)."""
    expected = _CONDUCTIVITY * (temperature_base - temperature_top) / (_R_OUTER - _R_INNER)
    assert _make_layer().calc_heat_flux_conductive(temperature_base, temperature_top) == pytest.approx(expected)


def test_heat_flux_conductive_zero_thickness():
    layer = SolidLiquidLayer("thin", 0, 1e6, 1e6, 1e20)
    assert layer.calc_heat_flux_conductive(3000.0, 1000.0) == pytest.approx(0.0)


def test_radiogenic_heating_without_model_is_zero():
    assert _make_layer().calc_radiogenic_heating(0.0, _MASS) == pytest.approx(0.0)


def test_attach_submodels_sets_flags():
    """set_cooling and set_radiogenics flip their flags, which start False."""
    layer = _make_layer()
    assert layer.cooling_set is False
    assert layer.radiogenics_set is False
    layer.set_cooling(ConvectiveCooling())
    layer.set_radiogenics(FixedRadiogenics(fixed_heat_production=1.0e-11))
    assert layer.cooling_set is True
    assert layer.radiogenics_set is True


def test_attach_radiogenics_changes_heating():
    layer = _make_layer()
    layer.set_radiogenics(FixedRadiogenics(fixed_heat_production=2.0e-11))
    assert layer.calc_radiogenic_heating(0.0, _MASS) == pytest.approx(2.0e-11 * _MASS, rel=1e-12)


def test_get_config_dict_has_all_keys():
    config = _make_layer().get_config_dict()
    for key in _ALL_KEYS:
        assert key in config, f"Missing key: {key}"
    for key in _MATERIAL_TABLE_KEYS:
        assert key in config["material"], f"Missing material key: {key}"


def test_get_config_dict_values():
    """Layer values sit at the top level and material constants in the material table."""
    config = _make_layer().get_config_dict()
    material = config["material"]
    assert config["name"] == "mantle"
    assert config["layer_index"] == 1
    assert material["shear_modulus_static_pa"] == pytest.approx(_SHEAR)
    assert material["shear_viscosity_static_pas"] == pytest.approx(_SHEAR_VISCOSITY)
    assert material["bulk_viscosity_static_pas"] == pytest.approx(_BULK_VISCOSITY)
    assert material["thermal_conductivity_w_mk"] == pytest.approx(_CONDUCTIVITY)
    assert material["thermal_expansion_1_k"] == pytest.approx(_EXPANSION)
    assert material["heat_capacity_j_kgk"] == pytest.approx(_HEAT_CAPACITY)
    assert material["reference_density_kg_m3"] == pytest.approx(_DENSITY)
    for moved in ("thermal_conductivity_ref_w_mk", "reference_density_kg_m3", "reference_temperature_k"):
        assert moved not in config
    assert config["love_number_k_re"] == pytest.approx(0.0)
    assert config["love_number_k_im"] == pytest.approx(0.0)


def test_get_config_dict_class_and_thermal_model_tables():
    """Cooling and radiogenics models appear as sub-tables once attached."""
    layer = _make_layer()
    config = layer.get_config_dict()
    assert config["class"] == "solidliquid"
    assert "cooling" not in config
    assert "radiogenics" not in config
    layer.set_cooling(ConvectiveCooling(convection_alpha=0.9))
    layer.set_radiogenics(FixedRadiogenics(fixed_heat_production=3.0e-11))
    config = layer.get_config_dict()
    assert config["cooling"]["convection_alpha"] == pytest.approx(0.9)
    assert config["radiogenics"]["fixed_heat_production_w_kg"] == pytest.approx(3.0e-11)


def test_save_config(tmp_path):
    """save_config writes the layer and material keys to TOML."""
    path = str(tmp_path / "layer.toml")
    _make_layer().save_config(path)
    data = _read_toml(path)
    assert data["name"] == "mantle"
    assert data["material"]["thermal_conductivity_w_mk"] == pytest.approx(_CONDUCTIVITY)
    assert data["material"]["heat_capacity_j_kgk"] == pytest.approx(_HEAT_CAPACITY)


def test_binary_roundtrip(tmp_path):
    """save_binary then load_binary restores every stored field."""
    loaded = _roundtrip(_make_layer(), tmp_path, _bare_layer())
    assert loaded.name == "mantle"
    assert loaded.layer_index == 1
    assert loaded.radius_outer == pytest.approx(_R_OUTER)
    assert loaded.radius_inner == pytest.approx(_R_INNER)
    assert loaded.mass == pytest.approx(_MASS)
    assert loaded.shear_modulus_static == pytest.approx(_SHEAR)
    assert loaded.shear_viscosity_static == pytest.approx(_SHEAR_VISCOSITY)
    assert loaded.bulk_viscosity_static == pytest.approx(_BULK_VISCOSITY)
    assert loaded.thermal_conductivity == pytest.approx(_CONDUCTIVITY)
    assert loaded.thermal_expansion == pytest.approx(_EXPANSION)
    assert loaded.heat_capacity == pytest.approx(_HEAT_CAPACITY)


def test_binary_roundtrip_derived_fields(tmp_path):
    """Derived geometry is recomputed after load_binary."""
    loaded = _roundtrip(_make_layer(), tmp_path, _bare_layer())
    assert loaded.thickness == pytest.approx(_R_OUTER - _R_INNER)


def test_binary_load_file_not_found():
    with pytest.raises(FileNotFoundError):
        _make_layer().load_binary("/nonexistent/path/xyz.tpyb")


def test_binary_roundtrip_with_all_submodels(tmp_path):
    """Rheology, cooling, and radiogenics models survive a binary round trip."""
    frequency = 1.0e-5
    layer = _make_layer()
    layer.set_shear_rheology(Maxwell())
    layer.set_cooling(ConvectiveCooling(convection_alpha=0.9, critical_rayleigh=1200.0))
    layer.set_radiogenics(FixedRadiogenics(fixed_heat_production=3.0e-11))
    modulus_before = layer.calc_complex_shear_modulus(frequency)
    heating_before = layer.calc_radiogenic_heating(0.0, _MASS)

    loaded = _roundtrip(layer, tmp_path, _make_layer(name="placeholder"))
    assert loaded.shear_rheology_set is True
    assert loaded.cooling_set is True
    assert loaded.radiogenics_set is True
    modulus_after = loaded.calc_complex_shear_modulus(frequency)
    assert modulus_after.real == pytest.approx(modulus_before.real, rel=1e-12)
    assert modulus_after.imag == pytest.approx(modulus_before.imag, rel=1e-12)
    assert loaded.calc_radiogenic_heating(0.0, _MASS) == pytest.approx(heating_before, rel=1e-12)


def test_binary_roundtrip_isotope_radiogenics(tmp_path):
    """Variable-length isotope radiogenics survives a binary round trip."""
    layer = _make_layer()
    layer.set_radiogenics(IsotopeRadiogenics.from_dataset("modern_day_chondritic"))
    heating_before = layer.calc_radiogenic_heating(0.0, _MASS)
    loaded = _roundtrip(layer, tmp_path, _make_layer(name="placeholder"))
    assert loaded.radiogenics_set is True
    assert loaded.calc_radiogenic_heating(0.0, _MASS) == pytest.approx(heating_before, rel=1e-12)


def test_binary_roundtrip_without_submodels(tmp_path):
    """A layer without sub-models loads back with every model flag False."""
    loaded = _roundtrip(_make_layer(), tmp_path, _make_layer(name="placeholder"))
    assert loaded.shear_rheology_set is False
    assert loaded.bulk_rheology_set is False
    assert loaded.cooling_set is False
    assert loaded.radiogenics_set is False


@pytest.mark.parametrize(
    "parent",
    [PhysicsLayer, BaseLayer, StructureBase, TidalPyBaseClass],
    ids=lambda cls: cls.__name__)
def test_is_instance_of_parents(parent):
    assert isinstance(_make_layer(), parent)
