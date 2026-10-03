"""The Layer class: construction, material and rheology forms, switches, state, cooling and radiogenics, the hand-set
profile, configuration, and binary round trips."""

import math

import numpy as np
import pytest

from TidalPy.Cooling.cooling import make_cooling
from TidalPy.Material import Material, Phase, load_material
from TidalPy.Radiogenics.radiogenics import make_radiogenics
from TidalPy.Rheology import Andrade, Elastic, Maxwell, make_rheology
from TidalPy.Structures.layers import Layer
from TidalPy.Utilities.classes.classes import StructureBase, TidalPyBaseClass
from TidalPy.Viscosity import ConstantViscosity

_RADIUS_INNER = 1.0e6
_RADIUS_OUTER = 2.0e6


def _material(rheology=None):
    return Material(solid=Phase(
        eos={"model": "constant", "reference_density_kg_m3": 3300.0, "bulk_modulus_pa": 1.3e11,
             "thermal_expansion_1_k": 3.0e-5, "reference_temperature_k": 300.0},
        shear_modulus={"model": "constant", "shear_modulus_pa": 6.0e10},
        shear_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e19},
        shear_rheology=rheology))


def _layer(**kwargs):
    return Layer("mantle", 1, _RADIUS_INNER, _RADIUS_OUTER, 5.0e22, **kwargs)


# =====================================================================================================================
# Construction
# =====================================================================================================================
def test_geometry_and_defaults():
    layer = _layer()
    assert layer.name == "mantle" and layer.layer_index == 1
    assert layer.radius == layer.radius_outer == _RADIUS_OUTER
    assert layer.thickness == pytest.approx(1.0e6)
    volume = (4.0 / 3.0) * math.pi * (_RADIUS_OUTER ** 3 - _RADIUS_INNER ** 3)
    assert layer.volume == pytest.approx(volume)
    assert layer.density_bulk == pytest.approx(5.0e22 / volume)
    assert layer.surface_area_outer == pytest.approx(4.0 * math.pi * _RADIUS_OUTER ** 2)
    # Every switch defaults to the simplest case.
    assert layer.use_tides and layer.is_volume_fixed and layer.is_static
    assert not (layer.use_thermal_expansion or layer.use_melting or layer.use_pressure_melting
                or layer.use_melt_density or layer.use_heating or layer.is_incompressible)
    assert layer.state == "auto" and layer.tidal_scale is None and layer.temperature == 0.0
    assert layer.material is None and not layer.material_set
    assert layer.shear_rheology is None and not layer.cooling_set and not layer.radiogenics_set
    assert isinstance(layer, StructureBase) and isinstance(layer, TidalPyBaseClass)


@pytest.mark.parametrize("radii", ((-1.0, 1.0e6), (5.0e5, 1.0e5), (0.0, math.inf)))
def test_bad_radii_are_refused(radii):
    with pytest.raises(ValueError, match="radius_inner"):
        Layer("bad", 0, radii[0], radii[1], 0.0)


def test_bad_mass_and_state_are_refused():
    with pytest.raises(ValueError, match="negative mass"):
        Layer("bad", 0, 0.0, 1.0, -1.0)
    with pytest.raises(ValueError, match="'auto', 'solid', or 'liquid'"):
        Layer("bad", 0, 0.0, 1.0, 0.0, state="gas")


# =====================================================================================================================
# Material
# =====================================================================================================================
def test_material_from_a_material_a_name_or_a_table():
    material = _material()
    assert _layer(material=material).material.get_config_dict() == material.get_config_dict()
    assert _layer(material="simple_rock").material.get_config_dict() == \
        load_material("simple_rock").get_config_dict()
    preset = {"preset": "simple_rock", "solid": {"shear_modulus": {"shear_modulus_pa": 4.0e10}}}
    assert _layer(material=preset).material.solid.shear_modulus.shear_modulus == 4.0e10


@pytest.mark.parametrize("material", [3.0, ConstantViscosity()])
def test_material_of_the_wrong_kind_is_refused(material):
    with pytest.raises(TypeError, match="Material, a MatPack name, or a material table"):
        _layer(material=material)


def test_material_is_shared_not_copied():
    material = _material()
    first, second = _layer(material=material), _layer(material=material)
    assert first.material.get_config_dict() == second.material.get_config_dict()
    first.material = None
    assert first.material is None and second.material is not None


def test_calc_state_uses_the_layer_switches_and_temperature():
    layer = _layer(material=_material(), temperature=1300.0)
    assert layer.calc_state(0.0)["density"] == 3300.0
    layer.use_thermal_expansion = True
    assert layer.calc_state(0.0)["density"] == pytest.approx(3300.0 * math.exp(-3.0e-5 * 1000.0))
    assert layer.calc_state(0.0, 300.0)["density"] == pytest.approx(3300.0)
    with pytest.raises(ValueError, match="has no material"):
        _layer().calc_state(0.0)


@pytest.mark.parametrize("switch", ["use_thermal_expansion", "use_melting", "use_pressure_melting",
                                    "use_melt_density", "use_heating", "use_tides", "is_static",
                                    "is_incompressible", "is_volume_fixed"])
def test_switches_set_and_read_back(switch):
    layer = _layer()
    before = getattr(layer, switch)
    setattr(layer, switch, not before)
    assert getattr(layer, switch) is (not before)


# =====================================================================================================================
# State
# =====================================================================================================================
def test_state_follows_the_material_unless_set():
    water = _layer(material="water")
    assert water.state == "auto" and water.is_liquid
    rock = _layer(material="simple_rock")
    assert not rock.is_liquid
    rock.state = "liquid"
    assert rock.is_liquid
    water.state = "solid"
    assert not water.is_liquid


def test_only_a_melting_layer_can_change_state():
    layer = _layer(material="peridotite")
    assert not layer.can_change_state
    layer.use_melting = True
    assert layer.can_change_state
    layer.state = "solid"
    assert not layer.can_change_state
    assert not _layer(material="simple_rock", use_melting=True).can_change_state


# =====================================================================================================================
# Rheology
# =====================================================================================================================
def test_rheology_override_and_material_default():
    layer = _layer(material=_material(rheology={"model": "andrade", "alpha": 0.2}))
    assert isinstance(layer.shear_rheology, Andrade) and layer.shear_rheology.alpha == 0.2
    layer.shear_rheology = "maxwell"
    assert isinstance(layer.shear_rheology, Maxwell)
    assert layer.get_config_dict()["shear_rheology"] == {"model": "maxwell"}
    layer.shear_rheology = None
    assert isinstance(layer.shear_rheology, Andrade)
    assert "shear_rheology" not in layer.get_config_dict()


def test_rheology_of_the_wrong_kind_is_refused():
    with pytest.raises(ValueError, match="takes a rheology model"):
        _layer(shear_rheology=ConstantViscosity())
    with pytest.raises(ValueError, match="needs a 'model' key"):
        _layer(bulk_rheology={"alpha": 0.3})


def test_layer_constant_complex_modulus():
    """The one-argument form applies the rheology to the material at zero pressure and the layer's temperature."""
    frequency = 1.0e-5
    layer = _layer(material=_material(), temperature=1500.0)
    assert layer.calc_complex_shear_modulus(frequency) == complex(6.0e10, 0.0)
    layer.shear_rheology = Maxwell()
    assert layer.calc_complex_shear_modulus(frequency) == pytest.approx(
        Maxwell().calc_complex_modulus(6.0e10, 1.0e19, frequency))
    layer.bulk_rheology = Elastic()
    assert layer.calc_complex_bulk_modulus(frequency) == complex(1.3e11, 0.0)
    assert np.isnan(_layer().calc_complex_shear_modulus(frequency).real)


# =====================================================================================================================
# Cooling and radiogenics
# =====================================================================================================================
def test_cooling_is_shared_and_radiogenics_move_in_once():
    cooling = make_cooling("convection")
    radiogenics = make_radiogenics("fixed", {"fixed_heat_production_w_kg": 1.0e-11, "ref_time_s": 0.0})
    layer = _layer(cooling=cooling, radiogenics=radiogenics)
    assert layer.cooling_set and layer.radiogenics_set
    assert layer.calc_radiogenic_heating(0.0, 2.0e22) == pytest.approx(2.0e11)
    # The cooling model is shared, like a rheology: the same model serves another layer, and the wrapper stays usable.
    other = _layer()
    other.cooling = cooling
    assert other.cooling.get_config_dict() == cooling.get_config_dict() == layer.cooling.get_config_dict()
    # A name or a table builds the model; None clears it.
    other.cooling = {"model": "convection", "critical_rayleigh": 1600.0}
    assert other.cooling.critical_rayleigh == 1600.0
    other.cooling = "conduction"
    assert other.cooling.model_name == "conduction"
    other.cooling = None
    assert other.cooling is None and not other.cooling_set
    # A model of another family is refused, naming the slot, rather than clearing it.
    with pytest.raises(ValueError, match="a layer's cooling model"):
        other.cooling = radiogenics
    with pytest.raises(ValueError, match="a layer's cooling model"):
        other.cooling = make_rheology("maxwell")
    assert _layer().calc_radiogenic_heating(0.0, 1.0) == 0.0


# =====================================================================================================================
# Hand-set profile
# =====================================================================================================================
def test_profile_getters_are_nan_until_populated_then_interpolate():
    layer = _layer()
    assert math.isnan(layer.get_density(1.5e6)) and not layer.eos_data_populated
    layer.update_eos_data([1.0e6, 2.0e6], [4000.0, 3000.0], [1.0, 2.0], [2.0e9, 0.0])
    assert layer.eos_data_populated and not layer.viscoelastic_populated
    assert layer.get_density(1.5e6) == pytest.approx(3500.0)
    assert layer.get_pressure(1.25e6) == pytest.approx(1.5e9)
    # A hand-set profile carries no material state.
    assert math.isnan(layer.get_shear_modulus(1.5e6))


@pytest.mark.parametrize("descending, density_size", [(False, 3), (True, 10)], ids=["length", "descending"])
def test_a_hand_set_profile_must_agree(descending, density_size):
    radius = np.linspace(1.0e6, 2.0e6, 10)
    if descending:
        radius = radius[::-1]
    with pytest.raises(ValueError):
        _layer().update_eos_data(radius, np.full(density_size, 3000.0), np.zeros(10), np.zeros(10))


# =====================================================================================================================
# Configuration and binary
# =====================================================================================================================
def _full_layer():
    return _layer(
        material="peridotite", temperature=1600.0, use_tides=False, tidal_scale=0.25, state="solid",
        is_static=False, is_incompressible=True, is_volume_fixed=False, use_thermal_expansion=True,
        use_melting=True, use_pressure_melting=True, use_melt_density=True, use_heating=True,
        shear_rheology="maxwell", bulk_rheology={"model": "zener", "relaxed_modulus_frac": 0.4},
        cooling=make_cooling("conduction"), radiogenics=make_radiogenics("isotope"))


def test_config_dict():
    config = _full_layer().get_config_dict()
    assert config["state"] == "solid" and config["tidal_scale"] == 0.25 and config["temperature_k"] == 1600.0
    assert config["use_tides"] is False and config["use_melt_density"] is True
    assert config["material"] == load_material("peridotite").get_config_dict()
    assert config["bulk_rheology"] == {"model": "zener", "relaxed_modulus_frac": 0.4}
    assert config["cooling"]["model"] == "conduction" and config["radiogenics"]["model"] == "isotope"
    assert "tidal_scale" not in _layer().get_config_dict()


def test_binary_round_trip(tmp_path):
    layer = _full_layer()
    path = str(tmp_path / "layer.tpyb")
    layer.save_binary(path)
    loaded = _layer()
    loaded.load_binary(path)
    assert loaded.get_config_dict() == layer.get_config_dict()
    loaded_state, state = loaded.calc_state(1.0e9), layer.calc_state(1.0e9)
    for key, value in state.items():
        assert loaded_state[key] == value or (value != value and loaded_state[key] != loaded_state[key]), key
