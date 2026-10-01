"""Cooling models (off, convection, conduction): physics against legacy values, factory, vectorization, and I/O.

Legacy values in ``frozen/test_cooling_01.npz`` come from the classic backend (0.8.0 snapshot 8b8e0b12).
"""

from pathlib import Path

import numpy as np
import pytest

import TidalPy.constants as tidalpy_constants
from TidalPy import Cooling
from TidalPy.Utilities.classes.classes import PhysicsBase, TidalPyBaseClass

# Each entry is the legacy (cooling_flux, boundary_layer_thickness, rayleigh, nusselt) tuple.
_FROZEN_PATH = Path(__file__).parent / "frozen" / "test_cooling_01.npz"
with np.load(_FROZEN_PATH, allow_pickle=False) as _frozen_file:
    _LEGACY_REFERENCE = {key: _frozen_file[key] for key in _frozen_file.files}

# A convecting silicate mantle (high Rayleigh number), MKS.
_INPUTS = {
    "delta_temp": 1000.0,
    "thickness": 1.0e6,
    "gravity": 9.8,
    "density": 3300.0,
    "viscosity": 1.0e21,
    "conductivity": 4.0,
    "diffusivity": 1.0e-6,
    "expansivity": 3.0e-5,
}
_RACR = 1100.0   # default critical Rayleigh number

_MODELS = [("OffCooling", "off"), ("ConvectiveCooling", "convection"), ("ConductiveCooling", "conduction")]


def _inputs(**overrides):
    """The eight calc_cooling inputs in order, with any overrides."""
    return tuple({**_INPUTS, **overrides}.values())


def _binary_copy(model, blank, tmp_path):
    """Save a model to a binary file and load it into a blank model."""
    path = str(tmp_path / "cooling.tpyb")
    model.save_binary(path)
    blank.load_binary(path)
    return blank


# =====================================================================================================================
# Model identity
# =====================================================================================================================
@pytest.mark.parametrize("cls_name,model_name", _MODELS)
def test_model_name(cls_name, model_name):
    assert getattr(Cooling, cls_name)().model_name == model_name


@pytest.mark.parametrize("cls_name,model_name", _MODELS)
def test_make_cooling_returns_rich_subclass(cls_name, model_name):
    assert type(Cooling.make_cooling(model_name)).__name__ == cls_name


@pytest.mark.parametrize("cls_name,model_name", _MODELS)
def test_isinstance_chain(cls_name, model_name):
    model = getattr(Cooling, cls_name)()
    assert isinstance(model, Cooling.CoolingBase)
    assert isinstance(model, PhysicsBase)
    assert isinstance(model, TidalPyBaseClass)


@pytest.mark.parametrize("call,error", [
    (lambda: Cooling.CoolingBase(), TypeError),
    (lambda: Cooling.make_cooling("not_a_real_model"), ValueError),
    (lambda: Cooling.ConvectiveCooling().calc_cooling_vectorize_all(
        *_inputs(delta_temp=np.ones(3), viscosity=np.ones(2))), ValueError),
    (lambda: Cooling.ConvectiveCooling().load_binary("does_not_exist_12345.tpyb"), FileNotFoundError),
], ids=["abstract_base", "unknown_model", "vectorize_size_mismatch", "missing_binary"])
def test_invalid_use_raises(call, error):
    with pytest.raises(error):
        call()


# =====================================================================================================================
# Cooling physics
# =====================================================================================================================
def test_off_cooling():
    """Off gives zero flux, a boundary layer of half the thickness, Ra = 0, and Nu = 1."""
    result = Cooling.OffCooling().calc_cooling(*_inputs())
    assert result.cooling_flux == 0.0
    assert result.boundary_layer_thickness == pytest.approx(0.5 * _INPUTS["thickness"])
    assert result.rayleigh == 0.0
    assert result.nusselt == 1.0


@pytest.mark.parametrize("key,cls_name,overrides", [
    ("conduction", "ConductiveCooling", {}),
    ("convection_high_rayleigh", "ConvectiveCooling", {}),
    # The legacy sub-critical case is not compared: the legacy model floored Nu at 2, which reported twice the
    # conductive flux; the floor is now 1, so a sub-critical layer conducts (test_cooling_fixes_01).
])
def test_matches_legacy(key, cls_name, overrides):
    """All four outputs match the legacy cooling model."""
    result = getattr(Cooling, cls_name)().calc_cooling(*_inputs(**overrides))
    flux, boundary_layer, rayleigh, nusselt = _LEGACY_REFERENCE[key]
    assert result.cooling_flux == pytest.approx(flux, rel=1e-12)
    assert result.boundary_layer_thickness == pytest.approx(boundary_layer, rel=1e-12)
    assert result.rayleigh == pytest.approx(rayleigh, rel=1e-12)
    assert result.nusselt == pytest.approx(nusselt, rel=1e-12)


def test_conduction_is_the_closed_form():
    result = Cooling.ConductiveCooling().calc_cooling(*_inputs())
    assert result.cooling_flux == pytest.approx(_INPUTS["conductivity"] * _INPUTS["delta_temp"] / _INPUTS["thickness"])


def test_convection_regimes():
    """The reference mantle convects (Nu > 2, Ra > Ra_cr); the stiff thin layer is floored at Nu = 1."""
    active = Cooling.ConvectiveCooling().calc_cooling(*_inputs())
    assert active.nusselt > 2.0
    assert active.rayleigh > _RACR
    floored = Cooling.ConvectiveCooling().calc_cooling(*_inputs(delta_temp=100.0, viscosity=1.0e24, thickness=1.0e4))
    assert floored.nusselt == pytest.approx(1.0)


def test_convection_guards_agree_at_the_minimum_thickness():
    """A layer exactly at the minimum thickness is too thin for every output, not for some of them."""
    thickness = tidalpy_constants.min_thickness
    assert thickness > 0.0
    # A viscosity low enough that the layer would convect hard if its thickness were accepted.
    result = Cooling.ConvectiveCooling().calc_cooling(*_inputs(thickness=thickness, viscosity=1.0e-6))
    assert result.rayleigh == 0.0
    assert result.nusselt == 1.0
    assert result.boundary_layer_thickness == thickness
    result = Cooling.ConvectiveCooling().calc_cooling(*_inputs(thickness=2.0 * thickness, viscosity=1.0e-6))
    assert result.rayleigh > 0.0
    assert result.nusselt > 2.0
    assert result.boundary_layer_thickness == pytest.approx(2.0 * thickness / result.nusselt)


def test_convection_without_contrast_keeps_the_floored_boundary_layer():
    """No temperature contrast gives no flux and a boundary layer of thickness / Nu_min."""
    # A whole-planet solve uses this boundary layer as the resistance to neighbors at other temperatures.
    result = Cooling.ConvectiveCooling().calc_cooling(*_inputs(delta_temp=0.0, thickness=5.0e5))
    assert result.cooling_flux == 0.0
    assert result.rayleigh == 0.0
    assert result.nusselt == 1.0
    assert result.boundary_layer_thickness == pytest.approx(5.0e5)


def test_convection_parameters_affect_result():
    """A larger critical Rayleigh number reduces the Nusselt number."""
    base = Cooling.ConvectiveCooling().calc_cooling(*_inputs())
    stiff = Cooling.ConvectiveCooling(critical_rayleigh=1.0e5).calc_cooling(*_inputs())
    assert stiff.nusselt < base.nusselt


def test_cooling_result_container():
    result = Cooling.ConductiveCooling().calc_cooling(*_inputs())
    assert set(result.to_dict()) == {"cooling_flux", "boundary_layer_thickness", "rayleigh", "nusselt"}
    assert tuple(result) == (result.cooling_flux, result.boundary_layer_thickness, result.rayleigh, result.nusselt)
    assert "CoolingResult" in repr(result)


# =====================================================================================================================
# Parameters and factory
# =====================================================================================================================
@pytest.mark.parametrize("use_factory,params", [
    (False, (0.5, 0.25, 2000.0)),
    (True, (0.6, 0.28, 1600.0)),
], ids=["constructor", "factory"])
def test_convective_parameters(use_factory, params):
    if use_factory:
        keys = ("convection_alpha", "convection_beta", "critical_rayleigh")
        convective = Cooling.make_cooling("convection", dict(zip(keys, params)))
    else:
        convective = Cooling.ConvectiveCooling(*params)
    assert convective.convection_alpha == pytest.approx(params[0])
    assert convective.convection_beta == pytest.approx(params[1])
    assert convective.critical_rayleigh == pytest.approx(params[2])


@pytest.mark.parametrize("alias,canonical", [
    ("off",         "off"),
    ("none",        "off"),
    ("OFF",         "off"),
    ("convection",  "convection"),
    ("convective",  "convection"),
    ("Convection",  "convection"),
    ("conduction",  "conduction"),
    ("conductive",  "conduction"),
    ("CONDUCTION",  "conduction"),
])
def test_make_cooling_aliases(alias, canonical):
    assert Cooling.make_cooling(alias).model_name == canonical


def test_make_cooling_adopted_object_is_usable(tmp_path):
    """An object built by the C++ enum factory is usable and round trips."""
    convective = Cooling.make_cooling("convection", {"critical_rayleigh": 1600.0})
    assert convective.critical_rayleigh == pytest.approx(1600.0)
    assert convective.calc_cooling(*_inputs()).cooling_flux > 0.0
    restored = _binary_copy(convective, Cooling.ConvectiveCooling(), tmp_path)
    assert restored.get_config_dict() == convective.get_config_dict()


# =====================================================================================================================
# Config dict and binary round trip
# =====================================================================================================================
@pytest.mark.parametrize("cls_name,keys", [
    ("OffCooling", {"model"}),
    ("ConductiveCooling", {"model"}),
    ("ConvectiveCooling", {"model", "convection_alpha", "convection_beta", "critical_rayleigh",
                           "liquid_convection_alpha", "liquid_convection_beta"}),
])
def test_config_dict_keys(cls_name, keys):
    assert set(getattr(Cooling, cls_name)().get_config_dict()) == keys


def test_save_config_writes_toml(tmp_path):
    import toml
    path = str(tmp_path / "cool.toml")
    Cooling.ConvectiveCooling(0.6, 0.28, 1600.0).save_config(path)
    loaded = toml.load(path)
    assert loaded["model"] == "convection"
    assert loaded["convection_alpha"] == pytest.approx(0.6)
    assert loaded["critical_rayleigh"] == pytest.approx(1600.0)


@pytest.mark.parametrize("name,cls,args", [
    ("off", "OffCooling", ()),
    ("conduction", "ConductiveCooling", ()),
    ("convection", "ConvectiveCooling", (0.6, 0.28, 1600.0)),
])
def test_binary_round_trip(name, cls, args, tmp_path):
    model_class = getattr(Cooling, cls)
    original = model_class(*args)
    restored = _binary_copy(original, model_class(), tmp_path)
    assert restored.model_name == name
    assert restored.get_config_dict() == original.get_config_dict()


# =====================================================================================================================
# Vectorized methods and direct convenience functions
# =====================================================================================================================
_ARRAY_INPUTS = [
    {"delta_temp": np.array([100.0, 500.0, 1000.0, 2000.0])},
    {"viscosity": np.array([1.0e19, 1.0e20, 1.0e21])},
    {"delta_temp": np.array([500.0, 1000.0, 2000.0]), "viscosity": np.array([1.0e20, 1.0e21, 1.0e22])},
]
_ARRAY_IDS = ["delta_temp", "viscosity", "both"]


def _scalar_results(model, array_inputs):
    """calc_cooling at each element of the array inputs."""
    size = len(next(iter(array_inputs.values())))
    return [model.calc_cooling(*_inputs(**{key: values[index] for key, values in array_inputs.items()}))
            for index in range(size)]


@pytest.mark.parametrize("array_inputs", _ARRAY_INPUTS, ids=_ARRAY_IDS)
def test_convenience_vectorizes(array_inputs):
    """Array inputs to the convenience function give float64 arrays that match scalar calls."""
    got = Cooling.convective(*_inputs(**array_inputs))
    singles = _scalar_results(Cooling.ConvectiveCooling(), array_inputs)
    assert isinstance(got.cooling_flux, np.ndarray)
    assert got.cooling_flux.shape == (len(singles),)
    assert got.cooling_flux.dtype == np.float64
    assert got.cooling_flux == pytest.approx(np.array([single.cooling_flux for single in singles]))
    assert got.nusselt == pytest.approx(np.array([single.nusselt for single in singles]))


@pytest.mark.parametrize("method,array_inputs", [
    ("calc_cooling_vectorize_temperature", {"delta_temp": np.array([100.0, 1000.0, 2000.0])}),
    ("calc_cooling_vectorize_viscosity", {"viscosity": np.array([1.0e20, 1.0e21, 1.0e22])}),
    ("calc_cooling_vectorize_all", {"delta_temp": np.array([100.0, 1000.0, 2000.0]),
                                    "viscosity": np.array([1.0e20, 1.0e21, 1.0e22])}),
], ids=["temperature", "viscosity", "all"])
def test_class_vectorize_methods_match_scalar(method, array_inputs):
    convective = Cooling.ConvectiveCooling()
    got = getattr(convective, method)(*_inputs(**array_inputs))
    singles = _scalar_results(convective, array_inputs)
    assert got.cooling_flux == pytest.approx(np.array([single.cooling_flux for single in singles]))


def test_convenience_scalar_matches_class():
    got = Cooling.convective(*_inputs())
    expected = Cooling.ConvectiveCooling().calc_cooling(*_inputs())
    assert isinstance(got.cooling_flux, float)
    assert got.cooling_flux == pytest.approx(expected.cooling_flux)
    assert got.nusselt == pytest.approx(expected.nusselt)


def test_convenience_preserves_2d_shape():
    got = Cooling.convective(*_inputs(delta_temp=np.array([[100.0, 500.0], [1000.0, 2000.0]])))
    assert got.cooling_flux.shape == (2, 2)
    assert got.nusselt.shape == (2, 2)


def test_conductive_off_convenience():
    delta_temp, thickness, conductivity = _INPUTS["delta_temp"], _INPUTS["thickness"], _INPUTS["conductivity"]
    assert Cooling.conductive(delta_temp, thickness, conductivity).cooling_flux == pytest.approx(
        conductivity * delta_temp / thickness)
    off = Cooling.cooling_off(delta_temp, thickness)
    assert off.cooling_flux == 0.0
    assert off.boundary_layer_thickness == pytest.approx(0.5 * thickness)
    assert Cooling.conductive(np.array([100.0, 200.0]), thickness, conductivity).cooling_flux.shape == (2,)
