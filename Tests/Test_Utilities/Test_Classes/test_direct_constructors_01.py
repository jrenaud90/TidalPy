"""Tests for direct construction of every concrete physics model class and agreement with its factory defaults."""
import gc
import inspect
import types

import pytest

from TidalPy.Rheology import Elastic, Viscous, Maxwell, Voigt, Burgers, Andrade, Sundberg
from TidalPy.Cooling.cooling import OffCooling, ConductiveCooling, ConvectiveCooling
from TidalPy.Radiogenics.radiogenics import OffRadiogenics, IsotopeRadiogenics, FixedRadiogenics
from TidalPy.Stellar.luminosity import FixedLuminosity, MassToLuminosity, PowerLawLuminosity
from TidalPy.Rheology.rheology import make_rheology
from TidalPy.Cooling.cooling import make_cooling
from TidalPy.Radiogenics.radiogenics import make_radiogenics
from TidalPy.Stellar.luminosity import make_luminosity


# Class -> the model name its C++ object must report.
_MODELS = [
    (Elastic,             "elastic"),
    (Viscous,             "viscous"),
    (Maxwell,             "maxwell"),
    (Voigt,               "voigt"),
    (Burgers,             "burgers"),
    (Andrade,             "andrade"),
    (Sundberg,            "sundberg"),
    (OffCooling,          "off"),
    (ConductiveCooling,   "conduction"),
    (ConvectiveCooling,   "convection"),
    (OffRadiogenics,      "off"),
    (IsotopeRadiogenics,  "isotope"),
    (FixedRadiogenics,    "fixed"),
    (FixedLuminosity,     "fixed"),
    (MassToLuminosity,    "mass_to_luminosity"),
    (PowerLawLuminosity,  "power_law"),
]


# Class -> readable property count, so the sweep cannot pass vacuously. Each count includes the two properties every
# model inherits, model_name and parameters.
_PROPERTY_COUNTS = {
    Elastic:            2,
    Viscous:            2,
    Maxwell:            2,
    Voigt:              4,
    Burgers:            4,
    Andrade:            4,
    Sundberg:           6,
    OffCooling:         2,
    ConductiveCooling:  2,
    ConvectiveCooling:  5,
    OffRadiogenics:     2,
    IsotopeRadiogenics: 9,
    FixedRadiogenics:   5,
    FixedLuminosity:    3,
    MassToLuminosity:   2,
    PowerLawLuminosity: 4,
}


def _cases():
    return [pytest.param(cls, model_name, id=cls.__name__) for cls, model_name in _MODELS]


def _property_names(model_class):
    """Readable property names of a Cython cdef class (its getset descriptors)."""
    return sorted(name for name, attribute in inspect.getmembers(model_class)
                  if isinstance(attribute, types.GetSetDescriptorType) and not name.startswith("_"))


@pytest.mark.parametrize("model_class,expected_model_name", _cases())
def test_default_constructor_builds_the_named_model(model_class, expected_model_name):
    """A no argument constructor builds its own model, not a sibling in the same family."""
    model = model_class()
    assert model.model_name == expected_model_name
    assert model.get_config_dict()["model"] == expected_model_name


@pytest.mark.parametrize("model_class,expected_model_name", _cases())
def test_every_property_is_readable(model_class, expected_model_name):
    """Every property reads without error through the class's typed pointer."""
    model = model_class()
    property_names = _property_names(model_class)
    assert len(property_names) == _PROPERTY_COUNTS[model_class], property_names
    for name in property_names:
        getattr(model, name)


@pytest.mark.parametrize("model_class,expected_model_name", _cases())
def test_object_outlives_a_collection(model_class, expected_model_name):
    """The C++ object survives a garbage collection."""
    model = model_class()
    config = model.get_config_dict()
    gc.collect()
    assert model.get_config_dict() == config


# Class -> the factory that builds the same model by name. These factories once repeated the C++ struct defaults
# as literals that could drift, which the factory tests below guard against.
_FACTORIES = {
    Elastic:            make_rheology,
    Viscous:            make_rheology,
    Maxwell:            make_rheology,
    Voigt:              make_rheology,
    Burgers:            make_rheology,
    Andrade:            make_rheology,
    Sundberg:           make_rheology,
    OffCooling:         make_cooling,
    ConductiveCooling:  make_cooling,
    ConvectiveCooling:  make_cooling,
    OffRadiogenics:     make_radiogenics,
    FixedRadiogenics:   make_radiogenics,
    FixedLuminosity:    make_luminosity,
    MassToLuminosity:   make_luminosity,
    PowerLawLuminosity: make_luminosity,
}


@pytest.mark.parametrize("model_class,expected_model_name",
                         [case for case in _cases() if case.values[0] in _FACTORIES])
def test_factory_default_matches_the_constructor(model_class, expected_model_name):
    """A factory given no config falls through to the same C++ struct defaults as the constructor."""
    factory = _FACTORIES[model_class]
    from_factory = factory(expected_model_name).get_config_dict()
    from_constructor = model_class().get_config_dict()
    assert from_factory == from_constructor


@pytest.mark.parametrize("model_class,expected_model_name",
                         [case for case in _cases() if case.values[0] in _FACTORIES])
def test_factory_empty_config_matches_no_config(model_class, expected_model_name):
    """A factory given an empty config builds the same model as one given no config."""
    factory = _FACTORIES[model_class]
    assert factory(expected_model_name, {}).get_config_dict() == factory(
        expected_model_name).get_config_dict()
