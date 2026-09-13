"""Direct construction of every concrete physics-model class.

The Cython constructors build their C++ object through the module's ``c_find_*`` factory, naming the
model with an enum member. A wrong enum would build the wrong model and cast the class's non-owning
typed pointer to the wrong type, which is silent rather than fatal, so each class is checked against
the model name it is supposed to produce. Reading every property afterwards exercises that typed
pointer, and ``get_config_dict`` exercises the owning ``unique_ptr``.

The last test pins where a default lives. A ``make_*`` factory given no config must agree with the
direct constructor, because both are supposed to fall through to the same C++ struct default. A factory
that re-introduces its own literal (``config.get(key, 0.3)``) passes until that literal drifts from the
struct, which is exactly the duplication this checks for.
"""
import gc
import inspect
import types

import pytest

from TidalPy.rheology_x import Elastic, Viscous, Maxwell, Voigt, Burgers, Andrade, Sundberg
from TidalPy.cooling_x.cooling import OffCooling, ConductiveCooling, ConvectiveCooling
from TidalPy.radiogenics_x.radiogenics import OffRadiogenics, IsotopeRadiogenics, FixedRadiogenics
from TidalPy.stellar_x.luminosity import FixedLuminosity, MassToLuminosity, PowerLawLuminosity
from TidalPy.rheology_x.rheology import make_rheology
from TidalPy.cooling_x.cooling import make_cooling
from TidalPy.radiogenics_x.radiogenics import make_radiogenics
from TidalPy.stellar_x.luminosity import make_luminosity


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


# Class -> number of readable properties, so the sweep below cannot pass vacuously if a class loses
# its parameter accessors. `model_name` comes from PhysicsBase and is counted here.
_PROPERTY_COUNTS = {
    Elastic:            1,
    Viscous:            1,
    Maxwell:            1,
    Voigt:              3,
    Burgers:            3,
    Andrade:            3,
    Sundberg:           5,
    OffCooling:         1,
    ConductiveCooling:  1,
    ConvectiveCooling:  4,
    OffRadiogenics:     1,
    IsotopeRadiogenics: 8,
    FixedRadiogenics:   4,
    FixedLuminosity:    2,
    MassToLuminosity:   1,
    PowerLawLuminosity: 3,
}


def _cases():
    return [pytest.param(cls, model_name, id=cls.__name__) for cls, model_name in _MODELS]


def _property_names(model_class):
    """Readable attribute names on a Cython cdef class (its properties are getset descriptors)."""
    return sorted(name for name, attribute in inspect.getmembers(model_class)
                  if isinstance(attribute, types.GetSetDescriptorType) and not name.startswith("_"))


@pytest.mark.parametrize("model_class,expected_model_name", _cases())
def test_default_constructor_builds_the_named_model(model_class, expected_model_name):
    """A no-argument constructor resolves to its own model, not a sibling in the same family."""
    model = model_class()
    assert model.model_name == expected_model_name
    assert model.get_config_dict()["model"] == expected_model_name


@pytest.mark.parametrize("model_class,expected_model_name", _cases())
def test_every_property_is_readable(model_class, expected_model_name):
    """Reading each property exercises the class's non-owning typed pointer into the C++ object."""
    model = model_class()
    property_names = _property_names(model_class)
    assert len(property_names) == _PROPERTY_COUNTS[model_class], property_names
    for name in property_names:
        getattr(model, name)   # a mis-cast or dangling typed pointer shows up here


@pytest.mark.parametrize("model_class,expected_model_name", _cases())
def test_object_outlives_a_collection(model_class, expected_model_name):
    """The owning unique_ptr keeps the C++ object alive; the typed pointer must not dangle."""
    model = model_class()
    config = model.get_config_dict()
    gc.collect()
    assert model.get_config_dict() == config


# Class -> the factory that builds the same model by name. Every one of these families had a factory
# that repeated its struct defaults as literals; see the module docstring.
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
    factory = _FACTORIES[model_class]
    assert factory(expected_model_name, {}).get_config_dict() == factory(
        expected_model_name).get_config_dict()
