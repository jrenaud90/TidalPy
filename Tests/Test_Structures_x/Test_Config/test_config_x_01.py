"""
Tests for the new ``_x`` configuration (``TidalPy_Configs_x.toml`` / ``TidalPy.config_x``)
built by ``TidalPy.defaultc_x`` and loaded during initialization.

Covers that the ``_x`` config is loaded with the expected numerical and per-material
layer sections, and that ``update_constants_x`` populates the shared C++ config
singleton (mirrored on the ``TidalPy.constants`` module globals) from it.

Requires the Cython extensions to be compiled first::

    uv pip install -v <repo_root>
"""

import math

import pytest

import TidalPy
from TidalPy.structures_x.configs.toml_loader import MATERIAL_TYPES


def test_config_x_loaded():
    assert TidalPy.config_x is not None
    assert TidalPy.config_x["schema_version"] == "0.2.0"


def test_config_x_has_numerical_section():
    numerical = TidalPy.config_x["numerical"]
    for key in ("minimum_frequency", "maximum_frequency", "min_spin_orbit_diff",
                "minimum_viscosity", "minimum_modulus", "minimum_layer_thickness",
                "test_constant"):
        assert key in numerical


@pytest.mark.parametrize("material_type", list(MATERIAL_TYPES))
def test_config_x_has_each_material_block(material_type):
    layers = TidalPy.config_x["layers"]
    assert material_type in layers
    # Every material block carries an EOS default.
    assert "eos" in layers[material_type]


def test_mantle_rock_defaults_present():
    mantle = TidalPy.config_x["layers"]["mantle_rock"]
    assert mantle["shear_rheology"]["model"] == "andrade"
    assert mantle["shear_rheology"]["zeta"] == 1.0
    assert mantle["cooling"]["model"] == "convection"
    assert mantle["radiogenics"]["model"] == "isotope"


def test_update_constants_x_populated_singleton():
    # update_constants_x runs during initialization and feeds the _x numerical
    # settings into the constants module globals.
    from TidalPy import constants
    numerical = TidalPy.config_x["numerical"]
    assert math.isclose(constants.min_viscosity, numerical["minimum_viscosity"])
    assert math.isclose(constants.min_modulus, numerical["minimum_modulus"])
    assert math.isclose(constants.min_thickness, numerical["minimum_layer_thickness"])
    assert math.isclose(constants.test_constant, numerical["test_constant"])
