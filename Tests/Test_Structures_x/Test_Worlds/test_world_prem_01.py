"""
Tests for building a PREM Earth from a PREM-like data file via the world builder
(``build_world("earth_prem")``, which carries ``data_file = "PREM.csv"``).

Confirms that the bundled PREM profile is loaded, the layers are auto-detected, an
interpolated EOS is built per layer, the whole-planet EOS solve reproduces Earth's
central pressure / surface gravity / mass, and the interpolated density and static
shear/bulk moduli are returned vs radius (zero shear in the liquid outer core).
"""

import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.structures_x import build_world
from TidalPy.structures_x.configs import prem, worldpack


def _prem_arrays():
    return prem.load_prem_arrays(worldpack.resolve_data_file("PREM.csv"))


# A radius well inside the (solid) lower mantle and one inside the (liquid) outer core.
_MANTLE_RADIUS_M = 5.0e6
_OUTER_CORE_RADIUS_M = 2.5e6


def test_earth_prem_builds_four_layers():
    world = build_world("earth_prem")
    assert world.name == "Earth-PREM"
    assert world.num_layers == 4
    assert world.all_eos_set is True


def test_earth_prem_eos_solve_reproduces_earth():
    world = build_world("earth_prem")
    result = world.solve_eos(G_to_use=G, slices_per_layer=120, verbose=False)
    assert result["success"] is True
    assert world.central_pressure > 0.0
    # Earth-like bulk properties from the PREM profile.
    assert math.isclose(world.surface_gravity_eos, 9.81, rel_tol=0.15)
    assert math.isclose(world.planet_mass_eos, 5.972e24, rel_tol=0.15)


def test_earth_prem_density_interpolated():
    arrays = _prem_arrays()
    world = build_world("earth_prem")
    world.solve_eos(G_to_use=G, verbose=False)
    expected = np.interp(_MANTLE_RADIUS_M, arrays["radius_m"], arrays["density_kg_m3"])
    assert math.isclose(world.get_density(_MANTLE_RADIUS_M), expected, rel_tol=0.05)


def test_earth_prem_shear_modulus_interpolated():
    arrays = _prem_arrays()
    world = build_world("earth_prem")
    world.solve_eos(G_to_use=G, verbose=False)
    # Solid mantle: shear modulus follows the interpolated PREM profile (rho*Vs^2).
    expected_shear = np.interp(_MANTLE_RADIUS_M, arrays["radius_m"], arrays["shear_modulus_pa"])
    assert expected_shear > 0.0
    assert math.isclose(world.get_shear_modulus(_MANTLE_RADIUS_M), expected_shear, rel_tol=0.10)
    # Liquid outer core: zero shear modulus.
    assert abs(world.get_shear_modulus(_OUTER_CORE_RADIUS_M)) < 1.0e3


def test_earth_prem_bulk_modulus_interpolated():
    arrays = _prem_arrays()
    world = build_world("earth_prem")
    world.solve_eos(G_to_use=G, verbose=False)
    expected_bulk = np.interp(_MANTLE_RADIUS_M, arrays["radius_m"], arrays["bulk_modulus_pa"])
    assert expected_bulk > 0.0
    assert math.isclose(world.get_bulk_modulus(_MANTLE_RADIUS_M), expected_bulk, rel_tol=0.10)


def test_earth_prem_toml_override_of_modulus(tmp_path):
    """A user layer table overriding a constant modulus replaces the PREM array."""
    import toml
    # earth_prem has 4 layers; override layer_2 (the solid mantle) bulk modulus to a
    # constant. Provide one table per detected layer (inner to outer).
    prem_path = worldpack.resolve_data_file("PREM.csv")
    config = {
        "schema_version": "0.2.0",
        "name": "Earth-PREM-Override",
        "type": "terrestrial",
        "radius_m": 6371000.0,
        "mass_kg": 5.972e24,
        "data_file": prem_path,
        "layers": {
            "layer_0": {"class": "solidliquid", "layer_index": 0},
            "layer_1": {"class": "physics", "layer_index": 1},
            "layer_2": {"class": "solidliquid", "layer_index": 2,
                        "bulk_modulus_static_pa": 1.0e11},
            "layer_3": {"class": "physics", "layer_index": 3},
        },
    }
    world = build_world(config)
    world.solve_eos(G_to_use=G, verbose=False)
    # The mantle bulk modulus is now the constant override (not the PREM value).
    assert math.isclose(world.get_bulk_modulus(_MANTLE_RADIUS_M), 1.0e11, rel_tol=1e-6)


def test_earth_prem_layer_count_mismatch_raises():
    prem_path = worldpack.resolve_data_file("PREM.csv")
    config = {
        "schema_version": "0.2.0",
        "name": "Earth-PREM-Bad",
        "type": "terrestrial",
        "radius_m": 6371000.0,
        "mass_kg": 5.972e24,
        "data_file": prem_path,
        "layers": {  # only 2 tables but 4 layers detected
            "layer_0": {"class": "solidliquid", "layer_index": 0},
            "layer_1": {"class": "physics", "layer_index": 1},
        },
    }
    with pytest.raises(ValueError, match="layer"):
        build_world(config)
