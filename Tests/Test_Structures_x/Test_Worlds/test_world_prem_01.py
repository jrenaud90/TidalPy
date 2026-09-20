"""
Tests for building a PREM Earth from a PREM-like data file via the world builder
(``build_world("earth_prem")``, which carries ``data_file = "PREM.csv"``).

Confirms that the bundled PREM profile is loaded, the layers are auto-detected with the
liquid outer core flagged for the radial solver, an interpolated EOS is built per layer,
the whole-planet EOS solve converges and reproduces Earth's central pressure / surface
gravity / mass, the interpolated density and static shear/bulk moduli are returned vs
radius (zero shear in the liquid outer core), and the degree-2 Love number is Earth's.
"""

import cmath
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

# Earth's mass and central pressure from the PREM profile itself.
_PREM_MASS = 5.9731e24
_PREM_CENTRAL_PRESSURE = 3.639e11
# Earth's solid-body degree-2 potential Love number (0.298 to 0.302 across studies).
_EARTH_K2 = 0.30


def test_earth_prem_builds_three_layers_with_a_liquid_outer_core():
    world = build_world("earth_prem")
    assert world.name == "Earth-PREM"
    assert world.num_layers == 3
    assert world.all_eos_set is True
    assert [layer.is_solid for layer in world] == [True, False, True]
    assert [layer.is_static for layer in world] == [True, True, True]
    assert [layer.is_tidal for layer in world] == [True, False, True]
    # The detected flags are config keys, so the expanded configuration carries them.
    assert [cfg["is_solid"] for cfg in world.source_config["layers"].values()] == [True, False, True]


def test_earth_prem_eos_solve_converges_and_reproduces_earth():
    world = build_world("earth_prem")
    result = world.solve_eos(G_to_use=G, slices_per_layer=120, verbose=False)
    assert result["success"] is True
    assert result["max_iters_hit"] is False
    assert result["iterations"] <= 10
    # The PREM profile integrates to Earth's mass, surface gravity, and central pressure.
    assert math.isclose(world.surface_gravity_eos, 9.81, rel_tol=0.01)
    assert math.isclose(world.planet_mass_eos, _PREM_MASS, rel_tol=1.0e-3)
    assert math.isclose(world.central_pressure, _PREM_CENTRAL_PRESSURE, rel_tol=0.01)


def test_earth_prem_love_number_is_earths():
    world = build_world("earth_prem")
    world.solve_eos(G_to_use=G, verbose=False)
    result = world.solve_love_numbers(frequency=2.0 * math.pi / 86400.0, degree_l=2)
    assert result["success"] is True, result["message"]
    k2 = world.love_number_k
    assert cmath.isclose(k2, _EARTH_K2, rel_tol=0.01)
    assert abs(k2.imag) < 1.0e-6   # elastic PREM: no dissipation


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
    # earth_prem has 3 layers; override layer_2 (the solid mantle) bulk modulus to a
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
            "layer_1": {"class": "physics", "layer_index": 1, "is_incompressible": True},
            "layer_2": {"class": "solidliquid", "layer_index": 2, "material": {"bulk_modulus_static_pa": 1.0e11}},
        },
    }
    world = build_world(config)
    world.solve_eos(G_to_use=G, verbose=False)
    # The mantle bulk modulus is now the constant override (not the PREM value).
    assert math.isclose(world.get_bulk_modulus(_MANTLE_RADIUS_M), 1.0e11, rel_tol=1e-6)
    # The user tables keep the detected liquid flag of the outer core and can add the other flags.
    assert [layer.is_solid for layer in world] == [True, False, True]
    assert [layer.is_incompressible for layer in world] == [False, True, False]


def test_earth_prem_layer_count_mismatch_raises():
    prem_path = worldpack.resolve_data_file("PREM.csv")
    config = {
        "schema_version": "0.2.0",
        "name": "Earth-PREM-Bad",
        "type": "terrestrial",
        "radius_m": 6371000.0,
        "mass_kg": 5.972e24,
        "data_file": prem_path,
        "layers": {  # only 2 tables but 3 layers detected
            "layer_0": {"class": "solidliquid", "layer_index": 0},
            "layer_1": {"class": "physics", "layer_index": 1},
        },
    }
    with pytest.raises(ValueError, match="layer"):
        build_world(config)
