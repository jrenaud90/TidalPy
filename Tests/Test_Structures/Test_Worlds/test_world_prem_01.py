"""PREM Earth built from a radial data file: detected layers, EOS solve, interpolated properties, Love number,
and layer tables that refine the detected layers.
"""

import cmath
import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Structures import build_world
from TidalPy.Structures.configs import data_file, worldpack


def _prem_arrays():
    return data_file.load_radial_data(worldpack.resolve_data_file("PREM.csv"))


# Radii inside the solid lower mantle and the liquid outer core.
_MANTLE_RADIUS_M = 5.0e6
_OUTER_CORE_RADIUS_M = 2.5e6

# Earth's mass and central pressure from the PREM profile itself.
_PREM_MASS = 5.9731e24
_PREM_CENTRAL_PRESSURE = 3.639e11
# Earth's solid-body degree-2 potential Love number (0.298 to 0.302 across studies).
_EARTH_K2 = 0.30


def _solved_prem():
    world = build_world("earth_prem")
    world.solve_eos(G_to_use=G, verbose=False)
    return world


def test_earth_prem_builds_three_layers_with_a_liquid_outer_core():
    """earth_prem detects three layers with a liquid, non-tidal outer core."""
    world = build_world("earth_prem")
    assert world.name == "Earth-PREM"
    assert world.num_layers == 3
    assert world.all_eos_set is True
    assert [layer.is_solid for layer in world] == [True, False, True]
    assert [layer.is_static for layer in world] == [True, True, True]
    assert [layer.is_tidal for layer in world] == [True, False, True]
    assert [cfg["is_solid"] for cfg in world.source_config["layers"].values()] == [True, False, True]


def test_earth_prem_eos_solve_converges_and_reproduces_earth():
    """The EOS solve converges to Earth's surface gravity, mass, and central pressure."""
    world = build_world("earth_prem")
    result = world.solve_eos(G_to_use=G, slices_per_layer=120, verbose=False)
    assert result["success"] is True
    assert result["max_iters_hit"] is False
    assert result["iterations"] <= 10
    assert math.isclose(world.surface_gravity_eos, 9.81, rel_tol=0.01)
    assert math.isclose(world.planet_mass_eos, _PREM_MASS, rel_tol=1.0e-3)
    assert math.isclose(world.central_pressure, _PREM_CENTRAL_PRESSURE, rel_tol=0.01)


def test_earth_prem_love_number_is_earths():
    """The elastic PREM k2 is Earth's, with no dissipation."""
    world = _solved_prem()
    result = world.solve_love_numbers(frequency=2.0 * math.pi / 86400.0, degree_l=2)
    assert result["success"] is True, result["message"]
    k2 = world.love_number_k
    assert cmath.isclose(k2, _EARTH_K2, rel_tol=0.01)
    assert abs(k2.imag) < 1.0e-6


@pytest.mark.parametrize(
    "getter_name, column, rel_tol",
    [
        ("get_density", "density_kg_m3", 0.05),
        ("get_shear_modulus", "shear_modulus_pa", 0.10),
        ("get_bulk_modulus", "bulk_modulus_pa", 0.10),
    ],
    ids=["density", "shear_modulus", "bulk_modulus"],
)
def test_earth_prem_mantle_property_interpolated(getter_name, column, rel_tol):
    """The mantle's density and static moduli follow the interpolated PREM profile."""
    arrays = _prem_arrays()
    world = _solved_prem()
    expected = np.interp(_MANTLE_RADIUS_M, arrays["radius_m"], arrays[column])
    assert expected > 0.0
    assert math.isclose(getattr(world, getter_name)(_MANTLE_RADIUS_M), expected, rel_tol=rel_tol)


def test_earth_prem_outer_core_has_no_shear_modulus():
    """The liquid outer core has zero shear modulus."""
    world = _solved_prem()
    assert abs(world.get_shear_modulus(_OUTER_CORE_RADIUS_M)) < 1.0e3


def _prem_config(**extra):
    config = {
        "schema_version": "0.2.0",
        "name": "Earth-PREM-Test",
        "type": "terrestrial",
        "radius_m": 6371000.0,
        "mass_kg": 5.972e24,
        "data_file": worldpack.resolve_data_file("PREM.csv"),
    }
    config.update(extra)
    return config


def test_earth_prem_toml_override_of_modulus():
    """A layer table's constant modulus replaces that layer's profile values; detected flags are kept."""
    world = build_world(_prem_config(layers={
        "layer_0": {"class": "solidliquid", "layer_index": 0},
        "layer_1": {"class": "base", "layer_index": 1, "is_incompressible": True},
        "layer_2": {"class": "solidliquid", "layer_index": 2, "material": {"bulk_modulus_static_pa": 1.0e11}},
    }))
    world.solve_eos(G_to_use=G, verbose=False)
    assert math.isclose(world.get_bulk_modulus(_MANTLE_RADIUS_M), 1.0e11, rel_tol=1e-6)
    assert [layer.is_solid for layer in world] == [True, False, True]
    assert [layer.is_incompressible for layer in world] == [False, True, False]


def test_one_layer_table_refines_one_layer():
    """One layer table adds a rheology to one layer; the other layers need no table."""
    world = build_world(_prem_config(layers={
        "mantle": {
            "layer_index": 2,
            "shear_rheology": {"model": "maxwell"},
            "material": {"shear_viscosity_static_pas": 1.0e21},
        },
    }))
    assert world.num_layers == 3
    assert [layer.name for layer in world] == ["layer_0", "layer_1", "mantle"]
    assert [layer.shear_rheology_set for layer in world] == [False, False, True]
    assert [layer.is_solid for layer in world] == [True, False, True]
    world.solve_eos(G_to_use=G, verbose=False)
    arrays = _prem_arrays()
    expected = np.interp(_MANTLE_RADIUS_M, arrays["radius_m"], arrays["density_kg_m3"])
    assert math.isclose(world.get_density(_MANTLE_RADIUS_M), expected, rel_tol=0.05)
    result = world.solve_love_numbers(frequency=2.0 * math.pi / 86400.0, degree_l=2)
    assert result["success"] is True, result["message"]
    assert abs(world.love_number_k.imag) > 0.0


def _prem_config_without_radius():
    config = _prem_config()
    del config["radius_m"]
    return config


@pytest.mark.parametrize(
    "make_config, match",
    [
        pytest.param(
            lambda: _prem_config(layers={"mantle": {"shear_rheology": {"model": "maxwell"}}}),
            "layer_index",
            id="table_without_layer_index",
        ),
        pytest.param(
            lambda: _prem_config(layers={"layer_7": {"class": "solidliquid"}}),
            "layer_index 7",
            id="table_out_of_range",
        ),
        pytest.param(
            lambda: _prem_config(layers={"core": {"layer_index": 0}, "inner": {"layer_index": 0}}),
            "both",
            id="two_tables_one_layer",
        ),
        # A table named layer_0 that refines layer 2 would displace the real layer_0.
        pytest.param(
            lambda: _prem_config(layers={"layer_0": {"layer_index": 2}}),
            "more than one",
            id="table_named_after_another_layer",
        ),
        pytest.param(
            lambda: _prem_config(data={"radius_km": [0.0, 1.0], "density": [1.0e3, 1.0e3],
                                       "vp": [1.0e4, 1.0e4], "vs": [0.0, 0.0]}),
            "one or the other",
            id="data_and_data_file",
        ),
        pytest.param(_prem_config_without_radius, "radius_m", id="profile_without_radius"),
    ],
)
def test_invalid_profile_config_raises(make_config, match):
    """An invalid profile configuration raises a ValueError naming the problem."""
    with pytest.raises(ValueError, match=match):
        build_world(make_config())


def test_a_world_can_be_built_from_arrays_in_memory():
    """A profile handed to build_world as in-memory arrays builds the same Earth."""
    arrays = _prem_arrays()
    world = build_world({
        "schema_version": "0.2.0",
        "name": "Earth-PREM-Arrays",
        "type": "terrestrial",
        "radius_m": 6371000.0,
        "mass_kg": 5.972e24,
        "data": {
            "radius_m":      arrays["radius_m"],
            "density_kg_m3": arrays["density_kg_m3"],
            "vp_m_s":        arrays["vp_m_s"],
            "vs_m_s":        arrays["vs_m_s"],
        },
    })
    assert world.num_layers == 3
    assert [layer.is_solid for layer in world] == [True, False, True]
    world.solve_eos(G_to_use=G, verbose=False)
    assert math.isclose(world.planet_mass_eos, _PREM_MASS, rel_tol=1.0e-3)
    # The arrays became the layers' materials, so the config does not carry them twice.
    assert "data" not in world.source_config
    result = world.solve_love_numbers(frequency=2.0 * math.pi / 86400.0, degree_l=2)
    assert result["success"] is True, result["message"]
    assert cmath.isclose(world.love_number_k, _EARTH_K2, rel_tol=0.01)


def test_a_profile_without_viscosities_is_elastic():
    """A profile with no viscosity column gives an elastic body: no rheology, viscosity, or melt."""
    world = build_world("earth_prem")
    for layer in world:
        assert layer.shear_rheology_set is False
        assert layer.bulk_rheology_set is False
    world.solve_eos(G_to_use=G, verbose=False)
    assert math.isnan(world.get_shear_viscosity(_MANTLE_RADIUS_M))
    assert math.isnan(world.get_bulk_viscosity(_MANTLE_RADIUS_M))
    assert world.get_melt_fraction(_MANTLE_RADIUS_M) == 0.0
    result = world.solve_love_numbers(frequency=2.0 * math.pi / 86400.0, degree_l=2)
    assert result["success"] is True, result["message"]
    assert abs(world.love_number_k.imag) < 1.0e-6
