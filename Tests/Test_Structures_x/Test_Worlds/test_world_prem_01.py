"""
Tests for building a PREM Earth from a radial data file via the world builder
(``build_world("earth_prem")``, which carries ``data_file = "PREM.csv"``).

Confirms that the bundled PREM profile is loaded, its layers are detected with the liquid outer core
flagged for the radial solver, each layer's slice of the profile becomes its material (an
interpolated EOS), the whole-planet EOS solve converges and reproduces Earth's central pressure /
surface gravity / mass, the interpolated density and static moduli are returned vs radius (zero
shear in the liquid outer core), and the degree-2 Love number is Earth's.

Also covers what a profile does not describe: it names no rheology, cooling model or radiogenics, so
those still come from layer tables, and a table refines one detected layer without the others
needing tables of their own.
"""

import cmath
import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.structures_x import build_world
from TidalPy.structures_x.configs import data_file, worldpack


def _prem_arrays():
    return data_file.load_radial_data(worldpack.resolve_data_file("PREM.csv"))


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
    """A layer table overriding a constant modulus replaces that layer's array from the profile."""
    world = build_world(_prem_config(layers={
        "layer_0": {"class": "solidliquid", "layer_index": 0},
        "layer_1": {"class": "physics", "layer_index": 1, "is_incompressible": True},
        "layer_2": {"class": "solidliquid", "layer_index": 2, "material": {"bulk_modulus_static_pa": 1.0e11}},
    }))
    world.solve_eos(G_to_use=G, verbose=False)
    # The mantle bulk modulus is now the constant override (not the PREM value).
    assert math.isclose(world.get_bulk_modulus(_MANTLE_RADIUS_M), 1.0e11, rel_tol=1e-6)
    # The tables keep the detected liquid flag of the outer core and can add the other flags.
    assert [layer.is_solid for layer in world] == [True, False, True]
    assert [layer.is_incompressible for layer in world] == [False, True, False]


def test_one_layer_table_refines_one_layer():
    """A profile holds no rheology, so it is named in a layer table; the other layers need no table."""
    world = build_world(_prem_config(layers={
        "mantle": {
            "layer_index": 2,
            "shear_rheology": {"model": "maxwell"},
            "material": {"shear_viscosity_static_pas": 1.0e21},
        },
    }))
    assert world.num_layers == 3
    # The refined layer takes the name of its table; the others keep the detected ones.
    assert [layer.name for layer in world] == ["layer_0", "layer_1", "mantle"]
    assert [layer.shear_rheology_set for layer in world] == [False, False, True]
    # The detected material and flags survive the refinement.
    assert [layer.is_solid for layer in world] == [True, False, True]
    world.solve_eos(G_to_use=G, verbose=False)
    arrays = _prem_arrays()
    expected = np.interp(_MANTLE_RADIUS_M, arrays["radius_m"], arrays["density_kg_m3"])
    assert math.isclose(world.get_density(_MANTLE_RADIUS_M), expected, rel_tol=0.05)
    # With a viscosity and a Maxwell rheology the mantle now dissipates.
    result = world.solve_love_numbers(frequency=2.0 * math.pi / 86400.0, degree_l=2)
    assert result["success"] is True, result["message"]
    assert abs(world.love_number_k.imag) > 0.0


def test_a_layer_table_must_say_which_layer_it_refines():
    with pytest.raises(ValueError, match="layer_index"):
        build_world(_prem_config(layers={"mantle": {"shear_rheology": {"model": "maxwell"}}}))


def test_a_layer_table_out_of_range_raises():
    with pytest.raises(ValueError, match="layer_index 7"):
        build_world(_prem_config(layers={"layer_7": {"class": "solidliquid"}}))


def test_two_layer_tables_may_not_refine_the_same_layer():
    with pytest.raises(ValueError, match="both"):
        build_world(_prem_config(layers={
            "core":  {"layer_index": 0},
            "inner": {"layer_index": 0},
        }))


def test_a_layer_table_may_not_be_named_after_another_layer():
    """Naming a table 'layer_0' while refining layer 2 would otherwise displace the real layer_0."""
    with pytest.raises(ValueError, match="more than one"):
        build_world(_prem_config(layers={"layer_0": {"layer_index": 2}}))


def test_a_world_takes_its_profile_from_one_source():
    with pytest.raises(ValueError, match="one or the other"):
        build_world(_prem_config(data={"radius_km": [0.0, 1.0], "density": [1.0e3, 1.0e3],
                                       "vp": [1.0e4, 1.0e4], "vs": [0.0, 0.0]}))


def test_a_profile_needs_the_world_radius():
    config = _prem_config()
    del config["radius_m"]
    with pytest.raises(ValueError, match="radius_m"):
        build_world(config)


def test_a_world_can_be_built_from_arrays_in_memory():
    """The build_world equivalent of a data file: the profile handed over as arrays."""
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
    # The arrays are the layers' materials now, so the config does not carry them a second time.
    assert "data" not in world.source_config
    result = world.solve_love_numbers(frequency=2.0 * math.pi / 86400.0, degree_l=2)
    assert result["success"] is True, result["message"]
    assert cmath.isclose(world.love_number_k, _EARTH_K2, rel_tol=0.01)


def test_a_profile_without_viscosities_is_elastic():
    """No viscosity column means an elastic body: no rheology, no viscosity, no melting."""
    world = build_world("earth_prem")
    for layer in world:
        assert layer.shear_rheology_set is False
        assert layer.bulk_rheology_set is False
    world.solve_eos(G_to_use=G, verbose=False)
    # Nothing supplies a viscosity: no column in the profile, no model built for it, no constant.
    assert math.isnan(world.get_shear_viscosity(_MANTLE_RADIUS_M))
    assert math.isnan(world.get_bulk_viscosity(_MANTLE_RADIUS_M))
    # No partial-melt model, so nothing is molten and the static moduli stand unreduced.
    assert world.get_melt_fraction(_MANTLE_RADIUS_M) == 0.0
    result = world.solve_love_numbers(frequency=2.0 * math.pi / 86400.0, degree_l=2)
    assert result["success"] is True, result["message"]
    assert abs(world.love_number_k.imag) < 1.0e-6   # elastic: no dissipation
