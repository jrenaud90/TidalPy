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


# PREM's discontinuities inside the mantle and crust [m], where the detected layers break.
_PREM_MANTLE_BREAKS = (3630.0e3, 5600.0e3, 5701.0e3, 5771.0e3, 5971.0e3, 6151.0e3, 6291.0e3, 6346.6e3, 6356.0e3)


def test_earth_prem_builds_its_layers_with_a_liquid_outer_core():
    """earth_prem detects the inner core, a liquid, non-tidal outer core, and ten mantle and crust layers between
    PREM's discontinuities."""
    world = build_world("earth_prem")
    assert world.name == "Earth-PREM"
    assert world.num_layers == 12
    assert world.all_materials_set is True
    assert [layer.radius_outer for layer in world] == pytest.approx(
        [1221.5e3, 3480.0e3, *_PREM_MANTLE_BREAKS, 6371.0e3])
    assert [not layer.is_liquid for layer in world] == [True, False] + [True] * 10
    assert all(layer.is_static for layer in world)
    assert [layer.use_tides for layer in world] == [True, False] + [True] * 10
    # A detected liquid layer gets a liquid-only material.
    phases = [list(cfg["material"]) for cfg in world.source_config["layers"].values()]
    assert phases == [["solid"], ["liquid"]] + [["solid"]] * 10


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


def test_earth_prem_solves_in_two_structure_integrations():
    """PREM's density is tabulated in radius, so the surface pressure falls one for one with the central pressure: the
    first unit-slope step lands on the root to rounding, and the integration after it keeps its output."""
    world = build_world("earth_prem")
    result = world.solve_eos(G_to_use=G)
    assert result["success"] is True
    assert result["structure_integrations"] == 2
    assert result["pressure_error"] < 1.0e-12 * result["central_pressure"]
    assert world.solve_eos(G_to_use=G)["structure_integrations"] == 1


@pytest.mark.parametrize("settings, tolerance", [
    ({}, 5.0e-8),                                                    # the defaults' error (1.5e-8 at the center)
    (dict(rtol=1.0e-12, atol=1.0e-16, pressure_tol=1.0e-9), 5.0e-10),
], ids=["defaults", "tight"])
def test_earth_prem_pressure_is_hydrostatic(settings, tolerance):
    """The pressure, left out of the step control for a density that does not depend on it, still integrates rho g:
    at every layer boundary it matches the surface pressure plus a fine trapezoid quadrature of the solved density
    times gravity, to the solve's accuracy on the central-pressure scale, and closer as the tolerances tighten."""
    world = build_world("earth_prem")
    world.solve_eos(G_to_use=G, **settings)
    pressure_above = 0.0
    for layer in reversed(list(world)):
        radii = np.linspace(layer.radius_inner, layer.radius_outer, 20001)
        integrand = layer.get_density(radii) * layer.get_gravity(radii)
        pressure_above += np.sum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(radii))
        assert abs(layer.get_pressure(layer.radius_inner) - pressure_above) < tolerance * world.central_pressure


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
    """A layer table's constant shear modulus replaces that layer's profile values; detected flags are kept."""
    world = build_world(_prem_config(layers={
        "layer_0": {"layer_index": 0},
        "layer_1": {"layer_index": 1, "is_incompressible": True},
        "mantle": {
            "radius_range_m": [3480.0e3, 6371.0e3],
            "material": {"solid": {"shear_modulus": {"model": "constant", "shear_modulus_pa": 1.0e11}}}},
    }))
    world.solve_eos(G_to_use=G, verbose=False)
    assert math.isclose(world.get_shear_modulus(_MANTLE_RADIUS_M), 1.0e11, rel_tol=1e-6)
    assert [not layer.is_liquid for layer in world] == [True, False] + [True] * 10
    assert [layer.is_incompressible for layer in world] == [False, True] + [False] * 10


@pytest.mark.parametrize("slot, key, value, getter", [
    ("eos", "bulk_modulus_pa", 1.5e11, "get_bulk_modulus"),
    ("eos", "density_kg_m3", 4500.0, "get_density"),
    ("shear_modulus", "shear_modulus_pa", 1.5e11, "get_shear_modulus"),
])
def test_a_single_value_holds_across_a_profile_layer(slot, key, value, getter):
    """One number in place of a profile column holds across the layer; the layer's other columns stay the file's."""
    reference = build_world(_prem_config())
    reference.solve_eos(G_to_use=G, verbose=False)
    world = build_world(_prem_config(layers={
        "mantle": {"radius_range_m": [3480.0e3, 6371.0e3], "material": {"solid": {slot: {key: value}}}}}))
    world.solve_eos(G_to_use=G, verbose=False)
    radii = np.array([4.0e6, _MANTLE_RADIUS_M, 6.2e6])
    np.testing.assert_array_equal(np.asarray(getattr(world, getter)(radii)), value)
    untouched = "get_shear_modulus" if getter != "get_shear_modulus" else "get_bulk_modulus"
    if key != "density_kg_m3":
        np.testing.assert_allclose(
            np.asarray(getattr(world, untouched)(radii)), np.asarray(getattr(reference, untouched)(radii)),
            rtol=1.0e-12)


def test_one_layer_table_refines_one_layer():
    """One layer table adds a rheology to one layer; the other layers need no table."""
    world = build_world(_prem_config(layers={
        "mantle": {
            "layer_index": 2,
            "shear_rheology": {"model": "maxwell"},
            "material": {"solid": {"shear_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e21}}},
        },
    }))
    assert world.num_layers == 12
    assert [layer.name for layer in world] == ["layer_0", "layer_1", "mantle"] + [f"layer_{i}" for i in range(3, 12)]
    assert [layer.shear_rheology is not None for layer in world] == [False, False, True] + [False] * 9
    assert [not layer.is_liquid for layer in world] == [True, False] + [True] * 10
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
            lambda: _prem_config(layers={"layer_17": {"is_static": True}}),
            "layer_index 17",
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
        pytest.param(
            lambda: _prem_config(layers={"mantle": {"radius_range_m": [3480.0e3, 6000.0e3]}}),
            "not a boundary",
            id="range_end_off_a_boundary"),
        pytest.param(
            lambda: _prem_config(layers={"mantle": {"radius_range_m": [3480.0e3, 6371.0e3], "layer_index": 2}}),
            "in a radius range or one detected layer",
            id="range_and_index"),
        pytest.param(
            lambda: _prem_config(layers={"mantle": {"radius_range_m": [6371.0e3, 3480.0e3]}}),
            "two increasing radii",
            id="range_not_increasing"),
        pytest.param(
            lambda: _prem_config(layers={"mantle": {"radius_range_m": [3480.0e3, 6371.0e3]},
                                         "layer_5": {"is_static": True}}),
            "both",
            id="range_overlapping_another_table"),
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
    assert world.num_layers == 12
    assert [not layer.is_liquid for layer in world] == [True, False] + [True] * 10
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
        assert layer.shear_rheology is None
        assert layer.bulk_rheology is None
    world.solve_eos(G_to_use=G, verbose=False)
    assert math.isnan(world.get_shear_viscosity(_MANTLE_RADIUS_M))
    assert math.isnan(world.get_bulk_viscosity(_MANTLE_RADIUS_M))
    assert world.get_melt_fraction(_MANTLE_RADIUS_M) == 0.0
    result = world.solve_love_numbers(frequency=2.0 * math.pi / 86400.0, degree_l=2)
    assert result["success"] is True, result["message"]
    assert abs(world.love_number_k.imag) < 1.0e-6


def test_a_radius_range_refines_every_layer_in_it():
    """One table with a radius range refines every detected layer between those boundaries, numbering their names
    from the inside out."""
    world = build_world(_prem_config(layers={
        "mantle": {
            "radius_range_m": [3480.0e3, 6371.0e3],
            "shear_rheology": {"model": "maxwell"},
            "material": {"solid": {"shear_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e21}}},
        },
    }))
    assert [layer.name for layer in world] == ["layer_0", "layer_1"] + [f"mantle_{i}" for i in range(10)]
    assert [layer.shear_rheology is not None for layer in world] == [False, False] + [True] * 10
    world.solve_eos(G_to_use=G, verbose=False)
    assert world.get_shear_viscosity(_MANTLE_RADIUS_M) == pytest.approx(1.0e21)
    result = world.solve_love_numbers(frequency=2.0 * math.pi / 86400.0, degree_l=2)
    assert result["success"] is True, result["message"]
    assert world.love_number_k.imag < 0.0


def _prem_rows_as_profile(split):
    """PREM's rows as standalone radial_solver arrays: three declared layers with the mantle's discontinuities inside
    one, or a layer declared at every repeated radius (the layers the data-file builder detects)."""
    arrays = _prem_arrays()
    columns = np.column_stack(
        [arrays[key] for key in ("radius_m", "density_kg_m3", "shear_modulus_pa", "bulk_modulus_pa")])
    # Midpoints on the existing segments give every stretch between repeated radii the five slices a layer of the
    # standalone solver needs, without changing the piecewise-linear profile.
    blocks = []
    edges = np.concatenate(([0], np.flatnonzero(np.diff(columns[:, 0]) == 0.0) + 1, [len(columns)]))
    for start, stop in zip(edges[:-1], edges[1:]):
        block = columns[start:stop]
        while len(block) < 8:
            merged = np.empty((2 * len(block) - 1, block.shape[1]))
            merged[0::2], merged[1::2] = block, 0.5 * (block[:-1] + block[1:])
            block = merged
        blocks.append(block)
    radius, density, shear, bulk = np.concatenate(blocks).T
    tops = [1221.5e3, 3480.0e3] + (list(_PREM_MANTLE_BREAKS) if split else []) + [float(radius[-1])]
    num_layers = len(tops)
    layer_types = tuple("liquid" if top == 3480.0e3 else "solid" for top in tops)
    return (np.ascontiguousarray(radius), np.ascontiguousarray(density), np.ascontiguousarray(bulk + 0j),
            np.ascontiguousarray(shear + 0j), 1.4052e-4, 5513.26, layer_types, (True,) * num_layers,
            (False,) * num_layers, np.asarray(tops))


def test_splitting_at_the_discontinuities_keeps_the_profile():
    """Declaring a layer at each discontinuity describes the same piecewise-linear profile, so the converged k2 agree;
    at the default tolerances the split profile is the more accurate one, and the standalone solver warns about the
    jumps left inside a declared layer."""
    from TidalPy.RadialSolver import radial_solver, warn_if_internal_discontinuity

    tight = dict(integration_rtol=1.0e-11, integration_atol=1.0e-15, eos_rtol=1.0e-11, eos_atol=1.0e-15,
                 eos_pressure_tol=1.0e-9, raise_on_fail=True, warnings=False)
    converged = {split: complex(radial_solver(*_prem_rows_as_profile(split), **tight).k) for split in (False, True)}
    assert converged[True] == pytest.approx(converged[False], rel=1.0e-8)
    k2_split = complex(radial_solver(*_prem_rows_as_profile(True), raise_on_fail=True, warnings=False).k)
    assert abs(k2_split - converged[True]) / abs(converged[True]) < 1.0e-5

    merged = _prem_rows_as_profile(False)
    assert warn_if_internal_discontinuity(merged[0], merged[9]) == pytest.approx(list(_PREM_MANTLE_BREAKS))
    split = _prem_rows_as_profile(True)
    assert warn_if_internal_discontinuity(split[0], split[9]) == []
