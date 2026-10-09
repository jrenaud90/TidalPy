"""A layer's static liquid zone solves like an equivalent static liquid layer.

The bundled Pluto's hydrosphere is one ice Ih layer, liquid below its melting point and solid above it. The radial
solver takes each zone as a layer, and the static-liquid interface conditions (Saito 1974) read the liquid's density at
the zone boundary. That boundary is a root of the rigidity margin, where the material read at the pressure and
temperature lies on whichever side of the melting curve the root finder stopped, so each zone must read its own side.
The same planet built as a water ocean under a solid ice shell is the reference.
"""
import copy
import math

import numpy as np
import pytest

from TidalPy.Structures.configs import build_world


DAY = 86400.0
PERIODS_DAYS = (1.0, 6.3872304, 100.0, 1000.0)


def _solved(source):
    world = build_world(source)
    world.solve_eos(raise_on_fail=True)
    return world


def _love_k(world, periods_days):
    frequencies = 2.0 * math.pi / (np.asarray(periods_days, dtype=float) * DAY)
    return world.calc_love_numbers(
        frequencies,
        degree_l=2,
        warnings=False)["k"]


@pytest.fixture(scope="module")
def pluto_and_split():
    """The bundled Pluto, and the same planet with its hydrosphere split at its solved zone boundary [m]."""
    pluto = _solved("pluto")
    boundary = next(zone for zone in pluto.zones if zone["state"] == "liquid")["radius_outer"]
    config = build_world("pluto").get_config_dict()
    hydrosphere = config["layers"].pop("hydrosphere")
    # A liquid-only water ocean (the melt of ice Ih) up to the boundary, under the same ice made solid.
    config["layers"]["ocean"] = {
        "layer_index": 1,
        "radius_outer_m": boundary,
        "material": "water",
        "temperature_k": hydrosphere["temperature_k"],
        "is_static": True}
    shell = copy.deepcopy(hydrosphere)
    shell["layer_index"] = 2
    shell["state"] = "solid"
    for key in ("use_melting", "use_pressure_melting", "use_melt_density"):
        shell.pop(key, None)
    config["layers"]["shell"] = shell
    split = _solved(config)
    return pluto, split, boundary


def test_split_world_has_the_same_profile(pluto_and_split):
    pluto, split, boundary = pluto_and_split
    assert [zone["state"] for zone in pluto.zones] == [zone["state"] for zone in split.zones]
    radii = np.linspace(pluto.hydrosphere.radius_inner + 1.0, pluto.radius - 1.0, 41)
    np.testing.assert_allclose(split.get_density(radii), pluto.get_density(radii), rtol=1.0e-7)
    np.testing.assert_allclose(split.get_gravity(radii), pluto.get_gravity(radii), rtol=1.0e-7)


def test_zone_boundary_reads_the_lower_zone(pluto_and_split):
    """A profile read on the boundary belongs to the lower zone, the ocean, and just above it to the ice."""
    pluto, _, boundary = pluto_and_split
    assert pluto.get_density(boundary) > 1000.0
    assert pluto.get_density(boundary + 1.0) < 950.0


def test_static_liquid_zone_matches_a_static_liquid_layer(pluto_and_split):
    """With the liquid zone's top read on the solid side, k2 was off by 3e-3 to 6e-3 at every period."""
    pluto, split, _ = pluto_and_split
    k_zone = _love_k(pluto, PERIODS_DAYS)
    k_layer = _love_k(split, PERIODS_DAYS)
    np.testing.assert_array_less(np.abs(k_zone - k_layer) / np.abs(k_layer), 1.0e-7)


def test_static_and_dynamic_liquid_zones_converge_at_long_periods(pluto_and_split):
    """pluto_dynamic's ocean is neutrally stratified, so its dynamic solve tends to the static one."""
    pluto, _, _ = pluto_and_split
    pluto_dynamic = _solved("pluto_dynamic")
    k_static = _love_k(pluto, [1000.0])
    k_dynamic = _love_k(pluto_dynamic, [1000.0])
    assert abs(k_dynamic[0] - k_static[0]) / abs(k_static[0]) < 1.0e-6
