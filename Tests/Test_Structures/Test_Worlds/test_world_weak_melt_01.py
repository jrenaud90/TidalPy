"""Partial melt stays a melt-weakened solid, and the tides change continuously as a layer melts.

With the default `[numerical] minimum_solid_rigidity` of 0 a layer that can change state is a liquid zone only where its
post-melt shear modulus vanishes. A band of partial melt keeps its own complex modulus, which the radial solver
integrates down to near-zero rigidity (with `minimum_complex_rigidity` 1e-15), so the Love numbers approach the
static-liquid answer continuously as the band's rigidity falls to zero.
"""
import math

import numpy as np
import pytest

import TidalPy
from TidalPy.configurations import find_invalid_config_values, get_packaged_config
from TidalPy.constants import update_constants
from TidalPy.Structures.configs import build_world


RADIUS = 2.0e6                # [m]
FORCING_FREQUENCY = 2.0e-5    # [rad s-1]
CORE_TOP = 0.4 * RADIUS
MIDDLE_TOP = 0.7 * RADIUS
DENSITIES = {"core": 7000.0, "middle": 3500.0, "shell": 3000.0}   # [kg m-3]
ROCK_THERMAL = {"thermal_conductivity_w_mk": 3.75, "heat_capacity_j_kgk": 1200.0}
SOLIDUS = 1400.0    # [K]
LIQUIDUS = 1600.0   # [K]
# Henning weakening with a breakdown band from melt fraction 0.3 to 0.5: the shear modulus falls continuously to zero
# across it, so the layer is a melt-weakened solid through the band and liquid past it.
BAND_MELT = {"solidus": {"model": "constant", "temperature_k": SOLIDUS},
             "liquidus": {"model": "constant", "temperature_k": LIQUIDUS},
             "weakening": {"model": "henning", "crit_melt_frac": 0.3, "crit_melt_frac_width": 0.2,
                           "hn_shear_param_1_k": 0.0}}


def material(density, shear_modulus, viscosity, melting=None):
    """A constant-density solid; with melting, also a melt phase of the same density and a 0.2 Pa s viscosity."""
    eos = {"model": "constant", "reference_density_kg_m3": density, "bulk_modulus_pa": 1.0e11,
           "thermal_expansion_1_k": 5.2e-5}
    out = {"solid": {**ROCK_THERMAL, "eos": eos,
                     "shear_modulus": {"model": "constant", "shear_modulus_pa": shear_modulus},
                     "shear_viscosity": {"model": "constant", "reference_viscosity_pas": viscosity}}}
    if melting is not None:
        out["liquid"] = {**ROCK_THERMAL, "eos": dict(eos),
                         "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 0.2}}
        out["melting"] = melting
    return out


def layer(index, radius_outer, temperature, layer_material, **flags):
    out = {"layer_index": index, "radius_outer_m": radius_outer, "use_tides": True, "temperature_k": temperature,
           "material": layer_material, "use_melting": "melting" in layer_material,
           "shear_rheology": {"model": "maxwell"}, "bulk_rheology": {"model": "elastic"},
           "cooling": {"model": "off"}, "radiogenics": {"model": "off"}}
    out.update(flags)
    return out


def world_with_middle(middle):
    """A solid core and shell around a middle layer, all constant density, its EOS solved."""
    layers = {"core": layer(0, CORE_TOP, 1000.0, material(DENSITIES["core"], 1.0e11, 1.0e22)),
              "middle": middle,
              "shell": layer(2, RADIUS, 1000.0, material(DENSITIES["shell"], 5.0e10, 1.0e22))}
    mass = (4.0 / 3.0) * math.pi * (
        DENSITIES["core"] * CORE_TOP ** 3
        + DENSITIES["middle"] * (MIDDLE_TOP ** 3 - CORE_TOP ** 3)
        + DENSITIES["shell"] * (RADIUS ** 3 - MIDDLE_TOP ** 3))
    world = build_world({"schema_version": "0.2.0", "name": "weak-melt", "type": "terrestrial",
                         "radius_m": RADIUS, "mass_kg": mass, "layers": layers})
    world.solve_eos()
    return world


def melting_world(temperature):
    """The middle layer, incompressible (neutrally stratified), melts uniformly at `temperature` [K]."""
    return world_with_middle(layer(1, MIDDLE_TOP, temperature,
                                   material(DENSITIES["middle"], 6.0e10, 1.0e19, BAND_MELT), is_incompressible=True))


def love_k(world):
    world.solve_love_numbers(frequency=FORCING_FREQUENCY, degree_l=2)
    assert world.love_success, world.love_message
    return complex(world.love_number_k)


def test_partial_melt_stays_solid_and_full_melt_is_liquid():
    """Inside the breakdown band the layer is one solid zone with a weakened modulus; past it, a liquid zone."""
    partial = melting_world(SOLIDUS + 0.4 * (LIQUIDUS - SOLIDUS))
    assert partial.molten_regions == []
    assert 0.0 < partial.middle.get_shear_modulus(0.5 * (CORE_TOP + MIDDLE_TOP)) < 6.0e10
    molten = melting_world(SOLIDUS + 0.8 * (LIQUIDUS - SOLIDUS))
    assert molten.molten_regions == [("middle", CORE_TOP, MIDDLE_TOP)]


def test_the_love_number_is_continuous_through_the_breakdown():
    """k2 changes smoothly as the band's rigidity falls to zero and meets the liquid layer's value without a step."""
    temperatures = np.linspace(SOLIDUS + 0.35 * (LIQUIDUS - SOLIDUS), SOLIDUS + 0.55 * (LIQUIDUS - SOLIDUS), 41)
    worlds = [melting_world(temperature) for temperature in temperatures]
    k = np.array([love_k(world) for world in worlds])
    liquid = np.array([bool(world.molten_regions) for world in worlds])
    assert (not liquid[0]) and liquid[-1]
    onset = int(np.argmax(liquid))
    steps = np.abs(np.diff(k))
    # The step into the liquid zone is no larger than the steps of the weak solid before it.
    assert steps[onset - 1] <= 2.0 * np.max(steps[:onset - 1]) + 1.0e-6 * abs(k[onset])
    assert np.all(np.isfinite(k))


def test_the_weak_solid_approaches_the_static_liquid():
    """A layer of vanishing rigidity solves to the static liquid's Love number."""
    liquid = world_with_middle(layer(1, MIDDLE_TOP, 1000.0, material(DENSITIES["middle"], 6.0e10, 1.0e19),
                                     state="liquid", is_static=True, is_incompressible=True))
    k_liquid = love_k(liquid)
    errors = []
    for shear in (1.0e3, 1.0e0, 1.0e-3):
        weak = world_with_middle(layer(1, MIDDLE_TOP, 1000.0, material(DENSITIES["middle"], shear, 1.0e19),
                                       is_incompressible=True))
        errors.append(abs(love_k(weak) / k_liquid - 1.0))
    # The weak solid's departure falls with its rigidity until it reaches the integration tolerance.
    assert errors[0] > 100.0 * max(errors[1:])
    assert max(errors[1:]) < 1.0e-8


def test_the_continuation_frequency_does_not_read_the_liquid_threshold():
    """calc_continuation_frequency uses its own near-fluid threshold, so a positive minimum_solid_rigidity leaves it
    unchanged."""
    numerical = TidalPy.config["numerical"]
    default = numerical["minimum_solid_rigidity"]
    world = build_world("charon")
    world["hydrosphere"].temperature = 260.0   # A warm shell raises the continuation frequency
    frequencies = []
    try:
        for value in (default, 1.0e-6):
            numerical["minimum_solid_rigidity"] = value
            update_constants()
            world.solve_eos()
            frequencies.append(world.calc_continuation_frequency())
    finally:
        numerical["minimum_solid_rigidity"] = default
        update_constants()
    assert frequencies[0] > 1.0e-16
    assert frequencies[0] == frequencies[1]


@pytest.mark.parametrize("value, problems", [(0.0, 0), (1.0e-6, 0), (-1.0, 1)])
def test_minimum_solid_rigidity_accepts_zero(value, problems):
    packaged = get_packaged_config()
    assert packaged["numerical"]["minimum_solid_rigidity"] == 0.0
    assert packaged["numerical"]["minimum_complex_rigidity"] == 1.0e-15
    assert len(find_invalid_config_values({"numerical": {"minimum_solid_rigidity": value}}, packaged)) == problems
