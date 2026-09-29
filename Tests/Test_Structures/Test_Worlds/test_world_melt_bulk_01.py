"""Melt in a solved world: the density the structure integrates, and the bulk response melt sets.

The mixture density enters the structure iteration and the dense readout through one function, so the mass the
solve integrated and the density it reports agree. The melt's compaction bulk viscosity reaches the Love number only
through a bulk rheology that is not elastic.
"""
import math

import numpy as np
import pytest

from TidalPy.PartialMelt import make_partial_melt
from TidalPy.Structures.configs.world_builder import construct_world

_SOLID_DENSITY = 3300.0
_FREQUENCY = 2.0 * math.pi / (1.77 * 86400.0)


def _world(temperature_k=1800.0, solve_temperature=False, density_melt_mixing=True, bulk_rheology=None,
           bulk_viscosity_melt_weakening=False):
    """A two-layer rocky world: an iron core under a mantle whose melt model mixes its melt into the density."""
    mantle = {
        "class": "solidliquid",
        "type": "mantle_rock",
        "radius_fraction": 1.0,
        "temperature_k": temperature_k,
        "material": {
            "model": "constant",
            "reference_density_kg_m3": _SOLID_DENSITY,
            "partial_melt": {"density_melt_mixing": density_melt_mixing,
                             "bulk_viscosity_melt_weakening": bulk_viscosity_melt_weakening},
        },
    }
    if bulk_rheology is not None:
        mantle["bulk_rheology"] = bulk_rheology
    config = {
        "schema_version": "0.2.0",
        "name": "melt_world",
        "type": "terrestrial",
        "radius_m": 1.8e6,
        "mass_kg": 8.9e22,
        "eos_solver": {"solve_temperature": solve_temperature},
        "layers": {
            "core": {"class": "solidliquid", "type": "iron", "radius_fraction": 0.5, "temperature_k": 1800.0},
            "mantle": mantle,
        },
    }
    world = construct_world(config)
    result = world.solve_eos()
    assert result["success"], result["message"]
    return world


def _mantle_radii(world, num=41):
    layer = world.mantle
    span = layer.radius_outer - layer.radius_inner
    return np.linspace(layer.radius_inner + 1.0e-3 * span, layer.radius_outer - 1.0e-3 * span, num)


def _integrated_mass(layer, num=2001):
    """The mass under a layer's reported density, read through the layer so an interface takes its own side."""
    radius = np.linspace(layer.radius_inner, layer.radius_outer, num)
    density = np.array([layer.get_density(r) for r in radius])
    return 4.0 * math.pi * np.trapezoid(radius ** 2 * density, radius)


@pytest.mark.parametrize("solve_temperature", [False, True])
def test_reported_density_is_the_mixture_the_structure_integrated(solve_temperature):
    world = _world(solve_temperature=solve_temperature)
    melt = make_partial_melt("henning", {"density_melt_mixing": True})
    molten = 0
    for radius in _mantle_radii(world):
        state = world.get_state(radius)
        phi = state["melt_fraction"]
        molten += phi > 0.0
        expected = (1.0 - phi) * _SOLID_DENSITY + phi * melt.calc_liquid_density(state["pressure"])
        assert state["density"] == pytest.approx(expected, rel=1.0e-12)
    assert molten > 0
    # The structure iteration used the same density the readout reports, so its mass is the readout's integral.
    for layer in world.layers:
        assert _integrated_mass(layer) == pytest.approx(layer.mass, rel=1.0e-6), layer.name


def test_mixing_is_off_by_default_and_changes_nothing_without_melt():
    melted = _world(density_melt_mixing=True)
    unmixed = _world(density_melt_mixing=False)
    assert unmixed.get_density(unmixed.mantle.radius_outer - 1.0) == pytest.approx(_SOLID_DENSITY)
    # Silicate melt is lighter at these pressures, so the mixed mantle holds less mass.
    assert melted.planet_mass_eos < unmixed.planet_mass_eos
    # Below the solidus the switch has nothing to mix.
    cold_mixed = _world(temperature_k=1500.0, density_melt_mixing=True)
    cold_unmixed = _world(temperature_k=1500.0, density_melt_mixing=False)
    assert cold_mixed.planet_mass_eos == cold_unmixed.planet_mass_eos
    assert make_partial_melt("henning").density_melt_mixing is False


def test_melt_bulk_viscosity_reaches_the_love_number_through_a_zener_bulk_rheology():
    """Melt compaction dissipates only when the bulk rheology lets the bulk modulus relax."""
    elastic = _world(bulk_viscosity_melt_weakening=True)
    zener_rheology = {"model": "zener", "relaxed_modulus_frac": 0.9}
    zener = _world(bulk_rheology=zener_rheology, bulk_viscosity_melt_weakening=True)
    zener_no_melt_viscosity = _world(bulk_rheology=zener_rheology, bulk_viscosity_melt_weakening=False)
    k2 = {}
    for name, world in (("elastic", elastic), ("zener", zener), ("zener_no_melt", zener_no_melt_viscosity)):
        result = world.solve_love_numbers(frequency=_FREQUENCY, degree_l=2)
        assert result["success"], result["message"]
        k2[name] = world.love_number_k
    # The mantle's bulk viscosity is 1e22 Pa s before melt: far too stiff to relax at this period, so a Zener
    # bulk rheology alone barely matters. Melt drops it to about c eta / phi, and the bulk response dissipates.
    assert abs(k2["zener_no_melt"] - k2["elastic"]) < 1.0e-6 * abs(k2["elastic"])
    assert -k2["zener"].imag > -k2["elastic"].imag
    assert k2["zener"].real != pytest.approx(k2["elastic"].real, rel=1.0e-9)
