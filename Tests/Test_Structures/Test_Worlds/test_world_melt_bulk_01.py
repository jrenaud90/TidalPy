"""Melt in a solved world: the density the structure integrates, and the bulk response melt sets.

The mixture density enters the structure iteration and the dense readout through one function, so the mass the
solve integrated and the density it reports agree. The melt's compaction bulk viscosity reaches the Love number only
through a bulk rheology that is not elastic.
"""
import math

import numpy as np
import pytest

from TidalPy.Cooling import make_cooling
from TidalPy.Material import Material, Phase
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds.terrestrial import TerrestrialWorld

_RADIUS = 1.8e6                 # [m]
_CORE_RADIUS = 0.5 * _RADIUS    # [m]
_MASS = 8.9e22                  # [kg]
_SOLID_DENSITY = 3300.0         # [kg m-3]
_MELT_DENSITY = 2900.0          # [kg m-3]
_SOLIDUS = 1600.0               # [K]
_LIQUIDUS = 2400.0              # [K]
_FREQUENCY = 2.0 * math.pi / (1.77 * 86400.0)


def _core_material():
    return Material(solid=Phase(
        eos={"model": "constant", "reference_density_kg_m3": 8000.0, "bulk_modulus_pa": 1.5e11},
        shear_modulus={"model": "constant", "shear_modulus_pa": 8.0e10}))


def _mantle_material(compaction):
    """A rock that melts between the solidus and liquidus with Henning weakening; with ``compaction`` its bulk
    viscosity follows the melt's compaction viscosity, otherwise it stays the solid's until fully molten."""
    solid = Phase(
        eos={"model": "constant", "reference_density_kg_m3": _SOLID_DENSITY, "bulk_modulus_pa": 1.0e11},
        shear_modulus={"model": "constant", "shear_modulus_pa": 6.0e10},
        shear_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e17},
        bulk_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e22})
    melt = Phase(
        eos={"model": "constant", "reference_density_kg_m3": _MELT_DENSITY, "bulk_modulus_pa": 2.0e10},
        shear_viscosity={"model": "constant", "reference_viscosity_pas": 1.0})
    return Material(
        solid=solid,
        liquid=melt,
        solidus={"model": "constant", "temperature_k": _SOLIDUS},
        liquidus={"model": "constant", "temperature_k": _LIQUIDUS},
        weakening="henning",
        bulk_viscosity_mixing="compaction" if compaction else None)


def _world(temperature=1800.0, solve_temperature=False, use_melt_density=True, bulk_rheology="elastic",
           compaction=False):
    """A two-layer rocky world: an iron core under a partially molten mantle."""
    world = TerrestrialWorld("melt_world", _RADIUS, _MASS)
    world.add_layer(Layer("core", 0, 0.0, _CORE_RADIUS, 0.0, _core_material(), temperature=1800.0))
    world.add_layer(Layer(
        "mantle",
        1,
        _CORE_RADIUS,
        _RADIUS,
        0.0,
        _mantle_material(compaction),
        temperature=temperature,
        use_melting=True,
        use_melt_density=use_melt_density,
        shear_rheology="maxwell",
        bulk_rheology=bulk_rheology,
        cooling=make_cooling("conduction") if solve_temperature else None))
    # A thermal solve conducts heat from the mantle out to a 300 K surface, so only the deep mantle melts.
    result = world.solve_eos(
        solve_temperature=solve_temperature,
        surface_temperature=300.0 if solve_temperature else None)
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
    material = world.mantle.material
    molten = 0
    for radius in _mantle_radii(world):
        state = world.get_state(radius)
        phi = state["melt_fraction"]
        molten += phi > 0.0
        expected = (1.0 - phi) * material.solid.calc_state(state["pressure"])["density"] \
            + phi * material.liquid.calc_state(state["pressure"])["density"]
        assert state["density"] == pytest.approx(expected, rel=1.0e-12)
    assert molten > 0
    # The structure iteration used the same density the readout reports, so its mass is the readout's integral.
    for layer in world.layers:
        assert _integrated_mass(layer) == pytest.approx(layer.mass, rel=1.0e-6), layer.name


def test_mixing_is_off_by_default_and_changes_nothing_without_melt():
    melted = _world(use_melt_density=True)
    unmixed = _world(use_melt_density=False)
    assert unmixed.get_density(unmixed.mantle.radius_outer - 1.0) == pytest.approx(_SOLID_DENSITY)
    # The melt is lighter, so the mixed mantle holds less mass.
    assert melted.planet_mass_eos < unmixed.planet_mass_eos
    # Below the solidus the switch has nothing to mix.
    cold_mixed = _world(temperature=1500.0, use_melt_density=True)
    cold_unmixed = _world(temperature=1500.0, use_melt_density=False)
    assert cold_mixed.planet_mass_eos == cold_unmixed.planet_mass_eos
    assert Layer("mantle", 0, 0.0, _RADIUS, 0.0).use_melt_density is False


def test_melt_bulk_viscosity_reaches_the_love_number_through_a_zener_bulk_rheology():
    """Melt compaction dissipates only when the bulk rheology lets the bulk modulus relax."""
    elastic = _world(compaction=True)
    zener_rheology = {"model": "zener", "relaxed_modulus_frac": 0.9}
    zener = _world(bulk_rheology=zener_rheology, compaction=True)
    zener_no_melt_viscosity = _world(bulk_rheology=zener_rheology, compaction=False)
    k2 = {}
    for name, world in (("elastic", elastic), ("zener", zener), ("zener_no_melt", zener_no_melt_viscosity)):
        result = world.solve_love_numbers(frequency=_FREQUENCY, degree_l=2)
        assert result["success"], result["message"]
        k2[name] = world.love_number_k
    # The mantle's bulk viscosity is 1e22 Pa s before melt: far too stiff to relax at this period, so a Zener
    # bulk rheology alone barely matters. Compaction drops it to about c eta / phi, and the bulk response dissipates.
    assert abs(k2["zener_no_melt"] - k2["elastic"]) < 1.0e-6 * abs(k2["elastic"])
    assert -k2["zener"].imag > -k2["elastic"].imag
    assert k2["zener"].real != pytest.approx(k2["elastic"].real, rel=1.0e-9)
