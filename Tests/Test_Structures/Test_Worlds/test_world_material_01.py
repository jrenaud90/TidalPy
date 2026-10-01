"""A world's EOS and Love solves read each layer's material with the layer's switches: closed forms, thermal expansion,
liquid-only and forced states, melting into molten stretches, the pressure-law range, and invalidation."""

import math
import os

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Cooling.cooling import make_cooling
from TidalPy.Material import Material, Phase
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds.base import BaseWorld

_RADIUS = 1.0e6
_DENSITY = 3300.0
_MASS = (4.0 / 3.0) * math.pi * _RADIUS ** 3 * _DENSITY
_FREQUENCY = 1.0e-5


def _solid(density=_DENSITY, alpha=0.0, gruneisen=0.0):
    return Phase(
        eos={"model": "constant", "reference_density_kg_m3": density, "bulk_modulus_pa": 1.0e11,
             "thermal_expansion_1_k": alpha, "reference_temperature_k": 300.0, "gruneisen_parameter": gruneisen},
        shear_modulus={"model": "constant", "shear_modulus_pa": 5.0e10},
        shear_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e20},
        thermal_conductivity=3.0)


def _liquid(density=_DENSITY):
    return Phase(
        eos={"model": "constant", "reference_density_kg_m3": density, "bulk_modulus_pa": 2.0e10},
        shear_viscosity={"model": "constant", "reference_viscosity_pas": 1.0},
        thermal_conductivity=3.0)


def _melting_material(solidus=1500.0, liquidus=1600.0):
    return Material(solid=_solid(), liquid=_liquid(),
                    solidus={"model": "constant", "temperature_k": solidus},
                    liquidus={"model": "constant", "temperature_k": liquidus})


def _one_layer_world(material, **layer_kwargs):
    world = BaseWorld("w", _RADIUS, _MASS)
    layer_kwargs.setdefault("shear_rheology", "maxwell")
    world.add_layer(Layer("body", 0, 0.0, _RADIUS, _MASS, material, **layer_kwargs))
    return world


def _two_layer_world(core_material, mantle_material=None, **core_kwargs):
    world = BaseWorld("w", _RADIUS, _MASS)
    world.add_layer(Layer("core", 0, 0.0, 0.5 * _RADIUS, 0.0, core_material, **core_kwargs))
    world.add_layer(Layer("mantle", 1, 0.5 * _RADIUS, _RADIUS, 0.0,
                          mantle_material or Material(solid=_solid()), shear_rheology="maxwell"))
    return world


def _k2(world):
    world.solve_love_numbers(frequency=_FREQUENCY)
    return complex(np.asarray(world.love_number_k).ravel()[0])


# =====================================================================================================================
# Structure
# =====================================================================================================================
def test_uniform_sphere_gravity_is_the_closed_form():
    """The EOS solve starts from the exact center limit of dg/dr, so a uniform sphere's gravity is exact."""
    world = _one_layer_world(Material(solid=_solid(4000.0)))
    world.solve_eos(G_to_use=G, rtol=1.0e-10, atol=1.0e-14)
    for radius in (1.0e3, 0.3 * _RADIUS, _RADIUS):
        assert world.get_gravity(radius) == pytest.approx((4.0 / 3.0) * math.pi * G * 4000.0 * radius, rel=1.0e-9)


def test_every_layer_needs_a_material():
    world = BaseWorld("w", _RADIUS, _MASS)
    world.add_layer(Layer("body", 0, 0.0, _RADIUS, _MASS))
    assert not world.all_materials_set
    with pytest.raises(ValueError, match="must have a material"):
        world.solve_eos()


def test_thermal_expansion_switch():
    material = Material(solid=_solid(alpha=3.0e-5))
    cold = _one_layer_world(material, temperature=1300.0)
    cold.solve_eos()
    assert cold.get_density(0.5 * _RADIUS) == pytest.approx(_DENSITY)
    hot = _one_layer_world(material, temperature=1300.0, use_thermal_expansion=True)
    hot.solve_eos()
    assert hot.get_density(0.5 * _RADIUS) == pytest.approx(_DENSITY * math.exp(-3.0e-5 * 1000.0), rel=1.0e-12)


def test_the_solved_bulk_modulus_is_adiabatic():
    """A tidal deformation is adiabatic: the layer reports K_S = K_T (1 + alpha gamma T)."""
    world = _one_layer_world(Material(solid=_solid(alpha=3.0e-5, gruneisen=1.2)), temperature=1500.0)
    world.solve_eos()
    assert world.get_bulk_modulus(0.5 * _RADIUS) == pytest.approx(1.0e11 * (1.0 + 3.0e-5 * 1.2 * 1500.0))


def test_the_tension_end_of_a_pressure_law_fails_the_solve():
    """Past the tension end of its law a hot layer cannot hold together, so the solve fails and names it."""
    material = Material(solid=Phase(
        eos={"model": "birch_murnaghan", "reference_density_kg_m3": _DENSITY, "reference_bulk_modulus_pa": 1.0e11,
             "thermal_expansion_1_k": 1.0e-3, "reference_temperature_k": 300.0},
        shear_modulus={"model": "constant", "shear_modulus_pa": 5.0e10}))
    world = _one_layer_world(material, temperature=5000.0, use_thermal_expansion=True, shear_rheology=None)
    result = world.solve_eos()
    assert not result["success"]
    assert "in tension past what its material's pressure law represents" in result["message"]


# =====================================================================================================================
# Liquids and melting
# =====================================================================================================================
def test_a_liquid_only_material_is_a_liquid_layer():
    liquid_core = _two_layer_world(Material(liquid=_liquid(8000.0)))
    forced_core = _two_layer_world(Material(solid=_solid(8000.0)), state="liquid")
    solid_core = _two_layer_world(Material(solid=_solid(8000.0)), shear_rheology="maxwell")
    for world in (liquid_core, forced_core, solid_core):
        world.solve_eos()
    assert liquid_core.core.is_liquid and forced_core.core.is_liquid and not solid_core.core.is_liquid
    assert _k2(liquid_core) == pytest.approx(_k2(forced_core), rel=1.0e-10)
    assert abs(_k2(liquid_core) - _k2(solid_core)) > 1.0e-3 * abs(_k2(solid_core))


@pytest.mark.parametrize("use_melting, temperature, molten", [
    (False, 2000.0, False), (True, 1400.0, False), (True, 1550.0, False), (True, 2000.0, True)])
def test_melting_needs_its_switch_and_a_fully_molten_state(use_melting, temperature, molten):
    """Without melt weakening a partially molten material keeps its solid modulus, so only a fully molten one is
    liquid for the radial solver; with melting off the material is its solid phase whatever the temperature."""
    world = _one_layer_world(_melting_material(), temperature=temperature, use_melting=use_melting)
    world.solve_eos()
    assert bool(world.molten_regions) == molten
    assert world.get_melt_fraction(0.5 * _RADIUS) == pytest.approx(
        0.0 if not use_melting else min(max((temperature - 1500.0) / 100.0, 0.0), 1.0))


def test_a_molten_interior_under_a_conducting_lid():
    """A hot conducting layer melts where its profile passes the liquidus: steady conduction from the mid-radius
    (2000 K) to the surface (300 K) is T = -1400 + 1700 R / r, which reaches the 1600 K liquidus at r = 17 R / 30."""
    world = _one_layer_world(_melting_material(), temperature=2000.0, use_melting=True,
                             cooling=make_cooling("conduction"))
    assert world.solve_eos(surface_temperature=300.0, solve_temperature=True)["success"]
    regions = world.molten_regions
    assert len(regions) == 1
    name, radius_inner, radius_outer = regions[0]
    assert name == "body" and radius_inner == 0.0
    assert radius_outer == pytest.approx(17.0 / 30.0 * _RADIUS, rel=1.0e-6)
    assert world.get_shear_modulus(0.3 * _RADIUS) == 0.0
    assert world.get_shear_modulus(0.8 * _RADIUS) == 5.0e10
    assert np.isfinite(_k2(world))


# =====================================================================================================================
# Invalidation and round trips
# =====================================================================================================================
@pytest.mark.parametrize("change", [
    lambda layer: setattr(layer, "material", Material(solid=_solid(3000.0))),
    lambda layer: setattr(layer, "use_thermal_expansion", True),
    lambda layer: setattr(layer, "use_melting", True),
    lambda layer: setattr(layer, "temperature", 1000.0),
], ids=["material", "thermal_expansion", "melting", "temperature"])
def test_a_layer_change_forgets_the_solve(change):
    world = _one_layer_world(Material(solid=_solid()))
    world.solve_eos()
    assert world.eos_solved
    change(world.body)
    assert not world.eos_solved


def test_a_material_swap_takes_effect_at_the_next_solve():
    world = _one_layer_world(Material(solid=_solid()))
    world.solve_eos()
    world.body.material = Material(solid=_solid(3000.0))
    world.solve_eos()
    assert world.get_density(0.5 * _RADIUS) == pytest.approx(3000.0)


def test_world_binary_round_trip_solves_alike(tmp_path):
    world = _two_layer_world(Material(liquid=_liquid(8000.0)), _melting_material())
    world.mantle.use_melting = True
    world.mantle.temperature = 1550.0
    world.solve_eos()
    path = os.path.join(tmp_path, "w.tpyb")
    world.save_binary(path)
    loaded = BaseWorld("other", 1.0, 1.0)
    loaded.load_binary(path)
    loaded.solve_eos()
    radii = np.linspace(0.0, _RADIUS, 9)
    np.testing.assert_array_equal(loaded.get_density(radii), world.get_density(radii))
    np.testing.assert_array_equal(loaded.get_melt_fraction(radii), world.get_melt_fraction(radii))
    assert loaded.core.is_liquid and loaded.mantle.use_melting
    assert _k2(loaded) == _k2(world)
