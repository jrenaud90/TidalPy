"""A world's EOS and Love solves read each layer's material with the layer's switches: closed forms, thermal expansion,
liquid-only and forced states, melting into molten stretches, the pressure-law range, and invalidation."""

import math
import os

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Cooling.cooling import make_cooling
from TidalPy.Material import Material, Phase, load_material
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


def _pressure_law_phase(liquid=False, **eos):
    """A Birch-Murnaghan phase (K0 = 1e11 Pa, K0' = 4, so a tension limit of -0.185 K0) that a hot layer with
    thermal expansion puts in tension; the liquid has no shear law."""
    eos = {"model": "birch_murnaghan", "reference_density_kg_m3": _DENSITY, "reference_bulk_modulus_pa": 1.0e11,
           "bulk_modulus_derivative": 4.0, "thermal_expansion_1_k": 1.0e-3, "reference_temperature_k": 300.0, **eos}
    if liquid:
        return Phase(eos=eos, shear_viscosity={"model": "constant", "reference_viscosity_pas": 1.0})
    return Phase(eos=eos, shear_modulus={"model": "constant", "shear_modulus_pa": 5.0e10})


def _two_phase_material(solid, liquid):
    return Material(solid=solid, liquid=liquid, solidus={"model": "constant", "temperature_k": 1500.0},
                    liquidus={"model": "constant", "temperature_k": 1600.0})


def test_a_liquid_layer_is_checked_against_its_liquid_s_law():
    """A forced-liquid ice Ih body at 1000 K, with thermal expansion and its melt density, is water: the ice law, in
    tension at its surface (a thermal pressure of 0.98 GPa), describes none of it and does not fail the solve. The
    law it does use is checked."""
    world = _one_layer_world(load_material("ice_ih"), state="liquid", temperature=1000.0,
                             use_thermal_expansion=True, use_melt_density=True, shear_rheology=None)
    result = world.solve_eos()
    assert result["success"], result["message"]
    water = _one_layer_world(load_material("water"), temperature=1000.0, use_thermal_expansion=True,
                             shear_rheology=None)
    assert water.solve_eos()["success"]
    radii = np.linspace(0.0, _RADIUS, 5)
    np.testing.assert_allclose(world.body.get_density(radii), water.body.get_density(radii), rtol=1.0e-10)

    # A liquid whose own law is in tension fails the solve, and the message names the liquid's law.
    hot_liquid = _two_phase_material(_solid(), _pressure_law_phase(liquid=True))
    world = _two_layer_world(hot_liquid, state="liquid", temperature=5000.0, use_thermal_expansion=True,
                             use_melt_density=True)
    result = world.solve_eos()
    assert not result["success"]
    assert "in tension past what its material's pressure law (its liquid phase's) represents" in result["message"]


@pytest.mark.parametrize("use_melt_density", [False, True])
def test_a_liquid_layer_s_density_law_is_the_one_checked(use_melt_density):
    """Without its melt density a liquid layer keeps its solid's density, so the solid's law is checked; with it,
    the solid's law describes nothing there."""
    hot_solid = _two_phase_material(_pressure_law_phase(), _liquid())
    world = _two_layer_world(hot_solid, state="liquid", temperature=5000.0, use_thermal_expansion=True,
                             use_melt_density=use_melt_density)
    result = world.solve_eos()
    assert result["success"] == use_melt_density
    if not use_melt_density:
        assert "pressure law (its solid phase's)" in result["message"]


@pytest.mark.parametrize("temperature, success", [(1550.0, False), (2000.0, True)])
def test_a_melting_layer_is_checked_against_the_phases_it_holds(temperature, success):
    """In a melting layer with its melt density the solid's law describes the partially molten points and none of
    the fully molten ones."""
    hot_solid = _two_phase_material(_pressure_law_phase(), _liquid())
    world = _one_layer_world(hot_solid, temperature=temperature, use_thermal_expansion=True, use_melting=True,
                             use_melt_density=True, shear_rheology=None)
    result = world.solve_eos()
    assert result["success"] == success
    if not success:
        assert "pressure law (its solid phase's)" in result["message"]


@pytest.mark.parametrize("derivative, success", [(float("nan"), True), (-7.4e7, False)])
def test_the_tension_end_follows_a_temperature_dependent_bulk_modulus(derivative, success):
    """A law whose K0 follows the temperature has its tension limit there. At 1800 K the thermal pressure,
    alpha0 K0 (T - T_ref) = 0.15 K0, is inside the law's -0.185 K0 at a fixed K0, but past it once K0 has fallen to
    half (a slope of K0' alpha0 K0 - 7.4e7 = -3.4e7 Pa K-1 over 1500 K): the layer fails only by its own, scaled
    range."""
    material = Material(solid=_pressure_law_phase(thermal_expansion_1_k=1.0e-4,
                                                  bulk_modulus_temperature_derivative_pa_k=derivative))
    world = _one_layer_world(material, temperature=1800.0, use_thermal_expansion=True, shear_rheology=None)
    result = world.solve_eos()
    assert result["success"] == success
    if not success:
        assert "in tension past what its material's pressure law represents" in result["message"]


@pytest.mark.parametrize("temperature, success", [(1100.0, True), (1300.0, False)])
def test_a_layer_past_its_law_s_range_in_temperature_fails(temperature, success):
    """A dK_T/dT of -1e8 Pa K-1 takes K0 = 1e11 Pa to a tenth of itself at 1200 K, where the law holds it. A layer
    hotter than that no longer follows its material, so the solve fails and says why."""
    material = Material(solid=_pressure_law_phase(thermal_expansion_1_k=0.0,
                                                  bulk_modulus_temperature_derivative_pa_k=-1.0e8))
    world = _one_layer_world(material, temperature=temperature, shear_rheology=None)
    result = world.solve_eos()
    assert result["success"] == success
    if not success:
        assert "reaches 1300 K, above the 1200 K past which its material's pressure law holds its bulk modulus" \
            in result["message"]


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


def _forced_liquid_pair(**core_kwargs):
    """A core set liquid with a material that melts but sits below its solidus with melting off, and the same core
    as a liquid-only material of that material's liquid phase. The melt density is mixed in, so both are its liquid.
    The solid phase carries a default Maxwell rheology, which a liquid layer does not take."""
    solid = Phase(config={**_solid().get_config_dict(), "shear_rheology": {"model": "maxwell"}})
    two_phase = Material(solid=solid, liquid=_liquid(8000.0),
                         solidus={"model": "constant", "temperature_k": 1500.0},
                         liquidus={"model": "constant", "temperature_k": 1600.0})
    forced = _two_layer_world(two_phase, state="liquid", temperature=1000.0, use_melt_density=True, **core_kwargs)
    liquid = _two_layer_world(Material(liquid=_liquid(8000.0)), temperature=1000.0, **core_kwargs)
    for world in (forced, liquid):
        assert world.solve_eos()["success"]
    return forced, liquid


def test_a_liquid_layer_reports_its_liquid():
    """A layer set liquid reports what the liquid equations use, the material fully molten: no rigidity, the liquid's
    viscosity and bulk modulus, in every getter, the complex moduli, and the radial solver's inputs."""
    forced, liquid = _forced_liquid_pair()
    core = forced.core
    radii = np.linspace(0.0, core.radius_outer, 5)
    np.testing.assert_array_equal(core.get_shear_modulus(radii), 0.0)
    np.testing.assert_array_equal(forced.get_shear_modulus(radii), 0.0)
    np.testing.assert_array_equal(core.get_shear_viscosity(radii), 1.0)
    np.testing.assert_array_equal(core.get_bulk_modulus(radii), 2.0e10)
    np.testing.assert_array_equal(core.get_melt_fraction(radii), 1.0)
    np.testing.assert_array_equal(core.get_density(radii), liquid.core.get_density(radii))
    np.testing.assert_array_equal(core.calc_complex_shear_modulus(radii, _FREQUENCY), 0.0)
    np.testing.assert_array_equal(forced.calc_complex_shear_modulus(radii, _FREQUENCY), 0.0)
    assert core.calc_complex_shear_modulus(_FREQUENCY) == 0.0
    assert core.shear_rheology is None
    state = core.calc_state(0.0)
    assert state["phase"] == "liquid" and state["shear_modulus"] == 0.0 and state["shear_viscosity"] == 1.0
    # The solid phase stays readable through the material.
    assert core.material.solid.calc_state(0.0, 1000.0)["shear_modulus"] == 5.0e10
    _k2(forced)
    solution = forced.release_radial_solution()
    assert solution.get_shear_modulus(0.25 * _RADIUS) == 0.0
    assert solution.get_complex_shear_modulus(0.25 * _RADIUS) == 0.0
    assert solution.get_shear_viscosity(0.25 * _RADIUS) == 1.0


@pytest.mark.parametrize("is_static", [True, False], ids=["static", "dynamic"])
def test_a_liquid_layer_solves_as_its_liquid(is_static):
    """The Love numbers of a layer set liquid are those of its liquid as a liquid-only material: a dynamic liquid
    reads the liquid's bulk modulus, and the quasi-homogeneous methods take the layer as a fluid."""
    forced, liquid = _forced_liquid_pair(is_static=is_static, use_tides=True)
    assert _k2(forced) == pytest.approx(_k2(liquid), rel=1.0e-12)
    for method, kwargs in (("homogeneous", {}), ("cpl", {"fixed_q": 30.0}), ("ctl", {"fixed_dt": 100.0})):
        loves = []
        for world in (forced, liquid):
            world.solve_love_numbers(frequency=_FREQUENCY, love_method=method, **kwargs)
            loves.append(complex(np.asarray(world.love_number_k).ravel()[0]))
        assert loves[0] == pytest.approx(loves[1], rel=1.0e-12)


def test_a_liquid_layer_without_a_liquid_phase_has_no_rigidity():
    """A layer set liquid whose material has no liquid phase keeps that phase, with a zero shear modulus."""
    world = _two_layer_world(Material(solid=_solid(8000.0)), state="liquid", temperature=1000.0)
    world.solve_eos()
    radius = 0.25 * _RADIUS
    assert world.core.get_shear_modulus(radius) == 0.0
    assert world.core.get_shear_viscosity(radius) == 1.0e20
    assert world.core.get_density(radius) == pytest.approx(8000.0)
    assert world.core.get_melt_fraction(radius) == 1.0


def test_a_liquid_layer_convects_as_its_liquid():
    """The magma-ocean scaling of a layer set liquid reads its liquid's viscosity."""
    nusselt_numbers = []
    for forced in (True, False):
        material = (Material(solid=_solid(), liquid=_liquid(),
                             solidus={"model": "constant", "temperature_k": 3000.0},
                             liquidus={"model": "constant", "temperature_k": 3500.0})
                    if forced else Material(liquid=_liquid()))
        world = _two_layer_world(material, state="liquid" if forced else "auto", temperature=2000.0,
                                 cooling=make_cooling("convection"))
        world.mantle.temperature = 1500.0
        world.mantle.cooling = make_cooling("conduction")
        result = world.solve_eos(surface_temperature=300.0, solve_temperature=True)
        assert result["success"]
        assert result["layer_magma_ocean"][0]
        assert result["layer_reference_viscosity"][0] == pytest.approx(1.0)
        nusselt_numbers.append(result["layer_nusselt_number"][0])
    assert nusselt_numbers[0] == pytest.approx(nusselt_numbers[1], rel=1.0e-12)


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
    lambda layer: setattr(layer, "state", "liquid"),
], ids=["material", "thermal_expansion", "melting", "temperature", "liquid_state"])
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
