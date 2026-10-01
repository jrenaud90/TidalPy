"""Solid and liquid zones found by the EOS integration, layers that hold their mass, and the latent heat a moving
melt boundary carries.

The EOS solve ends a piece of its integration where a layer's rigidity margin changes sign, so a layer whose material
melts is split into solid and liquid zones at the root itself, and a layer that holds its mass ends where it encloses
that mass. The cases here have closed forms: a uniform-density planet has P(r) = (2/3) pi G rho^2 (R^2 - r^2), so a
linear melting curve puts the edge of a solid inner core at a known radius, and a conducting shell carrying a constant
heat flow has T(r) = T_base - L (1/r_base - 1/r) / (4 pi k).
"""
import math

import numpy as np
import pytest

import TidalPy
from TidalPy.constants import G
from TidalPy.Material import Material
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Utilities.logging.logger import flush_logger, init_logger
from TidalPy.initialize import build_logging_config

_RADIUS = 2.0e6
_DENSITY = 4000.0
_MASS = 4.0 / 3.0 * math.pi * _DENSITY * _RADIUS**3
# A melting temperature linear in pressure, T_m = T0 (1 + P / a) (Simon-Glatzel with c = 1).
_MELT_T0 = 1500.0
_MELT_A = 2.0e10
_LATENT_HEAT = 4.0e5


def _central_pressure():
    """The center pressure of the uniform planet [Pa]."""
    return (2.0 / 3.0) * math.pi * G * _DENSITY**2 * _RADIUS**2


def _melting_material(conductivity=3.0):
    """Incompressible solid and liquid of one density, melting at one temperature T0 (1 + P / a)."""
    thermal = {"thermal_conductivity_w_mk": conductivity, "heat_capacity_j_kgk": 1000.0}
    eos = {"model": "constant", "reference_density_kg_m3": _DENSITY, "bulk_modulus_pa": 1.0e11}
    curve = {"model": "simon_glatzel", "temperature_k": _MELT_T0, "simon_a_pa": _MELT_A, "simon_c": 1.0}
    return Material(config={
        "solid": {**thermal, "eos": eos, "shear_modulus": {"model": "constant", "shear_modulus_pa": 5.0e10},
                  "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e20}},
        "liquid": {**thermal, "eos": eos, "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0}},
        "melting": {"solidus": curve, "liquidus": curve},
        "latent_heat_j_kg": _LATENT_HEAT})


def _uniform_planet(temperature):
    """One isothermal layer that melts at T0 (1 + P / a): solid below the radius where the melting curve meets its
    temperature, liquid above."""
    layer = Layer("body", 0, 0.0, _RADIUS, material=_melting_material(), temperature=temperature, use_melting=True,
                  use_pressure_melting=True)
    world = BaseWorld("uniform", _RADIUS, _MASS)
    world.add_layer(layer)
    result = world.solve_eos(G_to_use=G)
    assert result["success"], result["message"]
    return world, result


def _core_radius(temperature):
    """Where T0 (1 + P / a) = T in the uniform planet [m]."""
    boundary_pressure = _MELT_A * (temperature / _MELT_T0 - 1.0)
    return math.sqrt(_RADIUS**2 - boundary_pressure / ((2.0 / 3.0) * math.pi * G * _DENSITY**2))


@pytest.mark.parametrize("temperature", [1800.0, 2000.0, 2100.0])
def test_a_solid_core_ends_where_the_melting_curve_meets_the_temperature(temperature):
    world, result = _uniform_planet(temperature)
    zones = result["zones"]
    assert [zone["state"] for zone in zones] == ["solid", "liquid"]
    assert zones[0]["radius_inner"] == 0.0
    assert zones[-1]["radius_outer"] == _RADIUS
    assert zones[0]["radius_outer"] == pytest.approx(_core_radius(temperature), rel=1.0e-9)
    assert zones[1]["radius_inner"] == zones[0]["radius_outer"]
    # The enclosed mass at the boundary is the uniform sphere's.
    assert zones[0]["mass_outer"] == pytest.approx(
        4.0 / 3.0 * math.pi * _DENSITY * zones[0]["radius_outer"]**3, rel=1.0e-9)
    assert zones[1]["mass_outer"] == pytest.approx(_MASS, rel=1.0e-9)
    # The property reports the same zones, and the liquid zone is a molten region of the layer.
    assert world.zones == zones
    assert world.molten_regions == [("body", zones[1]["radius_inner"], _RADIUS)]


def test_the_core_grows_as_the_layer_cools():
    # The center melts at about 2170 K, so the whole planet is liquid above that.
    core_radii = [_uniform_planet(temperature)[1]["zones"][0]["radius_outer"]
                  for temperature in (2100.0, 2000.0, 1800.0)]
    assert core_radii[0] < core_radii[1] < core_radii[2]


def test_the_love_solve_takes_each_zone_as_a_layer():
    """The liquid zone over the solid core is solved with the liquid equations: its k2 lies between those of a body
    solid throughout and one liquid throughout, and it is the k2 of the same two layers built by hand."""
    world, result = _uniform_planet(2000.0)
    core_top = result["zones"][0]["radius_outer"]
    world.solve_love_numbers(frequency=1.0e-5)
    k2_zoned = world.love_number_k

    material = _melting_material()
    split = BaseWorld("split", _RADIUS, _MASS)
    split.add_layer(Layer("core", 0, 0.0, core_top, material=material, temperature=2000.0, state="solid"))
    split.add_layer(Layer("ocean", 1, core_top, _RADIUS, material=material, temperature=2000.0, state="liquid"))
    assert split.solve_eos(G_to_use=G)["success"]
    split.solve_love_numbers(frequency=1.0e-5)
    assert k2_zoned.real == pytest.approx(split.love_number_k.real, rel=1.0e-6)


def test_the_latent_heat_of_the_moving_boundary():
    """C_latent = L dM_liquid/dT: for an isothermal layer the boundary moves by 1 / |dT_m/dr| = a / (T0 rho g) per
    kelvin, so C_latent = L 4 pi r_b^2 rho a / (T0 rho g_b), with g_b = (4/3) pi G rho r_b. A centered difference of the
    liquid mass between two solves agrees."""
    temperature = 2000.0
    world, result = _uniform_planet(temperature)
    core_top = result["zones"][0]["radius_outer"]
    boundary_gravity = (4.0 / 3.0) * math.pi * G * _DENSITY * core_top
    expected = _LATENT_HEAT * 4.0 * math.pi * core_top**2 * _MELT_A / (_MELT_T0 * boundary_gravity)
    assert result["layer_latent_capacity"][0] == pytest.approx(expected, rel=1.0e-5)

    delta = 0.5
    liquid_mass = []
    for shifted in (temperature - delta, temperature + delta):
        zones = _uniform_planet(shifted)[1]["zones"]
        liquid_mass.append(zones[1]["mass_outer"] - zones[1]["mass_inner"])
    finite_difference = _LATENT_HEAT * (liquid_mass[1] - liquid_mass[0]) / (2.0 * delta)
    assert result["layer_latent_capacity"][0] == pytest.approx(finite_difference, rel=1.0e-5)


@pytest.fixture
def spdlog_text(tmp_path):
    """Route the C++ logger to a temporary file for the test and hand back a reader for its text."""
    log_path = tmp_path / "tidalpy.log"
    init_logger({"console_level": "off", "file_level": "warning", "log_to_file": True,
                 "log_file_path": str(log_path)})

    def read():
        flush_logger()
        return log_path.read_text(encoding="utf-8") if log_path.exists() else ""

    yield read
    init_logger(build_logging_config())


@pytest.mark.parametrize("is_static", [True, False])
def test_a_dynamic_liquid_zone_is_checked_for_its_instability(spdlog_text, is_static):
    """A liquid zone of a layer solved with the dynamic equations is a dynamic liquid like a liquid layer, and a
    constant-density one grows unstable at long periods: the Love solve warns about it, naming the layer."""
    world, _ = _uniform_planet(2000.0)
    world.body.is_static = is_static
    try:
        world.solve_love_numbers(frequency=1.0e-6)
    except Exception:  # noqa: BLE001
        pass
    text = spdlog_text()
    assert ("grows unstable at long periods" in text) == (not is_static)
    if not is_static:
        assert "'body'" in text


def test_a_state_change_that_lets_a_layer_melt_needs_a_new_solve():
    """The solve looks for zones only in a layer that can change state, so setting a forced layer back to auto (or
    the reverse) forgets the solve; a change between forced states does not."""
    layer = Layer("body", 0, 0.0, _RADIUS, material=_melting_material(), temperature=2000.0, use_melting=True,
                  use_pressure_melting=True, state="solid")
    world = BaseWorld("uniform", _RADIUS, _MASS)
    world.add_layer(layer)
    assert world.solve_eos(G_to_use=G)["success"]
    assert [zone["state"] for zone in world.zones] == ["solid"]
    world.body.state = "liquid"
    assert world.eos_solved
    assert world.molten_regions == []
    world.body.state = "auto"
    assert not world.eos_solved
    assert world.zones == []
    assert [zone["state"] for zone in world.solve_eos(G_to_use=G)["zones"]] == ["solid", "liquid"]


def test_a_layer_that_cannot_change_state_is_one_zone():
    layer = Layer("body", 0, 0.0, _RADIUS, material=_melting_material(), temperature=2000.0, use_melting=False)
    world = BaseWorld("uniform", _RADIUS, _MASS)
    world.add_layer(layer)
    result = world.solve_eos(G_to_use=G)
    assert [zone["state"] for zone in result["zones"]] == ["solid"]
    assert result["layer_latent_capacity"] == [0.0]
    # A layer forced liquid after the solve is one liquid zone to the Love solve, without a new EOS solve.
    world.body.state = "liquid"
    assert world.molten_regions == []
    assert world.zones[0]["state"] == "solid"


def test_a_thin_zone_takes_its_neighbors_state():
    """A liquid skin thinner than minimum_zone_fraction of the radius joins the solid below it."""
    numerical = TidalPy.config["numerical"]
    default_fraction = numerical["minimum_zone_fraction"]
    # The liquid shell at 1800 K is about a quarter of the radius thick.
    try:
        numerical["minimum_zone_fraction"] = 0.5
        TidalPy.constants.update_constants()
        world, result = _uniform_planet(1800.0)
        assert [zone["state"] for zone in result["zones"]] == ["solid"]
        assert world.molten_regions == []
    finally:
        numerical["minimum_zone_fraction"] = default_fraction
        TidalPy.constants.update_constants()
    world, result = _uniform_planet(1800.0)
    assert [zone["state"] for zone in result["zones"]] == ["solid", "liquid"]


# ======================================================================================================================
# A conducting shell
# ======================================================================================================================
_CORE_RADIUS = 1.0e6
_CORE_TEMPERATURE = 2500.0
_SOLIDUS = 1800.0
_CONDUCTIVITY = 4.0


def _shell_world():
    """A hot isothermal core under a conducting shell whose material melts at a constant 1800 K: the base of the shell
    is liquid, the rest solid."""
    curve = {"model": "constant", "temperature_k": _SOLIDUS}
    thermal = {"thermal_conductivity_w_mk": _CONDUCTIVITY, "heat_capacity_j_kgk": 1000.0}
    eos = {"model": "constant", "reference_density_kg_m3": 3300.0, "bulk_modulus_pa": 1.0e11}
    shell_material = Material(config={
        "solid": {**thermal, "eos": eos, "shear_modulus": {"model": "constant", "shear_modulus_pa": 5.0e10},
                  "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e20}},
        "liquid": {**thermal, "eos": eos, "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0}},
        "melting": {"solidus": curve, "liquidus": curve}})
    core = Layer("core", 0, 0.0, _CORE_RADIUS, material="simple_iron_core", temperature=_CORE_TEMPERATURE)
    shell = Layer("shell", 1, _CORE_RADIUS, _RADIUS, material=shell_material, temperature=1500.0, use_melting=True,
                  cooling="conduction")
    mass = 4.0 / 3.0 * math.pi * (8000.0 * _CORE_RADIUS**3 + 3300.0 * (_RADIUS**3 - _CORE_RADIUS**3))
    world = BaseWorld("shelled", _RADIUS, mass)
    world.add_layer(core)
    world.add_layer(shell)
    return world


def test_a_conducting_shell_melts_where_its_profile_crosses_the_solidus():
    """The base of the shell carries the heat flow L entering it, so T(r) = T_core - L (1/r_core - 1/r) / (4 pi k)
    there, and the liquid zone ends where that falls to the solidus."""
    world = _shell_world()
    result = world.solve_eos(G_to_use=G, solve_temperature=True, surface_temperature=300.0)
    assert result["success"], result["message"]
    assert result["thermal_converged"]
    zones = [zone for zone in result["zones"] if zone["layer"] == "shell"]
    assert [zone["state"] for zone in zones] == ["liquid", "solid"]
    heat_flow = result["layer_heat_flow_in"][1]
    expected = 1.0 / (1.0 / _CORE_RADIUS
                      - (_CORE_TEMPERATURE - _SOLIDUS) * 4.0 * math.pi * _CONDUCTIVITY / heat_flow)
    assert zones[0]["radius_outer"] == pytest.approx(expected, rel=1.0e-8)
    assert world.get_temperature(zones[0]["radius_outer"]) == pytest.approx(_SOLIDUS, rel=1.0e-9)


# ======================================================================================================================
# Layers that hold their mass
# ======================================================================================================================
def _expanding_world(mantle_holds_mass, crust=True):
    """An iron core under a mantle of fixed mass that expands with temperature (alpha 3e-5), under a thin crust that
    holds its volume."""
    mantle_material = Material(config={"solid": {
        "eos": {"model": "constant", "reference_density_kg_m3": 3300.0, "bulk_modulus_pa": 1.0e11,
                "thermal_expansion_1_k": 3.0e-5, "reference_temperature_k": 300.0},
        "shear_modulus": {"model": "constant", "shear_modulus_pa": 5.0e10}}})
    mantle_top = _RADIUS - 5.0e4 if crust else _RADIUS
    mantle_mass = 4.0 / 3.0 * math.pi * 3300.0 * (mantle_top**3 - _CORE_RADIUS**3)
    core = Layer("core", 0, 0.0, _CORE_RADIUS, material="simple_iron_core", temperature=300.0)
    mantle = Layer("mantle", 1, _CORE_RADIUS, mantle_top, mass=mantle_mass, material=mantle_material,
                   temperature=300.0, use_thermal_expansion=True, is_volume_fixed=not mantle_holds_mass)
    world = BaseWorld("expanding", _RADIUS, 6.0e22)
    world.add_layer(core)
    world.add_layer(mantle)
    if crust:
        world.add_layer(Layer("crust", 2, mantle_top, _RADIUS, material="simple_rock", temperature=300.0))
    return world, mantle_mass


def test_a_layer_holding_its_mass_expands_with_its_temperature():
    """At 300 K the mantle fills its given radii; heated by dT its density falls by exp(-alpha dT), so it holds the
    same mass in a volume larger by exp(alpha dT): its top moves out, about alpha dT / 3 in radius for a thin shell,
    and exactly where the solved structure encloses its mass. The crust above keeps its volume."""
    world, mantle_mass = _expanding_world(mantle_holds_mass=True)
    cold = world.solve_eos(G_to_use=G)
    assert cold["success"], cold["message"]
    assert world.mantle.mass == pytest.approx(mantle_mass, rel=1.0e-9)
    assert world.radius == pytest.approx(_RADIUS, rel=1.0e-9)

    heating = 500.0
    hot = world.solve_eos(G_to_use=G, temperature=300.0 + heating)
    assert hot["success"], hot["message"]
    mantle_top = world.mantle.radius_outer
    expected_top = (_CORE_RADIUS**3 + (_RADIUS - 5.0e4)**3 * math.exp(3.0e-5 * heating)
                    - _CORE_RADIUS**3 * math.exp(3.0e-5 * heating))**(1.0 / 3.0)
    assert mantle_top == pytest.approx(expected_top, rel=1.0e-9)
    assert world.mantle.mass == pytest.approx(mantle_mass, rel=1.0e-9)
    # The crust keeps its volume over the moved base, and the world ends where it does.
    crust_volume = (_RADIUS**3 - (_RADIUS - 5.0e4)**3)
    assert world.crust.radius_outer**3 - world.crust.radius_inner**3 == pytest.approx(crust_volume, rel=1.0e-9)
    assert world.crust.radius_inner == mantle_top
    assert world.radius == world.crust.radius_outer
    assert hot["layer_radius_outer"][-1] == world.radius
    # Its zones follow the moved radii.
    assert [zone["radius_outer"] for zone in hot["zones"]] == [layer.radius_outer for layer in world]


def test_a_convecting_layer_holding_its_mass_in_a_thermal_solve():
    """A convecting mantle holding its mass under a thermal profile: its segments (boundary layers around an adiabatic
    interior) end at fractions of its mass, so they move with it, and the solve converges with the mass held and the
    upper boundary layer ending at the mantle's top."""
    world, mantle_mass = _expanding_world(mantle_holds_mass=True, crust=False)
    world.mantle.temperature = 1600.0
    world.mantle.cooling = "convection"
    result = world.solve_eos(G_to_use=G, solve_temperature=True, surface_temperature=300.0)
    assert result["success"], result["message"]
    assert result["thermal_converged"]
    assert result["thermal_passes"] > 1
    assert world.mantle.mass == pytest.approx(mantle_mass, rel=1.0e-9)
    # The adiabatic interior expanded the mantle past its cold radii.
    assert world.radius > _RADIUS
    boundary = result["layer_boundary_thickness"][1]
    assert 0.0 < boundary < 0.5 * (world.mantle.radius_outer - world.mantle.radius_inner)
    # The layer's own temperature applies at the top of its interior, under the upper boundary layer.
    assert world.get_temperature(world.mantle.radius_outer - boundary) == pytest.approx(1600.0, rel=1.0e-6)
    assert world.get_temperature(world.radius) == pytest.approx(300.0, rel=1.0e-6)


def test_reads_above_a_layer_holding_its_mass_are_nan():
    """The integration toward a layer's mass ran past the top it found; nothing above that top is read from it."""
    world, _ = _expanding_world(mantle_holds_mass=True, crust=False)
    assert world.solve_eos(G_to_use=G, temperature=1300.0)["success"]
    assert math.isnan(world.get_pressure(1.1 * world.radius))
    assert math.isnan(world.mantle.get_pressure(1.1 * world.radius))
    assert math.isfinite(world.get_pressure(world.radius))

    world, _ = _expanding_world(mantle_holds_mass=True)
    assert world.solve_eos(G_to_use=G, temperature=1300.0)["success"]
    assert math.isnan(world.mantle.get_pressure(world.radius))
    assert world.crust.get_pressure(world.radius) == pytest.approx(0.0, abs=1.0e-6 * world.central_pressure)


def test_a_layer_holding_its_mass_that_melts():
    """A melting layer that holds its mass is split into zones and still ends where it encloses its mass."""
    core = Layer("core", 0, 0.0, 0.5 * _RADIUS, material="simple_iron_core", temperature=2000.0)
    outer_mass = 4.0 / 3.0 * math.pi * _DENSITY * (_RADIUS**3 - (0.5 * _RADIUS)**3)
    body = Layer("body", 1, 0.5 * _RADIUS, _RADIUS, mass=outer_mass, material=_melting_material(), temperature=2000.0,
                 use_melting=True, use_pressure_melting=True, is_volume_fixed=False)
    world = BaseWorld("melting", _RADIUS, _MASS)
    world.add_layer(core)
    world.add_layer(body)
    result = world.solve_eos(G_to_use=G)
    assert result["success"], result["message"]
    assert [zone["state"] for zone in result["zones"] if zone["layer"] == "body"] == ["solid", "liquid"]
    assert world.body.mass == pytest.approx(outer_mass, rel=1.0e-9)
    assert result["zones"][-1]["radius_outer"] == world.radius


def test_a_layer_holding_its_volume_keeps_its_radii():
    world, _ = _expanding_world(mantle_holds_mass=False)
    result = world.solve_eos(G_to_use=G, temperature=800.0)
    assert result["success"], result["message"]
    assert world.mantle.radius_outer == _RADIUS - 5.0e4
    assert world.radius == _RADIUS


def test_a_layer_holding_its_mass_on_top_sets_the_world_radius():
    world, mantle_mass = _expanding_world(mantle_holds_mass=True, crust=False)
    result = world.solve_eos(G_to_use=G, temperature=1300.0)
    assert result["success"], result["message"]
    assert world.radius > _RADIUS
    assert world.mantle.radius_outer == world.radius
    assert world.mantle.mass == pytest.approx(mantle_mass, rel=1.0e-9)
    assert np.isclose(result["surface_pressure"], 0.0, atol=1.0e-6 * result["central_pressure"])
