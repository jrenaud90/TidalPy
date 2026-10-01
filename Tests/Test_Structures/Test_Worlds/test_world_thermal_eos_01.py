"""Temperature profiles and heat flow in the whole-planet EOS solve against closed forms (isothermal, conductive,
adiabatic), plus the thermal EOS and bundled worlds.
"""
import math

import numpy as np
import pytest

from TidalPy.Structures.configs import build_world
from TidalPy.Structures.configs.world_builder import construct_world
from TidalPy.Viscosity import make_viscosity

_RADIUS = 2.0e6           # [m]
_MASS = 4.0e22            # [kg]
_CORE_FRACTION = 0.5
_CORE_DENSITY = 8000.0    # [kg/m^3]
_MANTLE_DENSITY = 3300.0  # [kg/m^3]
_CONDUCTIVITY = 4.0       # [W/(m K)]
_HEAT_CAPACITY = 1200.0   # [J/(kg K)]
_EXPANSION = 3.0e-5       # [1/K]
# Bundled-world silicate and iron thermal constants: conductivity [W m-1 K-1], heat capacity [J kg-1 K-1].
_ROCK_THERMAL = {"thermal_conductivity_w_mk": 3.75, "heat_capacity_j_kgk": 1200.0}
_ROCK_EXPANSION = 5.2e-5  # [1/K]
_IRON_THERMAL = {"thermal_conductivity_w_mk": 7.95, "heat_capacity_j_kgk": 840.0}
_IRON_EXPANSION = 1.2e-5  # [1/K]


def _config(core_temperature=1800.0, mantle_temperature=1600.0, cooling="conduction", **mantle_keys):
    """An isothermal iron core under a silicate mantle with the given cooling model."""
    mantle = {"layer_index": 1, "radius_fraction": 1.0, "temperature_k": mantle_temperature,
              # The expansivity drives the adiabat and convection; it drives density only with
              # use_thermal_expansion.
              "material": {"solid": {
                  "thermal_conductivity_w_mk": _CONDUCTIVITY, "heat_capacity_j_kgk": _HEAT_CAPACITY,
                  "eos": {"model": "constant", "reference_density_kg_m3": _MANTLE_DENSITY,
                          "thermal_expansion_1_k": _EXPANSION},
                  "shear_modulus": {"model": "constant", "shear_modulus_pa": 6.0e10}}},
              "cooling": {"model": cooling}}
    mantle["material"]["solid"].update(mantle_keys)
    return {
        "schema_version": "0.2.0",
        "name": "thermal_test",
        "type": "terrestrial",
        "radius_m": _RADIUS,
        "mass_kg": _MASS,
        "layers": {
            "core": {"layer_index": 0, "radius_fraction": _CORE_FRACTION, "temperature_k": core_temperature,
                     "material": {"solid": {
                         "thermal_conductivity_w_mk": _CONDUCTIVITY, "heat_capacity_j_kgk": _HEAT_CAPACITY,
                         "eos": {"model": "constant", "reference_density_kg_m3": _CORE_DENSITY}}},
                     "cooling": {"model": "off"}},
            "mantle": mantle,
        },
    }


def _solve(config=None, **kwargs):
    world = construct_world(config if config is not None else _config())
    result = world.solve_eos(**kwargs)
    assert result["success"], result["message"]
    return world, result


def _convecting_world():
    config = _config(
        cooling="convection",
        shear_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e18})
    config["layers"]["mantle"]["cooling"] = {"model": "convection"}
    return config


def _isothermal_config_with_expansion(temperature, use_thermal_expansion):
    """Both layers at one temperature with a thermal expansivity, optionally using the thermal EOS."""
    config = _config(core_temperature=temperature, mantle_temperature=temperature)
    for layer in config["layers"].values():
        eos = layer["material"]["solid"]["eos"]
        eos["thermal_expansion_1_k"] = _EXPANSION
        if use_thermal_expansion:
            layer["use_thermal_expansion"] = True
            eos["reference_temperature_k"] = 300.0
    return config


def _melting_rock(layer):
    """Give a silicate layer of a bundled world silicate thermal constants, convection, chondritic radiogenics, and a
    Henning-weakened melt between a 1600 K solidus and a 2000 K liquidus (a Murnaghan melt of 0.2 Pa s)."""
    solid = layer["material"]["solid"]
    solid.update(_ROCK_THERMAL)
    solid["eos"]["thermal_expansion_1_k"] = _ROCK_EXPANSION
    layer["material"]["liquid"] = {
        **_ROCK_THERMAL,
        "eos": {"model": "murnaghan", "reference_density_kg_m3": 2750.0, "reference_bulk_modulus_pa": 2.0e10,
                "bulk_modulus_derivative": 5.0, "thermal_expansion_1_k": _ROCK_EXPANSION},
        "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 0.2}}
    layer["material"]["melting"] = {"solidus": {"model": "constant", "temperature_k": 1600.0},
                                    "liquidus": {"model": "constant", "temperature_k": 2000.0},
                                    "weakening": {"model": "henning"}}
    layer["use_melting"] = True
    layer["cooling"] = {"model": "convection", "convection_alpha": 1.0, "convection_beta": 1.0 / 3.0,
                        "critical_rayleigh": 1100.0}
    layer["radiogenics"] = {"model": "isotope", "isotopes": "modern_day_chondritic"}


def _thermal_io():
    """Bundled Io whose silicate layers melt and convect over an iron core, ready for a thermal solve."""
    config = build_world("io").get_config_dict()
    for name in ("mantle", "asthenosphere"):
        _melting_rock(config["layers"][name])
    core = config["layers"]["core"]
    core["material"]["solid"].update(_IRON_THERMAL)
    core["material"]["solid"]["eos"]["thermal_expansion_1_k"] = _IRON_EXPANSION
    core["cooling"] = {"model": "off"}
    core["radiogenics"] = {"model": "off"}
    return build_world(config)


def test_uniform_world_matches_the_untemperatured_solve():
    """With one temperature everywhere the structure is bit-identical to a solve without temperature."""
    uniform, with_temperature = _solve(_config(core_temperature=1500.0, mantle_temperature=1500.0))
    without, no_temperature = _solve(_config(core_temperature=1500.0, mantle_temperature=1500.0),
                                     solve_temperature=False)
    for key in ("planet_mass", "planet_moi", "surface_gravity", "central_pressure"):
        assert with_temperature[key] == no_temperature[key], key
    assert np.array_equal(with_temperature["density"], no_temperature["density"])
    assert np.array_equal(with_temperature["pressure"], no_temperature["pressure"])
    assert uniform.get_temperature(0.5 * _RADIUS) == pytest.approx(1500.0)
    assert without.get_temperature(0.5 * _RADIUS) == pytest.approx(1500.0)
    assert with_temperature["thermal_passes"] == 0
    assert np.all(with_temperature["heat_flow"] == 0.0)


def test_switch_off_keeps_each_layer_at_its_own_temperature():
    """solve_temperature=False keeps each layer's temperature and reports no flow."""
    world, result = _solve(solve_temperature=False)
    assert world.get_temperature(0.25 * _RADIUS) == pytest.approx(1800.0)
    assert world.get_temperature(0.75 * _RADIUS) == pytest.approx(1600.0)
    assert result["thermal_passes"] == 0
    assert world.get_heat_flow(0.75 * _RADIUS) == 0.0


def test_uniform_override_replaces_every_layer_temperature():
    """A temperature argument overrides every layer's temperature."""
    world, _ = _solve(temperature=900.0)
    assert world.get_temperature(0.25 * _RADIUS) == pytest.approx(900.0)
    assert world.get_temperature(0.75 * _RADIUS) == pytest.approx(900.0)


def test_isothermal_layer_holds_its_temperature():
    """An isothermal layer holds its temperature and carries no heat flow."""
    world, _ = _solve()
    core_radius = _CORE_FRACTION * _RADIUS
    for radius in (0.1 * core_radius, 0.5 * core_radius, 0.99 * core_radius):
        assert world.get_temperature(radius) == pytest.approx(1800.0, rel=1e-10)
    assert world.get_heat_flow(0.5 * core_radius) == pytest.approx(0.0)


def test_conducting_half_follows_fouriers_law():
    """Above the mid-radius, T(r) = T_mid - (L / 4 pi k)(1/r_mid - 1/r), reaching the surface temperature."""
    world, result = _solve(_config(cooling="conduction"), surface_temperature=300.0)
    assert result["thermal_converged"]
    r_inner = _CORE_FRACTION * _RADIUS
    r_mid = 0.5 * (r_inner + _RADIUS)
    heat_flow = world.get_heat_flow(0.5 * (r_mid + _RADIUS))
    assert heat_flow > 0.0
    for radius in np.linspace(r_mid, _RADIUS, 7):
        expected = 1600.0 - (heat_flow / (4.0 * math.pi * _CONDUCTIVITY)) * (1.0 / r_mid - 1.0 / radius)
        assert world.get_temperature(radius) == pytest.approx(expected, rel=1e-8)
    assert world.get_temperature(_RADIUS) == pytest.approx(300.0, rel=1e-6)
    # The layer's own temperature applies at its mid-radius.
    assert world.get_temperature(r_mid) == pytest.approx(1600.0, rel=1e-9)


def test_surface_heat_flow_matches_the_shell_resistance():
    """L = (T_layer - T_surface) / R with R the resistance of the outer conducting half."""
    world, _ = _solve(_config(cooling="conduction"), surface_temperature=300.0)
    r_inner = _CORE_FRACTION * _RADIUS
    r_mid = 0.5 * (r_inner + _RADIUS)
    resistance = (1.0 / r_mid - 1.0 / _RADIUS) / (4.0 * math.pi * _CONDUCTIVITY)
    expected = (1600.0 - 300.0) / resistance
    assert world.get_heat_flow(_RADIUS) == pytest.approx(expected, rel=1e-8)


def test_heat_flow_is_continuous_across_an_interface():
    """The flow leaving one layer enters the next; nothing enters the isothermal core."""
    _, result = _solve(_config(cooling="conduction"), surface_temperature=300.0)
    assert result["layer_heat_flow_out"][0] == pytest.approx(result["layer_heat_flow_in"][1])
    assert result["layer_heat_flow_in"][0] == 0.0


def test_no_surface_temperature_leaves_no_flow():
    """Without a surface temperature no heat leaves the mantle."""
    world, result = _solve(_config(cooling="conduction"))
    assert result["layer_heat_flow_out"][1] == 0.0
    assert world.get_temperature(_RADIUS) == pytest.approx(1600.0, rel=1e-9)


def test_convecting_layer_has_an_adiabatic_interior_between_boundary_layers():
    """A convecting layer has a gently cooling interior and carries most of its drop in the boundary layers. Its
    temperature applies at the top of the interior, and the adiabat warms below it."""
    world, result = _solve(_convecting_world(), surface_temperature=300.0)
    assert result["thermal_converged"]
    r_inner = _CORE_FRACTION * _RADIUS
    boundary = world.mantle.radius_outer - world.mantle.radius_inner
    base_temperature = world.get_temperature(r_inner + 0.25 * boundary)
    top_temperature = world.get_temperature(_RADIUS - 0.25 * boundary)
    assert base_temperature > top_temperature > 300.0
    assert result["layer_top_temperature"][1] == 1600.0
    assert result["layer_base_temperature"][1] > 1600.0
    assert world.get_temperature(_RADIUS) == pytest.approx(300.0, rel=1e-6)
    assert (base_temperature - top_temperature) < 0.25 * (base_temperature - 300.0)


def test_adiabat_follows_its_closed_form():
    """The convecting interior follows T = T_base exp(-alpha / c_p int g dr)."""
    world, _ = _solve(_convecting_world(), surface_temperature=300.0)
    r_inner = world.mantle.radius_inner
    thickness = world.mantle.radius_outer - r_inner
    # Sample well away from both boundary layers.
    radii = np.linspace(r_inner + 0.42 * thickness, world.mantle.radius_outer - 0.42 * thickness, 9)
    temperatures = world.get_temperature(radii)
    gravity = world.get_gravity(radii)
    integral = np.concatenate([[0.0], np.cumsum(
        0.5 * (gravity[1:] + gravity[:-1]) * np.diff(radii) * _EXPANSION / _HEAT_CAPACITY)])
    expected = temperatures[0] * np.exp(-integral)
    assert temperatures == pytest.approx(expected, rel=1e-5)
    assert temperatures[-1] < temperatures[0]


def test_layer_temperature_is_the_top_of_the_adiabat():
    """The solved profile reaches the layer's own temperature at the top of its interior, under the upper boundary
    layer, and the reported base is the adiabat carried down to the interior's base."""
    world, result = _solve(_convecting_world(), surface_temperature=300.0)
    interior_top = world.mantle.radius_outer - result["layer_boundary_thickness"][1]
    interior_base = world.mantle.radius_inner + result["layer_boundary_thickness"][1]
    assert world.get_temperature(interior_top) == pytest.approx(1600.0, rel=1e-6)
    radii = np.linspace(interior_base, interior_top, 2001)
    gravity = world.get_gravity(radii)
    exponent = np.trapezoid(gravity, radii) * _EXPANSION / _HEAT_CAPACITY
    assert result["layer_base_temperature"][1] == pytest.approx(1600.0 * np.exp(exponent), rel=1e-6)
    assert world.get_temperature(interior_base) == pytest.approx(result["layer_base_temperature"][1], rel=1e-6)


@pytest.mark.parametrize("mantle_temperature", [1600.0, 2400.0])
def test_a_molten_convecting_layer_keeps_its_own_temperature(mantle_temperature):
    """A nearly inviscid convecting layer has millimetre boundary layers (resistances near 1e-18 K/W). They still
    conduct: the layer keeps its own temperature, and its heat loss follows its temperature drop."""
    config = _config(mantle_temperature=mantle_temperature, cooling="convection",
                     shear_viscosity={"model": "constant", "reference_viscosity_pas": 0.2})
    world, result = _solve(config, surface_temperature=300.0)
    assert result["layer_boundary_thickness"][1] < 1.0e-2
    interior_top = world.mantle.radius_outer - result["layer_boundary_thickness"][1]
    assert world.get_temperature(interior_top) == pytest.approx(mantle_temperature, rel=1e-6)
    resistance = (1.0 / interior_top - 1.0 / _RADIUS) / (4.0 * math.pi * _CONDUCTIVITY)
    assert result["layer_heat_flow_out"][1] == pytest.approx((mantle_temperature - 300.0) / resistance, rel=1e-6)


def test_convecting_layer_reports_its_rayleigh_and_boundary_layer():
    """A repeated convecting solve converges in at least one thermal pass with outward surface flow."""
    world, _ = _solve(_convecting_world(), surface_temperature=300.0)
    thermal = world.solve_eos(surface_temperature=300.0)
    assert thermal["thermal_passes"] >= 1
    assert thermal["thermal_converged"]
    assert world.get_heat_flow(_RADIUS) > 0.0


def test_layer_temperature_rate_is_the_heat_imbalance():
    """M c_p dT/dt = L_in - L_out, and a net-losing layer cools."""
    world, result = _solve(_config(cooling="conduction"), surface_temperature=300.0)
    mantle_mass = world.mantle.mass
    expected = ((result["layer_heat_flow_in"][1] - result["layer_heat_flow_out"][1])
                / (mantle_mass * _HEAT_CAPACITY))
    assert result["layer_temperature_rate"][1] == pytest.approx(expected, rel=1e-12)
    assert result["layer_temperature_rate"][1] < 0.0


def test_viscosity_follows_the_solved_profile():
    """An Arrhenius layer's viscosity follows the solved temperature and is stiffer where colder."""
    viscosity_config = {
        "reference_viscosity_pas": 1.0e20,
        "reference_temperature_k": 1600.0,
        "molar_activation_energy_j_mol": 3.0e5}
    config = _config(cooling="conduction")
    config["layers"]["mantle"]["material"]["solid"]["shear_viscosity"] = dict(model="reference", **viscosity_config)
    world, result = _solve(config, surface_temperature=300.0)
    viscosity_model = make_viscosity("reference", dict(viscosity_config))
    # Compare at the solve's slices: between slices the stored profile interpolates an exponential linearly.
    radii = result["radius"]
    for index in (len(radii) // 2 + 10, len(radii) - 10):
        radius = radii[index]
        expected = viscosity_model.calc_viscosity(world.get_temperature(radius), world.get_pressure(radius))
        assert world.get_shear_viscosity(radius) == pytest.approx(expected, rel=1e-9)
    r_mid = 0.5 * (_CORE_FRACTION * _RADIUS + _RADIUS)
    assert world.get_shear_viscosity(0.99 * _RADIUS) > world.get_shear_viscosity(r_mid)


def test_thermal_eos_expands_the_hot_interior():
    """With use_thermal_expansion a hot world is lighter than the same world cold, by exp(-alpha dT)."""
    _, cold_result = _solve(_isothermal_config_with_expansion(300.0, use_thermal_expansion=True))
    _, hot_result = _solve(_isothermal_config_with_expansion(2500.0, use_thermal_expansion=True))
    assert hot_result["planet_mass"] < cold_result["planet_mass"]
    assert hot_result["planet_mass"] / cold_result["planet_mass"] == pytest.approx(
        math.exp(-_EXPANSION * 2200.0), rel=5e-3)


def test_thermal_eos_is_off_by_default():
    """Without use_thermal_expansion the mass does not depend on temperature."""
    _, hot_result = _solve(_isothermal_config_with_expansion(2500.0, use_thermal_expansion=False))
    _, cold_result = _solve(_isothermal_config_with_expansion(300.0, use_thermal_expansion=False))
    assert hot_result["planet_mass"] == pytest.approx(cold_result["planet_mass"], rel=1e-12)


@pytest.mark.parametrize("world_name", ["io", "europa", "luna", "mercury", "earth_simple"])
def test_bundled_worlds_are_unchanged(world_name):
    """Isothermal bundled worlds solve identically with and without the temperature solve."""
    thermal = build_world(world_name)
    thermal_result = thermal.solve_eos()
    fast = build_world(world_name)
    fast_result = fast.solve_eos(solve_temperature=False)
    for key in ("planet_mass", "planet_moi", "surface_gravity", "central_pressure"):
        assert thermal_result[key] == fast_result[key], key
    assert thermal_result["thermal_passes"] == 0


@pytest.mark.parametrize("core_temperature", (1700.0, 1800.0, 1850.0, 1900.0))
def test_hot_core_under_an_isothermal_mantle_keeps_planetary_heat_flow(core_temperature):
    """Io (melting and convecting) with a core hotter than its mantle keeps a planetary core heat flow and a physical
    k2."""
    io = _thermal_io()
    io.core.temperature = core_temperature
    io.solve_eos(solve_temperature=True, surface_temperature=110.0)
    heat_flow_into_mantle = float(io.get_heat_flow(np.array([io.mantle.radius_inner + 10.0]))[0])
    assert abs(heat_flow_into_mantle) < 1.0e12   # [W]; Io's whole output is about 1e14 W
    io.solve_love_numbers(frequency=4.11e-5, degree_l=2)
    love_k2 = complex(io.love_number_k)
    assert io.love_success, io.love_message
    # Io's measured k2 is 0.125 +/- 0.047 (Park et al. 2024).
    assert 0.02 < love_k2.real < 0.15
    assert -0.05 < love_k2.imag < 0.0
    # A molten stretch raises the amplification (to about 10 here), far below the warning level.
    assert io.love_surface_amplification < 1.0e4
