"""Temperature and heat flow in the whole-planet EOS solve.

Each layer carries its own temperature and its cooling model says how heat moves inside it, so the solve returns a
profile, the heat flow across every interface, and the rate each layer's temperature changes at. The profiles have
closed forms, which is what these tests check: an isothermal layer holds its temperature, a conducting shell follows
Fourier's law, and an adiabat follows its own exponential. Interface flows must agree on both sides, and a world
with no temperature contrast must return exactly what it returned before any of this existed.

Requires the Cython extensions to be compiled first::

    uv pip install -v <repo_root>
"""
import math

import numpy as np
import pytest

from TidalPy.Material_x.eos import make_material_eos
from TidalPy.cooling_x import make_cooling
from TidalPy.structures_x.configs.world_builder import construct_world
from TidalPy.viscosity_x import make_viscosity

_RADIUS = 2.0e6           # [m]
_MASS = 4.0e22            # [kg]
_CORE_FRACTION = 0.5
_CORE_DENSITY = 8000.0    # [kg/m^3]
_MANTLE_DENSITY = 3300.0  # [kg/m^3]
_CONDUCTIVITY = 4.0       # [W/(m K)]
_HEAT_CAPACITY = 1200.0   # [J/(kg K)]
_EXPANSION = 3.0e-5       # [1/K]


def _config(core_temperature=1800.0, mantle_temperature=1600.0, cooling="conduction", **mantle_keys):
    """A two-layer world: an isothermal iron core under a silicate mantle with the given cooling model."""
    mantle = {"class": "solidliquid", "type": "none", "layer_index": 1, "radius_fraction": 1.0,
              "temperature_k": mantle_temperature,
              # The thermal constants are the material's. Its one expansivity drives the adiabat and convection,
              # and the density law too, but only for a layer that sets use_thermal_eos.
              "material": {"model": "constant", "reference_density_kg_m3": _MANTLE_DENSITY,
                           "shear_modulus_static_pa": 6.0e10, "thermal_conductivity_w_mk": _CONDUCTIVITY,
                           "heat_capacity_j_kgk": _HEAT_CAPACITY, "thermal_expansion_1_k": _EXPANSION},
              "cooling": {"model": cooling}}
    # Extra keys are the material's (a viscosity table, say), so they go into its table.
    mantle["material"].update(mantle_keys)
    return {
        "schema_version": "0.2.0",
        "name": "thermal_test",
        "type": "terrestrial",
        "radius_m": _RADIUS,
        "mass_kg": _MASS,
        "layers": {
            "core": {"class": "solidliquid", "type": "none", "layer_index": 0, "radius_fraction": _CORE_FRACTION,
                     "temperature_k": core_temperature,
                     "material": {"model": "constant", "reference_density_kg_m3": _CORE_DENSITY,
                                  "thermal_conductivity_w_mk": _CONDUCTIVITY, "heat_capacity_j_kgk": _HEAT_CAPACITY},
                     "cooling": {"model": "off"}},
            "mantle": mantle,
        },
    }


def _solve(config=None, **kwargs):
    world = construct_world(config if config is not None else _config())
    result = world.solve_eos(**kwargs)
    assert result["success"], result["message"]
    return world, result


# =====================================================================================================================
# No temperature contrast: the solve is what it was before temperature existed
# =====================================================================================================================
def test_uniform_world_matches_the_untemperatured_solve():
    """One temperature everywhere leaves no profile to integrate, so the structure must be bit-identical."""
    uniform, with_temperature = _solve(_config(core_temperature=1500.0, mantle_temperature=1500.0))
    without, no_temperature = _solve(_config(core_temperature=1500.0, mantle_temperature=1500.0),
                                     solve_temperature=False)
    for key in ("planet_mass", "planet_moi", "surface_gravity", "central_pressure"):
        assert with_temperature[key] == no_temperature[key], key
    assert np.array_equal(with_temperature["density"], no_temperature["density"])
    assert np.array_equal(with_temperature["pressure"], no_temperature["pressure"])
    # Both still report the temperature each layer carries.
    assert uniform.get_temperature(0.5 * _RADIUS) == pytest.approx(1500.0)
    assert without.get_temperature(0.5 * _RADIUS) == pytest.approx(1500.0)
    assert with_temperature["thermal_passes"] == 0
    assert np.all(with_temperature["heat_flow"] == 0.0)


def test_switch_off_keeps_each_layer_at_its_own_temperature():
    world, result = _solve(solve_temperature=False)
    assert world.get_temperature(0.25 * _RADIUS) == pytest.approx(1800.0)
    assert world.get_temperature(0.75 * _RADIUS) == pytest.approx(1600.0)
    assert result["thermal_passes"] == 0
    assert world.get_heat_flow(0.75 * _RADIUS) == 0.0


def test_uniform_override_replaces_every_layer_temperature():
    world, _ = _solve(temperature=900.0)
    assert world.get_temperature(0.25 * _RADIUS) == pytest.approx(900.0)
    assert world.get_temperature(0.75 * _RADIUS) == pytest.approx(900.0)


# =====================================================================================================================
# Conduction
# =====================================================================================================================
def test_isothermal_layer_holds_its_temperature():
    world, _ = _solve()
    core_radius = _CORE_FRACTION * _RADIUS
    for radius in (0.1 * core_radius, 0.5 * core_radius, 0.99 * core_radius):
        assert world.get_temperature(radius) == pytest.approx(1800.0, rel=1e-10)
    # A perfectly conducting layer carries no modeled gradient, so nothing flows through it.
    assert world.get_heat_flow(0.5 * core_radius) == pytest.approx(0.0)


def test_conducting_half_follows_fouriers_law():
    """Between the mid-radius and the surface, T(r) = T_mid - (L / 4 pi k)(1/r_mid - 1/r)."""
    world, result = _solve(_config(cooling="conduction"), surface_temperature=300.0)
    assert result["thermal_converged"]
    r_inner = _CORE_FRACTION * _RADIUS
    r_mid = 0.5 * (r_inner + _RADIUS)
    heat_flow = world.get_heat_flow(0.5 * (r_mid + _RADIUS))
    assert heat_flow > 0.0   # heat leaves a hot interior
    for radius in np.linspace(r_mid, _RADIUS, 7):
        expected = 1600.0 - (heat_flow / (4.0 * math.pi * _CONDUCTIVITY)) * (1.0 / r_mid - 1.0 / radius)
        assert world.get_temperature(radius) == pytest.approx(expected, rel=1e-8)
    # The profile arrives at the surface temperature it was given.
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
    """The flow leaving one layer is the flow entering the next."""
    _, result = _solve(_config(cooling="conduction"), surface_temperature=300.0)
    assert result["layer_heat_flow_out"][0] == pytest.approx(result["layer_heat_flow_in"][1])
    # The core is isothermal, so nothing crosses the core-mantle boundary.
    assert result["layer_heat_flow_in"][0] == 0.0


def test_no_surface_temperature_leaves_no_flow():
    world, result = _solve(_config(cooling="conduction"))
    assert result["layer_heat_flow_out"][1] == 0.0
    assert world.get_temperature(_RADIUS) == pytest.approx(1600.0, rel=1e-9)


# =====================================================================================================================
# Convection
# =====================================================================================================================
def _convecting_world():
    config = _config(
        cooling="convection",
        shear_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e18})
    config["layers"]["mantle"]["cooling"] = {"model": "convection"}
    return config


def test_convecting_layer_has_an_adiabatic_interior_between_boundary_layers():
    world, result = _solve(_convecting_world(), surface_temperature=300.0)
    assert result["thermal_converged"]
    r_inner = _CORE_FRACTION * _RADIUS
    boundary = world.mantle.radius_outer - world.mantle.radius_inner
    # The interior is adiabatic: T(r) = T_base exp(-int alpha g / c_p dr), which is nearly linear over a thin
    # shell, and always a far smaller drop than the boundary layers carry.
    base_temperature = world.get_temperature(r_inner + 0.25 * boundary)
    top_temperature = world.get_temperature(_RADIUS - 0.25 * boundary)
    assert base_temperature > top_temperature > 300.0
    assert base_temperature < 1600.0 + 1.0
    # The surface boundary layer carries the bulk of the drop.
    assert world.get_temperature(_RADIUS) == pytest.approx(300.0, rel=1e-6)
    assert (base_temperature - top_temperature) < 0.25 * (base_temperature - 300.0)


def test_adiabat_follows_its_closed_form():
    """Across the interior, dT/dr = -alpha g T / c_p integrates to T = T_base exp(-alpha/c_p int g dr)."""
    world, _ = _solve(_convecting_world(), surface_temperature=300.0)
    r_inner = world.mantle.radius_inner
    thickness = world.mantle.radius_outer - r_inner
    # Sample well inside the interior, away from both boundary layers.
    radii = np.linspace(r_inner + 0.42 * thickness, world.mantle.radius_outer - 0.42 * thickness, 9)
    temperatures = world.get_temperature(radii)
    gravity = world.get_gravity(radii)
    # Trapezoidal integral of alpha g / c_p from the first sample outward.
    integral = np.concatenate([[0.0], np.cumsum(
        0.5 * (gravity[1:] + gravity[:-1]) * np.diff(radii) * _EXPANSION / _HEAT_CAPACITY)])
    expected = temperatures[0] * np.exp(-integral)
    assert temperatures == pytest.approx(expected, rel=1e-5)
    assert temperatures[-1] < temperatures[0]


def test_convecting_layer_reports_its_rayleigh_and_boundary_layer():
    world, _ = _solve(_convecting_world(), surface_temperature=300.0)
    thermal = world.solve_eos(surface_temperature=300.0)
    assert thermal["thermal_passes"] >= 1
    assert thermal["thermal_converged"]
    assert world.get_heat_flow(_RADIUS) > 0.0


# =====================================================================================================================
# Secular terms and the viscosity the profile implies
# =====================================================================================================================
def test_layer_temperature_rate_is_the_heat_imbalance():
    """M c_p dT/dt = L_in - L_out."""
    world, result = _solve(_config(cooling="conduction"), surface_temperature=300.0)
    mantle_mass = world.mantle.mass
    expected = ((result["layer_heat_flow_in"][1] - result["layer_heat_flow_out"][1])
                / (mantle_mass * _HEAT_CAPACITY))
    assert result["layer_temperature_rate"][1] == pytest.approx(expected, rel=1e-12)
    # A layer losing more heat than it receives cools.
    assert result["layer_temperature_rate"][1] < 0.0


def test_viscosity_follows_the_solved_profile():
    """An Arrhenius layer is stiffer where the profile is colder."""
    config = _config(cooling="conduction")
    config["layers"]["mantle"]["material"]["shear_viscosity"] = {
        "model": "reference",
        "reference_viscosity_pas": 1.0e20,
        "reference_temperature_k": 1600.0,
        "molar_activation_energy_j_mol": 3.0e5}
    world, result = _solve(config, surface_temperature=300.0)
    viscosity_model = make_viscosity("reference", {
        "reference_viscosity_pas": 1.0e20,
        "reference_temperature_k": 1600.0,
        "molar_activation_energy_j_mol": 3.0e5})
    # The viscoelastic profile is sampled at the solve's slices, so compare there: away from a slice the stored
    # profile is a linear interpolation of an exponential and only agrees to a few parts in a thousand.
    radii = result["radius"]
    for index in (len(radii) // 2 + 10, len(radii) - 10):
        radius = radii[index]
        expected = viscosity_model.calc_viscosity(world.get_temperature(radius), world.get_pressure(radius))
        assert world.get_shear_viscosity(radius) == pytest.approx(expected, rel=1e-9)
    r_mid = 0.5 * (_CORE_FRACTION * _RADIUS + _RADIUS)
    assert world.get_shear_viscosity(0.99 * _RADIUS) > world.get_shear_viscosity(r_mid)


def test_thermal_eos_expands_the_hot_interior():
    """With use_thermal_eos the density follows the profile, so the same world is lighter when it is hot."""
    cold = _config(core_temperature=300.0, mantle_temperature=300.0)
    hot = _config(core_temperature=2500.0, mantle_temperature=2500.0)
    for config in (cold, hot):
        for layer in config["layers"].values():
            layer["use_thermal_eos"] = True
            layer["material"]["thermal_expansion_1_k"] = _EXPANSION
            layer["material"]["reference_temperature_k"] = 300.0
    _, cold_result = _solve(cold)
    _, hot_result = _solve(hot)
    assert hot_result["planet_mass"] < cold_result["planet_mass"]
    # Thermal expansion of alpha dT: about 6.6 percent less dense over 2200 K.
    assert hot_result["planet_mass"] / cold_result["planet_mass"] == pytest.approx(
        math.exp(-_EXPANSION * 2200.0), rel=5e-3)


def test_thermal_eos_is_off_by_default():
    hot = _config(core_temperature=2500.0, mantle_temperature=2500.0)
    for layer in hot.values() if False else hot["layers"].values():
        layer["material"]["thermal_expansion_1_k"] = _EXPANSION
    _, hot_result = _solve(hot)
    cold = _config(core_temperature=300.0, mantle_temperature=300.0)
    for layer in cold["layers"].values():
        layer["material"]["thermal_expansion_1_k"] = _EXPANSION
    _, cold_result = _solve(cold)
    assert hot_result["planet_mass"] == pytest.approx(cold_result["planet_mass"], rel=1e-12)


# =====================================================================================================================
# Bundled worlds
# =====================================================================================================================
@pytest.mark.parametrize("world_name", ["io", "europa", "luna", "mercury", "earth_simple"])
def test_bundled_worlds_are_unchanged(world_name):
    """Every bundled world is at one temperature, so its solve must be identical either way."""
    from TidalPy.structures_x.configs import build_world
    thermal = build_world(world_name)
    thermal_result = thermal.solve_eos()
    fast = build_world(world_name)
    fast_result = fast.solve_eos(solve_temperature=False)
    for key in ("planet_mass", "planet_moi", "surface_gravity", "central_pressure"):
        assert thermal_result[key] == fast_result[key], key
    assert thermal_result["thermal_passes"] == 0


@pytest.mark.parametrize("core_temperature", (1700.0, 1800.0, 1850.0, 1900.0))
def test_hot_core_under_an_isothermal_mantle_keeps_planetary_heat_flow(core_temperature):
    """Io with a core hotter than a mantle at the asthenosphere's temperature.

    The mantle has no contrast to the layer above, so its convection model once returned a fixed 1 m boundary layer.
    The core then drove about 1e16 W through that metre, melted a sliver at the mantle base, and the Love solve
    returned k2 near -0.58 + 0.25i while reporting success. The boundary layers now take the drop across both of
    them, including the one to the core, so the heat flow stays planetary. A core hot enough to melt the mantle base
    past the liquid limit of its partial-melt model leaves a molten stretch there, which the radial solver treats as
    a static liquid, so every case solves with a physical k2.
    """
    from TidalPy.structures_x.configs import build_world
    io = build_world("io")
    io.core.temperature = core_temperature
    io.solve_eos(solve_temperature=True, surface_temperature=110.0)
    heat_flow_into_mantle = float(io.get_heat_flow(np.array([io.mantle.radius_inner + 10.0]))[0])
    assert abs(heat_flow_into_mantle) < 1.0e12   # [W]; Io's whole output is about 1e14 W
    io.solve_love_numbers(frequency=4.11e-5, degree_l=2)
    love_k2 = complex(io.love_number_k)
    assert io.love_success, io.love_message
    assert 0.02 < love_k2.real < 0.1
    assert -0.05 < love_k2.imag < 0.0
    # A molten stretch is a liquid layer, which raises the surface-solve amplification from about 10 to a few
    # hundred; the solver warns only when the amplification times machine epsilon reaches the tolerance.
    assert io.love_surface_amplification < 1.0e4
