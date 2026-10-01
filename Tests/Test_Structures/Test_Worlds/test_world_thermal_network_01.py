"""Which layers the thermal network couples, the per-layer thermal detail it reports, and floating-layer grids."""
import math

import pytest

from TidalPy.Structures.configs.world_builder import construct_world

_RADIUS = 2.0e6           # [m]
_CORE_FRACTION = 0.5
_MANTLE_TEMPERATURE = 1600.0
_SURFACE_TEMPERATURE = 200.0


def _material(density):
    return {"solid": {"thermal_conductivity_w_mk": 4.0, "heat_capacity_j_kgk": 1200.0,
                      "eos": {"model": "constant", "reference_density_kg_m3": density, "thermal_expansion_1_k": 3.0e-5},
                      "shear_modulus": {"model": "constant", "shear_modulus_pa": 6.0e10}}}


def _config(core, core_density=8000.0, mantle_cooling="conduction", **mantle_keys):
    core = dict({"layer_index": 0, "radius_fraction": _CORE_FRACTION, "material": _material(core_density)}, **core)
    mantle = dict({"layer_index": 1, "radius_fraction": 1.0,
                   "temperature_k": _MANTLE_TEMPERATURE, "material": _material(3300.0),
                   "cooling": {"model": mantle_cooling}}, **mantle_keys)
    volume = (4.0 / 3.0) * math.pi
    mass = volume * (8000.0 * (_CORE_FRACTION * _RADIUS) ** 3
                     + 3300.0 * (_RADIUS ** 3 - (_CORE_FRACTION * _RADIUS) ** 3))
    return {"schema_version": "0.2.0", "name": "network", "type": "terrestrial", "radius_m": _RADIUS,
            "mass_kg": mass, "layers": {"core": core, "mantle": mantle}}


def _solve(config):
    world = construct_world(config)
    result = world.solve_eos(surface_temperature=_SURFACE_TEMPERATURE)
    assert result["success"], result["message"]
    return world, result


@pytest.mark.parametrize("core", [
    pytest.param({}, id="geometry_only"),
    pytest.param({"temperature_k": 0.0, "cooling": {"model": "off"}}, id="zero_kelvin"),
])
def test_a_layer_without_a_temperature_is_no_heat_sink(core):
    """A layer without a temperature is outside the network and exchanges no heat with the mantle."""
    world, result = _solve(_config(core))
    assert result["layer_in_thermal_network"] == [False, True]
    assert result["layer_heat_flow_in"][1] == 0.0
    assert result["layer_heat_flow_out"][0] == 0.0
    r_core = world.core.radius_outer
    r_mid = 0.5 * (r_core + world.radius)
    for radius in (r_core * (1.0 + 1.0e-6), 0.5 * (r_core + r_mid)):
        assert world.get_temperature(radius) == pytest.approx(_MANTLE_TEMPERATURE, rel=1.0e-9)
    assert result["layer_heat_flow_out"][1] > 0.0


def test_a_warm_core_still_couples():
    """A core with a temperature is in the network and feeds heat into the mantle."""
    _, result = _solve(_config({"temperature_k": 1800.0, "cooling": {"model": "off"}}))
    assert result["layer_in_thermal_network"] == [True, True]
    assert result["layer_heat_flow_in"][1] > 0.0


def test_the_solve_reports_the_convecting_detail():
    """The solve reports per-layer node, top, and base temperatures, boundary layers, Rayleigh and Nusselt numbers.
    A convecting layer's temperature is the top of its interior, and the adiabat warms below it."""
    _, result = _solve(_config({"temperature_k": 1800.0, "cooling": {"model": "off"}},
                               mantle_cooling="convection"))
    for key in ("layer_node_temperature", "layer_top_temperature", "layer_base_temperature", "layer_boundary_thickness",
                "layer_rayleigh_number", "layer_nusselt_number"):
        assert len(result[key]) == 2, key
    assert 0.0 < result["layer_boundary_thickness"][1] <= 0.4 * (1.0 - _CORE_FRACTION) * _RADIUS
    assert result["layer_node_temperature"][1] == _SURFACE_TEMPERATURE
    assert result["layer_top_temperature"][1] == _MANTLE_TEMPERATURE
    assert result["layer_base_temperature"][1] > _MANTLE_TEMPERATURE


def test_floating_layers_end_on_the_solved_grid():
    """Floating layer radii end on the grid of the last solve pass."""
    # A compressible core with less mass than its starting geometry, so its radius takes several passes.
    core = {"is_volume_fixed": False, "mass_kg": 2.0e22,
            "material": {"solid": {"eos": {"model": "birch_murnaghan", "reference_density_kg_m3": 8000.0,
                                           "reference_bulk_modulus_pa": 1.3e11, "bulk_modulus_derivative": 4.5}}}}
    world, result = _solve(_config(core))
    assert result["geometry_converged"]
    assert result["thermal_passes"] > 1
    slices = len(result["radius"]) // 2
    assert result["radius"][slices - 1] == world.core.radius_outer
    assert result["radius"][-1] == world.radius
    assert world.mantle.radius_inner == world.core.radius_outer
