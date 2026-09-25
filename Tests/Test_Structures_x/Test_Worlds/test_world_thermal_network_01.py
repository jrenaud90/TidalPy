"""Which layers the thermal network couples, and what the solve reports about it.

A layer with no temperature of its own (a geometry-only layer, or one left at the 0 K default) takes no part in
the network: its interfaces carry no heat, so it neither drains nor feeds its neighbors. Floating layers end on
the grid the last pass solved on, and the per-layer thermal detail (node and top temperatures, boundary layers,
Rayleigh and Nusselt numbers) comes back with the solve.

Requires the Cython extensions to be compiled first::

    uv pip install -v <repo_root>
"""
import math

import pytest

from TidalPy.structures_x.configs.world_builder import construct_world

_RADIUS = 2.0e6           # [m]
_CORE_FRACTION = 0.5
_MANTLE_TEMPERATURE = 1600.0
_SURFACE_TEMPERATURE = 200.0


def _material(density):
    return {"model": "constant", "reference_density_kg_m3": density, "shear_modulus_static_pa": 6.0e10,
            "thermal_conductivity_w_mk": 4.0, "heat_capacity_j_kgk": 1200.0, "thermal_expansion_1_k": 3.0e-5}


def _config(core, core_density=8000.0, mantle_cooling="conduction", **mantle_keys):
    core = dict({"layer_index": 0, "radius_fraction": _CORE_FRACTION, "type": "none",
                 "material": _material(core_density)}, **core)
    mantle = dict({"class": "solidliquid", "type": "none", "layer_index": 1, "radius_fraction": 1.0,
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
    {"class": "base"},                                                         # geometry only
    {"class": "solidliquid", "temperature_k": 0.0, "cooling": {"model": "off"}},  # the 0 K default
])
def test_a_layer_without_a_temperature_is_no_heat_sink(core):
    world, result = _solve(_config(core))
    assert result["layer_in_thermal_network"] == [False, True]
    # Nothing crosses the core's interface, so the mantle's lower half stays at the mantle's own temperature.
    assert result["layer_heat_flow_in"][1] == 0.0
    assert result["layer_heat_flow_out"][0] == 0.0
    r_core = world.core.radius_outer
    r_mid = 0.5 * (r_core + world.radius)
    for radius in (r_core * (1.0 + 1.0e-6), 0.5 * (r_core + r_mid)):
        assert world.get_temperature(radius) == pytest.approx(_MANTLE_TEMPERATURE, rel=1.0e-9)
    # The mantle still loses heat through the surface.
    assert result["layer_heat_flow_out"][1] > 0.0


def test_a_warm_core_still_couples():
    _, result = _solve(_config({"class": "solidliquid", "temperature_k": 1800.0, "cooling": {"model": "off"}}))
    assert result["layer_in_thermal_network"] == [True, True]
    assert result["layer_heat_flow_in"][1] > 0.0


def test_the_solve_reports_the_convecting_detail():
    _, result = _solve(_config({"class": "solidliquid", "temperature_k": 1800.0, "cooling": {"model": "off"}},
                               mantle_cooling="convection"))
    for key in ("layer_node_temperature", "layer_top_temperature", "layer_boundary_thickness",
                "layer_rayleigh_number", "layer_nusselt_number"):
        assert len(result[key]) == 2, key
    assert 0.0 < result["layer_boundary_thickness"][1] <= 0.4 * (1.0 - _CORE_FRACTION) * _RADIUS
    assert result["layer_node_temperature"][1] == _SURFACE_TEMPERATURE
    assert result["layer_top_temperature"][1] < _MANTLE_TEMPERATURE    # the adiabat cools on its way up


def test_floating_layers_end_on_the_solved_grid():
    # A compressible core that holds less mass than its starting geometry: its radius takes several passes.
    core = {"class": "physics", "is_volume_fixed": False, "mass_kg": 2.0e22,
            "material": {"model": "bm", "reference_density_kg_m3": 8000.0, "reference_bulk_modulus_pa": 1.3e11,
                         "bulk_modulus_derivative": 4.5}}
    world, result = _solve(_config(core))
    assert result["geometry_converged"]
    assert result["thermal_passes"] > 1
    slices = len(result["radius"]) // 2
    assert result["radius"][slices - 1] == world.core.radius_outer
    assert result["radius"][-1] == world.radius
    assert world.mantle.radius_inner == world.core.radius_outer
