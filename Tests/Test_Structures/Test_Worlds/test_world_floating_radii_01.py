"""Layers with ``is_volume_fixed = false`` hold their mass through the EOS solve while their radii float."""
import math

import pytest

from TidalPy.Structures.configs.world_builder import construct_world

_RADIUS = 2.0e6            # [m]
_CORE_FRACTION = 0.5
_CORE_DENSITY = 8000.0     # [kg/m^3]
_MANTLE_DENSITY = 3300.0   # [kg/m^3]


def _shell_mass(density, radius_inner, radius_outer):
    return (4.0 / 3.0) * math.pi * density * (radius_outer ** 3 - radius_inner ** 3)


_REFERENCE_CORE_MASS = _shell_mass(_CORE_DENSITY, 0.0, _CORE_FRACTION * _RADIUS)


def _config(core_floats=False, mantle_floats=False, core_density=_CORE_DENSITY, **core_keys):
    core = {"class": "physics", "type": "none", "layer_index": 0, "radius_fraction": _CORE_FRACTION,
            "is_volume_fixed": not core_floats,
            "material": {"model": "constant", "reference_density_kg_m3": core_density}}
    core.update(core_keys)
    return {
        "schema_version": "0.2.0",
        "name": "floating",
        "type": "terrestrial",
        "radius_m": _RADIUS,
        "mass_kg": (_shell_mass(_CORE_DENSITY, 0.0, _CORE_FRACTION * _RADIUS)
                    + _shell_mass(_MANTLE_DENSITY, _CORE_FRACTION * _RADIUS, _RADIUS)),
        "layers": {
            "core": core,
            "mantle": {"class": "physics", "type": "none", "layer_index": 1, "radius_fraction": 1.0,
                       "is_volume_fixed": not mantle_floats,
                       "material": {"model": "constant", "reference_density_kg_m3": _MANTLE_DENSITY}},
        },
    }


def _solve(config, **kwargs):
    world = construct_world(config)
    result = world.solve_eos(**kwargs)
    assert result["success"], result["message"]
    return world, result


def test_is_volume_fixed_defaults_to_true_and_round_trips():
    world, _ = _solve(_config())
    assert world.core.is_volume_fixed is True
    assert world.get_config_dict()["layers"]["core"]["is_volume_fixed"] is True
    floating, _ = _solve(_config(core_floats=True))
    assert floating.core.is_volume_fixed is False
    rebuilt = construct_world(floating.get_config_dict())
    assert rebuilt.core.is_volume_fixed is False


def test_fixed_volume_world_keeps_every_radius():
    world, result = _solve(_config())
    assert world.core.radius_outer == pytest.approx(_CORE_FRACTION * _RADIUS, rel=1e-15)
    assert world.radius == pytest.approx(_RADIUS, rel=1e-15)
    assert result["geometry_converged"]
    assert result["thermal_passes"] == 0


@pytest.mark.parametrize("core_density, core_mass", [
    pytest.param(_CORE_DENSITY, 0.5 * _REFERENCE_CORE_MASS, id="half_mass"),
    pytest.param(2.0 * _CORE_DENSITY, _REFERENCE_CORE_MASS, id="double_density"),
])
def test_a_floating_core_holds_its_mass(core_density, core_mass):
    """A floating core holding m / rho half its geometry's value shrinks by the cube root of two."""
    world, result = _solve(_config(core_floats=True, core_density=core_density, mass_kg=core_mass))
    assert result["geometry_converged"]
    expected_radius = _CORE_FRACTION * _RADIUS / 2.0 ** (1.0 / 3.0)
    assert world.core.radius_outer == pytest.approx(expected_radius, rel=1e-9)
    assert world.core.mass == pytest.approx(core_mass, rel=1e-9)


def test_the_layers_above_a_floating_layer_keep_their_volume():
    """The mantle keeps its volume and moves down with the shrunken core, and the world shrinks."""
    original_mantle_volume = (4.0 / 3.0) * math.pi * (_RADIUS ** 3 - (_CORE_FRACTION * _RADIUS) ** 3)
    world, _ = _solve(_config(core_floats=True, mass_kg=0.5 * _REFERENCE_CORE_MASS))
    assert world.mantle.volume == pytest.approx(original_mantle_volume, rel=1e-9)
    assert world.mantle.radius_inner == pytest.approx(world.core.radius_outer, rel=1e-15)
    assert world.radius < _RADIUS
    assert world.radius == pytest.approx(world.mantle.radius_outer, rel=1e-15)


def test_both_layers_can_float():
    core_mass = 0.8 * _REFERENCE_CORE_MASS
    config = _config(core_floats=True, mantle_floats=True, mass_kg=core_mass)
    mantle_mass = 1.2 * _shell_mass(_MANTLE_DENSITY, _CORE_FRACTION * _RADIUS, _RADIUS)
    config["layers"]["mantle"]["mass_kg"] = mantle_mass
    world, result = _solve(config)
    assert result["geometry_converged"]
    assert world.core.mass == pytest.approx(core_mass, rel=1e-9)
    assert world.mantle.mass == pytest.approx(mantle_mass, rel=1e-9)
    assert world.planet_mass_eos == pytest.approx(core_mass + mantle_mass, rel=1e-9)


def test_a_compressible_floating_layer_conserves_its_mass():
    """A compressible floating core holds its mass exactly in a smaller radius than the incompressible one."""
    config = _config(core_floats=True)
    config["layers"]["core"]["material"] = {
        "model": "bm",
        "reference_density_kg_m3": _CORE_DENSITY,
        "reference_bulk_modulus_pa": 1.3e11,
        "bulk_modulus_derivative": 4.5}
    reference = construct_world(_config())
    reference.solve_eos()
    target_mass = reference.core.mass
    config["layers"]["core"]["mass_kg"] = target_mass
    world, result = _solve(config)
    assert result["geometry_converged"]
    assert world.core.mass == pytest.approx(target_mass, rel=1e-6)
    assert world.core.radius_outer < reference.core.radius_outer


def test_reset_layer_masses_takes_the_current_geometry():
    """A floating layer without a configured mass holds its first solve's mass until reset_layer_masses redefines it."""
    world, _ = _solve(_config(core_floats=True))
    first_radius = world.core.radius_outer
    assert first_radius == pytest.approx(_CORE_FRACTION * _RADIUS, rel=1e-9)
    world.core.set_radii(0.0, 0.4 * _RADIUS)
    world.mantle.set_radii(0.4 * _RADIUS, _RADIUS)
    world.solve_eos(reset_layer_masses=True)
    assert world.core.radius_outer == pytest.approx(0.4 * _RADIUS, rel=1e-9)
    assert world.core.mass == pytest.approx(_shell_mass(_CORE_DENSITY, 0.0, 0.4 * _RADIUS), rel=1e-6)


def test_a_tiny_floating_layer_shrinks_to_almost_nothing():
    """A layer holding far less mass than its geometry implies collapses toward the center and holds it."""
    world, result = _solve(_config(core_floats=True, mass_kg=1.0e16))
    assert result["geometry_converged"]
    assert world.core.radius_outer < 0.01 * _RADIUS
    assert world.core.mass == pytest.approx(1.0e16, rel=1e-6)
    assert world.mantle.radius_inner == pytest.approx(world.core.radius_outer, rel=1e-15)


@pytest.mark.parametrize("world_name", ["io", "europa", "luna", "mercury", "earth_simple", "earth_prem"])
def test_bundled_worlds_hold_their_volume(world_name):
    from TidalPy.Structures.configs import build_world
    world = build_world(world_name)
    assert all(layer.is_volume_fixed for layer in world)
    radii_before = [layer.radius_outer for layer in world]
    result = world.solve_eos()
    assert [layer.radius_outer for layer in world] == radii_before
    assert result["geometry_converged"]
    assert result["thermal_passes"] == 0
