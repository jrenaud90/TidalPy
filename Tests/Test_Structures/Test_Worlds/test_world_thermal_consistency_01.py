"""The thermal network and the integrated profile agree where the material's thermal properties vary.

The network relates each layer's temperature to its neighbors through conduction resistances and, for a convecting
layer, the warming of its adiabat; the structure solve then integrates the profile itself. The two agree only when the
network's resistances and adiabats follow the material as the integration does: a conductivity that follows the
temperature (ice, k ~ T^-0.84), and an adiabat and conductivity that change through a melting range. A converged solve
then lands on the surface temperature, and each layer holds its own temperature where its cooling model applies it.
"""
from pathlib import Path

import numpy as np
import pytest

import TidalPy
from TidalPy.Structures import build_world
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds.base import BaseWorld

_CORE_RADIUS = 1.0e6       # [m]
_RADIUS = 1.5e6            # [m]
_CORE_TEMPERATURE = 300.0  # [K]
_SHELL_TEMPERATURE = 200.0  # [K]
_HOT_MANTLE_TEMPERATURE = 1900.0   # [K] the bundled Earth's mantle 300 K warmer


def _icy_world():
    """A rock core under a conducting shell of ice Ih, whose conductivity falls with temperature."""
    core = Layer("core", 0, 0.0, _CORE_RADIUS, material="simple_rock", temperature=_CORE_TEMPERATURE)
    shell = Layer("shell", 1, _CORE_RADIUS, _RADIUS, material="ice_ih", temperature=_SHELL_TEMPERATURE,
                  cooling="conduction")
    mass = 4.0 / 3.0 * np.pi * (3300.0 * _CORE_RADIUS**3 + 920.0 * (_RADIUS**3 - _CORE_RADIUS**3))
    world = BaseWorld("icy", _RADIUS, mass)
    world.add_layer(core)
    world.add_layer(shell)
    return world


@pytest.mark.parametrize("surface_temperature", [40.0, 100.0])
def test_a_shell_whose_conductivity_follows_the_temperature(surface_temperature):
    """The shell's conductivity changes by about a factor of two across it; the solve still reaches the surface
    temperature, and the shell holds its own temperature at its mid-radius."""
    world = _icy_world()
    result = world.solve_eos(solve_temperature=True, surface_temperature=surface_temperature)
    assert result["success"], result["message"]
    assert result["thermal_converged"]
    conductivity = world.shell.calc_state(
        world.get_pressure(np.array([_CORE_RADIUS, _RADIUS])),
        world.get_temperature(np.array([_CORE_RADIUS, _RADIUS])))["thermal_conductivity"]
    assert conductivity[1] / conductivity[0] > 1.4
    assert world.get_temperature(_RADIUS) == pytest.approx(surface_temperature, abs=1.0e-2)
    assert world.get_temperature(0.5 * (_CORE_RADIUS + _RADIUS)) == pytest.approx(_SHELL_TEMPERATURE, abs=1.0e-2)


@pytest.mark.parametrize("surface_temperature", [300.0, 800.0])
def test_a_mantle_partially_molten_at_depth_reaches_the_surface_temperature(surface_temperature):
    """The bundled Earth's convecting mantle, 300 K warmer than bundled, is partially molten through much of its
    interior, where both its adiabat and its conductivity change; the converged solve still lands on the surface
    temperature."""
    world = build_world(str(Path(TidalPy.__file__).parent / "WorldPack" / "earth_simple.toml"))
    world.mantle.temperature = _HOT_MANTLE_TEMPERATURE
    result = world.solve_eos(solve_temperature=True, surface_temperature=surface_temperature)
    assert result["success"], result["message"]
    assert result["thermal_converged"]
    mantle = world.mantle
    boundary = result["layer_boundary_thickness"][-1]
    interior = np.linspace(mantle.radius_inner + 2.0 * boundary, mantle.radius_outer - 2.0 * boundary, 200)
    melt = mantle.calc_state(world.get_pressure(interior), world.get_temperature(interior))["melt_fraction"]
    assert np.count_nonzero((melt > 0.0) & (melt < 1.0)) > 20
    assert world.get_temperature(world.radius) == pytest.approx(surface_temperature, abs=0.1)
