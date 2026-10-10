"""One-line reprs of worlds, layers, materials, phases, physics models, tide models, and radial solutions, and the
multi-line world summary."""
import re

import pytest

from TidalPy.Material import Material, Phase, load_material
from TidalPy.Rheology import Andrade, Maxwell
from TidalPy.Structures import build_world
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds import BaseWorld
from TidalPy.Tides.classes.tide import make_tide

# A repr is one line and shows no machine address.
ADDRESS = re.compile(r"0x[0-9a-fA-F]+")


def check_one_line(text):
    assert "\n" not in text
    assert ADDRESS.search(text) is None
    return text


@pytest.mark.parametrize("name, class_name", [
    ("io", "TerrestrialWorld"), ("jupiter_simple", "GasGiantWorld"), ("sol", "StarWorld")])
def test_world_repr(name, class_name):
    world = build_world(name)
    text = check_one_line(repr(world))
    assert text.startswith(f"{class_name}('{world.name}', radius_km=")
    assert f"num_layers={world.num_layers}" in text
    assert check_one_line(repr(BaseWorld.__new__(BaseWorld))) == "BaseWorld(no world)"


def test_layer_repr():
    layer = Layer("mantle", radius_inner=1.0e6, radius_outer=1.5e6, temperature=1600.0)
    assert check_one_line(repr(layer)) == (
        "Layer('mantle', index=0, radius_inner_km=1000, radius_outer_km=1500, state='auto', temperature_k=1600)")
    world = build_world("io")
    assert check_one_line(repr(world[1])).startswith(f"Layer('{world[1].name}', index=1,")


def test_material_and_phase_reprs():
    assert check_one_line(repr(Phase(eos="vinet", shear_modulus="constant"))) == (
        "Phase(eos='vinet', shear_modulus='constant')")
    peridotite = load_material("peridotite")
    text = check_one_line(repr(peridotite))
    assert text.startswith("Material(solid_eos=")
    assert "liquid_eos=" in text and "solidus=" in text
    assert check_one_line(repr(Material(liquid=Phase()))) == "Material(liquid_eos='constant')"


def test_physics_model_reprs_share_one_form():
    assert check_one_line(repr(Maxwell())) == "Maxwell('maxwell')"
    assert check_one_line(repr(Andrade(alpha=0.25))).startswith("Andrade('andrade', alpha=0.25, ")
    assert check_one_line(repr(make_tide("cpl", fixed_k=[0.3], fixed_q=[50.0]))) == (
        "FixedQTide('fixed_q', fixed_k=[0.3], fixed_q=[50])")
    assert check_one_line(repr(make_tide("rheology"))) == "RheologyTide('rheology')"


def test_radial_solution_repr(io):
    world = io.copy()
    world.solve_eos()
    world.solve_love_numbers(frequency=4.1e-5)
    solution = world.release_radial_solution()
    text = check_one_line(repr(solution))
    assert text.startswith("RadialSolverSolution(success=True, degree_l=2, num_ytypes=1, k=0.0")


def test_world_summary():
    world = build_world("io")
    unsolved = world.summary()
    lines = unsolved.splitlines()
    assert lines[0].startswith("TerrestrialWorld 'Io': radius 1821.49 km")
    assert "EOS not solved" in lines[0]
    assert len(lines) == 2 + world.num_layers
    assert lines[1].split()[:3] == ["layer", "radii", "[km]"]
    # Each row: name, inner "to" outer radius, state, material, density ("-" unsolved), temperature, rheology.
    for layer, line in zip(world, lines[2:]):
        assert line.split()[0] == layer.name
        assert line.split()[6] == "-"
    world.solve_eos()
    solved_lines = world.summary().splitlines()
    assert "EOS solved" in solved_lines[0]
    density_low, density_high = (float(value) for value in solved_lines[2].split()[6:9:2])
    assert 0.0 < density_low <= density_high
