"""Saved world and system files carry a version header, and a data-file world saves its file reference."""
import math

import pytest

from TidalPy.constants import G
from TidalPy.Structures import build_system, build_world


def _read_saved(path):
    with open(path, encoding="utf-8") as file:
        text = file.read()
    assert "\r" not in text
    return text.splitlines(), text


@pytest.mark.parametrize("build, source, title, same_attribute", [
    pytest.param(build_world, "io", "TidalPy world configuration: Io", "name", id="world"),
    pytest.param(build_system, "sol_system", "TidalPy system configuration:", "num_worlds", id="system"),
])
def test_a_saved_file_starts_with_the_version_header(tmp_path, build, source, title, same_attribute):
    original = build(source)
    path = str(tmp_path / "copy.toml")
    original.save_to_toml(path)
    lines, _ = _read_saved(path)
    assert lines[0].startswith("# ===")
    assert title in lines[1]
    assert any(line.startswith("#  TidalPy version:") for line in lines[:6])
    assert any(line.startswith("#  CyRK version:") for line in lines[:6])
    # The header is a comment, so the file loads as before.
    assert getattr(build(path), same_attribute) == getattr(original, same_attribute)


def test_a_data_file_world_saves_its_file_reference(tmp_path):
    world = build_world("earth_prem")
    assert world.portable_config is not None
    assert world.portable_config["data_file"] == "PREM.csv"
    assert "layers" not in world.portable_config or "interpolate" not in str(world.portable_config["layers"])
    # The inner core, the outer core, and ten mantle and crust layers between PREM's discontinuities.
    assert len(world.source_config["layers"]) == 12

    path = str(tmp_path / "earth_prem_copy.toml")
    world.save_to_toml(path)
    _, text = _read_saved(path)
    assert 'data_file = "PREM.csv"' in text
    assert "interpolate" not in text and "density_kg_m3 = [" not in text

    rebuilt = build_world(path)
    assert rebuilt.num_layers == 12
    world.solve_eos(G_to_use=G)
    rebuilt.solve_eos(G_to_use=G)
    assert math.isclose(rebuilt.planet_mass_eos, world.planet_mass_eos, rel_tol=1.0e-12)


def test_a_layered_world_has_no_portable_config():
    assert build_world("io").portable_config is None
