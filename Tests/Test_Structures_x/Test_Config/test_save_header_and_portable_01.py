"""Saved world and system files carry a version header, and a data-file world saves its file reference."""
import math

from TidalPy.constants import G
from TidalPy.structures_x import build_system, build_world


def _header_lines(path):
    with open(path, encoding="utf-8") as file:
        text = file.read()
    assert "\r" not in text
    lines = text.splitlines()
    return lines, text


def test_a_saved_world_starts_with_the_version_header(tmp_path):
    world = build_world("io")
    path = str(tmp_path / "io_copy.toml")
    world.save_to_toml(path)
    lines, _ = _header_lines(path)
    assert lines[0].startswith("# ===")
    assert "TidalPy world configuration: Io" in lines[1]
    assert any(line.startswith("#  TidalPy version:") for line in lines[:6])
    assert any(line.startswith("#  CyRK version:") for line in lines[:6])
    # The header is a comment: the file loads as before.
    assert build_world(path).name == world.name


def test_a_saved_system_starts_with_the_version_header(tmp_path):
    system = build_system("sol_system")
    path = str(tmp_path / "system_copy.toml")
    system.save_to_toml(path)
    lines, _ = _header_lines(path)
    assert "TidalPy system configuration:" in lines[1]
    assert build_system(path).num_worlds == system.num_worlds


def test_a_data_file_world_saves_its_file_reference(tmp_path):
    world = build_world("earth_prem")
    assert world.portable_config is not None
    assert world.portable_config["data_file"] == "PREM.csv"
    assert "layers" not in world.portable_config or "interpolate" not in str(world.portable_config["layers"])
    # The expanded form is still the normalized configuration.
    assert len(world.source_config["layers"]) == 3

    path = str(tmp_path / "earth_prem_copy.toml")
    world.save_to_toml(path)
    _, text = _header_lines(path)
    assert 'data_file = "PREM.csv"' in text
    assert "interpolate" not in text and "density_kg_m3 = [" not in text

    rebuilt = build_world(path)
    assert rebuilt.num_layers == 3
    world.solve_eos(G_to_use=G)
    rebuilt.solve_eos(G_to_use=G)
    assert math.isclose(rebuilt.planet_mass_eos, world.planet_mass_eos, rel_tol=1.0e-12)


def test_a_layered_world_has_no_portable_config():
    assert build_world("io").portable_config is None
