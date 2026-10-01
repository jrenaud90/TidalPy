"""Tests that a failed binary load leaves the object it loads into as it was.

The state compared is the object's own binary record, which holds everything the object saves.
"""
import math
import struct

import pytest

from TidalPy.constants import G
from TidalPy.Rheology import make_rheology
from TidalPy.Structures import build_world
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.system import System
from TidalPy.Structures.worlds.stellar import StarWorld
from TidalPy.Utilities.classes import StructureBase


# Header layout: magic (4), schema major, minor, patch (1 each), byte order (1), class id (4), payload size (8).
HEADER_BYTES = 20
PAYLOAD_SIZE_OFFSET = 12
TRUNCATION_FRACTIONS = (0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)


def saved_bytes(obj, path):
    """The object's binary record, as save_binary writes it."""
    obj.save_binary(str(path))
    return path.read_bytes()


def make_system(name, companion_radius):
    """A star hosting one companion on an eccentric orbit."""
    star = StarWorld("star", 7.0e8, 1.9e30)
    companion = StarWorld("companion", companion_radius, 1.898e27)
    companion.set_spin_frequency(math.sqrt(G * (1.9e30 + 1.898e27) / 1.0e10 ** 3))
    system = System(name)
    system.add_world(star, is_star=True)
    system.add_world(companion, tidal_host="star", semi_major_axis=1.0e10, eccentricity=0.05)
    return system


# (factory for the saved object, factory for the load target); the two differ in every saved field that can differ.
CASES = {
    "structure": (lambda: StructureBase(7.0e6, 8.0e23), lambda: StructureBase(1.0e6, 2.0e22)),
    "rheology": (lambda: make_rheology("sundberg", {"alpha": 0.2}), lambda: make_rheology("sundberg")),
    "layer": (lambda: Layer("probe", 0, 0.0, 1.0e6, 1.0e20), lambda: Layer("other", 0, 0.0, 2.0e6, 3.0e20)),
    "world": (lambda: build_world("pluto"), lambda: build_world("io")),
    "system": (lambda: make_system("saved", 7.0e7), lambda: make_system("existing", 6.0e7)),
}


def corrupt_versions(data):
    """Truncations at several fractions, trailing bytes, and a record whose header claims only part of its payload."""
    versions = {f"truncated_{fraction}": data[:int(len(data) * fraction)] for fraction in TRUNCATION_FRACTIONS}
    versions["trailing_bytes"] = data + b"trailing garbage"
    # A root payload size that fits the truncated file lets the reader start on the record before it runs out.
    cut = HEADER_BYTES + (len(data) - HEADER_BYTES) // 2
    short_claim = bytearray(data[:cut])
    struct.pack_into("<Q", short_claim, PAYLOAD_SIZE_OFFSET, cut - HEADER_BYTES)
    versions["short_payload_claim"] = bytes(short_claim)
    return versions


@pytest.mark.parametrize("case", list(CASES))
def test_failed_load_leaves_the_target_unchanged(tmp_path, case):
    """Every corrupt version of a good file raises IOError and leaves the target's saved state as it was."""
    make_saved, make_target = CASES[case]
    data = saved_bytes(make_saved(), tmp_path / "good.tpyb")
    target = make_target()
    before = saved_bytes(target, tmp_path / "before.tpyb")
    assert before != data
    for label, corrupt in corrupt_versions(data).items():
        corrupt_path = tmp_path / f"{label}.tpyb"
        corrupt_path.write_bytes(corrupt)
        with pytest.raises(IOError):
            target.load_binary(str(corrupt_path))
        assert saved_bytes(target, tmp_path / "after.tpyb") == before, f"{case}: {label} changed the target"


@pytest.mark.parametrize("case", list(CASES))
def test_a_good_load_after_failed_loads_still_works(tmp_path, case):
    """A target that refused corrupt files still loads a good one completely."""
    make_saved, make_target = CASES[case]
    good_path = tmp_path / "good.tpyb"
    data = saved_bytes(make_saved(), good_path)
    target = make_target()
    corrupt_path = tmp_path / "corrupt.tpyb"
    corrupt_path.write_bytes(data[:len(data) // 2])
    with pytest.raises(IOError):
        target.load_binary(str(corrupt_path))
    target.load_binary(str(good_path))
    assert saved_bytes(target, tmp_path / "after.tpyb") == data


def test_audit_world_keeps_its_radius_and_layers(tmp_path):
    """Loading a truncated Pluto into Io keeps Io's radius and its own layers."""
    data = saved_bytes(build_world("pluto"), tmp_path / "pluto.tpyb")
    target = build_world("io")
    radius_before = target.radius
    layers_before = [(layer.name, type(layer).__name__, layer.radius_outer) for layer in target]
    for fraction in (0.1, 0.5, 0.9):
        truncated_path = tmp_path / "truncated.tpyb"
        truncated_path.write_bytes(data[:int(len(data) * fraction)])
        with pytest.raises(IOError):
            target.load_binary(str(truncated_path))
        assert target.radius == radius_before
        assert [(layer.name, type(layer).__name__, layer.radius_outer) for layer in target] == layers_before


def test_system_keeps_its_worlds_after_trailing_bytes(tmp_path):
    """Bytes after a system record, found only once the worlds are read, leave the system's worlds as they were."""
    data = saved_bytes(make_system("saved", 7.0e7), tmp_path / "good.tpyb")
    trailing_path = tmp_path / "trailing.tpyb"
    trailing_path.write_bytes(data + b"\x00")
    target = make_system("existing", 6.0e7)
    with pytest.raises(IOError, match="after the end of its record"):
        target.load_binary(str(trailing_path))
    assert target.name == "existing"
    assert [world.name for world in target] == ["star", "companion"]
    assert target.companion.radius == 6.0e7
