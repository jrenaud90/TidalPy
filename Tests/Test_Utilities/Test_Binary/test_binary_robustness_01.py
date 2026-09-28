"""Tests that binary loads refuse another class's record or a corrupt file and leave the object unchanged."""
import os
import struct

import pytest

from TidalPy.Rheology import Andrade, Maxwell, Sundberg
from TidalPy.Structures import build_system, build_world
from TidalPy.Structures.layers.base import BaseLayer


@pytest.mark.parametrize("target_class", [Maxwell, Andrade])
def test_a_record_of_another_class_is_refused(tmp_path, target_class):
    """A rheology record of another class raises IOError and leaves the target unchanged."""
    path = os.path.join(str(tmp_path), "sundberg.tpyb")
    Sundberg(alpha=0.2).save_binary(path)
    target = target_class()
    before = target.get_config_dict()
    with pytest.raises(IOError, match="class id"):
        target.load_binary(path)
    assert target.get_config_dict() == before
    reloaded = Sundberg()
    reloaded.load_binary(path)
    assert reloaded.get_config_dict()["alpha"] == 0.2


def test_a_layer_of_another_class_is_refused(tmp_path):
    """A layer record of another class raises IOError and leaves the target unchanged."""
    world = build_world("io")
    path = os.path.join(str(tmp_path), "mantle.tpyb")
    world.mantle.save_binary(path)
    standalone = BaseLayer("probe", 0, 0.0, 1.0e6, 1.0e20)
    with pytest.raises(IOError, match="class id"):
        standalone.load_binary(path)
    assert standalone.name == "probe"


def test_a_corrupt_string_length_raises_instead_of_allocating(tmp_path):
    """A huge string length raises IOError before any allocation."""
    path = os.path.join(str(tmp_path), "maxwell.tpyb")
    Maxwell().save_binary(path)
    with open(path, "rb") as file:
        data = bytearray(file.read())
    # The model name follows the 20-byte header.
    name_offset = 4 + 4 + 4 + 8
    data[name_offset:name_offset + 4] = struct.pack("<I", 0xFFFFFFF0)
    with open(path, "wb") as file:
        file.write(bytes(data))
    with pytest.raises(IOError, match="corrupt or truncated"):
        Maxwell().load_binary(path)


def test_a_truncated_system_file_leaves_the_system_unchanged(tmp_path):
    """A truncated System file raises IOError and leaves the System unchanged."""
    system = build_system("sol_system")
    path = os.path.join(str(tmp_path), "system.tpyb")
    system.save_binary(path)
    size = os.path.getsize(path)
    with open(path, "rb") as file:
        data = file.read()
    with open(path, "wb") as file:
        file.write(data[:int(0.9 * size)])

    target = build_system("sol_system")
    names_before = [world.name for world in target]
    with pytest.raises(IOError):
        target.load_binary(path)
    assert [world.name for world in target] == names_before
    assert len(target) == len(names_before)
    for index in range(len(target)):
        assert target.get_tidal_host_index(index) < len(target)
