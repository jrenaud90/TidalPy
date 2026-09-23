"""Binary loads refuse a file of another class and a corrupt or truncated file, and leave the object as it was.

``load_binary`` first compares the file's class id with the object's own, so a record of another class is never read
field by field into the wrong layout. Counts read from a file are checked against the bytes left in it before
anything is allocated, and a System reads its whole record before committing any of it.
"""
import os
import struct

import pytest

from TidalPy.rheology_x import Andrade, Maxwell, Sundberg
from TidalPy.structures_x import build_system, build_world
from TidalPy.structures_x.layers.base import BaseLayer


def test_a_record_of_another_class_is_refused(tmp_path):
    path = os.path.join(str(tmp_path), "sundberg.tpyb")
    Sundberg(alpha=0.2).save_binary(path)
    for target in (Maxwell(), Andrade()):
        before = target.get_config_dict()
        with pytest.raises(IOError, match="class id"):
            target.load_binary(path)
        assert target.get_config_dict() == before
    reloaded = Sundberg()
    reloaded.load_binary(path)
    assert reloaded.get_config_dict()["alpha"] == 0.2


def test_a_layer_of_another_class_is_refused(tmp_path):
    world = build_world("io")
    path = os.path.join(str(tmp_path), "mantle.tpyb")
    world.mantle.save_binary(path)
    standalone = BaseLayer("probe", 0, 0.0, 1.0e6, 1.0e20)
    with pytest.raises(IOError, match="class id"):
        standalone.load_binary(path)
    assert standalone.name == "probe"


def test_a_corrupt_string_length_raises_instead_of_allocating(tmp_path):
    path = os.path.join(str(tmp_path), "maxwell.tpyb")
    Maxwell().save_binary(path)
    with open(path, "rb") as file:
        data = bytearray(file.read())
    # The model name follows the 20-byte header: magic, three schema bytes and a reserved one, class id, payload size.
    name_offset = 4 + 4 + 4 + 8
    data[name_offset:name_offset + 4] = struct.pack("<I", 0xFFFFFFF0)
    with open(path, "wb") as file:
        file.write(bytes(data))
    with pytest.raises(IOError, match="corrupt or truncated"):
        Maxwell().load_binary(path)


def test_a_truncated_system_file_leaves_the_system_unchanged(tmp_path):
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
