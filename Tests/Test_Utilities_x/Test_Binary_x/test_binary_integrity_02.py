"""Binary loads check a file's record sizes and byte order, and saves never destroy the file they replace.

``load_binary`` refuses a record whose header claims more payload than the file holds, a physics model record whose
payload size is not what this build reads, a file written in another byte order, and a file with bytes left over after
its record. ``save_binary`` writes a temporary file beside the target and renames it over the target once complete.
"""
import os
import stat
import struct
import sys

import pytest

from TidalPy.rheology_x import make_rheology
from TidalPy.structures_x import build_world
from TidalPy.structures_x.layers.base import BaseLayer
from TidalPy.Utilities_x.binary_x import check_binary_file


# =====================================================================================================================
# Helpers
# =====================================================================================================================
# Header layout: magic (4), schema major, minor, patch (1 each), byte order (1), class id (4), payload size (8).
BYTE_ORDER_OFFSET = 7
PAYLOAD_SIZE_OFFSET = 12
HOST_ORDER = "<" if sys.byteorder == "little" else ">"
HOST_BYTE_ORDER_FLAG = 0 if sys.byteorder == "little" else 1


def make_sundberg():
    return make_rheology("sundberg", {"alpha": 0.2})


def make_sundberg_target():
    return make_rheology("sundberg")


def make_layer():
    return BaseLayer("probe", 0, 0.0, 1.0e6, 1.0e20)


def make_layer_target():
    return BaseLayer("other", 0, 0.0, 2.0e6, 3.0e20)


def make_world():
    return build_world("io")


def sundberg_state(model):
    return model.get_config_dict()


def layer_state(layer):
    return (layer.name, layer.radius_outer, layer.mass)


def world_state(world):
    return [(layer.name, type(layer).__name__, layer.radius_outer) for layer in world]


# (saved object factory, load target factory, state reader)
CASES = {
    "rheology": (make_sundberg, make_sundberg_target, sundberg_state),
    "layer": (make_layer, make_layer_target, layer_state),
    "world": (make_world, make_world, world_state),
}


def read_bytes(path):
    with open(path, "rb") as file:
        return bytearray(file.read())


def write_bytes(path, data):
    with open(path, "wb") as file:
        file.write(bytes(data))


def leftover_temporary_files(directory):
    return [name for name in os.listdir(directory) if name.endswith(".partial")]


# =====================================================================================================================
# Loading
# =====================================================================================================================
@pytest.mark.parametrize("case", CASES)
def test_round_trip_still_works(tmp_path, case):
    make_saved, make_target, state = CASES[case]
    saved = make_saved()
    path = os.path.join(str(tmp_path), f"{case}.tpyb")
    saved.save_binary(path)
    target = make_target()
    target.load_binary(path)
    assert state(target) == state(saved)


@pytest.mark.parametrize("case", CASES)
def test_trailing_bytes_raise(tmp_path, case):
    make_saved, make_target, _ = CASES[case]
    path = os.path.join(str(tmp_path), f"{case}.tpyb")
    make_saved().save_binary(path)
    write_bytes(path, read_bytes(path) + b"trailing garbage")
    with pytest.raises(IOError, match="after the end of its record"):
        make_target().load_binary(path)


@pytest.mark.parametrize("case", CASES)
@pytest.mark.parametrize("force", (False, True))
def test_payload_larger_than_file_raises(tmp_path, case, force):
    make_saved, make_target, state = CASES[case]
    path = os.path.join(str(tmp_path), f"{case}.tpyb")
    make_saved().save_binary(path)
    data = read_bytes(path)
    data[PAYLOAD_SIZE_OFFSET:PAYLOAD_SIZE_OFFSET + 8] = struct.pack(f"{HOST_ORDER}Q", 10**15)
    write_bytes(path, data)
    target = make_target()
    before = state(target)
    with pytest.raises(IOError, match="corrupt or truncated"):
        target.load_binary(path, force)
    # The size is checked before the record is read, so nothing was loaded.
    assert state(target) == before


@pytest.mark.parametrize("force", (False, True))
def test_physics_record_with_an_extra_parameter_raises(tmp_path, force):
    path = os.path.join(str(tmp_path), "maxwell.tpyb")
    make_rheology("maxwell").save_binary(path)
    data = read_bytes(path)
    payload_size = struct.unpack(f"{HOST_ORDER}Q", bytes(data[PAYLOAD_SIZE_OFFSET:PAYLOAD_SIZE_OFFSET + 8]))[0]
    # A record written with one more parameter than this build reads: the payload and the file both grow by a double.
    data[PAYLOAD_SIZE_OFFSET:PAYLOAD_SIZE_OFFSET + 8] = struct.pack(f"{HOST_ORDER}Q", payload_size + 8)
    data += struct.pack(f"{HOST_ORDER}d", 1.0)
    write_bytes(path, data)
    with pytest.raises(IOError, match="payload bytes, but this TidalPy build reads"):
        make_rheology("maxwell").load_binary(path, force)


@pytest.mark.parametrize("wrong_size", (0, 3))
def test_physics_record_with_a_wrong_payload_size_raises(tmp_path, wrong_size):
    path = os.path.join(str(tmp_path), "sundberg.tpyb")
    make_sundberg().save_binary(path)
    data = read_bytes(path)
    data[PAYLOAD_SIZE_OFFSET:PAYLOAD_SIZE_OFFSET + 8] = struct.pack(f"{HOST_ORDER}Q", wrong_size)
    write_bytes(path, data)
    target = make_sundberg_target()
    before = target.get_config_dict()
    with pytest.raises(IOError, match="payload bytes"):
        target.load_binary(path)
    assert target.get_config_dict() == before


def test_the_byte_order_is_recorded_and_checked(tmp_path):
    path = os.path.join(str(tmp_path), "maxwell.tpyb")
    make_rheology("maxwell").save_binary(path)
    data = read_bytes(path)
    assert data[BYTE_ORDER_OFFSET] == HOST_BYTE_ORDER_FLAG

    data[BYTE_ORDER_OFFSET] = 1 - HOST_BYTE_ORDER_FLAG
    write_bytes(path, data)
    with pytest.raises(IOError, match="byte order"):
        make_rheology("maxwell").load_binary(path)
    with pytest.raises(IOError, match="byte order"):
        check_binary_file(path)

    data[BYTE_ORDER_OFFSET] = 7
    write_bytes(path, data)
    with pytest.raises(IOError, match="unknown"):
        make_rheology("maxwell").load_binary(path, True)


# =====================================================================================================================
# Saving
# =====================================================================================================================
def test_save_replaces_an_existing_file(tmp_path):
    path = os.path.join(str(tmp_path), "sundberg.tpyb")
    make_rheology("sundberg", {"alpha": 0.2}).save_binary(path)
    make_rheology("sundberg", {"alpha": 0.4}).save_binary(path)
    reloaded = make_sundberg_target()
    reloaded.load_binary(path)
    assert reloaded.get_config_dict()["alpha"] == 0.4
    assert leftover_temporary_files(str(tmp_path)) == []


def test_save_into_a_missing_directory_raises_and_writes_nothing(tmp_path):
    path = os.path.join(str(tmp_path), "missing", "maxwell.tpyb")
    with pytest.raises(IOError, match="cannot open file for writing"):
        make_rheology("maxwell").save_binary(path)
    assert os.listdir(str(tmp_path)) == []


def test_a_failed_save_leaves_the_previous_file_intact(tmp_path):
    directory = str(tmp_path)
    path = os.path.join(directory, "world.tpyb")
    make_world().save_binary(path)
    good_data = read_bytes(path)
    replacement = make_sundberg()

    if sys.platform == "win32":
        # A file held open without delete sharing (as Python opens files) cannot be replaced on Windows, so the save
        # fails at the rename, after the new record is fully written.
        with open(path, "rb"):
            with pytest.raises(IOError, match="cannot replace"):
                replacement.save_binary(path)
    else:
        # A read-only directory refuses the temporary file. Root ignores the permission bits, so skip there.
        original_mode = os.stat(directory).st_mode
        os.chmod(directory, stat.S_IRUSR | stat.S_IXUSR)
        try:
            probe_path = os.path.join(directory, "probe")
            try:
                with open(probe_path, "wb"):
                    pass
            except OSError:
                pass
            else:
                os.remove(probe_path)
                pytest.skip("the directory stays writable after chmod (running as root)")
            with pytest.raises(IOError, match="cannot open file for writing"):
                replacement.save_binary(path)
        finally:
            os.chmod(directory, original_mode)

    assert read_bytes(path) == good_data
    assert leftover_temporary_files(directory) == []
    reloaded = make_world()
    reloaded.load_binary(path)
    assert world_state(reloaded) == world_state(make_world())
