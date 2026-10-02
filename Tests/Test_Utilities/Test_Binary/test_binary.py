"""Tests for ``check_binary_file`` and ``get_current_schema_version`` in ``TidalPy.Utilities.binary.binary``."""
import struct
import sys

import pytest

from TidalPy.Utilities.binary import binary

_ENDIAN = "<" if sys.byteorder == "little" else ">"


def _header_bytes(major, minor, patch, class_id, payload_size=0, magic=b"TPYB"):
    """Return a TidalPy binary header: magic, 3 version bytes + reserved, class id (uint32), payload size (uint64)."""
    return (
        magic
        + bytes([major, minor, patch, 0])
        + struct.pack(f"{_ENDIAN}I", class_id)
        + struct.pack(f"{_ENDIAN}Q", payload_size))


def _write_file(tmp_path, data):
    path = tmp_path / "file.tpyb"
    path.write_bytes(data)
    return str(path)


def test_schema_version_format():
    """The schema version is a 3 component dot separated string of integers."""
    parts = binary.get_current_schema_version().split(".")
    assert len(parts) == 3
    assert all(part.isdigit() for part in parts)


def test_schema_version_value():
    """The current schema version is 0.2.0."""
    assert binary.get_current_schema_version() == "0.2.0"


@pytest.mark.parametrize(
    "patch, class_id, payload_size",
    [
        (0, 100, 256),
        (0, 1, 0),
        (99, 1, 0),
        (0, 201, 2**32 + 12345),
    ],
    ids=["basic", "empty_payload", "different_patch", "payload_above_uint32"])
def test_check_binary_file_reads_header(tmp_path, patch, class_id, payload_size):
    """A well formed header is reported as a dict with exactly the header fields."""
    path = _write_file(tmp_path, _header_bytes(0, 2, patch, class_id, payload_size))
    assert binary.check_binary_file(path) == {
        "schema_major": 0,
        "schema_minor": 2,
        "schema_patch": patch,
        "schema_version": f"0.2.{patch}",
        "class_id": class_id,
        "payload_size": payload_size,
    }


def test_check_binary_file_not_found():
    """A missing path raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        binary.check_binary_file("/nonexistent/path/to/file_abc123.tpyb")


@pytest.mark.parametrize(
    "data",
    [_header_bytes(0, 2, 0, 1, magic=b"FAKE"), b"TPYB", b""],
    ids=["wrong_magic", "truncated", "empty"])
def test_check_binary_file_bad_header_raises(tmp_path, data):
    """Wrong magic bytes or a header shorter than 20 bytes raises IOError."""
    path = _write_file(tmp_path, data)
    with pytest.raises(IOError):
        binary.check_binary_file(path)


@pytest.mark.parametrize("path", [123, b"/some/path.tpyb"], ids=["int", "bytes"])
def test_check_binary_file_type_error(path):
    """A path that is not a str raises TypeError."""
    with pytest.raises(TypeError):
        binary.check_binary_file(path)


def test_the_binary_and_toml_schema_versions_agree():
    """The schema version has two homes, the C++ binary header (binary_.hpp) and the TOML schema (schema.py, which
    imports nothing from TidalPy so the configuration can load first); a bump to one must reach the other."""
    from TidalPy.schema import SCHEMA_VERSION
    assert binary.get_current_schema_version() == SCHEMA_VERSION
