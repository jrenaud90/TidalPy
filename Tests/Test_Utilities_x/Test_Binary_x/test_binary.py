"""
Tests for TidalPy.Utilities_x.binary_x.binary

Covers check_binary_file and get_current_schema_version.
These tests require the Cython extension to be compiled first::

    uv pip install -v <repo_root>
"""

import os
import struct
import sys
import tempfile
import pytest


# =====================================================================================================================
# Helpers
# =====================================================================================================================
def _import_binary():
    """Import the compiled binary module, skipping if not built."""
    try:
        from TidalPy.Utilities_x.binary_x import binary as _mod
        return _mod
    except ImportError:
        raise ImportError(
            "TidalPy.Utilities_x.binary_x.binary not compiled — run uv pip install first."
        )


_ENDIAN = "<" if sys.byteorder == "little" else ">"


def _write_tpyb_header(f, major, minor, patch, class_id, payload_size=0):
    """Write a well-formed TidalPy binary header to file object f."""
    f.write(b"TPYB")
    f.write(bytes([major, minor, patch, 0]))           # 3 version bytes + reserved
    f.write(struct.pack(f"{_ENDIAN}I", class_id))      # class_id (uint32, host order)
    f.write(struct.pack(f"{_ENDIAN}Q", payload_size))  # payload_size (uint64, host order)


# =====================================================================================================================
# get_current_schema_version
# =====================================================================================================================
def test_schema_version_format():
    """get_current_schema_version returns a 3-component dot-separated string."""
    mod = _import_binary()
    ver = mod.get_current_schema_version()
    parts = ver.split(".")
    assert len(parts) == 3
    assert all(p.isdigit() for p in parts)


def test_schema_version_value():
    """Current schema version is 0.2.0."""
    mod = _import_binary()
    assert mod.get_current_schema_version() == "0.2.0"


# =====================================================================================================================
# check_binary_file — valid inputs
# =====================================================================================================================
def test_check_binary_file_returns_expected_dict():
    """check_binary_file returns correct version and class_id for a well-formed header."""
    mod = _import_binary()
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        _write_tpyb_header(f, major=0, minor=2, patch=0, class_id=100, payload_size=256)
        path = f.name
    try:
        result = mod.check_binary_file(path)
        assert result["schema_major"] == 0
        assert result["schema_minor"] == 2
        assert result["schema_patch"] == 0
        assert result["schema_version"] == "0.2.0"
        assert result["class_id"] == 100
        assert result["payload_size"] == 256
    finally:
        os.unlink(path)


def test_check_binary_file_dict_has_all_keys():
    """check_binary_file dict contains exactly the expected keys."""
    mod = _import_binary()
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        _write_tpyb_header(f, major=0, minor=2, patch=0, class_id=1, payload_size=0)
        path = f.name
    try:
        result = mod.check_binary_file(path)
        expected_keys = {
            "schema_major", "schema_minor", "schema_patch",
            "schema_version", "class_id", "payload_size",
        }
        assert set(result.keys()) == expected_keys
    finally:
        os.unlink(path)


def test_check_binary_file_different_patch():
    """check_binary_file reports the file's schema version even when patch differs."""
    mod = _import_binary()
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        _write_tpyb_header(f, major=0, minor=2, patch=99, class_id=1, payload_size=0)
        path = f.name
    try:
        result = mod.check_binary_file(path)
        assert result["schema_patch"] == 99
        assert result["schema_version"] == "0.2.99"
    finally:
        os.unlink(path)


def test_check_binary_file_large_payload_size():
    """check_binary_file correctly reads a large payload_size (uint64)."""
    mod = _import_binary()
    large_size = 2**32 + 12345  # larger than uint32 max
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        _write_tpyb_header(f, major=0, minor=2, patch=0, class_id=201, payload_size=large_size)
        path = f.name
    try:
        result = mod.check_binary_file(path)
        assert result["payload_size"] == large_size
        assert result["class_id"] == 201
    finally:
        os.unlink(path)


# =====================================================================================================================
# check_binary_file — error inputs
# =====================================================================================================================
def test_check_binary_file_not_found():
    """check_binary_file raises FileNotFoundError for a missing path."""
    mod = _import_binary()
    with pytest.raises(FileNotFoundError):
        mod.check_binary_file("/nonexistent/path/to/file_abc123.tpyb")


def test_check_binary_file_wrong_magic():
    """check_binary_file raises IOError when magic bytes are wrong."""
    mod = _import_binary()
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        f.write(b"FAKE")
        f.write(bytes([0, 2, 0, 0]))
        f.write(struct.pack(f"{_ENDIAN}I", 1))
        f.write(struct.pack(f"{_ENDIAN}Q", 0))
        path = f.name
    try:
        with pytest.raises(IOError):
            mod.check_binary_file(path)
    finally:
        os.unlink(path)


def test_check_binary_file_truncated():
    """check_binary_file raises IOError for a file shorter than 20 bytes."""
    mod = _import_binary()
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        f.write(b"TPYB")  # only magic — header is incomplete
        path = f.name
    try:
        with pytest.raises(IOError):
            mod.check_binary_file(path)
    finally:
        os.unlink(path)


def test_check_binary_file_empty_file():
    """check_binary_file raises IOError for an empty file."""
    mod = _import_binary()
    with tempfile.NamedTemporaryFile(suffix=".tpyb", delete=False) as f:
        path = f.name  # write nothing
    try:
        with pytest.raises(IOError):
            mod.check_binary_file(path)
    finally:
        os.unlink(path)


def test_check_binary_file_type_error():
    """check_binary_file raises TypeError when path is not a str."""
    mod = _import_binary()
    with pytest.raises(TypeError):
        mod.check_binary_file(123)


def test_check_binary_file_bytes_type_error():
    """check_binary_file raises TypeError for bytes path (must be str)."""
    mod = _import_binary()
    with pytest.raises(TypeError):
        mod.check_binary_file(b"/some/path.tpyb")
