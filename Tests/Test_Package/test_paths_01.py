"""
Tests for TidalPy data-directory path helpers (``TidalPy.paths``).

The data/config directories are scoped to the package's ``<major>.<minor>.X``
version (not the full version), so every patch release of a given major.minor
shares the same directory.
"""

import os

import pytest

from TidalPy import paths


def test_data_version_is_major_minor_x():
    data_version = paths.get_data_version()
    parts = data_version.split(".")
    assert len(parts) == 3
    assert parts[2] == "X"
    # Major and minor are integers.
    assert parts[0].isdigit()
    assert parts[1].isdigit()


def test_data_version_matches_package_version():
    import TidalPy
    pkg_parts = str(TidalPy.version).split(".")
    data_parts = paths.get_data_version().split(".")
    # Leading-integer major/minor agree with the package version.
    assert data_parts[0] == "".join(c for c in pkg_parts[0] if c.isdigit())
    assert data_parts[1] == "".join(c for c in pkg_parts[1] if c.isdigit())


@pytest.mark.parametrize("getter", [
    "get_config_dir", "get_log_dir", "get_worlds_dir", "get_worlds_x_dir"])
def test_data_dirs_contain_data_version(getter):
    data_version = paths.get_data_version()
    directory = getattr(paths, getter)()
    assert os.path.isdir(directory)
    # The version-scoped segment appears in the path; the full patch version does not.
    assert data_version in directory.replace(os.sep, "/").split("/")
