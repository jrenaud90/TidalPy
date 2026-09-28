"""Tests for the data directory path helpers in ``TidalPy.paths`` (scoped to ``<major>.<minor>.X``)."""
import os

import pytest

import TidalPy
from TidalPy import paths


def test_data_version_is_major_minor_x():
    """The data version has integer major and minor parts and an ``X`` patch."""
    data_version = paths.get_data_version()
    parts = data_version.split(".")
    assert len(parts) == 3
    assert parts[2] == "X"
    assert parts[0].isdigit()
    assert parts[1].isdigit()


def test_data_version_matches_package_version():
    """The data version's major and minor match the leading integers of the package version."""
    pkg_parts = str(TidalPy.version).split(".")
    data_parts = paths.get_data_version().split(".")
    assert data_parts[0] == "".join(char for char in pkg_parts[0] if char.isdigit())
    assert data_parts[1] == "".join(char for char in pkg_parts[1] if char.isdigit())


@pytest.mark.parametrize("getter", ["get_config_dir", "get_log_dir", "get_worlds_dir"])
def test_data_dirs_contain_data_version(getter):
    """Each data directory exists and has the data version as one of its path segments."""
    data_version = paths.get_data_version()
    directory = getattr(paths, getter)()
    assert os.path.isdir(directory)
    assert data_version in directory.replace(os.sep, "/").split("/")
