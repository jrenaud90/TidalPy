"""Shared pytest fixtures for the TidalPy test suite.

The world pack is redirected to a temporary directory so tests read the packaged worlds of the build under test
(the user's copy is copy-if-absent and may be stale) and never write into the user's data directory.
"""
import pytest


@pytest.fixture(scope="session", autouse=True)
def isolated_worlds_dir(tmp_path_factory):
    """Redirect the world pack data directory to a temporary folder for the whole session."""
    from TidalPy.Structures.configs import worldpack

    directory = tmp_path_factory.mktemp("Worlds")
    original = worldpack.get_worlds_dir
    worldpack.get_worlds_dir = lambda: str(directory)
    yield str(directory)
    worldpack.get_worlds_dir = original
