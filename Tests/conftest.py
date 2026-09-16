"""Shared pytest fixtures for the TidalPy test suite.

Bundled worlds are built by name from a copy of the packaged world pack in the user's data directory
(copy-if-absent), so an edited packaged world or data file never reaches a machine that already holds an older
copy. The session-wide fixture below points the world pack at a temporary directory, so every test reads the
packaged files of the build under test whatever the user's data directory holds, and no test writes into it.
"""
import pytest


@pytest.fixture(scope="session", autouse=True)
def isolated_worlds_x_dir(tmp_path_factory):
    """Redirect the new-backend world pack data directory to a temporary folder for the whole session."""
    from TidalPy.structures_x.configs import worldpack

    directory = tmp_path_factory.mktemp("Worlds_x")
    original = worldpack.get_worlds_x_dir
    worldpack.get_worlds_x_dir = lambda: str(directory)
    yield str(directory)
    worldpack.get_worlds_x_dir = original
