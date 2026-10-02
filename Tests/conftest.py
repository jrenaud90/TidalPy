"""Shared pytest setup for the TidalPy test suite.

The suite runs against a fresh, temporary TidalPy data directory (``TIDALPY_DATA_DIR``), set here before anything
imports TidalPy. So every test reads the packaged configuration, worlds, and materials of the build under test, not the
user's ``TidalPy_Configs.toml`` (whose overrides of tolerances, truncations, or threads would change results) or their
copy-if-absent world and material files (which may be stale), and nothing is written into the user's data directory.
Subprocesses and notebook kernels the tests start inherit the variable. Each pytest process (each xdist worker) gets its
own directory, removed when the process exits.
"""
import atexit
import os
import shutil
import tempfile

import pytest

# TidalPy.paths.DATA_DIR_ENVIRONMENT_VARIABLE, spelled out: importing it would initialize TidalPy before it is set.
DATA_DIR_ENVIRONMENT_VARIABLE = "TIDALPY_DATA_DIR"

_TEST_DATA_DIR = tempfile.mkdtemp(prefix="tidalpy_tests_")
atexit.register(shutil.rmtree, _TEST_DATA_DIR, ignore_errors=True)
os.environ[DATA_DIR_ENVIRONMENT_VARIABLE] = _TEST_DATA_DIR
os.environ["TIDALPY_TEST_MODE"] = "1"


@pytest.fixture
def spdlog_text(tmp_path):
    """Route the C++ logger to a temporary file for the test and hand back a reader for its text; the logger returns to
    the configured settings afterward."""
    from TidalPy.initialize import build_logging_config
    from TidalPy.Utilities.logging.logger import flush_logger, init_logger

    log_path = tmp_path / "tidalpy.log"
    init_logger({"console_level": "off", "file_level": "warning", "log_to_file": True,
                 "log_file_path": str(log_path)})

    def read():
        flush_logger()
        return log_path.read_text(encoding="utf-8") if log_path.exists() else ""

    yield read
    init_logger(build_logging_config())


@pytest.fixture
def restore_config():
    """Restore ``TidalPy.config``, the path it was loaded from, and the C++ settings fed from it after a test changes
    them."""
    import copy

    import TidalPy
    from TidalPy.constants import update_constants

    original = copy.deepcopy(TidalPy.config)
    original_path = TidalPy._config_path
    yield
    TidalPy.config = original
    TidalPy._config_path = original_path
    update_constants()


@pytest.fixture(scope="module")
def io():
    """The bundled Io with its interior solved, built once for each test module that asks for it."""
    from TidalPy.Structures import build_world

    world = build_world("io")
    world.solve_eos()
    return world
