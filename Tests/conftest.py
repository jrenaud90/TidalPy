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

# TidalPy.paths.DATA_DIR_ENVIRONMENT_VARIABLE, spelled out: importing it would initialize TidalPy before it is set.
DATA_DIR_ENVIRONMENT_VARIABLE = "TIDALPY_DATA_DIR"

_TEST_DATA_DIR = tempfile.mkdtemp(prefix="tidalpy_tests_")
atexit.register(shutil.rmtree, _TEST_DATA_DIR, ignore_errors=True)
os.environ[DATA_DIR_ENVIRONMENT_VARIABLE] = _TEST_DATA_DIR
os.environ["TIDALPY_TEST_MODE"] = "1"
