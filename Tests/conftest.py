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
# One BLAS thread per test process: TidalPy's solvers thread in C++, and with pytest-xdist each worker would otherwise
# reserve a BLAS buffer per logical processor.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

# The reduced suite (the default) skips tests marked `slow` and runs a subset of each heavily parametrized test: cases
# are taken in order while they bring a parameter value not yet run, so every value of every parameter runs at least
# once. `--full-suite` or TIDALPY_FULL_SUITE=1 runs everything; the GitHub workflow sets it for its "full" suite.
FULL_SUITE_ENVIRONMENT_VARIABLE = "TIDALPY_FULL_SUITE"
REDUCED_MAX_CASES = 8


def pytest_addoption(parser):
    parser.addoption("--full-suite", action="store_true", default=False,
                     help=f"Run every test and every parametrized case (also {FULL_SUITE_ENVIRONMENT_VARIABLE}=1).")


@pytest.hookimpl(tryfirst=True)
def pytest_cmdline_main(config):
    # With pytest-xdist's -n and no --dist, distribute by work stealing rather than xdist's default `load`, which hands
    # each worker a contiguous block of the collection up front, so the demo notebooks land on one worker in a row.
    if getattr(config.option, "numprocesses", None) and (config.option.dist == "no"):
        config.option.dist = "worksteal"


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: runs only in the full suite (--full-suite or TIDALPY_FULL_SUITE=1)")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--full-suite") or os.environ.get(FULL_SUITE_ENVIRONMENT_VARIABLE, "0") not in ("", "0"):
        return
    cases = {}
    reduced = set()
    for item in items:
        if (not item.get_closest_marker("slow")) and hasattr(item, "callspec"):
            cases.setdefault((item.path, item.cls, item.originalname), []).append(item)
    for group in cases.values():
        if len(group) <= REDUCED_MAX_CASES:
            reduced.update(group)
            continue
        seen = set()
        for item in group:
            values = {(name, repr(value)) for name, value in item.callspec.params.items()}
            if not values <= seen:
                reduced.add(item)
                seen |= values
    kept = [item for item in items
            if (not item.get_closest_marker("slow")) and ((not hasattr(item, "callspec")) or (item in reduced))]
    if len(kept) < len(items):
        kept_ids = set(map(id, kept))
        config.hook.pytest_deselected(items=[item for item in items if id(item) not in kept_ids])
        items[:] = kept


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
