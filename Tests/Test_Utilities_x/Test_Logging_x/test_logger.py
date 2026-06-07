"""
Tests for TidalPy.Utilities_x.logging_x.logger

Covers init_logger, set_log_level, shutdown_logger, and bad-input handling.
These tests require the Cython extension to be compiled first::

    uv pip install -v <repo_root>
"""

import pytest


# =====================================================================================================================
# Helpers
# =====================================================================================================================
def _import_logger():
    """Import the compiled logger module, skipping if not built."""
    try:
        from TidalPy.Utilities_x.logging_x import logger as _logger_mod
        return _logger_mod
    except ImportError:
        pytest.skip("TidalPy.Utilities_x.logging_x.logger not compiled — run uv pip install first.")


# =====================================================================================================================
# Init / Shutdown
# =====================================================================================================================
def test_init_logger_default():
    """init_logger with no arguments should not raise."""
    mod = _import_logger()
    mod.shutdown_logger()   # ensure clean state
    mod.init_logger()
    mod.shutdown_logger()


def test_init_logger_with_config():
    """init_logger accepts a dict of recognized keys."""
    mod = _import_logger()
    mod.shutdown_logger()
    mod.init_logger({
        "console_level": "debug",
        "file_level":    "info",
        "log_to_file":   False,
        "log_file_path": "",
    })
    mod.shutdown_logger()


def test_init_logger_idempotent():
    """Calling init_logger twice does not raise (second call is a no-op)."""
    mod = _import_logger()
    mod.shutdown_logger()
    mod.init_logger({"console_level": "info"})
    mod.init_logger({"console_level": "debug"})  # second call → no-op
    mod.shutdown_logger()


def test_shutdown_logger_idempotent():
    """Calling shutdown_logger when not initialised should not raise."""
    mod = _import_logger()
    mod.shutdown_logger()
    mod.shutdown_logger()  # second shutdown → no-op


# =====================================================================================================================
# set_log_level — valid inputs
# =====================================================================================================================

@pytest.mark.parametrize("level", [
    "trace", "debug", "info", "warning", "warn", "error", "critical", "off",
    "TRACE", "DEBUG", "INFO",          # case-insensitive strings
    0, 1, 2, 3, 4, 5, 6,               # integer values
])
def test_set_log_level_valid(level):
    """set_log_level accepts all valid string and integer levels."""
    mod = _import_logger()
    mod.shutdown_logger()
    mod.init_logger()
    mod.set_log_level(level)
    mod.shutdown_logger()


# =====================================================================================================================
# set_log_level — invalid inputs
# =====================================================================================================================
@pytest.mark.parametrize("bad_level", [
    "verbose",      # unknown name
    "INFO2",
    -1,             # out-of-range integer
    7,
    3.5,            # wrong type
    None,
])
def test_set_log_level_invalid(bad_level):
    """set_log_level raises for unknown or out-of-range levels."""
    mod = _import_logger()
    mod.shutdown_logger()
    mod.init_logger()
    with pytest.raises((ValueError, TypeError)):
        mod.set_log_level(bad_level)
    mod.shutdown_logger()


# =====================================================================================================================
# init_logger — invalid config values
# =====================================================================================================================
@pytest.mark.parametrize("bad_level", ["verbose", "INFO2", -1, 7])
def test_init_logger_bad_console_level(bad_level):
    """init_logger raises when console_level is invalid."""
    mod = _import_logger()
    mod.shutdown_logger()
    with pytest.raises((ValueError, TypeError)):
        mod.init_logger({"console_level": bad_level})
    mod.shutdown_logger()
