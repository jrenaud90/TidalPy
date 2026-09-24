"""Level handling, flushing, and exit behavior of the new backend's C++ (spdlog) logger.

- ``set_log_level`` sets the logger-level threshold only; the sinks keep the levels ``init_logger`` gave them, and
  ``set_console_level`` / ``set_file_level`` change one sink each.
- Warnings and above are flushed as they are written, so they reach the file without ``flush_logger`` and survive
  a process that is killed outright.
- ``log_message("off", ...)`` emits nothing.
- ``flush_logger`` is registered with ``atexit`` exactly once, however often ``TidalPy.reinit()`` runs.

spdlog writes to the console directly, so neither ``caplog`` nor ``capfd`` sees its output reliably on every
platform; the file sink is the observable stand-in for sink levels. The console sink is set to ``"off"`` throughout.
"""
import os
import subprocess
import sys
from types import SimpleNamespace

import pytest

import TidalPy
from TidalPy.initialize import build_logging_x_config
from TidalPy.Utilities_x.logging_x import logger as logger_module
from TidalPy.Utilities_x.logging_x import (
    flush_logger,
    init_logger,
    log_critical,
    log_debug,
    log_error,
    log_info,
    log_message,
    log_warning,
    set_console_level,
    set_file_level,
    set_log_level,
)


# =====================================================================================================================
# Fixtures
# =====================================================================================================================
@pytest.fixture
def file_log(tmp_path):
    """Route the C++ logger to a temporary file; restore the package configuration afterward."""
    log_path = tmp_path / "tidalpy_x_levels.log"

    def configure(file_level="trace"):
        init_logger({
            "console_level": "off",
            "file_level":    file_level,
            "log_to_file":   True,
            "log_file_path": str(log_path),
        })

    def read(flush=True):
        if flush:
            flush_logger()
        return log_path.read_text(encoding="utf-8") if log_path.exists() else ""

    yield SimpleNamespace(configure=configure, read=read, path=log_path)
    init_logger(build_logging_x_config())


# =====================================================================================================================
# Logger-Level Threshold and Sink Levels
# =====================================================================================================================
def test_set_log_level_keeps_the_file_sink_level(file_log):
    """Lowering the logger threshold does not lower the file sink's level."""
    file_log.configure(file_level="warning")
    set_log_level("trace")
    log_info("below the file level")
    log_error("at or above the file level")
    text = file_log.read()
    assert "below the file level" not in text
    assert "at or above the file level" in text


def test_set_log_level_filters_every_sink(file_log):
    """Raising the logger threshold drops messages the file sink would otherwise write."""
    file_log.configure(file_level="trace")
    set_log_level("error")
    log_warning("below the threshold")
    log_critical("above the threshold")
    text = file_log.read()
    assert "below the threshold" not in text
    assert "above the threshold" in text


def test_init_logger_resets_the_threshold(file_log):
    """The next init_logger call resets the logger threshold to trace, so the sink levels apply again."""
    file_log.configure(file_level="trace")
    set_log_level("off")
    file_log.configure(file_level="trace")
    log_debug("after reinitialization")
    assert "after reinitialization" in file_log.read()


def test_set_file_level_changes_the_file_sink(file_log):
    """set_file_level lowers and raises the file sink's level on its own."""
    file_log.configure(file_level="error")
    log_info("before lowering")
    assert set_file_level("debug") is True
    log_debug("after lowering")
    assert set_file_level("critical") is True
    log_error("after raising")
    text = file_log.read()
    assert "before lowering" not in text
    assert "after lowering" in text
    assert "after raising" not in text


def test_set_console_level_leaves_the_file_sink(file_log):
    """Changing the console level does not change what the file sink writes."""
    file_log.configure(file_level="warning")
    assert set_console_level("critical") is True
    log_info("info stays out of the file")
    log_warning("warning still reaches the file")
    text = file_log.read()
    assert "info stays out of the file" not in text
    assert "warning still reaches the file" in text


def test_set_file_level_without_a_file_sink():
    """set_file_level reports False and changes nothing when no file is written."""
    init_logger({"console_level": "off", "log_to_file": False})
    try:
        assert set_file_level("debug") is False
    finally:
        init_logger(build_logging_x_config())


@pytest.mark.parametrize("setter", [set_console_level, set_file_level])
@pytest.mark.parametrize("bad_level", ["verbose", 7, -1, 3.5, None])
def test_sink_level_setters_reject_bad_levels(setter, bad_level):
    """The sink level setters validate their level as set_log_level does."""
    with pytest.raises((ValueError, TypeError)):
        setter(bad_level)


# =====================================================================================================================
# Flushing
# =====================================================================================================================
@pytest.mark.parametrize("emit", [log_warning, log_error, log_critical])
def test_warning_and_above_reach_the_file_without_flush(file_log, emit):
    """A warning or worse is on disk as soon as it is logged, along with the buffered lines before it."""
    file_log.configure(file_level="trace")
    log_info("buffered line before the warning")
    emit("flushed without flush_logger")
    text = file_log.read(flush=False)
    assert "buffered line before the warning" in text
    assert "flushed without flush_logger" in text


def test_reinit_does_not_register_the_exit_flush_again():
    """The atexit registration runs once, at the logger module's import; reinit re-configures the logger without
    importing the module again. (atexit offers no count of one function's registrations: unregister empties slots
    without shrinking _ncallbacks, so the module identity is what is checked. The exit flush itself is covered by
    the subprocess tests below.)"""
    module_before = sys.modules[logger_module.__name__]
    TidalPy.reinit()
    TidalPy.reinit()
    assert sys.modules[logger_module.__name__] is module_before


# Runs in a fresh interpreter from a neutral working directory, so the installed package is imported and the exit
# path is real. argv: log file path, then "exit" (normal interpreter exit) or "kill" (SIGTERM to itself: the process
# ends at once, as in a crash, without atexit, static destructors, or a stdio flush; TerminateProcess on Windows).
EXIT_SCRIPT = """
import os
import signal
import sys
from TidalPy.Utilities_x.logging_x import init_logger, log_info, log_warning
init_logger({
    "console_level": "off",
    "file_level": "trace",
    "log_to_file": True,
    "log_file_path": sys.argv[1],
})
log_info("info line before exit")
if sys.argv[2] == "kill":
    log_warning("warning line before kill")
    os.kill(os.getpid(), signal.SIGTERM)
"""


@pytest.mark.parametrize("exit_mode, expected", [
    ("exit", "info line before exit"),     # buffered info lines are flushed at a normal exit
    ("kill", "warning line before kill"),  # the warning flushed itself; nothing ran at exit
])
def test_log_file_is_complete_at_exit(tmp_path, exit_mode, expected):
    """A normal exit writes buffered lines; a process killed outright still keeps every warning."""
    log_path = tmp_path / "tidalpy_x_exit.log"
    environment = dict(os.environ)
    environment["TIDALPY_TEST_MODE"] = "1"
    result = subprocess.run(
        [sys.executable, "-c", EXIT_SCRIPT, str(log_path), exit_mode],
        capture_output=True,
        text=True,
        timeout=300,
        cwd=tmp_path,
        env=environment,
    )
    # The killed process has a nonzero return code by design; a traceback means the script itself failed.
    assert "Traceback" not in result.stderr, f"subprocess failed:\n{result.stderr}"
    if exit_mode == "exit":
        assert result.returncode == 0, f"subprocess failed:\n{result.stderr}"
    text = log_path.read_text(encoding="utf-8") if log_path.exists() else ""
    assert expected in text
    if exit_mode == "kill":
        # The warning's flush also wrote the info line buffered before it.
        assert "info line before exit" in text


# =====================================================================================================================
# The "off" Level
# =====================================================================================================================
@pytest.mark.parametrize("off_level", ["off", "OFF", 6])
def test_log_message_off_writes_nothing(file_log, off_level):
    """log_message at the off level emits nothing, even with every sink at trace."""
    file_log.configure(file_level="trace")
    log_message(off_level, "message logged at off")
    log_message("error", "marker after the off message")
    text = file_log.read()
    assert "message logged at off" not in text
    assert "marker after the off message" in text
