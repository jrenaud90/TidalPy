"""Tests that the C++ logger writes its file at a non-ASCII path, on Windows as on Linux and macOS.

spdlog writes to the console directly, out of reach of ``caplog`` and ``capfd``, so the file sink is observed instead.
"""
import pytest

from TidalPy.initialize import build_logging_config
from TidalPy.Utilities.logging import flush_logger, init_logger, log_warning


@pytest.fixture
def restore_logging():
    """Restore the package logging configuration after the test."""
    yield
    init_logger(build_logging_config())


@pytest.mark.parametrize("directory_name", ["logdir_José_漢", "éèê", "漢字"])
def test_log_file_opens_at_a_non_ascii_path(tmp_path, restore_logging, directory_name):
    """The log file lands at the path given, and no directory with a garbled name is created beside it."""
    log_directory = tmp_path / directory_name
    log_path = log_directory / "tidalpy_unicode.log"
    init_logger({
        "console_level": "off",
        "file_level":    "info",
        "log_to_file":   True,
        "log_file_path": str(log_path),
    })
    log_warning("message at a non-ASCII path")
    flush_logger()
    assert log_path.is_file()
    assert "message at a non-ASCII path" in log_path.read_text(encoding="utf-8")
    assert [entry.name for entry in tmp_path.iterdir()] == [directory_name]
