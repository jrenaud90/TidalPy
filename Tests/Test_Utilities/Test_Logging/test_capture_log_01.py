"""TidalPy.capture_log collects a block's log messages and restores the logger; a notebook prints warnings by
default, without timestamps and with LF line endings."""
import pathlib

import pytest

import TidalPy
from TidalPy import initialize
from TidalPy.Utilities.logging import (
    capture_log,
    get_logger_config,
    init_logger,
    log_error,
    log_info,
    log_warning,
    resolve_log_level,
    set_console_level,
)
from TidalPy.Utilities.logging.capture import split_log_records
from TidalPy.Utilities.logging.logger import print_pending_messages


@pytest.fixture
def configured_logger():
    """The logger as the configuration sets it, before and after the test."""
    init_logger(initialize.build_logging_config())
    yield
    init_logger(initialize.build_logging_config())


def test_capture_log_collects_warnings_and_above(configured_logger):
    before = get_logger_config()
    with capture_log() as records:
        log_info("an info line")
        log_warning("a warning line")
        log_error("an error line")
        assert records == []
    assert len(records) == 2
    assert "[warning] a warning line" in records[0]
    assert "[error] an error line" in records[1]
    assert get_logger_config() == before


def test_capture_log_level_and_file(configured_logger, tmp_path):
    path = pathlib.Path(tmp_path) / "captured.log"
    with capture_log(path, level="info") as records:
        log_info("kept at info")
    assert len(records) == 1 and "kept at info" in records[0]
    assert "kept at info" in path.read_text(encoding="utf-8")
    # A second capture into the same file appends and returns only its own messages.
    with capture_log(path) as records:
        log_warning("second block")
    assert len(records) == 1 and "second block" in records[0]
    assert "kept at info" in path.read_text(encoding="utf-8")


def test_capture_log_restores_the_logger_after_an_exception(configured_logger):
    set_console_level("error")
    before = get_logger_config()
    assert before["console_level"] == resolve_log_level("error")
    with pytest.raises(RuntimeError):
        with capture_log() as records:
            log_warning("before the error")
            raise RuntimeError("raised inside the block")
    assert get_logger_config() == before
    assert len(records) == 1


def test_capture_log_refuses_an_unknown_level():
    with pytest.raises(ValueError, match="Unknown log level"):
        with capture_log(level="verbose"):
            pass


def test_capture_log_is_exported_at_the_top_level():
    assert TidalPy.capture_log is capture_log


def test_split_log_records_joins_continuation_lines():
    text = ("[2026-10-02 12:00:00.000] [TidalPy] [warning] first\nsecond line\n"
            "[2026-10-02 12:00:01.000] [TidalPy] [error] x\n")
    assert split_log_records(text) == [
        "[2026-10-02 12:00:00.000] [TidalPy] [warning] first\nsecond line",
        "[2026-10-02 12:00:01.000] [TidalPy] [error] x"]


@pytest.mark.parametrize("print_log_notebook, console_level, expected", [
    (False, "info", "warning"),
    (False, "error", "error"),
    (True, "info", "info"),
])
def test_a_notebook_prints_warnings_by_default(monkeypatch, restore_config, print_log_notebook, console_level,
                                               expected):
    monkeypatch.setattr(initialize, "is_notebook", lambda: True)
    monkeypatch.setitem(TidalPy.config["logging"], "print_log_notebook", print_log_notebook)
    monkeypatch.setitem(TidalPy.config["logging"], "console_level", console_level)
    assert resolve_log_level(initialize.build_logging_config()["console_level"]) == resolve_log_level(expected)


def test_the_notebook_sink_prints_without_a_timestamp_and_with_lf(configured_logger, capsys):
    """A notebook message prints as `[TidalPy] [warning] text` and an LF on every platform, so a re-run notebook's
    output does not change."""
    # Drain anything an earlier notebook configuration kept.
    print_pending_messages()
    capsys.readouterr()
    init_logger({"console_level": "warning", "console_pending": True})
    log_info("below the notebook level")
    log_warning("a notebook warning")
    print_pending_messages()
    assert capsys.readouterr().err == "[TidalPy] [warning] a notebook warning\n"
