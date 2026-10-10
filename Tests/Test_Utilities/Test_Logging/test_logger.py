"""Tests for init_logger, set_log_level, shutdown_logger, and message emission in the C++ logger."""
import pytest

from TidalPy.Utilities.logging import logger


@pytest.mark.parametrize(
    "init_args",
    [
        (),
        ({"console_level": "debug", "file_level": "info", "log_to_file": False, "log_file_path": ""},),
    ],
    ids=["default", "with_config"])
def test_init_logger(init_args):
    """init_logger accepts no arguments or a dict of recognized keys."""
    logger.shutdown_logger()
    logger.init_logger(*init_args)
    logger.shutdown_logger()


def test_init_logger_idempotent():
    """A second init_logger call reconfigures the sinks without raising."""
    logger.shutdown_logger()
    logger.init_logger({"console_level": "info"})
    logger.init_logger({"console_level": "debug"})
    logger.shutdown_logger()


def test_shutdown_logger_idempotent():
    """shutdown_logger on an uninitialized logger does not raise."""
    logger.shutdown_logger()
    logger.shutdown_logger()


@pytest.mark.parametrize("level", [
    "trace", "debug", "info", "warning", "warn", "error", "critical", "off",
    "TRACE", "DEBUG", "INFO",
    0, 1, 2, 3, 4, 5, 6,
])
def test_set_log_level_valid(level):
    """set_log_level accepts every valid string (any case) and integer level."""
    logger.shutdown_logger()
    logger.init_logger()
    logger.set_log_level(level)
    logger.shutdown_logger()


@pytest.mark.parametrize("bad_level", ["verbose", "INFO2", -1, 7, 3.5, None])
def test_set_log_level_invalid(bad_level):
    """set_log_level raises for unknown names, out of range integers, and wrong types."""
    logger.shutdown_logger()
    logger.init_logger()
    with pytest.raises((ValueError, TypeError)):
        logger.set_log_level(bad_level)
    logger.shutdown_logger()


@pytest.mark.parametrize("bad_level", ["verbose", "INFO2", -1, 7])
def test_init_logger_bad_console_level(bad_level):
    """init_logger raises for an invalid console_level."""
    logger.shutdown_logger()
    with pytest.raises((ValueError, TypeError)):
        logger.init_logger({"console_level": bad_level})
    logger.shutdown_logger()


def test_log_message_reaches_the_file_sink(tmp_path):
    """log_message and the level helpers write through a file sink."""
    log_path = tmp_path / "logging.log"
    logger.init_logger({
        "console_level": "off",
        "file_level": "trace",
        "log_to_file": True,
        "log_file_path": str(log_path),
    })
    try:
        logger.log_message("warning", "message one")
        logger.log_info("message two")
        logger.log_debug("message three")
        logger.flush_logger()
        text = log_path.read_text(encoding="utf-8")
    finally:
        # Restore the default console-only sinks.
        logger.init_logger()
    assert "message one" in text
    assert "[warning]" in text
    assert "message two" in text
    assert "message three" in text


def test_file_sink_honors_its_level(tmp_path):
    """Messages below the file level are dropped."""
    log_path = tmp_path / "logging_level.log"
    logger.init_logger({
        "console_level": "off",
        "file_level": "warning",
        "log_to_file": True,
        "log_file_path": str(log_path),
    })
    try:
        logger.log_info("quiet")
        logger.log_error("loud")
        logger.flush_logger()
        text = log_path.read_text(encoding="utf-8")
    finally:
        logger.init_logger()
    assert "quiet" not in text
    assert "loud" in text


@pytest.mark.parametrize("bad_level", ["loud", "7", 7, -1])
def test_log_message_bad_level(bad_level):
    """log_message rejects unknown level names and out of range integers."""
    with pytest.raises((ValueError, TypeError)):
        logger.log_message(bad_level, "text")
