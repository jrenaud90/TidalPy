"""Tests that the C++ logger is configured at startup from the ``[logging]`` configuration section."""
import TidalPy
from TidalPy.initialize import build_logging_config


def test_logging_config_mirrors_the_logging_section():
    """The mapped logger config matches the ``[logging]`` section."""
    config = build_logging_config()
    logging_config = TidalPy.config["logging"]
    assert config["file_level"] == logging_config["file_level"]
    # Test mode never writes a log file.
    assert config["log_to_file"] is False
    assert config["log_file_path"] == ""
    if TidalPy._in_jupyter and not logging_config["print_log_notebook"]:
        assert config["console_level"] == "off"
    else:
        assert config["console_level"] == logging_config["console_level"]


def test_reinit_reapplies_the_logging_config():
    """reinit re-applies the logger config without raising."""
    TidalPy.reinit()
    assert build_logging_config()["file_level"] == TidalPy.config["logging"]["file_level"]
