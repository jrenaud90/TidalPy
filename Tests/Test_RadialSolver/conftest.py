"""Shared fixtures for the radial solver tests."""
import pytest

from TidalPy.initialize import build_logging_config
from TidalPy.Utilities.logging.logger import flush_logger, init_logger


@pytest.fixture
def spdlog_text(tmp_path):
    """Route the C++ logger to a temporary file and return a reader for its text."""
    log_path = tmp_path / "tidalpy.log"
    init_logger({"console_level": "off", "file_level": "warning", "log_to_file": True,
                 "log_file_path": str(log_path)})

    def read():
        flush_logger()
        return log_path.read_text(encoding="utf-8") if log_path.exists() else ""

    yield read
    init_logger(build_logging_config())
