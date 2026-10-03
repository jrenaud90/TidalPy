"""TidalPy C++ logging interface (spdlog wrapper)."""

from TidalPy.Utilities.logging.logger import (
    init_logger,
    get_logger_config,
    resolve_log_level,
    set_log_level,
    set_console_level,
    set_file_level,
    shutdown_logger,
    flush_logger,
    log_message,
    log_trace,
    log_debug,
    log_info,
    log_warning,
    log_error,
    log_critical,
)
from TidalPy.Utilities.logging.capture import capture_log

__all__ = [
    "init_logger",
    "get_logger_config",
    "resolve_log_level",
    "set_log_level",
    "set_console_level",
    "set_file_level",
    "shutdown_logger",
    "flush_logger",
    "log_message",
    "log_trace",
    "log_debug",
    "log_info",
    "log_warning",
    "log_error",
    "log_critical",
    "capture_log",
]
