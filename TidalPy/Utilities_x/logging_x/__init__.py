"""TidalPy C++ logging interface (spdlog wrapper)."""

from TidalPy.Utilities_x.logging_x.logger import (
    init_logger,
    set_log_level,
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

__all__ = [
    "init_logger",
    "set_log_level",
    "shutdown_logger",
    "flush_logger",
    "log_message",
    "log_trace",
    "log_debug",
    "log_info",
    "log_warning",
    "log_error",
    "log_critical",
]
