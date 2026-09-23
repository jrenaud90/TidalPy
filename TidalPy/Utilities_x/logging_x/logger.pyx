# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Python interface to TidalPy's C++ spdlog logger.

Importing this extension creates the TidalPy spdlog logger with a default console sink and stores a
stable raw pointer to it; other Cython extensions wire their own DLL-local pointer from
``get_tidalpy_logger_address()``. ``init_logger`` runs once at TidalPy startup, after the config loads.
"""

from TidalPy.Utilities_x.logging_x.logger cimport (
    c_LoggerConfig,
    cy_create_default_logger,
    cy_init_logger,
    cy_set_log_level,
    cy_log_message,
    cy_flush_logger,
    cy_shutdown_logger,
    cy_get_logger_ptr,
)

# =====================================================================================================================
# Module-Init: create logger immediately for a stable pointer address
# =====================================================================================================================

# Runs at import, before init_logger: creates the logger (stdout, info) and sets tidalpy_logger_ptr.
cy_create_default_logger()


# =====================================================================================================================
# Cross-DLL Pointer Export
# =====================================================================================================================

# Non-owning raw address of the TidalPy spdlog logger, owned by this extension's spdlog registry and
# valid for the lifetime of the process. Other extensions pass it to set_tidalpy_logger_ptr_void at
# their own module-init level.
cdef api void* get_tidalpy_logger_address():
    return cy_get_logger_ptr()


# =====================================================================================================================
# Level Name :: spdlog Integer Mapping
# =====================================================================================================================

# spdlog level enum: trace=0, debug=1, info=2, warn=3, error=4, critical=5, off=6
_LEVEL_MAP: dict = {
    "trace":    0,
    "debug":    1,
    "info":     2,
    "warning":  3,
    "warn":     3,
    "error":    4,
    "critical": 5,
    "off":      6,
}


cdef int cy_resolve_level(object level) except -1:
    """Convert a level name (case-insensitive) or an integer 0 to 6 into the spdlog level integer."""
    cdef str key
    if isinstance(level, str):
        key = level.lower()
        if key not in _LEVEL_MAP:
            raise ValueError(
                f"Unknown log level '{level}'. "
                f"Valid levels: {list(_LEVEL_MAP.keys())}"
            )
        return _LEVEL_MAP[key]
    elif isinstance(level, int):
        if not (0 <= level <= 6):
            raise ValueError(
                f"Integer log level must be in range [0, 6], got {level}."
            )
        return level
    else:
        raise TypeError(
            f"Log level must be a str or int, got {type(level).__name__}."
        )


# =====================================================================================================================
# Public Python Functions
# =====================================================================================================================

def init_logger(dict config = None):
    """Initialize the TidalPy C++ logger from a configuration dictionary.

    Replaces the sinks on the logger created at import, so every DLL holding the shared pointer sees the
    new configuration immediately.

    Parameters
    ----------
    config : dict, optional
        Logging settings: ``console_level`` and ``file_level`` (name or integer, default ``"info"``),
        ``log_to_file`` (bool, default False), and ``log_file_path`` (str, used only when writing a file).

    Assumptions
    -----------
    Called at startup before any C++ code emits TIDALPY_LOG_* messages; not thread-safe against
    concurrent logging.
    """
    cdef c_LoggerConfig c_config

    if config is None:
        config = {}

    c_config.console_level = cy_resolve_level(config.get("console_level", "info"))
    c_config.file_level    = cy_resolve_level(config.get("file_level", "info"))
    c_config.log_to_file   = True if config.get("log_to_file", False) else False

    # `object`, not `str`: a config can hold a non-string here and the isinstance check below is what falls
    # back to an empty path. A `cdef str` would raise on the assignment instead.
    cdef object log_path = config.get("log_file_path", "")
    c_config.log_file_path = (log_path.encode("utf-8") if isinstance(log_path, str)
                              else b"")

    cy_init_logger(c_config)


def set_log_level(level):
    """Set the active log level on the TidalPy logger and all its sinks.

    Parameters
    ----------
    level : str or int
        Level name (``"trace"``, ``"debug"``, ``"info"``, ``"warning"``/``"warn"``, ``"error"``,
        ``"critical"``, ``"off"``; case-insensitive) or the equivalent integer 0 to 6.
    """
    cdef int int_level = cy_resolve_level(level)
    cy_set_log_level(int_level)


def flush_logger():
    """Flush every sink of the TidalPy C++ logger (file sinks buffer their output)."""
    cy_flush_logger()


def log_message(level, str message):
    """Emit ``message`` through the TidalPy C++ logger at ``level``.

    Cython and Python code in the new backend log through this function (or the level helpers below) so
    their messages reach the same sinks as the C++ ``TIDALPY_LOG_*`` macros.

    Parameters
    ----------
    level : str or int
        Log level name (case-insensitive) or integer 0-6; see :func:`set_log_level`.
    message : str
        Text to log.
    """
    cdef int c_level = cy_resolve_level(level)
    cy_log_message(c_level, message.encode("utf-8"))


def log_trace(str message):
    """Emit ``message`` at the trace level."""
    cy_log_message(0, message.encode("utf-8"))


def log_debug(str message):
    """Emit ``message`` at the debug level."""
    cy_log_message(1, message.encode("utf-8"))


def log_info(str message):
    """Emit ``message`` at the info level."""
    cy_log_message(2, message.encode("utf-8"))


def log_warning(str message):
    """Emit ``message`` at the warning level."""
    cy_log_message(3, message.encode("utf-8"))


def log_error(str message):
    """Emit ``message`` at the error level."""
    cy_log_message(4, message.encode("utf-8"))


def log_critical(str message):
    """Emit ``message`` at the critical level."""
    cy_log_message(5, message.encode("utf-8"))


def shutdown_logger():
    """Flush pending log messages and make all TIDALPY_LOG_* macros no-ops.

    Call on TidalPy shutdown, for example through atexit. The logger object stays in spdlog's registry so
    that raw pointers held by other DLLs cannot dangle; it is released at process exit.
    """
    cy_shutdown_logger()
