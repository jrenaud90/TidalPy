# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
logger.pyx
Python-facing interface to TidalPy's C++ spdlog logger.

At module-init time this extension creates the TidalPy spdlog logger with a
default console sink and stores a stable raw pointer to it.  Other Cython
extensions obtain that address via get_tidalpy_logger_address() and set their
own DLL-local tidalpy_logger_ptr via set_tidalpy_logger_ptr_void() — mirroring
the constants.pyx / constants_.hpp pattern.

Called once at TidalPy startup from the package __init__ after config loads::

    from TidalPy.Utilities_x.logging_x.logger import init_logger
    init_logger(config.get("logging", {}))
"""

from TidalPy.Utilities_x.logging_x.logger cimport (
    c_LoggerConfig,
    cy_create_default_logger,
    cy_init_logger,
    cy_set_log_level,
    cy_shutdown_logger,
    cy_get_logger_ptr,
)

# =====================================================================================================================
# Module-Init: create logger immediately for a stable pointer address
# =====================================================================================================================

# This runs when the extension is first imported, before init_logger() is called.
# Creates the TidalPy spdlog logger (stdout, info level) and sets tidalpy_logger_ptr.
cy_create_default_logger()


# =====================================================================================================================
# Cross-DLL Pointer Export  (mirrors constants.pyx get_shared_config_address)
# =====================================================================================================================

# Returns the raw address of the TidalPy spdlog logger for cross-DLL sharing.
# Other Cython extensions call at their module-init level:
#
#   from TidalPy.Utilities_x.logging_x.logger cimport (
#       set_tidalpy_logger_ptr_void, get_tidalpy_logger_address)
#   set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
#
# The returned void* is the non-owning raw pointer to the spdlog::logger
# instance owned by this extension's spdlog registry.
# Valid for the lifetime of the process.
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


cdef int _resolve_level(object level) except -1:
    """Convert a level name or integer to the spdlog level integer.

    Parameters
    ----------
    level : str or int
        Level name (case-insensitive) or integer 0–6.

    Returns
    -------
    int
        spdlog level enum value.

    Raises
    ------
    ValueError
        If level is a string not in the level map, or an integer outside 0–6.
    TypeError
        If level is not a str or int.
    """
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

    Replaces the sinks on the existing logger (created at module-import time)
    with properly configured sinks based on the user's config.  Because all
    Cython extensions share the same raw pointer, every DLL immediately sees
    the updated sink configuration.

    Parameters
    ----------
    config : dict, optional
        Sub-dict from TidalPy's global config containing logging settings.
        Recognized keys:

        ``console_level`` : str or int, default ``"info"``
            Log level for console output.
        ``file_level`` : str or int, default ``"info"``
            Log level for file output (only used when log_to_file is true).
        ``log_to_file`` : bool, default ``False``
            Whether to write logs to a file.
        ``log_file_path`` : str, default ``""``
            Absolute path for the log file (UTF-8 encoded).

    Assumptions
    -----------
    - Called before any C++ code emits TIDALPY_LOG_* messages.
    - Not thread-safe with concurrent logging; safe in practice at startup.

    Examples
    --------
    >>> from TidalPy.Utilities_x.logging_x.logger import init_logger
    >>> init_logger({"console_level": "debug", "log_to_file": False})
    """
    cdef c_LoggerConfig c_config

    if config is None:
        config = {}

    c_config.console_level = _resolve_level(config.get("console_level", "info"))
    c_config.file_level    = _resolve_level(config.get("file_level", "info"))
    c_config.log_to_file   = True if config.get("log_to_file", False) else False

    log_path = config.get("log_file_path", "")
    c_config.log_file_path = (log_path.encode("utf-8") if isinstance(log_path, str)
                              else b"")

    cy_init_logger(c_config)


def set_log_level(level):
    """Set the active log level on the TidalPy logger and all its sinks.

    Can be called at any time after import to adjust verbosity at runtime.

    Parameters
    ----------
    level : str or int
        Log level name (case-insensitive) or integer 0–6.
        Valid names: ``"trace"``, ``"debug"``, ``"info"``, ``"warning"``/``"warn"``,
        ``"error"``, ``"critical"``, ``"off"``.

    Examples
    --------
    >>> from TidalPy.Utilities_x.logging_x.logger import set_log_level
    >>> set_log_level("debug")
    >>> set_log_level(1)
    """
    cdef int int_level = _resolve_level(level)
    cy_set_log_level(int_level)


def shutdown_logger():
    """Flush pending log messages and make all TIDALPY_LOG_* macros no-ops.

    Should be called on TidalPy shutdown (e.g. via atexit) to ensure all
    messages are flushed and file handles are released.  The underlying logger
    object is retained in spdlog's registry to prevent dangling pointers in
    other DLLs; it is released automatically at process exit.

    Examples
    --------
    >>> from TidalPy.Utilities_x.logging_x.logger import shutdown_logger
    >>> shutdown_logger()
    """
    cy_shutdown_logger()
