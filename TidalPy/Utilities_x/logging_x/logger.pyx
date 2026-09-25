# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Python interface to TidalPy's C++ spdlog logger.

Importing this extension creates the logger with a default console sink and a stable raw pointer to it;
other extensions wire their own DLL-local pointer from ``get_tidalpy_logger_address()``. Importing it also
registers ``flush_logger`` with ``atexit`` so buffered lines reach the log file at interpreter exit. The module is
imported once per process, so ``TidalPy.reinit()`` never registers the hook a second time.
"""

import atexit

from TidalPy.Utilities_x.logging_x.logger cimport (
    c_LoggerConfig,
    cy_create_default_logger,
    cy_init_logger,
    cy_set_log_level,
    cy_set_console_level,
    cy_set_file_level,
    cy_log_message,
    cy_flush_logger,
    cy_shutdown_logger,
    cy_get_logger_ptr,
)

# At import, before init_logger, so the pointer address is stable: stdout sink at info level.
cy_create_default_logger()


# Non-owning address of the logger, owned by this extension's spdlog registry for the life of the
# process. Other extensions pass it to set_tidalpy_logger_ptr_void at their own module-init.
cdef api void* get_tidalpy_logger_address():
    return cy_get_logger_ptr()


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
    """A level name (case-insensitive) or an integer 0 to 6, as the spdlog level integer."""
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


def init_logger(dict config = None):
    """Initialize the TidalPy C++ logger from a configuration dictionary.

    Replaces the console and file sinks behind the logger created at import, so every DLL holding the shared
    pointer sees the new configuration at once, and resets the logger-level threshold (``set_log_level``) to
    ``"trace"`` so only the sink levels filter.

    Parameters
    ----------
    config : dict, optional
        ``console_level`` and ``file_level`` (name or integer, default ``"info"``), ``log_to_file``
        (default False), and ``log_file_path``, read only when writing a file.

    Notes
    -----
    Safe while C++ threads are logging: the sinks are swapped under the lock the logger holds while it writes, so
    each message goes entirely to the old sinks or entirely to the new ones.
    """
    cdef c_LoggerConfig c_config

    if config is None:
        config = {}

    c_config.console_level = cy_resolve_level(config.get("console_level", "info"))
    c_config.file_level    = cy_resolve_level(config.get("file_level", "info"))
    c_config.log_to_file   = True if config.get("log_to_file", False) else False

    # `object`, not `str`: a config can hold a non-string here, which the isinstance check below turns
    # into an empty path. `cdef str` would raise on the assignment first.
    cdef object log_path = config.get("log_file_path", "")
    c_config.log_file_path = (log_path.encode("utf-8") if isinstance(log_path, str)
                              else b"")

    cy_init_logger(c_config)


def set_log_level(level):
    """Set the logger-level threshold: messages below it reach no sink.

    The console and file sinks keep their own levels (``set_console_level``, ``set_file_level``), so this only
    ever narrows what they write: ``set_log_level("info")`` in a notebook, where the console sink is off, leaves the
    console silent. The next ``init_logger`` call (including ``TidalPy.reinit()``) resets the threshold to
    ``"trace"``.

    Parameters
    ----------
    level : str or int
        ``"trace"``, ``"debug"``, ``"info"``, ``"warning"``/``"warn"``, ``"error"``, ``"critical"``, or
        ``"off"`` (case-insensitive), or the equivalent integer 0 to 6.
    """
    cdef int int_level = cy_resolve_level(level)
    cy_set_log_level(int_level)


def set_console_level(level):
    """Set the console sink's level, for example to turn console output back on in a notebook.

    The logger-level threshold (``set_log_level``) still applies on top of it. The next ``init_logger`` call
    (including ``TidalPy.reinit()``) restores the configured console level.

    Parameters
    ----------
    level : str or int
        A level name (case-insensitive) or the equivalent integer 0 to 6, as for ``set_log_level``.

    Returns
    -------
    bool
        True when the console sink was updated.
    """
    cdef int int_level = cy_resolve_level(level)
    return cy_set_console_level(int_level)


def set_file_level(level):
    """Set the file sink's level.

    The logger-level threshold (``set_log_level``) still applies on top of it. The next ``init_logger`` call
    (including ``TidalPy.reinit()``) restores the configured file level.

    Parameters
    ----------
    level : str or int
        A level name (case-insensitive) or the equivalent integer 0 to 6, as for ``set_log_level``.

    Returns
    -------
    bool
        True when the file sink was updated; False, changing nothing, when no log file is being written.
    """
    cdef int int_level = cy_resolve_level(level)
    return cy_set_file_level(int_level)


def flush_logger():
    """Flush every sink.

    File sinks buffer lines below the warning level; warnings and above are flushed as they are written. Registered
    with ``atexit`` at import, so a normal interpreter exit flushes the file.
    """
    cy_flush_logger()


def log_message(level, str message):
    """Emit ``message`` through the TidalPy C++ logger at ``level``.

    Cython and Python code logs through this (or the level helpers below) so its messages reach the same
    sinks as the C++ ``TIDALPY_LOG_*`` macros. A level of ``"off"`` (6) emits nothing.
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

    Not needed at interpreter exit, where the ``atexit`` hook flushes. The logger stays in spdlog's registry so raw
    pointers held by other DLLs cannot dangle; it is released at process exit. A later ``init_logger`` call turns it
    back on.
    """
    cy_shutdown_logger()


# Flush buffered info and debug lines at a normal interpreter exit. A crash skips atexit, which is why warnings and
# above are flushed as they are written.
atexit.register(flush_logger)
