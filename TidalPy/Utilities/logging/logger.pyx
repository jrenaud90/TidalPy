# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Python interface to TidalPy's C++ spdlog logger.

Importing this extension creates the logger with a default console sink and a stable raw pointer to it;
other extensions wire their own DLL-local pointer from ``get_tidalpy_logger_address()``. Importing it also
registers ``flush_logger`` with ``atexit`` so buffered lines reach the log file at interpreter exit. The module is
imported once per process, so ``TidalPy.reinit()`` never registers the hook a second time.
"""

import atexit
import numbers
import sys

from TidalPy.schema import LOG_LEVELS, LOG_LEVEL_RANGE

from TidalPy.Utilities.logging.logger cimport (
    c_LoggerConfig,
    cy_create_default_logger,
    cy_init_logger,
    cy_set_log_level,
    cy_set_console_level,
    cy_set_file_level,
    cy_log_message,
    cy_flush_logger,
    cy_drain_pending_messages,
    cy_shutdown_logger,
    cy_get_logger_ptr,
)

# At import, before init_logger, so the pointer address is stable: stdout sink at info level.
cy_create_default_logger()

# The logger's configuration as init_logger takes it: the last init_logger call's, with the level changes made since
# (set_console_level, set_file_level). get_logger_config hands out copies, and capture_log restores it.
p_logger_config = {
    "console_level": "info", "file_level": "info", "log_to_file": False, "log_file_path": "", "console_pending": False}


# Non-owning address of the logger, owned by this extension's spdlog registry for the life of the
# process. Other extensions pass it to set_tidalpy_logger_ptr_void at their own module-init.
cdef api void* get_tidalpy_logger_address():
    return cy_get_logger_ptr()


cdef int cy_resolve_level(object level) except -1:
    """A level name (case-insensitive, TidalPy.schema.LOG_LEVELS) or an integer 0 to 6, as the spdlog level integer
    (trace = 0 to off = 6). A bool is neither."""
    cdef str key
    if isinstance(level, str):
        key = level.lower()
        if key not in LOG_LEVELS:
            raise ValueError(
                f"Unknown log level '{level}'. "
                f"Valid levels: {list(LOG_LEVELS.keys())}"
            )
        return LOG_LEVELS[key]
    elif isinstance(level, numbers.Integral) and not isinstance(level, bool):
        if not (LOG_LEVEL_RANGE[0] <= level <= LOG_LEVEL_RANGE[1]):
            raise ValueError(
                f"Integer log level must be in range [{LOG_LEVEL_RANGE[0]}, {LOG_LEVEL_RANGE[1]}], got {level}."
            )
        return int(level)
    else:
        raise TypeError(
            f"Log level must be a str or int, got {type(level).__name__}."
        )


def resolve_log_level(level) -> int:
    """A log level as the spdlog level integer: trace 0, debug 1, info 2, warning 3, error 4, critical 5, off 6.

    Parameters
    ----------
    level : str or int
        A level name (case-insensitive, ``"warn"`` included) or the integer itself.

    Returns
    -------
    int

    Raises
    ------
    ValueError
        An unknown name or an integer outside 0 to 6.
    TypeError
        Neither a str nor an int.
    """
    return cy_resolve_level(level)


def get_logger_config() -> dict:
    """The logger's current configuration, in the form :func:`init_logger` takes, so passing it back restores it.

    Returns
    -------
    dict
        ``console_level``, ``file_level``, ``log_to_file``, and ``log_file_path``: the last :func:`init_logger` call's
        values, with the levels set since by :func:`set_console_level` and :func:`set_file_level`. A copy.
    """
    return dict(p_logger_config)


def init_logger(dict config = None):
    """Initialize the TidalPy C++ logger from a configuration dictionary.

    Replaces the console and file sinks behind the logger created at import, so every DLL holding the shared
    pointer sees the new configuration at once, and resets the logger-level threshold (``set_log_level``) to
    ``"trace"`` so only the sink levels filter.

    Parameters
    ----------
    config : dict, optional
        ``console_level`` and ``file_level`` (name or integer, default ``"info"``), ``log_to_file``
        (default False), ``log_file_path``, read only when writing a file, and ``console_pending`` (default False):
        keep the console's messages for :func:`print_pending_messages` instead of writing them to stdout, as a
        notebook does.

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
    c_config.console_pending = True if config.get("console_pending", False) else False

    # `object`, not `str`: a config can hold a non-string here, which the isinstance check below turns
    # into an empty path. `cdef str` would raise on the assignment first.
    cdef object log_path = config.get("log_file_path", "")
    c_config.log_file_path = (log_path.encode("utf-8") if isinstance(log_path, str)
                              else b"")

    cy_init_logger(c_config)
    p_logger_config.update({
        "console_level": config.get("console_level", "info"),
        "file_level": config.get("file_level", "info"),
        "log_to_file": True if c_config.log_to_file else False,
        "log_file_path": log_path if isinstance(log_path, str) else "",
        "console_pending": True if c_config.console_pending else False,
    })


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
    cdef cpp_bool updated = cy_set_console_level(int_level)
    if updated:
        p_logger_config["console_level"] = int_level
    return updated


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
    cdef cpp_bool updated = cy_set_file_level(int_level)
    if updated:
        p_logger_config["file_level"] = int_level
    return updated


def flush_logger():
    """Flush every sink.

    File sinks buffer lines below the warning level; warnings and above are flushed as they are written. Registered
    with ``atexit`` at import, so a normal interpreter exit flushes the file.
    """
    cy_flush_logger()


def print_pending_messages(*args):
    """Print the messages the notebook console sink kept since the last call to ``sys.stderr``.

    In a Jupyter notebook the console sink keeps its messages (``console_pending``), because the kernel does not show
    what C++ writes to the process's stdout on every platform; TidalPy registers this function as an IPython
    ``post_run_cell`` hook, so a cell's warnings print below it once it finishes. Each message prints as
    ``[TidalPy] [<level>] <text>`` with an LF line ending on every platform and no timestamp, so a re-run notebook's
    output stays the same; the console and log file keep their timestamps. The arguments IPython passes are ignored.
    """
    cdef vector[string] messages = cy_drain_pending_messages()
    cdef size_t i
    for i in range(messages.size()):
        sys.stderr.write(messages[i].decode("utf-8", errors="replace"))
    if messages.size() > 0:
        sys.stderr.flush()


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
    """Flush pending log messages and turn the logger off, so no ``TIDALPY_LOG_*`` macro or ``log_*`` call writes.

    Not needed at interpreter exit, where the ``atexit`` hook flushes. The logger stays in spdlog's registry, at the
    same address, so raw pointers held by other DLLs cannot dangle and an extension imported while it is off still
    reaches it. A later ``init_logger`` call (or ``set_log_level``) turns it back on for every extension.
    """
    cy_shutdown_logger()


# Flush buffered info and debug lines at a normal interpreter exit. A crash skips atexit, which is why warnings and
# above are flushed as they are written.
atexit.register(flush_logger)
