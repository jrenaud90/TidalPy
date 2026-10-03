# Logging (`Utilities.logging`)

_Updated: 2026-10-02_

TidalPy's compiled code logs through [spdlog](https://github.com/gabime/spdlog), with a thin Cython wrapper so Python can configure and write to the same logger. A single named logger, `"TidalPy"`, is created at package startup and shared by every compiled extension.

Most of the code producing the messages runs in C++, often with the interpreter lock released, so a warning raised inside a radial solve cannot call back into Python's `logging`. The logging therefore lives on the C++ side, and the Python entry points let Cython and Python code reach the same sinks.

```
TidalPy.__init__
  -> init_logger(config)          Python, in logger.pyx
       -> cy_init_logger(config)  C++, in logger_.hpp
            -> spdlog logger -> distribution sink -> console, and optionally a file
```

## Python API

```python
from TidalPy.Utilities.logging import (
    init_logger, get_logger_config, set_log_level, set_console_level, set_file_level, shutdown_logger,
    flush_logger, log_message, log_trace, log_debug, log_info, log_warning, log_error, log_critical)

# Called automatically at startup; each call reconfigures the sinks.
init_logger({
    "console_level": "info",
    "file_level":    "debug",
    "log_to_file":   True,
    "log_file_path": "/tmp/tidalpy.log",
})

set_log_level("warning")   # global threshold: drop everything below warning from every sink
set_console_level("info")  # console sink only, for example to see output in a notebook
set_file_level("trace")    # file sink only

log_warning("Surface solve is poorly conditioned; results may be unreliable.")
log_message("debug", "level by name or by integer from 0 to 6")

flush_logger()             # info and debug lines are buffered, so flush before reading the file
shutdown_logger()          # flush and turn the logger off
init_logger(get_logger_config())   # the current configuration, as init_logger takes it: turns the logger back on
```

### `capture_log(path=None, level="warning")`

`TidalPy.capture_log` (also `TidalPy.Utilities.logging.capture_log`) is a context manager that collects the log messages of a block of code. spdlog writes to the console and the log file directly, so Python's `logging` handlers and pytest's `caplog` never see TidalPy's messages; `capture_log` is how a script or a notebook inspects them.

```python
import TidalPy
from TidalPy.Structures import build_world

io = build_world("io")
io.solve_eos()
with TidalPy.capture_log() as records:           # Warnings and above, by default
    io.calc_tides(
        4.11e-5,
        4.11e-5,
        0.5,                                     # Past the range of the default eccentricity truncation (0.395)
        0.0,
        4.217e8,
        1.898e27)
print(len(records), records[0])                  # Each record as the log file holds it: time, logger, level, text
```

The block's messages at `level` and above go to a log file (`path`, appended to, or a temporary file removed afterward), and the list the block receives is filled with them when it ends, one entry per message. The console keeps printing as before, and a log file the logger was writing receives nothing during the block. On exit, through an exception included, the logger returns to the configuration it had on entry (`get_logger_config()`).

### `init_logger(config=None)`

Initialize or reconfigure the logger. It is safe to call repeatedly. Each call replaces the console and file sinks and resets the global threshold of `set_log_level` to `trace`. `TidalPy.reinit()` re-applies the package settings.

| Config key | Type | Default | Description |
|---|---|---|---|
| `console_level` | str or int | `"info"` | Level for console output. |
| `file_level` | str or int | `"info"` | Level for file output. |
| `log_to_file` | bool | `False` | Whether to write a log file at all. |
| `log_file_path` | str | `""` | Absolute path to that file. Non-ASCII paths work on every platform. |
| `console_pending` | bool | `False` | Keep the console's messages for `print_pending_messages` instead of writing them to stdout, as a notebook does (see [Notebooks](#notebooks)). |

`get_logger_config()` returns the current configuration in this form (the last `init_logger` call's, with the levels set since by `set_console_level` and `set_file_level`), so `init_logger(get_logger_config())` restores it. `resolve_log_level(level)` converts a level name or integer to the integer.

Level names are case-insensitive: `trace`, `debug`, `info`, `warning` or `warn`, `error`, `critical`, and `off`. Integers 0 through 6 are accepted in their place.

### Changing Levels at Runtime

Levels filter in two places. The logger level is a global threshold. A message below it reaches no sink. Each sink then applies its own level.

- `set_log_level(level)` sets the global threshold only. The sinks keep their levels, so it can only narrow what they write. In a notebook, where the console sink starts at warnings, `set_log_level("info")` leaves info messages silent.
- `set_console_level(level)` sets the console sink's level, for example `set_console_level("info")` to see info messages in a notebook. It returns `True`.
- `set_file_level(level)` sets the file sink's level. It returns `False` and changes nothing when no log file is being written.

Each takes a level name or integer, as in the table above. The next `init_logger` call, including `TidalPy.reinit()`, restores the configured sink levels and resets the global threshold to `trace`.

### `log_message(level, message)` and Level Helpers

`log_trace`, `log_debug`, `log_info`, `log_warning`, `log_error`, and `log_critical` each take just the message. `log_message` takes the level first, named or as an integer. `log_message("off", ...)` emits nothing. These are what TidalPy's Cython and Python code should use, so their output interleaves correctly with the messages the C++ macros emit.

### `flush_logger()` and `shutdown_logger()`

Warnings, errors, and critical messages are flushed to every sink as they are written, so they reach the log file even if the process later crashes. Trace, debug, and info lines are buffered, so a log file read while the interpreter is still running (e.g., while using a Jupyter notebook) may be missing them until `flush_logger` runs.

Importing the logging module registers `flush_logger` with `atexit`, so a normal interpreter exit writes the buffered lines. A crash skips `atexit`, and buffered lines below the warning level are lost.

`shutdown_logger` flushes and turns the logger off, making every `TIDALPY_LOG_*` macro a no-op. The logger keeps its address, so an extension imported while it is off is still wired to it, and a later `init_logger` call (or `set_log_level`) turns it back on for every extension. It is not needed at interpreter exit.

### Notebooks

A Jupyter kernel does not show what compiled code writes to the process's stdout on every platform (a Windows kernel never does), while Python's `sys.stderr` reaches the running cell. In a notebook the console sink therefore keeps its messages (`console_pending`), and TidalPy registers `print_pending_messages` as an IPython `post_run_cell` hook that writes them below the cell once it finishes. By default only warnings and above print there (`[logging] notebook_console_level`); `print_log_notebook = true` prints everything at `console_level` instead. A notebook message prints as `[TidalPy] [warning] <text>`, without the timestamp the console and the log file carry, and ends with an LF on every platform, so re-running a notebook leaves its saved output unchanged.

## C++ API

```cpp
#include "TidalPy/Utilities/logging/logger_.hpp"

TIDALPY_LOG_DEBUG("Initializing layer: {}", layer_name);
TIDALPY_LOG_WARN("EOS data not populated for layer {}", index);
TIDALPY_LOG_ERROR("Failed to open binary file: {}", path);
```

The macros are `TIDALPY_LOG_TRACE`, `TIDALPY_LOG_DEBUG`, `TIDALPY_LOG_INFO`, `TIDALPY_LOG_WARN`, `TIDALPY_LOG_ERROR`, and `TIDALPY_LOG_CRITICAL`. Each checks for a null logger and does nothing if one has not been created, so they are safe to call from any translation unit at any point in startup. Format strings use the Python-style braces that spdlog's bundled `fmt` provides.

## Sharing the Logger Pointer Across Extensions

`logger.pyx` creates the logger at import time and stores a raw pointer to it. The macros use that pointer directly instead of looking the logger up in spdlog's registry on every call, which keeps a debug-level log statement cheap enough to leave inside a solver loop.

As a result, every Cython extension using C++ logging has to wire itself to that pointer at module-init level, outside any function:

```cython
from TidalPy.Utilities.logging.logger cimport (
    set_tidalpy_logger_ptr_void, get_tidalpy_logger_address)
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
```

`get_tidalpy_logger_address` is a `cdef api` function, the same mechanism the shared config pointer uses. Because another module must `cimport` it, Cython imports the logger module first, which guarantees the logger exists before its address is taken.

On Linux and macOS the inline variable is process-wide and the pointer is already shared, so the call is a harmless no-op. On Windows each compiled extension is a separate DLL with its own copy of the variable, and the call is what connects that DLL to the shared logger. Writing it unconditionally makes the same source work on all three.

## Thread Safety

The logger object, and so the address every extension holds, never changes after it is created. It writes to a single spdlog distribution sink (`dist_sink_mt`), whose children are the console sink and the optional file sink. Reconfiguring while C++ threads are logging is safe. Each message goes entirely to the old sinks or entirely to the new ones. The configuration functions (`init_logger`, the level setters, `shutdown_logger`) are called from Python with the interpreter lock held, so they never race each other.

## Configuration

Package initialization (and `TidalPy.reinit`) configures this logger from the `[logging]` section of `TidalPy_Configs.toml`:

```toml
[logging]
use_cwd = true               # log directory: <run output dir>/Logs, or the TidalPy data directory's Logs folder
write_log_to_disk = false    # write a timestamped log file
file_level = "debug"
console_level = "info"
print_log_notebook = false   # in a notebook, print every message at console_level when true
notebook_console_level = "warning"   # otherwise a notebook prints this level and above (and no lower than console_level)
write_log_notebook = false   # no log file from a notebook unless true
```

The levels are `"trace"`, `"debug"`, `"info"`, `"warning"`, `"error"`, `"critical"`, and `"off"`. The log file is named `TidalPy_<YYYYMMDD-HHMMSS>.log`. Test mode (the `TIDALPY_TEST_MODE` environment variable) disables the file sink entirely.

spdlog writes to the console directly rather than through Python's `logging`, so pytest's `caplog` does not see these messages. A test that needs to assert on log output enables the file sink and reads the file, which is what `capture_log` and the `spdlog_text` fixture of `Tests/conftest.py` do.

## Dependencies

[spdlog v1.15.3](https://github.com/gabime/spdlog/releases/tag/v1.15.3), a git submodule at `Dependencies/spdlog`. It is header only, so there is no separate compilation step.
