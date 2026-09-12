# Logging (`Utilities_x.logging_x`)

_Updated: 2026-09-12_

TidalPy's compiled code logs through [spdlog](https://github.com/gabime/spdlog), wrapped thinly in Cython so Python can configure and write to the same logger. A single named logger, `"TidalPy"`, is created at package startup and shared by every compiled extension.

The reason this needs a module of its own rather than Python's `logging` is that most of the code producing the messages runs in C++, often with the interpreter lock released. A warning raised from inside a radial solve cannot call back into Python to log itself, so the logging has to live on the C++ side, and the Python entry points exist so that Cython and Python code in the new backend reach the same sinks rather than a parallel set.

```
TidalPy.__init__
  -> init_logger(config)          Python, in logger.pyx
       -> cy_init_logger(config)  C++, in logger_.hpp
            -> spdlog registry -> console, and optionally a file
```

## Python API

```python
from TidalPy.Utilities_x.logging_x import (
    init_logger, set_log_level, shutdown_logger, flush_logger,
    log_message, log_trace, log_debug, log_info, log_warning, log_error, log_critical)

# Called automatically at startup; each call reconfigures the sinks.
init_logger({
    "console_level": "info",
    "file_level":    "debug",
    "log_to_file":   True,
    "log_file_path": "/tmp/tidalpy.log",
})

set_log_level("debug")     # change verbosity at runtime

log_warning("Surface solve is poorly conditioned; results may be unreliable.")
log_message("debug", "level by name or by integer from 0 to 6")

flush_logger()             # file sinks buffer, so flush before reading the file
shutdown_logger()          # flush and release handles; also runs at interpreter exit
```

### `init_logger(config=None)`

Initialize or reconfigure the logger. Safe to call repeatedly; each call replaces the sinks, and `TidalPy.reinit()` re-applies the package settings.

| Config key | Type | Default | Description |
|---|---|---|---|
| `console_level` | str or int | `"info"` | Level for console output. |
| `file_level` | str or int | `"info"` | Level for file output. |
| `log_to_file` | bool | `False` | Whether to write a log file at all. |
| `log_file_path` | str | `""` | Absolute path to that file, UTF-8. |

Level names are case-insensitive: `trace`, `debug`, `info`, `warning` or `warn`, `error`, `critical`, and `off`. Integers 0 through 6 are accepted in their place.

### `set_log_level(level)`

Adjust the active level at runtime, updating both the logger and all of its sinks. A later `init_logger` call resets the logger-level filter, after which the configured per-sink levels apply again.

### `log_message(level, message)` and the level helpers

`log_trace`, `log_debug`, `log_info`, `log_warning`, `log_error`, and `log_critical` each take just the message. `log_message` takes the level first, named or as an integer. These are what Cython and Python code in the new backend should use, so their output interleaves correctly with the messages the C++ macros emit.

### `flush_logger()` and `shutdown_logger()`

`flush_logger` flushes every sink, which matters because file sinks buffer: a log file read immediately after a solve may be missing its last lines without it. `shutdown_logger` flushes and releases file handles, and runs automatically at interpreter exit.

## C++ API

```cpp
#include "TidalPy/Utilities_x/logging_x/logger_.hpp"

TIDALPY_LOG_DEBUG("Initializing layer: {}", layer_name);
TIDALPY_LOG_WARN("EOS data not populated for layer {}", index);
TIDALPY_LOG_ERROR("Failed to open binary file: {}", path);
```

The macros are `TIDALPY_LOG_TRACE`, `TIDALPY_LOG_DEBUG`, `TIDALPY_LOG_INFO`, `TIDALPY_LOG_WARN`, `TIDALPY_LOG_ERROR`, and `TIDALPY_LOG_CRITICAL`. Each checks for a null logger and does nothing if one has not been created, so they are safe to call from any translation unit at any point in startup. Format strings use the Python-style braces that spdlog's bundled `fmt` provides.

## Sharing the logger pointer across extensions

`logger.pyx` creates the logger at import time and stores a raw pointer to it. The macros use that pointer directly instead of looking the logger up in spdlog's registry on every call, which is what keeps a debug-level log statement cheap enough to leave inside a solver loop.

The consequence is that every Cython extension using C++ logging has to wire itself to that pointer at module-init level, outside any function:

```cython
from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void, get_tidalpy_logger_address)
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
```

`get_tidalpy_logger_address` is a `cdef api` function, the same mechanism the shared config pointer uses. Because another module must `cimport` it, Cython imports the logger module first, which guarantees the logger exists before its address is taken.

On Linux and macOS the inline variable is process-wide and the pointer is already shared, so the call is a harmless no-op. On Windows each compiled extension is a separate DLL with its own copy of the variable, and the call is what connects that DLL to the shared logger. Writing it unconditionally is what makes the same source work on all three.

## Configuration

Package initialization maps the classic `[logging]` section onto this logger, so one set of settings drives both logging systems:

```toml
[logging]
use_cwd = true               # log directory: <output dir>/Logs, or the TidalPy data path
write_log_to_disk = false    # enables this logger's own timestamped file
file_level = "DEBUG"
console_level = "INFO"
print_log_notebook = false   # console output is silenced in a notebook unless true
write_log_notebook = false   # no log file from a notebook unless true
```

The two loggers never share a file: this one writes its own timestamped `TidalPy_x` log in the same directory. Test mode disables the file sink entirely.

One practical consequence for test authors: spdlog writes to the console directly rather than through Python's `logging`, so `caplog` and friends do not see these messages. A test that needs to assert on log output has to enable the file sink and read the file.

## Dependencies

[spdlog v1.15.3](https://github.com/gabime/spdlog/releases/tag/v1.15.3), a git submodule at `Dependencies/spdlog`. Header only, so no separate compilation step.
