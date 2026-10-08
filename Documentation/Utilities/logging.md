# Logging (`Utilities.logging`)

_Updated: 2026-10-06_

TidalPy logs through one C++ logger, `"TidalPy"`, built on [spdlog](https://github.com/gabime/spdlog) and shared by every compiled extension. Most messages come from C++ code running without the interpreter lock, where Python's `logging` cannot be called, so Python and Cython code write to the same C++ logger through the functions below. Package startup and `TidalPy.reinit` configure it from the `[logging]` section of `TidalPy_Configs.toml` (see [Configurations](../Overview/2_TidalPy_Configurations.md#logging)).

spdlog writes to the console and the log file directly, so Python's `logging` handlers and pytest's `caplog` never see TidalPy's messages. Use [`capture_log`](logging.md#capture_logpathnone-levelwarning) to inspect them.

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

Level names are case-insensitive: `trace`, `debug`, `info`, `warning` or `warn`, `error`, `critical`, and `off`. Integers 0 through 6 are accepted in their place.

### `capture_log(path=None, level="warning")`

`TidalPy.capture_log` (also `TidalPy.Utilities.logging.capture_log`) is a context manager that collects the log messages of a block of code:

```python
import TidalPy
from TidalPy.Structures import build_world

io = build_world("io")
io.solve_eos()
with TidalPy.capture_log() as records:           # Warnings and above, by default
    io.calc_tides(
        4.11e-5,
        4.11e-5,
        0.5,                                     # Past the range of the default eccentricity truncation (0.4 for a synchronous spin)
        0.0,
        4.217e8,
        1.898e27)
print(len(records), records[0])                  # Each record as the log file holds it: time, logger, level, text
```

The block's messages at `level` and above go to a log file (`path`, appended to, or a temporary file removed afterward), and the list fills with them, one entry per message, when the block ends. The console keeps printing, and a log file the logger was writing receives nothing during the block. On exit, through an exception included, the logger returns to its configuration on entry.

### `init_logger(config=None)`

Initialize or reconfigure the logger; safe to call repeatedly. Each call replaces the console and file sinks and resets the global threshold of `set_log_level` to `trace`. `TidalPy.reinit()` re-applies the package settings.

| Config key | Type | Default | Description |
|---|---|---|---|
| `console_level` | str or int | `"info"` | Level for console output. |
| `file_level` | str or int | `"info"` | Level for file output. |
| `log_to_file` | bool | `False` | Whether to write a log file at all. |
| `log_file_path` | str | `""` | Absolute path to that file. Non-ASCII paths work on every platform. |
| `console_pending` | bool | `False` | Keep the console's messages for `print_pending_messages` instead of writing them to stdout, as a notebook does (see [Notebooks](#notebooks)). |

`get_logger_config()` returns the current configuration in this form (the last `init_logger` call's, with any levels set since), so `init_logger(get_logger_config())` restores it. `resolve_log_level(level)` converts a level name or integer to the integer.

### Changing Levels at Runtime

Levels filter twice: the logger level is a global threshold, below which a message reaches no sink, and each sink then applies its own level.

- `set_log_level(level)` sets the global threshold only, so it can only narrow what the sinks write. In a notebook, where the console sink starts at warnings, `set_log_level("info")` leaves info messages silent.
- `set_console_level(level)` sets the console sink's level (`set_console_level("info")` shows info messages in a notebook) and returns `True`.
- `set_file_level(level)` sets the file sink's level. It returns `False` and changes nothing when no log file is being written.

The next `init_logger` call, including `TidalPy.reinit()`, restores the configured sink levels and resets the global threshold to `trace`.

### `log_message(level, message)` and Level Helpers

`log_trace`, `log_debug`, `log_info`, `log_warning`, `log_error`, and `log_critical` take just the message; `log_message` takes the level first, by name or integer (`"off"` emits nothing). TidalPy's Cython and Python code use these, so their output interleaves correctly with the C++ messages.

### `flush_logger()` and `shutdown_logger()`

Warnings, errors, and critical messages are flushed as they are written, so they reach the log file even if the process later crashes. Trace, debug, and info lines are buffered: a log file read while the interpreter runs (in a Jupyter notebook, say) may lack them until `flush_logger` runs. A normal interpreter exit flushes them; a crash loses them.

`shutdown_logger` flushes and turns the logger off for every extension; a later `init_logger` (or `set_log_level`) turns it back on. It is not needed at interpreter exit.

### Notebooks

A Jupyter kernel does not show what compiled code writes to stdout on every platform (a Windows kernel never does). In a notebook the console sink therefore holds its messages (`console_pending`), and TidalPy prints them below the cell once it finishes. By default only warnings and above print there (`[logging] notebook_console_level`); `print_log_notebook = true` prints everything at `console_level`. A notebook message prints as `[TidalPy] [warning] <text>`, without the timestamp of the console and log file, and ends with an LF on every platform, so re-running a notebook leaves its saved output unchanged.

## Configuration

The `[logging]` keys, their defaults, and the log file location are listed on the [Configurations](../Overview/2_TidalPy_Configurations.md#logging) page. The log file is named `TidalPy_<YYYYMMDD-HHMMSS>.log`, and test mode (the `TIDALPY_TEST_MODE` environment variable) never writes one. A test that asserts on log output reads a log file, as `capture_log` and the `spdlog_text` fixture of `Tests/conftest.py` do.

## C++ API

```cpp
#include "TidalPy/Utilities/logging/logger_.hpp"

TIDALPY_LOG_DEBUG("Initializing layer: {}", layer_name);
TIDALPY_LOG_WARN("EOS data not populated for layer {}", index);
TIDALPY_LOG_ERROR("Failed to open binary file: {}", path);
```

The macros are `TIDALPY_LOG_TRACE`, `TIDALPY_LOG_DEBUG`, `TIDALPY_LOG_INFO`, `TIDALPY_LOG_WARN`, `TIDALPY_LOG_ERROR`, and `TIDALPY_LOG_CRITICAL`. Each does nothing before the logger exists, so they are safe anywhere in startup, and cheap enough to leave at debug level inside a solver loop. Format strings use spdlog's Python-style braces. Logging from several C++ threads, and reconfiguring while they log, is safe: each message goes wholly to the old sinks or wholly to the new ones.

## Sharing the Logger Pointer Across Extensions

Every Cython extension that logs from C++ must wire itself to the shared logger at module-init level, outside any function:

```cython
from TidalPy.Utilities.logging.logger cimport (
    set_tidalpy_logger_ptr_void, get_tidalpy_logger_address)
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
```

On Windows each compiled extension is a separate DLL with its own copy of the header-only logger pointer, and this call connects it to the shared logger. On Linux and macOS the pointer is already shared and the call does nothing, so the same source works everywhere. The `cimport` also makes Cython import the logger module first, so the logger exists before its address is taken.

spdlog [v1.15.3](https://github.com/gabime/spdlog/releases/tag/v1.15.3) is header only and ships as a git submodule at `Dependencies/spdlog`.
