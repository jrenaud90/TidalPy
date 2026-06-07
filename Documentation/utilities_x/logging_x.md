# Utilities_x: Logging (`logging_x`)

TidalPy uses [spdlog](https://github.com/gabime/spdlog) for C++ logging, with a thin
Cython wrapper that lets Python code initialize and control the same logger.

---

## Overview

All C++ code in TidalPy logs through a single named spdlog logger (`"TidalPy"`).
The logger is initialized once at package startup and is accessible from both
C++ (via macros) and Python (via `TidalPy.Utilities_x.logging_x`).

```
TidalPy.__init__
    └─ init_logger(config["logging"])     # Python: logger.pyx
           └─ cy_init_logger(config)      # C++: logger_.hpp
                  └─ spdlog registry → stdout (+ optional file)
```

---

## Python API

```python
from TidalPy.Utilities_x.logging_x import init_logger, set_log_level, shutdown_logger

# Called automatically at package startup; can be called again to reconfigure.
init_logger({
    "console_level": "info",   # or "trace", "debug", "warning", "error", "critical", "off"
    "file_level":    "debug",
    "log_to_file":   True,
    "log_file_path": "/tmp/tidalpy.log",
})

set_log_level("debug")   # adjust verbosity at runtime

shutdown_logger()        # flush + release file handles (called via atexit)
```

### `init_logger(config: dict = None)`

Initialize the TidalPy C++ logger. Safe to call multiple times; only the first
call has any effect.

| Config key | Type | Default | Description |
|---|---|---|---|
| `console_level` | str or int | `"info"` | Log level for console output |
| `file_level` | str or int | `"info"` | Log level for file output |
| `log_to_file` | bool | `False` | Enable file logging |
| `log_file_path` | str | `""` | Absolute path to log file (UTF-8) |

Level names (case-insensitive): `trace`, `debug`, `info`, `warning`/`warn`,
`error`, `critical`, `off`. Integer values 0–6 are also accepted.

### `set_log_level(level: str | int)`

Adjust the active log level at runtime. Updates both the logger and all sinks.

### `shutdown_logger()`

Flush all pending log messages and release file handles. Called automatically
on interpreter exit.

---

## C++ API

Include `logger_.hpp` and use the macros directly:

```cpp
#include "TidalPy/Utilities_x/logging_x/logger_.hpp"

// Inside any C++ function:
TIDALPY_LOG_DEBUG("Initializing layer: {}", layer_name);
TIDALPY_LOG_WARN("EOS data not populated for layer {}", index);
TIDALPY_LOG_ERROR("Failed to open binary file: {}", path);
```

Macros available: `TIDALPY_LOG_TRACE`, `TIDALPY_LOG_DEBUG`, `TIDALPY_LOG_INFO`,
`TIDALPY_LOG_WARN`, `TIDALPY_LOG_ERROR`, `TIDALPY_LOG_CRITICAL`.

All macros check for a nullptr logger and are no-ops if the logger has not been
initialized — safe to call from any translation unit.

spdlog supports Python-style `{}` format strings via the bundled `fmtlib`.

---

## Cross-Platform Notes — Pointer Sharing Pattern

`logger.pyx` creates the spdlog logger at import time and stores a raw pointer
to it (`tidalpy_logger_ptr`). The `TIDALPY_LOG_*` macros use that pointer
directly rather than calling `spdlog::get()` on every invocation.

This mirrors the `constants.pyx` / `constants_.hpp` pattern used for TidalPy's
global config pointer. Any new Cython extension that uses C++ logging **must**
add the following two lines at its module-init level (the `.pyx` file, outside
any function):

```cython
from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void, get_tidalpy_logger_address)
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
```

`get_tidalpy_logger_address()` is a `cdef api` function (the same mechanism as
`get_shared_config_address()` in `constants.pyx`).  When another module
`cimport`s it, Cython automatically imports `TidalPy.Utilities_x.logging_x.logger`
first, which guarantees the logger exists before the pointer is fetched.

- **Linux/macOS**: The inline variable is process-wide; the pointer is already
  shared automatically.  The set call is a harmless no-op.
- **Windows**: Each `.pyd` DLL holds its own copy of `tidalpy_logger_ptr`.
  The set call is required to wire the DLL into the shared logger.

---

## Configuration in `TidalPy_Configs.toml`

```toml
[logging]
console_level = "info"
file_level    = "info"
log_to_file   = false
log_file_path = ""
```

---

## Dependencies

- [spdlog v1.15.3](https://github.com/gabime/spdlog/releases/tag/v1.15.3) —
  git submodule at `Dependencies/spdlog`. Header-only; no separate compilation
  step required.
