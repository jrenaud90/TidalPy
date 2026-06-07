# Utilities_x: Binary Serialization (`binary_x`)

TidalPy uses a custom binary format for fast, compact storage of world and layer
configurations. Every TidalPy binary file begins with a fixed 20-byte header that
identifies the format version, the serialized class type, and the payload size.

---

## File Format

| Offset | Size | Field | Description |
|--------|------|-------|-------------|
| 0 | 4 | `magic` | ASCII bytes `TPYB` |
| 4 | 1 | `schema_major` | Schema major version |
| 5 | 1 | `schema_minor` | Schema minor version |
| 6 | 1 | `schema_patch` | Schema patch version |
| 7 | 1 | `reserved` | Always `0` |
| 8 | 4 | `class_id` | Class type ID (`uint32_t`, host byte order) |
| 12 | 8 | `payload_size` | Bytes of payload after this header (`uint64_t`, host byte order) |

**Total header size:** 20 bytes.

**Byte order:** Host (native) byte order. All supported TidalPy platforms
(Windows/Linux/macOS on x64 and ARM64) are little-endian, so files are
cross-platform in practice.

Fields are written individually with explicit `stream.write` calls, no implicit
struct padding, so the byte layout is identical on all platforms.

---

## Schema Version

**Current schema version: `0.2.0`**.

Compatibility rules:

| Condition | Result |
|-----------|--------|
| Same `major.minor`, any `patch` | Compatible; informational log if `patch` differs |
| Different `major` or `minor` | Incompatible; `read_binary` raises unless `force=True` |

---

## Python API

```python
from TidalPy.Utilities_x.binary_x import check_binary_file, get_current_schema_version

# Inspect a binary file without loading it:
info = check_binary_file("/path/to/world.tpyb")
# {
#   "schema_major": 0,
#   "schema_minor": 2,
#   "schema_patch": 0,
#   "schema_version": "0.2.0",
#   "class_id": 201,       # BinaryClassID::LayeredWorld = 201
#   "payload_size": 4096,
# }

ver = get_current_schema_version()  # "0.2.0"
```

### `check_binary_file(path: str) -> dict`

Read the header of a TidalPy binary file and return its contents as a dict.

Raises `FileNotFoundError` if the file does not exist. Raises `IOError` if the
magic bytes are wrong or the file is too short (less than 20 bytes).

### `get_current_schema_version() -> str`

Return the schema version compiled into this TidalPy build (`"0.2.0"`).

---

## C++ API

Include `binary_.hpp` in any C++ code that reads or writes binary files.
This header also includes `logger_.hpp` so that version-mismatch warnings
are emitted through the TidalPy logger.

```cpp
#include "binary_.hpp"

// Write a binary file:
std::ofstream out("world.tpyb", std::ios::binary);
tidalpy::write_binary_header(
    out,
    static_cast<uint32_t>(tidalpy::BinaryClassID::LayeredWorld),
    payload_size
);
// ... write payload bytes ...

// Read a binary file:
std::ifstream in("world.tpyb", std::ios::binary);
tidalpy::c_BinaryHeader header = tidalpy::read_binary_header(in);
if (!tidalpy::check_binary_schema_version(header)) {
    throw std::runtime_error("Incompatible binary schema version");
}
// ... read payload bytes ...

// Convenience: read header from a path string:
tidalpy::c_BinaryHeader header = tidalpy::read_binary_header_from_file(path);
```

### Functions

| Function | Description |
|----------|-------------|
| `write_binary_header(out, class_id, payload_size)` | Write 20-byte header to `std::ostream&` |
| `read_binary_header(in)` | Read 20-byte header from `std::istream&` |
| `read_binary_header_from_file(path)` | Open file by path and read header |
| `check_binary_schema_version(header, force)` | Validate version; logs warning on mismatch |

### Constants

| Name | Value | Description |
|------|-------|-------------|
| `TIDALPY_SCHEMA_MAJOR` | `0` | Schema major version |
| `TIDALPY_SCHEMA_MINOR` | `2` | Schema minor version |
| `TIDALPY_SCHEMA_PATCH` | `0` | Schema patch version |
| `TIDALPY_BINARY_MAGIC` | `"TPYB"` | File magic bytes |
| `TIDALPY_BINARY_HEADER_BYTES` | `20` | Header size in bytes |

---

## Class Type IDs (`BinaryClassID`)

| Class | ID |
|-------|----|
| `TidalPyBase` | 1 |
| `StructureBase` | 2 |
| `PhysicsBase` | 3 |
| `BaseLayer` | 100 |
| `PhysicsLayer` | 101 |
| `SolidLiquidLayer` | 102 |
| `GasLayer` | 103 |
| `BaseWorld` | 200 |
| `LayeredWorld` | 201 |
| `GasGiantWorld` | 202 |
| `StarWorld` | 203 |
| `RheologyBase` | 300 |
| `CoolingBase` | 400 |
| `RadiogenicsBase` | 500 |

---

## Cross-Platform Notes

- All field sizes use `<cstdint>` fixed-width types (`uint8_t`, `uint32_t`, `uint64_t`).
- Fields are written one at a time (`stream.write`), not as a packed struct, so
  compiler padding never affects the layout.
- File paths are passed as UTF-8-encoded `std::string`. Non-ASCII paths on Windows
  may not work on all systems; prefer ASCII paths for maximum portability.
- The header-only nature of `binary_.hpp` means no separate compilation step is needed.

---

## Dependencies

- [spdlog v1.15.3](https://github.com/gabime/spdlog/releases/tag/v1.15.3) via `logger_.hpp`
  (for version-mismatch warning messages).
