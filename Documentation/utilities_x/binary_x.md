# Binary Serialization (`Utilities_x.binary_x`)

_Updated: 2026-09-12_

TidalPy writes worlds, layers, systems, and physics models to a compact binary format. The purpose is not archival: a TOML configuration is the readable, editable, portable way to describe a world, and it is what you should commit to a repository. The binary format exists for the cases TOML handles badly, chiefly saving and restoring an object graph exactly as it stands, including every attached sub-model, without going back through the builders.

Every file starts with a fixed 20-byte header naming the format version, the class that wrote it, and the payload size. That is what makes a file self-describing: a reader can identify what it is holding before deciding whether it can read it.

## File format

| Offset | Size | Field | Description |
|---|---|---|---|
| 0 | 4 | `magic` | The ASCII bytes `TPYB`. |
| 4 | 1 | `schema_major` | Schema major version. |
| 5 | 1 | `schema_minor` | Schema minor version. |
| 6 | 1 | `schema_patch` | Schema patch version. |
| 7 | 1 | `reserved` | Always zero. |
| 8 | 4 | `class_id` | Class type id, `uint32_t` in host byte order. |
| 12 | 8 | `payload_size` | Payload bytes following the header, `uint64_t` in host byte order. |

Fields are written one at a time with explicit stream writes rather than as a packed struct, so compiler padding never affects the layout. The byte order is the host's. Every platform TidalPy supports, Windows, Linux, and macOS on x64 and ARM64, is little-endian, so files move between them in practice, but the format does not promise it.

## Schema version

The current schema version is `0.2.0`.

| Condition | Result |
|---|---|
| Same major and minor, any patch | Compatible. A differing patch produces an informational log line. |
| Different major or minor | Incompatible. Reading raises unless `force=True` is given. |

The rule is strict for a reason. A minor version bump can change a class's member layout, and reading an old payload into a new layout produces an object that looks valid and is not. `force=True` exists for the case where you know the layout did not change, and it still warns.

## Python API

```python
from TidalPy.Utilities_x.binary_x import check_binary_file, get_current_schema_version

info = check_binary_file("world.tpyb")
# {'schema_major': 0, 'schema_minor': 2, 'schema_patch': 0,
#  'schema_version': '0.2.0', 'class_id': 201, 'payload_size': 87}

get_current_schema_version()   # '0.2.0'
```

`check_binary_file(path)` reads only the header, so it is cheap and safe on a file of unknown provenance. It raises `FileNotFoundError` if the path does not exist and `IOError` if the magic bytes are wrong or the file is shorter than the header.

`get_current_schema_version()` returns the version compiled into this build, which is what a file would be written with now.

## C++ API

Include `binary_.hpp` in any code that reads or writes these files. It pulls in `logger_.hpp` so version-mismatch warnings go through the shared logger.

```cpp
#include "binary_.hpp"

// Writing.
std::ofstream out("world.tpyb", std::ios::binary);
tidalpy::write_binary_header(
    out,
    static_cast<uint32_t>(tidalpy::BinaryClassID::LayeredWorld),
    payload_size);
// ... payload bytes ...

// Reading.
std::ifstream in("world.tpyb", std::ios::binary);
const tidalpy::c_BinaryHeader header = tidalpy::read_binary_header(in);
if (!tidalpy::check_binary_schema_version(header)) {
    throw std::runtime_error("Incompatible binary schema version");
}
```

| Function | Description |
|---|---|
| `write_binary_header(out, class_id, payload_size)` | Write the 20-byte header. |
| `read_binary_header(in)` | Read the 20-byte header. |
| `read_binary_header_from_file(path)` | Open a path and read its header. |
| `check_binary_schema_version(header, force)` | Validate the version, logging a warning on mismatch. |
| `write_binary_string(out, text)` and `read_binary_string(in)` | Length-prefixed string input and output. |
| `binary_string_bytes(text)` | The payload bytes a length-prefixed string contributes, for sizing a header. |
| `write_optional_binary(out, unique_ptr)` | Write an optional owned sub-object: a presence flag, then its record if present. |
| `read_optional_binary<T>(in, force, factory)` | Read an optional sub-object, rebuilding it through the given factory. |
| `optional_binary_flag_bytes()` | The bytes one presence flag contributes, which is one. |

| Constant | Value |
|---|---|
| `TIDALPY_SCHEMA_MAJOR`, `TIDALPY_SCHEMA_MINOR`, `TIDALPY_SCHEMA_PATCH` | `0`, `2`, `0` |
| `TIDALPY_BINARY_MAGIC` | `"TPYB"` |
| `TIDALPY_BINARY_HEADER_BYTES` | `20` |

## Variable-length strings

Model names, layer names, and material names are written as a `uint32_t` length followed by the raw UTF-8 bytes:

```
[uint32_t length][length bytes of UTF-8 text]
```

Use `binary_string_bytes(text)` when computing a record's payload size, so the header and the payload cannot disagree.

## Nested and recursive serialization

Containers own sub-objects that have to round-trip with them: a layer owns its physics models, a world owns its layers, a system owns its worlds. The encoding is uniform at every level.

An optional owned sub-object, held in a `unique_ptr`, is written as a one-byte presence flag, zero for absent and one for present. When present, the sub-object's own complete record follows immediately, header and all. The presence flag counts toward the owning record's payload size, while the nested record is a separate self-describing record appended after it, so a file is a sequence of concatenated records.

On read, the owning class reads the flag and, when set, calls a binary-dispatch factory. The factory peeks the upcoming record's class id, default-constructs the matching concrete subclass, and delegates to its `read_binary`. That is what allows a layer to restore an Andrade rheology it never knew it had.

| Module | Dispatch factory |
|---|---|
| `rheology_x` | `c_rheology_from_binary` |
| `viscosity_x` | `c_viscosity_from_binary` |
| `partial_melt_x` | `c_partial_melt_from_binary` |
| `cooling_x` | `c_cooling_from_binary` |
| `radiogenics_x` | `c_radiogenics_from_binary` |
| `Material_x.eos` | `c_material_eos_from_binary` |
| `stellar_x` | `c_luminosity_from_binary` |
| `Tides_x` | `c_tide_from_binary` |
| `structures_x.layers` | `c_layer_from_binary` |

```cpp
// Writing, as c_PhysicsLayer does it.
write_optional_binary(out, this->p_shear_rheology);
write_optional_binary(out, this->p_bulk_rheology);

// Reading.
this->p_shear_rheology =
    read_optional_binary<c_RheologyBase>(in, force, c_rheology_from_binary);
this->p_bulk_rheology =
    read_optional_binary<c_RheologyBase>(in, force, c_rheology_from_binary);
```

What each layer class carries recursively, after its own scalar payload:

| Layer | Recursively serialized sub-objects |
|---|---|
| `c_BaseLayer` | material EOS model |
| `c_PhysicsLayer` | material EOS model, shear rheology, bulk rheology, shear viscosity, bulk viscosity, partial melt |
| `c_GasLayer` | the same six, inherited |
| `c_SolidLiquidLayer` | the same six, plus cooling and radiogenics |

> [!NOTE]
> The equation-of-state profile data is never serialized, because it is derived from the attached model: `solve_eos` runs directly on a loaded world and regenerates it.

## Class type ids

Each concrete class needs a unique id so the dispatch factories can reconstruct the right subclass. The ranges are grouped by family, which leaves room to add models without renumbering.

| Range | Family | Members |
|---|---|---|
| 1-3 | Base classes | `TidalPyBase` 1, `StructureBase` 2, `PhysicsBase` 3 |
| 100-199 | Layers | `BaseLayer` 100, `PhysicsLayer` 101, `SolidLiquidLayer` 102, `GasLayer` 103 |
| 200-299 | Worlds and systems | `BaseWorld` 200, `LayeredWorld` 201, `GasGiantWorld` 202, `StarWorld` 203, `System` 210 |
| 300-399 | Rheology | `RheologyBase` 300, `Elastic` 301, `Viscous` 302, `Voigt` 303, `Maxwell` 304, `Burgers` 305, `Andrade` 306, `Sundberg` 307 |
| 400-499 | Cooling | `CoolingBase` 400, `OffCooling` 401, `ConvectiveCooling` 402, `ConductiveCooling` 403 |
| 500-599 | Radiogenics | `RadiogenicsBase` 500, `OffRadiogenics` 501, `IsotopeRadiogenics` 502, `FixedRadiogenics` 503 |
| 600-699 | Material EOS | `MaterialEOSBase` 600, `ConstantDensityEOS` 601, `BirchMurnaghanEOS` 602, `VinetEOS` 603, `InterpolatedEOS` 604 |
| 700-799 | Partial melt | `PartialMeltBase` 700, `OffPartialMelt` 701, `SpohnPartialMelt` 702, `HenningPartialMelt` 703 |
| 800-899 | Viscosity | `ViscosityBase` 800, `ArrheniusViscosity` 801, `ReferenceViscosity` 802, `ConstantViscosity` 803 |
| 900-999 | Tide models | `TideBase` 900, `RheologyTide` 901, `FixedQTide` 902, `FixedLagTide` 903, `CTLQTide` 904 |
| 1000-1099 | Luminosity | `LuminosityBase` 1000, `FixedLuminosity` 1001, `MassToLuminosity` 1002, `PowerLawLuminosity` 1003 |

Id zero is `Unknown` and is never written.

## Portability notes

Every field uses a fixed-width type from `<cstdint>`, and fields are written individually rather than as a struct, so the layout does not depend on the compiler. File paths are passed as UTF-8 `std::string`; non-ASCII paths on Windows are not guaranteed to work everywhere, so prefer ASCII paths when portability matters. The header is header-only, so no separate compilation step is involved.

The only external dependency is [spdlog](https://github.com/gabime/spdlog/releases/tag/v1.15.3), reached through `logger_.hpp`, and only for the version-mismatch warnings.
