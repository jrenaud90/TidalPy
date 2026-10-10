# Binary Serialization (`Utilities.binary`)

_Updated: 2026-10-02_

TidalPy writes worlds, layers, systems, and physics models to a compact binary format. A TOML configuration is the readable, editable way to describe a world; the binary format saves and restores an object as it stands, with every attached sub-model, without going back through the builders. This page describes the format and how a class takes part in it; the other pages only list each class's id.

Every file starts with a fixed 20-byte header naming the format version, the class that wrote it, and the payload size, so a reader can tell what a file holds before reading it.

## Python API

```python
from TidalPy.Structures import TerrestrialWorld, load_world
from TidalPy.Utilities.binary import binary_file_class, check_binary_file, get_current_schema_version

TerrestrialWorld("world", 1.0e6, 1.0e22).save_binary("world.tpyb")   # A world with no layers yet
info = check_binary_file("world.tpyb")
# {'schema_major': 0, 'schema_minor': 2, 'schema_patch': 0,
#  'schema_version': '0.2.0', 'class_id': 201, 'payload_size': 148}

binary_file_class("world.tpyb")   # 'TerrestrialWorld'
world = load_world("world.tpyb")  # A new TerrestrialWorld, no placeholder object needed
get_current_schema_version()      # '0.2.0'
```

- `check_binary_file(path)` reads only the header, so it is cheap and safe on a file of unknown origin. It raises `FileNotFoundError` for a missing path and `IOError` for wrong magic bytes, a foreign byte order, or a file shorter than the header. It does not check the payload size; `load_binary` does.
- `binary_file_class(path)` reads only the header and returns the name of the class that wrote the file, or `None` for a file without the TidalPy magic bytes (a TOML file, say). `TidalPy.Structures.load_world(path)` and `load_system(path)` use it to load a file into a new object of the right class, and `build_world` and `build_system` to accept a binary file in place of TOML.
- `get_current_schema_version()` returns the version this build writes.

Every object also has `save_binary(path)` and `load_binary(path, force=False)` (see [Base Classes](classes.md)).

## Saving

`save_binary` writes to a temporary file beside the target (`<target>.<random>.partial`) and renames it over the target only once the whole record is written and closed. A failed save (a full disk, an object that cannot be written, a target that cannot be replaced) raises `IOError`, removes the temporary file, and leaves any previous file unchanged. On Windows, a target another program holds open cannot be replaced, so the save raises until that program closes it. The saved file takes the directory's default permissions, and a target that is a symbolic link is replaced by the new file, not written through.

## Integrity Checks

`load_binary` refuses a file, raising `IOError`, when:

- The root record's class id is not that of the object being loaded into. The message names both classes ("it is a GasGiantWorld file, not a TerrestrialWorld one").
- The magic bytes are wrong, or the byte order is not this machine's.
- The schema version is incompatible and `force=True` was not given (see [Schema Version](#schema-version)).
- Any record's header claims a payload larger than what is left in the file.
- A nested record's class id is not the one its owner expects.
- A record's payload ends before this build has read its fields, or bytes are left in it after them: it was written with a different layout, or it is corrupt.
- A count read from the file (a string length, a table length, a number of layers or worlds) is larger than what is left in the record.
- Bytes are left over after the root record.

A load that raises leaves every saved setting of the object as it was. The header checks run before anything is read into the object, and the object's own record is restored if a later check fails. State the object does not save, such as a solved structure, also survives a failed load for a class that reads the file into a scratch object first (`make_binary_scratch` in C++).

## Schema Version

The current schema version is `0.2.0`.

| Condition | Result |
|---|---|
| Same major and minor, any patch | Compatible. A differing patch logs an informational line. |
| Different major or minor | Incompatible. Reading raises unless `force=True` is given. |

A minor bump can change a class's member layout, and reading an old payload into a new layout gives an object that looks valid but is not. Use `force=True` only when the layout is known not to have changed; it still warns, and it relaxes only the version check, never the other integrity checks.

## File Format

| Offset | Size | Field | Description |
|---|---|---|---|
| 0 | 4 | `magic` | The ASCII bytes `TPYB`. |
| 4 | 1 | `schema_major` | Schema major version. |
| 5 | 1 | `schema_minor` | Schema minor version. |
| 6 | 1 | `schema_patch` | Schema patch version. |
| 7 | 1 | `byte_order` | The writer's byte order: `0` little-endian, `1` big-endian. |
| 8 | 4 | `class_id` | Class type id, `uint32_t` in the writer's byte order. |
| 12 | 8 | `payload_size` | Every byte of the record after the header, nested sub-object records included, `uint64_t` in the writer's byte order. |

Fields are written one at a time with fixed-width types, never as a packed struct, so the layout does not depend on the compiler. A reader on a machine of the other byte order refuses the file rather than converting it. Every platform TidalPy supports (Windows, Linux, and macOS on x64 and ARM64) is little-endian, so files move freely between them.

A record's payload holds its own fields and the complete records of the sub-objects it owns, so a record spans exactly its header plus `payload_size` bytes. Model, layer, and material names are written as a `uint32_t` length followed by the raw UTF-8 bytes:

```text
[uint32_t length][length bytes of UTF-8 text]
```

A physics model writes its model name, the number of parameters (`uint32_t`), and each parameter by key: the key, a kind byte (`c_ParamKind`), a value count (`uint64_t`), and the values (doubles). Reading validates as construction does, and a key the record lacks reads at its default, so adding a parameter never breaks older records. State outside the parameter table follows (the isotope radiogenics model writes its label count, `uint64_t`, then each label), and a composite (a phase or a material) then writes one optional record per slot.

```text
[string model_name][uint32_t num_params]
    num_params x [string key][uint8_t kind][uint64_t num_values][num_values x double]
[extra state, if any]
```

## Nested and Recursive Serialization

A layer owns its physics models, a world its layers, and a system its worlds, and all round-trip together. An optional owned sub-object is written as a one-byte presence flag (zero absent, one present), followed when present by the sub-object's complete record, header and all. Both belong to the owner's payload, so a file is one root record with every sub-object record nested inside it. On read, the owner reads the flag and, when set, calls its family's dispatch factory, which peeks at the record's class id, constructs the matching subclass, and has it read the record.

Each family has its factory: `c_rheology_from_binary`, `c_viscosity_from_binary`, `c_cooling_from_binary`, `c_radiogenics_from_binary`, `c_luminosity_from_binary`, `c_tide_from_binary`, and `c_layer_from_binary`; `c_melting_curve_from_binary`, `c_melt_weakening_from_binary`, `c_bulk_modulus_mixing_from_binary`, and `c_bulk_viscosity_mixing_from_binary` (`PartialMelt`); and `c_eos_from_binary`, `c_shear_modulus_from_binary`, `c_phase_from_binary`, and `c_material_from_binary` (`Material`).

```cpp
// Writing, in c_Layer::p_write_payload.
write_optional_binary(out, this->p_shear_rheology);

// Reading, in c_Layer::p_read_payload.
this->p_shear_rheology = read_optional_binary<c_RheologyBase>(in, force, c_rheology_from_binary);
```

| Class | Recursively serialized sub-objects |
|---|---|
| `c_Layer` | material, shear rheology override, bulk rheology override, cooling, radiogenics (after its geometry, switches, and flags) |
| `c_Material` | solid phase, liquid phase, solidus, liquidus, melt weakening, bulk-modulus mixing, bulk-viscosity mixing (after its parameters) |
| `c_Phase` | equation of state, shear modulus, shear viscosity, bulk viscosity, shear rheology, bulk rheology (after its thermal parameters) |
| `c_BaseWorld`, `c_TerrestrialWorld`, `c_GasGiantWorld` | tide model (after the tide configuration scalars), then every layer in index order (the spin model's moment-of-inertia factor and the pinned solver settings are scalars in the payload) |
| `c_StarWorld` | the `c_BaseWorld` sub-objects, then the luminosity model |

> [!NOTE]
> The equation-of-state profile is never saved, because it is derived from the layers' materials. Run `solve_eos` on a loaded world to regenerate it.

## Class Type IDs

Each concrete class has a unique id so the dispatch factories can rebuild the right subclass. The ranges are grouped by family, leaving room to add models without renumbering. ID zero is `Unknown` and is never written.

| Range | Family | Members |
|---|---|---|
| 1-3 | Base classes | `TidalPyBase` 1, `StructureBase` 2, `PhysicsBase` 3 |
| 100-199 | Layers | `Layer` 100 |
| 200-299 | Worlds and systems | `BaseWorld` 200, `TerrestrialWorld` 201, `GasGiantWorld` 202, `StarWorld` 203, `System` 210 |
| 300-399 | Rheology | `RheologyBase` 300, `Elastic` 301, `Viscous` 302, `Voigt` 303, `Maxwell` 304, `Burgers` 305, `Andrade` 306, `Sundberg` 307, `Zener` 308, `SeismicQ` 309 |
| 400-499 | Cooling | `CoolingBase` 400, `OffCooling` 401, `ConvectiveCooling` 402, `ConductiveCooling` 403 |
| 500-599 | Radiogenics | `RadiogenicsBase` 500, `OffRadiogenics` 501, `IsotopeRadiogenics` 502, `FixedRadiogenics` 503 |
| 600-699 | Materials | Equation-of-state laws: `ConstantEOSLaw` 610, `BirchMurnaghanEOSLaw` 611, `VinetEOSLaw` 612, `MurnaghanEOSLaw` 613, `PolytropeEOSLaw` 614, `ModifiedPolytropeEOSLaw` 615, `InterpolatedEOSLaw` 616; shear-modulus laws: `ConstantShearModulus` 620, `LinearShearModulus` 621, `InterpolatedShearModulus` 622; `Phase` 630, `Material` 631 |
| 700-799 | Melting laws | Melting curves: `ConstantMeltingCurve` 710, `SimonGlatzelCurve` 711, `SimonGlatzel2Curve` 712, `InterpolatedMeltingCurve` 713; weakening: `NoMeltWeakening` 720, `SpohnMeltWeakening` 721, `HenningMeltWeakening` 722; `HashinShtrikmanMixing` 730; `CompactionViscosity` 740 |
| 800-899 | Viscosity | `ViscosityBase` 800, `ArrheniusViscosity` 801, `ReferenceViscosity` 802, `ConstantViscosity` 803, `InterpolatedViscosity` 804, `CompositeViscosity` 805 |
| 900-999 | Tide models | `TideBase` 900, `RheologyTide` 901, `FixedQTide` 902, `FixedLagTide` 903, `CTLQTide` 904 |
| 1000-1099 | Luminosity | `LuminosityBase` 1000, `FixedLuminosity` 1001, `MassToLuminosity` 1002, `PowerLawLuminosity` 1003 |

## C++ API

`c_TidalPyBaseClass` (`tidalpy_base_.hpp`) holds the only `write_binary` and `read_binary` in the package. A class takes part through a few overrides:

| Member | Role |
|---|---|
| `write_binary(out)` | Writes the header, with the measured payload size, then the payload from `p_write_payload`. |
| `read_binary(in, force)` | Reads and checks the header and class id, passes `p_read_payload` a stream holding only the payload, and raises unless it read every byte. |
| `get_binary_class_id()` | The `BinaryClassID` of the class's records. Every concrete class returns its own. |
| `write_binary_bytes()` | The whole record as a byte string, the in-memory form of `save_binary`; a world's and a system's `copy()` and pickle use it. |
| `load_binary_bytes(record_bytes, source, force)` | `load_binary` from a record in memory; `source` names it in error messages. |
| `p_write_payload(out)`, `p_read_payload(in, force)` | The class's payload. An override calls its parent's first, then writes or reads its own fields and owned sub-object records. |
| `make_binary_scratch()` | Optional: a new object of the class that `load_binary` reads a file into first (see [Integrity Checks](#integrity-checks)). |

A parameter-table model (`c_SpecModel`) needs only its table and its `C_CLASS_ID`, from which it gets `get_binary_class_id` and its payload; state outside the table goes through `p_write_extra` and `p_read_extra`.

Include `binary_.hpp` (header only) in code that reads or writes these files directly; it routes version-mismatch warnings through the shared logger.

| Function | Description |
|---|---|
| `write_binary_header(out, class_id, payload_size)` | Write the 20-byte header. |
| `read_binary_header(in)` | Read the 20-byte header, validating the magic bytes and byte order. |
| `read_binary_header_from_file(path)` | Open a path and read its header. |
| `c_binary_file_class_id(path)` | The class id of a file's record, or `BinaryClassID::Unknown` (0) for a file without the magic bytes. |
| `c_peek_binary_header(in)` | Read the header at the read position and rewind. |
| `check_binary_schema_version(header, force)` | Validate the version, logging a warning on mismatch. |
| `c_read_binary_record_header(in, force)` | Read and fully validate a record header, refusing a payload larger than the rest of the stream. |
| `binary_bytes_remaining(in)` | The bytes left in a seekable stream after its read position. |
| `check_binary_count(in, count, element_bytes, what)` | Throw when a count read from a file needs more bytes than are left. |
| `c_host_binary_byte_order()` | The `byte_order` value this machine writes. |
| `c_binary_class_name(class_id)` | The Python class name of a class id (`"TerrestrialWorld"`). |
| `c_binary_temporary_path(target)` | The temporary sibling path `save_binary` writes before renaming. |
| `write_binary_string(out, text)` and `read_binary_string(in)` | Length-prefixed string input and output. |
| `write_optional_binary(out, pointer)`, `read_optional_binary<T>(in, force, factory)` | Write or read an optional owned sub-object (see [Nested and Recursive Serialization](#nested-and-recursive-serialization)). |

| Constant | Value |
|---|---|
| `TIDALPY_SCHEMA_MAJOR`, `TIDALPY_SCHEMA_MINOR`, `TIDALPY_SCHEMA_PATCH` | `0`, `2`, `0` |
| `TIDALPY_BINARY_MAGIC` | `"TPYB"` |
| `TIDALPY_BINARY_HEADER_BYTES` | `20` |
| `TIDALPY_BINARY_LITTLE_ENDIAN`, `TIDALPY_BINARY_BIG_ENDIAN` | `0`, `1` |

## Portability Notes

File paths are passed as UTF-8 `std::string`. Non-ASCII paths on Windows are not guaranteed to work everywhere, so prefer ASCII paths when portability matters. The only external dependency is [spdlog](https://github.com/gabime/spdlog/releases/tag/v1.15.3), for the version-mismatch warnings.
