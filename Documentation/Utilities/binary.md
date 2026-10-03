# Binary Serialization (`Utilities.binary`)

_Updated: 2026-10-02_

TidalPy writes worlds, layers, systems, and physics models to a compact binary format. This page describes that format and how a class takes part in it. The other pages only list each class's id. A TOML configuration is the readable, editable way to describe a world. The binary format saves and restores an object structure as it stands, including every attached sub-model, without going back through the builders.

Every file starts with a fixed 20-byte header naming the format version, the class that wrote it, and the payload size, so a reader can identify what a file holds before deciding whether it can read it.

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

Fields are written one at a time with explicit stream writes rather than as a packed struct, so compiler padding never affects the layout. Multi-byte fields use the writer's byte order, which the header records. A reader on a machine of the other byte order refuses the file instead of converting it. Every platform TidalPy supports (Windows, Linux, and macOS on x64 and ARM64) is little-endian, so files move freely between them.

A record's payload holds its own fields and the complete records of the sub-objects it owns (see [Nested and Recursive Serialization](#nested-and-recursive-serialization)), so a record spans exactly its header plus `payload_size` bytes.

## Schema Version

The current schema version is `0.2.0`.

| Condition | Result |
|---|---|
| Same major and minor, any patch | Compatible. A differing patch produces an informational log line. |
| Different major or minor | Incompatible. Reading raises unless `force=True` is given. |

A minor version bump can change a class's member layout, and reading an old payload into a new layout produces an object that looks valid but is not. `force=True` is for the case where the layout is known not to have changed, and it still warns. It relaxes only the version check, never the integrity checks below.

## Integrity Checks

`load_binary` refuses a file, raising `IOError`, when:

- The root record's class id is not the class id of the object being loaded into. The message names both classes, as in "it is a GasGiantWorld file, not a TerrestrialWorld one".
- The magic bytes are wrong, or the byte order is not this machine's.
- The schema version is incompatible and `force=True` was not given.
- Any record's header claims a payload larger than what is left in the file.
- A nested record's class id is not the one its owner expects.
- A record's payload ends before this build has read its fields, or bytes are left in it after them. The record was written with a different layout, which would otherwise misalign every record after it, or it is corrupt.
- A count read from the file (a string length, a table length, a number of layers or worlds) is larger than what is left in the record.
- Bytes are left over after the root record. The file was written with a layout this build does not read, or it is corrupt.

A load that raises an error leaves every setting the object had before the call. The whole file is read into memory first, and the initial checks are performed before anything is read into the object. For the other checks, which can only run during or after the record is read, `load_binary` keeps a copy of the object's own record in memory and reads it back into the object when a failed load is caught. A class that provides a scratch object (`make_binary_scratch` in C++) is protected more fully. The record is read into the scratch first and reaches the object only once the whole file has passed every check, so state the object does not save, such as a solved structure, survives a failed load too.

## Saving

`save_binary` writes the record to a temporary file beside the target (`<target>.<random>.partial`) and renames it over the target only after the whole record is written and the file is closed. A save that fails, whether the disk fills, the object cannot be written, or the target cannot be replaced, raises `IOError`, removes the temporary file, and leaves any previous file at the target unchanged. On Windows, a target that another program holds open cannot be replaced, so the save raises until the program closes it. The saved file takes the directory's default permissions, and a target that is a symbolic link is replaced by the new file rather than written through.

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

`check_binary_file(path)` reads only the header, so it is cheap and safe on a file of unknown provenance. It raises `FileNotFoundError` if the path does not exist and `IOError` if the magic bytes are wrong, the byte order is not this machine's, or the file is shorter than the header. It does not compare the payload size with the file. That check is done by `load_binary`.

`binary_file_class(path)` also reads only the header and returns the name of the class whose record the file holds, or `None` for a file that does not start with the TidalPy magic bytes (a TOML file, say). `TidalPy.Structures.load_world(path)` and `load_system(path)` use it to read a file into a new object of the class that saved it, and `build_world` and `build_system` use it to load a binary file they are given instead of reading it as TOML.

`get_current_schema_version()` returns the version compiled into this build, which is what a file would be written with now.

## C++ API

`c_TidalPyBaseClass` (`Utilities/classes/tidalpy_base_.hpp`) holds the only `write_binary` and `read_binary` in the package, and every class takes part in the format through a few overrides:

| Member | Role |
|---|---|
| `write_binary(out)` | Writes the class's payload to memory through `p_write_payload`, then the header, with the measured payload size, and the payload. |
| `read_binary(in, force)` | Reads and checks the header (`c_read_binary_record_header` and the class id), then passes `p_read_payload` a stream that holds only the payload, and raises unless it read every byte. |
| `get_binary_class_id()` | The `BinaryClassID` of the class's records. Every concrete class returns its own. |
| `write_binary_bytes()` | The object's whole record (header and payload) as a byte string, the in-memory form of `save_binary`; a world's and a system's `copy()` and pickle go through it, read back by `c_world_from_binary_bytes` (`Structures/worlds/factory_.hpp`) or `load_binary_bytes`. |
| `load_binary_bytes(record_bytes, source, force)` | `load_binary` from a record held in memory (`load_binary` reads the file and calls it); `source` names the record in the error messages. |
| `p_write_payload(out)`, `p_read_payload(in, force)` | The class's payload. An override calls its parent's first, then writes or reads its own fields and the records of the sub-objects it owns, so a record holds its parent's payload followed by its own additions. |
| `make_binary_scratch()` | Optional: a new object of the class that `load_binary` reads a file into first (see [Integrity Checks](#integrity-checks)). |

Every concrete physics model is declared through a parameter table (`c_SpecModel`, `Utilities/classes/spec_model_.hpp`: rheology, viscosity, cooling, radiogenics, tides, luminosity, the material laws, melting laws, phases, and materials). Its payload is its model name, the number of parameters (`uint32_t`), and then each parameter by key: the key, a kind byte (`c_ParamKind`), a value count (`uint64_t`), and the values (doubles). Reading goes through the same validation as construction, and a key the record does not hold reads at its default, so adding a parameter never changes how older records are read. A model with state outside its table writes it after the parameters through `p_write_extra` and reads it back through `p_read_extra`: the isotope radiogenics model writes its label count (`uint64_t`) and then each label as a string. A composite (a phase or a material) then writes one optional record per slot. A model therefore needs only its table and its `C_CLASS_ID`, from which `c_SpecModel` provides `get_binary_class_id`. The payload of a bare `c_PhysicsBase` is its model name alone.

```
[string model_name][uint32_t num_params]
    num_params x [string key][uint8_t kind][uint64_t num_values][num_values x double]
[extra state, if any]
```

Include `binary_.hpp` in any code that reads or writes these files directly. It pulls in `logger_.hpp` so version-mismatch warnings go through the shared logger.

| Function | Description |
|---|---|
| `write_binary_header(out, class_id, payload_size)` | Write the 20-byte header. |
| `read_binary_header(in)` | Read the 20-byte header, validating the magic bytes and byte order. |
| `read_binary_header_from_file(path)` | Open a path and read its header. |
| `c_binary_file_class_id(path)` | The class id of the record a file holds, or `BinaryClassID::Unknown` (0) for a file that does not start with the magic bytes. |
| `c_peek_binary_header(in)` | Read the header at the read position and rewind to it, for a factory that picks the class to read the record. |
| `check_binary_schema_version(header, force)` | Validate the version, logging a warning on mismatch. |
| `c_read_binary_record_header(in, force)` | Read a header before loading its record: validates the magic bytes, byte order, and schema version, and refuses a payload larger than the rest of the stream. |
| `binary_bytes_remaining(in)` | The bytes left in a seekable stream after its read position. |
| `check_binary_count(in, count, element_bytes, what)` | Throw when a count read from a file needs more bytes than are left. |
| `c_host_binary_byte_order()` | The `byte_order` value this machine writes. |
| `c_binary_class_name(class_id)` | The Python class a record of that class id loads into (`"TerrestrialWorld"`, `"Maxwell"`), which the load errors name. |
| `c_binary_temporary_path(target)` | The temporary sibling path `save_binary` writes before renaming over the target. |
| `write_binary_string(out, text)` and `read_binary_string(in)` | Length-prefixed string input and output. |
| `write_optional_binary(out, pointer)` | Write an optional owned sub-object: a presence flag, then its record if present. |
| `read_optional_binary<T>(in, force, factory)` | Read an optional sub-object, rebuilding it through the given factory. |

| Constant | Value |
|---|---|
| `TIDALPY_SCHEMA_MAJOR`, `TIDALPY_SCHEMA_MINOR`, `TIDALPY_SCHEMA_PATCH` | `0`, `2`, `0` |
| `TIDALPY_BINARY_MAGIC` | `"TPYB"` |
| `TIDALPY_BINARY_HEADER_BYTES` | `20` |
| `TIDALPY_BINARY_LITTLE_ENDIAN`, `TIDALPY_BINARY_BIG_ENDIAN` | `0`, `1` |

## Variable-length Strings

Model names, layer names, and material names are written as a `uint32_t` length followed by the raw UTF-8 bytes:

```
[uint32_t length][length bytes of UTF-8 text]
```

## Nested and Recursive Serialization

Containers own sub-objects that have to round-trip with them: a layer owns its physics models, a world owns its layers, a system owns its worlds. The encoding is the same at every level.

An optional owned sub-object is written as a one-byte presence flag, zero for absent and one for present. When present, the sub-object's own complete record follows immediately, header and all. Both belong to the owning record's payload, so a file is one root record with every sub-object record nested inside it.

On read, the owning class reads the flag and, when set, calls a binary-dispatch factory. The factory peeks at the upcoming record's class id, default-constructs the matching concrete subclass, and delegates to its `read_binary`.

| Module | Dispatch factory |
|---|---|
| `Rheology` | `c_rheology_from_binary` |
| `Viscosity` | `c_viscosity_from_binary` |
| `PartialMelt` | `c_melting_curve_from_binary`, `c_melt_weakening_from_binary`, `c_bulk_modulus_mixing_from_binary`, `c_bulk_viscosity_mixing_from_binary` |
| `Cooling` | `c_cooling_from_binary` |
| `Radiogenics` | `c_radiogenics_from_binary` |
| `Material` | `c_eos_from_binary`, `c_shear_modulus_from_binary`, `c_phase_from_binary`, `c_material_from_binary` |
| `Stellar` | `c_luminosity_from_binary` |
| `Tides` | `c_tide_from_binary` |
| `Structures.layers` | `c_layer_from_binary` |

```cpp
// Writing, in c_Layer::p_write_payload.
write_optional_binary(out, this->p_shear_rheology);

// Reading, in c_Layer::p_read_payload.
this->p_shear_rheology = read_optional_binary<c_RheologyBase>(in, force, c_rheology_from_binary);
```

The sub-objects a layer and its material carry:

| Class | Recursively serialized sub-objects |
|---|---|
| `c_Layer` | material, shear rheology override, bulk rheology override, cooling, radiogenics (after its geometry, switches, and flags) |
| `c_Material` | solid phase, liquid phase, solidus, liquidus, melt weakening, bulk-modulus mixing, bulk-viscosity mixing (after its parameters) |
| `c_Phase` | equation of state, shear modulus, shear viscosity, bulk viscosity, shear rheology, bulk rheology (after its thermal parameters) |

Worlds carry their layers and world-scale models the same way:

| World | Recursively serialized sub-objects |
|---|---|
| `c_BaseWorld`, `c_TerrestrialWorld`, `c_GasGiantWorld` | tide model (after the tide configuration scalars), then every layer in index order (the spin model's moment-of-inertia factor and the pinned solver settings are scalars in the payload) |
| `c_StarWorld` | the `c_BaseWorld` sub-objects, then the luminosity model |

> [!NOTE]
> The equation-of-state profile data is never serialized, because it is derived from the layers' materials. Run `solve_eos` on a loaded world to regenerate it.

## Class Type IDs

Each concrete class needs a unique id so the dispatch factories can reconstruct the right subclass. The ranges are grouped by family, which leaves room to add models without renumbering.

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

ID zero is `Unknown` and is never written.

## Portability Notes

Every field uses a fixed-width type from `<cstdint>`, and fields are written individually rather than as a struct, so the layout does not depend on the compiler. File paths are passed as UTF-8 `std::string`. Non-ASCII paths on Windows are not guaranteed to work everywhere, so prefer ASCII paths when portability matters. `binary_.hpp` is header-only.

The only external dependency is [spdlog](https://github.com/gabime/spdlog/releases/tag/v1.15.3), reached through `logger_.hpp`, and only for the version-mismatch warnings.
