# Base Classes (`Utilities_x.classes_x`)

_Updated: 2026-09-12_

Three C++ base classes sit underneath every object TidalPy builds. They are the reason a rheology model, a cooling model, a layer, and a whole world all answer to the same four methods for saving and restoring themselves, and the reason adding a new physics model does not mean writing serialization for the fourth time.

If you plan to add a model of any kind, read this page first. The contract described here is what a new class has to satisfy, and the per-module "adding a new model" sections assume it.

## Inheritance chain

```
c_TidalPyBaseClass          (abstract; binary input and output, schema version)
    ├── c_StructureBase     (spherical geometry: radius, mass, calc_* helpers)
    └── c_PhysicsBase       (physics models: model name, layer observer pointer)
```

Every layer inherits `c_StructureBase`, and every physics model inherits `c_PhysicsBase`. The Cython wrappers `TidalPyBaseClass`, `StructureBase`, and `PhysicsBase` expose the same three levels to Python.

## Python API

```python
from TidalPy.Utilities_x.classes_x import StructureBase, PhysicsBase

body = StructureBase(radius=6.371e6, mass=5.972e24)
body.radius                       # 6371000.0
body.mass                         # 5.972e+24
body.get_schema_version_str()     # '0.2.0'

# Geometry helpers take explicit arguments rather than reading the stored state.
body.calc_surface_area(body.radius)                    # [m^2]
body.calc_volume_sphere(body.radius)                   # [m^3]
body.calc_surface_gravity(body.mass, body.radius)      # [m s-2]
body.calc_escape_velocity(body.mass, body.radius)      # [m s-1]

# Serialization, inherited by everything.
body.save_binary("body.tpyb")
body.get_config_dict()            # {'radius_m': 6371000.0, 'mass_kg': 5.972e+24}
body.save_config("body.toml")

model = PhysicsBase(model_name="maxwell")
model.model_name                  # 'maxwell', readable and writable
model.get_config_dict()           # {'model': 'maxwell'}
```

## `TidalPyBaseClass`

Abstract; instantiate a concrete subclass. It provides the file and version surface.

| Method | Returns | Description |
|---|---|---|
| `get_schema_version_str()` | `str` | The schema version, for example `"0.2.0"`. |
| `save_binary(path)` | | Serialize to a binary file. |
| `load_binary(path, force=False)` | | Load from a binary file. |
| `get_config_dict()` | `dict` | Empty at this level; subclasses fill it. |
| `save_config(path)` | | Write `get_config_dict()` as TOML. |

The version check is the part worth knowing about. A file written by a different minor version is refused, because the class layout it encodes may no longer match. `force=True` bypasses the refusal with a warning, which is occasionally what you want and is never safe to assume.

## `StructureBase`

```python
StructureBase(radius: float, mass: float)
```

| Property or method | Returns | Description |
|---|---|---|
| `.radius` | `float` | Stored radius [m]. |
| `.mass` | `float` | Stored mass [kg]. |
| `calc_surface_area(radius)` | `float` | $4 \pi r^2$ [m$^2$]. |
| `calc_volume_sphere(radius)` | `float` | $\tfrac{4}{3} \pi r^3$ [m$^3$]. |
| `calc_volume_shell(radius_outer, radius_inner)` | `float` | Shell volume [m$^3$]. |
| `calc_surface_gravity(mass, radius)` | `float` | $G m / r^2$ [m s$^{-2}$]. |
| `calc_mean_density(mass, volume)` | `float` | $m / V$ [kg m$^{-3}$]. |
| `calc_escape_velocity(mass, radius)` | `float` | $\sqrt{2 G m / r}$ [m s$^{-1}$]. |

Every `calc_` method is const and takes its inputs explicitly rather than reading the object's stored radius and mass. That is deliberate: a layer needs the volume of a shell between two radii that are not its own, and a world needs the surface area at an arbitrary radius, so binding these helpers to the object's own state would make them useless in exactly the cases they are called for.

The binary record is 36 bytes: the 20-byte header, then the radius and mass as doubles in host byte order.

## `PhysicsBase`

```python
PhysicsBase(model_name: str)
```

| Property or method | Returns | Description |
|---|---|---|
| `.model_name` | `str` | The physics model's resolved name, readable and writable. |
| `get_config_dict()` | `dict` | `{"model": ...}` plus the model's own parameters. |

Every physics model's configuration comes from one place. The C++ base declares the virtual `append_config_entries(std::vector<c_ConfigEntry>&)`, which pushes the model name; each concrete model calls its parent and then appends its own parameters using the builders in `config_entry_.hpp`. The Cython `get_config_dict` converts the entries to a dict, so the wrapper classes never override it, and a layer or world writer can read the configuration of any attached model through its raw pointer.

The payoff is that the keys are exactly what the matching factory accepts, so `make_<family>(config["model"], config)` rebuilds the model. That is what makes a world's configuration round-trip: the world asks each layer, each layer asks each attached model, and every answer is valid builder input.

The config entries are not part of the binary format; they are a separate, human-readable view. The layer observer pointer is a C++ only field that the owning layer sets after construction, and it is neither serialized nor exposed to Python.

The binary record is 24 bytes plus the model name: the 20-byte header, the name length as a `uint32_t`, then the UTF-8 name bytes.

## Checking physics-model config keys

`check_config_keys(config, accepted_keys, family)` is the guard every `make_*` factory runs before building a model. It raises `ValueError` for any key that no model in the family reads, always accepts `model` so a `get_config_dict()` result can be passed straight back, and names the closest accepted key for each rejected one. That last part matters because the most common mistake is a missing unit suffix, such as `solidus` for `solidus_k`.

```python
from TidalPy.Utilities_x.classes_x import check_config_keys

check_config_keys({"model": "henning", "solidus_k": 1500.0}, {"solidus_k", "liquidus_k"}, "partial-melt")
check_config_keys({"solidus": 1500.0}, {"solidus_k", "liquidus_k"}, "partial-melt")
# ValueError: TidalPy: unrecognized partial-melt config key(s): 'solidus' (did you mean 'solidus_k'?). ...
```

The check is per family rather than per model on purpose. The world builder merges material defaults beneath a user's table, so a table can legitimately carry a key that belongs to a different model of the same family; only a key that no model reads is an error. The world builder adds the table name to the message, for example `[layers.mantle.partial_melt]`, so the offending line can be found in the TOML file.

## C++ API

### `tidalpy_base_.hpp`

```cpp
#include "tidalpy_base_.hpp"

c_StructureBase body(6.371e6, 5.972e24);
body.save_binary("body.tpyb");

c_StructureBase restored;
restored.load_binary("body.tpyb");
```

| Method | Description |
|---|---|
| `get_schema_version_str() const` | The schema version string. |
| `check_schema_compatibility(major, minor) const` | Version check; logs a warning on mismatch. |
| `write_binary(ostream&) const` | Pure virtual; every subclass implements it. |
| `read_binary(istream&, force = false)` | Virtual; the base reads and validates the header. |
| `save_binary(path) const` and `load_binary(path, force = false)` | Open the file and delegate to the two above. |

### `config_entry_.hpp`

```cpp
// A concrete model reports its parameters by extending its parent's entries.
void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
    c_RheologyBase::append_config_entries(out);   // pushes {"model": "andrade"}
    out.push_back(c_config_double("alpha", this->p_alpha));
    out.push_back(c_config_double("zeta",  this->p_zeta));
}
```

`c_ConfigEntry` carries a key, a kind tag, and one payload: a double, a 64-bit integer, a bool, a string, a list of doubles, or a list of strings. The builders `c_config_double`, `c_config_int`, `c_config_bool`, `c_config_string`, `c_config_doubles`, and `c_config_strings` construct them. `c_PhysicsBase::get_config_entries()` returns the filled vector.

### `structure_base_.hpp` and `physics_base_.hpp`

```cpp
#include "structure_base_.hpp"
#include "physics_base_.hpp"

c_StructureBase body(1.0e6, 1.0e22);
const double gravity = body.calc_surface_gravity(body.get_mass(), body.get_radius());

c_PhysicsBase model("maxwell");
model.set_layer_ptr(layer_ptr);           // called by the owning layer
const std::string& name = model.get_model_name();
```

## Logger wiring across extensions

`classes.pyx` calls `set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())` when the module initializes, so the logging macros inside `tidalpy_base_.hpp` and `binary_.hpp` reach the shared logger. Every compiled extension in the new backend does the same thing, for the reason explained under [Logging](logging_x.md): each extension is a separate dynamic library with its own copy of the header-only state, and the pointer has to be handed across explicitly.

## Include chain

```
classes.pyx
    -> classes.pxd
        -> tidalpy_base_.hpp   -> binary_.hpp -> logger_.hpp -> spdlog
        -> structure_base_.hpp -> tidalpy_base_.hpp
        -> physics_base_.hpp   -> tidalpy_base_.hpp
```

`cython_extensions.json` lists all three directories, `classes_x`, `binary_x`, and `logging_x`, so each header resolves its dependencies.
