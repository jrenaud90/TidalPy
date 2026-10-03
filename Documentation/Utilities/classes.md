# Base Classes (`Utilities.classes`)

_Updated: 2026-10-02_

Three C++ base classes underlie every object TidalPy builds. They give a rheology model, a cooling model, a layer, and a world the same methods for saving and restoring themselves, so a new physics model needs no serialization code of its own.

A new model class has to satisfy the contract described here. The per-module "adding a new model" sections build on it.

## Inheritance

```
c_TidalPyBaseClass          (abstract; binary input and output, schema version)
    ├── c_StructureBase     (spherical geometry: radius, mass, calc_* helpers)
    └── c_PhysicsBase       (physics models: model name, generic parameter interface)
```

Every layer inherits `c_StructureBase`, and every physics model inherits `c_PhysicsBase`. The Cython wrappers `TidalPyBaseClass`, `StructureBase`, and `PhysicsBase` expose the same three levels to Python.

## Python API

```python
from TidalPy.Utilities.classes import StructureBase, PhysicsBase

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
model.model_name                  # 'maxwell', read-only
model.get_config_dict()           # {'model': 'maxwell'}
```

## `TidalPyBaseClass`

This class is abstract, so instantiate a concrete subclass instead. It provides the file and version methods.

| Method | Returns | Description |
|---|---|---|
| `get_schema_version_str()` | `str` | The schema version, for example `"0.2.0"`. |
| `save_binary(path)` | | Serialize to a binary file (`path` a `str` or `os.PathLike`). |
| `load_binary(path, force=False)` | | Load from a binary file. A file of another class is refused with both classes named ("it is a Sundberg file, not a Maxwell one"). |
| `get_config_dict()` | `dict` | Empty at this level; subclasses fill it. |
| `save_config(path)` | | Write `get_config_dict()` as TOML, under the comment header naming the TidalPy, SciPy, and CyRK versions, with LF newlines. |

A file written by a different minor version is refused, because the class layout it encodes may no longer match. `force=True` bypasses the refusal with a warning.

## `StructureBase`

```python
StructureBase(radius: float, mass: float)
```

| Property or method | Returns | Description |
|---|---|---|
| `.radius` | `float` | Stored radius \[m\]. |
| `.mass` | `float` | Stored mass \[kg\]. |
| `calc_surface_area(radius)` | `float` | $4 \pi r^2$ \[m$^2$\]. |
| `calc_volume_sphere(radius)` | `float` | $\tfrac{4}{3} \pi r^3$ \[m$^3$\]. |
| `calc_volume_shell(radius_outer, radius_inner)` | `float` | Shell volume \[m$^3$\]. |
| `calc_surface_gravity(mass, radius)` | `float` | $G m / r^2$ \[m s$^{-2}$\]. |
| `calc_mean_density(mass, volume)` | `float` | $m / V$ \[kg m$^{-3}$\]. |
| `calc_escape_velocity(mass, radius)` | `float` | $\sqrt{2 G m / r}$ \[m s$^{-1}$\]. |

Every `calc_` method is const and takes its inputs explicitly rather than reading the object's stored radius and mass, because a layer needs the volume of a shell between two radii that are not its own and a world needs the surface area at an arbitrary radius.

## `PhysicsBase`

```python
PhysicsBase(model_name: str)
```

| Property or method | Returns | Description |
|---|---|---|
| `.model_name` | `str` | The physics model's canonical name. Read-only: a different model is a new object. |
| `get_config_dict()` | `dict` | `{"model": ...}` plus the model's own parameters. |
| `parameters` | `dict` | Every parameter by argument name (no unit suffix); a table as a list. |
| `get_parameter(name)` | value | One parameter by argument name or config key. |
| `get_parameter_info()` | `list` | One dict per parameter: `name`, `key` (the config key), `kind`, `default`, `bounds`, and `doc`. |
| `with_parameters(**changes)` | model | A new model with some parameters changed, validated like a new one; this model is unchanged. |

Every concrete physics model is declared through a parameter table (`c_SpecModel`, `spec_model_.hpp`), which gives it the methods above, and its parameters also read as attributes (`model.alpha`). A bare `PhysicsBase` is a name alone: it reports no parameters, and the methods that copy it raise. Models are not changed in place, so one model can be shared by several layers or worlds, and `copy.copy` returns the model itself. `repr(model)` is one line with the class, the model name, and the first three parameters in table order (`Andrade('andrade', alpha=0.3, zeta=1)`; `...` follows when there are more), the form every family shares; `Phase` and `Material` list their components instead.

Every physics model's configuration comes from one place. The C++ base declares the virtual `append_config_entries(std::vector<c_ConfigEntry>&)`, which pushes the model name, and `c_SpecModel` appends each parameter from its table. A model that holds state outside its table (the isotope labels of a radiogenics model) extends the override using the builders in `config_entry_.hpp`. The Cython `get_config_dict` converts the entries to a dict, so the wrapper classes never override it, and a layer or world writer can read the configuration of any attached model through its raw pointer.

The keys are the ones the matching factory accepts, so `make_<family>(config["model"], config)` rebuilds the model. A world's configuration therefore round-trips. The world collects the configuration of each layer, each layer collects the configuration of each attached model, and every result is valid builder input.

The config entries are not part of the binary format. They are a separate, human-readable view.

## Physics-Model Config Keys

Every physics model checks its config keys when it is built, whether through its class or its family's `make_<family>` factory. A key the model does not read raises `ValueError` naming the model, the key, and the closest key the model does read. The check is per model, so a key that another model of the same family reads is refused too: a `fixed` radiogenics table that carries `isotopes`, or a `fixed_q` tide table that carries `fixed_dt_s`, raises. The key `model` is always accepted, so a `get_config_dict()` result can be passed straight back. A parameter is accepted under either of its spellings, the config key (`fixed_heat_production_w_kg`) or the argument name (`fixed_heat_production`), but not under both at once.

```python
from TidalPy.Radiogenics import make_radiogenics, radiogenics_config_keys

print(sorted(radiogenics_config_keys("fixed")))  # ['average_half_life_s', 'fixed_heat_production_w_kg', 'ref_time_s']
try:
    make_radiogenics(
        "fixed",
        {"heat_production_w_kg": 1.0e-11})  # A key of the isotope model
except ValueError as error:
    print(error)  # ... has no parameter 'heat_production_w_kg' (did you mean 'fixed_heat_production_w_kg'?) ...
```

Each family lists its models' keys: `<family>_config_keys(name)` gives the keys one model reads and `<family>_model_names()` the canonical names, and the module constant `<FAMILY>_CONFIG_KEYS` holds every key some model of the family reads. The world builder adds the table name to the message, for example `[layers.mantle.radiogenics]`, so the offending line can be found in the TOML file.

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
| `write_binary(ostream&) const`, `read_binary(istream&, force = false)` | The record of this object. A subclass supplies `get_binary_class_id` and its payload (see [Binary Serialization](binary.md#c-api)). |
| `save_binary(path) const` and `load_binary(path, force = false)` | Delegate to the two above. A save writes a temporary file beside the target and renames it over the target, so a failed save leaves the old file intact. A load raises if bytes remain after the root record, and a load that raises an error leaves the object's saved state as it was (see [Binary Format](binary.md)). |

### `config_entry_.hpp`

```cpp
// A model with state outside its parameter table (the isotope radiogenics model) writes its own entries.
void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
    c_RadiogenicsBase::append_config_entries(out);  // pushes {"model": "isotope"}
    out.push_back(c_config_doubles("heat_production_w_kg", this->p_heat_production));
    // ... the other three tables and ref_time_s ...
    if (!this->p_isotope_names.empty()) {
        out.push_back(c_config_strings("isotope_names", this->p_isotope_names));
    }
}
```

A model declared through a parameter table gets this override from the table, and a composite (a phase or a material) adds one nested table per filled slot. `c_ConfigEntry` carries a key, a kind tag, and one payload: a double, a 64-bit integer, a bool, a string, a list of doubles, a list of strings, a nested table, or a list of tables. The builders `c_config_double`, `c_config_int`, `c_config_bool`, `c_config_string`, `c_config_doubles`, `c_config_strings`, `c_config_table`, and `c_config_table_list` construct them. `c_PhysicsBase::get_config_entries()` returns the filled vector.

### `structure_base_.hpp` and `physics_base_.hpp`

```cpp
#include "structure_base_.hpp"
#include "physics_base_.hpp"

c_StructureBase body(1.0e6, 1.0e22);
const double gravity = body.calc_surface_gravity(body.get_mass(), body.get_radius());

c_PhysicsBase model("maxwell");  // A name alone, with no parameters
const std::string& name = model.get_model_name();
```

`c_PhysicsBase` declares the generic parameter interface that `c_SpecModel` implements from a model's table: `get_parameter_info()`, `get_parameter(name_or_key)`, `clone_physics()`, and `with_parameters(changes)`, plus `get_family_name()`, which the Cython layer uses to wrap a model in its family's class. Its binary payload is the model name alone; a spec model follows it with its parameters by key (see [Binary Serialization](binary.md#c-api)).

## Logger Wiring Across Extensions

`classes.pyx` calls `set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())` when the module initializes, so the logging macros inside `tidalpy_base_.hpp` and `binary_.hpp` reach the shared logger. Every compiled extension does the same thing, for the reason explained under [Logging](logging.md): each extension is a separate dynamic library with its own copy of the header-only state, and the pointer has to be passed across explicitly.

## Include Chain

```
classes.pyx
    -> classes.pxd
        -> tidalpy_base_.hpp   -> binary_.hpp -> logger_.hpp -> spdlog
        -> structure_base_.hpp -> tidalpy_base_.hpp
        -> physics_base_.hpp   -> tidalpy_base_.hpp
```

`cython_extensions.json` lists all three directories, `classes`, `binary`, and `logging`, so each header resolves its dependencies.
