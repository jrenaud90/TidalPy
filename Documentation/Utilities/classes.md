# Base Classes (`Utilities.classes`)

_Updated: 2026-10-02_

Three C++ base classes underlie every object TidalPy builds. They give a rheology model, a cooling model, a layer, and a world the same methods for saving and restoring themselves, so a new physics model needs no serialization code of its own. A new model class has to satisfy the contract described here; the per-module "adding a new model" sections build on it.

```text
c_TidalPyBaseClass          (abstract; binary input and output, schema version)
    ├── c_StructureBase     (spherical geometry: radius, mass, calc_* helpers)
    └── c_PhysicsBase       (physics models: model name, generic parameter interface)
```

Every layer inherits `c_StructureBase` and every physics model `c_PhysicsBase`. The Python classes `TidalPyBaseClass`, `StructureBase`, and `PhysicsBase` expose the same three levels.

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

### Configuration Round Trip

A physics model's `get_config_dict()` holds the keys its family factory accepts, so `make_<family>(config["model"], config)` rebuilds the model. A world collects the configuration of each layer, and each layer that of each attached model, so a world's configuration is valid builder input and round-trips. The configuration is a human-readable view, separate from the binary format.

## Physics-Model Config Keys

Every physics model checks its config keys when it is built, by class or by `make_<family>` factory. A key the model does not read raises `ValueError` naming the model, the key, and the closest key it does read. The check is per model, so a key of another model in the same family is refused too (`isotopes` in a `fixed` radiogenics table, `fixed_dt_s` in a `fixed_q` tide table). The key `model` is always accepted, so a `get_config_dict()` result can be passed straight back. A parameter is accepted under its config key (`fixed_heat_production_w_kg`) or its argument name (`fixed_heat_production`), not both.

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

`<family>_config_keys(name)` lists the keys one model reads, `<family>_model_names()` the canonical names, and the constant `<FAMILY>_CONFIG_KEYS` every key of the family. The world builder adds the table name (for example `[layers.mantle.radiogenics]`) to the message, to locate the line in the TOML file.

## `TidalPyBaseClass`

Abstract; it provides the file and version methods.

| Method | Returns | Description |
|---|---|---|
| `get_schema_version_str()` | `str` | The schema version, for example `"0.2.0"`. |
| `save_binary(path)` | | Serialize to a binary file (`path` a `str` or `os.PathLike`). |
| `load_binary(path, force=False)` | | Load from a binary file. A file of another class is refused with both classes named ("it is a Sundberg file, not a Maxwell one"). |
| `get_config_dict()` | `dict` | Empty at this level; subclasses fill it. |
| `save_config(path)` | | Write `get_config_dict()` as TOML, under the comment header naming the TidalPy, SciPy, and CyRK versions, with LF newlines. |

A file written by a different minor version is refused, since its class layout may no longer match; `force=True` bypasses this with a warning (see [Binary Serialization](binary.md)).

## `StructureBase`

`StructureBase(radius, mass)`, with `radius` \[m\] and `mass` \[kg\] as floats.

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

The `calc_` methods take their inputs explicitly, not from the stored radius and mass.

## `PhysicsBase`

`PhysicsBase(model_name)`, with `model_name` a string.

| Property or method | Returns | Description |
|---|---|---|
| `.model_name` | `str` | The model's canonical name. Read-only: a different model is a new object. |
| `get_config_dict()` | `dict` | `{"model": ...}` plus the model's own parameters. |
| `parameters` | `dict` | Every parameter by argument name (no unit suffix); a table as a list. |
| `get_parameter(name)` | value | One parameter by argument name or config key. |
| `get_parameter_info()` | `list` | One dict per parameter: `name`, `key` (the config key), `kind`, `default`, `bounds`, and `doc`. |
| `with_parameters(**changes)` | model | A new model with some parameters changed, validated like a new one; this model is unchanged. |

Every concrete physics model is declared through a parameter table, which gives it these methods and attribute access to its parameters (`model.alpha`). A bare `PhysicsBase` has no parameters, and its copying methods raise. Models are never changed in place, so one model can be shared by several layers or worlds, and `copy.copy` returns the model itself. `repr(model)` shows the class, the model name, and the first three parameters (`Andrade('andrade', alpha=0.3, zeta=1)`, with `...` when there are more); `Phase` and `Material` list their components instead.

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
| `save_binary(path) const` and `load_binary(path, force = false)` | Delegate to the two above, with the file safeguards of [Binary Serialization](binary.md#saving). |

### `config_entry_.hpp`

A model's configuration comes from the virtual `append_config_entries(std::vector<c_ConfigEntry>&)`: the base pushes the model name and a parameter-table model (`c_SpecModel`) appends each parameter. A model with state outside its table extends the override:

```cpp
// The isotope radiogenics model writes its own entries.
void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
    c_RadiogenicsBase::append_config_entries(out);  // pushes {"model": "isotope"}
    out.push_back(c_config_doubles("heat_production_w_kg", this->p_heat_production));
    // ... the other three tables and ref_time_s ...
    if (!this->p_isotope_names.empty()) {
        out.push_back(c_config_strings("isotope_names", this->p_isotope_names));
    }
}
```

A composite (a phase or a material) adds one nested table per filled slot. `c_ConfigEntry` holds a key, a kind tag, and one payload (a double, 64-bit integer, bool, string, list of doubles or strings, nested table, or list of tables), built by `c_config_double`, `c_config_int`, `c_config_bool`, `c_config_string`, `c_config_doubles`, `c_config_strings`, `c_config_table`, and `c_config_table_list`. `c_PhysicsBase::get_config_entries()` returns the filled vector.

### `structure_base_.hpp` and `physics_base_.hpp`

`c_StructureBase` has `get_radius()`, `get_mass()`, and the `calc_` methods above. `c_PhysicsBase` has `get_model_name()` and declares the generic parameter interface that `c_SpecModel` implements from a model's table: `get_parameter_info()`, `get_parameter(name_or_key)`, `clone_physics()`, and `with_parameters(changes)`, plus `get_family_name()`, which picks the Python family class that wraps a model. Its binary payload is the model name alone; a spec model follows it with its parameters by key (see [Binary Serialization](binary.md#c-api)).

A compiled extension that logs from these headers must wire itself to the shared logger at module init; see [Logging](logging.md#sharing-the-logger-pointer-across-extensions).
