# Utilities_x: Base Class Hierarchy (`classes_x`)

`classes_x` provides the three C++ base classes that the entire TidalPy class
system is built on, plus Cython/Python wrappers exposing them to Python.

---

## Inheritance Chain

```
c_TidalPyBaseClass          (abstract; binary I/O, schema version)
    ├── c_StructureBase     (spherical geometry: radius, mass, calc_* methods)
    └── c_PhysicsBase       (physics models: model_name, layer observer ptr)
```

In later phases every layer (`c_BaseLayer`, `c_PhysicsLayer`, ...) inherits
`c_StructureBase`, and every physics model (`c_RheologyBase`, `c_CoolingBase`,
...) inherits `c_PhysicsBase`.

---

## Python API

```python
from TidalPy.Utilities_x.classes_x import StructureBase, PhysicsBase

# StructureBase — radius and mass in MKS
earth = StructureBase(radius_m=6.371e6, mass_kg=5.972e24)
print(earth.radius)   # 6371000.0
print(earth.mass)     # 5.972e+24
print(earth.get_schema_version_str())  # '0.2.0'

# Geometry helpers (explicit args, not implicit object state)
area  = earth.calc_surface_area(earth.radius)       # m^2
vol   = earth.calc_volume_sphere(earth.radius)       # m^3
g     = earth.calc_surface_gravity(earth.mass, earth.radius)   # m/s^2
rho   = earth.calc_mean_density(earth.mass, vol)     # kg/m^3
v_esc = earth.calc_escape_velocity(earth.mass, earth.radius)   # m/s

# Binary serialization
earth.save_binary("earth.tpyb")
earth2 = StructureBase(0.0, 0.0)
earth2.load_binary("earth.tpyb")
assert earth2.radius == earth.radius

# TOML config
earth.save_config("earth.toml")  # writes {"radius_m": ..., "mass_kg": ...}

# Config dict
cfg = earth.get_config_dict()  # {"radius_m": 6371000.0, "mass_kg": 5.972e24}

# PhysicsBase — model name
ph = PhysicsBase(model_name="maxwell")
print(ph.model_name)   # 'maxwell'
ph.model_name = "andrade"
ph.save_binary("maxwell.tpyb")
```

---

## `TidalPyBaseClass`

Abstract base. Not directly instantiable — only use concrete subclasses.

| Method | Returns | Description |
|--------|---------|-------------|
| `get_schema_version_str()` | `str` | Schema version, e.g. `"0.2.0"` |
| `save_binary(path)` | — | Serialize to binary file |
| `load_binary(path, force=False)` | — | Load from binary file |
| `get_config_dict()` | `dict` | Returns `{}` in base; overridden in subclasses |
| `save_config(path)` | — | Write `get_config_dict()` to TOML file |

---

## `StructureBase`

```python
StructureBase(radius_m: float, mass_kg: float)
```

| Property/Method | Returns | Description |
|----------------|---------|-------------|
| `.radius` | `float` | Stored radius [m] |
| `.mass` | `float` | Stored mass [kg] |
| `calc_surface_area(radius)` | `float` | `4π r²` [m²] |
| `calc_volume_sphere(radius)` | `float` | `(4/3)π r³` [m³] |
| `calc_volume_shell(r_outer, r_inner)` | `float` | Shell volume [m³] |
| `calc_surface_gravity(mass, radius)` | `float` | `G m / r²` [m/s²] |
| `calc_mean_density(mass, volume)` | `float` | `m / V` [kg/m³] |
| `calc_escape_velocity(mass, radius)` | `float` | `√(2Gm/r)` [m/s] |

All `calc_*` methods are const and take explicit parameters — they do not read
the object's stored `radius`/`mass` implicitly.

**Binary format** (36 bytes total):
- 20-byte header (`BinaryClassID::StructureBase = 2`, payload 16)
- `p_radius` (double, 8 bytes, host byte order)
- `p_mass` (double, 8 bytes, host byte order)

---

## `PhysicsBase`

```python
PhysicsBase(model_name: str)
```

| Property/Method | Returns | Description |
|----------------|---------|-------------|
| `.model_name` | `str` | Physics model name (read/write) |
| `get_config_dict()` | `dict` | `{"model_name": "..."}` |

The layer observer pointer (`p_layer_ptr`) is a C++ only field set by the owning
layer after construction. It is not serialized and not exposed to Python in Phase 0d.

**Binary format** (24 + n bytes total):
- 20-byte header (`BinaryClassID::PhysicsBase = 3`, payload `4 + n`)
- model name length (uint32_t, 4 bytes)
- model name bytes (UTF-8, n bytes)

---

## C++ API

### `tidalpy_base_.hpp`

```cpp
#include "tidalpy_base_.hpp"

// Create a concrete subclass and use its file methods:
c_StructureBase earth(6.371e6, 5.972e24);
earth.save_binary("earth.tpyb");

// Load back into another instance:
c_StructureBase earth2;
earth2.load_binary("earth.tpyb");
```

Key methods on `c_TidalPyBaseClass`:

| Method | Description |
|--------|-------------|
| `get_schema_version_str() const` | Returns `"0.2.0"` |
| `check_schema_compatibility(major, minor) const` | Checks version; logs warning |
| `write_binary(ostream&) const = 0` | Pure virtual; subclasses implement |
| `read_binary(istream&, force=false)` | Virtual; base reads/validates header |
| `save_binary(path) const` | Opens file, calls `write_binary` |
| `load_binary(path, force=false)` | Opens file, calls `read_binary` |

### `structure_base_.hpp`

```cpp
#include "structure_base_.hpp"

c_StructureBase body(1.0e6, 1.0e22);
double area = body.calc_surface_area(body.get_radius());
double g    = body.calc_surface_gravity(body.get_mass(), body.get_radius());
```

### `physics_base_.hpp`

```cpp
#include "physics_base_.hpp"

c_PhysicsBase model("maxwell");
model.set_layer_ptr(layer_ptr);  // called by owning layer
const std::string& name = model.get_model_name();
```

---

## Cross-DLL Logger Wiring

`classes.pyx` calls `set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())`
at module-init time so that `TIDALPY_LOG_*` calls inside `tidalpy_base_.hpp`
and `binary_.hpp` are routed to the shared spdlog logger. This follows the same
pattern as all other `_x` Cython extensions.

---

## Include Chain

```
classes.pyx
    ↓ cimport
classes.pxd
    ↓ cdef extern from
tidalpy_base_.hpp  ─→  binary_.hpp  ─→  logger_.hpp  ─→  spdlog
structure_base_.hpp ─→ tidalpy_base_.hpp
physics_base_.hpp   ─→ tidalpy_base_.hpp
```

`cython_extensions.json` includes all three directories (`classes_x`, `binary_x`,
`logging_x`) so each header can find its dependency via `#include "..."`.
