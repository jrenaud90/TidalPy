# World Configuration & TOML Schema (`structures_x.configs`)

The `structures_x` configuration system builds a fully wired world (the world
object, its inner-to-outer stack of layers, and each layer's attached physics
models) from a single TOML file or an equivalent Python `dict`, and writes a world
back out to TOML. It is the user-facing entry point to TidalPy's class system.

Following the `structures_x` design, **C++ never touches TOML**. Files are read and
written at the Python/Cython level with the `toml` package, converted to a `dict`,
validated against the schema, and handed to the builder, which calls the layer/world
constructors and the physics-model factories.

The current schema is version **`0.2.0`**. A config's `schema_version` is checked
against it with a graded policy: a **patch** difference (`0.0.X`) is allowed
silently, a **minor** difference (`0.X.0`) is allowed with a warning that some
functionality may break, and a **major** difference (`X.0.0`) is refused with a
`ValueError`. A missing `schema_version` is allowed with a warning. Pass
`force=True` (to `build_world` / `BaseWorld.build` / `validate_schema_version`) to
bypass these checks entirely.

---

## Quick start

```python
from TidalPy.structures_x import build_world, available_worlds

# Build one of the bundled example worlds by name.
print(available_worlds())            # ['earth_prem', 'earth_simple', 'jupiter_simple', 'sol', ...]
earth = build_world("earth_simple")  # returns the Cython world (a BaseWorld subclass)

# build_world returns the world object directly, so its methods are immediate.
earth.solve_eos()
print(earth.surface_gravity_eos, earth.planet_mass_eos)

# Build from a file path or an in-memory dict instead.
world = build_world("/path/to/my_world.toml")
world = build_world(my_config_dict)

# Save a (possibly modified) world back out; schema_version is stamped on write.
earth.save_to_toml("earth_copy.toml")
```

`build_world(source)` returns the underlying Cython world directly (a `BaseWorld`
subclass: `LayeredWorld`, `GasGiantWorld`, or `StarWorld`). It is a thin wrapper over
`BaseWorld.build(source)` (the build logic lives on the world class). `source` may be
a bundled world name, a path to a `.toml` file, or a configuration `dict`.

---

## Bundled worlds (`WorldPack_x`)

A set of example worlds (schema `0.2.0`) ships in the package directory
`TidalPy/WorldPack_x/`. The install/resolution mechanism is documented 
in [`worldpack.md`](worldpack.md); in brief, these are copied into a version-scoped,
user-editable data directory on first use:

```
<user documents>/TidalPy/<version>/Worlds_x/
```

When a world is requested by bare name, the **data-directory copy is preferred**
over the packaged copy, so editing the installed TOML (e.g.
`.../TidalPy/<version>/Worlds_x/earth_simple.toml`) changes the world a user gets
from `build_world("earth_simple")` without touching the installed package. The copy
is per-file and only happens when the data directory does not already hold a file of
that name, so user edits are never overwritten; worlds newly added to the package
appear on the next run. Pass `force=True` to `install_worldpack_x` to re-copy the
packaged versions and discard local edits.

```python
from TidalPy.structures_x import available_worlds, install_worldpack_x
install_worldpack_x()          # copy packaged worlds into the data dir (copy-if-absent)
print(available_worlds())      # data-dir worlds unioned with packaged worlds
```

---

## World-level schema

| Key | Required | Applies to | Description |
|-----|----------|------------|-------------|
| `schema_version` | optional | all | Schema marker (e.g. `"0.2.0"`). Validated when present; stamped on save. |
| `name` | **yes** | all | World name. |
| `type` | **yes** | all | `star`, `gasgiant`, `terrestrial`, or `layered`. |
| `radius_m` | **yes** | all | World radius [m]. |
| `mass_kg` | **yes** | all | World mass [kg]. |
| `albedo` | optional | all | Bond albedo. |
| `emissivity` | optional | all | Surface emissivity. |
| `obliquity_rad` | optional | all | Axial obliquity [rad]. |
| `spin_frequency_rad_s` | optional | all | Rotation rate [rad/s]. |
| `effective_temperature_k` | optional | `star` | Effective temperature [K]. |
| `luminosity_w` | optional | `star` | Luminosity [W]. |
| `[layers.<name>]` | **yes** (non-star) | layered families | One table per layer (see below). |

World `type` maps to a class as follows:

| `type` | Class |
|--------|-------|
| `terrestrial`, `layered` | `LayeredWorld` |
| `gasgiant` | `GasGiantWorld` |
| `star` | `StarWorld` (no layers) |

Optional physical keys that are omitted are **not** passed to the constructor, so
the constructor's own default applies. Defaults live in exactly one place (the C++
class or the physics-model factory) and are never duplicated in the loader.

---

## Layer-level schema

Each non-star world declares one or more `[layers.<layer_name>]` tables. The table
key is the layer's name. Layers are ordered inner-to-outer by their `layer_index`
when given, otherwise by declaration order; each layer's inner radius is assumed to be equal to 
the previous layer's outer radius (the innermost layer has a inner radius of 0).

A layer carries two distinct keys: `class` selects which Cython layer class to
build, and the optional `type` names a material whose default parameters are pulled
from the `_x` config (see "Default resolution" below).

_Most layers for rocky or icy planets and moons should use the `solidliquid` class._

| Key | Required | Layer classes | Description |
|-----|----------|---------------|-------------|
| `class` | **yes** | all | `base`, `physics`, `solidliquid`, or `gas`. Selects the layer class. |
| `type` | optional | all | Material type for default lookup: `gas`, `mantle_rock`, `ice`, `hp_ice`, or `iron`. |
| `layer_index` | optional | all | Inner-to-outer position (0 = innermost). Falls back to declaration order. |
| `radius_outer_m` | one-of | all | Outer radius [m] (absolute). |
| `radius_fraction` | one-of | all | Outer radius as a fraction of the world radius (`radius_outer_m = radius_fraction * world radius_m`). |
| `volume_fraction` | one-of | all | Layer shell volume as a fraction of the whole-world volume; the outer radius is solved from it. |
| `mass_kg` | optional | all | Layer mass [kg]. Defaults to 0.0; the EOS solve recomputes the structure. |
| `material_name` | optional | all | Free-form material label. |
| `is_tidal` | optional | all | Whether the layer participates in tides. |
| `tidal_scale` | optional | all | Tidal scaling factor. |
| `shear_modulus_static_pa` | optional | physics, solidliquid, gas | Static shear modulus [Pa]. |
| `bulk_modulus_static_pa` | optional | physics, solidliquid, gas | Static bulk modulus [Pa]. |
| `shear_viscosity_static_pas` | optional | physics, solidliquid, gas | Static shear viscosity [Pa s]. |
| `bulk_viscosity_static_pas` | optional | physics, solidliquid, gas | Static bulk viscosity [Pa s]. |
| solidliquid thermal/melt params | optional | solidliquid | See below. |
| gas params | optional | gas | See below. |

**Solid-liquid thermal and melting parameters:**
- `thermal_conductivity_ref_w_mk`
- `thermal_expansion_ref_1_k`
- `heat_capacity_ref_j_kgk`
- `activation_energy_j_mol`
- `activation_volume_m3_mol`
- `solidus_temperature_k`
- `liquidus_temperature_k`
- `melt_fraction_exponent`
- `reference_density_kg_m3`
- `reference_temperature_k`
- `melt_viscosity_reduction`.

**Gas layer parameters:**
- `mean_molecular_weight_kg_mol`
- `adiabatic_index`
- `reference_temperature_k`
- `reference_density_kg_m3`.

An unrecognized scalar key (or a model table not allowed for the layer's class) is a
validation error, which protects against typos.

#### Geometry: inner radius is derived; outer radius is one specifier

Layers are always built **inner-to-outer**, so a layer's inner radius is never
written by the user: it is the previous layer's outer radius (0 for the innermost).
Supplying `radius_inner_m`, for example, will raise an error. Each layer must specify its
outer radius with **exactly one** of `radius_outer_m`, `radius_fraction`, or `volume_fraction`
(supplying more than one, or none, is an error). For `volume_fraction`, the layer's
spherical-shell volume equals that fraction of the whole-world volume, i.e.
`r_out = (r_in^3 + volume_fraction * R_world^3)^(1/3)`.

### Attached physics models

A layer attaches a physics model through a nested table carrying a `model` key plus
that model's parameters. Every other key in the table is forwarded verbatim to the
matching `make_*` factory as its parameter dict (omitted keys keep their factory
defaults). The model tables and the layer types that may hold them are:

| Model table | Factory | Allowed layer classes |
|-------------|---------|-----------------------|
| `[layers.<name>.eos]` | `make_material_eos` | base, physics, solidliquid, gas |
| `[layers.<name>.shear_rheology]` | `make_rheology` | physics, solidliquid, gas |
| `[layers.<name>.bulk_rheology]` | `make_rheology` | physics, solidliquid, gas |
| `[layers.<name>.shear_viscosity]` | `make_viscosity` | physics, solidliquid, gas |
| `[layers.<name>.bulk_viscosity]` | `make_viscosity` | physics, solidliquid, gas |
| `[layers.<name>.partial_melt]` | `make_partial_melt` | physics, solidliquid, gas |
| `[layers.<name>.cooling]` | `make_cooling` | solidliquid only |
| `[layers.<name>.radiogenics]` | `make_radiogenics` | solidliquid only |

See each module's documentation for the available model names and parameters.

---

## Default resolution (three tiers)

Any layer parameter or physics-model table is resolved through three tiers, in order:

1. **The user world (dict / TOML).** Anything the user writes wins.
2. **The `_x` config (`TidalPy_Configs_x.toml`), keyed by material `type`.** The
   builder reads `TidalPy.config_x['layers'][<type>]` and fills in anything the user
   omitted. Material-block keys / model tables that the layer's `class` cannot hold
   are ignored (so an `ice` block applied to a `physics` layer simply drops its
   cooling/radiogenics sections). With no `type` set, this tier is skipped.
3. **The constructor / factory default.** Anything still unset falls through to the
   hardcoded C++/Cython default.

For example, an Andrade shear rheology's `zeta` for a `solidliquid` / `mantle_rock`
layer resolves as: `layers.<name>.shear_rheology.zeta` in the user world; else
`[layers.mantle_rock.shear_rheology].zeta` in `TidalPy_Configs_x.toml`; else the
Cython class' factory default.

`TidalPy_Configs_x.toml` is the main configuration file for TidalPy's new `_x` system.
It is generated from `TidalPy.defaultc_x` into the user's TidalPy `Config` directory
(next to the legacy `TidalPy_Configs.toml`) on first use and is then user-editable.
**Any new default configuration for the `_x` system belongs in `TidalPy_Configs_x.toml`
(via `defaultc_x.py`), not the legacy config.** Its `[numerical]` section also feeds
the shared C++ config singleton used by all `_x` modules (frequency / viscosity /
modulus / thickness floors).

Because the per-material defaults supply the EOS and physics models, a world can be
specified very compactly by naming only `class`, `type`, and geometry (this is how
the bundled `earth_simple` world is written).

---

## Example: two-layer terrestrial world

This world relies on per-material defaults: each layer names only its `class`,
material `type`, and geometry, and the EOS / rheology / viscosity / melt / cooling /
radiogenics come from the matching `[layers.<type>]` blocks of `TidalPy_Configs_x.toml`.

```toml
schema_version = "0.2.0"
name = "Earth-Simple"
type = "terrestrial"
radius_m = 6371000.0
mass_kg = 5.972e24
spin_frequency_rad_s = 7.292e-5

[layers.core]
class = "physics"
type = "iron"
layer_index = 0
radius_outer_m = 3480000.0   # inner radius is derived (0 for the innermost)
is_tidal = false

[layers.mantle]
class = "solidliquid"
type = "mantle_rock"
layer_index = 1
radius_fraction = 1.0        # outer radius = full world radius; inner = core's outer
is_tidal = true
```

Any default can be overridden by adding the key or sub-table. For example, to give
the mantle a specific shear viscosity and override its EOS density:

```toml
[layers.mantle.shear_viscosity]
model = "constant"
reference_viscosity = 1.0e21

[layers.mantle.eos]
model = "constant"
reference_density_kg_m3 = 4500.0
```

A star is far simpler (no layers):

```toml
schema_version = "0.2.0"
name = "Sol"
type = "star"
radius_m = 695700000.0
mass_kg = 1.988435e30
effective_temperature_k = 5772.0
```

---

## Building from a PREM-like data file (`data_file`)

Instead of writing layer tables by hand, a world can be built from a **PREM-like
radial data file** by giving a top-level `data_file` key. The layers are then
auto-detected from the data.

```toml
schema_version = "0.2.0"
name = "Earth-PREM"
type = "terrestrial"
radius_m = 6371000.0
mass_kg = 5.972e24
data_file = "PREM.csv"
```

### Data file format

A delimited table (comma, tab, or whitespace; `#` comment lines ignored) with
columns:

| Column | Quantity | Units |
|--------|----------|-------|
| 1 | radius | km |
| 2 | density | kg/m³ |
| 3 | P-wave velocity `Vp` | m/s |
| 4 | S-wave velocity `Vs` | m/s |
| 5 *(optional)* | shear viscosity | Pa s |
| 6 *(optional)* | bulk viscosity | Pa s |

The file may be ordered surface-first or center-first (it is sorted internally).
Static moduli are derived per row: shear `μ = ρ·Vs²`, bulk `K = ρ·(Vp² − 4/3·Vs²)`.

### Automatic layer detection

The profile is scanned from the center outward and split into layers by shear
modulus: `Vs = 0` (zero shear) is **liquid**, non-zero is **solid**, and every
solid↔liquid transition starts a new layer. Layers are named `layer_0`, `layer_1`,
… inner to outer. (Duplicate-radius boundary points are absorbed so no
zero-thickness layers are produced.) For PREM this yields four layers: inner core
(solid), outer core (liquid), mantle (solid), ocean (liquid).

Each detected layer gets an **interpolated EOS** carrying that layer's
radius-varying density and static shear/bulk moduli (and viscosities, if the file
has those columns). During `solve_eos` the structure ODE integrates using the
interpolated density, and the world's viscoelastic profile is taken from the
interpolated moduli/viscosities (rather than a per-layer constant).

### Refining auto-detected layers

Add `[layers.layer_N]` tables to refine the auto-detected layers (e.g. attach a
shear rheology, or override a modulus). When layer tables are provided there must be
**one per detected layer** (matched in `layer_index` order, inner to outer); a
count mismatch raises. A provided outer radius (`radius_outer_m` / `radius_fraction`)
must match the detected boundary or an error is raised. A user-provided **constant**
modulus or viscosity (e.g. `bulk_modulus_static_pa = 1.0e11`) replaces that layer's
interpolated array with the constant ("TOML overrides the data file"), while other
keys (`class`, rheology sub-tables, …) override the auto values.

```toml
[layers.layer_2]            # the solid mantle
class = "solidliquid"
layer_index = 2
bulk_modulus_static_pa = 1.0e11   # override the PREM bulk profile with a constant
[layers.layer_2.shear_rheology]
model = "maxwell"
```

The `data_file` path is resolved relative to the world TOML's directory, then the
worlds data directory, then the packaged `WorldPack_x` (see
[`worldpack.md`](worldpack.md)).

---

## System schema (`[worlds.<name>]`)

A **system** groups several worlds and the orbits that connect them into one TOML,
built with `build_system` (the system analogue of `build_world`). The file carries a
top-level `schema_version` and `name`, then one `[worlds.<key>]` table per member
world. The table key becomes that world's name within the system.

| Key | Required | Description |
|-----|----------|-------------|
| `world` | **yes** | The member world: a bundled world name, a path to a world TOML, or an inline `[worlds.<key>.world]` table. |
| `is_host` | optional | Marks the gravitational host the others orbit (at most one per system). |
| `is_star` | optional | Marks the star that lights the system (at most one; often the same body as the host). |
| `semi_major_axis_m` | optional | Orbital semi-major axis about the host [m]. |
| `eccentricity` | optional | Orbital eccentricity about the host. |
| `stellar_semi_major_axis_m` | optional | Distance from the star [m], tracked separately from the host distance so a moon can orbit a non-star host. |
| `stellar_eccentricity` | optional | Orbital eccentricity about the star. |

A system needs at least one world, at most one host, and at most one star.

```toml
schema_version = "0.2.0"
name = "Sol System"

[worlds.sun]
world = "sol"          # a bundled world name (also accepts a path or an inline table)
is_host = true
is_star = true

[worlds.earth]
world = "earth_simple"
semi_major_axis_m = 1.495978707e11
eccentricity = 0.0167
stellar_semi_major_axis_m = 1.495978707e11
stellar_eccentricity = 0.0167
```

```python
from TidalPy.structures_x.configs import build_system
system = build_system("sol_system")     # or a path / a config dict
```

The full system API (evolution, insolation, save/load) is documented in
[`../system/system.md`](../system/system.md).

---

## Python API

All entry points are re-exported from `TidalPy.structures_x` and from
`TidalPy.structures_x.configs`.

### High level

* `build_world(source, force=False) -> BaseWorld`: resolve `source` (bundled name /
  file path / dict), validate, and return the built Cython world directly.
  `force=True` bypasses the schema-version warning. Thin wrapper over
  `BaseWorld.build(source, force=False)`, which holds the build logic and returns the
  type-appropriate subclass.
* The returned world exposes its methods directly: `world.solve_eos(...)`,
  `world.solve_love_numbers(...)`, `world.get_density(r)`, etc.
* `world.save_to_toml(path, overwrite=True)`: write the retained build configuration
  (stamped with the current `schema_version`); falls back to `get_config_dict()` if
  the world was constructed directly rather than via `build_world`.
* `world.config` (alias of `world.source_config`): the normalized configuration dict
  the world was built from (`None` if constructed directly).
* `available_worlds() -> list[str]`: names of the bundled example worlds (data dir
  unioned with packaged `WorldPack_x`).
* `install_worldpack_x(force=False) -> str`: copy the packaged `WorldPack_x` worlds
  into the user data directory (copy-if-absent unless `force`); returns that
  directory.

### Low level

* `construct_world(config) -> LayeredWorld | GasGiantWorld | StarWorld`: validate a
  dict and build the underlying Cython world (and its layers).
* `construct_layer(name, layer_cfg, layer_index) -> BaseLayer`: build a single
  layer and attach its physics models.
* `save_world_to_toml(config, path, overwrite=True)`: serialize a config dict.

### Loader / validation

From `TidalPy.structures_x.configs.toml_loader`:

* `SCHEMA_VERSION`: the current schema version string (`"0.2.0"`).
* `load_toml(source)`: parse a file path or pass through a dict.
* `validate_schema_version(config, force=False)`: graded schema check (patch =
  silent, minor = warn, major = raise `ValueError`); `force=True` bypasses it.
* `validate_world_config(config)` / `validate_layer_config(name, cfg)`: structural
  validation.
* `merge_with_defaults(config)`: apply structural (non-physical) defaults.

---

## Round-trip

`build_world` retains the exact normalized configuration it built from, so a
`build -> save_to_toml -> build` cycle reproduces the same world. The saved file
always carries the current `schema_version`.

```python
world = build_world("earth_simple")
world.save_to_toml("earth_copy.toml")
reloaded = build_world("earth_copy.toml")   # identical structure
```

---
