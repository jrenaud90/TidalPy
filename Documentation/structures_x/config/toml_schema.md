# World Configuration & TOML Schema (`structures_x.configs`)

_Updated: 2026-09-21_

Schema version `0.2.0`.

The `structures_x` configuration system builds a fully wired world (the world object, its inner-to-outer stack of layers, and each layer's attached physics models) from a single TOML file or an equivalent Python `dict`, and writes a world back out to TOML. It is the user-facing entry point to TidalPy's class system.

C++ never touches TOML. Files are read and written at the Python/Cython level with the `toml` package, converted to a `dict`, validated against the schema, and handed to the builder, which calls the layer and world constructors and the physics-model factories.

A config's `schema_version` is checked with a graded policy: a patch difference (`0.0.X`) is allowed silently, a minor difference (`0.X.0`) is allowed with a warning that some functionality may break, and a major difference (`X.0.0`) is refused with a `ValueError`. A missing `schema_version` is allowed with a warning. Pass `force=True` (to `build_world`, `BaseWorld.build`, or `validate_schema_version`) to bypass these checks entirely.

## Example

```python
from TidalPy.structures_x import build_world, available_worlds

# Build one of the bundled example worlds by name.
print(available_worlds())            # ['earth_prem', 'earth_simple', 'jupiter_simple', 'sol']
earth = build_world("earth_simple")  # returns the Cython world (a BaseWorld subclass)

# build_world returns the world object directly, so its methods are immediate.
earth.solve_eos()
print(earth.surface_gravity_eos, earth.planet_mass_eos)

# Build from a file path or an in-memory dict instead.
world = build_world("/path/to/my_world.toml")
world = build_world(my_config_dict)

# Save a (possibly modified) world back out; schema_version is included on write.
earth.save_to_toml("earth_copy.toml")
```

`build_world(source)` returns the underlying Cython world directly (a `BaseWorld` subclass: `LayeredWorld`, `GasGiantWorld`, or `StarWorld`). It is a thin wrapper over `BaseWorld.build(source)` (the build logic lives on the world class). `source` may be a bundled world name, a path to a `.toml` file, or a configuration `dict`.

## Bundled Worlds (`WorldPack_x`)

A set of example worlds ships in the package directory `TidalPy/WorldPack_x/`. The install/resolution mechanism is documented in [`worldpack.md`](worldpack.md); in brief, these are copied into a version-scoped, user-editable data directory on first use:

```
<user documents>/TidalPy/<TidalPy version>/Worlds_x/
```

When a world is requested by bare name, the data-directory copy is preferred over the packaged copy, so editing the installed TOML (e.g. `.../TidalPy/<TidalPy version>/Worlds_x/earth_simple.toml`) changes the world a user gets from `build_world("earth_simple")` without touching the installed package. The copy is per-file and only happens when the data directory does not already hold a file of that name, so user edits are never overwritten; worlds newly added to the package appear on the next run. Pass `force=True` to `install_worldpack_x` to re-copy the packaged versions and discard local edits.

```python
from TidalPy.structures_x import available_worlds, install_worldpack_x
install_worldpack_x()          # copy packaged worlds into the data dir (copy-if-absent)
print(available_worlds())      # data-dir worlds unioned with packaged worlds
```

## World-Level Schema

| Key | Required | Applies to | Description |
|-----|----------|------------|-------------|
| `schema_version` | optional | all | Schema marker (e.g. `"0.2.0"`). Validated when present; included on save. |
| `name` | **yes** | all | World name. |
| `type` | **yes** | all | `star`, `gasgiant`, `terrestrial`, or `layered`. |
| `radius_m` | **yes** | all | World radius \[m\]. |
| `mass_kg` | **yes** | all | World mass \[kg\]. |
| `albedo` | optional | all | Bond albedo. |
| `emissivity` | optional | all | Surface emissivity. |
| `obliquity_rad` | optional | all | Axial obliquity \[rad\]. |
| `spin_frequency_rad_s` | optional | all | Rotation rate \[rad/s\]. |
| `effective_temperature_k` | optional | `star` | Effective temperature \[K\]. |
| `luminosity_w` | optional | `star` | Luminosity \[W\]. |
| `[luminosity]` | optional | `star` | The star's mass-to-luminosity model: `model` (`fixed`, `mass_to_luminosity`, or `power_law`) plus that model's parameters, as `stellar_x.make_luminosity` takes them. Attaching it does not change the stored `luminosity_w`. |
| `moment_of_inertia_factor` | optional | layered families | $C/(MR^2)$ of the world's spin model, within $(0, 2/3]$ (0.4, a uniform sphere, when left out). It gives the moment of inertia until the EOS is solved. |
| `[layers.<name>]` | **yes** (non-star) | layered families | One table per layer (see below). |
| `[tides]` | optional | all | Tidal dissipation settings (see below). Omitted entirely, the world still gets a dissipation model from the `_x` config defaults. |

World `type` maps to a class as follows:

| `type` | Class |
|--------|-------|
| `terrestrial`, `layered` | `LayeredWorld` |
| `gasgiant` | `GasGiantWorld` |
| `star` | `StarWorld` (no layers) |

An omitted optional key resolves through the same three tiers the layers use: the world's own table, then the `[worlds]` block of `TidalPy_Configs_x.toml` (with a `[worlds.star]` sub-table for the two star-only keys), then the C++ class default. Whatever neither tier supplies is not passed to the constructor at all, so the class default applies. Defaults live in either the C++ class or the physics-model factory, and are never duplicated in the loader.

## Layer-Level Schema

Each non-star world declares one or more `[layers.<layer_name>]` tables. The table key is the layer's name. Layers are ordered inner-to-outer by their `layer_index` when given, otherwise by declaration order; each layer's inner radius is assumed to be equal to the previous layer's outer radius (the innermost layer has a inner radius of 0).

A layer carries two distinct keys: `class` selects which Cython layer class to build, and the optional `type` names a material whose default parameters are pulled from the `_x` config (see "Default resolution" below).

_Most layers for rocky or icy planets and moons should use the `solidliquid` class._

| Key | Required | Layer classes | Description |
|-----|----------|---------------|-------------|
| `class` | **yes** | all | `base`, `physics`, `solidliquid`, or `gas`. Selects the layer class. |
| `type` | optional | all | Material type for default lookup: `gas`, `mantle_rock`, `ice`, `hp_ice`, `iron`, `default` (the block a layer without a type takes, a copy of `mantle_rock`), or `none` (no material defaults at all; `get_config_dict` writes this because it lists every model explicitly). |
| `layer_index` | optional | all | Inner-to-outer position (0 = innermost). Falls back to declaration order. |
| `radius_outer_m` | one-of | all | Outer radius \[m\] (absolute). |
| `radius_fraction` | one-of | all | Outer radius as a fraction of the world radius (`radius_outer_m = radius_fraction * world radius_m`). |
| `volume_fraction` | one-of | all | Layer shell volume as a fraction of the whole-world volume; the outer radius is solved from it. |
| `mass_kg` | optional | all | Layer mass \[kg\]. Defaults to 0.0; every successful EOS solve overwrites it with the solved layer mass. |
| `material_name` | optional | all | Free-form material label. |
| `is_tidal` | optional | all | Whether the layer participates in tides. |
| `is_volume_fixed` | optional | all | `false` lets the layer grow or shrink to hold its mass while the EOS solve redistributes the interior; the layers above it move with it. Default `true`. |
| `tidal_scale` | optional | all | Tidal scaling factor, used for homogeneous tidal solvers. |
| `is_solid` | optional | physics, solidliquid, gas | `false` makes the layer a liquid in the radial Love-number solve. Default `true` (`false` for `gas`). |
| `is_static` | optional | physics, solidliquid, gas | Static approximation (no inertia) in the radial solve. Default `true`, so a liquid layer is a static liquid unless this is `false`. |
| `is_incompressible` | optional | physics, solidliquid, gas | Incompressible approximation in the radial solve. Default `false`. |
| `temperature_k` | optional | physics, solidliquid, gas | Layer temperature \[K\] at which the material's viscosity and melt models are evaluated. Default `0.0`, the cold rigid limit of the viscosity laws. |
| `use_thermal_eos` | optional | physics, solidliquid, gas | Let the density law of the layer's material see the temperature, so its density and bulk modulus depend on it (set `thermal_expansion_1_k` in the `material` table). Default `false`. |
| `use_heating` | optional | physics, solidliquid, gas | Let the world's heat sources act inside the layer during a thermal EOS solve: its `radiogenics` model then heats it, as a specific rate times the local density. Default `false`. |
| gas params | optional | gas | See below. |

The static moduli, the shear law, the thermal constants, and the viscosity and melting parameters are not layer keys. They are
properties of the material, so they live in the layer's `material` table (below), and the layer and the solve read
the same numbers. Setting one of them on the layer is a validation error whose message says where it moved.

**Gas layer parameters:**
- `mean_molecular_weight_kg_mol`
- `adiabatic_index`
- `reference_temperature_k`
- `reference_density_kg_m3`.

An unrecognized scalar key (or a model table not allowed for the layer's class) is a validation error, which protects against typos.

### Geometry

Layers are always built inner-to-outer, so a layer's inner radius is never written by the user: it is the previous layer's outer radius (0 for the innermost). Supplying `radius_inner_m`, for example, raises an error. Each layer must specify its outer radius with exactly one of `radius_outer_m`, `radius_fraction`, or `volume_fraction` (supplying more than one, or none, is an error). For `volume_fraction`, the layer's spherical-shell volume equals that fraction of the whole-world volume, so with $f_V$ the volume fraction and $R$ the world radius, $r_\mathrm{out} = \left(r_\mathrm{in}^3 + f_V R^3\right)^{1/3}$.

### Physical Values

After the keys are checked, so are their values (`validate_physical_values`, run by `validate_world_config` and therefore by every build). A configuration is refused with a `ValueError` that names the world or layer, the key, the value found, and the range allowed when:

* the world's `radius_m` or `mass_kg` is not a positive finite number, `albedo` lies outside $[0, 1]$, `emissivity` outside $(0, 1]$, or the obliquity, spin rate, stellar temperature, or luminosity is not finite (the last two may be zero, which asks for the value to be derived);
* a `radius_fraction` or `volume_fraction` lies outside $(0, 1]$, or a `radius_outer_m` is not positive;
* a layer ends at or below the top of the layer under it, or above the world's radius;
* the outermost layer stops short of the world's radius, since the layers have to fill the world (to a relative $10^{-6}$, which absorbs the roundoff of stacked fractions);
* a layer's `mass_kg`, `tidal_scale`, or `temperature_k` is negative or not finite;
* a `layer_index` is not a non-negative integer, or two layers resolve to the same index (a layer with none takes its position in the file).

### Attached Physics Models

A layer attaches a physics model through a nested table carrying a `model` key plus that model's parameters. Every other key in the table is forwarded verbatim to the matching `make_*` factory as its parameter dict (omitted keys keep their factory defaults, and a key that no model in that family reads raises `ValueError` naming the table). The model tables and the layer types that may hold them are:

| Model table | Factory | Allowed layer classes |
|-------------|---------|-----------------------|
| `[layers.<name>.material]` | `make_material_eos` | base, physics, solidliquid, gas |
| `[layers.<name>.shear_rheology]` | `make_rheology` | physics, solidliquid, gas |
| `[layers.<name>.bulk_rheology]` | `make_rheology` | physics, solidliquid, gas |
| `[layers.<name>.cooling]` | `make_cooling` | solidliquid only |
| `[layers.<name>.radiogenics]` | `make_radiogenics` | solidliquid only |

See each module's documentation for the available model names and parameters.

### The Material Table

`[layers.<name>.material]` is the layer's EOS model, and the EOS model is the layer's material:

- The density law: `model` (`"constant"`, `"bm"`, `"vinet"`, `"interpolate"`) and its parameters (`reference_density_kg_m3`, `reference_bulk_modulus_pa`, `thermal_expansion_1_k`, ...);
- The static constants `shear_modulus_static_pa`, `bulk_modulus_static_pa`, `shear_viscosity_static_pas`, and `bulk_viscosity_static_pas` (a viscosity left out is unset);
- The thermal constants `thermal_conductivity_w_mk` (default `4.0`), `heat_capacity_j_kgk` (default `1200.0`), and `thermal_expansion_1_k` (default `0.0`). There is one expansivity: it sets the adiabat and convection of a cooling layer, and the density law uses the same number, but only on a layer that sets `use_thermal_eos`;
- The static shear law $\mu = \mu_0 + \mu'_P P + \mu'_T (T - T_\mathrm{ref})$ through `shear_modulus_pressure_derivative`, `shear_modulus_temperature_derivative_pa_k` \[Pa K$^{-1}$\], and `shear_modulus_reference_temperature_k` \[K\] (defaults `0.0`, `0.0`, `300.0`);
- Three optional nested model tables, each with its own `model` key: `[layers.<name>.material.shear_viscosity]` and `[layers.<name>.material.bulk_viscosity]` (built by `make_viscosity`) and `[layers.<name>.material.partial_melt]` (built by `make_partial_melt`).

The whole table is handed to `make_material_eos`; see [Material EOS Models](../../material_x/material_eos.md). The rheology tables stay on the layer, because the rheology is the one thing that needs a frequency.

Unlike the other model tables, `material` may leave out `model` when the layer's material `type` supplies one, so a fitted number can be overridden without restating the rest:

```toml
[layers.mantle]
class = "solidliquid"
type = "mantle_rock"
radius_fraction = 1.0

[layers.mantle.material]
shear_modulus_static_pa = 4.17e10     # everything else comes from the mantle_rock defaults

[layers.mantle.material.shear_viscosity]
reference_viscosity_pas = 3.0e21      # one key of a nested default table
```

## Tidal Dissipation (`[tides]`)

An optional world-level `[tides]` table sets how the world dissipates tidal energy. It applies to every world family, stars and gas giants included, and every key is optional: what the table omits falls back to the `[tides]` block of `TidalPy_Configs_x.toml`, and what that omits falls back to a built-in default. A world with no `[tides]` table is still given a dissipation model.

| Key | Applies to | Description |
|-----|------------|-------------|
| `global_tidal_model` | all | Dissipation model: `rheology`, `cpl` (`fixed_q`), `ctl` (`fixed_dt`), or `ctl_q` (`fixed_dt_q`). Defaults per world family from `[tides.default_model]` in the `_x` config: `rheology` for terrestrial and layered worlds, `fixed_dt` for gas giants, `fixed_q` for stars. |
| `fixed_k` | all | Per-degree static potential Love numbers $k_l$, a list indexed from $l = 2$ (nine slots, $l = 2 \ldots 10$). Read only by the analytic models. A list shorter than nine zero-fills the remaining degrees, and a zero $k_l$ is no dissipation at that degree, so the list has to reach `max_degree_l`. |
| `fixed_q` | all | Per-degree tidal quality factors $Q_l$, same indexing. Read by `cpl` and `ctl_q`. |
| `fixed_dt` | all | Per-degree tidal time lags $\Delta t_l$ \[s\], same indexing. Read by `ctl` and `ctl_q`. |
| `min_degree_l` | all | Lowest harmonic degree in the mode sum. Default `2`. |
| `max_degree_l` | all | Highest harmonic degree in the mode sum. Default `2`. |
| `eccentricity_trunc_lvl` | all | Eccentricity truncation order $e^n$. Tabulated at 1, 2, 3, 4, 5, 10, 15, and 20; default `3`. An untabulated level is promoted to the next tabulated one with a once-per-session warning, so accuracy never drops silently. `eccentricity_truncation` is accepted as an alias. |
| `obliquity_trunc_lvl` | all | Obliquity truncation order $I^n$. Tabulated at 0, 1, 2, and 10; default `"off"`. `"off"` means 0 (no obliquity terms), 1 and 2 keep every term through $I^1$ and $I^2$, and `"gen"` or `"general"` means 10 (the exact, untruncated form). Untabulated integers are promoted like the eccentricity levels. `obliquity_truncation` is accepted as an alias. |
| `tidal_timescale_width_decades` | layered families | Width \[decades\] of the log-Gaussian bell a layer's `tidal_timescale` scale method uses. The bell peaks where the layer's Maxwell time equals the forcing period. Default `1.0`. |
| `love_method` | layered families | How the Love numbers are obtained: `radial_solver` (aliases `shooting`, `rs`; the default), `propagation_matrix` (`prop_matrix`, `pm`, `prop`), `homogeneous` (`homogen`), `cpl`, `ctl`, or `laterally_inhomogeneous` (`3d`, `lat_inhom`, reserved for the 3D solver). The three homogeneous methods use the analytic homogeneous-sphere formulas instead of a radial solve, so they have no depth-resolved solution and the 3D stress, strain, and heating path raises `RuntimeError` while one of them is configured. |
| `love_fixed_q` | layered families | Scalar $Q$ the `cpl` Love method applies to the static Love numbers. Unset by default, in which case the tide model's own `fixed_q` is used. |
| `love_fixed_dt` | layered families | Scalar time lag \[s\] the `ctl` Love method applies. Unset by default, falling back to the tide model's `fixed_dt`. |

`global_tidal_model` selects the model that turns Love numbers into dissipation; `love_method` selects how the Love numbers themselves are computed. They are independent: a world can solve its Love numbers with the radial solver and still collapse them with an analytic tide model.

```toml
[tides]
global_tidal_model = "rheology"      # dissipation from the layers' complex moduli
min_degree_l = 2
max_degree_l = 3                     # include the degree-3 tide
eccentricity_trunc_lvl = 5           # e^5; 3 is the default
obliquity_trunc_lvl = "gen"          # exact obliquity terms
love_method = "radial_solver"
```

A star or gas giant has no interior to solve, so it uses an analytic model and its per-degree parameters:

```toml
[tides]
global_tidal_model = "fixed_q"
fixed_k = [0.3, 0.15, 0.1]           # l = 2, 3, 4; higher degrees would be zero-filled
fixed_q = [1.0e5, 1.0e5, 1.0e5]
```

`world.get_tide_config()` returns the degree and truncation settings under these same key names, and `get_config_dict()` puts them in a `[tides]` table together with the tide model's own parameters (the model's name is emitted as `global_tidal_model`), so a world's tidal configuration survives a save and rebuild. Note that the round trip writes the resolved integer for `obliquity_trunc_lvl`, so a world written with `"off"` reads back as `0`.

## Default Configuration Resolution

Any layer parameter or physics-model table is resolved through three tiers, in order:

1. The user world (dict or TOML): anything the user writes is used over anything else.
2. The TidalPy `_x` config (`TidalPy_Configs_x.toml`), keyed by material `type`: the builder reads `TidalPy.config_x['layers'][<type>]` and fills in anything the user omitted. Material-block keys and model tables that the layer's `class` cannot hold are ignored (so an `ice` block applied to a `physics` layer drops its cooling and radiogenics sections). With no `type` set, this tier is skipped.
3. The constructor or factory default: anything still unset falls through to the C++ or Cython default.

For example, an Andrade shear rheology's `zeta` for a `solidliquid` / `mantle_rock` layer resolves as: `layers.<name>.shear_rheology.zeta` in the user world; else `[layers.mantle_rock.shear_rheology].zeta` in `TidalPy_Configs_x.toml`; else the Cython class' factory default.

World-level properties resolve the same way through the `[worlds]` block instead of `[layers.<type>]`: a world's `albedo` comes from its own table, else `[worlds].albedo`, else the class default. A `[worlds.<type>]` sub-table specializes the block for one world type, and is how the star-only `effective_temperature_k` and `luminosity_w` are kept off every other world.

`TidalPy_Configs_x.toml` is the main configuration file for TidalPy's new `_x` system. It is generated from `TidalPy.defaultc_x` into the user's TidalPy `Config` directory (next to the legacy `TidalPy_Configs.toml`) on first use and is then user-editable. Any new default configuration for the `_x` system belongs in `TidalPy_Configs_x.toml` (through `defaultc_x.py`), not the legacy config. Its `[numerical]` section also feeds the shared C++ config singleton used by all `_x` modules (frequency / viscosity / modulus / thickness floors, plus `numerical_floor`, the magnitude a guarded denominator is raised to, and `layer_continuity_rtol`, how closely a layer's inner radius must match the previous layer's outer radius; see [Constants](../../utilities_x/constants.md)).

Because the per-material defaults supply the EOS and physics models, a world can be specified very compactly by naming only `class`, `type`, and geometry (this is how the bundled `earth_simple` world is written).

## Example: Two-Layer Terrestrial World

This world relies on per-material defaults: each layer names only its `class`, material `type`, and geometry, and the EOS / rheology / viscosity / melt / cooling / radiogenics come from the matching `[layers.<type>]` blocks of `TidalPy_Configs_x.toml`.

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

Any default can be overridden by adding the key or sub-table. For example, to give the mantle a specific shear viscosity and override its EOS density:

```toml
[layers.mantle.material]
model = "constant"
reference_density_kg_m3 = 4500.0

[layers.mantle.material.shear_viscosity]
model = "constant"
reference_viscosity_pas = 1.0e21
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

## Building a World From a Radial Profile

Instead of writing layer tables by hand, a world can describe its interior with a PREM-like radial profile. This is useful if comparing to published data or if you use a more sophisticated EOS solver and want to feed those results into TidalPy for, _e.g._, Love number calculations. These files are defined via a top-level `data_file` key naming a delimited file.

```toml
schema_version = "0.2.0"
name = "Earth-PREM"
type = "terrestrial"
radius_m = 6371000.0
mass_kg = 5.972e24
data_file = "PREM.csv"  # Can be a path to a file too.
```

From Python the same world is built by handing `build_world` the arrays themselves under a `data` key, which is the `build_world` equivalent of a data file. A world gives its profile one way or the other, never both.

```python
world = build_world({
    "schema_version": "0.2.0", "name": "My-Earth", "type": "terrestrial",
    "radius_m": 6371000.0, "mass_kg": 5.972e24,
    "data": {"radius_km": radius, "density": density, "vp": vp, "vs": vs},
})
```

### Profile Format

A delimited table (comma, semicolon, tab, or whitespace; `#` comment lines ignored), or a mapping of arrays with the same names. A profile must give a radius and a density, and must give either both seismic velocities or both static moduli:

| Quantity | Accepted names | Units | Required |
|----------|----------------|-------|----------|
| radius | `radius`, `r`, `rad`, `radii` | m or km | yes, or a depth |
| depth | `depth`, `z` | m or km | converted with the world's `radius_m` |
| density | `density`, `rho`, `dens` | kg/m³ | yes |
| P-wave velocity | `vp`, `v_p`, `p_velocity`, … | m/s or km/s | yes, unless moduli are given |
| S-wave velocity | `vs`, `v_s`, `shear_velocity`, … | m/s or km/s | yes, unless moduli are given |
| shear modulus | `shear_modulus`, `mu`, `rigidity` | Pa | instead of the velocities |
| bulk modulus | `bulk_modulus`, `k`, `incompressibility` | Pa | instead of the velocities |
| shear viscosity | `shear_viscosity`, `eta`, `viscosity` | Pa s | no |
| bulk viscosity | `bulk_viscosity`, `eta_bulk`, `zeta` | Pa s | no |

Names are matched ignoring case and punctuation, so `Vp`, `V_P` and `vp` are one name, and a name may state its unit: `radius_km`, `Vp [km/s]`, `rho_kg_m3`. A unit this reader does not convert is taken to be MKS already. Columns are found by name, so their **order does not matter**; the names come from a header row or from the last `#` comment line before the data (which must name every column, so prose about the data is not mistaken for a header). A file with no header at all is read positionally as radius, density, `Vp`, `Vs`, shear viscosity, bulk viscosity.

A radius or depth with no stated unit is read as kilometers below 100 km and as meters above it, ranges that cannot overlap for a real body. The file may be ordered surface-first or center-first (it is sorted internally). Where the velocities are given, the static moduli are derived per row: shear $\mu = \rho V_s^2$, bulk $K = \rho \left(V_p^2 - \tfrac{4}{3} V_s^2\right)$.

### Layer Detection

The profile is scanned from the center outward and split into layers by shear modulus: `Vs = 0` (zero shear) is liquid, non-zero is solid, and every solid-liquid transition starts a new layer. Layers are named `layer_0`, `layer_1`, and so on, inner to outer. Duplicate-radius boundary points are absorbed so no zero-thickness layers are produced, and a duplicated boundary radius keeps the lower layer's row first whichever way the file is ordered. A liquid layer gets `is_solid = false` and `is_static = true`, so the radial solver treats it as a static liquid; a layer table can override either flag. For the bundled `PREM.csv`, which replaces PREM's 3 km ocean with the upper crust, this yields three layers: inner core (solid), outer core (liquid), mantle plus crust (solid).

Each layer's slice of the profile becomes its material: an interpolated EOS carrying that layer's density, static shear and bulk moduli, and viscosities when the profile gave any. That EOS owns those arrays, and is the only place a radial grid persists. A layer built this way takes **no defaults from a material type**, so a profile that names no viscosity produces an elastic layer: no viscosity model, no partial-melt model, and no rheology it did not ask for.

### Refining Detected Layers

A profile fixes how many layers the world has, where their boundaries are, and what each is made of. It does not hold complex moduli, so a rheology is still named in a `[layers.<name>]` table, as are a cooling model and radiogenics.

Such a table says which detected layer it refines with `layer_index` (or by being named `layer_N`), and **only the layers being refined need a table**. The table lends the layer its name. A provided outer radius (`radius_outer_m` or `radius_fraction`) must match the detected boundary, two tables may not claim the same layer, and an index outside the detected range raises. A constant modulus or viscosity (e.g. `bulk_modulus_static_pa = 1.0e11`) replaces that layer's array with the constant ("TOML overrides the data file"); other keys (`class`, `type`, the model sub-tables, …) override the detected values. Naming a `type` brings that material block's defaults back.

```toml
[layers.mantle]             # refines the outermost detected layer; the other two need no table
layer_index = 2
[layers.mantle.shear_rheology]
model = "maxwell"
[layers.mantle.material]
shear_viscosity_static_pas = 1.0e21   # the profile named no viscosity, so give one here
```

The `data_file` path is resolved relative to the world TOML's directory, then the worlds data directory, then the packaged `WorldPack_x` (see [`worldpack.md`](worldpack.md)).

## System Schema

A system groups several worlds and the orbits that connect them into one TOML, built with `build_system` (the system analogue of `build_world`). The file carries a top-level `schema_version` and `name`, then one `[worlds.<key>]` table per member world. The table key becomes that world's name within the system.

| Key | Required | Description |
|-----|----------|-------------|
| `world` | **yes** | The member world: a bundled world name, a path to a world TOML, or an inline `[worlds.<key>.world]` table. |
| `tidal_host` | optional | The table key of the world that raises this world's tides. It must be another world of the system, and may be declared later in the file. Left out, the world is not tidally forced. Two worlds may name each other; they then share one orbit, which only one of them needs to state. |
| `is_star` | optional | Marks the star that provides insolation to the system (at most one; does not have to be anyone's tidal host). |
| `semi_major_axis_m` | optional | Orbital semi-major axis about the tidal host [m]. Requires `tidal_host`. |
| `eccentricity` | optional | Orbital eccentricity about the tidal host. Requires `tidal_host`. |
| `stellar_semi_major_axis_m` | optional | Distance from the star [m], tracked separately from the host distance so a moon can orbit a non-star host. |
| `stellar_eccentricity` | optional | Orbital eccentricity about the star. |

A system needs at least one world and at most one star. There is no system-wide host: each tidally forced world names its own.

```toml
schema_version = "0.2.0"
name = "Sol System"

[worlds.sun]
world = "sol"          # a bundled world name (also accepts a path or an inline table)
is_star = true

[worlds.earth]
world = "earth_simple"
tidal_host = "sun"     # the world that raises this one's tides
semi_major_axis_m = 1.495978707e11
eccentricity = 0.0167
stellar_semi_major_axis_m = 1.495978707e11
stellar_eccentricity = 0.0167
```

```python
from TidalPy.structures_x.configs import build_system
system = build_system("sol_system")     # or a path / a config dict
```

The full system API (evolution, insolation, save/load) is documented in [`../system/system.md`](../system/system.md).

## Python API

All entry points are re-exported from `TidalPy.structures_x` and from `TidalPy.structures_x.configs`.

### High Level

* `build_world(source, force=False) -> BaseWorld`: resolve `source` (bundled name / file path / dict), validate, and return the built Cython world directly. `force=True` bypasses the schema-version warning. Thin wrapper over `BaseWorld.build(source, force=False)`, which holds the build logic and returns the type-appropriate subclass.
* `load_radial_data(source, surface_radius=None) -> dict`: read a radial profile (a data-file path or a mapping of arrays) into MKS arrays ascending in radius, the same reader `data_file` and `data` worlds use. `detect_layer_boundaries(radius, shear_modulus)` returns the `(start, end, is_solid)` runs it splits into.
* The returned world exposes its methods directly: `world.solve_eos(...)`, `world.solve_love_numbers(...)`, `world.get_density(r)`, etc.
* `world.save_to_toml(path, overwrite=True)`: write the retained build configuration (stamped with the current `schema_version`); falls back to `get_config_dict()` if the world was constructed directly rather than via `build_world`. The fallback is validated against this schema first, so it writes a buildable file or raises `ValueError`.
* `world.get_config_dict()`: the live world as a builder-valid table (`type`, name-keyed `layers` with `class` and attached-model sub-tables, `tides`, `schema_version`).
* `build_world_from_dict(config, force=False) -> BaseWorld`, `build_layer_from_dict(config) -> BaseLayer`, and `build_system_from_dict(config, force=False) -> System`: rebuild an object from the dictionary its `get_config_dict()` returns (see [Round Trip](#round-trip)). Each takes only a `dict`, leaves it unmodified, and raises `TypeError` for anything else.
* `world.config` (alias of `world.source_config`): the normalized configuration dict the world was built from (`None` if constructed directly).
* `available_worlds() -> list[str]`: names of the bundled example worlds (data dir unioned with packaged `WorldPack_x`). Bundled system files share that directory and are listed by `available_systems()` instead; `build_world` on a system config raises a `ValueError` naming `build_system`, and the reverse holds too.
* `install_worldpack_x(force=False) -> str`: copy the packaged `WorldPack_x` worlds into the user data directory (copy-if-absent unless `force`); returns that directory.

### Low Level

* `construct_world(config) -> LayeredWorld | GasGiantWorld | StarWorld`: validate a dict and build the underlying Cython world (and its layers).
* `construct_layer(name, layer_cfg, layer_index, radius_inner, radius_outer) -> BaseLayer`: build a single layer and attach its physics models.
* `save_world_to_toml(config, path, overwrite=True)`: serialize a config dict.

### Loader / Validation

From `TidalPy.structures_x.configs.toml_loader`:

* `SCHEMA_VERSION`: the current schema version string (`"0.2.0"`).
* `load_toml(source)`: parse a file path or pass through a dict.
* `validate_schema_version(config, force=False)`: graded schema check (patch = silent, minor = warn, major = raise `ValueError`); `force=True` bypasses it.
* `validate_world_config(config)` / `validate_layer_config(name, cfg)`: structural validation. `validate_world_config` ends with `validate_physical_values(config)`, the value checks of [Physical Values](#physical-values).
* `merge_with_defaults(config)`: apply structural (non-physical) defaults.

## Round Trip

`build_world` retains the exact normalized configuration it built from, so a `build -> save_to_toml -> build` cycle reproduces the same world. The saved file always carries the current `schema_version`.

```python
world = build_world("earth_simple")
world.save_to_toml("earth_copy.toml")
reloaded = build_world("earth_copy.toml")   # identical structure
```

The retained configuration is the world as its file described it. `get_config_dict()` is the world as it stands now, with every change made since the build, and `build_*_from_dict` rebuilds the same class with the same parameters and attached models from it. This holds for a layer, a world, and a system:

```python
from TidalPy.structures_x import build_world_from_dict, build_layer_from_dict, build_system_from_dict

world.set_spin_frequency(2.0e-5)                 # A change made after the build
config = world.get_config_dict()                 # Builder-valid: type, layers by name, tides, schema_version
twin = build_world_from_dict(config)             # Same class, same parameters, same models
assert twin.get_config_dict() == config

mantle = world.mantle                            # The layer named "mantle"
mantle_twin = build_layer_from_dict(mantle.get_config_dict())    # A standalone layer, owned by no world

system_twin = build_system_from_dict(system.get_config_dict())   # Worlds, tidal hosts, star, and orbits
```

The live dict carries every scalar and every attached model explicitly, so no material defaults are needed and each layer is written with `type = "none"`, which keeps a rebuild from adding the `[layers.default]` models a typeless layer would otherwise take; the dict is a frozen snapshot of the object as configured. A standalone layer's dict adds the keys a world would supply from the layer's place in its `layers` table (`name` and `radius_inner_m`) and, for a physics layer, its six Love number components (`love_number_k_re` through `love_number_l_im`, since TOML has no complex type); a world drops those keys when it nests the layer. A system's dict inlines each member world's own live dict under `world`. Solved state (the EOS, Love numbers, tides) is not configuration, so run the solves again on a rebuilt object.
