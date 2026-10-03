# World Configuration & TOML Schema (`Structures.configs`)

_Updated: 2026-10-02_

Schema version `0.2.0`.

The `Structures` configuration system builds a world (the world object, its inner-to-outer stack of layers, and each layer's material and attached physics models) from a TOML file or an equivalent Python `dict`, and writes a world back out to TOML. It is the user-facing entry point to TidalPy's class system.

C++ never touches TOML. The `toml` package reads and writes files at the Python/Cython level. The resulting `dict` is validated against the schema and passed to the builder, which calls the layer and world constructors and the physics-model factories.

A config's `schema_version` is checked with a graded policy:

- Patch difference (`0.0.X`): allowed silently.
- Minor difference (`0.X.0`): allowed with a warning that some functionality may break.
- Major difference (`X.0.0`): refused with a `ValueError`.
- Missing `schema_version` in a file: allowed with a warning. A dict built in Python without one targets the current schema and builds silently.

Pass `force=True` (to `build_world`, `BaseWorld.build`, or `validate_schema_version`) to bypass these checks entirely.

## Example

```python
from TidalPy.Structures import build_world, available_worlds

# Build one of the bundled example worlds by name.
print(available_worlds())            # 26 names, from 'charon' to 'triton' (see the WorldPack page)
earth = build_world("earth_simple")  # returns the Cython world (a BaseWorld subclass)

# build_world returns the world object directly, so its methods are immediate.
earth.solve_eos()
print(earth.surface_gravity_eos, earth.planet_mass_eos)

# Save a (possibly modified) world back out; schema_version is included on write.
earth.save_to_toml("earth_copy.toml")

# Build from a file path or an in-memory dict instead.
world = build_world("earth_copy.toml")
world = build_world(earth.get_config_dict())

# Change a few values of a source as it is built; nested tables merge key by key.
hot_earth = build_world(
    "earth_simple",
    overrides={"layers": {"mantle": {"temperature_k": 2000.0}}})
```

`build_world(source)` returns a `BaseWorld`: a `TerrestrialWorld`, `GasGiantWorld`, `StarWorld`, or, for the `layered` type, the `BaseWorld` class itself (see [Python API](#python-api)).

## Bundled Worlds (`WorldPack`)

The [Schema Examples](schema_examples.md) page has six commented files that together use every key of this schema. Example worlds ship in the package directory `TidalPy/WorldPack/`. On first use they are copied into a version-scoped, user-editable data directory ([`worldpack.md`](worldpack.md) covers the install and resolution mechanism):

```
<user documents>/TidalPy/<major>.<minor>.X/Worlds/
```

A world requested by bare name is read from the data directory before the package. Editing the installed TOML (_e.g._, `.../TidalPy/<major>.<minor>.X/Worlds/earth_simple.toml`) therefore changes what `build_world("earth_simple")` returns. A file is copied only if the data directory lacks a file of that name, so user edits are never overwritten and new packaged worlds appear on the next run. `install_worldpack(force=True)` re-copies the packaged versions and discards local edits.

```python
from TidalPy.Structures import available_worlds, install_worldpack
install_worldpack()          # copy packaged worlds into the data dir (copy-if-absent)
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
| `[prescribed_heating]` | optional | all | A layer's prescribed internal heating, keyed by layer name: `power_w` \[W\] spread over the layer by mass, or `specific_rate_w_kg` \[W kg$^{-1}$\], exactly one per layer (`BaseWorld.set_prescribed_heating`). It acts in a layer with `use_heating` during a thermal EOS solve. |
| `[luminosity]` | optional | `star` | The star's mass-to-luminosity model: `model` (`fixed`, `mass_to_luminosity`, or `power_law`) plus that model's parameters, as `Stellar.make_luminosity` takes them. Attaching it does not change the stored `luminosity_w`. |
| `moment_of_inertia_factor` | optional | all | $C/(MR^2)$ of the world's spin model, within $(0, 2/3]$. It gives the moment of inertia until the EOS is solved. Left out, it comes from `[worlds]` in `TidalPy_Configs.toml`: 0.4 (a uniform sphere), or 0.0754 for a star (an $n = 3$ polytrope). |
| `[layers.<name>]` | **yes** (optional for a star) | all | One table per layer (see below). |
| `[tides]` | optional | all | Tidal dissipation settings (see below). Omitted entirely, the world still gets a dissipation model from the `[tides]` defaults of `TidalPy_Configs.toml`. |

World `type` maps to a class as follows:

| `type` | Class |
|--------|-------|
| `terrestrial` | `TerrestrialWorld` |
| `gasgiant` | `GasGiantWorld` |
| `star` | `StarWorld` (layers optional) |
| `layered` | `BaseWorld` |

An omitted optional key falls back to the `[worlds]` block of `TidalPy_Configs.toml`, then to the C++ class default (see [Default Configuration Resolution](#default-configuration-resolution)). The loader never duplicates a default.

## Layer-Level Schema

Each world other than a star declares one or more `[layers.<layer_name>]` tables (a star may declare them too), keyed by the layer's name. Layers are ordered inner-to-outer by `layer_index` when given, otherwise by declaration order ([Geometry](#geometry) covers the radii). Each table builds one [`Layer`](../layers/layer.md), and every key the table leaves out takes the layer's default, which is the simplest case.

| Key | Required | Description |
|-----|----------|-------------|
| `layer_index` | optional | Inner-to-outer position (0 = innermost). Falls back to declaration order. |
| `radius_outer_m` | one-of | Outer radius \[m\] (absolute). |
| `radius_fraction` | one-of | Outer radius as a fraction of the world radius (`radius_outer_m = radius_fraction * world radius_m`). |
| `volume_fraction` | one-of | Layer shell volume as a fraction of the whole-world volume; the outer radius is solved from it. |
| `material` | optional | The layer's material: a MatPack name (`material = "peridotite"`) or a `[layers.<name>.material]` table (see [The Material Table](#the-material-table)). Absent takes `[layers] material` of `TidalPy_Configs.toml` (`simple_rock`). |
| `mass_kg` | optional | Layer mass \[kg\]. Defaults to 0.0; every successful EOS solve overwrites it with the solved layer mass. A layer that holds its mass holds this one. |
| `use_tides` | optional | Whether the layer dissipates tidal energy: off, it has no share in the quasi-homogeneous Love methods, and the radial solver treats it as elastic (its static moduli, no rheology), so it adds no heating. Default `true`. |
| `tidal_scale` | optional | The layer's share of the planet in the quasi-homogeneous Love methods (`homogeneous`, `cpl`, `ctl`) and of an analytic tide model's heating. Absent takes the layer's volume over the planet's. The radial solver resolves the layers directly and does not use it. |
| `is_volume_fixed` | optional | `false` makes the layer hold its mass rather than its volume: the EOS solve ends it where it encloses that mass, and the layers above it move with it. Default `true`. |
| `state` | optional | `"auto"` (the material decides, and a melting layer splits into solid and liquid zones), `"solid"`, or `"liquid"`. Default `"auto"`. |
| `is_static` | optional | Static approximation (no inertia) in the radial solve. Default `true`, so a liquid layer or liquid zone is a static liquid unless this is `false`. |
| `is_incompressible` | optional | Incompressible approximation in the radial solve. Default `false`. |
| `temperature_k` | optional | Layer temperature \[K\]. Default `0.0`, the cold rigid limit of the viscosity laws; a layer at 0 K takes no part in a thermal solve, and a solve warns once when it leaves the layer rigid. |
| `use_thermal_expansion` | optional | Let the material's density follow the temperature. Default `false`. |
| `use_melting` | optional | Let the material melt: its liquid phase, melt fraction, and weakening. Default `false`. |
| `use_pressure_melting` | optional | Let the material's melting curves follow the pressure (off: their zero-pressure values). Default `false`. |
| `use_melt_density` | optional | Mix the melt's density in by melt fraction (off: the solid's density). Default `false`. |
| `use_heating` | optional | Let the world's heat sources (the layer's radiogenics, the tides, and prescribed heating) act inside the layer during a thermal EOS solve. Default `false`. |
| `[layers.<name>.shear_rheology]`, `[layers.<name>.bulk_rheology]` | optional | Override the rheology the material's phase carries (see [Attached Physics Models](#attached-physics-models)). |
| `[layers.<name>.cooling]` | optional | The layer's cooling model. Absent leaves the layer isothermal. |
| `[layers.<name>.radiogenics]` | optional | The layer's radiogenics model. Absent gives it none. |

The [Layer](../layers/layer.md#physics-switches-and-flags) page describes what each switch does. The static moduli, the density laws, the thermal constants, and the viscosity and melting laws are material properties, not layer keys. They live in the layer's material, so the layer and the solve read the same numbers.

An unrecognized key, a model table without a `model` key, a switch that is not `true` or `false`, or a `state` outside the three values is a validation error. A retired key (`class`, `type`, `material_name`, `is_tidal`, `is_solid`, `use_thermal_eos`, the gas-layer scalars, and layer-level `eos`, `shear_viscosity`, `bulk_viscosity`, and `partial_melt` tables) is refused with a message naming what replaced it (`TidalPy.schema.RETIRED_LAYER_KEYS`).

### Geometry

Layers are built inner-to-outer, so the user never writes a layer's inner radius: it is the previous layer's outer radius (0 for the innermost). Supplying `radius_inner_m` raises an error. Each layer sets its outer radius with exactly one of `radius_outer_m`, `radius_fraction`, or `volume_fraction`. Supplying more than one, or none, is an error. With `volume_fraction` $f_V$ and world radius $R$, the shell volume is $f_V$ times the world volume, so $r_\mathrm{out} = \left(r_\mathrm{in}^3 + f_V R^3\right)^{1/3}$.

### Physical Values

After the keys, `validate_physical_values` checks the values. Every build runs it through `validate_world_config`. It raises a `ValueError` naming the world or layer, the key, the value found, and the allowed range when:

* the world's `radius_m` or `mass_kg` is not a positive finite number, `albedo` lies outside $[0, 1]$, `emissivity` outside $(0, 1]$, or the obliquity, spin rate, stellar temperature, or luminosity is not finite (the last two may be zero, which asks for the value to be derived);
* a `radius_fraction` or `volume_fraction` lies outside $(0, 1]$, or a `radius_outer_m` is not positive;
* a layer ends at or below the top of the layer under it, or above the world's radius;
* the outermost layer stops short of the world's radius, since the layers must fill the world (to a relative $10^{-6}$, which absorbs the roundoff of stacked fractions);
* a layer's `mass_kg`, `tidal_scale`, or `temperature_k` is negative or not finite;
* a `layer_index` is not a non-negative integer, or two layers resolve to the same index (a layer with none takes its position in the file).

### Attached Physics Models

A layer attaches a model through a nested table holding a `model` key and that model's parameters. The other keys are forwarded to the matching `make_*` factory as its parameter dict. Omitted keys keep the model's own defaults. A key the model does not read raises a `ValueError` naming the table and the closest accepted key.

| Model table | Factory |
|-------------|---------|
| `[layers.<name>.shear_rheology]` | `make_rheology` |
| `[layers.<name>.bulk_rheology]` | `make_rheology` |
| `[layers.<name>.cooling]` | `make_cooling` |
| `[layers.<name>.radiogenics]` | `make_radiogenics` |

An isotope radiogenics table that names neither a dataset (`isotopes`) nor its own isotope arrays takes the dataset of `[radiogenics] isotopes` in `TidalPy_Configs.toml`. See each module's documentation for the available model names and parameters.

### The Material Table

A layer's `material` is a MatPack name or a table, built by `TidalPy.Material.load_material` (see [MatPack](../../Material/matpack.md) and [Phases and Materials](../../Material/materials.md)). A table takes one of two forms:

- A preset with overrides: `preset` names a MatPack material, and the rest of the table merges over it. A law table that names the same model changes that law's keys, one that names another model replaces the law, and a law the preset lacks is added.
- A full definition: a `solid` phase table, a `liquid` phase table, or both, a `melting` table (`solidus`, `liquidus`, `weakening`, `bulk_modulus_mixing`, `bulk_viscosity_mixing`) joining them, and the material's `latent_heat_j_kg`.

A phase table holds its laws as sub-tables, each with its own `model` key: `eos` (the density law and its thermal terms), `shear_modulus`, `shear_viscosity`, `bulk_viscosity`, and the phase's default `shear_rheology` and `bulk_rheology`, plus its thermal constants (`thermal_conductivity_w_mk`, `heat_capacity_j_kgk`, their temperature exponents, and `thermal_reference_temperature_k`). A phase table may also name a `preset`, taking that material's phase in the same slot. A material with only a liquid phase is liquid everywhere (an ocean, a gas envelope). With both phases it melts between its solidus and liquidus when the layer has `use_melting` on, and equal curves give a single melting temperature.

```toml
[layers.mantle]
radius_fraction = 1.0
material = "peridotite"                 # a MatPack name

[layers.crust.material]                 # a preset with one law replaced
preset = "basalt"
[layers.crust.material.solid.shear_viscosity]
model = "constant"
reference_viscosity_pas = 1.0e23

[layers.ocean.material.liquid.eos]      # a full definition: a liquid-only material
model = "constant"
reference_density_kg_m3 = 1000.0
```

A radius-tabulated law (`model = "interpolate"`, with a `radius_m` table \[m\]) must span its layer to within 0.1% of the layer's outer radius at each end, or the build fails: the law holds its end values beyond the table, so a table in km would otherwise build a uniform layer. An error in the material names the layer's table (`[layers.<name>.material]`). The [Schema Examples](schema_examples.md) show every law with its keys.

## Tidal Dissipation (`[tides]`)

The optional world-level `[tides]` table sets how a world of any family, stars and gas giants included, dissipates tidal energy. Every key is optional: an omitted key falls back to the world family's `[tides.<type>]` table of `TidalPy_Configs.toml` when there is one, then to that file's `[tides]` block, then to a built-in default.

The per-degree lists of the `[tides]` block (`fixed_k = 0.3`, `fixed_q = 100` at degree 2) describe a planet. Stars take theirs from `[tides.star]`: the fluid Love numbers of an $n = 3$ polytrope ($k_2 = 0.0289$, a Sun-like star) with $Q = 1.93 \times 10^4$. This gives a modified quality factor $Q' = 3Q / (2k_2) = 10^6$, the usual assumption for stellar tides. The bundled `sol` and `trappist1` files state their tides explicitly. TRAPPIST-1 is fully convective and uses $n = 1.5$ values ($k_2 = 0.287$) at the same $Q'$.

| Key | Applies to | Description |
|-----|------------|-------------|
| `global_tidal_model` | all | Dissipation model: `rheology`, `cpl` (`fixed_q`), `ctl` (`fixed_dt`), or `ctl_q` (`fixed_dt_q`). Defaults per world family from `[tides.default_model]` in `TidalPy_Configs.toml`: `rheology` for terrestrial and layered worlds, `fixed_dt` for gas giants, `fixed_q` for stars. |
| `fixed_k` | all | Per-degree static potential Love numbers $k_l$, a list indexed from $l = 2$ (nine slots, $l = 2 \ldots 10$). Read only by the analytic models. A list shorter than nine zero-fills the remaining degrees, and a zero $k_l$ is no dissipation at that degree, so the list has to reach `max_degree_l`; the builder warns when a list the model reads stops short (the `short_degree_list` switch of `[warnings]`). |
| `fixed_q` | all | Per-degree tidal quality factors $Q_l$, same indexing. Read by `cpl` and `ctl_q`. |
| `fixed_dt_s` | all | Per-degree tidal time lags $\Delta t_l$ \[s\], same indexing. Read by `ctl` and `ctl_q`. The key carries the unit; the model alias stays `fixed_dt`. |
| `min_degree_l` | all | Lowest harmonic degree in the mode sum. Default `2`. |
| `max_degree_l` | all | Highest harmonic degree in the mode sum. Default `2`. |
| `eccentricity_trunc_lvl` | all | Eccentricity truncation level $N$: every product of two eccentricity functions, and so the heating, is kept through $e^N$ (see [Eccentricity Functions](../../Tides/Eccentricity.md)). Tabulated at 2, 4, 6, 8, 10, 20, and 50, or `"exact"` for the functions from the exact orbit (any $e < 1$); default `10`. An untabulated level is promoted to the next tabulated one with a once-per-session warning, so accuracy never drops silently. `eccentricity_truncation` is accepted as an alias. |
| `eccentricity_exact_tolerance` | all | For `eccentricity_trunc_lvl = "exact"`: the modes kept leave a $q^2$-weighted tail of the squared eccentricity functions below this fraction of the total, which bounds the relative error of the heating. In (0, 1); default `1e-4`. Ignored by the tabulated levels. |
| `obliquity_trunc_lvl` | all | Obliquity truncation level $N$: every product of two obliquity functions (the heating) is kept through $I^N$. Tabulated at 0, 2, and 4; default `"off"`. `"off"` means 0 (no obliquity terms), and `"gen"` or `"general"` the general functions (exact at any obliquity). Level 2 stays within 1% of the general heating to $I \approx 8^\circ$, level 4 to $27^\circ$ ([Obliquity Functions](../../Tides/Obliquity.md)). Untabulated integers are promoted like the eccentricity levels (1 to 2, and anything past 4 to `"gen"`). `obliquity_truncation` is accepted as an alias. |
| `layer_tidal_heating` | all | Whether `calc_tides` also resolves each layer's heating when the Love numbers come from the radial solver, a volume integral of the radial solution that costs about as much as the global solve again. The other paths share out the heating at no extra cost. Default `true`. |
| `love_method` | all | How the Love numbers are obtained: `radial_solver` (aliases `shooting`, `rs`; the default), `propagation_matrix` (`prop_matrix`, `pm`, `prop`), `homogeneous` (`homogen`), `cpl`, `ctl`, or `laterally_inhomogeneous` (`3d`, `lat_inhom`, reserved for the 3D solver: setting it raises `NotImplementedError`). The three homogeneous methods use the analytic homogeneous-sphere formulas instead of a radial solve, so they have no depth-resolved solution and the 3D stress, strain, and heating path raises `RuntimeError` while one of them is configured. |
| `love_fixed_q` | all | Scalar $Q$ the `cpl` Love method applies to the static Love numbers. Unset by default, in which case the tide model's own `fixed_q` is used. |
| `love_fixed_dt_s` | all | Scalar time lag \[s\] the `ctl` Love method applies. Unset by default, falling back to the tide model's `fixed_dt_s`. |

`global_tidal_model` (Love numbers to dissipation) and `love_method` (how the Love numbers are computed) are independent: a world can solve its Love numbers with the radial solver and still collapse them with an analytic tide model.

```toml
[tides]
global_tidal_model = "rheology"      # dissipation from the layers' complex moduli
min_degree_l = 2
max_degree_l = 3                     # include the degree-3 tide
eccentricity_trunc_lvl = 20          # heating through e^20; 10 is the default
obliquity_trunc_lvl = "gen"          # exact obliquity terms
love_method = "radial_solver"
```

A star or gas giant has no interior to solve and uses an analytic model with per-degree parameters:

```toml
[tides]
global_tidal_model = "fixed_q"
fixed_k = [0.3, 0.15, 0.1]           # l = 2, 3, 4; higher degrees would be zero-filled
fixed_q = [1.0e5, 1.0e5, 1.0e5]
```

`world.get_tide_config()` returns the degree and truncation settings under these key names. `get_config_dict()` writes them in a `[tides]` table with the tide model's own parameters (the model's name as `global_tidal_model`), so a world's tidal configuration survives a save and rebuild. The round trip writes the resolved `obliquity_trunc_lvl`: a world written with `"off"` reads back as `0`, and `"general"` as `"gen"`.

## Solver Settings (`[eos_solver]`, `[radial_solver]`)

A world's file may pin the solver settings its results depend on, so the file and a TidalPy configuration file reproduce a run on another machine. The two tables take the keys of the same-named sections of `TidalPy_Configs.toml` (see [Configurations](../../Overview/2_TidalPy_Configurations.md)):

- `[eos_solver]`: `integration_method`, `rtol`, `atol`, `pressure_tol`, `max_iters`, `slices_per_layer`, `nondimensionalize`, `solve_temperature`, `max_thermal_passes`, and `thermal_tol`.
- `[radial_solver]`: `integration_method`, `rtol`, `atol`, `use_kamata`, `start_radius_tolerance`, `scale_rtols`, `max_num_steps`, `expected_size`, `max_ram_mb`, and `nondimensionalize`.

A pinned key overrides the configuration for every solve the world runs (`solve_eos`, `solve_love_numbers`, `calc_tides`, and the 3D paths). A call's own argument still overrides the pinned key. A key left out follows the configuration.

```toml
[eos_solver]
integration_method = "RK45"
rtol = 1.0e-8

[radial_solver]
use_kamata = true
```

`world.set_solver_defaults(eos_solver=..., radial_solver=...)` pins the same keys on a built world. `world.get_solver_defaults()` returns the pinned tables. `get_config_dict()` carries them, so they survive a save and rebuild. The tables hold no physical parameters: they change how a result is computed, not what is computed. Every bundled layered world except `earth_prem` and `earth_prem_q` pins `solve_temperature = false` in its `[eos_solver]` table. The two PREM worlds are built from a profile and take `integration_method = "RK45"` (see [Building a World From a Radial Profile](#building-a-world-from-a-radial-profile)). `luna_dynamic` also pins `rtol = 1.0e-8` in a `[radial_solver]` table so that its dynamic liquid core keeps its yearly Love number.

## Default Configuration Resolution

A value a world file leaves out is resolved in this order:

1. The user world (dict or TOML): a value the user writes takes precedence over every other tier.
2. The TidalPy configuration (`TidalPy_Configs.toml`): a world-level property comes from its `[worlds]` block, a `[tides]` value from its `[tides]` block, and the material of a layer that names none from `[layers] material`.
3. The constructor or factory default: anything still unset falls through to the C++ or Cython default. A layer's switches and flags take the `Layer` defaults, a model table's omitted keys take the model's own defaults, and a layer gets no cooling or radiogenics model unless its table names one.

A layer's material values come from the material it names: a MatPack material's values are edited in its file in the `Materials` folder of the TidalPy data directory (`TidalPy.Material.material_info(name)["path"]`), or overridden in the layer's own material table.

World-level properties resolve through the `[worlds]` block: a world's `albedo` comes from its own table, else `[worlds].albedo`, else the class default. A `[worlds.<type>]` sub-table specializes the block for one world type. It keeps the star-only `effective_temperature_k` and `luminosity_w` off every other world.

`TidalPy_Configs.toml` is TidalPy's configuration file (see [TidalPy Configurations](../../Overview/2_TidalPy_Configurations.md)). It is generated from `TidalPy.defaultc` into the user's TidalPy `Config` directory on first use and is then user-editable. A new default belongs in `defaultc.py`. Its `[numerical]` section also feeds the shared C++ config singleton read by every compiled module (see [Constants](../../Utilities/constants.md)):

- the frequency, modulus, and thickness floors;
- `numerical_floor`: the magnitude a guarded denominator is raised to;
- `layer_continuity_rtol`: how closely a layer's inner radius must match the previous layer's outer radius;
- `minimum_solid_rigidity` and `minimum_zone_fraction`: where a melting layer turns liquid, and the thinnest zone the radial solver takes as a layer of its own.

## Example: Three-Layer Terrestrial World

Each layer of this world names its geometry, its material, and its temperature, and the outer core is a liquid-only MatPack material, so it is a liquid layer. It is the skeleton of the bundled `earth_simple`, which gives full material tables with fitted densities, viscosity laws, and mantle rigidity, and explains them in its file comments.

```toml
schema_version = "0.2.0"
name = "Earth-Simple"
type = "terrestrial"
radius_m = 6371000.0
mass_kg = 5.972e24
spin_frequency_rad_s = 7.292e-5

[layers.inner_core]
layer_index = 0
radius_outer_m = 1221500.0   # inner radius is derived (0 for the innermost)
temperature_k = 5500.0
material = "iron"

[layers.outer_core]
layer_index = 1
radius_outer_m = 3480000.0
use_tides = false            # a liquid that dissipates nothing in this model
temperature_k = 4500.0
material = "liquid_iron"     # liquid-only, so a static liquid layer

[layers.mantle]
layer_index = 2
radius_fraction = 1.0        # outer radius = full world radius; inner = the outer core's outer
temperature_k = 1600.0
material = "lower_mantle"
```

Adding a key or sub-table overrides a default. For example, this replaces the mantle's material with a MatPack preset whose shear viscosity is replaced, and turns its melting on:

```toml
[layers.mantle]
layer_index = 2
radius_fraction = 1.0
temperature_k = 1600.0
use_melting = true
use_pressure_melting = true

[layers.mantle.material]
preset = "peridotite"

[layers.mantle.material.solid.shear_viscosity]
model = "constant"
reference_viscosity_pas = 1.0e21
```

A world whose tide model is analytic (`fixed_q`, `fixed_dt`, or `ctl_q`, set by `global_tidal_model` or the type's default) needs no layers, which is how a star or a simple gas giant is usually given; the `rheology` model solves the interior's Love numbers, so a world on it needs at least one. A star:

```toml
schema_version = "0.2.0"
name = "Sol"
type = "star"
radius_m = 695700000.0
mass_kg = 1.988435e30
effective_temperature_k = 5772.0
```

## Building a World From a Radial Profile

A world can describe its interior with a PREM-like radial profile instead of layer tables, _e.g._, to compare with published data or to pass the output of a more sophisticated EOS solver to TidalPy's Love number calculations. A top-level `data_file` key names the delimited profile file.

```toml
schema_version = "0.2.0"
name = "Earth-PREM"
type = "terrestrial"
radius_m = 6371000.0
mass_kg = 5.972e24
data_file = "PREM.csv"  # Can be a path to a file too.
```

From Python, pass the arrays to `build_world` under a `data` key instead of a data file. A world gives its profile one way or the other, never both.

```python
import numpy as np

radius_km = np.array([0.0, 3480.0, 3480.0, 6371.0])        # The boundary row repeated: a liquid core, a solid mantle
density = np.array([12000.0, 10000.0, 5500.0, 3300.0])     # [kg m-3]
vp = np.array([11000.0, 8000.0, 13700.0, 8000.0])          # [m s-1]
vs = np.array([0.0, 0.0, 7300.0, 4500.0])                  # [m s-1]; zero is liquid

world = build_world({
    "schema_version": "0.2.0",
    "name": "My-Earth",
    "type": "terrestrial",
    "radius_m": 6371000.0,
    "mass_kg": 5.972e24,
    "data": {"radius_km": radius_km, "density": density, "vp": vp, "vs": vs},
})                                                         # Two layers, detected at the solid-liquid transition
print([layer.name for layer in world], world.layer_0.is_liquid)   # ['layer_0', 'layer_1'] True
```

### Profile Format

A profile is a delimited table (comma, semicolon, tab, or whitespace; `#` comment lines ignored) or a mapping of arrays with the same names. It must give a radius, a density, and either both seismic velocities or both static moduli:

| Quantity | Accepted names | Units | Required |
|----------|----------------|-------|----------|
| radius | `radius`, `r`, `rad`, `radii` | m or km | yes, or a depth |
| depth | `depth`, `z` | m or km | converted with the world's `radius_m` |
| density | `density`, `rho`, `dens` | kg/m³ | yes |
| P-wave velocity | `vp`, `v_p`, `p_velocity`, … | m/s or km/s | yes, unless moduli are given |
| S-wave velocity | `vs`, `v_s`, `shear_velocity`, … | m/s or km/s | yes, unless moduli are given |
| shear modulus | `shear_modulus`, `mu`, `rigidity` | Pa | instead of the velocities |
| bulk modulus | `bulk_modulus`, `k`, `incompressibility` | Pa | instead of the velocities |
| shear viscosity | `shear_viscosity`, `eta_shear`, `viscosity`, `visc` | Pa s | no |
| bulk viscosity | `bulk_viscosity`, `eta_bulk`, `zeta` | Pa s | no |
| shear quality factor $Q_\mu$ | `q_mu`, `qmu`, `q_shear`, `q_s`, `q_beta` | none | no; used only with `q_provided = true` |
| bulk quality factor $Q_\kappa$ | `q_kappa`, `qkappa`, `q_bulk`, `q_k` | none | no; used only with `q_provided = true` |

Names are matched ignoring case and punctuation, so `Vp`, `V_P`, and `vp` are one name. A name may state its unit: `radius_km`, `Vp [km/s]`, `rho_kg_m3`. A name whose unit the reader does not convert for its quantity is refused with a `ValueError` listing the units it accepts. Columns are found by name, so their order does not matter, and a column whose name matches none of these is not read. A bare `eta` is not a viscosity: in PREM and the IRIS tables it is the dimensionless anisotropy parameter, so an `eta` column is not read. The names come from a header row or from the last `#` comment line before the data. A header row must name every column, or the file is refused. A `#` comment line is taken as the header only when it names every column, so prose about the data is not mistaken for a header. In a whitespace-delimited header a bracketed unit stays with its name (`depth (km)`). A file with no header is read positionally as radius, density, `Vp`, `Vs`, shear viscosity, bulk viscosity, so quality factors need a header.

A radius or depth column with no stated unit is read as kilometers when its largest value is below 100000, and as meters otherwise. A density below 100 kg/m³, a velocity below 100 m/s, and a bulk modulus or a solid's shear modulus below 10⁶ Pa are refused as a column in g/cm³, km/s, or GPa that did not state its unit (name it `density_g_cm3`, `vp_km_s`, or `shear_modulus_gpa`). The file may be ordered surface-first or center-first (it is sorted internally). Where the velocities are given, the static moduli are derived per row: shear $\mu = \rho V_s^2$, bulk $K = \rho \left(V_p^2 - \tfrac{4}{3} V_s^2\right)$.

### Layer Detection

The profile is scanned from the center outward and split into layers by shear modulus: `Vs = 0` is liquid, non-zero is solid, and every solid-liquid transition starts a new layer. Layers are named `layer_0`, `layer_1`, and so on, inner to outer. Duplicate-radius boundary points are absorbed, so no zero-thickness layer is produced. At a duplicated boundary radius the lower layer's row comes first, whichever way the file is ordered. The bundled `PREM.csv` replaces PREM's 3 km ocean with the upper crust and yields three layers: inner core (solid), outer core (liquid), and mantle plus crust (solid).

Each layer's slice of the profile becomes its material, a phase of radius-tabulated (`interpolate`) laws: the density and bulk modulus in its `eos`, the shear modulus of a solid layer, and any viscosities the profile gave. A solid layer's material is a solid phase, and a liquid layer's is a liquid-only material, so the layer is a liquid (static, unless its table sets `is_static = false`). The material holds those arrays and is the only place a radial grid persists. A layer starts where the one below it ends. When the profile repeats no row at that boundary, the layer's first row is repeated at the boundary radius, so its table spans the layer. A solid layer takes part in the tides (`use_tides`) and a liquid one does not. No other law and no melting is added: a profile that names no viscosity produces an elastic layer, with no viscosity law, no melting, and no rheology it did not ask for.

### Refining Detected Layers

A profile fixes the number of layers, their boundaries, and their materials. It holds no complex moduli, so a rheology, a cooling model, and radiogenics are still named in a `[layers.<name>]` table.

A table picks the detected layer it refines with `layer_index`, or by being named `layer_N`. Only the refined layers need a table, and the table gives the layer its name.

- A given outer radius (`radius_outer_m` or `radius_fraction`) must match the detected boundary.
- Two tables may not claim the same layer.
- An index outside the detected range raises an error.
- A `material` table merges over the layer's slice of the profile: a law naming another model replaces the profile's (a constant viscosity, say), and a law the profile lacks is added. A MatPack name is refused, since the profile is the layer's material.
- A single number given for a tabulated quantity of the profile's laws (_e.g._, `bulk_modulus_pa = 1.0e11` in the `solid.eos` table) holds across the layer: TOML overrides the data file.
- Every other key (a switch, a flag, a model table) is set on the layer.
- In a world that sets `q_provided = true`, a table may not give a viscosity law or a `preset`, and names no rheology but `seismic_q` (or `elastic` for the bulk). A solid layer's table may not turn `use_melting` on or name the `convection` cooling model. See below.

```toml
[layers.mantle]             # refines the outermost detected layer; the other two need no table
layer_index = 2
[layers.mantle.shear_rheology]
model = "maxwell"
[layers.mantle.material.solid.shear_viscosity]
model = "constant"          # the profile named no viscosity, so give one here
reference_viscosity_pas = 1.0e21
```

### Quality Factors in Place of Viscosities

Seismic profiles give quality factors, not viscosities: PREM tabulates $Q_\mu$ and $Q_\kappa$ at its 1 s reference period. Setting `q_provided = true` builds each solid layer's loss from them with the [`seismic_q` rheology](../../Rheology/rheology_models.md), with no viscosity involved. The complex modulus is the profile's static modulus with $Q(\omega) = Q_\mathrm{ref}(\omega/\omega_\mathrm{ref})^{a}$ and the dispersion that loss implies.

```toml
data_file = "PREM.csv"                             # carries q_mu and q_kappa columns
q_provided = true                                  # default false: the Q columns are read but not used
q_reference_frequency_rad_s = 6.283185307179586    # where the Q and the moduli were measured; default 2 pi (1 s)
q_frequency_exponent = 0.0                         # a in Q ~ omega^a, in [0, 1); default 0
```

| Key | Default | Meaning |
|-----|---------|---------|
| `q_provided` | `false` | Use the profile's quality factors. The profile must then give `q_mu` and must not give viscosities. |
| `q_reference_frequency_rad_s` | $2\pi$ | Frequency \[rad s⁻¹\] at which the profile's $Q$ and moduli were measured. |
| `q_frequency_exponent` | 0.0 | Exponent $a$ of $Q \propto \omega^{a}$. $a = 0$ keeps the seismic $Q$ at tidal periods; $a > 0$ lowers it there. |

These keys belong only to a world with a profile. The last two require `q_provided = true`, so a setting cannot be silently ignored.

- **Solid layers.** Each gets `shear_rheology = {model = "seismic_q", ...}` with the world's two settings, and a `bulk_rheology` of the same kind when the profile gives `q_kappa`. The quality factors ride in the layer's viscosity arrays, which is where `seismic_q` reads them. They must be positive in every solid row.
- **Liquid layers.** These (a liquid's $Q_\mu$ is conventionally 0) take no quality factor and no rheology.
- **Refinement tables.** A `[layers.<name>]` table may name `seismic_q` with its own `reference_frequency_rad_s` or `q_frequency_exponent` to override the world's for that layer, or set an `elastic` bulk rheology to ignore $Q_\kappa$. It may not name another rheology, give a viscosity law, or name a material `preset`: each would put a viscosity where `seismic_q` reads a quality factor. On a solid layer it may not turn `use_melting` on or name the `convection` cooling model: the melt weakening and the Rayleigh number both read the layer's viscosity, which there holds $Q_\mu$.

The bundled `earth_prem_q` world is PREM with its own quality factors. At 1 s its $k_2$ about equals the elastic PREM's (0.298351 against 0.298362). At the M2 tide the dispersion raises it from 0.298 to 0.302, with $|k_2|/|\mathrm{Im}\,k_2| \approx 500$; with `q_frequency_exponent = 0.15` that ratio falls to about 100. See [`example_profile_q_world.toml`](examples/example_profile_q_world.toml).

The `data_file` path is resolved relative to the world TOML's directory, then the worlds data directory, then the packaged `WorldPack` (see [`worldpack.md`](worldpack.md)). A world built this way pins `integration_method = "RK45"` in its `[eos_solver]` table unless the file sets that key, and `get_solver_defaults()` reports it. On an interpolated profile RK45 is about 2.8 times faster than DOP853 at equal accuracy, because the profile's kinks defeat the higher order.

The world keeps the configuration as given (the `data_file` reference as written and the refining tables) on `portable_config`. While the world is unchanged since its build, `save_to_toml` writes that form, with the reference rewritten to find the same file from the folder saved into, so the saved file holds no expanded profile. A world changed after its build is saved as its live `get_config_dict()`, the expanded profile included, so the change is kept. `source_config` holds the expanded form as built.

## System Schema

A system TOML, built with `build_system`, groups several worlds and the orbits that connect them. The file carries a top-level `schema_version` and `name`, then one `[worlds.<key>]` table per member world. The table key is that world's name within the system.

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
from TidalPy.Structures import build_system
system = build_system("sol_system")     # or a path / a config dict
```

[`../system/system.md`](../system/system.md) documents the full system API (evolution, insolation, save/load).

## Python API

Every entry point below is exported from `TidalPy.Structures.configs`, and the main ones (`build_world`, `build_layer_from_dict`, `build_system`, `load_world`, `load_system`, `available_worlds`, `available_systems`, `save_world_to_toml`, `install_worldpack`, `SCHEMA_VERSION`) from `TidalPy.Structures` too, together with the classes (`Layer`, `BaseWorld`, `TerrestrialWorld`, `GasGiantWorld`, `StarWorld`, `System`) and `make_tide`.

### High Level

* `build_world(source, overrides=None, force=False) -> BaseWorld`: resolve `source` (bundled name, file path, or dict), validate it, and return the built Cython world. `overrides` is a nested dict merged over the configuration before the build (`TidalPy.configurations.merge_configs`: tables merge key by key, any other value replaces the source's), so it only needs the keys it changes. `force=True` bypasses the schema-version warning. A path to a binary file (`save_binary`) is loaded with `load_world` instead, and refuses `overrides`. A thin wrapper over `BaseWorld.build(source, overrides=None, force=False)`, which returns the type-appropriate subclass. `build_system(source, overrides=None, force=False)` is the same for a system.
* `load_world(path, force=False) -> BaseWorld` and `load_system(path, force=False) -> System`: read a binary file into a new object of the class that saved it, with no placeholder object. Each raises `IOError` naming what a file holds when it is not of its kind (`"it is a System file; load it with load_system"`, or not a TidalPy binary file). `TidalPy.Utilities.binary.binary_file_class(path)` reads only the header and returns the class name, or `None` for a file that is not a TidalPy binary file.
* An unknown bundled name raises `FileNotFoundError` naming the closest bundled name and listing them all, and an unknown world, `[tides]`, layer, solver-table, or system key raises `ValueError` naming the closest accepted key (`"did you mean 'temperature_k'?"`).
* `load_radial_data(source, surface_radius=None) -> dict`: read a radial profile (a data-file path or a mapping of arrays) into MKS arrays ascending in radius, the same reader `data_file` and `data` worlds use. `detect_layer_boundaries(radius, shear_modulus)` returns the `(start, end, is_solid)` runs it splits into.
* `world.save_to_toml(path, overwrite=True)`: write the world as it is now (`world.get_save_config(destination_dir)`), stamped with the current `schema_version` under a comment header naming the TidalPy, SciPy, and CyRK versions that wrote it. That is the live `get_config_dict()`, validated against this schema first so it writes a buildable file or raises `ValueError`, except for a world built from a `data_file` and unchanged since its build, which writes its `portable_config` with the file reference rewritten for the destination folder. The layer masses an EOS solve sets and the world's name do not count as changes. `path` may be a string or a `pathlib.Path`.
* `world.get_config_dict()`: the live world as a builder-valid table (`type`, name-keyed `layers` with each layer's scalars, `material` table, and model sub-tables, `tides`, `schema_version`).
* `build_world(world.get_config_dict())`, `build_system(system.get_config_dict())`, and `build_layer_from_dict(layer.get_config_dict()) -> Layer`: rebuild an object from the dictionary its `get_config_dict()` returns (see [Round Trip](#round-trip)), leaving the dict unmodified. `build_layer_from_dict` takes only a `dict` and raises `TypeError` for anything else.
* `world.config` (alias of `world.source_config`): the normalized configuration dict the world was built from (`None` if constructed directly). `world.portable_config`: for a world built from a `data_file`, the configuration as given (`None` otherwise). `world.built_config`: the world's `get_config_dict()` at the end of its build, which `save_to_toml` compares the live state against. A successful `load_binary` clears all three, since they describe the world before the load, and `save_to_toml` then writes `get_config_dict()`.
* `EOS_SOLVER_KEYS`, `RADIAL_SOLVER_KEYS`, and `validate_solver_table(section, table, where)`: the keys a world's `[eos_solver]` and `[radial_solver]` tables may pin, and the check the loader and `set_solver_defaults` apply to them.
* `available_worlds() -> list[str]`: names of the bundled example worlds (data dir unioned with packaged `WorldPack`). Bundled system files share that directory and are listed by `available_systems()` instead. `build_world` on a system config raises a `ValueError` naming `build_system`, and the reverse holds too.
* `install_worldpack(force=False) -> str`: copy the packaged `WorldPack` worlds into the user data directory (copy-if-absent unless `force`); returns that directory.

### Low Level

* `load_build_config(source, overrides=None, force=False) -> dict`: the configuration a build uses, the step `build_world` and `build_system` share: the source loaded, its schema version graded, the overrides merged, and the structural defaults filled in.
* `save_world_to_toml(config, path, overwrite=True)`: serialize a config dict.

### Loader / Validation

The schema's key sets (`WORLD_TYPES`, `LAYER_SCALAR_KEYS`, `LAYER_MODEL_SECTIONS`, `LAYER_STATES`, `RETIRED_LAYER_KEYS`, `ALLOWED_TIDES_KEYS`, `EOS_SOLVER_KEYS`, and the rest) are defined in `TidalPy.schema`. That module has no TidalPy imports, so the configuration loader can check `TidalPy_Configs.toml` against them while `import TidalPy` runs. The loader below re-exports them.

From `TidalPy.Structures.configs.toml_loader`:

* `SCHEMA_VERSION`: the current schema version string (`"0.2.0"`).
* `load_toml(source)`: parse a file path or pass through a dict.
* `validate_schema_version(config, force=False)`: the graded schema check described at the top of this page. `force=True` bypasses it.
* `validate_world_config(config)` / `validate_layer_config(name, cfg)`: structural validation. `validate_world_config` ends with `validate_physical_values(config)`, the value checks of [Physical Values](#physical-values).
* `merge_with_defaults(config)`: apply structural (non-physical) defaults.

## Round Trip

`save_to_toml` writes the world as it is now, so a `build -> change -> save_to_toml -> build` cycle reproduces the changed world.

```python
world = build_world("earth_simple")
world.save_to_toml("earth_copy.toml")
reloaded = build_world("earth_copy.toml")   # identical structure
```

The retained configuration is the world as its file described it. `get_config_dict()` is the world as it stands now, including every change made since the build. `build_world` and `build_system` rebuild the same class with the same parameters and attached models from that dict, and `build_layer_from_dict` does the same for a standalone layer. A world constructed directly in Python with no tide model attached writes no `tides` table, so its rebuild takes the builder's default tide model for its type (see [Tidal Dissipation](#tidal-dissipation-tides)):

```python
from TidalPy.Structures import build_layer_from_dict

world.set_spin_frequency(2.0e-5)                 # A change made after the build
config = world.get_config_dict()                 # Builder-valid: type, layers by name, tides, schema_version
twin = build_world(config)                       # Same class, same parameters, same models
assert twin.get_config_dict() == config

mantle = world.mantle                            # The layer named "mantle"
mantle_twin = build_layer_from_dict(mantle.get_config_dict())    # A standalone layer, owned by no world

system_twin = build_system(system.get_config_dict())             # Worlds, tidal hosts, star, and orbits
```

The live dict carries every layer scalar, the layer's full material table, its rheology overrides, and its cooling and radiogenics tables explicitly, so a rebuild takes nothing from `[layers] material` or a MatPack file that may have changed since. A standalone layer's dict adds the keys a world would supply from the layer's place in its `layers` table (`name` and `radius_inner_m`). A world drops those keys when it nests the layer. A system's dict inlines each member world's live dict under `world`. Solved state (the EOS, Love numbers, tides) is not configuration, so run the solves again on a rebuilt object.
