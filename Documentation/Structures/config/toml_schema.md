# World Configuration & TOML Schema (`Structures.configs`)

_Updated: 2026-10-08_

Schema version `0.2.0`.

A world file (TOML) or the equivalent `dict` describes a world: its type, its layers, their materials, and their physics models. `build_world` builds it and `save_to_toml` writes a world back out. This page lists every key of a world file and a system file. The [Schema Examples](schema_examples.md) use every key, and [WorldPack](worldpack.md) covers the bundled worlds.

## Example

```python
from TidalPy.Structures import build_world, available_worlds

print(available_worlds())            # 26 bundled names, 'charon' to 'triton'
earth = build_world("earth_simple")  # a BaseWorld subclass
earth.solve_eos()
print(earth.surface_gravity_eos, earth.planet_mass_eos)

earth.save_to_toml("earth_copy.toml")
world = build_world("earth_copy.toml")         # from a file
world = build_world(earth.get_config_dict())   # from a dict

# Change values as it is built; nested tables merge key by key.
hot_earth = build_world(
    "earth_simple",
    overrides={"layers": {"mantle": {"temperature_k": 2000.0}}})
```

A bundled name is read from its editable copy in `<user documents>/TidalPy/<major>.<minor>.X/Worlds/` first, so editing that copy changes what `build_world("earth_simple")` returns ([Name Resolution](worldpack.md#name-resolution)).

## Example: Three-Layer Terrestrial World

The skeleton of the bundled `earth_simple`, whose file adds fitted material tables:

```toml
schema_version = "0.2.0"
name = "Earth-Simple"
type = "terrestrial"
radius_m = 6371000.0
mass_kg = 5.972e24
spin_frequency_rad_s = 7.292e-5

[layers.inner_core]
layer_index = 0
radius_outer_m = 1221500.0   # the inner radius is never given
temperature_k = 5500.0
material = "iron"

[layers.outer_core]
layer_index = 1
radius_outer_m = 3480000.0
use_tides = false
temperature_k = 4500.0
material = "liquid_iron"     # liquid-only: a static liquid layer

[layers.mantle]
layer_index = 2
radius_fraction = 1.0        # reaches the surface
temperature_k = 1600.0
material = "lower_mantle"
```

Any other key or sub-table overrides a default (_e.g._, `use_melting = true`, or a [material table](#the-material-table)). A world on an analytic tide model (`fixed_q`, `fixed_dt`, or `ctl_q`, given or by default) needs no layers, as for most stars and simple gas giants; the `rheology` model needs at least one. A star:

```toml
schema_version = "0.2.0"
name = "Sol"
type = "star"
radius_m = 695700000.0
mass_kg = 1.988435e30
effective_temperature_k = 5772.0
```

## World-Level Schema

| Key | Required | Description |
|-----|----------|-------------|
| `schema_version` | optional | Schema marker (e.g. `"0.2.0"`). Checked when present (see [Schema Version](#schema-version)); written on save. |
| `name` | **yes** | World name. |
| `type` | **yes** | `terrestrial` (`TerrestrialWorld`), `gasgiant` (`GasGiantWorld`), `star` (`StarWorld`), or `layered` (`BaseWorld`). |
| `radius_m` | **yes** | World radius \[m\]. |
| `mass_kg` | **yes** | World mass \[kg\]. |
| `albedo` | optional | Bond albedo. |
| `emissivity` | optional | Surface emissivity. |
| `obliquity_rad` | optional | Axial obliquity \[rad\]. |
| `spin_frequency_rad_s` | optional | Rotation rate \[rad/s\]. |
| `effective_temperature_k` | optional | Star only: effective temperature \[K\]. |
| `luminosity_w` | optional | Star only: luminosity \[W\]. |
| `[luminosity]` | optional | Star only: the mass-to-luminosity model, `model` (`fixed`, `mass_to_luminosity`, or `power_law`) plus its parameters, as `Stellar.make_luminosity` takes them. It does not change the stored `luminosity_w`. |
| `[prescribed_heating]` | optional | Heating keyed by layer name, exactly one of `power_w` \[W\] (spread by mass) or `specific_rate_w_kg` \[W kg$^{-1}$\] per layer (`BaseWorld.set_prescribed_heating`). It acts in a layer with `use_heating` during a thermal EOS solve. |
| `moment_of_inertia_factor` | optional | $C/(MR^2)$ of the spin model, in $(0, 2/3]$; the moment of inertia until the EOS is solved. Default from `[worlds]` in `TidalPy_Configs.toml`: 0.4 (uniform sphere), or 0.0754 for a star ($n = 3$ polytrope). |
| `[layers.<name>]` | **yes**, unless the tide model is analytic | One table per layer (see [Layer-Level Schema](#layer-level-schema)). |
| `[tides]` | optional | Tidal dissipation (see [Tidal Dissipation](#tidal-dissipation-tides)). Omitted, the world gets its type's default model. |
| `[eos_solver]`, `[radial_solver]` | optional | Pinned solver settings (see [Solver Settings](toml_schema.md#solver-settings-eos_solver-radial_solver)). |
| `data_file`, `data` | optional | A radial profile in place of layer tables (see [Building a World From a Radial Profile](#building-a-world-from-a-radial-profile)). |

An omitted optional key falls back to `TidalPy_Configs.toml`, then to the class default (see [Default Configuration Resolution](#default-configuration-resolution)).

## Layer-Level Schema

Each `[layers.<layer_name>]` table builds one [`Layer`](../layers/layer.md), named by its key. Layers are ordered inner-to-outer by `layer_index`, or by declaration order. An omitted key takes the layer's default.

| Key | Required | Description |
|-----|----------|-------------|
| `layer_index` | optional | Position from the center (0 = innermost). Default: declaration order. |
| `radius_outer_m` | one-of | Outer radius \[m\]. |
| `radius_fraction` | one-of | Outer radius over the world radius. |
| `volume_fraction` | one-of | Shell volume over the world volume (see [Geometry](#geometry)). |
| `material` | optional | A MatPack name (`material = "peridotite"`) or a table (see [The Material Table](#the-material-table)). Default: `[layers] material` of `TidalPy_Configs.toml` (`simple_rock`). |
| `mass_kg` | optional | Layer mass \[kg\]. Default 0.0; each successful EOS solve overwrites it. A layer that holds its mass holds this one. |
| `use_tides` | optional | Whether the layer dissipates. Off, it has no share in the quasi-homogeneous Love methods and is elastic (static moduli) in the radial solver, so it adds no heating. Default `true`. |
| `tidal_scale` | optional | The layer's share in the quasi-homogeneous Love methods (`homogeneous`, `cpl`, `ctl`) and of an analytic tide model's heating. Default: its volume over the planet's. Unused by the radial solver. |
| `is_volume_fixed` | optional | `false` makes the layer hold its mass, not its volume: the EOS solve ends it where it encloses that mass, and the layers above move with it. Default `true`. |
| `state` | optional | `"auto"` (the material decides; a melting layer splits into solid and liquid zones), `"solid"`, or `"liquid"`. Default `"auto"`. |
| `is_static` | optional | No inertia in the radial solve. Default `true`, so a liquid layer or zone is static unless this is `false`. |
| `is_incompressible` | optional | Incompressible radial solve. Default `false`. |
| `temperature_k` | optional | Temperature \[K\]. Default `0.0`, the cold rigid limit of the viscosity laws; a 0 K layer takes no part in a thermal solve, and a solve warns once when it leaves the layer rigid. |
| `use_thermal_expansion` | optional | Density follows temperature. Default `false`. |
| `use_melting` | optional | The material can melt (liquid phase, melt fraction, weakening). Default `false`. |
| `use_pressure_melting` | optional | Melting curves follow pressure (off: zero-pressure values). Default `false`. |
| `use_melt_density` | optional | Mix the melt's density in by melt fraction (off: the solid's). Default `false`. |
| `use_heating` | optional | The world's heat sources (radiogenics, tides, prescribed heating) act in the layer during a thermal EOS solve. Default `false`. |
| `[layers.<name>.shear_rheology]`, `[layers.<name>.bulk_rheology]` | optional | Override the material phase's rheology (see [Attached Physics Models](#attached-physics-models)). |
| `[layers.<name>.cooling]` | optional | Cooling model. Absent: isothermal. |
| `[layers.<name>.radiogenics]` | optional | Radiogenics model. Absent: none. |

[Layer](../layers/layer.md#physics-switches-and-flags) describes each switch. Moduli, density, thermal, viscosity, and melting laws belong to the material.

An unknown key, a model table without `model`, a non-boolean switch, or another `state` is an error. A retired key (`class`, `type`, `material_name`, `is_tidal`, `is_solid`, `use_thermal_eos`, the gas-layer scalars, and layer-level `eos`, `shear_viscosity`, `bulk_viscosity`, and `partial_melt` tables) is refused with a message naming its replacement (`TidalPy.schema.RETIRED_LAYER_KEYS`).

### Geometry

A layer's inner radius is the outer radius of the layer below (0 for the innermost), so `radius_inner_m` is an error. Each layer gives exactly one of `radius_outer_m`, `radius_fraction`, or `volume_fraction` $f_V$. With world radius $R$, the shell volume is $f_V$ times the world volume, so $r_\mathrm{out} = \left(r_\mathrm{in}^3 + f_V R^3\right)^{1/3}$.

### The Material Table

A layer's `material` is a MatPack name or a table, built by `TidalPy.Material.load_material` (see [MatPack](../../Material/matpack.md) and [Phases and Materials](../../Material/materials.md)). A table is either:

- A preset with overrides: `preset` names a MatPack material and the rest merges over it. A law naming the same model changes its keys, one naming another model replaces it, and a new law is added.
- A full definition: a `solid` and/or `liquid` phase table, a `melting` table (`solidus`, `liquidus`, `weakening`, `bulk_modulus_mixing`, `bulk_viscosity_mixing`), and `latent_heat_j_kg`.

A phase table holds law sub-tables, each with a `model` key (`eos` for density and its thermal terms, `shear_modulus`, `shear_viscosity`, `bulk_viscosity`, and the default `shear_rheology` and `bulk_rheology`), thermal constants (`thermal_conductivity_w_mk`, `heat_capacity_j_kgk`, their temperature exponents, `thermal_reference_temperature_k`), and optionally a `preset` to take that material's phase. A liquid-only material is always liquid (an ocean, a gas envelope). A two-phase material melts between solidus and liquidus when `use_melting` is on; equal curves give one melting point.

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

A radius-tabulated law (`model = "interpolate"`, `radius_m` in \[m\]) must span its layer to within 0.1% of the outer radius at each end, or the build fails (this catches a table in km). Material errors name the layer's table. The [Schema Examples](schema_examples.md) show every law.

### Attached Physics Models

A model table holds a `model` key and that model's parameters, which go to the factory: `make_rheology` for `shear_rheology` and `bulk_rheology`, `make_cooling` for `cooling`, and `make_radiogenics` for `radiogenics`. Omitted parameters keep the model's defaults. An unknown parameter raises a `ValueError` naming the table and the closest accepted key. An isotope radiogenics table with neither a dataset (`isotopes`) nor its own isotope arrays takes `[radiogenics] isotopes` from `TidalPy_Configs.toml`. Each module's page lists its models and parameters.

### Physical Values

Every build checks the values (`validate_physical_values`), raising a `ValueError` that names the world or layer, key, value, and allowed range when:

* `radius_m` or `mass_kg` is not positive and finite, `albedo` is outside $[0, 1]$, `emissivity` outside $(0, 1]$, or the obliquity, spin, stellar temperature, or luminosity is not finite (the last two may be zero: derive them);
* a `radius_fraction` or `volume_fraction` is outside $(0, 1]$, or `radius_outer_m` is not positive;
* a layer ends at or below the layer under it, above the world's radius, or (the outermost) short of it by more than a relative $10^{-6}$;
* a layer's `mass_kg`, `tidal_scale`, or `temperature_k` is negative or not finite;
* a `layer_index` is not a non-negative integer, or two layers share one.

## Tidal Dissipation (`[tides]`)

An omitted `[tides]` key falls back to the type's `[tides.<type>]` table of `TidalPy_Configs.toml`, then its `[tides]` block (planet values, `fixed_k = 0.3` and `fixed_q = 100` at degree 2), then a built-in default. `[tides.star]` holds an $n = 3$ polytrope's fluid $k_2 = 0.0289$ (Sun-like) with $Q = 1.93 \times 10^4$, the usual stellar $Q' = 3Q / (2k_2) = 10^6$. The bundled `trappist1`, fully convective, states $n = 1.5$ values ($k_2 = 0.287$) at the same $Q'$.

| Key | Description |
|-----|-------------|
| `global_tidal_model` | `rheology`, `cpl` (`fixed_q`), `ctl` (`fixed_dt`), or `ctl_q` (`fixed_dt_q`). Default per type (`[tides.default_model]`): `rheology` for terrestrial and layered worlds, `fixed_dt` for gas giants, `fixed_q` for stars. |
| `fixed_k` | Static potential Love numbers $k_l$ for $l = 2 \ldots 10$, read only by the analytic models. A shorter list zero-fills, and a zero $k_l$ means no dissipation, so reach `max_degree_l` (a short list warns: `short_degree_list` in `[warnings]`). |
| `fixed_q` | Quality factors $Q_l$, same indexing. Read by `cpl` and `ctl_q`. |
| `fixed_dt_s` | Time lags $\Delta t_l$ \[s\], same indexing. Read by `ctl` and `ctl_q`. |
| `min_degree_l`, `max_degree_l` | Lowest and highest harmonic degree in the mode sum. Default `2` and `2`. |
| `eccentricity_trunc_lvl` | Heating kept through $e^N$ ([Eccentricity Functions](../../Tides/Eccentricity.md)). 2, 4, 6, 8, 10, 20, 50, or `"exact"` (any $e < 1$); default `10`. Other levels are promoted to the next tabulated one with a once-per-session warning. Alias `eccentricity_truncation`. |
| `eccentricity_exact_tolerance` | For `"exact"`: the dropped modes' $q^2$-weighted share of the squared eccentricity functions, which bounds the heating's relative error. In (0, 1); default `1e-4`. |
| `obliquity_trunc_lvl` | Heating kept through $I^N$: 0, 2, 4, `"off"` (0, the default), or `"gen"`/`"general"` (exact at any obliquity). Level 2 is within 1% of `"gen"` to $I \approx 8^\circ$, level 4 to $27^\circ$ ([Obliquity Functions](../../Tides/Obliquity.md)). 1 is promoted to 2, and above 4 to `"gen"`. Alias `obliquity_truncation`. |
| `layer_tidal_heating` | Whether `calc_tides` resolves each layer's heating; with the radial solver this costs about one more global solve (free otherwise). Default `true`. |
| `love_method` | `radial_solver` (`shooting`, `rs`; default), `propagation_matrix` (`prop_matrix`, `pm`, `prop`), `homogeneous` (`homogen`), `cpl`, `ctl`, or `laterally_inhomogeneous` (`3d`, `lat_inhom`; reserved, raises `NotImplementedError`). The three homogeneous-sphere methods have no depth-resolved solution, so the 3D stress, strain, and heating path raises `RuntimeError` with them. |
| `love_fixed_q` | Scalar $Q$ the `cpl` Love method applies. Unset: the tide model's `fixed_q`. |
| `love_fixed_dt_s` | Scalar lag \[s\] the `ctl` Love method applies. Unset: the tide model's `fixed_dt_s`. |

`global_tidal_model` (Love numbers to dissipation) and `love_method` (how the Love numbers are computed) are independent.

```toml
[tides]
global_tidal_model = "rheology"      # dissipation from the layers' complex moduli
min_degree_l = 2
max_degree_l = 3                     # include the degree-3 tide
eccentricity_trunc_lvl = 20          # heating through e^20; 10 is the default
obliquity_trunc_lvl = "gen"          # exact obliquity terms
love_method = "radial_solver"
```

```toml
[tides]
global_tidal_model = "fixed_q"       # a star or gas giant: per-degree parameters
fixed_k = [0.3, 0.15, 0.1]           # l = 2, 3, 4; higher degrees would be zero-filled
fixed_q = [1.0e5, 1.0e5, 1.0e5]
```

`world.get_tide_config()` returns these settings, and `get_config_dict()` saves them with the tide model's parameters. A saved `obliquity_trunc_lvl` of `"off"` reads back as `0`, and `"general"` as `"gen"`.

## Solver Settings (`[eos_solver]`, `[radial_solver]`)

A world file may pin its solver settings so it reproduces a run elsewhere, with the keys of the same-named [configuration](../../Overview/2_TidalPy_Configurations.md) sections:

- `[eos_solver]`: `integration_method`, `rtol`, `atol`, `pressure_tol`, `max_iters`, `slices_per_layer`, `nondimensionalize`, `solve_temperature`, `max_thermal_passes`, `thermal_tol`.
- `[radial_solver]`: `integration_method`, `rtol`, `atol`, `starting_method`, `degree1_frame`, `start_radius_tolerance`, `scale_rtols`, `max_num_steps`, `expected_size`, `max_ram_mb`, `nondimensionalize`.

```toml
[eos_solver]
integration_method = "RK45"
rtol = 1.0e-8

[radial_solver]
starting_method = "kamata"
```

A pinned key applies to every solve the world runs (`solve_eos`, `solve_love_numbers`, `calc_tides`, the 3D paths); a call's own argument still wins. They change how a result is computed, not what, except `degree1_frame`, the frame of degree-1 load Love numbers. `world.set_solver_defaults(eos_solver=..., radial_solver=...)` pins them on a built world and `get_solver_defaults()` returns them. Several bundled worlds pin keys ([WorldPack](worldpack.md#bundled-worlds)).

## Default Configuration Resolution

An omitted value is resolved in this order:

1. The world file or dict.
2. `TidalPy_Configs.toml`: world properties from `[worlds]` (a `[worlds.<type>]` sub-table specializes it, keeping the star-only keys off other worlds), tide values from `[tides]`, and an unnamed material from `[layers] material`.
3. The class or factory default: `Layer` defaults for switches, the model's defaults for its omitted keys, and no cooling or radiogenics model unless named.

A MatPack material's values are edited in its file (`TidalPy.Material.material_info(name)["path"]`) or overridden in the layer's material table. `TidalPy_Configs.toml` is created in the user's TidalPy `Config` directory on first use and is user-editable; its `[numerical]` floors and tolerances are listed under [Numerical Settings](../../Overview/2_TidalPy_Configurations.md#numerical-settings).

## Building a World From a Radial Profile

A world can take its interior from a PREM-like radial profile instead of layer tables, _e.g._, to compare with published data or to use another EOS solver's output. A top-level `data_file` names the delimited file:

```toml
schema_version = "0.2.0"
name = "Earth-PREM"
type = "terrestrial"
radius_m = 6371000.0
mass_kg = 5.972e24
data_file = "PREM.csv"  # Can be a path to a file too.
```

From Python, pass arrays under `data` instead (never both):

```python
import numpy as np

radius_km = np.array([0.0, 3480.0, 3480.0, 6371.0])        # Repeated boundary row
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
})                                                         # A liquid core and a solid mantle
print([layer.name for layer in world], world.layer_0.is_liquid)   # ['layer_0', 'layer_1'] True
```

`data_file` resolves as in [WorldPack](worldpack.md#locations). A profile world pins `integration_method = "RK45"` in `[eos_solver]` unless the file sets it: on an interpolated profile RK45 is about 2.8 times faster than DOP853 at equal accuracy.

### Profile Format

A profile is a delimited table (comma, semicolon, tab, or whitespace; `#` lines are comments) or a mapping of arrays. It needs a radius, a density, and both seismic velocities or both static moduli:

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

Velocities give the moduli per row: $\mu = \rho V_s^2$, $K = \rho \left(V_p^2 - \tfrac{4}{3} V_s^2\right)$. Rows may run surface-first or center-first. Reading rules:

- Columns are matched by name, ignoring case and punctuation (`Vp`, `V_P`, `vp`), in any order; other columns are skipped. A bare `eta` (PREM's anisotropy parameter) is not a viscosity and is skipped.
- A name may carry its unit (`radius_km`, `Vp [km/s]`, `rho_kg_m3`). An unconvertible unit raises a `ValueError` listing the accepted ones. In a whitespace-delimited header a bracketed unit stays with its name (`depth (km)`).
- Names come from a header row, which must name every column, or from the last `#` line before the data if it names every column. Without a header the columns are radius, density, `Vp`, `Vs`, shear viscosity, bulk viscosity, so quality factors need a header.
- An unlabeled radius or depth is km when its largest value is below 100000, else m.
- A density below 100 kg/m³, a velocity below 100 m/s, or a bulk or solid shear modulus below 10⁶ Pa is refused as an unlabeled g/cm³, km/s, or GPa column (name it `density_g_cm3`, `vp_km_s`, or `shear_modulus_gpa`).

### Layer Detection

From the center outward, every change between `Vs = 0` (liquid) and nonzero (solid) starts a new layer, and so does a radius given twice (a discontinuity). Of the two rows, the one listed first belongs to the lower layer, whichever way the file is ordered. With the jump on a layer boundary, where the integrations restart, PREM runs in half the time with Love numbers about 20 times closer to converged. Layers are named `layer_0`, `layer_1`, ..., from the center; a duplicate that would make a zero-thickness layer joins a neighbor. The bundled `PREM.csv` (PREM's 3 km ocean replaced by upper crust) gives twelve layers: the solid inner core, the liquid outer core, and ten solid layers.

Each layer's slice becomes its material: radius-tabulated (`interpolate`) laws for the density and bulk modulus (`eos`), a solid's shear modulus, and any viscosities given. A solid layer gets a solid phase and takes part in the tides. A liquid layer gets a liquid-only material (static unless `is_static = false`) and has `use_tides` off. A layer's first row is repeated at its lower boundary when the profile does not repeat it there. Nothing else is added: without viscosities the layer is elastic, with no melting or rheology.

### Refining Detected Layers

A profile holds no complex moduli, so rheology, cooling, and radiogenics go in `[layers.<name>]` tables. A table refines one detected layer (by `layer_index` or the name `layer_N`), or every layer between two detected boundaries with `radius_range_m = [inner, outer]` \[m\]. It names them: its own name for one, `mantle_0`, `mantle_1`, ... from the inside for several. Unrefined layers need no table.

- An outer radius given (`radius_outer_m` or `radius_fraction`) must match the detected boundary, and both ends of `radius_range_m` must be detected boundaries. A range excludes `layer_index` and an outer radius.
- Two tables may not claim one layer, and an index past the detected layers is an error.
- A `material` table merges over the profile's: another model replaces a law (a constant viscosity, say), and a new law is added. A MatPack name is refused.
- A single number for a tabulated quantity (_e.g._, `bulk_modulus_pa = 1.0e11` in `solid.eos`) holds across the layer: TOML overrides the data file.
- Every other key (switch, flag, model table) is set on the layer.

```toml
[layers.mantle]             # refines every detected layer from the core-mantle boundary to the surface
radius_range_m = [3480.0e3, 6371.0e3]
[layers.mantle.shear_rheology]
model = "maxwell"
[layers.mantle.material.solid.shear_viscosity]
model = "constant"          # the profile named no viscosity, so give one here
reference_viscosity_pas = 1.0e21
```

### Quality Factors in Place of Viscosities

Seismic profiles such as PREM give $Q_\mu$ and $Q_\kappa$ (at 1 s), not viscosities. `q_provided = true` builds each solid layer's loss from them with the [`seismic_q` rheology](../../Rheology/rheology_models.md): the profile's static modulus with $Q(\omega) = Q_\mathrm{ref}(\omega/\omega_\mathrm{ref})^{a}$ and the dispersion that loss implies.

```toml
data_file = "PREM.csv"                             # has q_mu and q_kappa columns
q_provided = true
q_reference_frequency_rad_s = 6.283185307179586    # 1 s period
q_frequency_exponent = 0.0                         # a in [0, 1)
```

| Key | Default | Meaning |
|-----|---------|---------|
| `q_provided` | `false` | Use the profile's quality factors (off, its Q columns are read but unused). The profile must give `q_mu` and no viscosities. |
| `q_reference_frequency_rad_s` | $2\pi$ | Frequency \[rad s⁻¹\] of the profile's $Q$ and moduli. |
| `q_frequency_exponent` | 0.0 | $a$ in $Q \propto \omega^{a}$. $a = 0$ keeps the seismic $Q$ at tidal periods; $a > 0$ lowers it there. |

These keys need a profile, and the last two need `q_provided = true`, so none is silently ignored. Each solid layer gets a `seismic_q` shear rheology, plus a `seismic_q` bulk rheology when the profile gives `q_kappa`; its quality factors (held in the layer's viscosity arrays) must be positive. Liquid layers get no rheology. A refinement table may name `seismic_q` with its own `reference_frequency_rad_s` or `q_frequency_exponent`, or an `elastic` bulk rheology to ignore $Q_\kappa$. It may not name another rheology, a viscosity law, or a `preset`, or, on a solid layer, turn on `use_melting` or name the `convection` cooling model, since each would use a viscosity where the layer holds $Q_\mu$.

The bundled `earth_prem_q` ([`example_profile_q_world.toml`](examples/example_profile_q_world.toml) shows the keys) has, at 1 s, the elastic PREM's $k_2$ (0.298363 against 0.298365). At the M2 tide dispersion raises it to 0.302 with $|k_2|/|\mathrm{Im}\,k_2| \approx 500$, or about 100 with `q_frequency_exponent = 0.15`.

### Saving a Profile World

`portable_config` keeps the configuration as given (the `data_file` reference and the refining tables). While the world is unchanged since its build, `save_to_toml` writes that, with the reference rewritten for the destination folder, so the file holds no expanded profile. A changed world is saved as its live `get_config_dict()`, profile included. `source_config` holds the expanded form as built.

## System Schema

A system file has a top-level `schema_version` and `name`, then one `[worlds.<key>]` table per member, whose key is the member's name.

| Key | Required | Description |
|-----|----------|-------------|
| `world` | **yes** | A bundled world name, a path to a world file (relative paths try the system file's folder, then the working directory), or an inline `[worlds.<key>.world]` table. |
| `tidal_host` | optional | The key of the world that raises this one's tides, possibly declared later. Absent: not tidally forced. Two worlds may host each other; they share one orbit, which only one needs to state. |
| `is_star` | optional | Marks the insolation source (at most one; need not be a tidal host). |
| `semi_major_axis_m` | optional | Semi-major axis about the tidal host \[m\]. Requires `tidal_host`. |
| `eccentricity` | optional | Eccentricity about the tidal host. Requires `tidal_host`. |
| `synchronous` | optional | `true` sets the world's spin to its mean motion about the tidal host once the system is built (`System.set_synchronous_rotation`), replacing the world file's spin. Requires `tidal_host` and a semi-major axis (its own or its mutual partner's). |
| `stellar_semi_major_axis_m` | optional | Distance from the star \[m\], for a world whose host is not the star. |
| `stellar_eccentricity` | optional | Eccentricity about the star. |

A system has at least one world and at most one star, and no system-wide host.

```toml
schema_version = "0.2.0"
name = "Sol System"

[worlds.sun]
world = "sol"          # a bundled name, a path, or an inline table
is_star = true

[worlds.earth]
world = "earth_simple"
tidal_host = "sun"
semi_major_axis_m = 1.495978707e11
eccentricity = 0.0167  # the host is the star, so no stellar orbit is needed
```

[System](../system/system.md#building-a-system-from-toml-build_system) covers `build_system` and the system API.

## Round Trip

`get_config_dict()` is the world as it stands now, every change since the build included. `build_world` and `build_system` rebuild the same class, parameters, and models from it, and `build_layer_from_dict` a standalone layer. A world constructed in Python with no tide model writes no `tides` table, so its rebuild takes its type's default tide model.

```python
from TidalPy.Structures import build_layer_from_dict

world.set_spin_frequency(2.0e-5)                 # A change after the build
config = world.get_config_dict()
twin = build_world(config)                       # Same class, parameters, and models
assert twin.get_config_dict() == config

mantle_twin = build_layer_from_dict(world.mantle.get_config_dict())   # A standalone layer
system_twin = build_system(system.get_config_dict())                  # Worlds, hosts, star, orbits
```

The dict states every layer value and model table explicitly, so a rebuild does not depend on `[layers] material` or a MatPack file that may have changed. A standalone layer's dict adds `name` and `radius_inner_m`; a system's inlines each member's dict under `world`. Solved state (EOS, Love numbers, tides) is not saved; solve again after a rebuild.

## Schema Version

`schema_version` is graded against the current one:

- Patch difference (`0.0.X`): allowed silently.
- Minor difference (`0.X.0`): allowed with a warning that some functionality may break.
- Major difference (`X.0.0`): refused with a `ValueError`.
- Missing in a file: allowed with a warning. A Python dict without one builds silently as the current schema.

`force=True` (to `build_world`, `BaseWorld.build`, or `validate_schema_version`) skips these checks.

## Python API

All are in `TidalPy.Structures.configs`. The main ones (`build_world`, `build_layer_from_dict`, `build_system`, `load_world`, `load_system`, `available_worlds`, `available_systems`, `save_world_to_toml`, `install_worldpack`, `SCHEMA_VERSION`) are also in `TidalPy.Structures`, with the classes and `make_tide`.

### High Level

* `build_world(source, overrides=None, force=False) -> BaseWorld` (wraps `BaseWorld.build`) and `build_system(source, overrides=None, force=False)`: build from a bundled name, file path, or dict. `overrides` merges over the source (`TidalPy.configurations.merge_configs`: tables key by key, other values replace). `force=True` skips the schema-version check. A binary file is loaded with `load_world` or `load_system` and refuses `overrides`.
* `load_world(path, force=False) -> BaseWorld`, `load_system(path, force=False) -> System`: read a binary file into a new object of its saved class. A file of the other kind or not a TidalPy file raises `IOError` saying so (`"it is a System file; load it with load_system"`). `TidalPy.Utilities.binary.binary_file_class(path)` returns a file's class name, or `None`.
* `world.save_to_toml(path, overwrite=True)`: write `world.get_save_config(destination_dir)` with the current `schema_version`, under a header naming the TidalPy, SciPy, and CyRK versions. That is the live `get_config_dict()` (`type`, `layers`, `tides`, `schema_version`), validated first (else `ValueError`), except for an unchanged profile world ([Saving a Profile World](#saving-a-profile-world)). EOS-solved layer masses and the name do not count as changes. `path` is a `str` or `pathlib.Path`.
* `build_layer_from_dict(config) -> Layer`: a standalone layer from a layer's `get_config_dict()` ([Round Trip](#round-trip)), leaving the dict unmodified; a non-`dict` raises `TypeError`.
* `world.config` (alias `source_config`, `None` if constructed directly), `world.portable_config` (profile worlds only), `world.built_config` (`get_config_dict()` after the build, the reference for "unchanged"): the configurations kept from the build. A successful `load_binary` clears all three.
* `load_radial_data(source, surface_radius=None) -> dict`: read a profile (path or mapping) into MKS arrays ascending in radius. `detect_layer_boundaries(radius, shear_modulus)` returns its `(start, end, is_solid)` runs.
* `EOS_SOLVER_KEYS`, `RADIAL_SOLVER_KEYS`, `validate_solver_table(section, table, where)`: the solver-table keys and their check.
* `available_worlds() -> list[str]`, `available_systems() -> list[str]`: the bundled world and system names. `build_world` on a system config raises a `ValueError` naming `build_system`, and the reverse.
* `install_worldpack(force=False) -> str`: copy the packaged worlds into the user data directory (copy-if-absent unless `force`); returns it.

An unknown bundled name raises `FileNotFoundError` naming the closest name and listing all. An unknown world, `[tides]`, layer, solver-table, or system key raises `ValueError` naming the closest accepted key (`"did you mean 'temperature_k'?"`).

### Low Level

* `load_build_config(source, overrides=None, force=False) -> dict`: the configuration a build uses (loaded, version-graded, overrides merged, structural defaults filled).
* `save_world_to_toml(config, path, overwrite=True)`: serialize a config dict.

### Loader / Validation

The key sets (`WORLD_TYPES`, `LAYER_SCALAR_KEYS`, `LAYER_MODEL_SECTIONS`, `LAYER_STATES`, `RETIRED_LAYER_KEYS`, `ALLOWED_TIDES_KEYS`, `EOS_SOLVER_KEYS`, ...) live in `TidalPy.schema` and are re-exported by `TidalPy.Structures.configs.toml_loader`, which also has:

* `SCHEMA_VERSION`: `"0.2.0"`.
* `load_toml(source)`: parse a path or pass a dict through.
* `validate_schema_version(config, force=False)`: see [Schema Version](#schema-version).
* `validate_world_config(config)`, `validate_layer_config(name, cfg)`: structural checks; the first ends with `validate_physical_values(config)` ([Physical Values](#physical-values)).
* `merge_with_defaults(config)`: apply structural (non-physical) defaults.
