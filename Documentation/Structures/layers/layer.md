# Layer (`Structures.layers`)

_Updated: 2026-10-07_

`TidalPy.Structures.layers.Layer` is one spherically symmetric shell of a world: its radii \[m\] and mass \[kg\], its material, the physics switches that say how much of that material it uses, its temperature \[K\], the assumptions the radial (Love number) solver makes inside it, and optional rheology overrides and cooling and radiogenics models. Rock, ice, iron, an ocean, or a gas envelope differ only in material and switches. A world's EOS solve evaluates the layer's material at every radius it integrates, and the layer's profile getters read those values back. `Layer` is a `StructureBase` (C++ `c_Layer`).

```python
from TidalPy.Structures.layers import Layer

mantle = Layer(
    "mantle",                       # Name
    1,                              # Layer index (0 = innermost)
    3.485e6,                        # Inner radius [m]
    6.371e6,                        # Outer radius [m]
    material="peridotite",          # A MatPack material
    temperature=1600.0,             # [K]
    use_melting=True,
    use_pressure_melting=True,
    cooling="convection",
)                                   # A melting, convecting upper mantle
print(mantle.thickness, mantle.state, mantle.can_change_state)   # 2886000.0 auto True
print(mantle)                       # Layer('mantle', index=1, radius_inner_km=3485, radius_outer_km=6371, ...)
```

Inside a world the index and inner radius follow from the stack, so a layer for `add_layer` needs only `radius_outer` (see [Worlds: Layers](../worlds/worlds.md#layers)).

## Construction

```text
Layer(name, layer_index=None, radius_inner=None, radius_outer=None, mass=0.0, material=None, *,
      use_tides=True, is_volume_fixed=True, tidal_scale=None, state="auto", is_static=True,
      is_incompressible=False, temperature=0.0, use_thermal_expansion=False, use_melting=False,
      use_pressure_melting=False, use_melt_density=False, use_heating=False,
      shear_rheology=None, bulk_rheology=None, cooling=None, radiogenics=None)
```

Everything after `material` is keyword-only, and `material` and every keyword-only argument are also read-write properties. The switches and flags are in [Physics Switches and Flags](#physics-switches-and-flags); the rheologies and models (each a model object, a model name, or a config table) in [Rheology Overrides](#rheology-overrides) and [Cooling and Radiogenics](#cooling-and-radiogenics).

| Parameter | Units | Default | Description |
|---|---|---|---|
| `name` | - | | Unique within a world, which reaches the layer by it (`world.mantle`). |
| `layer_index` | - | `None` | Position in the world, 0 for the innermost. Left out, `add_layer` gives the layer its place in the stack; given, `add_layer` refuses any other place. A standalone layer without one reports 0. |
| `radius_inner` | m | `None` | Left out, `add_layer` starts the layer at the top of the stack; given, `add_layer` refuses a mismatch. A standalone layer without one starts at 0. |
| `radius_outer` | m | | Required (by keyword when the two before it are left out), with `0 <= radius_inner <= radius_outer`; otherwise `ValueError`. |
| `mass` | kg | `0.0` | Replaced by each successful world EOS solve with the mass the solved profile places in the layer. A layer that holds its mass (`is_volume_fixed = False`) holds this one when positive. |
| `material` | - | `None` | A `Material`, a MatPack name, or a material table (see [Material](#material)); a world's EOS solve needs one in every layer. |

The geometry is read-only: `name`, `layer_index`, `radius` (the outer radius), `radius_inner`, `radius_outer`, `thickness`, `volume` \[m$^3$\], `surface_area_inner` and `surface_area_outer` \[m$^2$\], `mass`, and `density_bulk` (`mass / volume` \[kg m$^{-3}$\], NaN at zero volume). `set_radii(radius_inner, radius_outer)` moves both boundaries; on a layer of a world it leaves the world's radius and other layers alone (keep the stack continuous), and the world forgets its solve.

## Material

A `Material` ([Phases and Materials](../../Material/materials.md)) has a solid phase, a liquid phase, or both, with melting curves, melt weakening, bulk-mixing laws, and a latent heat. The `material` argument and property take:

- A `Material` object.
- A MatPack name, loaded with `load_material` (see [MatPack](../../Material/matpack.md); `TidalPy.Material.available_materials()` lists them).
- A material table: a `preset` naming a MatPack material plus overrides, or a full `solid`, `liquid`, and `melting` definition, as under `[layers.<name>.material]` in a world file.

```python
from TidalPy.Material import Material, Phase, load_material

mantle.material = load_material("peridotite")   # A Material object

mantle.material = {    # A MatPack preset with one law replaced
    "preset": "peridotite",
    "solid": {"shear_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e20}},
}

mantle.material = Material(    # A full definition: one solid phase, no melting
    solid=Phase(
        eos={"model": "constant", "reference_density_kg_m3": 3300.0},
        shear_modulus={"model": "constant", "shear_modulus_pa": 6.0e10},
        shear_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e21},
    ),
)
print(mantle.material_set)   # True
```

A material is immutable and shared: one can serve any number of layers and solves, and `mantle.material` returns the layer's material, not an editable copy. To change a value, set a new material (`load_material` takes overrides, and `Material.replace` swaps a slot). Setting the material of a layer in a solved world makes the world forget its solve.

`calc_state(pressure, temperature=None, radius=None)` evaluates the material at a point with the layer's four material switches and returns the dict of `Material.calc_state`: `phase`, `density`, `bulk_modulus`, `adiabatic_bulk_modulus`, `thermal_expansion`, `heat_capacity`, `thermal_conductivity`, `shear_modulus`, `shear_viscosity`, `bulk_viscosity`, `melt_fraction`, `solidus`, `liquidus`, and `latent_expansion`. The temperature defaults to the layer's own, and the radius (read only by tabulated laws) to the mid-radius. Any argument may be an array. A layer without a material raises `ValueError`.

```python
mantle.material = "peridotite"
state = mantle.calc_state(1.0e9)                          # At 1 GPa and the layer's 1600 K
print(state["phase"], state["melt_fraction"])             # solid 0.0
state = mantle.calc_state(1.0e9, temperature=1900.0)
print(state["phase"], round(state["melt_fraction"], 3))   # partial 0.447
```

## Physics Switches and Flags

The four material switches and `use_heating` default to off, the simplest and fastest case. The material switches apply wherever the material is evaluated: the world's EOS solve, the thermal network, the radial solver, and `calc_state`.

| Switch | Default | Off | On |
|---|---|---|---|
| `use_thermal_expansion` | `False` | The density ignores the temperature. | The density and the thermal pressure follow the temperature. The expansivity shapes an adiabat either way. |
| `use_melting` | `False` | The material's base phase only (its solid, or the liquid of a liquid-only material), so the state is fixed. | The liquid phase, melt fraction, and melt weakening take part, and the layer can split into solid and liquid zones. |
| `use_pressure_melting` | `False` | The solidus and liquidus are read at zero pressure. | They follow the local pressure, and inside a melting range the latent heat also steepens an adiabat (`latent_expansion`). |
| `use_melt_density` | `False` | The density is the solid phase's. | The density mixes the two phases' by melt fraction. |
| `use_heating` | `False` | The layer generates no heat, whatever models it carries. | The world's heat sources (radiogenic, tidal, prescribed) act inside it during a thermal EOS solve. |
| `use_tides` | `True` | No part in tidal dissipation: its tidal scale is 0 in the quasi-homogeneous Love methods, and the radial solver deforms it with its static moduli but applies no rheology, so its moduli (and `calc_complex_shear_modulus`) are purely real and it adds nothing to Im(k) or to the heating. | The layer takes part in the tides. |

| Flag | Default | Description |
|---|---|---|
| `state` | `"auto"` | How the radial solver treats the layer: `"auto"`, `"solid"`, or `"liquid"` (case-insensitive; anything else raises `ValueError`). See [State and Zones](#state-and-zones). |
| `is_static` | `True` | The static equations, without inertia (Beuthe 2015). `False` takes the dynamic form, which a liquid needs at short forcing periods. The liquid zones of a layer that can change state take the layer's choice. |
| `is_incompressible` | `False` | The incompressible equations, which the propagation-matrix Love method needs. |
| `temperature` | `0.0` | \[K\]: the profile of a solve without a temperature contrast, and the layer's lumped temperature in a thermal solve. 0 K is the cold, rigid limit of the viscosity laws, and a layer whose temperature is not positive takes no part in the thermal network; a solve warns once when that leaves the layer rigid. |
| `is_volume_fixed` | `True` | `False` makes the layer hold its mass instead of its volume: the EOS solve ends it where it encloses that mass, and the layers above move with it (see [Layer Size](../worlds/worlds.md#layer-size)). |
| `tidal_scale` | `None` | The layer's share of the planet in the quasi-homogeneous Love methods (`homogeneous`, `cpl`, `ctl`) and of an analytic tide model's heating. `None` takes its volume over the planet's. |

`get_tidal_heating()` returns the heating \[W\] the world's last `calc_tides` put in the layer, NaN before one (see [Tidal Heating of Each Layer](../worlds/worlds.md#tidal-heating-of-each-layer)).

### State and Zones

With `state = "auto"` the material decides: a liquid-only material (`water`, `simple_liquid_iron`, `h2_he_molecular`) makes a liquid layer (`is_liquid`), any other a solid one. A solid layer whose material melts and that has `use_melting` on can change state (`can_change_state`): the world's EOS solve splits it into solid and liquid zones where its post-melt rigidity $\mu / (\bar{\rho} g R)$ crosses `[numerical] minimum_solid_rigidity`, and the Love solve takes each zone as a layer (see [Pieces and Zones](../worlds/worlds.md#pieces-and-zones)). `"solid"` or `"liquid"` forces the state for an idealized problem (`"liquid"` sets `is_liquid`), and the layer is never split.

### Changes That Clear a Solve

When something the EOS solve reads changes, the layer's world forgets its solved structure (see [Solved State](../worlds/worlds.md#solved-state)):

- The material, the temperature, any of the four material switches, `use_heating`, or `is_volume_fixed`
- The cooling or radiogenics model
- The radii, moved with `set_radii`
- A `state` change that changes `can_change_state` (the EOS solve looks for zones only in a layer that can change state)

The rheologies, `use_tides`, `tidal_scale`, `is_static`, `is_incompressible`, and other `state` changes are read afresh by each Love or tidal solve and leave the solved structure standing.

## Cooling and Radiogenics

A [cooling model](../../Cooling/cooling_models.md) sets how heat moves through the layer in a thermal EOS solve, and so its temperature profile (see [Temperature and Heat Flow](../worlds/worlds.md#temperature-and-heat-flow)). A [radiogenics model](../../Radiogenics/radiogenics_models.md) heats the layer during a solve when `use_heating` is on.

The `cooling` and `radiogenics` arguments and properties take a model, a name (`"off"`, `"conduction"`, `"convection"`; `"off"`, `"isotope"`, `"fixed"`), or a config table with a `model` key; `None` clears it (a layer without a cooling model is held at one temperature). The layer shares a model rather than consuming it, so one model can serve several layers. `cooling_set` and `radiogenics_set` report whether one is attached, and `calc_radiogenic_heating(time, mass)` returns the heating \[W\] at a time \[s\] for a mass \[kg\], 0.0 without a model.

```python
mantle.cooling = {
    "model": "convection",
    "critical_rayleigh": 1100.0,
}  # Boundary-layer convection
mantle.radiogenics = {
    "model": "isotope",
    "isotopes": "bulk_silicate_earth",
}  # A decaying isotope dataset
mantle.use_heating = True  # Let the radiogenics heat the layer in a thermal solve
print(mantle.calc_radiogenic_heating(1.4516496e17, 4.0e24))   # [W] at the dataset's reference time, about 2e13
```

## Rheology Overrides

A rheology maps the static modulus and viscosity onto a complex modulus at a forcing frequency (see [Rheology Models](../../Rheology/rheology_models.md)). A material's phase may carry default shear and bulk rheologies (MatPack rocks carry Andrade), and the layer may override either. `shear_rheology` and `bulk_rheology` return the rheology in effect: the override, else the material's default, else `None`. Setting one sets the override (a model, name, or config table, shared rather than consumed); `None` clears it. With no rheology the complex modulus is the real static modulus: perfectly elastic.

```python
from TidalPy.Rheology import Maxwell

print(type(mantle.shear_rheology).__name__)   # Andrade, from the material
mantle.shear_rheology = Maxwell()             # Override it
mantle.bulk_rheology = {
    "model": "voigt",
    "voigt_modulus_frac": 0.5,
}                                             # A table with a model key
mantle.shear_rheology = None                  # Back to the material's Andrade
```

## Profile Getters

A world's EOS solve populates every layer's profile, and the getters read it back at any radius \[m\] (float or `np.ndarray` in, the same shape out, NaN where unsolved). Each value comes from the solve's dense output at that radius, so a getter and the solve agree exactly.

| Getter | Returns |
|---|---|
| `get_density`, `get_gravity`, `get_pressure` | Density \[kg m$^{-3}$\], gravity \[m s$^{-2}$\], and pressure \[Pa\]. |
| `get_shear_modulus`, `get_bulk_modulus` | Post-melt static shear modulus and adiabatic bulk modulus \[Pa\]. |
| `get_shear_viscosity`, `get_bulk_viscosity` | Post-melt viscosities \[Pa s\]. |
| `get_melt_fraction` | Melt fraction: 0 for a layer that does not melt, 1 for a liquid-only material. |
| `get_static_viscoelastics(radius)` | `(shear_modulus, shear_viscosity, bulk_modulus, bulk_viscosity)` from one evaluation per radius. |
| `get_state(radius)` | Every profile above as a dict. |
| `calc_complex_shear_modulus(radius, frequency)`, `calc_complex_bulk_modulus(radius, frequency)` | The rheology in effect applied to the solved post-melt modulus and viscosity at `radius`, at a frequency \[rad s$^{-1}$\]. |

The one-argument forms `calc_complex_shear_modulus(frequency)` and `calc_complex_bulk_modulus(frequency)` apply the rheology to the material at zero pressure and the layer's temperature, with no solve needed. `eos_data_populated` and `viscoelastic_populated` report whether a profile and its material state are present. `update_eos_data(radius, density, gravity, pressure)` sets the density, gravity, and pressure profile from arrays (radius ascending), for tests and hand-built profiles; such a profile has no material state, so the other getters read NaN.

```python
import numpy as np

from TidalPy.Structures.worlds import TerrestrialWorld

world = TerrestrialWorld("Earth", 6.371e6, 5.972e24)
world.add_layer(Layer(
    "core", 0, 0.0, 3.485e6,
    material="simple_iron_core",
    temperature=4000.0))                       # A uniform solid core
world.add_layer(Layer(
    "mantle", 1, 3.485e6, 6.371e6,
    material="simple_rock",
    temperature=1600.0))                       # A uniform rock mantle
world.solve_eos()                              # Populates both layers' profiles

radii = np.linspace(3.6e6, 6.3e6, 4)
print(world.mantle.get_pressure(radii))        # [Pa] from the solve's dense output
print(world.mantle.calc_complex_shear_modulus(
    5.0e6,
    2.0e-5))                                  # Andrade at 5000 km and 2e-5 rad/s [Pa]
```

A layer of a world takes turns with the world's solves on other threads, and one call on an array reads every value from one solve (see [Solved State](../worlds/worlds.md#solved-state)).

## Serialization

`get_config_dict()` returns the layer as the builder's layer table: the constructor's values under their config names (`radius_inner_m`, `radius_outer_m`, `temperature_k`, ...), `mass_kg` once the layer has a mass (given or solved; an unset 0.0 is left out), `tidal_scale` when set, the `material` table, and the rheology-override, `cooling`, and `radiogenics` tables when set. `name` and `radius_inner_m` belong to a standalone layer only (`LAYER_STANDALONE_CONFIG_KEYS`); a world drops them when it nests the layer under its name. `save_config(path)` writes the dict to TOML, and `TidalPy.Structures.build_layer_from_dict` builds a standalone layer from it.

`save_binary(path)` and `load_binary(path)` (a `str` or `os.PathLike` path) write and read the layer with its material and models (binary class id 100; see [Binary Serialization](../../Utilities/binary.md)), but not its solved profile or tidal heating. A layer that belongs to a world cannot be loaded in place: load the world, or a standalone layer.

```python
from TidalPy.Structures import build_layer_from_dict

config = world.mantle.get_config_dict()  # Builder-valid layer table, material included
twin = build_layer_from_dict(config)     # A standalone layer, owned by no world
print(twin.get_config_dict() == config)  # True
```

## TOML Table

A world file builds each layer from a `[layers.<name>]` table ([TOML Schema](../config/toml_schema.md#layer-level-schema) has every key, [Schema Examples](../config/schema_examples.md) every form). The key is the layer's name, the inner radius comes from the layer below, and exactly one of `radius_outer_m`, `radius_fraction`, or `volume_fraction` sets the outer radius. Scalars take their config spelling (`temperature_k`, `mass_kg`), and models are sub-tables:

```toml
[layers.mantle]
layer_index = 1
radius_fraction = 1.0
temperature_k = 1600.0
use_melting = true
use_pressure_melting = true
use_heating = true
material = "peridotite"             # a MatPack name; or a [layers.mantle.material] table

[layers.mantle.shear_rheology]
model = "maxwell"                   # overrides the material's Andrade

[layers.mantle.cooling]
model = "convection"

[layers.mantle.radiogenics]
model = "isotope"
isotopes = "modern_day_chondritic"
```

A layer that names no material takes `[layers] material` of `TidalPy_Configs.toml` (`simple_rock`), and a layer gets no cooling or radiogenics model unless its table names one.

## C++ API

`c_Layer : c_StructureBase` is built from a `c_LayerConfig` holding the constructor's scalars (`state` as a `c_LayerState`, the switches as `c_MaterialSwitches`). The material is a `shared_ptr<const c_Material>` (`set_material`, `get_material`), and the rheology overrides and models are shared pointers too (`set_shear_rheology`, `set_bulk_rheology`, `set_cooling`, `set_radiogenics`). `get_shear_rheology()` returns the rheology in effect and `apply_shear_rheology(static_modulus, viscosity, frequency)` makes a complex modulus. `calc_state(point, out)` fills a `c_MaterialState` with the layer's switches. `get_is_liquid()`, `get_can_change_state()`, the vectorized reads `get_eos_fields(field_indices, num_fields, radii, num_radii, values_out)` and `calc_complex_moduli`, and `c_layer_state_from_name` / `c_layer_state_name` complete it.

## References

- Beuthe, M. (2015). Tidal Love numbers of membrane worlds: Europa, Titan, and Co. *Icarus*, 258, 239-266.
