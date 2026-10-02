# Layer (`Structures.layers`)

_Updated: 2026-10-01_

`TidalPy.Structures.layers.Layer` contains one spherically symmetric shell of a world: its radii \[m\] and mass \[kg\], the material it is made of, the physics switches that say how much of that material it uses, its temperature \[K\], the assumptions the radial (Love number) solver makes inside it, optional rheology overrides, and optional cooling and radiogenics models. A world's EOS solve evaluates the layer's material at every radius it integrates, and the layer's profile getters read those values back. There is one layer class for every kind of shell: rock, ice, iron, an ocean, or a gas envelope differ only in their material and switches.

## Inheritance

```
TidalPyBaseClass
  └── StructureBase
        └── Layer
```

In C++ the class is `c_Layer` (`Structures/layers/layer_.hpp`).

## Construction

```python
Layer(
    name:                  str,
    layer_index:           int,
    radius_inner:          float,
    radius_outer:          float,
    mass:                  float = 0.0,
    material:              Material | str | dict = None,
    *,
    use_tides:             bool  = True,
    is_volume_fixed:       bool  = True,
    tidal_scale:           float = None,
    state:                 str   = "auto",
    is_static:             bool  = True,
    is_incompressible:     bool  = False,
    temperature:           float = 0.0,
    use_thermal_expansion: bool  = False,
    use_melting:           bool  = False,
    use_pressure_melting:  bool  = False,
    use_melt_density:      bool  = False,
    use_heating:           bool  = False,
    shear_rheology:        RheologyBase | str | dict = None,
    bulk_rheology:         RheologyBase | str | dict = None,
    cooling:               CoolingBase | str | dict = None,
    radiogenics:           RadiogenicsBase = None,
)
```

Everything after `material` is keyword-only. `material` and every keyword-only argument are also properties of the same name, read and written after construction, except `radiogenics`, which is attached with `set_radiogenics`. The geometry is read-only (see below).

**Geometry and material**

| Parameter | Units | Default | Description |
|---|---|---|---|
| `name` | - | | The layer's name (`"mantle"`). A world reaches the layer by it (`world.mantle`), so names are unique within a world. |
| `layer_index` | - | | Position in the world, 0 for the innermost layer. |
| `radius_inner`, `radius_outer` | m | | Boundary radii, with `0 <= radius_inner <= radius_outer`; otherwise `ValueError`. |
| `mass` | kg | `0.0` | Each successful world EOS solve replaces it with the mass the solved profile places in the layer. A layer that holds its mass (`is_volume_fixed = False`) holds this one when it is positive. |
| `material` | - | `None` | A `Material`, a MatPack name (`"peridotite"`), or a material table (see [Material](#material)). A world's EOS solve needs every layer to have one. |

The switches and flags are described in [Physics Switches and Flags](#physics-switches-and-flags), and the models in [Cooling and Radiogenics](#cooling-and-radiogenics) and [Rheology Overrides](#rheology-overrides).

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
```

The geometry is read through properties: `name`, `layer_index`, `radius` (the outer radius), `radius_inner`, `radius_outer`, `thickness`, `volume` \[m$^3$\], `surface_area_inner` and `surface_area_outer` \[m$^2$\], `mass`, and `density_bulk` (`mass / volume` \[kg m$^{-3}$\], NaN for a zero-volume layer). `set_radii(radius_inner, radius_outer)` moves both boundaries and keeps the derived geometry in step. On a layer of a world it leaves the world's radius and its other layers alone, so keep the stack continuous, and the world forgets its solved structure.

## Material

The material is a `Material` from `TidalPy.Material`: a solid phase, a liquid phase, or both, with melting curves, melt weakening, bulk-mixing laws, and a latent heat (see [Phases and Materials](../../Material/materials.md)). The `material` argument and property take it in three forms:

- A `Material` object.
- A MatPack name, loaded with `load_material` (see [MatPack](../../Material/matpack.md); `TidalPy.Material.available_materials()` lists them).
- A material table: a `preset` naming a MatPack material plus overrides, or a full `solid`, `liquid`, and `melting` definition. These are the tables a world TOML file gives under `[layers.<name>.material]`.

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

A material is immutable and shared: one material can serve any number of layers and solves, and `mantle.material` returns the layer's material, not a copy that could be edited. To change a value, set a new material (`load_material` takes overrides, and `Material.replace` swaps a slot). Setting the material of a layer in a solved world makes the world forget its solved structure.

`calc_state(pressure, temperature=None, radius=None)` evaluates the layer's material at a point with the layer's four material switches, and returns the dict `Material.calc_state` returns: `phase`, `density`, `bulk_modulus`, `adiabatic_bulk_modulus`, `thermal_expansion`, `heat_capacity`, `thermal_conductivity`, `shear_modulus`, `shear_viscosity`, `bulk_viscosity`, `melt_fraction`, `solidus`, `liquidus`, and `latent_expansion`. The temperature defaults to the layer's own, and the radius (read only by tabulated laws) to the layer's mid-radius. Any argument may be an array. It raises `ValueError` for a layer without a material.

```python
mantle.material = "peridotite"
state = mantle.calc_state(1.0e9)                          # At 1 GPa and the layer's 1600 K
print(state["phase"], state["melt_fraction"])             # solid 0.0
state = mantle.calc_state(1.0e9, temperature=1900.0)
print(state["phase"], round(state["melt_fraction"], 3))   # partial 0.447
```

## Physics Switches and Flags

The four material switches and `use_heating` default to off, the simplest and fastest case. The material switches decide how much of the material the layer uses wherever it is evaluated: the world's EOS solve, the thermal network, the radial solver, and `calc_state`.

| Switch | Default | Off | On |
|---|---|---|---|
| `use_thermal_expansion` | `False` | The density ignores the temperature. | The density and the thermal pressure follow the temperature. The expansivity shapes an adiabat either way. |
| `use_melting` | `False` | The material's base phase only (its solid phase, or the liquid of a liquid-only material), so the state is fixed. | The liquid phase, the melt fraction, and the melt weakening take part, and the layer can split into solid and liquid zones. |
| `use_pressure_melting` | `False` | The solidus and liquidus are read at zero pressure. | The solidus and liquidus follow the local pressure. Inside a melting range the latent heat then also steepens an adiabat (`latent_expansion`). |
| `use_melt_density` | `False` | The density is the solid phase's. | The density mixes the two phases' by melt fraction. |
| `use_heating` | `False` | The layer generates no heat, whatever models it carries. | The world's heat sources (radiogenic, tidal, and prescribed) act inside the layer during a thermal EOS solve. |
| `use_tides` | `True` | The layer takes no share of the tidal heating: its tidal scale is 0, and `calc_tides` gives it no heating, leaving the total to the tidal layers. | The layer takes part in the tides. |

The flags set the radial solver's assumptions, the temperature, and how the layer sizes itself:

| Flag | Default | Description |
|---|---|---|
| `state` | `"auto"` | How the radial solver treats the layer: `"auto"`, `"solid"`, or `"liquid"` (case-insensitive; anything else raises `ValueError`). See [State and Zones](#state-and-zones). |
| `is_static` | `True` | The static equations, without inertia (see Beuthe 2015). `False` takes the dynamic form, which a liquid needs at short forcing periods. The liquid zones of a layer that can change state take the layer's choice. |
| `is_incompressible` | `False` | The incompressible equations. The propagation-matrix Love method needs it. |
| `temperature` | `0.0` | The layer temperature \[K\]: the profile of a solve without a temperature contrast, and the layer's lumped temperature in a thermal solve. The default 0 K is the cold, rigid limit of the viscosity laws, and a layer whose temperature is not positive takes no part in the thermal network. A solve warns once when that leaves the layer rigid. |
| `is_volume_fixed` | `True` | `False` makes the layer hold its mass instead of its volume: the EOS solve ends it where it encloses that mass, and the layers above it move with it (see [Worlds](../worlds/worlds.md#layer-size)). |
| `tidal_scale` | `None` | The layer's share of the planet in the quasi-homogeneous Love methods (`homogeneous`, `cpl`, `ctl`) and of an analytic tide model's heating. `None` takes the layer's volume over the planet's. |

`get_tidal_heating()` returns the heating \[W\] the world's last `calc_tides` put in the layer, NaN before one (see [Worlds](../worlds/worlds.md#tidal-heating-of-each-layer)).

### State and Zones

With `state = "auto"` the material decides: a liquid-only material (`water`, `simple_liquid_iron`, `h2_he_molecular`) makes a liquid layer, and any other material a solid one. A solid layer whose material melts and that has `use_melting` on can change state: the world's EOS solve splits it into solid and liquid zones where its post-melt rigidity $\mu / (\bar{\rho} g R)$ crosses `[numerical] minimum_solid_rigidity`, and the Love solve takes each zone as a layer of its own (see [Pieces and Zones](../worlds/worlds.md#pieces-and-zones)). `"solid"` or `"liquid"` forces the state for an idealized problem, and the layer is then never split.

| Property | Description |
|---|---|
| `is_liquid` | Whether the radial solver treats the whole layer as a liquid: `state = "liquid"`, or `"auto"` with a liquid-only material. |
| `can_change_state` | Whether part of the layer can turn liquid in a solve: `state = "auto"`, `use_melting` on, and a material with both phases. |

The radial solver reads `state`, `is_static`, and `is_incompressible` afresh at every Love solve, so changing them keeps a world's solved structure. The exception is a `state` change that changes `can_change_state`, since the EOS solve looks for zones only in a layer that can change state: the world then forgets its solved structure.

### Changes That Clear a Solve

A layer in a world tells the world when something its EOS solve reads changes, and the world forgets its solved structure (see [Solved State](../worlds/worlds.md#solved-state)):

- The material, the temperature, any of the four material switches, `use_heating`, or `is_volume_fixed`
- The cooling or radiogenics model
- The radii, moved with `set_radii`
- A `state` change that changes `can_change_state`

The rheologies, `use_tides`, `tidal_scale`, `is_static`, `is_incompressible`, and other `state` changes are read by each Love or tidal solve and leave the solved structure standing.

## Cooling and Radiogenics

A cooling model says how heat moves through the layer in a thermal EOS solve, and so the shape of its temperature profile (see [Cooling Models](../../Cooling/cooling_models.md) and [Temperature and Heat Flow](../worlds/worlds.md#temperature-and-heat-flow)). The `cooling` argument and property take a model, a model name (`"off"`, `"conduction"`, `"convection"`), or a config table with a `model` key. The layer shares the model rather than consuming it, so one model can serve several layers. `None` clears it, which holds the layer at one temperature. `cooling_set` reports whether one is attached.

A radiogenics model gives the layer's radiogenic heating (see [Radiogenic Models](../../Radiogenics/radiogenics_models.md)). `set_radiogenics(model)` (or the `radiogenics` argument) moves the C++ model into the layer, so the Python model is left empty and attaching it again raises `ValueError`. `radiogenics_set` reports whether one is attached, and `calc_radiogenic_heating(time, mass)` returns its heating \[W\] at a time \[s\] for a mass \[kg\], 0.0 without one. The model heats the layer during a solve only when `use_heating` is on.

```python
from TidalPy.Radiogenics import make_radiogenics

mantle.cooling = {
    "model": "convection",
    "critical_rayleigh": 1100.0,
}                                                        # Boundary-layer convection
mantle.set_radiogenics(
    make_radiogenics(
        "isotope",
        {"isotopes": "bulk_silicate_earth"}))            # A decaying isotope dataset
mantle.use_heating = True                                # Let the radiogenics heat the layer in a thermal solve
print(mantle.calc_radiogenic_heating(1.4516496e17, 4.0e24))   # [W] at the dataset's reference time, about 2e13
```

## Rheology Overrides

A rheology maps the static modulus and viscosity onto a complex modulus at a forcing frequency (see [Rheology Models](../../Rheology/rheology_models.md)). A layer's material may carry default shear and bulk rheologies on its phase (MatPack rocks carry an Andrade shear rheology), and the layer may override either. The `shear_rheology` and `bulk_rheology` properties return the rheology in effect: the layer's override, else the material's default, else `None`. Setting one sets the override from a model, a model name, or a config table; setting `None` clears the override, so the material's default applies again. With no rheology in effect the complex modulus is the static modulus as a purely real number, which is perfectly elastic. The layer shares the model, so it stays usable.

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

A world's EOS solve populates every layer's profile, and the getters read it back at any radius \[m\]. Each takes a float or an `np.ndarray` and returns a float or an array of the same shape, NaN where nothing is solved. Nothing is stored on a grid: each value comes from the solve's dense output at the radius asked for, evaluated with the layer's material, so a getter and the solve agree exactly.

| Getter | Returns |
|---|---|
| `get_density`, `get_gravity`, `get_pressure` | Density \[kg m$^{-3}$\], gravity \[m s$^{-2}$\], and pressure \[Pa\]. |
| `get_shear_modulus`, `get_bulk_modulus` | Post-melt static shear modulus and adiabatic bulk modulus \[Pa\]. |
| `get_shear_viscosity`, `get_bulk_viscosity` | Post-melt viscosities \[Pa s\]. |
| `get_melt_fraction` | Melt fraction: 0 for a layer that does not melt, 1 for a liquid-only material. |
| `get_static_viscoelastics(radius)` | `(shear_modulus, shear_viscosity, bulk_modulus, bulk_viscosity)` from one evaluation per radius. |
| `get_state(radius)` | Every profile above as a dict. |
| `calc_complex_shear_modulus(radius, frequency)`, `calc_complex_bulk_modulus(radius, frequency)` | The rheology in effect applied to the solved post-melt modulus and viscosity at `radius`, at a frequency \[rad s$^{-1}$\]. |

The one-argument forms `calc_complex_shear_modulus(frequency)` and `calc_complex_bulk_modulus(frequency)` apply the rheology to the layer's material at zero pressure and the layer's own temperature, with no solve needed. `eos_data_populated` and `viscoelastic_populated` report whether a profile and its material state are present. `update_eos_data(radius, density, gravity, pressure)` populates the density, gravity, and pressure profile directly from arrays (radius ascending), for tests and hand-built profiles; such a profile carries no material state, so the other getters read NaN.

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

A layer that belongs to a world reads its profile under the world's call lock, so a getter takes turns with the world's `solve_eos` on another thread (see [Solved State](../worlds/worlds.md#solved-state)), and one call on an array reads every value from one solve. A standalone layer has no lock to take.

## Serialization

`get_config_dict()` returns the layer as the world builder's layer table: `name`, `layer_index`, `radius_inner_m`, `radius_outer_m`, `mass_kg`, `use_tides`, `is_volume_fixed`, `tidal_scale` when one is set, `state`, `is_static`, `is_incompressible`, `temperature_k`, the four material switches, `use_heating`, and the `material` table, the `shear_rheology` and `bulk_rheology` overrides, and the `cooling` and `radiogenics` tables when set. `name` and `radius_inner_m` belong to a standalone layer only (`LAYER_STANDALONE_CONFIG_KEYS`): a world drops them when it nests the layer under its name. `save_config(path)` writes the dict to TOML, and `TidalPy.Structures.build_layer_from_dict` builds a standalone layer from it.

`save_binary(path)` and `load_binary(path)` write and read the layer with its material, rheology overrides, cooling model, and radiogenics model (binary class id 100; see [Binary Serialization](../../Utilities/binary.md)). The solved profile and the tidal heating are not saved. A layer that belongs to a world cannot be loaded in place: load the world, or load into a standalone layer.

```python
from TidalPy.Structures import build_layer_from_dict

config = world.mantle.get_config_dict()  # Builder-valid layer table, material included
twin = build_layer_from_dict(config)     # A standalone layer, owned by no world
print(twin.get_config_dict() == config)  # True
```

## TOML Table

A world file builds each layer from a `[layers.<name>]` table (see the [TOML Schema](../config/toml_schema.md#layer-level-schema) for every key and the [Schema Examples](../config/schema_examples.md) for every form). The table key is the layer's name, the inner radius comes from the layer below, and exactly one of `radius_outer_m`, `radius_fraction`, or `volume_fraction` sets the outer radius. The scalars are the constructor's keyword arguments under their config spelling (`temperature_k` for `temperature`, `mass_kg` for `mass`), and the models are sub-tables:

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

`c_Layer : c_StructureBase` is built from a `c_LayerConfig` (name, index, radii, mass, `use_tides`, `is_volume_fixed`, `tidal_scale`, `state` as a `c_LayerState`, `is_static`, `is_incompressible`, `temperature`, the `c_MaterialSwitches`, and `use_heating`). It holds the material as a `shared_ptr<const c_Material>` (`set_material`, `get_material`), the rheology overrides and cooling model as shared pointers (`set_shear_rheology`, `set_bulk_rheology`, `set_cooling`), and owns its radiogenics model (`set_radiogenics(unique_ptr)`). `get_shear_rheology()` returns the rheology in effect, and `apply_shear_rheology(static_modulus, viscosity, frequency)` is the one place a complex modulus is made. `calc_state(point, out)` evaluates the material with the layer's switches into a `c_MaterialState`. `get_is_liquid()` and `get_can_change_state()` give the state predicates, and `get_eos_fields(field_indices, num_fields, radii, num_radii, values_out)` and `calc_complex_moduli` are the vectorized profile reads. `c_layer_state_from_name` and `c_layer_state_name` convert the state to and from its name.

## References

- Beuthe, M. (2015). Tidal Love numbers of membrane worlds: Europa, Titan, and Co. *Icarus*, 258, 239-266.
