# Phases and Materials (`Material`)

_Updated: 2026-10-02_

A `Phase` is one phase of a material, a solid or its melt: an [equation-of-state law](material_eos.md), an optional shear-modulus law, optional shear and bulk [viscosity laws](../Viscosity/viscosity_models.md), optional default [rheologies](../Rheology/rheology_models.md), and its thermal conductivity and heat capacity. A `Material` is a solid phase, a liquid phase, or both. With both it melts between a solidus and a liquidus, through a melt-weakening law and optional bulk-mixing laws ([`PartialMelt`](../PartialMelt/partial_melt_models.md)), and with one it is that phase everywhere. `Material.calc_state(pressure, temperature, radius)` maps a point onto every property of the material, as a layer with a given set of physics switches sees it.

`calc_state` is the one place a point is mapped onto a material's properties. The structure and thermal integrations, the thermal network, the radial solver, and every getter of a solved world read it, so they cannot disagree about a property. Phases and materials, like every law, are not changed once built: a change returns a new object, and several layers (and a solve) may hold the same one.

## Inheritance

```
c_TidalPyBaseClass
  └── c_PhysicsBase
        ├── c_PhaseFamily
        │     └── c_Phase      model "phase"
        └── c_MaterialFamily
              └── c_Material   model "material"
```

Each composite holds its sub-models in slots (`c_PhaseComponents`, `c_MaterialComponents`) as shared pointers to immutable laws, and its own scalar parameters through `c_SpecModel`. The Cython classes are `Phase` and `Material`.

## Phases

| Slot | Family | Without one |
|---|---|---|
| `eos` | [equation of state](material_eos.md) | Required; a `Phase` built without one takes a constant-density law at its defaults. |
| `shear_modulus` | [shear modulus](material_eos.md#shear-modulus-laws) | The phase is a fluid: a shear modulus of 0. |
| `shear_viscosity`, `bulk_viscosity` | [viscosity](../Viscosity/viscosity_models.md) | A NaN viscosity. |
| `shear_rheology`, `bulk_rheology` | [rheology](../Rheology/rheology_models.md) | No default; a layer of the phase uses its own rheology or none (elastic). |

The phase's own parameters are its thermal constants. The conductivity and heat capacity follow power laws in temperature about a reference temperature, $k = k_0 (T / T_\mathrm{ref})^{n_k}$ and $c_p = c_{p0} (T / T_\mathrm{ref})^{n_c}$ (MatPack's `ice_ih`, for example, has $k \propto T^{-0.84}$):

| Parameter | Config key | Default | Units |
|---|---|---|---|
| `thermal_conductivity` | `thermal_conductivity_w_mk` | 4.0 | W m$^{-1}$ K$^{-1}$ |
| `conductivity_temperature_exponent` | `conductivity_temperature_exponent` | 0.0 | - |
| `heat_capacity` | `heat_capacity_j_kgk` | 1200.0 | J kg$^{-1}$ K$^{-1}$ |
| `heat_capacity_temperature_exponent` | `heat_capacity_temperature_exponent` | 0.0 | - |
| `thermal_reference_temperature` | `thermal_reference_temperature_k` | 300.0 | K |

The expansivity belongs to the equation-of-state law (see [Thermal Terms](material_eos.md#thermal-terms)), so the density and the adiabat use the same one. `Phase.calc_state(pressure, temperature, radius, use_thermal_expansion=False)` returns the phase's own properties at a point, with no melting.

## Materials

| Slot | Family | Meaning |
|---|---|---|
| `solid` | `Phase` | The solid phase. Its equation of state and default rheologies are the material's own. |
| `liquid` | `Phase` | The liquid phase. With a solid phase it is the melt and needs a `shear_viscosity` law (the melt's viscosity); without one the material is liquid everywhere. |
| `solidus`, `liquidus` | [melting curve](../PartialMelt/partial_melt_models.md#melting-curves) | Required with both phases. Equal curves give a single melting temperature. |
| `weakening` | [melt weakening](../PartialMelt/partial_melt_models.md#melt-weakening) | How the shear modulus and viscosity fall with melt; without one, the solid's values until fully molten. |
| `bulk_modulus_mixing`, `bulk_viscosity_mixing` | [bulk mixing](../PartialMelt/partial_melt_models.md#bulk-mixing) | How melt changes the bulk modulus and bulk viscosity. Without a bulk modulus law, the solid's bulk modulus blends linearly into the liquid's across the weakening law's breakdown band (or steps at full melt with no weakening law); without a bulk viscosity law, the solid's bulk viscosity until fully molten. |

The material's one parameter is `latent_heat` (`latent_heat_j_kg`, default 0.0) \[J kg$^{-1}$\], the latent heat of melting (see [Latent Heat](#latent-heat)).

A material needs at least one phase, and a `Material` built with neither takes a default solid phase. The melting slots need both phases to melt between, and building a material that breaks either rule raises `ValueError`. `can_melt` is true for a material with both phases, and `is_liquid_only` for one with only a liquid phase (a water ocean, a liquid-iron core, or a gas envelope): such a material is liquid, with a melt fraction of 1, whatever the switches.

## State

`calc_state` returns a dict:

| Key | Units | Meaning |
|---|---|---|
| `phase` | - | `"solid"`, `"partial"`, or `"liquid"`. |
| `density` | kg m$^{-3}$ | Density (see [Physics Switches](#physics-switches) for when it mixes the phases). |
| `bulk_modulus` | Pa | Isothermal bulk modulus $K_T$. |
| `adiabatic_bulk_modulus` | Pa | Adiabatic bulk modulus $K_S$, the one the radial solver reads. |
| `thermal_expansion` | K$^{-1}$ | Expansivity $\alpha$. |
| `heat_capacity` | J kg$^{-1}$ K$^{-1}$ | Effective isobaric heat capacity, latent heat included. |
| `thermal_conductivity` | W m$^{-1}$ K$^{-1}$ | Thermal conductivity. |
| `shear_modulus` | Pa | Static (unrelaxed) shear modulus after melt weakening. |
| `shear_viscosity`, `bulk_viscosity` | Pa s | Viscosities after melt weakening and mixing. |
| `melt_fraction` | m$^3$ m$^{-3}$ | Volumetric melt fraction $\phi$; NaN without a finite temperature. |
| `solidus`, `liquidus` | K | The melting curves at the point (at zero pressure without pressure melting); NaN for a material that cannot melt or a layer that does not use melting. |
| `latent_heat_capacity` | J kg$^{-1}$ K$^{-1}$ | The latent heat's share of `heat_capacity` inside a melting range, $L / (T_\mathrm{liq} - T_\mathrm{sol})$; zero elsewhere. Heat diffusion (the thermal diffusivity of a convecting layer's Rayleigh number) uses `heat_capacity` without it. |
| `latent_expansion` | K$^{-1}$ | The latent heat's share of an adiabat's expansivity inside a melting range whose curves follow the pressure (see [Latent Heat](#latent-heat)); zero elsewhere. |

Every value is a float (and `phase` a string) for float inputs, and an array of the broadcast shape otherwise.

## Melting

A material melts only when the layer reading it uses melting (`use_melting`) and the material has both phases. Below the solidus it is its solid phase.

### Melt Fraction

Between the solidus $T_\mathrm{sol}(P)$ and the liquidus $T_\mathrm{liq}(P)$ the melt fraction runs linearly,

$$\phi = \mathrm{clip}\left( \frac{T - T_\mathrm{sol}}{T_\mathrm{liq} - T_\mathrm{sol}},\; 0,\; 1 \right),$$

and above the liquidus the material is its liquid phase, $\phi = 1$. When the two curves coincide (or the liquidus falls below the solidus) the material melts as a step at the solidus: solid at or below it, liquid above it. The ices, `olivine`, `iron`, and `nitrogen_ice` in [MatPack](matpack.md) melt this way. Without a finite temperature there is no melt state: the material is its solid phase with a NaN melt fraction.

### Properties in the Melting Range

Inside the range ($0 < \phi < 1$) the material evaluates both phases at the point and combines them:

- Shear modulus and shear viscosity: the [weakening law](../PartialMelt/partial_melt_models.md#melt-weakening), between the solid's and the liquid's values. Without one, the solid's values until fully molten.
- Bulk modulus and bulk viscosity: the [mixing laws](../PartialMelt/partial_melt_models.md#bulk-mixing) when present. Otherwise the bulk modulus (isothermal and adiabatic) blends linearly from the solid's into the liquid's across the weakening law's breakdown band, as the shear modulus does, so the aggregate reaches the liquid's moduli together (with no weakening law, at full melt), and the bulk viscosity is the solid's until fully molten.
- Density: mixed by volume, $(1 - \phi) \rho_s + \phi \rho_l$, with `use_melt_density`; otherwise the solid's.
- Expansivity, heat capacity, and conductivity: linear in $\phi$ between the two phases, with the latent heat added to the heat capacity.

A fully molten material ($\phi = 1$) takes the liquid phase's shear modulus, shear viscosity, and thermal properties. Its bulk modulus and bulk viscosity are the liquid's, or the mixing laws' values at $\phi = 1$ when it has them, and its density stays the solid's without `use_melt_density`.

### Latent Heat

The latent heat $L$ \[J kg$^{-1}$\] is spread over the melting range as an effective heat capacity,

$$c_{p,\mathrm{eff}} = (1 - \phi) c_{p,s} + \phi c_{p,l} + \frac{L}{T_\mathrm{liq} - T_\mathrm{sol}},$$

so a layer's temperature rises more slowly while it melts. Where the curves follow the pressure (`use_pressure_melting`), the melt fraction also changes with pressure, and an isentrope through the range, $(c_p + L \, \partial\phi/\partial T) \, dT = (\alpha T / \rho - L \, \partial\phi/\partial P) \, dP$, takes the latent heat of that change as an extra expansivity:

$$\alpha_L = \frac{\rho L \left[ (1 - \phi) \, dT_\mathrm{sol}/dP + \phi \, dT_\mathrm{liq}/dP \right]}{(T_\mathrm{liq} - T_\mathrm{sol}) \, T}.$$

This is the reported `latent_expansion`. Adiabats (in the thermal integration and in a convecting layer's profile) run at $dT/dr = -(\alpha + \alpha_L) g T / c_p$. Without it, an adiabat inside the range would be shallower than the melting curve it nears while one outside is steeper, and the two would meet on the curve. $\alpha_L$ is zero outside the range and when `use_pressure_melting` is off.

A material that melts as a step has no range to spread its latent heat over. The boundary between a layer's solid and liquid zones carries it instead (a Stefan condition): the world's thermal network adds the latent heat that boundary absorbs per kelvin of the layer's temperature to the layer's heat capacity, and the solve reports it as `layer_latent_capacity` (see [Worlds](../Structures/worlds/worlds.md)).

## Physics Switches

A material is evaluated with four switches. A layer passes its own, so one material can serve layers that use more or less of it. Each switch off, the default, is the simpler and faster case.

| Switch | Off | On |
|---|---|---|
| `use_thermal_expansion` | The density ignores the temperature. | The density (and a pressure law's thermal pressure) follows the temperature. The expansivity shapes the adiabat either way. |
| `use_melting` | The solid phase only; the melt fraction is 0 and the solidus and liquidus are NaN. | The liquid phase, the melt fraction, weakening, and mixing take part. |
| `use_pressure_melting` | The melting curves are read at zero pressure. | The melting curves follow the local pressure, and `latent_expansion` can be nonzero. |
| `use_melt_density` | The solid's density, even when fully molten. | The phases' densities mixed by volume. |

A layer has two more switches that a material does not see: `use_heating` (the world's heat sources act in the layer) and `use_tides` (the layer dissipates tidal energy; off, it responds elastically). Its `state` option (`"auto"`, `"solid"`, `"liquid"`) says how the radial solver treats it. With `"auto"`, `use_melting` on, and a material that melts, the world's EOS solve splits the layer into solid and liquid zones, liquid wherever the post-melt shear modulus falls to `[numerical] minimum_solid_rigidity` (default 1e-6) times the world's $\rho g R$ or below. `Layer.calc_state(pressure, temperature=None, radius=None)` evaluates the layer's material with the layer's switches, at its own temperature and mid-radius unless given.

### Behavior at the Limits

- A phase with no shear-modulus law has a shear modulus of 0, and one with no viscosity law a NaN viscosity. A viscoelastic rheology applied to a NaN viscosity returns NaN, so a phase that a layer tides through needs a viscosity law unless its rheology is elastic.
- A shear-modulus law's value is floored at `[numerical] minimum_modulus`.
- A NaN melting temperature leaves the material solid.
- A melting curve fitted over a limited pressure range can fall to zero past its end (ice Ih's above 415 MPa), and the material then reads as fully molten there.
- The switches do nothing for a material with one phase.

## Python API

```python
from TidalPy.Material import Material, Phase, make_eos, make_shear_modulus
from TidalPy.PartialMelt import make_melt_weakening, make_melting_curve
from TidalPy.Rheology import Maxwell
from TidalPy.Viscosity import make_viscosity

# The solid phase, from models
rock = Phase(
    eos=make_eos(
        "birch_murnaghan",
        {"reference_density_kg_m3": 3300.0,
         "reference_bulk_modulus_pa": 1.3e11,
         "bulk_modulus_derivative": 4.2,
         "thermal_expansion_1_k": 3.0e-5}),
    shear_modulus=make_shear_modulus(
        "linear",
        {"shear_modulus_pa": 6.0e10,
         "pressure_derivative": 1.4,
         "temperature_derivative_pa_k": -8.0e6}),
    shear_viscosity=make_viscosity(
        "reference",
        {"reference_viscosity_pas": 1.0e19,
         "reference_temperature_k": 1400.0}),
    shear_rheology=Maxwell(),
    thermal_conductivity=3.3)

# The melt, from config tables
melt = Phase(
    eos={"model": "murnaghan", "reference_density_kg_m3": 2750.0},
    shear_viscosity={"model": "constant", "reference_viscosity_pas": 0.1})

mantle_rock = Material(
    solid=rock,
    liquid=melt,
    solidus=make_melting_curve("constant", {"temperature_k": 1600.0}),
    liquidus=make_melting_curve("constant", {"temperature_k": 2000.0}),
    weakening=make_melt_weakening("henning"),
    latent_heat=4.0e5)

# A quarter molten at 2 GPa and 1700 K
state = mantle_rock.calc_state(
    2.0e9,
    1700.0,
    use_melting=True)
print(state["phase"], state["melt_fraction"], state["shear_modulus"], state["heat_capacity"])
```

Each slot takes a model, a config table with a `model` key, or a model name for that model's defaults. `Material`'s phase slots take a `Phase` or a phase table. The scalar parameters are keywords, by argument name or config key. `config=` takes a whole phase or material as one nested table, the form `get_config_dict()` returns and a TOML material table holds, and `make_phase(config)` and `make_material(config)` do the same. A misnamed slot or parameter raises `ValueError` naming the closest accepted one.

| Member | Returns | Description |
|---|---|---|
| `calc_state(pressure, temperature=nan, radius=nan, *, use_thermal_expansion=False, use_melting=False, use_pressure_melting=False, use_melt_density=False)` | `dict` | The [state](#state) at a point, with the given switches (`Material`). |
| `calc_state(pressure, temperature=nan, radius=nan, use_thermal_expansion=False)` | `dict` | The phase's own state at one point, without melting (`Phase`). |
| `calc_melting_range(pressure, use_pressure_melting=True)` | `tuple` | `(solidus, liquidus)` \[K\] at a pressure \[Pa\]; NaN for a material that cannot melt. |
| `solid`, `liquid`, `solidus`, `liquidus`, `weakening`, `bulk_modulus_mixing`, `bulk_viscosity_mixing` | model or `None` | The material's slots, each wrapped in its class. |
| `eos`, `shear_modulus`, `shear_viscosity`, `bulk_viscosity`, `shear_rheology`, `bulk_rheology` | model or `None` | The phase's slots. |
| `can_melt`, `is_liquid_only` | `bool` | Whether the material has both phases, or only a liquid one. |
| `with_parameters(**changes)` | composite | A copy with some scalar parameters changed (`latent_heat`, or a phase's thermal constants); the slots are kept. |
| `replace(**changes)` | `Material` | A copy with some slots or parameters replaced; `None` removes a slot. |
| `parameters`, `get_parameter(name)`, `get_parameter_info()` | | The scalar parameters, as for any law. |
| `get_config_dict()`, `save_config(path)` | `dict`, - | The nested table, and that table written as TOML. |

### Changing a Material

A material is not changed in place. `with_parameters` changes its scalar parameters and `replace` its slots, each returning a new material and leaving this one as it was. To change a law inside a phase, build a new phase and replace the old one, or load the material with an override (see [MatPack](matpack.md#overrides-and-presets)).

```python
from TidalPy.Material import load_material

peridotite = load_material("peridotite")

# No latent heat
no_latent = peridotite.with_parameters(latent_heat=0.0)

# Fischer and Spohn weakening in place of Henning's
spohn_peridotite = peridotite.replace(weakening="spohn")

# A lower conductivity in the solid phase
cooler_solid = peridotite.solid.with_parameters(thermal_conductivity=2.5)
low_conductivity = peridotite.replace(solid=cooler_solid)

# A law inside a phase changed through an override
constant_viscosity = load_material(
    "peridotite",
    solid={"shear_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e18}})

print(no_latent.latent_heat, spohn_peridotite.weakening.model_name, low_conductivity.solid.thermal_conductivity)
```

### Vectorized Evaluation

`Material.calc_state` takes a float or an array for each of the pressure, temperature, and radius, and broadcasts them together. The switches are scalars.

```python
import numpy as np

from TidalPy.Material import load_material

peridotite = load_material("peridotite")

# The melt fraction and latent expansion along a 2100 K isotherm
pressure = np.linspace(0.0, 1.0e10, 50)   # [Pa]
state = peridotite.calc_state(
    pressure,
    2100.0,
    use_melting=True,
    use_pressure_melting=True)
melt_fraction    = state["melt_fraction"]
latent_expansion = state["latent_expansion"]   # [K-1]
phases           = state["phase"]              # an object array of "solid", "partial", "liquid"
```

At the C++ level `calc_state_vectorize(pressure, temperature, radius, switches, out_states)` fills one `c_MaterialState` per point.

### Attaching a Material to a `Layer`

```python
from TidalPy.Material import load_material
from TidalPy.Structures.layers import Layer

# A MatPack name, a material table, or a Material
mantle = Layer(
    "mantle",
    1,
    9.0e5,
    1.82e6,
    material="peridotite",
    temperature=1650.0,
    use_melting=True,
    use_pressure_melting=True)

mantle.material = load_material("peridotite", latent_heat_j_kg=3.0e5)   # Replace it later
print(mantle.calc_state(1.0e9)["solidus"])   # [K] at 1 GPa, the layer's switches and temperature
```

The layer shares the material rather than copying it. A world built from a TOML file reads each layer's `material` key or `[layers.<name>.material]` table and the switches beside it. See the [TOML schema](../Structures/config/toml_schema.md) and [Worlds](../Structures/worlds/worlds.md).

## C++ API

`TidalPy/Material/material_.hpp` (namespace `tidalpy`, header only) holds both composites.

| Type | Contents |
|---|---|
| `c_MaterialSwitches` | `use_thermal_expansion`, `use_melting`, `use_pressure_melting`, `use_melt_density` (all false by default). |
| `c_MaterialPhase` | `Solid`, `Partial`, `Liquid`. |
| `c_PhaseState` | One phase's properties at a point. |
| `c_MaterialState` | The material's state at a point, every field of [State](#state). |
| `c_PhaseComponents`, `c_MaterialComponents` | The slots as shared pointers; `set(slot, model)` checks the model's family and throws `std::invalid_argument` for a wrong one. |

`c_Material` provides:

- `calc_state(point, switches, out)`, the full state.
- `calc_thermal(point, switches, out)`, only what the thermal integration needs (no moduli, viscosities, or weakening).
- `calc_density(point, switches)`, what the structure integration asks for.
- `calc_rigidity_margin(point, switches, minimum_shear_modulus)`, $\ln(\mu / \mu_\mathrm{min})$ for the solve's state events.
- `calc_melting_range(pressure, switches, solidus, liquidus)` and `calc_state_vectorize`.
- The accessors `get_components`, `get_solid`, `get_liquid`, `get_base_phase` (the solid, or the liquid of a liquid-only material), `get_can_melt`, `get_is_liquid_only`, and `get_latent_heat`.

Every evaluation is `const`, `noexcept`, and allocation-free. `c_Phase` provides `calc_phase_state(point, thermal, out)`, `calc_phase_thermal`, `calc_density`, and its slot accessors. `c_make_phase(params, components)` and `c_make_material(params, components)` build them from slots (filling a missing equation of state or phase with the defaults), and `c_phase_from_binary` and `c_material_from_binary` read them back.

## Serialization

| Call | Result |
|---|---|
| `get_config_dict()` | A phase: its thermal parameters plus one table per filled slot. A material: `latent_heat_j_kg`, a `solid` and a `liquid` table, and a `melting` table holding `solidus`, `liquidus`, `weakening`, `bulk_modulus_mixing`, and `bulk_viscosity_mixing`. `make_material` (or `load_material`) accepts it, so a material round-trips. |
| `save_config(path)` | That table written as TOML; it is a [MatPack](matpack.md) file without the metadata keys. |
| `save_binary(path)` / `load_binary(path)` | The parameters, then one optional record per slot in slot order. |

A layer's material is saved and restored with the layer, in its binary record and in its config dict (as its `material` table), so a world round-trip needs no separate handling of materials.

## References

- Solomatov, V. S. (2000). Fluid dynamics of a terrestrial magma ocean. In *Origin of the Earth and Moon*, University of Arizona Press, 323-338. Melting ranges and latent heat in planetary interiors.
- Monteux, J., Andrault, D., and Samuel, H. (2016). On the cooling of a deep terrestrial magma ocean. *Earth and Planetary Science Letters*, 448, 140-149. Latent heat across a pressure-dependent melting range.
- Turcotte, D. L., and Schubert, G. (2002). *Geodynamics*, second edition. The Stefan problem of a moving solidification front.
