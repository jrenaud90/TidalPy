# Material EOS Models (`Material.eos`)

_Updated: 2026-10-01_

A material equation-of-state model returns a mass density [kg m$^{-3}$]. The analytic models calculate it from the local pressure [Pa]. The interpolated model looks it up by radius [m]. All four use the same call, `calc_density(pressure, temperature=None, radius=0.0)`, so the whole-planet solve does not need to know which model a layer uses.

An EOS model is also the layer's material. It contains, and is the only object that calculates, every frequency-independent property of the layer: the static moduli, the shear law, the static viscosities, and the optional viscosity and partial-melt models (see [The Material](#the-material)).

## Inheritance

```
c_TidalPyBaseClass
  └── c_PhysicsBase
        └── c_MaterialEOSBase  (abstract)
              ├── c_ConstantDensityEOS      aliases "constant", "uniform", "constant_density"
              ├── c_PressureLawEOS          (shared by the two analytic pressure laws)
              │     ├── c_BirchMurnaghanEOS aliases "bm", "birch_murnaghan", "birch-murnaghan"
              │     └── c_VinetEOS          alias "vinet"
              └── c_InterpolatedEOS         aliases "interp", "interpolate", "interpolated"
```

The base declares `calc_density(pressure, temperature, radius)` pure virtual, holds the two thermal parameters, provides `calc_bulk_modulus(pressure, temperature, radius)`, and adds four optional radius-varying getters (see [Interpolated](#interpolated)). The Cython classes mirror the hierarchy: `MaterialEOSBase`, `ConstantDensityEOS`, `BirchMurnaghanEOS`, `VinetEOS`, `InterpolatedEOS`.

## Models

| Model (aliases) | Density from | Parameters |
|---|---|---|
| `constant` (`uniform`, `constant_density`) | nothing; incompressible | `reference_density` |
| `birch_murnaghan` (`bm`) | pressure | `reference_density`, `reference_bulk_modulus`, `bulk_modulus_derivative` |
| `vinet` | pressure | same as Birch-Murnaghan |
| `interpolated` (`interp`, `interpolate`) | radius, by table lookup | `radius`, `density`, and optional viscoelastic tables |

Every model also accepts `thermal_expansion` [K$^{-1}$] and `reference_temperature` [K] so that density can depend on temperature (see [Thermal Terms](#thermal-terms)).

### Constant Density

Returns the same density everywhere. An incompressible body is not realistic but can be a useful diagnostic or applicable to small moons.

### Birch-Murnaghan, Third Order

A finite-strain expansion around a reference state. With the compression $\eta = \rho / \rho_0 = V_0 / V$,

$$P(\eta) = \frac{3}{2} K_0 \left( \eta^{7/3} - \eta^{5/3} \right) \left[ 1 + \frac{3}{4} \left( K_0' - 4 \right) \left( \eta^{2/3} - 1 \right) \right]$$

where $K_0$ is the reference bulk modulus and $K_0'$ its pressure derivative. It is the standard equation of state in mineral physics, fitted to compression experiments across the mantle pressure range. It is the default choice for a silicate or iron layer.

### Vinet

Is derived from a scaled interatomic potential rather than a strain expansion. With $x = (V / V_0)^{1/3} = \eta^{-1/3}$,

$$P(x) = 3 K_0 \frac{1 - x}{x^2} \exp\left[ \frac{3}{2} \left( K_0' - 1 \right) \left( 1 - x \right) \right]$$

The two forms agree closely at modest compression and diverge at high compression, where Vinet is generally the better extrapolation. Fitted $K_0$ and $K_0'$ values are specific to their form: do not use a Birch-Murnaghan fit in the Vinet model.

### Interpolated

Linear interpolation of a sorted radius-to-density table, clamped at both ends. TidalPy uses it for a PREM-like profiles provided by a user for a planet.

The model optionally carries four more radius-varying tables: static shear modulus, static bulk modulus, shear viscosity, and bulk viscosity. They are read with `get_tabulated_shear_modulus(radius)`, `get_tabulated_bulk_modulus(radius)`, `get_tabulated_shear_viscosity(radius)`, and `get_tabulated_bulk_viscosity(radius)`. The analytic models return NaN from these base-class lookups. A tabulated value takes precedence over the material's law or constant (see [The Material](#the-material)). A world TOML that names a `data_file` (in a PREM-like format) fills these tables automatically (see the [TOML schema](../Structures/config/toml_schema.md)).

### Thermal Terms

Each model has a thermal expansivity $\alpha_0$ [K$^{-1}$] and a reference temperature $T_\mathrm{ref}$ [K], the temperature at which $\rho_0$ and $K_0$ apply (default 300 K, where mineral-physics parameters are usually quoted).

Birch-Murnaghan and Vinet add a thermal pressure to their cold pressure law:

$$P(\eta, T) = P_\mathrm{cold}(\eta) + \alpha_0 K_0 \left( T - T_\mathrm{ref} \right)$$

The product $\alpha K_T$ is taken as constant, which is its high-temperature limit (Anderson 1995). The density then comes from inverting the cold law at $P - \alpha_0 K_0 (T - T_\mathrm{ref})$.

The constant and interpolated models have no pressure law, so they scale their density by the thermal expansion directly:

$$\rho(T) = \rho \, \exp\left[ -\alpha_0 \left( T - T_\mathrm{ref} \right) \right]$$

A model is athermal when $\alpha_0 = 0$ (the default) or when no temperature is passed (`None` in Python, a non-finite value in C++).

#### Expansivity Under Compression

The thermal expansivity of rock falls with compression, by several times across the Earth's mantle. A constant $\alpha_0$ therefore makes the adiabat of a thick convecting layer far too steep. The Anderson-Gruneisen parameter $\delta_T = -(\partial \ln \alpha / \partial \ln \rho)_T$ describes this fall. $\delta_T$ itself decreases with compression, $\delta_T = \delta_{T0} (\rho_0 / \rho)^{\kappa}$ (Chopelas and Boehler 1992), and integrating gives

$$\alpha(\rho) = \alpha_0 \exp\left[ \frac{\delta_{T0}}{\kappa} \left( \left( \frac{\rho_0}{\rho} \right)^{\kappa} - 1 \right) \right].$$

For $\kappa = 0$ this is the single power law $\alpha_0 (\rho_0 / \rho)^{\delta_{T0}}$ (Anderson 1967). $\rho_0$ is the model's reference density. The interpolated model has no reference density, so its expansivity stays $\alpha_0$.

| Parameter | Symbol | Default | Mantle silicates |
|---|---|---|---|
| `anderson_gruneisen_parameter` | $\delta_{T0}$ | `0.0` (a constant $\alpha_0$) | 5 to 6 |
| `anderson_gruneisen_exponent` | $\kappa$ | `0.0` | About 1.4 |

With $\delta_{T0}$ = 5.5 and $\kappa$ = 1.4, a rock with $\rho_0$ = 3300 kg m$^{-3}$ keeps 40 percent of $\alpha_0$ at 4000 kg m$^{-3}$ and 13 percent at 5500 kg m$^{-3}$.

```python
from TidalPy.Material.eos import BirchMurnaghanEOS

rock = BirchMurnaghanEOS(
    reference_density=3300.0,
    reference_bulk_modulus=1.3e11,
    bulk_modulus_derivative=4.2,
    thermal_expansion=3.0e-5,
    anderson_gruneisen_parameter=5.5,
    anderson_gruneisen_exponent=1.4
)

rock.calc_thermal_expansion(3300.0)   # [K-1] 3.0e-5, alpha_0 at the reference density
rock.calc_thermal_expansion(5500.0)   # [K-1] about 4.0e-6 under compression
```

A world's thermal solve uses the compressed expansivity for a convecting layer's adiabat, $dT/dr = -\alpha(\rho) g T / c_p$ along the solved density, and for that layer's Rayleigh number. The [Choosing Melting Curves](../PartialMelt/partial_melt_models.md#choosing-melting-curves) table shows the effect on an Earth-like mantle. The density law keeps its thermal pressure $\alpha_0 K_0 (T - T_\mathrm{ref})$, since the product $\alpha K_T$ stays nearly constant as $\alpha$ falls and $K_T$ rises (Anderson 1995).

### Bulk Modulus

`calc_bulk_modulus(pressure, temperature=None, radius=0.0)` returns the isothermal bulk modulus $K_T = \rho \, \partial P / \partial \rho$ [Pa]. Birch-Murnaghan and Vinet evaluate the analytic derivative of their pressure law at the solved compression, so the modulus is consistent with the density and equals $K_0$ at the reference state. The interpolated model returns its bulk table at `radius`. A model with neither returns NaN, and `calc_material_state` then uses the material's `bulk_modulus_static`.

### The Material

The EOS model contains the layer's material parameters:

| Parameter | Units | Default | Meaning |
|---|---|---|---|
| `shear_modulus_static` | Pa | `0.0` | $\mu_0$ of the shear law below. |
| `bulk_modulus_static` | Pa | `0.0` | Bulk modulus of a model with no pressure law or bulk table of its own. |
| `shear_viscosity_static`, `bulk_viscosity_static` | Pa s | NaN (unset) | Used when no viscosity model is attached. |
| `shear_modulus_pressure_derivative` | - | `0.0` | $\mu'_P$ of the shear law. |
| `shear_modulus_temperature_derivative` | Pa K$^{-1}$ | `0.0` | $\mu'_T$ of the shear law. |
| `shear_modulus_reference_temperature` | K | `300.0` | $T_\mathrm{ref}$ of the shear law. |
| `thermal_conductivity` | W m$^{-1}$ K$^{-1}$ | `4.0` | Conductivity $k$ of a conducting layer and of a convecting layer's boundary layers. |
| `heat_capacity` | J kg$^{-1}$ K$^{-1}$ | `1200.0` | Specific heat $c_p$: the adiabat, the diffusivity, and the secular cooling rate. |
| `anderson_gruneisen_parameter` | - | `0.0` | $\delta_{T0}$ of the expansivity's fall with compression. 0 keeps $\alpha_0$. |
| `anderson_gruneisen_exponent` | - | `0.0` | $\kappa$ in $\delta_T = \delta_{T0} (\rho_0 / \rho)^\kappa$. |

The thermal expansivity is the model's `thermal_expansion`, $\alpha_0$ at the reference density, which falls with compression when `anderson_gruneisen_parameter` is set (see [Expansivity Under Compression](#expansivity-under-compression)). The expansivity at the local density sets the adiabatic gradient $\alpha T g / c_p$ and the Rayleigh number of a convecting layer. The density law uses the same $\alpha$ only when it receives a temperature (see step 3 below). `calc_thermal_diffusivity(density)` returns $\kappa = k / (\rho c_p)$.

The material also holds three optional models, attached with `set_shear_viscosity(model)`, `set_bulk_viscosity(model)` (from [`Viscosity`](../Viscosity/viscosity_models.md)) and `set_partial_melt(model)` (from [`PartialMelt`](../PartialMelt/partial_melt_models.md)). A layer has the same three methods, which pass the model to its EOS.

`calc_material_state(pressure, temperature=None, radius=0.0, thermal_density=True)` maps a point onto all of it, in this order:

1. The static shear modulus from the linear law

   $$\mu = \mu_0 + \mu'_P P + \mu'_T \left( T - T_\mathrm{ref} \right)$$

   floored at the config's `minimum_modulus`, and the bulk modulus from `bulk_modulus_static`.
2. The viscosities: the attached viscosity model at the temperature and pressure, else the static viscosity.
3. The density and, where the model defines one, the bulk modulus of the density law. The law receives the temperature only when `thermal_density` is set, which a layer controls with its `use_thermal_eos` switch. The viscosity and partial-melt models always receive it. When the partial-melt model's `density_melt_mixing` switch is on, the density becomes the solid and melt phases mixed by volume, $(1 - \phi)\rho + \phi\rho_l(P)$; the structure iteration uses the same density, so the solved structure and this state agree.
4. A table of an interpolated model replaces the law or constant, and the bulk modulus of a Birch-Murnaghan or Vinet law replaces `bulk_modulus_static`.
5. The partial-melt model, applied to the shear modulus and shear viscosity and, behind their own switches, to the bulk modulus (`bulk_melt_weakening`) and the bulk viscosity (`bulk_viscosity_melt_weakening`), which come after the shear pair because they read the post-melt shear modulus and viscosity. Without a finite temperature this step is skipped, the density is the law's, and `melt_fraction` is NaN. See [Density and Bulk Response](../PartialMelt/partial_melt_models.md#density-and-bulk-response).

It returns a dict of `density`, `melt_fraction`, `shear_modulus`, `bulk_modulus`, `shear_viscosity`, and `bulk_viscosity`, all after the partial-melt step.

The whole-planet EOS solve evaluates this same method as it integrates.

> [!WARNING]
> After a solve, do not call `calc_material_state` to find the state of the planet. Read the world and layer getters (`get_shear_modulus(radius)`, `get_melt_fraction(radius)`, `get_state(radius)`, ...), which report what the solve used. `calc_material_state` evaluates the material at a pressure and temperature of your choosing.

```python
from TidalPy.Material.eos import BirchMurnaghanEOS
from TidalPy.Viscosity import make_viscosity
from TidalPy.PartialMelt import make_partial_melt

rock = BirchMurnaghanEOS(
    reference_density=3300.0, reference_bulk_modulus=1.3e11, bulk_modulus_derivative=4.2,
    shear_modulus_static=6.0e10, shear_modulus_pressure_derivative=1.4,
    shear_modulus_temperature_derivative=-8.0e6)
rock.set_shear_viscosity(make_viscosity("reference", {
    "reference_viscosity_pas": 1.0e19, "reference_temperature_k": 1400.0}))
rock.set_partial_melt(make_partial_melt("henning"))

state = rock.calc_material_state(pressure=2.0e10, temperature=1700.0)
print(state["shear_modulus"], state["shear_viscosity"], state["melt_fraction"])
```

### Pressure Inversion

Both compressible models invert their pressure law for the compression with one shared Newton's method implementation. The slope is exact because the bulk modulus of each law is $K = \eta \, dP/d\eta$, and one evaluation returns both the pressure and $K$. The first guess is the Murnaghan law, $\eta = (1 + K_0' P / K_0)^{1/K_0'}$, which inverts in closed form and follows both laws closely over planetary compressions. Convergence to `invert_rtol` usually takes four or five evaluations. Each evaluation tightens a bracket around the root, and a step that leaves the bracket is replaced by the bracket midpoint.

The laws are monotonic in $\eta$ only over a finite range. Every law turns over in tension. When $K_0' < 4$, the third-order Birch-Murnaghan correction term changes sign at large compression, so $P(\eta)$ also turns over there. The range depends only on $K_0$ and $K_0'$. Each model finds it once, when it is built or loaded, by stepping outward from $\eta = 1$ until $K$ is no longer positive and then bisecting that sign change. A pressure outside the range returns the compression at that end of the range, so the density is continuous in pressure everywhere. The structure solve relies on this: while its central pressure is still a guess, its outer radii can sit far into tension.

Two settings control the iteration. Their defaults come from `[numerical]` in `TidalPy_Configs.toml` (`eos_invert_rtol`, `eos_invert_max_iters`), and each material can override them:

| Setting | Default | Meaning |
|---|---|---|
| `invert_rtol` | `1e-13` | Relative convergence tolerance on the compression. |
| `invert_max_iters` | `60` | Hard iteration cap. A termination safeguard only; convergence normally takes well under ten steps. |

```python
bm = BirchMurnaghanEOS(3500.0, 1.3e11, 4.5, invert_rtol=1e-9, invert_max_iters=80)
bm.invert_rtol, bm.invert_max_iters      # (1e-09, 80)
```

A model built without its own values copies the configured defaults at construction, so a later config change does not affect an existing model. Both settings appear in `get_config_dict()` and survive the binary round trip, so a saved world records the values it was solved with.

### Choosing a Model

- `constant`: an analytic check, a regression test with a known closed-form answer, or a layer thin enough that compression is negligible.
- `birch_murnaghan`: reproducing published mineral-physics parameters, which are most often quoted in this form.
- `vinet`: a fit made in that form, or a layer that reaches compressions where the two forms visibly disagree.
- `interpolated`: a profile that already exists, from a seismic reference model, a mineral-physics or thermal-evolution code, or a previous TidalPy run. It is also the only model that can vary the moduli and viscosities with radius inside a single layer.

## Python API

```python
from TidalPy.Material.eos import (
    ConstantDensityEOS, BirchMurnaghanEOS, VinetEOS, InterpolatedEOS,
    make_material_eos, birch_murnaghan_pressure, vinet_pressure)

bm = BirchMurnaghanEOS(reference_density=3500.0,
                       reference_bulk_modulus=1.3e11,
                       bulk_modulus_derivative=4.5)
density = bm.calc_density(5.0e10)       # [kg/m^3] at 50 GPa
bm.reference_bulk_modulus               # 1.3e11
bm.calc_bulk_modulus(5.0e10)            # [Pa] isothermal bulk modulus at 50 GPa

# A thermal model: hotter material is less dense at the same pressure.
hot = BirchMurnaghanEOS(reference_density=3500.0,
                        reference_bulk_modulus=1.3e11,
                        bulk_modulus_derivative=4.5,
                        thermal_expansion=3.0e-5,
                        reference_temperature=300.0)
hot.calc_density(5.0e10, 2000.0)        # [kg/m^3] at 50 GPa and 2000 K
hot.calc_density(5.0e10)                # no temperature: the athermal density

# Name factory: case-insensitive, aliases accepted.
eos = make_material_eos("vinet", {"reference_density_kg_m3":      3500.0,
                                  "reference_bulk_modulus_pa":    1.3e11,
                                  "bulk_modulus_derivative":      4.5})

# A tabulated profile; density is looked up by radius, not pressure.
prem = InterpolatedEOS(radius=[0.0, 1.0e6, 2.0e6],
                       density=[5000.0, 4000.0, 3000.0])
prem.calc_density(0.0, 0.0, 0.5e6)      # 4500.0

# The forward pressure laws, useful as an inversion cross-check.
birch_murnaghan_pressure(1.2, 1.3e11, 4.5)   # [Pa] at compression eta = 1.2
vinet_pressure(1.2, 1.3e11, 4.5)
```

`make_material_eos(model_name, config=None)` raises `ValueError` for an unrecognized name. Its config keys match those from `get_config_dict()` and, unlike the constructor arguments, carry unit suffixes.

| Config key | Model | Constructor argument |
|---|---|---|
| `reference_density_kg_m3` | all analytic | `reference_density` |
| `reference_bulk_modulus_pa` | Birch-Murnaghan, Vinet | `reference_bulk_modulus` |
| `bulk_modulus_derivative` | Birch-Murnaghan, Vinet | `bulk_modulus_derivative` |
| `invert_rtol`, `invert_max_iters` | Birch-Murnaghan, Vinet | same |
| `thermal_expansion_1_k`, `reference_temperature_k` | all | `thermal_expansion`, `reference_temperature` |
| `anderson_gruneisen_parameter`, `anderson_gruneisen_exponent` | all (a constant $\alpha$ for interpolated) | same |
| `radius_m`, `density_kg_m3` | interpolated | `radius`, `density` |
| `shear_modulus_pa`, `bulk_modulus_pa`, `shear_viscosity_pas`, `bulk_viscosity_pas` | interpolated | `shear_modulus`, `bulk_modulus`, `shear_viscosity`, `bulk_viscosity` |

A key that no material EOS model reads raises `ValueError` naming the closest accepted key, so a misspelling or a missing unit suffix does not silently build a default model.

### Attaching a Model to a `Layer`

```python
from TidalPy.Material.eos import make_material_eos

core.set_eos(make_material_eos("constant", {"reference_density_kg_m3": 11000.0}))
mantle.set_eos(make_material_eos("bm", {"reference_density_kg_m3":   3500.0,
                                        "reference_bulk_modulus_pa": 1.3e11,
                                        "bulk_modulus_derivative":   4.5}))
world.solve_eos(surface_pressure=0.0)
```

`set_eos` moves ownership of the C++ model into the layer and leaves the Python wrapper empty. Every layer class accepts an EOS, and `solve_eos` raises `ValueError` if any layer lacks one. See [Worlds](../Structures/worlds/worlds.md) for the solve itself and its results.

## Serialization

Every model supports the standard interfaces inherited from the base class.

- `get_config_dict()` returns the model name under the key `model`, its parameters (interpolated tables as lists), the material keys (`shear_modulus_static_pa`, `bulk_modulus_static_pa`, the two `*_viscosity_static_pas` when set, and the three shear-law keys), and a sub-table for each attached model (`shear_viscosity`, `bulk_viscosity`, `partial_melt`). `make_material_eos` accepts the dict and builds and attaches the nested models, so a material round-trips through it. The dict is the `[layers.<name>.material]` table of a world TOML.
- `save_config(path)` writes the same content as TOML.
- A model attached to a layer is saved and restored with that layer, in its binary file and in its `get_config_dict()` (under `material`). A world or layer round trip needs no separate handling of materials.

## C++ API

The models are in `material_eos_.hpp` (namespace `tidalpy`, header only). The Cython classes above wrap them, and every C++ consumer (layers attaching an EOS, the whole-planet solve, and binary reconstruction) uses these types directly.

```cpp
#include "material_eos_.hpp"

using namespace tidalpy;

// Defaults come from the struct's own member initializers.
c_MaterialEOSConfig config;
config.reference_density       = 3500.0;
config.reference_bulk_modulus  = 1.3e11;
config.bulk_modulus_derivative = 4.5;

const c_BirchMurnaghanEOS bm(config);
const double density = bm.calc_density(5.0e10, 0.0, 0.0);

// Or through the enum factory, which returns an owning pointer to the base.
std::unique_ptr<c_MaterialEOSBase> eos =
    c_find_material_eos(c_MaterialEOSModel::Vinet, config);
```

### `c_MaterialEOSConfig`

One config struct is shared by every model, and each model reads only the fields it needs. Its member initializers are the single source of default values for both C++ and Python, except the two inversion settings, which take the `[numerical]` defaults when a model is built.

| Field | Used by | Default |
|---|---|---|
| `reference_density` | all | `3500.0` |
| `reference_bulk_modulus` | Birch-Murnaghan, Vinet | `1.0e11` |
| `bulk_modulus_derivative` | Birch-Murnaghan, Vinet | `4.0` |
| `invert_rtol` | Birch-Murnaghan, Vinet | NaN: `[numerical]` `eos_invert_rtol` (`1e-13`) |
| `invert_max_iters` | Birch-Murnaghan, Vinet | -1: `[numerical]` `eos_invert_max_iters` (`60`) |
| `thermal_expansion` | all | `0.0` (athermal) |
| `reference_temperature` | all | `d_EOS_REFERENCE_TEMPERATURE` (`300.0`) |
| `shear_modulus_static`, `bulk_modulus_static` | all | `0.0` |
| `shear_viscosity_static`, `bulk_viscosity_static` | all | NaN (unset) |
| `shear_modulus_pressure_derivative`, `shear_modulus_temperature_derivative` | all | `0.0` |
| `shear_modulus_reference_temperature` | all | `d_EOS_REFERENCE_TEMPERATURE` (`300.0`) |
| `thermal_conductivity`, `heat_capacity` | all | `4.0`, `1200.0` |
| `anderson_gruneisen_parameter`, `anderson_gruneisen_exponent` | all | `0.0`, `0.0` (constant $\alpha$) |
| `radius`, `density` and the four viscoelastic tables | interpolated | empty vectors |

### Classes and Free Functions

All models derive from `c_MaterialEOSBase` and override `calc_density(pressure, temperature, radius)` (see [Inheritance](#inheritance)). The base also has the virtual `calc_density_and_bulk_modulus(pressure, temperature, radius, density, bulk_modulus)`, which returns both from one pressure inversion. Each model has a default constructor and one taking the config. Accessors:

- All models: `get_thermal_expansion()`, `get_reference_temperature()`, `get_reference_density()`, `get_anderson_gruneisen_parameter()`, `get_anderson_gruneisen_exponent()`, `get_expansion_reference_density()` (NaN for interpolated), and `calc_thermal_expansion(density)`.
- Birch-Murnaghan and Vinet: `get_reference_bulk_modulus()`, `get_bulk_modulus_derivative()`, `get_invert_rtol()`, and `get_invert_max_iters()`.
- Interpolated: `get_num_points()`.

| Function | Description |
|---|---|
| `eos_bm_pressure(eta, K0, K0_prime)` | Third-order Birch-Murnaghan pressure [Pa] at compression $\eta$. |
| `eos_vinet_pressure(eta, K0, K0_prime)` | Vinet pressure [Pa] at $\eta$. |
| `eos_bm_bulk_modulus(eta, K0, K0_prime)`, `eos_vinet_bulk_modulus(eta, K0, K0_prime)` | Isothermal bulk modulus $\eta \, dP/d\eta$ [Pa] of each law at $\eta$. |
| `eos_bm_pressure_and_bulk_modulus(eta, K0, K0_prime, pressure, bulk_modulus)`, `eos_vinet_pressure_and_bulk_modulus(...)` | Both values of a law from one evaluation. The four functions above call these. |
| `eos_find_monotonic_range(K0, K0_prime, law_fn, rtol)` | Returns a `c_PressureLawRange`: the compressions between which the law rises, and the pressures there. An end the search does not reach stays unbounded. |
| `eos_invert_eta(pressure_target, K0, K0_prime, law_fn, range, rtol, max_iters)` | Inverts a pressure law for $\eta$. Shared by both compressible models; pass one of the two combined law functions and the model's range. |
| `c_material_eos_model_from_name(name)` | Name or alias to enum, throwing `std::invalid_argument` on an unknown name. |
| `c_find_material_eos(model, config)` | Heap-allocates the model as a `unique_ptr`. A name-string overload is also provided. |
| `c_material_eos_from_binary(stream, force)` | Peeks the binary class id and reconstructs the matching model, used when a layer with an attached EOS is loaded. |

## Adding a New Model

**C++ (`TidalPy/Material/eos/material_eos_.hpp`)**

1. Add any new parameters to `c_MaterialEOSConfig` with sensible defaults.
2. If the law is analytic, add a free function that fills the pressure and the bulk modulus $\eta \, dP/d\eta$ at a compression so the shared `eos_find_monotonic_range` and `eos_invert_eta` can be reused. Give the model an `update_law_range` called from its constructors and from `set_binary_params`. Otherwise, compute the density directly.
3. Add the model class deriving from `c_MaterialEOSBase`: constructors (pass the config to the base so the thermal parameters are stored), `get_*` accessors, the `calc_density` override, a `calc_density_and_bulk_modulus` override if the law defines a bulk modulus, `get_binary_class_id`, and `get_binary_params` / `set_binary_params` appending the law's parameters to the base's material values ([Binary Serialization](../Utilities/binary.md)).
4. Add the enum value, the name and alias branch in `c_material_eos_model_from_name`, and the cases in `c_find_material_eos` and `c_material_eos_from_binary`.

**C++ (`TidalPy/Utilities/binary/binary_.hpp`)**

5. Add a unique `BinaryClassID` in the 60X block.

**Cython (`material_eos.pxd` and `material_eos.pyx`)**

6. Declare the C++ class and the new enum value in the `.pxd`.
7. Add the `cdef class` wrapper with parameter properties and the adoption branch in `make_material_eos`. The config dict comes from the C++ `append_config_entries` override.

**Package, tests, and docs**

8. Export the class from `__init__.py`.
9. Extend `Tests/Test_Material/test_material_eos_01.py`: density evaluation, an inversion cross-check for an analytic law, factory and aliases, config dict, and binary round trip.
10. Document the model here with its formula and references, and add a changelog entry.

No build-system change is needed: `Material.eos.material_eos` is already registered in `cython_extensions.json`.

## References

- Birch, F. (1947). Finite elastic strain of cubic crystals. *Physical Review*, 71(11), 809-824.
- Vinet, P., Ferrante, J., Rose, J. H., and Smith, J. R. (1987). Compressibility of solids. *Journal of Geophysical Research*, 92(B9), 9319-9325.
- Anderson, O. L. (1967). Equation for thermal expansivity in planetary interiors. *Journal of Geophysical Research*, 72(14), 3661-3668. The Anderson-Gruneisen relation for the expansivity.
- Anderson, O. L. (1995). *Equations of State of Solids for Geophysics and Ceramic Science*. Oxford University Press. Thermal pressure and the near constancy of $\alpha K_T$ at high temperature.
- Chopelas, A., and Boehler, R. (1992). Thermal expansivity in the lower mantle. *Geophysical Research Letters*, 19(19), 1983-1986. The decrease of the Anderson-Gruneisen parameter with compression.
- Poirier, J.-P. (2000). *Introduction to the Physics of the Earth's Interior*, second edition. Comparison of the finite-strain and universal forms.
- Dziewonski, A. M., and Anderson, D. L. (1981). Preliminary reference Earth model. *Physics of the Earth and Planetary Interiors*, 25(4), 297-356.
