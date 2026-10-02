# Equation-of-State and Shear-Modulus Laws (`Material.laws`)

_Updated: 2026-10-02_

An equation-of-state law maps a pressure \[Pa\], temperature \[K\], and radius \[m\] onto a phase's density \[kg m$^{-3}$\], its isothermal and adiabatic bulk moduli \[Pa\], and its thermal expansivity \[K$^{-1}$\]. The analytic laws calculate the density from the pressure, and the interpolated law looks it up by radius. A shear-modulus law maps the same point onto the phase's static (unrelaxed) shear modulus \[Pa\]. Every law takes the same point, so a phase, and the whole-planet solve that reads it, never needs to know which law it holds.

Each law is one slot of a [phase](materials.md). The frequency dependence of the moduli is the [rheology's](../Rheology/rheology_models.md) job, and melt weakening is the material's.

## Inheritance

```
c_TidalPyBaseClass
  └── c_PhysicsBase
        ├── c_EOSBase  (abstract)
        │     ├── c_ConstantEOS            aliases "constant", "uniform", "constant_density"
        │     ├── c_BirchMurnaghanEOS      aliases "birch_murnaghan", "bm", "birch-murnaghan"
        │     ├── c_VinetEOS               alias "vinet"
        │     ├── c_MurnaghanEOS           alias "murnaghan"
        │     ├── c_PolytropeEOS           alias "polytrope"
        │     ├── c_ModifiedPolytropeEOS   aliases "modified_polytrope", "seager"
        │     └── c_InterpolatedEOS        aliases "interpolate", "interp", "interpolated"
        └── c_ShearModulusBase  (abstract)
              ├── c_ConstantShearModulus       aliases "constant", "const"
              ├── c_LinearShearModulus         alias "linear"
              └── c_InterpolatedShearModulus   aliases "interpolate", "interp", "interpolated"
```

Birch-Murnaghan and Vinet share one template, `c_PressureLawEOSModel`, which holds the pressure-law inversion. Every concrete law derives from `c_SpecModel` (see [C++ API](#c-api)). The Cython classes mirror the hierarchy: `EOSBase` with `ConstantEOS`, `BirchMurnaghanEOS`, `VinetEOS`, `MurnaghanEOS`, `PolytropeEOS`, `ModifiedPolytropeEOS`, and `InterpolatedEOS`, and `ShearModulusBase` with `ConstantShearModulus`, `LinearShearModulus`, and `InterpolatedShearModulus`.

## Inputs and Result

Every law takes a `c_ThermoPoint`: the pressure \[Pa\], the temperature \[K\] (NaN for a law evaluated without one), and the radius \[m\] (read only by the laws tabulated in radius). An equation-of-state law also takes a `thermal` flag, which says whether the density sees the temperature. A layer passes its `use_thermal_expansion` switch.

An equation-of-state law returns a `c_EOSPoint`:

| Field | Units | Meaning |
|---|---|---|
| `density` | kg m$^{-3}$ | Density at the point. |
| `bulk_modulus` | Pa | Isothermal bulk modulus $K_T = \rho \, \partial P / \partial \rho$; NaN for a law that gives none. |
| `adiabatic_bulk_modulus` | Pa | $K_S$, the modulus a tidal (adiabatic) deformation sees (see [Adiabatic Bulk Modulus](#adiabatic-bulk-modulus)). |
| `thermal_expansion` | K$^{-1}$ | Expansivity $\alpha$ at the point's density, reported whether or not the density sees the temperature, since the adiabat needs it either way. |

The Python `calc_eos` returns the same four fields as a dict.

```python
from TidalPy.Material.laws import BirchMurnaghanEOS

rock = BirchMurnaghanEOS(
    reference_density=3300.0,
    reference_bulk_modulus=1.3e11,
    bulk_modulus_derivative=4.2,
    thermal_expansion=3.0e-5)

# 50 GPa and 2000 K, with the density following the temperature
state = rock.calc_eos(
    5.0e10,
    2000.0,
    thermal=True)
print(state["density"], state["bulk_modulus"], state["thermal_expansion"])
```

## Models

| Model (aliases) | Density from | Parameters |
|---|---|---|
| `constant` (`uniform`, `constant_density`) | nothing; incompressible | `reference_density`, `bulk_modulus` |
| `birch_murnaghan` (`bm`) | pressure, by inversion | `reference_density`, `reference_bulk_modulus`, `bulk_modulus_derivative`, `invert_rtol`, `invert_max_iters` |
| `vinet` | pressure, by inversion | as Birch-Murnaghan |
| `murnaghan` | pressure, closed form | `reference_density`, `reference_bulk_modulus`, `bulk_modulus_derivative` |
| `polytrope` | pressure, closed form | `polytropic_constant`, `polytropic_index` |
| `modified_polytrope` (`seager`) | pressure, closed form | `reference_density`, `polytrope_coefficient`, `polytrope_exponent`, `anderson_gruneisen_parameter`, `anderson_gruneisen_exponent` |
| `interpolate` (`interp`, `interpolated`) | radius, by table lookup | `radius`, `density`, optional `bulk_modulus` |

Every law also carries the three thermal parameters of [Thermal Terms](#thermal-terms).

### Constant

Returns the same density everywhere, $\rho_0 \exp[-\alpha_0 (T - T_\mathrm{ref})]$, with a constant bulk modulus `bulk_modulus` for the radial solver. An incompressible body is not realistic but can be a useful diagnostic or applicable to small moons.

### Birch-Murnaghan, Third Order

A finite-strain expansion around a reference state. With the compression $\eta = \rho / \rho_0 = V_0 / V$,

$$P(\eta) = \frac{3}{2} K_0 \left( \eta^{7/3} - \eta^{5/3} \right) \left[ 1 + \frac{3}{4} \left( K_0' - 4 \right) \left( \eta^{2/3} - 1 \right) \right]$$

where $K_0$ is the reference bulk modulus and $K_0'$ its pressure derivative. It is the standard equation of state in mineral physics, fitted to compression experiments across the mantle pressure range, and the usual choice for a silicate or iron phase.

### Vinet

Is derived from a scaled interatomic potential rather than a strain expansion. With $x = (V / V_0)^{1/3} = \eta^{-1/3}$,

$$P(x) = 3 K_0 \frac{1 - x}{x^2} \exp\left[ \frac{3}{2} \left( K_0' - 1 \right) \left( 1 - x \right) \right]$$

The two forms agree closely at modest compression and diverge at high compression, where Vinet is generally the better extrapolation. Fitted $K_0$ and $K_0'$ values are specific to their form: do not use a Birch-Murnaghan fit in the Vinet law.

### Murnaghan

The bulk modulus rises linearly with pressure, $K = K_0 + K_0' P$, which inverts in closed form:

$$\rho = \rho_0 \left( 1 + \frac{K_0' P}{K_0} \right)^{1/K_0'}$$

$K_0' = 0$, and any pressure in tension, gives $\rho_0 e^{P / K_0}$, which joins the law smoothly at $P = 0$. Because $\rho / (d\rho/dP) = K$, a fully liquid layer of a Murnaghan phase is neutrally stratified under the bulk modulus the tidal equations see. It suits melts and liquids over modest pressures.

### Polytrope

$P = K \rho^{1 + 1/n}$ with the polytropic constant $K$ and index $n$, so $\rho = (P / K)^{n/(n+1)}$ and $K_T = (1 + 1/n) P$. A polytrope is a barotrope: it takes no temperature, and its thermal parameters only shape the adiabat. The density is zero where the pressure is not positive, the surface of a gas envelope. $K \approx 2 \times 10^5$ in SI units at $n = 1$ fits Jupiter.

### Modified Polytrope

$\rho = \rho_0 + c P^n$ for $P > 0$ and $\rho_0$ otherwise (Seager et al. 2007), a fit to the cold compression of planetary materials to TPa pressures, with $K_T = \rho / (c \, n \, P^{n-1})$. It is isothermal: its thermal parameters only shape the adiabat, with an expansivity that can fall with compression (see [Expansivity Under Compression](#expansivity-under-compression)). The defaults are Seager et al.'s fit for iron.

### Interpolated

Linear interpolation of a radius-to-density table, held at the end values beyond it, and scaled by $\exp[-\alpha_0 (T - T_\mathrm{ref})]$. An optional `bulk_modulus` table gives the bulk modulus the same way. Without one, the law reports NaN. A layer holding an interpolated law must keep its volume, since the table is in radius. A world TOML that names a `data_file` (in a PREM-like format) gives each of its layers interpolated laws built from the file's columns (see the [TOML schema](../Structures/config/toml_schema.md)).

### Thermal Terms

Every law carries the same three thermal parameters, and reports the expansivity its own density has, $\alpha = -(1/\rho)(\partial \rho / \partial T)_P$, so the density, the adiabat, and the Rayleigh number agree:

| Parameter | Config key | Default | Meaning |
|---|---|---|---|
| `thermal_expansion` | `thermal_expansion_1_k` | `0.0` | $\alpha_0$ \[K$^{-1}$\] at the reference state. |
| `reference_temperature` | `reference_temperature_k` | `300.0` | $T_\mathrm{ref}$ \[K\], where $\rho_0$ and $K_0$ apply (where mineral-physics parameters are usually quoted). |
| `gruneisen_parameter` | `gruneisen_parameter` | `0.0` | $\gamma$ in $K_S = K_T (1 + \alpha \gamma T)$; 0 makes $K_S$ equal $K_T$. |

Birch-Murnaghan, Vinet, and Murnaghan add a thermal pressure to their cold law:

$$P(\eta, T) = P_\mathrm{cold}(\eta) + \alpha_0 K_0 \left( T - T_\mathrm{ref} \right)$$

The product $\alpha K_T$ is taken as constant, which is its high-temperature limit (Anderson 1995). The density then comes from the cold law at $P - \alpha_0 K_0 (T - T_\mathrm{ref})$, and its expansivity is $\alpha = \alpha_0 K_0 / K_T$, with $K_T$ the law's bulk modulus at that density. The constant and interpolated laws have no pressure law, so they scale their density by $\exp[-\alpha_0 (T - T_\mathrm{ref})]$ instead, an expansivity of $\alpha_0$. The polytrope and the modified polytrope ignore the temperature.

A law gives its athermal density when $\alpha_0 = 0$ (the default), when the `thermal` flag is off (a layer with `use_thermal_expansion` off), or when the temperature is not finite. Its expansivity is still reported.

### Expansivity Under Compression

The thermal expansivity of rock falls with compression, by several times across the Earth's mantle. A constant $\alpha_0$ would make the adiabat of a thick convecting layer far too steep. For the pressure laws the fall follows from the thermal pressure: $\alpha = \alpha_0 K_0 / K_T$ falls as compression stiffens the law, the Anderson-Gruneisen behavior with $\delta_T = -(\partial \ln \alpha / \partial \ln \rho)_T$ equal to the law's $\partial \ln K_T / \partial \ln \rho$ (about $K_0'$). The bundled `earth_simple` mantle, a Birch-Murnaghan law, has $\alpha$ = 5.2e-5 K$^{-1}$ at its top and 8.8e-6 K$^{-1}$ at its base, where $K_T$ is about six times larger.

The modified polytrope is a barotrope, so its density gives no expansivity. Its expansivity shapes the adiabat alone and falls with compression through an Anderson-Gruneisen parameter $\delta_T$, which itself decreases with compression, $\delta_T = \delta_{T0} (\rho_0 / \rho)^{\kappa}$ (Chopelas and Boehler 1992). Integrating gives

$$\alpha(\rho) = \alpha_0 \exp\left[ \frac{\delta_{T0}}{\kappa} \left( \left( \frac{\rho_0}{\rho} \right)^{\kappa} - 1 \right) \right],$$

the single power law $\alpha_0 (\rho_0 / \rho)^{\delta_{T0}}$ for $\kappa = 0$ (Anderson 1967), with $\delta_{T0}$ from `anderson_gruneisen_parameter` and $\kappa$ from `anderson_gruneisen_exponent` (both 0 by default, which keeps $\alpha_0$). The other laws take no such parameters, since their expansivity is their density's.

A world's thermal solve uses the expansivity for a convecting layer's adiabat, $dT/dr = -\alpha g T / c_p$ along the solved structure (with the latent term of a [melting range](materials.md#latent-heat) added), and for that layer's Rayleigh number.

### Adiabatic Bulk Modulus

A tide deforms a planet faster than heat can diffuse, so the deformation is adiabatic and sees $K_S$, not $K_T$. With a Gruneisen parameter $\gamma$,

$$K_S = K_T \left( 1 + \alpha \gamma T \right).$$

For a silicate at 2000 K with $\alpha$ = 3e-5 K$^{-1}$ and $\gamma$ = 1.2, $K_S / K_T$ is about 1.07. A material reports both, and the radial solver reads the adiabatic one.

### Pressure Inversion

Birch-Murnaghan and Vinet invert their pressure law for the compression with one shared Newton's method implementation. The slope is exact because the bulk modulus of each law is $K = \eta \, dP/d\eta$, and one evaluation returns both the pressure and $K$. The first guess is the Murnaghan law, $\eta = (1 + K_0' P / K_0)^{1/K_0'}$, which inverts in closed form and follows both laws closely over planetary compressions. Convergence to `invert_rtol` usually takes four or five evaluations. Each evaluation tightens a bracket around the root, and a step that leaves the bracket is replaced by the bracket midpoint.

The laws are monotonic in $\eta$ only over a finite range. Every law turns over in tension. When $K_0' < 4$, the third-order Birch-Murnaghan correction term changes sign at large compression, so $P(\eta)$ also turns over there. The range depends only on $K_0$ and $K_0'$. Each law finds it once, when it is built or loaded, by stepping outward from $\eta = 1$ until $K$ is no longer positive and then bisecting that sign change. A pressure outside the range returns the compression at that end of the range, so the density is continuous in pressure everywhere. The structure solve relies on this: while its central pressure is still a guess, its outer radii can sit far into tension.

Two parameters control the iteration. Left unset, each takes its value from `[numerical]` in `TidalPy_Configs.toml` (`eos_invert_rtol`, `eos_invert_max_iters`) when the law is built, so a later config change does not affect an existing law:

| Parameter | Unset value | Configured default | Meaning |
|---|---|---|---|
| `invert_rtol` | NaN | `1e-13` | Relative convergence tolerance on the compression. |
| `invert_max_iters` | `-1` | `60` | Hard iteration cap. A termination safeguard only; convergence normally takes well under ten steps. |

```python
from TidalPy.Material.laws import BirchMurnaghanEOS

# A looser inversion for a quick survey
survey_rock = BirchMurnaghanEOS(
    reference_density=3500.0,
    reference_bulk_modulus=1.3e11,
    bulk_modulus_derivative=4.5,
    invert_rtol=1.0e-9,
    invert_max_iters=80)
print(survey_rock.invert_rtol, survey_rock.invert_max_iters)   # 1e-09 80
```

### Shear-Modulus Laws

| Model (aliases) | Shear modulus $\mu$ \[Pa\] | Parameters |
|---|---|---|
| `constant` (`const`) | $\mu_0$ | `shear_modulus` |
| `linear` | $\mu_0 + \mu'_P (P - P_\mathrm{ref}) + \mu'_T (T - T_\mathrm{ref})$ | `shear_modulus`, `pressure_derivative`, `temperature_derivative`, `reference_pressure`, `reference_temperature` |
| `interpolate` (`interp`, `interpolated`) | linear in radius between the points of a `radius` and `shear_modulus` table, held at the end values beyond it | `radius`, `shear_modulus` |

The linear law leaves its temperature term out when the temperature is not finite. A phase floors any law's value at `[numerical] minimum_modulus`, since a steep temperature derivative can take the law negative far from its fit. A phase with no shear-modulus law is a fluid (a shear modulus of 0).

### Parameters

Each parameter carries two names: the constructor keyword, which also reads as an attribute, and the config key used in a TOML table, a `make_eos` or `make_shear_modulus` config dict, and `get_config_dict()`. A dimensional config key ends in its unit while the code name does not. Either name is accepted wherever a parameter is given.

**Equation-of-state laws**

| Parameter | Config key | Default | Units | Used by |
|---|---|---|---|---|
| `reference_density` | `reference_density_kg_m3` | 3500.0 (Murnaghan 2750.0, modified polytrope 8300.0) | kg m$^{-3}$ | All but polytrope and interpolated |
| `bulk_modulus` | `bulk_modulus_pa` | 1.0e11 | Pa | Constant |
| `reference_bulk_modulus` | `reference_bulk_modulus_pa` | 1.0e11 (Murnaghan 2.0e10) | Pa | Birch-Murnaghan, Vinet, Murnaghan |
| `bulk_modulus_derivative` | `bulk_modulus_derivative` | 4.0 (Murnaghan 5.0) | - | Birch-Murnaghan, Vinet, Murnaghan |
| `invert_rtol`, `invert_max_iters` | same | unset | - | Birch-Murnaghan, Vinet |
| `polytropic_constant` | `polytropic_constant` | 2.0e5 | Pa (m$^3$ kg$^{-1}$)$^{1 + 1/n}$ | Polytrope |
| `polytropic_index` | `polytropic_index` | 1.0 | - | Polytrope |
| `polytrope_coefficient` | `polytrope_coefficient` | 0.00349 | kg m$^{-3}$ Pa$^{-n}$ | Modified polytrope |
| `polytrope_exponent` | `polytrope_exponent` | 0.528 | - | Modified polytrope |
| `radius` | `radius_m` | `[0.0]` | m | Interpolated |
| `density` | `density_kg_m3` | `[3500.0]` | kg m$^{-3}$ | Interpolated |
| `bulk_modulus` | `bulk_modulus_pa` | `[]` (none) | Pa | Interpolated |
| thermal parameters | see [Thermal Terms](#thermal-terms) | | | All |

**Shear-modulus laws**

| Parameter | Config key | Default | Units | Used by |
|---|---|---|---|---|
| `shear_modulus` | `shear_modulus_pa` | 5.0e10 (`[5.0e10]` for interpolated) | Pa | All |
| `pressure_derivative` | `pressure_derivative` | 0.0 | - | Linear |
| `temperature_derivative` | `temperature_derivative_pa_k` | 0.0 | Pa K$^{-1}$ | Linear |
| `reference_pressure` | `reference_pressure_pa` | 0.0 | Pa | Linear |
| `reference_temperature` | `reference_temperature_k` | 300.0 | K | Linear |
| `radius` | `radius_m` | `[0.0]` | m | Interpolated |

`get_parameter_info()` on any law lists its parameters with their keys, defaults, bounds, and descriptions.

### Behavior at the Limits

- Birch-Murnaghan and Vinet hold the compression at the end of their monotonic range for a pressure beyond it (see [Pressure Inversion](#pressure-inversion)). Murnaghan continues into tension as $\rho_0 e^{P / K_0}$.
- The polytrope gives zero density, and zero bulk modulus, where the pressure is not positive. The modified polytrope gives $\rho_0$ there, with a bulk modulus of 0 for $n < 1$.
- The tabulated laws hold their end values beyond their tables.
- A non-finite temperature gives the athermal density.
- A parameter outside its bounds (a non-positive density or reference temperature, for example) raises `ValueError` when the law is built.

### Choosing a Model

- `constant`: an analytic check, a regression test with a known closed-form answer, or a layer thin enough that compression is negligible.
- `birch_murnaghan`: reproducing published mineral-physics parameters, which are most often quoted in this form.
- `vinet`: a fit made in that form, or a phase that reaches compressions where the two forms visibly disagree.
- `murnaghan`: a melt or a liquid over modest pressures, the usual liquid phase of a melting material.
- `polytrope`: a gas-giant envelope.
- `modified_polytrope`: a cold core or a super-Earth interior compressed to TPa pressures, beyond the range of the finite-strain fits.
- `interpolate`: a profile that already exists, from a seismic reference model, a mineral-physics or thermal-evolution code, or a previous TidalPy run. With the interpolated shear-modulus and viscosity laws, it is the only way to vary the properties with radius inside a single layer independently of pressure and temperature.

## Python API

```python
import numpy as np

from TidalPy.Material import make_eos, make_shear_modulus
from TidalPy.Material.laws import InterpolatedEOS, LinearShearModulus

# Name factory: case-insensitive, aliases accepted, parameters by config key.
melt = make_eos(
    "murnaghan",
    {"reference_density_kg_m3": 2750.0,
     "reference_bulk_modulus_pa": 2.0e10,
     "bulk_modulus_derivative": 5.0})
print(melt.calc_density(5.0e9))           # [kg m-3] at 5 GPa

# A tabulated profile; the density is looked up by radius, not pressure.
profile = InterpolatedEOS(
    radius=[0.0, 1.0e6, 2.0e6],
    density=[5000.0, 4000.0, 3000.0])
print(profile.calc_density(0.0, radius=0.5e6))   # 4500.0

# A shear modulus that rises with pressure and falls with temperature
shear_law = LinearShearModulus(
    shear_modulus=6.0e10,
    pressure_derivative=1.4,
    temperature_derivative=-8.0e6)
print(shear_law.calc_shear_modulus(2.0e9, 1500.0))   # [Pa] 5.32e10
same_law = make_shear_modulus(
    "linear",
    {"shear_modulus_pa": 6.0e10,
     "pressure_derivative": 1.4,
     "temperature_derivative_pa_k": -8.0e6})
```

Constructors take every parameter their law uses, by its argument name or config key, as keywords or positionally in the order `get_parameter_info()` lists them, with the defaults from the tables above. `make_eos(model_name, config=None)` and `make_shear_modulus(model_name, config=None)` resolve a name or alias case-insensitively and build that law from `config`. Absent keys take the law's defaults. An unknown name, a key the law does not read, or a value outside a parameter's bounds raises `ValueError` naming the closest accepted name or key.

| Member | Returns | Description |
|---|---|---|
| `calc_eos(pressure, temperature=nan, radius=nan, thermal=False)` | `dict` | `density`, `bulk_modulus`, `adiabatic_bulk_modulus`, and `thermal_expansion` (equation-of-state laws). |
| `calc_density(pressure, temperature=nan, radius=nan, thermal=False)` | `float` or `np.ndarray` | The density alone \[kg m$^{-3}$\]. |
| `calc_shear_modulus(pressure, temperature=nan, radius=nan)` | `float` or `np.ndarray` | The static shear modulus \[Pa\] (shear-modulus laws). |
| `model_name` | `str` | The canonical model name. |
| `parameters` | `dict` | Every parameter by argument name. |
| `get_parameter(name)` | value | One parameter by argument name or config key. |
| `get_parameter_info()` | `list` of `dict` | Each parameter's name, config key, kind, default, bounds, and description. |
| `with_parameters(**changes)` | law | A new law with some parameters changed; this one is unchanged. |
| `get_config_dict()` | `dict` | `model` plus every parameter under its config key, ready for the factory. An unset (NaN) `invert_rtol` is left out. |
| `save_config(path)` | - | Writes that dict to a TOML file. |

Parameters also read as attributes under either name (`law.reference_density`, `law.reference_density_kg_m3`). Laws are not changed once built, so one law can serve several phases.

`eos_model_names()` and `shear_modulus_model_names()` list the canonical names, `canonical_eos_name(name)` and `canonical_shear_modulus_name(name)` resolve an alias, and `eos_config_keys(name)` and `shear_modulus_config_keys(name)` list the config keys one law reads. These live in `TidalPy.Material.laws`.

### Factory Internals

At the C++ level each family has one registry, `c_eos_registry()` and `c_shear_modulus_registry()`, listing each law's names (canonical first, then aliases), its binary class id, and its constructor. `c_find_eos(name, params)` and `c_find_shear_modulus(name, params)` build a law from a `c_ParamMap` (config key to a list of values) and throw `std::invalid_argument` for an unknown name or parameter. The Python factories pass the config dict through to them and wrap the result in the matching class.

### Vectorized Evaluation

`calc_eos`, `calc_density`, and `calc_shear_modulus` take a float or an array for each of the pressure, temperature, and radius, and broadcast them together. Floats give floats, and anything else gives arrays of the broadcast shape.

```python
import numpy as np

from TidalPy.Material.laws import BirchMurnaghanEOS

rock = BirchMurnaghanEOS(
    reference_density=3300.0,
    reference_bulk_modulus=1.3e11,
    bulk_modulus_derivative=4.2,
    thermal_expansion=3.0e-5,
    gruneisen_parameter=1.2)

# A pressure sweep at one temperature
sweep = rock.calc_eos(
    np.linspace(0.0, 1.0e11, 50),
    2000.0,
    thermal=True)
expansivity = sweep["thermal_expansion"]          # [K-1], alpha0 K0 / K_T, falls with compression
modulus_ratio = sweep["adiabatic_bulk_modulus"] / sweep["bulk_modulus"]
```

At the C++ level `calc_eos_vectorize(pressure, temperature, radius, thermal, ...)` and `calc_shear_modulus_vectorize(pressure, temperature, radius, out)` loop over vectors, each holding one value per point or a single value used at every point, and fill caller-supplied output vectors. Mismatched lengths raise `ValueError` from Python.

### Attaching a Law to a `Phase`

A law reaches a layer through a phase of the layer's material.

```python
from TidalPy.Material import Material, Phase, make_eos
from TidalPy.Structures.layers import Layer

rock_phase = Phase(
    eos=make_eos(
        "vinet",
        {"reference_density_kg_m3": 3300.0,
         "reference_bulk_modulus_pa": 1.3e11,
         "bulk_modulus_derivative": 4.2}),
    shear_modulus={"model": "constant", "shear_modulus_pa": 6.0e10})   # A config table works too

mantle = Layer(
    "mantle",
    0,
    0.0,
    1.0e6,
    material=Material(solid=rock_phase))
```

The phase shares the law rather than copying it. The declarative form is a `[layers.<name>.material.solid.eos]` or `[layers.<name>.material.solid.shear_modulus]` table in a world's TOML. See [Phases and Materials](materials.md) and the [TOML schema](../Structures/config/toml_schema.md).

## C++ API

The laws are in `TidalPy/Material/laws/eos_law_.hpp` and `shear_modulus_law_.hpp` (namespace `tidalpy`, header only), with the Birch-Murnaghan and Vinet pressure laws and their inversion in `pressure_laws_.hpp`. The Cython classes wrap them, and every C++ consumer (the phase, the whole-planet solve, and binary reconstruction) uses these types directly.

```cpp
#include "eos_law_.hpp"

using namespace tidalpy;

// Parameters by config key; absent keys take the law's defaults.
c_ParamMap params;
params["reference_density_kg_m3"]   = {3300.0};
params["reference_bulk_modulus_pa"] = {1.3e11};
params["bulk_modulus_derivative"]   = {4.2};
std::unique_ptr<c_EOSBase> rock = c_find_eos("vinet", params);

c_ThermoPoint point;
point.pressure    = 5.0e10;
point.temperature = 2000.0;
c_EOSPoint result;
rock->calc_eos(point, true, result);   // true: the density sees the temperature
```

`c_EOSBase` provides `calc_eos(point, thermal, out)`, `calc_density(point, thermal)`, `calc_thermal_expansion(density)`, `calc_thermal_pressure(temperature, thermal)`, `get_pressure_law_range()` (unbounded for a law without one), `get_reference_temperature()`, and `calc_eos_vectorize`. A law implements the protected `p_calc_law(point, temperature_offset, out)`, which fills the density and isothermal bulk modulus, and overrides `p_expansion_reference_density()` and `p_thermal_pressure(temperature_offset)` where it has them. The base adds the expansivity and the adiabatic modulus. `c_ShearModulusBase` declares `calc_shear_modulus(point)` pure virtual and provides `calc_shear_modulus_vectorize`.

| Function | Description |
|---|---|
| `eos_bm_pressure(eta, K0, K0_prime)`, `eos_vinet_pressure(eta, K0, K0_prime)` | Birch-Murnaghan and Vinet pressure \[Pa\] at compression $\eta$. |
| `eos_bm_bulk_modulus(eta, K0, K0_prime)`, `eos_vinet_bulk_modulus(eta, K0, K0_prime)` | Isothermal bulk modulus $\eta \, dP/d\eta$ \[Pa\] of each law at $\eta$. |
| `eos_bm_pressure_and_bulk_modulus(eta, K0, K0_prime, pressure, bulk_modulus)`, `eos_vinet_pressure_and_bulk_modulus(...)` | Both values of a law from one evaluation. The four functions above call these. |
| `eos_find_monotonic_range(K0, K0_prime, law_fn, rtol)` | Returns a `c_PressureLawRange`: the compressions between which the law rises, and the pressures there. An end the search does not reach stays unbounded. |
| `eos_invert_eta(pressure_target, K0, K0_prime, law_fn, range, rtol, max_iters)` | Inverts a pressure law for $\eta$; pass one of the two combined law functions and the law's range. |
| `c_find_eos(name, params)`, `c_find_shear_modulus(name, params)` | Build a law by name or alias as a `std::unique_ptr`. |
| `c_eos_from_binary(stream, force)`, `c_shear_modulus_from_binary(stream, force)` | Peek a record's class id and rebuild the matching law, used when a phase is loaded. |
| `c_eos_canonical_name(name)`, `c_eos_model_names()` and the shear-modulus pair | Name resolution. |

## Serialization

| Call | Result |
|---|---|
| `get_config_dict()` | `model` plus every parameter under its config key (tables as lists), except an unset (NaN) `invert_rtol`, which is left out. The factory accepts it, so a law round-trips through it. |
| `save_config(path)` | That dict written as TOML. |
| `save_binary(path)` / `load_binary(path)` | The law's TidalPy binary record, its parameters written by key; a key missing from a record reads at its default. |

A law in a phase is saved and restored with that phase, its material, and its layer, in both the binary record and the config dict.

## Adding a New Model

To add an equation-of-state law named `Foo` (a shear-modulus law is the same in `shear_modulus_law_.hpp`):

**C++ (`TidalPy/Material/laws/eos_law_.hpp`)**

1. Add `c_FooEOS : public c_SpecModel<c_FooEOS, c_EOSBase>`: its `parameter_specs()` table (argument name, config key, member, default, bounds, one-line description per parameter) with `p_append_thermal_specs(rows)` at the end, `C_CLASS_ID`, two constructors that call `p_initialize`, and `p_calc_law`. Override `p_expansion_reference_density` and `p_thermal_pressure` if the law has a reference density or adds a thermal pressure, `p_validate` for checks across parameters, and `p_update_derived` for values cached from them (the pressure-law range, say).
2. Add one row to `c_eos_registry()` with the law's names and aliases.

**C++ (`TidalPy/Utilities/binary/binary_.hpp`)**

3. Reserve a unique `BinaryClassID::FooEOSLaw`. Equation-of-state laws occupy the 61X range and shear-modulus laws the 62X range.

**Cython (`laws.pyx`)**

4. Add `cdef class FooEOS(EOSBase)` with a docstring and `MODEL_NAME = "foo"`, and include it in the `ModelFamily` list.

**Package, tests, and documentation**

5. Export `FooEOS` from `TidalPy/Material/laws/__init__.py`.
6. Add its physics tests to `Tests/Test_Material/test_eos_laws_01.py` (the density against a closed form, and an inversion cross-check for a law that inverts). The generic tests in `Tests/Test_Utilities/Test_Classes/test_spec_models_01.py` cover its parameters, config, binary record, and errors without changes.
7. Document the law here with its formula and references.

No build-system change is needed: `Material.laws.laws` is already registered in `cython_extensions.json`.

## References

- Birch, F. (1947). Finite elastic strain of cubic crystals. *Physical Review*, 71(11), 809-824.
- Vinet, P., Ferrante, J., Rose, J. H., and Smith, J. R. (1987). Compressibility of solids. *Journal of Geophysical Research*, 92(B9), 9319-9325.
- Murnaghan, F. D. (1944). The compressibility of media under extreme pressures. *Proceedings of the National Academy of Sciences*, 30(9), 244-247.
- Chandrasekhar, S. (1939). *An Introduction to the Study of Stellar Structure*. University of Chicago Press. Polytropes.
- Seager, S., Kuchner, M., Hier-Majumder, C. A., and Militzer, B. (2007). Mass-radius relationships for solid exoplanets. *The Astrophysical Journal*, 669, 1279-1297. The modified polytrope.
- Anderson, O. L. (1967). Equation for thermal expansivity in planetary interiors. *Journal of Geophysical Research*, 72(14), 3661-3668. The Anderson-Gruneisen relation for the expansivity.
- Anderson, O. L. (1995). *Equations of State of Solids for Geophysics and Ceramic Science*. Oxford University Press. Thermal pressure, the near constancy of $\alpha K_T$ at high temperature, and $K_S / K_T$.
- Chopelas, A., and Boehler, R. (1992). Thermal expansivity in the lower mantle. *Geophysical Research Letters*, 19(19), 1983-1986. The decrease of the Anderson-Gruneisen parameter with compression.
- Poirier, J.-P. (2000). *Introduction to the Physics of the Earth's Interior*, second edition. Comparison of the finite-strain and universal forms.
- Dziewonski, A. M., and Anderson, D. L. (1981). Preliminary reference Earth model. *Physics of the Earth and Planetary Interiors*, 25(4), 297-356.
