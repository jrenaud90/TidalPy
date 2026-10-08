# Equation-of-State and Shear-Modulus Laws (`Material.laws`)

_Updated: 2026-10-07_

An equation-of-state law maps a pressure \[Pa\], temperature \[K\], and radius \[m\] onto a phase's density \[kg m$^{-3}$\], its isothermal and adiabatic bulk moduli \[Pa\], and its thermal expansivity \[K$^{-1}$\]. The analytic laws calculate the density from the pressure; the interpolated law looks it up by radius. A shear-modulus law maps the same point onto the phase's static (unrelaxed) shear modulus \[Pa\]. Each law fills one slot of a [phase](materials.md). The frequency dependence of the moduli belongs to the [rheology](../Rheology/rheology_models.md), and melt weakening to the material.

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

The `thermal` flag says whether the density sees the temperature; a layer passes its `use_thermal_expansion` switch.

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

Constructors take parameters by argument name or config key, as keywords or positionally in `get_parameter_info()` order, with the defaults of [Parameters](#parameters). `make_eos(model_name, config=None)` and `make_shear_modulus(model_name, config=None)` resolve a name or alias case-insensitively; absent keys take the defaults. An unknown name or key, or a value outside its bounds, raises `ValueError` naming the closest accepted one.

| Member | Returns | Description |
|---|---|---|
| `calc_eos(pressure, temperature=nan, radius=nan, thermal=False)` | `dict` | `density` \[kg m$^{-3}$\]; `bulk_modulus` $K_T = \rho \, \partial P / \partial \rho$ \[Pa\] (NaN for a law that gives none); `adiabatic_bulk_modulus` $K_S$ \[Pa\], the modulus a tidal deformation sees; and `thermal_expansion` $\alpha$ \[K$^{-1}$\], reported whether or not the density sees the temperature, since the adiabat needs it either way. |
| `calc_density(pressure, temperature=nan, radius=nan, thermal=False)` | `float` or `np.ndarray` | The density alone \[kg m$^{-3}$\]. |
| `calc_shear_modulus(pressure, temperature=nan, radius=nan)` | `float` or `np.ndarray` | The static shear modulus \[Pa\] (shear-modulus laws). |
| `save_config(path)` | - | Writes `get_config_dict()` to a TOML file. |

Every law also has the members shared by all physics models (`model_name`, `parameters`, `get_parameter`, `get_parameter_info`, `with_parameters`, `get_config_dict`; see [`PhysicsBase`](../Utilities/classes.md#physicsbase)), and its parameters read as attributes under either name (`law.reference_density`, `law.reference_density_kg_m3`). Laws are not changed once built, so one law can serve several phases. In `TidalPy.Material.laws`, `eos_model_names()`, `canonical_eos_name(name)`, and `eos_config_keys(name)` (and their `shear_modulus_` versions) list the names, resolve an alias, and list a law's config keys.

### Vectorized Evaluation

`calc_eos`, `calc_density`, and `calc_shear_modulus` broadcast array pressure, temperature, and radius together. Floats give floats; mismatched lengths raise `ValueError`.

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

### Attaching a Law to a `Phase`

A law reaches a layer through a phase of the layer's material. `Phase(eos=..., shear_modulus=...)` takes a law, a config table with a `model` key, or a model name, and shares the law rather than copying it (see [Building a Material](materials.md#building-a-material)). In a world's TOML the law is a `[layers.<name>.material.solid.eos]` or `[layers.<name>.material.solid.shear_modulus]` table (see the [TOML schema](../Structures/config/toml_schema.md)).

## Models

| Model (aliases) | Python class | Density from | Use for |
|---|---|---|---|
| `constant` (`uniform`, `constant_density`) | `ConstantEOS` | nothing; incompressible | Analytic checks, regression tests, layers too thin to compress, small moons. |
| `birch_murnaghan` (`bm`, `birch-murnaghan`) | `BirchMurnaghanEOS` | pressure, by inversion | Published mineral-physics parameters, most often quoted in this form. |
| `vinet` | `VinetEOS` | pressure, by inversion | A fit made in this form, or compressions where the two forms visibly disagree. |
| `murnaghan` | `MurnaghanEOS` | pressure, closed form | A melt or liquid over modest pressures. |
| `polytrope` | `PolytropeEOS` | pressure, closed form | A gas-giant envelope. |
| `modified_polytrope` (`seager`) | `ModifiedPolytropeEOS` | pressure, closed form | A cold core or super-Earth interior at TPa pressures, beyond the finite-strain fits. |
| `interpolate` (`interp`, `interpolated`) | `InterpolatedEOS` | radius, by table lookup | An existing profile: a seismic reference model, another code's output, or a previous TidalPy run. |

Every law also carries the three parameters of [Thermal Terms](#thermal-terms). Each Python class wraps the C++ class of the same name with a `c_` prefix.

### Constant

The density is $\rho_0 \exp[-\alpha_0 (T - T_\mathrm{ref})]$ everywhere, with a constant `bulk_modulus` for the radial solver. An incompressible body is not realistic but is a useful diagnostic.

### Birch-Murnaghan, Third Order

A finite-strain expansion around a reference state. With the compression $\eta = \rho / \rho_0 = V_0 / V$,

$$P(\eta) = \frac{3}{2} K_0 \left( \eta^{7/3} - \eta^{5/3} \right) \left[ 1 + \frac{3}{4} \left( K_0' - 4 \right) \left( \eta^{2/3} - 1 \right) \right]$$

where $K_0$ is the reference bulk modulus and $K_0'$ its pressure derivative. It is the standard equation of state in mineral physics, fitted to compression experiments across the mantle pressure range, and the usual choice for a silicate or iron phase.

### Vinet

Derived from a scaled interatomic potential rather than a strain expansion. With $x = (V / V_0)^{1/3} = \eta^{-1/3}$,

$$P(x) = 3 K_0 \frac{1 - x}{x^2} \exp\left[ \frac{3}{2} \left( K_0' - 1 \right) \left( 1 - x \right) \right]$$

The two forms agree closely at modest compression and diverge at high compression, where Vinet is generally the better extrapolation. Fitted $K_0$ and $K_0'$ are specific to their form: do not use a Birch-Murnaghan fit in the Vinet law.

### Murnaghan

The bulk modulus rises linearly with pressure, $K = K_0 + K_0' P$, which inverts in closed form:

$$\rho = \rho_0 \left( 1 + \frac{K_0' P}{K_0} \right)^{1/K_0'}$$

$K_0' = 0$, and any pressure in tension, gives $\rho_0 e^{P / K_0}$, which joins the law smoothly at $P = 0$. Because $\rho / (d\rho/dP) = K$, a fully liquid Murnaghan layer is neutrally stratified under the bulk modulus the tidal equations see.

### Polytrope

$P = K \rho^{1 + 1/n}$ with the polytropic constant $K$ and index $n$, so $\rho = (P / K)^{n/(n+1)}$ and $K_T = (1 + 1/n) P$. A polytrope is a barotrope: it takes no temperature, and its thermal parameters only shape the adiabat. The density is zero where the pressure is not positive, the surface of a gas envelope. $K \approx 2 \times 10^5$ in SI units at $n = 1$ fits Jupiter.

### Modified Polytrope

$\rho = \rho_0 + c P^n$ for $P > 0$ and $\rho_0$ otherwise (Seager et al. 2007), a fit to the cold compression of planetary materials to TPa pressures, with $K_T = \rho / (c \, n \, P^{n-1})$. It is isothermal: its thermal parameters only shape the adiabat (see [Expansivity Under Compression](#expansivity-under-compression)). The defaults are Seager et al.'s fit for iron.

### Interpolated

Linear interpolation of a radius-to-density table, held at the end values beyond it, and scaled by $\exp[-\alpha_0 (T - T_\mathrm{ref})]$. An optional `bulk_modulus` table gives the bulk modulus the same way; without one the law reports NaN. Its layer must keep its volume, since the table is in radius. A world TOML that names a PREM-like `data_file` gives each layer interpolated laws built from the file's columns (see the [TOML schema](../Structures/config/toml_schema.md)). With the interpolated shear-modulus and viscosity laws, this is the only way to vary the properties with radius inside one layer independently of pressure and temperature.

### Thermal Terms

Every law carries three thermal parameters and reports the expansivity its own density has, $\alpha = -(1/\rho)(\partial \rho / \partial T)_P$, so the density, the adiabat, and the Rayleigh number agree:

| Parameter | Config key | Default | Meaning |
|---|---|---|---|
| `thermal_expansion` | `thermal_expansion_1_k` | `0.0` | $\alpha_0$ \[K$^{-1}$\] at the reference state. |
| `reference_temperature` | `reference_temperature_k` | `300.0` | $T_\mathrm{ref}$ \[K\], where $\rho_0$ and $K_0$ apply (where mineral-physics parameters are usually quoted). |
| `gruneisen_parameter` | `gruneisen_parameter` | `0.0` | $\gamma$ in $K_S = K_T (1 + \alpha \gamma T)$; 0 makes $K_S$ equal $K_T$. |

Birch-Murnaghan, Vinet, and Murnaghan add a thermal pressure to their cold law:

$$P(\eta, T) = P_\mathrm{cold}(\eta) + \alpha_0 K_0 \left( T - T_\mathrm{ref} \right)$$

taking $\alpha K_T$ as constant, its high-temperature limit (Anderson 1995). The density comes from the cold law at $P - \alpha_0 K_0 (T - T_\mathrm{ref})$, with expansivity $\alpha = \alpha_0 K_0 / K_T$ at that density. The constant and interpolated laws scale their density by $\exp[-\alpha_0 (T - T_\mathrm{ref})]$ instead, an expansivity of $\alpha_0$. The polytrope and the modified polytrope ignore the temperature.

A law gives its athermal density when $\alpha_0 = 0$ (the default), when the `thermal` flag is off, or when the temperature is not finite. Its expansivity is still reported.

### Expansivity Under Compression

The thermal expansivity of rock falls several times across the Earth's mantle; a constant $\alpha_0$ would make a thick convecting layer's adiabat far too steep. For the pressure laws the fall follows from the thermal pressure: $\alpha = \alpha_0 K_0 / K_T$ falls as compression stiffens the law. This is Anderson-Gruneisen behavior, with $\delta_T = -(\partial \ln \alpha / \partial \ln \rho)_T$ equal to the law's $\partial \ln K_T / \partial \ln \rho$ (about $K_0'$). The bundled `earth_simple` mantle, a Birch-Murnaghan law, has $\alpha$ = 5.2e-5 K$^{-1}$ at its top and 8.8e-6 K$^{-1}$ at its base, where $K_T$ is about six times larger.

The modified polytrope's density gives no expansivity. Its expansivity shapes the adiabat alone and falls with compression through an Anderson-Gruneisen parameter that itself decreases with compression, $\delta_T = \delta_{T0} (\rho_0 / \rho)^{\kappa}$ (Chopelas and Boehler 1992). Integrating gives

$$\alpha(\rho) = \alpha_0 \exp\left[ \frac{\delta_{T0}}{\kappa} \left( \left( \frac{\rho_0}{\rho} \right)^{\kappa} - 1 \right) \right],$$

the power law $\alpha_0 (\rho_0 / \rho)^{\delta_{T0}}$ for $\kappa = 0$ (Anderson 1967), with $\delta_{T0}$ from `anderson_gruneisen_parameter` and $\kappa$ from `anderson_gruneisen_exponent` (both 0 by default, which keeps $\alpha_0$).

A world's thermal solve uses the expansivity for a convecting layer's adiabat, $dT/dr = -\alpha g T / c_p$ (plus the latent term of a [melting range](materials.md#latent-heat)), and for its Rayleigh number.

### Adiabatic Bulk Modulus

A tide deforms a planet faster than heat can diffuse, so the deformation is adiabatic and sees $K_S$, not $K_T$. With a Gruneisen parameter $\gamma$,

$$K_S = K_T \left( 1 + \alpha \gamma T \right).$$

For a silicate at 2000 K with $\alpha$ = 3e-5 K$^{-1}$ and $\gamma$ = 1.2, $K_S / K_T$ is about 1.07. The radial solver reads $K_S$.

### Pressure Inversion

Birch-Murnaghan and Vinet invert their pressure law for the compression with a bracketed Newton's method started from the Murnaghan law. Convergence to `invert_rtol` usually takes four or five evaluations.

The laws rise monotonically in $\eta$ only over a finite range, which depends only on $K_0$ and $K_0'$. Every law turns over in tension, and when $K_0' < 4$ the third-order Birch-Murnaghan correction changes sign at large compression, so $P(\eta)$ turns over there too. A pressure outside the range returns the compression at that end, so the density is continuous in pressure everywhere. The structure solve relies on this: while its central pressure is still a guess, its outer radii can sit far into tension.

Left unset, the two iteration parameters take their values from `[numerical]` in `TidalPy_Configs.toml` (`eos_invert_rtol`, `eos_invert_max_iters`) when the law is built; a later config change does not affect an existing law.

| Parameter | Unset value | Configured default | Meaning |
|---|---|---|---|
| `invert_rtol` | NaN | `1e-13` | Relative convergence tolerance on the compression. |
| `invert_max_iters` | `-1` | `60` | Hard iteration cap, a safeguard only. |

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

| Model (aliases) | Python class | Shear modulus $\mu$ \[Pa\] |
|---|---|---|
| `constant` (`const`) | `ConstantShearModulus` | $\mu_0$ |
| `linear` | `LinearShearModulus` | $\mu_0 + \mu'_P (P - P_\mathrm{ref}) + \mu'_T (T - T_\mathrm{ref})$ |
| `interpolate` (`interp`, `interpolated`) | `InterpolatedShearModulus` | linear in radius between the points of a `radius` and `shear_modulus` table |

A phase floors any law's value at `[numerical] minimum_modulus`, since a steep temperature derivative can take the law negative far from its fit. A phase with no shear-modulus law is a fluid (a shear modulus of 0).

## Parameters

Each parameter has a constructor keyword (also an attribute) and a config key (TOML tables, factory configs, `get_config_dict()`); a dimensional config key ends in its unit. Either name is accepted.

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
| `anderson_gruneisen_parameter`, `anderson_gruneisen_exponent` | same | 0.0 | - | Modified polytrope |
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

## Behavior at the Limits

- Pressure beyond a Birch-Murnaghan or Vinet law's monotonic range: the compression at that end of the range (see [Pressure Inversion](#pressure-inversion)). Murnaghan continues into tension as $\rho_0 e^{P / K_0}$.
- Non-positive pressure: zero density and bulk modulus for the polytrope; $\rho_0$ for the modified polytrope, with a bulk modulus of 0 for $n < 1$.
- Outside a table: the end values.
- Non-finite temperature: the athermal density; the linear shear law drops its temperature term.
- Parameter outside its bounds (a non-positive density or reference temperature, for example): `ValueError` when the law is built.

## Serialization

`get_config_dict()` returns `model` plus every parameter under its config key (tables as lists, an unset `invert_rtol` left out), which the factory accepts, so a law round-trips. `save_binary(path)` and `load_binary(path)` write and read the law's TidalPy binary record, parameters by key; a key missing from a record reads at its default. A law in a phase is saved and restored with that phase, its material, and its layer.

## C++ API

The laws are header only, in namespace `tidalpy`, in `TidalPy/Material/laws`: `eos_law_.hpp`, `shear_modulus_law_.hpp`, and `pressure_laws_.hpp` (the Birch-Murnaghan and Vinet pressure laws and their inversion).

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

A `c_ThermoPoint` holds the pressure, temperature (NaN for none), and radius; a `c_EOSPoint` the four fields of `calc_eos`. `c_EOSBase` provides `calc_eos(point, thermal, out)`, `calc_density(point, thermal)`, `calc_thermal_pressure(temperature, thermal)`, `get_pressure_law_range()` (unbounded for a law without one), `get_reference_temperature()`, and `calc_eos_vectorize(pressure, temperature, radius, thermal, ...)`. `c_ShearModulusBase` provides `calc_shear_modulus(point)` and `calc_shear_modulus_vectorize(pressure, temperature, radius, out)`. The vectorized calls take vectors of one value per point or a single value for every point, and fill caller-supplied output vectors.

`c_find_eos(name, params)` and `c_find_shear_modulus(name, params)` build a law by name or alias as a `std::unique_ptr` and throw `std::invalid_argument` for an unknown name or parameter. `c_eos_from_binary(stream, force)` and `c_shear_modulus_from_binary(stream, force)` rebuild the law a binary record names, and `c_eos_canonical_name(name)`, `c_eos_model_names()`, and the shear-modulus pair resolve names.

The pressure laws are free functions of `(eta, K0, K0_prime)`: `eos_bm_pressure` and `eos_vinet_pressure` give $P$ \[Pa\], `eos_bm_bulk_modulus` and `eos_vinet_bulk_modulus` give $\eta \, dP/d\eta$ \[Pa\], and `eos_bm_pressure_and_bulk_modulus` and `eos_vinet_pressure_and_bulk_modulus` fill both from one evaluation. `eos_find_monotonic_range(K0, K0_prime, law_fn, rtol)` returns a `c_PressureLawRange` (the compressions between which the law rises, and the pressures there; an end the search does not reach stays unbounded), and `eos_invert_eta(pressure_target, K0, K0_prime, law_fn, range, rtol, max_iters)` inverts a law for $\eta$ given a combined law function and its range.

## Adding a New Model

An equation-of-state law `Foo` (a shear-modulus law is the same, in `shear_modulus_law_.hpp`) needs:

1. `c_FooEOS : public c_SpecModel<c_FooEOS, c_EOSBase>` in `eos_law_.hpp`: a `parameter_specs()` table (argument name, config key, member, default, bounds, description) ending in `p_append_thermal_specs(rows)`, `C_CLASS_ID`, two constructors that call `p_initialize`, and `p_calc_law(point, temperature_offset, out)`, which fills the density and $K_T$ (the base adds the expansivity and $K_S$). Override `p_thermal_pressure(temperature_offset)` for a thermal pressure and `p_calc_thermal_expansion(law_point)` when the expansivity is not the constant $\alpha_0$ (a pressure law returns `p_thermal_pressure_expansion(K0, law_point)`, which is $\alpha_0 K_0 / K_T$); `p_validate` and `p_update_derived` take cross-parameter checks and cached values.
2. A row in `c_eos_registry()`, and a `BinaryClassID::FooEOSLaw` in `Utilities/binary/binary_.hpp` (EOS laws 61X, shear-modulus laws 62X).
3. `cdef class FooEOS(EOSBase)` in `laws.pyx` with a docstring and `MODEL_NAME = "foo"`, added to the `ModelFamily` list and exported from `TidalPy/Material/laws/__init__.py`.
4. Physics tests in `Tests/Test_Material/test_eos_laws_01.py` (the density against a closed form, and an inversion cross-check); `Tests/Test_Utilities/Test_Classes/test_spec_models_01.py` covers parameters, config, binary record, and errors unchanged. Then document the law here with its formula and references.

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
