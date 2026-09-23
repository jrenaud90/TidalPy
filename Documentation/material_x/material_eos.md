# Material EOS Models (`Material_x.eos`)

_Updated: 2026-09-23_

A material equation-of-state model returns a mass density [kg m$^{-3}$]. The analytic models return it as a function of the local pressure [Pa]; the interpolated model returns it as a function of radius [m]. All four are evaluated through the same call, `calc_density(pressure, temperature=None, radius=0.0)`, so the whole-planet solve does not need to know which kind it is holding.

All of the models are built on an abstract base class deriving from `PhysicsBase`. Every model also takes a thermal expansivity [K$^{-1}$] and a reference temperature [K], so the density can depend on temperature; with the default expansivity of zero a model is athermal.

An EOS model is the layer's **material**. It owns every frequency-independent property of the layer (the static moduli, the shear law, the static viscosities, and the optional viscosity and partial-melt models), and it is the only thing that calculates them; see [The Material](#the-material).

## Inheritance

```
c_TidalPyBaseClass
  └── c_PhysicsBase
        └── c_MaterialEOSBase  (abstract)
              ├── c_ConstantDensityEOS  aliases "constant", "uniform", "constant_density"
              ├── c_BirchMurnaghanEOS   aliases "bm", "birch_murnaghan", "birch-murnaghan"
              ├── c_VinetEOS            alias "vinet"
              └── c_InterpolatedEOS     aliases "interp", "interpolate", "interpolated"
```

The base declares `calc_density(pressure, temperature, radius)` pure virtual, holds the two thermal parameters, provides `calc_bulk_modulus(pressure, temperature, radius)`, and adds four optional radius-varying getters, described below. The Cython classes mirror the hierarchy: `MaterialEOSBase`, `ConstantDensityEOS`, `BirchMurnaghanEOS`, `VinetEOS`, `InterpolatedEOS`.

## Models

| Model (aliases) | Density from | Parameters |
|---|---|---|
| `constant` (`uniform`, `constant_density`) | nothing; incompressible | `reference_density` |
| `birch_murnaghan` (`bm`) | pressure | `reference_density`, `reference_bulk_modulus`, `bulk_modulus_derivative` |
| `vinet` | pressure | same as Birch-Murnaghan |
| `interpolated` (`interp`, `interpolate`) | radius, by table lookup | `radius`, `density`, and optional viscoelastic tables |

Every model also accepts `thermal_expansion` and `reference_temperature`; see Thermal Terms below.

### Constant Density

Returns the same density everywhere. An incompressible body is not realistic but can be a useful diagnostic or applicable to small moons. It is also used to check and test the analytic models.

### Birch-Murnaghan, Third Order

A finite-strain expansion around a reference state. With the compression $\eta = \rho / \rho_0 = V_0 / V$,

$$P(\eta) = \frac{3}{2} K_0 \left( \eta^{7/3} - \eta^{5/3} \right) \left[ 1 + \frac{3}{4} \left( K_0' - 4 \right) \left( \eta^{2/3} - 1 \right) \right]$$

where $K_0$ is the reference bulk modulus and $K_0'$ its pressure derivative. This is the standard equation of state of mineral physics, fitted to compression experiments across the mantle pressure range, and it is the default choice for a silicate or iron layer.

### Vinet

An alternative functional form, also called the universal equation of state, derived from a scaled interatomic potential rather than from a strain expansion. With $x = (V / V_0)^{1/3} = \eta^{-1/3}$,

$$P(x) = 3 K_0 \frac{1 - x}{x^2} \exp\left[ \frac{3}{2} \left( K_0' - 1 \right) \left( 1 - x \right) \right]$$

Vinet and Birch-Murnaghan agree closely at modest compression and diverge at high compression, where Vinet is generally the better extrapolation. Fitted $K_0$ and $K_0'$ values are specific to the form they were fitted with, so do not mix a Birch-Murnaghan fit into the Vinet model.

### Interpolated

Linear interpolation of a sorted radius-to-density table, clamped at both ends. This is the route for any profile computed elsewhere: run a full mineral-physics or thermal-evolution code, export the result as arrays, and load them here. For example, TidalPy uses a PREM profile of the Earth via this model.

The interpolated model optionally carries four more radius-varying tables alongside density: static shear modulus, static bulk modulus, shear viscosity, and bulk viscosity. When present they are read back with `get_tabulated_shear_modulus(radius)`, `get_tabulated_bulk_modulus(radius)`, `get_tabulated_shear_viscosity(radius)`, and `get_tabulated_bulk_viscosity(radius)`. These four lookups exist on the base class, and the analytic models return NaN from all of them. A tabulated value takes precedence over the material's law or constant (see [The Material](#the-material)), so a tabulated layer's moduli and viscosities vary with radius the way its density does. A world TOML that names a `data_file` gets these tables populated automatically (using a PREM-like data file format); see the [TOML schema](../structures_x/config/toml_schema.md).

### Thermal Terms

Each model holds a thermal expansivity $\alpha_0$ [K$^{-1}$] and a reference temperature $T_\mathrm{ref}$ [K], the temperature at which $\rho_0$ and $K_0$ apply (default 300 K, where mineral-physics parameters are usually quoted).

Birch-Murnaghan and Vinet add a thermal pressure to their cold pressure law:

$$P(\eta, T) = P_\mathrm{cold}(\eta) + \alpha_0 K_0 \left( T - T_\mathrm{ref} \right)$$

The product $\alpha K_T$ is taken as constant, which is its high-temperature limit (Anderson 1995). The density then comes from inverting the cold law at $P - \alpha_0 K_0 (T - T_\mathrm{ref})$.

The constant and interpolated models have no pressure law to carry a thermal pressure, so they scale their density by the thermal expansion directly:

$$\rho(T) = \rho \, \exp\left[ -\alpha_0 \left( T - T_\mathrm{ref} \right) \right]$$

A model is athermal when $\alpha_0 = 0$ (the default) or when no temperature is passed (`None` in Python, a non-finite value in C++).

### Bulk Modulus

`calc_bulk_modulus(pressure, temperature=None, radius=0.0)` returns the isothermal bulk modulus $K_T = \rho \, \partial P / \partial \rho$ [Pa]. Birch-Murnaghan and Vinet evaluate the analytic derivative of their pressure law at the solved compression, so the modulus is consistent with the density and equals $K_0$ at the reference state. The interpolated model returns its bulk table at `radius`, and a model with neither returns NaN, which tells a layer to use its own constant.

### The Material

The EOS model carries the layer's material information:

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

The thermal expansivity is the `thermal_expansion` every model already takes. There is one $\alpha$ per material: it sets the adiabatic gradient $\alpha T g / c_p$ and the Rayleigh number of a convecting layer, and the density law uses the same number, but only when it is handed a temperature, which a layer controls with `use_thermal_eos`. `calc_thermal_diffusivity(density)` returns $\kappa = k / (\rho c_p)$.

The material also holds three optional models, attached with `set_shear_viscosity(model)`, `set_bulk_viscosity(model)` (from [`viscosity_x`](../viscosity_x/viscosity_models.md)) and `set_partial_melt(model)` (from [`partial_melt_x`](../partial_melt_x/partial_melt_models.md)). A layer has the same three methods as helpers; they hand the model to its EOS.

`calc_material_state(pressure, temperature=None, radius=0.0, thermal_density=True)` maps a point onto all of it, in this order:

1. The static shear modulus from the linear law

   $$\mu = \mu_0 + \mu'_P P + \mu'_T \left( T - T_\mathrm{ref} \right)$$

   floored at the config's `minimum_modulus`, and the bulk modulus from `bulk_modulus_static`.
2. The viscosities: the attached viscosity model at the temperature and pressure, else the static viscosity.
3. The density and, where the model defines one, the bulk modulus of the density law. The law sees the temperature only when `thermal_density` is set, which is what a layer's `use_thermal_eos` switch controls; the viscosity and partial-melt models always see it.
4. A table of an interpolated model replaces the law or constant, and the bulk modulus of a Birch-Murnaghan or Vinet law replaces `bulk_modulus_static`.
5. The partial-melt model, applied to the shear modulus and viscosity and then to the bulk pair.

It returns a dict of `density`, `melt_fraction`, `shear_modulus`, `bulk_modulus`, `shear_viscosity`, and `bulk_viscosity`, all after the partial-melt step.

This is what the whole-planet EOS solve evaluates as it integrates, so there is one path to these numbers. **After a solve, do not call it to find the state of the planet**: read the world and layer getters (`get_shear_modulus(radius)`, `get_melt_fraction(radius)`, `get_state(radius)`, ...), which report what the solve used. `calc_material_state` is for asking the material about a pressure and temperature of your choosing.

```python
from TidalPy.Material_x.eos import BirchMurnaghanEOS
from TidalPy.viscosity_x import make_viscosity
from TidalPy.partial_melt_x import make_partial_melt

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

The analytic laws give pressure as a function of compression, but the solve needs the inverse. Both compressible models invert their own law with Newton's method and share one implementation. The slope is exact, because the bulk modulus of each law is already $K = \eta \, dP/d\eta$, and one evaluation returns the pressure and $K$ together. The first guess is the Murnaghan law, $\eta = (1 + K_0' P / K_0)^{1/K_0'}$, which inverts in closed form and follows both laws closely over planetary compressions, so convergence to `invert_rtol` usually takes four or five evaluations. Every evaluation tightens a bracket around the root, and a step that leaves the bracket is replaced by its midpoint.

The laws are monotonic in $\eta$ only over a finite range. Every law turns over in tension, and the third-order Birch-Murnaghan correction term changes sign at large compression when $K_0' < 4$, so $P(\eta)$ turns over there too. This range depends only on $K_0$ and $K_0'$, so each model finds it once, when it is built or loaded, by stepping outward from $\eta = 1$ until $K$ stops being positive and bisecting that sign change. A pressure outside the range has no compression to find and returns the compression at that end of the range. The density is therefore continuous in pressure everywhere, which the structure solve relies on: while its central pressure is still a guess, its outer radii can sit far into tension.

Two numerical knobs control the iteration. Their defaults come from `[numerical]` in `TidalPy_Configs_x.toml` (`eos_invert_rtol`, `eos_invert_max_iters`), and each material can set its own:

| Setting | Default | Meaning |
|---|---|---|
| `invert_rtol` | `1e-13` | Relative convergence tolerance on the compression. |
| `invert_max_iters` | `60` | Hard iteration cap. A termination safeguard only; convergence normally takes well under ten steps. |

```python
bm = BirchMurnaghanEOS(3500.0, 1.3e11, 4.5, invert_rtol=1e-9, invert_max_iters=80)
bm.invert_rtol, bm.invert_max_iters      # (1e-09, 80)
```

A model built without its own values takes the configured defaults when it is built, so a later change to the config does not reach an existing model. Both appear in `get_config_dict()` and survive the binary round trip, so a saved world records the values it was solved with.

### Choosing a Model

Use `constant` for an analytic check, for a regression test with a known closed-form answer, or for a layer thin enough that compression is negligible.

Use `birch_murnaghan` when reproducing published mineral-physics parameters, which are most often quoted in this form.

Use `vinet` when the fit was made in that form, or when the layer reaches compressions where the two forms visibly disagree.

Use `interpolated` whenever a profile already exists, whether from a seismic reference model, a mineral-physics package, or a previous TidalPy run. It is also the only model that can vary the moduli and viscosities with radius inside a single layer.

## Python API

```python
from TidalPy.Material_x.eos import (
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

`make_material_eos(model_name, config=None)` resolves the name or alias case-insensitively and raises `ValueError` for an unrecognized one. Its config keys are the ones `get_config_dict()` emits and carry unit suffixes, unlike the constructor arguments.

| Config key | Model | Constructor argument |
|---|---|---|
| `reference_density_kg_m3` | all analytic | `reference_density` |
| `reference_bulk_modulus_pa` | Birch-Murnaghan, Vinet | `reference_bulk_modulus` |
| `bulk_modulus_derivative` | Birch-Murnaghan, Vinet | `bulk_modulus_derivative` |
| `invert_rtol`, `invert_max_iters` | Birch-Murnaghan, Vinet | same |
| `thermal_expansion_1_k`, `reference_temperature_k` | all | `thermal_expansion`, `reference_temperature` |
| `radius_m`, `density_kg_m3` | interpolated | `radius`, `density` |
| `shear_modulus_pa`, `bulk_modulus_pa`, `shear_viscosity_pas`, `bulk_viscosity_pas` | interpolated | `shear_modulus`, `bulk_modulus`, `shear_viscosity`, `bulk_viscosity` |

A key that no material EOS model reads raises `ValueError` naming the closest accepted key, so a misspelling or a missing unit suffix fails loudly instead of silently building a default model.

### Attaching a Model to a `Layer`

```python
from TidalPy.Material_x.eos import make_material_eos

core.set_eos(make_material_eos("constant", {"reference_density_kg_m3": 11000.0}))
mantle.set_eos(make_material_eos("bm", {"reference_density_kg_m3":   3500.0,
                                        "reference_bulk_modulus_pa": 1.3e11,
                                        "bulk_modulus_derivative":   4.5}))
world.solve_eos(surface_pressure=0.0)
```

`set_eos` moves ownership of the C++ model into the layer, leaving the Python wrapper an empty shell. Every layer class accepts one, and `solve_eos` raises `ValueError` if any layer is missing it. See [Worlds](../structures_x/worlds/worlds.md) for the solve itself and its results.

## Serialization

Every model supports the standard interfaces inherited from the base class.

- `get_config_dict()` returns the model name under the key `model` plus its parameters, with the interpolated tables as lists, the material keys (`shear_modulus_static_pa`, `bulk_modulus_static_pa`, the two `*_viscosity_static_pas` when set, and the three shear-law keys), and a sub-table for each attached model (`shear_viscosity`, `bulk_viscosity`, `partial_melt`). The dict is accepted by `make_material_eos`, which builds and attaches the nested models, so a material round-trips through it. This is the `[layers.<name>.material]` table of a world TOML.
- `save_config(path)` writes the same content as TOML.
- A model attached to a layer is saved and restored with that layer: the layer's binary record and its `get_config_dict()` (under `material`) both carry it, so a world or a layer round trip needs no separate handling of its materials.
- `save_binary(path)` and `load_binary(path, force=False)` use the TidalPy binary format, through the shared `c_PhysicsBase` helpers. Each model's record is followed by the material section: nine doubles (the four static constants, the three shear-law parameters, the conductivity, and the heat capacity), then a presence flag and nested record for each of the three optional models.

Binary class ids: 601 constant, 602 Birch-Murnaghan, 603 Vinet, 604 interpolated.

## C++ API

The models live in `material_eos_.hpp` (namespace `tidalpy`, header only). The C++ layer is canonical: the Cython classes above are wrappers over it, and every C++ consumer (layers attaching an EOS, the whole-planet solve, and binary reconstruction) uses these types directly.

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

One combined config shared by every model; each reads only the fields it needs. Its member initializers are the single source of default values for both C++ and Python, except the two inversion settings, which are unset (NaN and -1) and take the `[numerical]` defaults when a model is built.

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
| `radius`, `density` and the four viscoelastic tables | interpolated | empty vectors |

### Classes and Free Functions

All models derive from `c_MaterialEOSBase : c_PhysicsBase` and override `calc_density(pressure, temperature, radius)`. The base also has `calc_bulk_modulus(pressure, temperature, radius)` and the virtual `calc_density_and_bulk_modulus(pressure, temperature, radius, density, bulk_modulus)`, which returns both from one pressure inversion. Each model has a default constructor and one taking the config. Accessors are `get_thermal_expansion()`, `get_reference_temperature()`, and `get_reference_density()` on all of them, plus `get_reference_bulk_modulus()`, `get_bulk_modulus_derivative()`, `get_invert_rtol()`, and `get_invert_max_iters()` on the two compressible models, and `get_num_points()` on the interpolated one.

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

**C++ (`TidalPy/Material_x/eos/material_eos_.hpp`)**

1. Add any new parameters to `c_MaterialEOSConfig` with sensible defaults.
2. If the law is analytic, add a free function that fills the pressure and the bulk modulus $\eta \, dP/d\eta$ at a compression, so the shared `eos_find_monotonic_range` and `eos_invert_eta` can be reused, and give the model an `update_law_range` called from its constructors and from `read_binary`. Otherwise compute the density directly.
3. Add the model class deriving from `c_MaterialEOSBase`: constructors (pass the config to the base so the thermal parameters are stored), `get_*` accessors, the `calc_density` override, a `calc_density_and_bulk_modulus` override if the law defines a bulk modulus, and `write_binary` / `read_binary` through the `c_PhysicsBase` helpers, including the two thermal parameters.
4. Add the enum value, the name and alias branch in `c_material_eos_model_from_name`, and the cases in `c_find_material_eos` and `c_material_eos_from_binary`.

**C++ (`TidalPy/Utilities_x/binary_x/binary_.hpp`)**

5. Add a unique `BinaryClassID` in the 60X block.

**Cython (`material_eos.pxd` and `material_eos.pyx`)**

6. Declare the C++ class and the new enum value in the `.pxd`.
7. Add the `cdef class` wrapper with parameter properties and the adoption branch in `make_material_eos`. The config dict comes from the C++ `append_config_entries` override.

**Package, tests, and docs**

8. Export the class from `__init__.py`.
9. Extend `Tests/Test_Material_x/test_material_eos_01.py`: density evaluation, an inversion cross-check for an analytic law, factory and aliases, config dict, and binary round trip.
10. Document the model here with its formula and references, and add a changelog entry.

No build-system change is needed; `Material_x.eos.material_eos` is already registered in `cython_extensions.json`.

## References

- Birch, F. (1947). Finite elastic strain of cubic crystals. *Physical Review*, 71(11), 809-824.
- Vinet, P., Ferrante, J., Rose, J. H., and Smith, J. R. (1987). Compressibility of solids. *Journal of Geophysical Research*, 92(B9), 9319-9325.
- Anderson, O. L. (1995). *Equations of State of Solids for Geophysics and Ceramic Science*. Oxford University Press. Thermal pressure and the near constancy of $\alpha K_T$ at high temperature.
- Poirier, J.-P. (2000). *Introduction to the Physics of the Earth's Interior*, second edition. Comparison of the finite-strain and universal forms.
- Dziewonski, A. M., and Anderson, D. L. (1981). Preliminary reference Earth model. *Physics of the Earth and Planetary Interiors*, 25(4), 297-356.
