# Material EOS Models (`Material_x.eos`)

_Updated: 2026-09-12_

A material equation-of-state model returns a mass density [kg m$^{-3}$]. The analytic models return it as a function of the local pressure [Pa]; the interpolated model returns it as a function of radius [m]. All four answer through the same call, `calc_density(pressure, temperature=0.0, radius=0.0)`, so the whole-planet solve does not need to know which kind it is holding.

The models follow the same pattern as the rheology, viscosity, cooling, and radiogenics hierarchies: an abstract base deriving from `PhysicsBase`, concrete subclasses, a name-based factory, and shared serialization. The analytic models are isothermal: the temperature argument exists for interface uniformity and is not currently used.

## Inheritance

```
c_TidalPyBaseClass
  └── c_PhysicsBase
        └── c_MaterialEOSBase          (abstract)
              ├── c_ConstantDensityEOS   aliases "constant", "uniform", "constant_density"
              ├── c_BirchMurnaghanEOS    aliases "bm", "birch_murnaghan", "birch-murnaghan"
              ├── c_VinetEOS             alias "vinet"
              └── c_InterpolatedEOS      aliases "interp", "interpolate", "interpolated"
```

The base declares `calc_density(pressure, temperature, radius)` pure virtual and adds four optional radius-varying getters, described below. The Cython classes mirror the hierarchy: `MaterialEOSBase`, `ConstantDensityEOS`, `BirchMurnaghanEOS`, `VinetEOS`, `InterpolatedEOS`.

## The four models

| Model (aliases) | Density from | Parameters |
|---|---|---|
| `constant` (`uniform`, `constant_density`) | nothing; incompressible | `reference_density` |
| `birch_murnaghan` (`bm`) | pressure | `reference_density`, `reference_bulk_modulus`, `bulk_modulus_derivative` |
| `vinet` | pressure | same as Birch-Murnaghan |
| `interpolated` (`interp`, `interpolate`) | radius, by table lookup | `radius`, `density`, and optional viscoelastic tables |

### Constant density

Returns the same density everywhere. An incompressible body is not a realistic planet, but it is the case with closed-form Love numbers, so it is what analytic checks and regression tests are written against. It is also the honest choice for a thin layer where the pressure range is too small for compression to matter.

### Birch-Murnaghan, third order

A finite-strain expansion around a reference state. With the compression $\eta = \rho / \rho_0 = V_0 / V$,

$$P(\eta) = \frac{3}{2} K_0 \left( \eta^{7/3} - \eta^{5/3} \right) \left[ 1 + \frac{3}{4} \left( K_0' - 4 \right) \left( \eta^{2/3} - 1 \right) \right]$$

where $K_0$ is the reference bulk modulus and $K_0'$ its pressure derivative. This is the standard equation of state of mineral physics, fitted to compression experiments across the mantle pressure range, and it is the default choice for a silicate or iron layer.

### Vinet

An alternative functional form, also called the universal equation of state, derived from a scaled interatomic potential rather than from a strain expansion. With $x = (V / V_0)^{1/3} = \eta^{-1/3}$,

$$P(x) = 3 K_0 \frac{1 - x}{x^2} \exp\left[ \frac{3}{2} \left( K_0' - 1 \right) \left( 1 - x \right) \right]$$

Vinet and Birch-Murnaghan agree closely at modest compression and diverge at high compression, where Vinet is generally the better extrapolation. Fitted $K_0$ and $K_0'$ values are specific to the form they were fitted with, so do not mix a Birch-Murnaghan fit into the Vinet model.

### Interpolated

Linear interpolation of a sorted radius-to-density table, clamped at both ends. This is the route for any profile computed elsewhere: run a full mineral-physics or thermal-evolution code, export the result as arrays, and load them here. It is also how the shipped PREM profile of the Earth is used.

The interpolated model optionally carries four more radius-varying tables alongside density: static shear modulus, static bulk modulus, shear viscosity, and bulk viscosity. When present they are read back with `calc_static_shear_modulus(radius)`, `calc_static_bulk_modulus(radius)`, `calc_shear_viscosity(radius)`, and `calc_bulk_viscosity(radius)`. These four getters exist on the base class, and the analytic models return NaN from all of them. During a whole-planet solve an interpolated layer's tabulated values take precedence over the layer's own constants, so a tabulated layer's moduli and viscosities vary with radius the way its density does. A world TOML that names a `data_file` gets these tables populated automatically; see the [TOML schema](../structures_x/config/toml_schema.md).

### The pressure inversion

The analytic laws give pressure as a function of compression, but the solve needs the inverse. Both compressible models invert their own law with a safeguarded Newton iteration that falls back to bisection, and they share one implementation.

The laws are monotonic in $\eta$ only over a finite range. The third-order Birch-Murnaghan correction term changes sign at large compression when $K_0' \ne 4$, so $P(\eta)$ turns over and can even go negative beyond that point. The inverter therefore brackets the root by expanding outward from $\eta = 1$, where the pressure is zero by construction, and stops at the turning point rather than assuming monotonicity across a fixed wide interval.

Two numerical knobs control the iteration, both carried in the config so they can be set per material:

| Setting | Default | Meaning |
|---|---|---|
| `invert_rtol` | `1e-13` | Relative convergence tolerance on the compression. |
| `invert_max_iters` | `60` | Hard iteration cap. A termination safeguard only; convergence normally takes well under ten steps. |

```python
bm = BirchMurnaghanEOS(3500.0, 1.3e11, 4.5, invert_rtol=1e-9, invert_max_iters=80)
bm.invert_rtol, bm.invert_max_iters      # (1e-09, 80)
```

Both appear in `get_config_dict()` and survive the binary round trip. Their defaults live in one place, the C++ `c_MaterialEOSConfig` member initializers, and the Python wrappers override them only when a value is supplied explicitly.

## Choosing a model

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
| `radius_m`, `density_kg_m3` | interpolated | `radius`, `density` |
| `shear_modulus_pa`, `bulk_modulus_pa`, `shear_viscosity_pas`, `bulk_viscosity_pas` | interpolated | `shear_modulus`, `bulk_modulus`, `shear_viscosity`, `bulk_viscosity` |

An unrecognized key is ignored without warning, so a misspelling produces a model built from defaults rather than an error.

### Attaching a model to a layer

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

- `get_config_dict()` returns the model name under the key `model` plus its parameters, with the interpolated tables as lists. The dict is accepted by `make_material_eos`, so a model round-trips through it.
- `save_config(path)` writes the same content as TOML.
- `save_binary(path)` and `load_binary(path, force=False)` use the TidalPy binary format, through the shared `c_PhysicsBase` helpers.

Binary class ids: 601 constant, 602 Birch-Murnaghan, 603 Vinet, 604 interpolated.

## C++ API

The models live in `material_eos_.hpp` (namespace `tidalpy`, header only). The C++ layer is the canonical one; the Cython classes above are wrappers over it, and every other C++ consumer, including layers attaching an EOS, the whole-planet solve, and binary reconstruction, uses these types directly.

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

One combined config shared by every model; each reads only the fields it needs. Its member initializers are the single source of default values for both C++ and Python.

| Field | Used by | Default |
|---|---|---|
| `reference_density` | all | `3500.0` |
| `reference_bulk_modulus` | Birch-Murnaghan, Vinet | `1.0e11` |
| `bulk_modulus_derivative` | Birch-Murnaghan, Vinet | `4.0` |
| `invert_rtol` | Birch-Murnaghan, Vinet | `d_EOS_INVERT_RTOL` (`1e-13`) |
| `invert_max_iters` | Birch-Murnaghan, Vinet | `d_EOS_INVERT_MAX_ITERS` (`60`) |
| `radius`, `density` and the four viscoelastic tables | interpolated | empty vectors |

### Classes and free functions

All models derive from `c_MaterialEOSBase : c_PhysicsBase` and override `calc_density(pressure, temperature, radius)`. Each has a default constructor and one taking the config. Accessors are `get_reference_density()` on all of them, plus `get_reference_bulk_modulus()`, `get_bulk_modulus_derivative()`, `get_invert_rtol()`, and `get_invert_max_iters()` on the two compressible models, and `get_num_points()` on the interpolated one.

| Function | Description |
|---|---|
| `eos_bm_pressure(eta, K0, K0_prime)` | Third-order Birch-Murnaghan pressure [Pa] at compression $\eta$. |
| `eos_vinet_pressure(eta, K0, K0_prime)` | Vinet pressure [Pa] at $\eta$. |
| `eos_invert_eta(pressure_target, K0, K0_prime, pressure_fn, rtol, max_iters)` | Inverts a monotonic pressure law for $\eta$. Shared by both compressible models; pass one of the two pressure functions. |
| `c_material_eos_model_from_name(name)` | Name or alias to enum, throwing `std::invalid_argument` on an unknown name. |
| `c_find_material_eos(model, config)` | Heap-allocates the model as a `unique_ptr`. A name-string overload is also provided. |
| `c_material_eos_from_binary(stream, force)` | Peeks the binary class id and reconstructs the matching model, used when a layer with an attached EOS is loaded. |

## Adding a new model

**C++ (`TidalPy/Material_x/eos/material_eos_.hpp`)**

1. Add any new parameters to `c_MaterialEOSConfig` with sensible defaults.
2. If the law is analytic, add a free pressure function monotonic in $\eta$ so the shared `eos_invert_eta` inverter can be reused. Otherwise compute the density directly.
3. Add the model class deriving from `c_MaterialEOSBase`: constructors, `get_*` accessors, the `calc_density` override, and `write_binary` / `read_binary` through the `c_PhysicsBase` helpers.
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
- Poirier, J.-P. (2000). *Introduction to the Physics of the Earth's Interior*, second edition. Comparison of the finite-strain and universal forms.
- Dziewonski, A. M., and Anderson, D. L. (1981). Preliminary reference Earth model. *Physics of the Earth and Planetary Interiors*, 25(4), 297-356.
