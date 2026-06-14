# Cooling (`cooling_x`)

`TidalPy.cooling_x` provides the C++ cooling (heat-transport) model hierarchy that
maps a layer's physical state to a **cooling result**: the surface heat flux `q`
[W/m²], the thermal boundary-layer thickness [m], and the Rayleigh and Nusselt
numbers. This is what the models compute and expose (`calc_cooling`); the heat
flux drives a layer's thermal evolution.

Each `SolidLiquidLayer` can hold one cooling model.

## Inheritance

```
c_TidalPyBaseClass
  └── c_PhysicsBase
        └── c_CoolingBase          (abstract)
              ├── c_OffCooling         alias "none"
              ├── c_ConvectiveCooling  alias "convective"
              └── c_ConductiveCooling  alias "conductive"
```

The abstract base declares `calc_cooling(inputs)` as pure virtual; every model
returns a `CoolingResult`.

## Inputs and result

A cooling evaluation depends on eight physical quantities, bundled (per the style
guide's >5-argument rule) in a `c_CoolingInputs` struct (the Python `calc_cooling`
takes them as explicit MKS arguments):

| Input | Meaning |
|-------|---------|
| `delta_temp_k` | temperature drop across the layer [K] |
| `thickness_m` | layer (or sub-layer) thickness [m] |
| `gravity_m_s2` | gravitational acceleration [m/s²] |
| `density_kg_m3` | bulk density [kg/m³] |
| `viscosity_pas` | dynamic viscosity [Pa·s] |
| `thermal_conductivity_w_mk` | thermal conductivity [W/m/K] |
| `thermal_diffusivity_m2_s` | thermal diffusivity [m²/s] |
| `thermal_expansion_1_k` | thermal expansivity [1/K] |

The result is a `CoolingResult` with `cooling_flux` [W/m²],
`boundary_layer_thickness` [m], `rayleigh`, and `nusselt`. Each field is a Python
`float` for scalar evaluations or a `float64` ndarray for vectorized ones.
`CoolingResult` supports `to_dict()` and iteration
(`flux, blt, ray, nu = result`).

## The three models

All quantities are MKS. The math mirrors TidalPy's validated legacy
`cooling.cooling_models` functions.

| Model | Cooling flux `q` [W/m²] | Boundary layer | Ra / Nu |
|-------|-------------------------|----------------|---------|
| `OffCooling` (`none`) | `0` | `0.5 · thickness` | `0` / `1` |
| `ConductiveCooling` (`conductive`) | `k · ΔT / thickness` | `thickness` | `0` / `1` |
| `ConvectiveCooling` (`convective`) | `k · ΔT / blt` | `thickness / Nu` | see below |

For convection the Rayleigh number is

```
Ra = expansion · density · gravity · ΔT · thickness³ / (viscosity · diffusivity)
Nu = max(alpha · (Ra / Ra_crit)^beta, 2)
blt = thickness / Nu
q   = k · ΔT / blt
```

with parameters `convection_alpha` (default `1.0`), `convection_beta`
(default `1/3`), and `critical_rayleigh` (default `1100.0`). The Nusselt number
has a floor of 2 (stagnant-lid limit). Degenerate inputs (`ΔT ≤ 0`, or a
thickness below the shared `minimum_layer_thickness` config floor) collapse to the
legacy edge behaviour (`Ra = 0`, `Nu = 2`).

## Usage

```python
from TidalPy.cooling_x import ConvectiveCooling, make_cooling

# delta_temp, thickness, gravity, density, viscosity, conductivity, diffusivity, expansion
cv = ConvectiveCooling(convection_alpha=1.0, convection_beta=1/3, critical_rayleigh=1100.0)
result = cv.calc_cooling(1000.0, 1.0e6, 9.8, 3300.0, 1.0e21, 4.0, 1.0e-6, 3.0e-5)
print(result.cooling_flux, result.nusselt)   # [W/m^2], [dimensionless]

# Name/alias-based factory (case-insensitive).
g = make_cooling("conductive")
```

`make_cooling(model_name, config=None)` resolves aliases case-insensitively and
forwards the convection parameters (`convection_alpha`, `convection_beta`,
`critical_rayleigh`) from `config`. Unknown names raise `ValueError`.

### Factory internals

At the C++ level the factory is enum-based. `c_CoolingModel` names one value per
model, `c_cooling_model_from_name(name)` maps a (case-insensitive) name or alias
to that enum (throwing `std::invalid_argument` on an unknown name), and
`c_find_cooling(model, config)` returns a `std::unique_ptr<c_CoolingBase>` to a
freshly heap-allocated model. The Python `make_cooling` wraps it and adopts the
returned `unique_ptr` into the matching rich Python wrapper.

## Vectorized evaluation

Each model also provides vectorized methods (defined once on `c_CoolingBase`, so
all models inherit them via virtual dispatch). The two "live" inputs that change
during thermal evolution are the temperature drop and the viscosity; the
remaining inputs are held constant:

- `calc_cooling_vectorize_temperature(delta_temp[], <other 7 scalars>)` — a
  temperature-drop sweep.
- `calc_cooling_vectorize_viscosity(<delta_temp>, ..., viscosity[], ...)` — a
  viscosity sweep.
- `calc_cooling_vectorize_all(delta_temp[], ..., viscosity[], ...)` — element-wise
  over two equal-length arrays.

At the C++ level each fills a caller-supplied `std::vector<c_CoolingResult>&`
(copying the base `c_CoolingInputs` and overriding the swept field per element);
mismatched input lengths throw `std::invalid_argument`. The Cython wrappers return
a `CoolingResult` whose fields are `float64` NumPy arrays.

## Direct convenience functions

For one-shot evaluation without explicitly constructing a model object:

```python
from TidalPy.cooling_x import convective, conductive, cooling_off
import numpy as np

# Scalar in -> CoolingResult of floats.
r = convective(1000.0, 1.0e21, 1.0e6, 9.8, 3300.0, 4.0, 1.0e-6, 3.0e-5)

# delta_temp and viscosity may be floats or arrays (broadcast together);
# the remaining inputs are scalar constants.
sweep = convective(np.linspace(100.0, 2000.0, 50), 1.0e21,
                   1.0e6, 9.8, 3300.0, 4.0, 1.0e-6, 3.0e-5)
```

Signatures:
`cooling_off(delta_temp_k, thickness_m)`,
`conductive(delta_temp_k, thickness_m, thermal_conductivity_w_mk)`,
`convective(delta_temp_k, viscosity_pas, thickness_m, gravity_m_s2, density_kg_m3, thermal_conductivity_w_mk, thermal_diffusivity_m2_s, thermal_expansion_1_k, convection_alpha=1.0, convection_beta=1/3, critical_rayleigh=1100.0)`.

Each builds a *stack-allocated* C++ model, solves (picking the most specific
vectorized routine for the input pattern), and returns a `CoolingResult`. The off
and conduction functions only take the inputs they actually use.

> **Note on argument order:** the class method `calc_cooling` takes all eight
> physical inputs in a fixed order (`delta_temp, thickness, gravity, density,
> viscosity, conductivity, diffusivity, expansion`); the convenience functions
> lead with the live inputs (`delta_temp`, then `viscosity` for convection).

## Serialization

Every model supports the standard TidalPy interfaces:

- `get_config_dict()` → dict of `model_name` plus any model parameters.
- `save_config(path)` → TOML config file.
- `save_binary(path)` / `load_binary(path, force=False)` → TidalPy binary format
  (binary class IDs 401–403; Off and Conduction write no parameters, Convection
  writes its three scalars). The layer observer pointer is not serialized.

## Adding a new cooling model

The hierarchy is designed so a new model is a small, local addition. To add a new
cooling model named `Foo`:

**C++ (`TidalPy/cooling_x/cooling_.hpp`)**

1. If `Foo` needs new parameters, add them to `c_CoolingConfig` (with defaults).
2. Add a `cool_foo(const c_CoolingInputs&[, const c_CoolingConfig&])` free function
   implementing the heat-transport law (use the `cool_guard` floor on any
   denominator).
3. Add the model class `c_Foo : public c_CoolingBase`:
   - constructors `c_Foo()` and `explicit c_Foo(const c_CoolingConfig&)` that pass
     a model-name string to the base and copy any params into `p_*` members;
   - `get_*` accessors for each param;
   - override `calc_cooling(const c_CoolingInputs&)` returning a `c_CoolingResult`;
   - override `write_binary` / `read_binary` using the shared base helpers.

**C++ (`TidalPy/Utilities_x/binary_x/binary_.hpp`)**

4. Add a unique `BinaryClassID::Foo` value (cooling models occupy 40X).

**C++ factory (`cooling_.hpp`)**

5. Add `Foo` to the `c_CoolingModel` enum, map its name/aliases in
   `c_cooling_model_from_name`, and add a `case` to `c_find_cooling`.

**Cython (`cooling.pxd` / `cooling.pyx`)**

6. Declare `c_Foo` (constructors + getters) in `cooling.pxd` and add the enum
   value to the `c_CoolingModel` cimport.
7. Add the `cdef class Foo(CoolingBase)` wrapper in `cooling.pyx` (with param
   properties and a `get_config_dict` override), the adoption branch in
   `make_cooling`, and a lower-case `foo(...)` convenience function.

**Package + tests + docs**

8. Export `Foo` (and `foo`) from `TidalPy/cooling_x/__init__.py` (and the C++
   names from `__init__.pxd`).
9. Add `Foo` to the test lists in `Tests/Test_Cooling_x/test_cooling_01.py`
   (model name, cooling cross-check against an independent reference,
   factory/alias, vectorization, config dict, binary round-trip, isinstance).
10. Document the model here (formula, parameters, references).

No build-system change is needed — `cooling_x.cooling` is already registered in
`cython_extensions.json`.

## References

- Turcotte and Schubert (2002), *Geodynamics* — Rayleigh/Nusselt convection scaling and conduction.
- Solomatov (1995), Physics of Fluids — stagnant-lid convection scaling.
- Schubert, Turcotte, and Olson (2001), *Mantle Convection in the Earth and Planets* — boundary-layer theory.
