# Radiogenics (`radiogenics_x`)

`TidalPy.radiogenics_x` provides the C++ radiogenics model hierarchy that maps a
layer's mass and elapsed time to the **radiogenic heating** `Q` [W] produced by
the decay of radioactive isotopes. This is the single quantity the models compute
and expose (`calc_heating`); it feeds a layer's thermal budget.

Each `SolidLiquidLayer` can hold one radiogenics model.
`SolidLiquidLayer.calc_radiogenic_heating(time, mass)` delegates to the attached
model's `calc_heating` (returning `0` when no model is attached).

## Inheritance

```
c_TidalPyBaseClass
  └── c_PhysicsBase
        └── c_RadiogenicsBase        (abstract)
              ├── c_OffRadiogenics      alias "none"
              ├── c_IsotopeRadiogenics
              └── c_FixedRadiogenics    alias "constant"
```

The abstract base declares `calc_heating(time, mass)` as pure virtual; every
model returns the heating `Q` [W] directly.

## Radiogenic Models

All inputs are MKS: time and half-lives in seconds `[s]`, mass in `[kg]`, heat
production rates in `[W/kg]`, and heating in `[W]`. Each model decays from a
shared reference time so times can be expressed relative to any fixed epoch.

| Model | Heating `Q` [W] | Parameters |
|-------|-----------------|------------|
| `OffRadiogenics` (`none`) | `0` | — |
| `IsotopeRadiogenics` | `m · Σᵢ hprᵢ · fracᵢ · concᵢ · exp(γᵢ(t − t_ref))` | a list of isotopes, `ref_time_s` |
| `FixedRadiogenics` (`constant`) | `m · q · exp(γ(t − t_ref))` | `fixed_heat_production_w_kg`, `average_half_life_s`, `ref_time_s` |

where the decay constant for a half life `t½` is `γ = ln(0.5) / t½` (`d_LN_HALF`
is a module-level `constexpr`). The Fixed model treats `average_half_life_s ≤ 0`
as "no decay" (a constant heating rate). This math mirrors TidalPy's validated
legacy `radiogenics.radiogenic_models` functions.

### The `c_Isotope` value type

A single radioactive isotope is described by the lightweight `c_Isotope` struct
(C++) / per-isotope fields (Python). It is a plain value type (no base class, no
virtuals) carrying:

| Field | Meaning |
|-------|---------|
| `name` | isotope label (e.g. `"U238"`) |
| `heat_production_w_kg` | specific heat production of the pure isotope [W/kg] |
| `half_life_s` | half life [s] |
| `mass_frac` | isotopic mass fraction within its element [kg/kg] |
| `concentration` | element concentration in the layer material [kg/kg] |

It provides `decay_constant()` [1/s] and `specific_heating(time_s, ref_time_s)`
[W/kg]. `c_IsotopeRadiogenics` holds a `std::vector<c_Isotope>` and sums each
isotope's specific heating before scaling by the layer mass. At the Python level,
`IsotopeRadiogenics` accepts parallel arrays (`heat_production_w_kg`,
`half_lives_s`, `mass_fracs`, `concentrations`) plus optional `names`, and exposes
them back through the `heat_production`, `half_lives`, `mass_fracs`,
`concentrations`, and `isotope_names` properties.

## Usage

```python
from TidalPy.radiogenics_x import IsotopeRadiogenics, FixedRadiogenics, make_radiogenics

mass = 1.0e22       # kg
time = 1.0e17       # s (relative to the reference time)

# Fixed lumped rate, no decay.
f = FixedRadiogenics(fixed_heat_production_w_kg=1.0e-11)
Q = f.calc_heating(time, mass)   # heating [W]

# Explicit isotope set (all MKS).
iso = IsotopeRadiogenics(
    heat_production_w_kg=[9.48e-5, 2.69e-5],
    half_lives_s=[4.47e17, 1.40e18],
    mass_fracs=[0.9928, 0.9998],
    concentrations=[0.012e-6, 0.04e-6],
    ref_time_s=0.0,
    names=["U238", "Th232"])

# Built-in literature dataset (no hand-entered abundances needed).
chondritic = IsotopeRadiogenics.from_dataset("modern_day_chondritic")

# Name/alias-based factory (case-insensitive).
g = make_radiogenics("fixed", {"fixed_heat_production_w_kg": 2.0e-11})
```

`make_radiogenics(model_name, config=None)` resolves aliases case-insensitively
and forwards the relevant keys from `config`. Unknown names raise `ValueError`.

### Factory

At the C++ level the factory is enum-based. `c_RadiogenicsModel` names one value
per model, `c_radiogenics_model_from_name(name)` maps a (case-insensitive) name
or alias to that enum (throwing `std::invalid_argument` on an unknown name), and
`c_find_radiogenics(model, config)` returns a
`std::unique_ptr<c_RadiogenicsBase>` to a freshly heap-allocated model. This is
the canonical C++ factory used by all C++ consumers (layers attaching
radiogenics, future binary reconstruction). The Python `make_radiogenics` wraps
it: it builds the `c_RadiogenicsConfig`, calls `c_radiogenics_model_from_name`
then `c_find_radiogenics`, and adopts the returned `unique_ptr` into the matching
rich Python wrapper (`OffRadiogenics`, `IsotopeRadiogenics`, `FixedRadiogenics`).

### Built-in isotope datasets

The module ships a small catalog of literature-sourced isotope sets so a caller
can build a realistic model without entering abundances by hand. List them with
`available_isotope_datasets()`, inspect one (as an MKS dict) with
`isotope_dataset(name)`, and build a model with
`IsotopeRadiogenics.from_dataset(name)` or
`make_radiogenics("isotope", {"isotopes": name})`.

| Dataset | Isotopes | Applicability | Source |
|---------|----------|---------------|--------|
| `modern_day_chondritic` | U238, U235, Th232, K40 | Present-day rocky/icy bodies of broadly chondritic composition | Hussmann & Spohn (2004); Turcotte & Schubert (2001) |
| `llri_and_slri` | U238, U235, Th232, K40, Mn53, Fe60, Al26 | Early-solar-system thermal evolution (first few Myr), where short-lived isotopes dominate | Castillo-Rogez et al. (2007) |
| `bulk_silicate_earth` | U238, U235, Th232, K40 | Present-day Earth-like silicate mantles (U = 20.3 ppb, Th = 79.5 ppb, K = 240 ppm) | McDonough & Sun (1995) concentrations; Turcotte & Schubert (2002) rates |

These built-in datasets are defined in C++ (`c_get_isotope_dataset`) directly in
MKS (half lives and reference times converted from the literature's Myr using
`d_SECONDS_PER_MYR`); all use a reference time of 4600 Myr (concentrations quoted
at the present epoch). They take precedence in `make_radiogenics` over any
same-named entry in the global config.

`make_radiogenics("isotope", {"isotopes": "<name>"})` also falls back to a dataset
under `TidalPy.config['physics']['radiogenics']['known_isotope_data']` if the name
is not a built-in; those config datasets (and inline dataset dicts) store half
lives and reference times in **mega-years (Myr)** and are converted to seconds by
the Python factory. Alternatively, supply explicit MKS arrays
(`heat_production_w_kg`, `half_lives_s`, `mass_fracs`, `concentrations`,
`ref_time_s`, optional `isotope_names`) directly.

## Vectorized evaluation

Each model also provides vectorized heating methods (defined once on
`c_RadiogenicsBase`, so all models inherit them via virtual dispatch):

- `calc_heating_vectorize_time(time[], mass)` — a time sweep at constant mass.
- `calc_heating_vectorize_mass(time, mass[])` — a mass sweep at constant time.
- `calc_heating_vectorize_all(time[], mass[])` — element-wise over two
  equal-length arrays.

At the C++ level each fills a caller-supplied `std::vector<double>&` output;
mismatched input lengths throw `std::invalid_argument`. The Cython wrappers
accept array-likes and return a `float64` NumPy array.

## Direct convenience functions

For one-shot evaluation without explicitly constructing a model object, each
model has a lower-case module function:

```python
from TidalPy.radiogenics_x import fixed, isotope
import numpy as np

# Scalar in -> float out.
Q = fixed(time, mass, fixed_heat_production_w_kg=1.0e-11)

# Arrays in -> float64 ndarray out (time and mass may each be a float or an
# ndarray; they are broadcast together).
Q_decay = fixed(np.linspace(0.0, 1.0e18, 50), mass,
                fixed_heat_production_w_kg=1.0e-11, average_half_life_s=4.47e17)
Q_iso   = isotope(np.array([0.0, 5.0e17]), mass,
                  heat_production_w_kg=[9.48e-5], half_lives_s=[4.47e17],
                  mass_fracs=[0.9928], concentrations=[0.012e-6])
```

Signatures: `off(time, mass)`,
`isotope(time, mass, heat_production_w_kg, half_lives_s, mass_fracs, concentrations, ref_time_s=0.0, names=None)`,
`fixed(time, mass, fixed_heat_production_w_kg=0.0, average_half_life_s=0.0, ref_time_s=0.0)`.

Each builds a *stack-allocated* C++ model, solves (picking the most specific
vectorized routine for the input pattern), and returns. `time` and `mass` accept
floats or NumPy arrays; the model parameters (isotope arrays, fixed rate) are
always constants.

> **Note on argument order:** both the class methods (`calc_heating(time, mass)`)
> and the convenience functions (`func(time, mass, ...)`) take `time` first, then
> `mass`, then any model parameters. The order is consistent across the module.

## Serialization

Every model supports the standard TidalPy interfaces:

- `get_config_dict()` → dict of `model_name` plus any model parameters (isotope
  arrays are returned as lists).
- `save_config(path)` → TOML config file.
- `save_binary(path)` / `load_binary(path, force=False)` → TidalPy binary format
  (preserves the model name and parameters; the layer observer pointer is not
  serialized). The Isotope model writes its variable-length isotope list (each
  isotope's name plus its four doubles) after the shared header and model name (a
  justified exception to the scalar-only `write_physics_binary` helper used by the
  Off and Fixed models).

## Adding a new radiogenics model

The hierarchy is designed so a new model is a small, local addition. To add a
new radiogenics model named `Foo`:

**C++ (`TidalPy/radiogenics_x/radiogenics_.hpp`)**

1. If `Foo` needs new parameters, add them to `c_RadiogenicsConfig` (with
   sensible defaults). The single combined config is shared by all models.
2. Add a `rad_heating_foo(...)` free function implementing the heating law (use
   the `rad_guard` floor on any half-life denominator).
3. Add the model class `c_Foo : public c_RadiogenicsBase`:
   - constructors `c_Foo()` and `explicit c_Foo(const c_RadiogenicsConfig&)` that
     pass a model-name string to the base and copy any params into `p_*` members;
   - `get_*` accessors for each param;
   - override `calc_heating(time, mass)` returning the heating [W];
   - override `write_binary` / `read_binary` using the shared base helpers
     (`this->write_physics_binary(out, class_id, {params...})` and
     `this->read_physics_binary(in, force, n_params)` from `c_PhysicsBase`). If
     the model has variable-length data, write the header + model name + payload
     directly (see `c_IsotopeRadiogenics`).

**C++ (`TidalPy/Utilities_x/binary_x/binary_.hpp`)**

4. Add a unique `BinaryClassID::Foo` value (radiogenics models occupy 50X).

**C++ factory (`radiogenics_.hpp`)**

5. Add `Foo` to the `c_RadiogenicsModel` enum, map its name/aliases in
   `c_radiogenics_model_from_name`, and add a `case` to `c_find_radiogenics`. Also
   add a `case BinaryClassID::Foo` to `c_radiogenics_from_binary` (the
   binary-dispatch factory) so a `Foo` attached to a `SolidLiquidLayer` can be
   reconstructed when the layer is loaded from a binary file.

**Cython (`radiogenics.pxd` / `radiogenics.pyx`)**

6. Declare `c_Foo` (constructors + getters) in `radiogenics.pxd` and add the enum
   value to the `c_RadiogenicsModel` cimport.
7. Add the `cdef class Foo(RadiogenicsBase)` wrapper in `radiogenics.pyx` (with
   param properties and a `get_config_dict` override), the adoption branch in
   `make_radiogenics`, and the lower-case `foo(time, mass, ...)` convenience
   function.

**Package + tests + docs**

8. Export `Foo` (and `foo`) from `TidalPy/radiogenics_x/__init__.py` (and the C++
   names from `__init__.pxd`).
9. Add `Foo` to the test lists in
   `Tests/Test_Radiogenics_x/test_radiogenics_01.py` (model name, heating
   cross-check against an independent reference, factory/alias, vectorization,
   config dict, binary round-trip, isinstance).
10. Document the model here (formula, parameters, references).

No build-system change is needed — `radiogenics_x.radiogenics` is already
registered in `cython_extensions.json`.

## References

- Hussmann and Spohn (2004), Icarus — chondritic radiogenic isotope data.
- Turcotte and Schubert (2001, 2002), *Geodynamics* — radiogenic heat production rates.
- Castillo-Rogez et al. (2007), Icarus, [DOI: 10.1016/j.icarus.2007.02.018](https://doi.org/10.1016/j.icarus.2007.02.018) — long- and short-lived radiogenic isotopes.
- McDonough and Sun (1995), Chemical Geology, [DOI: 10.1016/0009-2541(94)00140-4](https://doi.org/10.1016/0009-2541(94)00140-4) — bulk silicate Earth elemental abundances.
