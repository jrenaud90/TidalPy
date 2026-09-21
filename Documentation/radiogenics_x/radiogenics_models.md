# Radiogenic Models (`radiogenics_x`)

_Updated: 2026-09-20_

A radiogenics model utilizes a layer of mass $m$ at time $t$ to find how much power is being released inside it by radioactive decay. The heating $Q$ \[W\] is returned by `calc_heating(time, mass)`.

Time is measured in seconds from an epoch the caller chooses, and each model carries the reference time `ref_time` at which its rates or concentrations were quoted. Only the difference $t - t_{\text{ref}}$ enters the physics, so a model built from present-day abundances with a reference time of 4600 Myr is evaluated at $t = 0$ to get the heating at the birth of the solar system, or at $t = t_{\text{ref}}$ to get today's.

## Inheritance

```
c_TidalPyBaseClass
  └── c_PhysicsBase
        └── c_RadiogenicsBase  (abstract)
              ├── c_OffRadiogenics       alias "none"
              ├── c_IsotopeRadiogenics
              └── c_FixedRadiogenics     alias "constant"
```

The abstract base declares `calc_heating(time, mass)` and supplies the three vectorized wrappers, so a new model only has to implement the heating law.

## Models

For a half life $t_{1/2}$ the decay constant is $\gamma = \ln(0.5) / t_{1/2}$, a negative number whose magnitude grows as the half life shortens.

| Model (aliases) | Heating $Q$ [W] | Parameters |
|---|---|---|
| `off` (`none`) | $0$ | none |
| `isotope` (`isotopes`) | $m \sum_i q_i f_i c_i \exp[\gamma_i (t - t_{\text{ref}})]$ | a list of isotopes, `ref_time` |
| `fixed` (`constant`) | $m\, q \exp[\gamma (t - t_{\text{ref}})]$ | `fixed_heat_production`, `average_half_life`, `ref_time` |

Here $q_i$ is the specific heat production of pure isotope $i$, $f_i$ its mass fraction within its parent element, and $c_i$ that element's concentration in the layer material. The product $q_i f_i c_i$ is the specific heating the isotope contributed at the reference time, per kilogram of layer, and the exponential carries it forward or backward in time.

### Off

Returns zero for any time and mass.

### Isotope

Sums the decay of an arbitrary list of isotopes, each with its own half life.

### Fixed

Applies one lumped specific rate to the whole layer, optionally with a single effective half life. Setting `average_half_life` to zero (the default) or to any non-positive value means no decay at all, and the heating is then constant for all time.

### Behavior at the Limits

A half life at or below zero is treated as infinite rather than as an error, which is what makes the constant-rate case fall out of the same formula. A half life that is finite but smaller than the module's floor is clamped to that floor, so no decay constant is ever divided by zero.

Evaluating a model far before its reference time asks for an exponential that would overflow. Both decaying models guard against this and return NaN, so a bad epoch shows up as NaN heating.

## Isotope Value Type

One isotope is described by the `c_Isotope` C++ struct, a plain value type with no base class and no virtual functions.

| Field | Meaning |
|---|---|
| `name` | Isotope label, for example `"U238"`. |
| `heat_production` | Specific heat production of the pure isotope [W kg$^{-1}$]. |
| `half_life` | Half life [s]. |
| `mass_frac` | Mass fraction of the isotope within its element [kg kg$^{-1}$]. |
| `concentration` | Concentration of the parent element in the layer material [kg kg$^{-1}$]. |

It provides `decay_constant()` [s$^{-1}$] and `specific_heating(time, ref_time)` [W kg$^{-1}$]. `c_IsotopeRadiogenics` holds a `std::vector<c_Isotope>` and sums the specific heating of each member before scaling by the layer mass.

`IsotopeRadiogenics` is the Python wrapper to this struct. It takes parallel arrays and reads them back through properties of the same names.

```python
from TidalPy.radiogenics_x import IsotopeRadiogenics

model = IsotopeRadiogenics(
    heat_production=[9.48e-5, 2.69e-5],   # [W/kg] of the pure isotope
    half_lives=[4.47e17, 1.40e18],        # [s]
    mass_fracs=[0.9928, 0.9998],          # [kg/kg] within the element
    concentrations=[0.012e-6, 0.04e-6],   # [kg/kg] of the element in the material
    ref_time=0.0,                         # [s]
    names=["U238", "Th232"])

model.num_isotopes      # 2
model.isotope_names     # ['U238', 'Th232']
model.heat_production   # ndarray, and likewise half_lives, mass_fracs, concentrations
```

The four numeric arrays must be the same length. Labels are optional and are auto-generated as `isotope_0`, `isotope_1`, and so on when omitted.

## Built-in Isotope Datasets

TidalPy provides some sets of isotopes popular in the literature. List them with `available_isotope_datasets()`, inspect one with `isotope_dataset(name)`, and build a model from one with `IsotopeRadiogenics.from_dataset(name)` or `make_radiogenics("isotope", {"isotopes": name})`.

| Dataset | Isotopes | Reference time | Applicability | Source |
|---|---|---|---|---|
| `modern_day_chondritic` | U238, U235, Th232, K40 | 4600 Myr | Present-day rocky and icy bodies of broadly chondritic composition. | Hussmann and Spohn (2004); Turcotte and Schubert (2001) |
| `llri_and_slri` | U238, U235, Th232, K40, Mn53, Fe60, Al26 | 0 Myr | Early solar system thermal evolution, where the short-lived isotopes dominate the first few million years. | Castillo-Rogez et al. (2007) |
| `bulk_silicate_earth` | U238, U235, Th232, K40 | 4600 Myr | Present-day Earth-like silicate mantles (U 20.3 ppb, Th 79.5 ppb, K 240 ppm). | McDonough and Sun (1995) concentrations; Turcotte and Schubert (2002) rates |

> [!NOTE]
> Notice that the reference times differ. The two present-day sets quote concentrations at 4600 Myr, so evaluating them at $t = 0$ gives the heating at solar system formation and evaluating at $t = $ `ref_time` gives today's. The short-lived set quotes its concentrations at formation instead, so it is evaluated at the body's age directly.

`isotope_dataset(name)` returns the dataset as a dict with the keys `heat_production_w_kg`, `half_lives_s`, `mass_fracs`, `concentrations`, `isotope_names`, and `ref_time_s`.

A name that is not one of the three built-ins is looked up in the global config under `TidalPy.config['physics']['radiogenics']['known_isotope_data']`, and an inline dict is accepted in the same place. Those two sources follow the older convention of storing half lives and reference times in mega-years, with the per-isotope keys `hpr`, `half_life`, `iso_mass_fraction`, and `element_concentration`; the Python factory converts them to seconds. The built-in catalog is always preferred over a same-named config entry.

## Python API

```python
from TidalPy.radiogenics_x import (
    OffRadiogenics, IsotopeRadiogenics, FixedRadiogenics,
    make_radiogenics, available_isotope_datasets, isotope_dataset)

mass = 1.0e22   # [kg]
time = 1.0e17   # [s] from the caller's epoch

# A lumped rate with no decay.
model = FixedRadiogenics(fixed_heat_production=1.0e-11)
heating = model.calc_heating(time, mass)                    # [W]

# A literature isotope set, evaluated at its own reference epoch (present day).
chondritic = IsotopeRadiogenics.from_dataset("modern_day_chondritic")
heating_today = chondritic.calc_heating(chondritic.ref_time, mass)

# Name factory: case-insensitive, aliases accepted.
model = make_radiogenics("constant", {"fixed_heat_production_w_kg": 2.0e-11})
model = make_radiogenics("isotope", {"isotopes": "bulk_silicate_earth"})
```

`make_radiogenics(model_name, config=None)` resolves the name or alias case-insensitively and raises `ValueError` for an unrecognized one. Its config keys carry unit suffixes, matching what `get_config_dict()` emits, and are not the same spellings as the constructor arguments.

| Config key | Model | Meaning |
|---|---|---|
| `fixed_heat_production_w_kg` | fixed | Lumped specific rate [W kg$^{-1}$]. |
| `average_half_life_s` | fixed | Effective half life [s]; non-positive means no decay. |
| `ref_time_s` | fixed, isotope | Reference time [s]. |
| `isotopes` | isotope | A built-in dataset name, a config dataset name, or an inline dict. |
| `heat_production_w_kg`, `half_lives_s`, `mass_fracs`, `concentrations`, `isotope_names` | isotope | Explicit parallel arrays, in MKS. |

### Attaching a Model to a `Layer`

```python
layer.set_radiogenics(IsotopeRadiogenics.from_dataset("modern_day_chondritic"))
layer.radiogenics_set                         # True
layer.calc_radiogenic_heating(time, mass)     # [W] for the mass supplied
world.calc_internal_heating(time)             # [W] summed over all layers
```

`set_radiogenics` moves ownership of the C++ model into the layer, leaving the Python wrapper an empty shell, so build a fresh model if the same parameters are needed elsewhere. Only `SolidLiquidLayer` accepts one. A layer without a model reports zero heating rather than raising, and a world sums only the layers that carry one. A layer with `use_heating` set also feeds its model to the world's thermal EOS solve, which heats the layer at the model's specific rate times the local density and reports the total as `layer_heating` (see [Worlds](../structures_x/worlds/worlds.md)).

The mass is an argument so the caller can choose which mass is radiogenic: usually the layer's own mass, but possibly one differentiated component of it. The world-level sum uses each layer's `mass` attribute, which the equation-of-state solve sets, so solve the world's structure first: a layer built without a mass reports zero until then.

### Vectorized Evaluation

Three vectorized entry points are defined once on the base class, so every model inherits them.

- `calc_heating_vectorize_time(time[], mass)`: a time sweep at constant mass.
- `calc_heating_vectorize_mass(time, mass[])`: a mass sweep at constant time.
- `calc_heating_vectorize_all(time[], mass[])`: element-wise over two equal-length arrays.

The Cython wrappers accept any array-like and return a `float64` NumPy array. At the C++ level each fills a caller-supplied `std::vector<double>&`, and mismatched lengths throw `std::invalid_argument`.

```python
import numpy as np
from TidalPy.radiogenics_x import IsotopeRadiogenics

ages = np.linspace(0.0, 1.45e17, 200)                      # [s] 0 to 4.6 Gyr
model = IsotopeRadiogenics.from_dataset("llri_and_slri")
curve = model.calc_heating_vectorize_time(ages, 1.0e22)    # [W] at each age
```

### Convenience Functions

For a single number without keeping a model around, each model has a lower-case module function that builds a stack-allocated C++ model, evaluates it, and discards it.

```python
import numpy as np
from TidalPy.radiogenics_x import off, isotope, fixed

heating = fixed(time, mass, fixed_heat_production=1.0e-11)

# time and mass may each be a float or an array and are broadcast together.
curve = fixed(np.linspace(0.0, 1.0e18, 50), mass,
              fixed_heat_production=1.0e-11, average_half_life=4.47e17)
```

Signatures: 

- `off(time, mass)`;
- `isotope(time, mass, heat_production, half_lives, mass_fracs, concentrations, ref_time=0.0, names=None)`;
- `fixed(time, mass, fixed_heat_production=0.0, average_half_life=0.0, ref_time=0.0)`.

All-scalar input returns a float, anything else a `float64` array. The model parameters themselves are always constants.

Both the class methods and the convenience functions take `time` first, then `mass`, then the model parameters, matching the argument order used across the physics modules.

## Serialization

Every model supports the standard interfaces inherited from the TidalPy base class.

- `get_config_dict()` returns the model name under the key `model` plus its parameters, with isotope arrays as lists. The dict is accepted by `make_radiogenics`, so a model round-trips through it.
- `save_config(path)` writes the same content as TOML.
- `save_binary(path)` and `load_binary(path, force=False)` use the TidalPy binary format. The layer back-pointer is not serialized, so re-attach the model after loading.

The Off and Fixed models write their scalars through the shared `write_physics_binary` helper. The Isotope model writes its variable-length list, each isotope's name plus its four doubles, directly after the shared header and model name, which is the one place in the module that bypasses the scalar-only helper.

## C++ API

```cpp
#include "radiogenics_.hpp"   // pulls in radiogenics_base_.hpp

using namespace tidalpy;

c_RadiogenicsConfig config;
config.fixed_heat_production = 1.0e-11;
config.average_half_life     = 4.47e17;

const c_RadiogenicsModel model_id = c_radiogenics_model_from_name("fixed");
std::unique_ptr<c_RadiogenicsBase> model = c_find_radiogenics(model_id, config);

const double heating = model->calc_heating(time, mass);   // [W]
```

- `c_RadiogenicsBase : c_PhysicsBase`: abstract, with `calc_heating(time, mass)` pure virtual plus the three vectorized wrappers.
- Concrete models `c_OffRadiogenics`, `c_IsotopeRadiogenics`, and `c_FixedRadiogenics`, along with the `c_Isotope` value type.
- `enum class c_RadiogenicsModel { Off, Isotope, Fixed }`.
- `c_radiogenics_model_from_name(name)`: name or alias to enum, throwing `std::invalid_argument` on an unknown name.
- `c_find_radiogenics(model, config)`: heap-allocates the model as a `unique_ptr`.
- `c_radiogenics_from_binary(stream, force)`: reconstructs from a binary record, so a model attached to a layer is restored when the layer is loaded.
- `c_isotope_dataset_names()` and `c_get_isotope_dataset(name)`: the built-in catalog, already in MKS.

Binary class ids 500 through 503 are reserved for this module.

## Adding a New Model

**C++ (`TidalPy/radiogenics_x/radiogenics_.hpp`)**

1. Add any new parameters to `c_RadiogenicsConfig` with sensible defaults. The single combined config is shared by all models, and each reads only the fields it needs.
2. Add a free function implementing the heating law, guarding any half-life denominator with `rad_guard` and any growth term with `c_safe_exp` or `c_safe_pow`, so an overflow returns NaN like the existing models.
3. Add the model class deriving from `c_RadiogenicsBase`: a default constructor and one taking the config, `get_*` accessors, the `calc_heating` override, and `write_binary` / `read_binary` through the `c_PhysicsBase` helpers. Variable-length data is written directly after the header and model name, as `c_IsotopeRadiogenics` does.
4. Add the enum value, the name and alias branch in `c_radiogenics_model_from_name`, and the cases in `c_find_radiogenics` and `c_radiogenics_from_binary`.

**C++ (`TidalPy/Utilities_x/binary_x/binary_.hpp`)**

5. Add a unique `BinaryClassID` in the 50X block.

**Cython (`radiogenics.pxd` and `radiogenics.pyx`)**

6. Declare the C++ class and the new enum value in the `.pxd`.
7. Add the `cdef class` wrapper with parameter properties, the adoption branch in `make_radiogenics`, and the lower-case convenience function. The config dict comes from the C++ `append_config_entries` override, not from Python.

**Package, tests, and docs**

8. Export the class and the function from `__init__.py`, and the C++ names from `__init__.pxd`.
9. Extend `Tests/Test_Radiogenics_x/test_radiogenics_01.py`: model name, heating against an independent reference, factory and aliases, vectorization, config dict, and binary round trip.
10. Document the model here and add a changelog entry.

No build-system change is needed; `radiogenics_x.radiogenics` is already registered in `cython_extensions.json`.

## References

- Turcotte, D. L., and Schubert, G. (2001, 2002). *Geodynamics*. Radiogenic heat production rates.
- Hussmann, H., and Spohn, T. (2004). Thermal-orbital evolution of Io and Europa. *Icarus*, 171(2), 391-410. Chondritic isotope data.
- Castillo-Rogez, J. C., et al. (2007). Iapetus geophysics: Rotation rate, shape, and equatorial ridge. *Icarus*, 190(1), 179-202, [doi:10.1016/j.icarus.2007.02.018](https://doi.org/10.1016/j.icarus.2007.02.018). Long- and short-lived isotope inventories.
- McDonough, W. F., and Sun, S.-s. (1995). The composition of the Earth. *Chemical Geology*, 120(3-4), 223-253, [doi:10.1016/0009-2541(94)00140-4](https://doi.org/10.1016/0009-2541(94)00140-4). Bulk silicate Earth abundances.
