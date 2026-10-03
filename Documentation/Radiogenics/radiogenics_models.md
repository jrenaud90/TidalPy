# Radiogenic Models (`Radiogenics`)

_Updated: 2026-10-02_

A radiogenics model utilizes a layer of mass $m$ at time $t$ to find how much power is being released inside it by radioactive decay. The heating $Q$ \[W\] is returned by `calc_heating(time, mass)`.

Time is measured in seconds from an epoch the caller chooses, and each model carries the reference time `ref_time` at which its rates or concentrations were quoted. Only the difference $t - t_{\text{ref}}$ is used in the calculations, so a model built from present-day abundances with a reference time of 4600 Myr is evaluated at $t = 0$ to get the heating at the birth of the Solar System, or at $t = t_{\text{ref}}$ to get today's.

## Inheritance

```
c_TidalPyBaseClass
  └── c_PhysicsBase
        └── c_RadiogenicsBase  (abstract)
              ├── c_OffRadiogenics      aliases "off", "none"
              ├── c_IsotopeRadiogenics  aliases "isotope", "isotopes"
              └── c_FixedRadiogenics    aliases "fixed", "constant"
```

`c_RadiogenicsBase` declares `calc_heating(time, mass)` pure virtual and supplies one vectorized method, `calc_heating_vectorize(time, mass, out)`, that loops over it, so a new model only has to implement the heating law. A model may override the vectorized method; the isotope model does, summing one isotope at a time over every point. Every concrete model derives from `c_SpecModel`, which gives it its parameters, config dict, and binary record from one table. The canonical names are `off`, `isotope`, and `fixed`, and those are the names written to a config dict.

## Models

For a half life $t_{1/2}$ the decay constant is $\gamma = \ln(0.5) / t_{1/2}$, a negative number whose magnitude grows as the half life shortens.

| Model (aliases) | Heating $Q$ \[W\] | Parameters |
|---|---|---|
| `off` (`none`) | $0$ | none |
| `isotope` (`isotopes`) | $m \sum_i q_i f_i c_i \exp[\gamma_i (t - t_{\text{ref}})]$ | `heat_production`, `half_lives`, `mass_fracs`, `concentrations`, `ref_time` |
| `fixed` (`constant`) | $m\, q \exp[\gamma (t - t_{\text{ref}})]$ | `fixed_heat_production`, `average_half_life`, `ref_time` |

Here $q_i$ is the specific heat production of pure isotope $i$, $f_i$ its mass fraction within its parent element, and $c_i$ that element's concentration in the layer material, both at the reference time. The product $q_i f_i c_i$ is the specific heating the isotope contributed at the reference time, per kilogram of layer, and the exponential carries it forward or backward in time. A source that quotes the isotope's own concentration rather than its element's is entered with $f_i = 1$.

### Off

Returns zero for any time and mass.

### Isotope

Sums the decay of any number of isotopes, each with its own half life. Its first four parameters are tables with one value per isotope (see [Isotope Tables](#isotope-tables)).

### Fixed

Applies one lumped specific rate to the whole layer, optionally with a single effective half life. Setting `average_half_life` to zero (the default) or to any non-positive value means no decay at all, and the heating is then constant for all time.

### Behavior at the Limits

A `fixed` model's average half life at or below zero is treated as infinite, not zero, so the constant-rate case uses the same formula. A half life that is finite but smaller than the module's floor is clamped to that floor, so no decay constant is ever divided by zero.

The models refuse values that are not physical with `ValueError`, which names the model and the parameter:

- `heat_production`, `concentrations`, and `fixed_heat_production`: finite and non-negative.
- `half_lives`: positive, or infinite for a stable isotope.
- `mass_fracs`: from 0 to 1.
- `average_half_life` and `ref_time`: finite.

A negative half life or abundance would otherwise make the heating grow without bound or turn negative. The isotope model also refuses tables of different lengths and labels that are not one per isotope.

Evaluating a model far before its reference time asks for an exponential that would overflow. Both decaying models guard against this and return NaN, so a bad epoch shows up as NaN heating.

## Isotope Tables

The isotope model takes its isotopes as four tables of equal length, one value per isotope:

| Parameter | Config key | Meaning |
|---|---|---|
| `heat_production` | `heat_production_w_kg` | Specific heat production of the pure isotope \[W kg$^{-1}$\]. |
| `half_lives` | `half_lives_s` | Half life \[s\]; infinite for a stable isotope. |
| `mass_fracs` | `mass_fracs` | Mass fraction of the isotope within its element at the reference time \[kg kg$^{-1}$\]; 1 when `concentrations` holds the isotope's own concentration. |
| `concentrations` | `concentrations` | Concentration of the parent element in the layer material at the reference time \[kg kg$^{-1}$\]. |

`isotope_names` labels the isotopes, one label per isotope, or none. The labels are written to the config dict and the binary record and take no part in the heating.

```python
from TidalPy.Radiogenics import IsotopeRadiogenics

model = IsotopeRadiogenics(
    heat_production=[9.48e-5, 2.69e-5],  # [W/kg] of the pure isotope
    half_lives=[4.47e17, 1.40e18],  # [s]
    mass_fracs=[0.9928, 0.9998],  # [kg/kg] within the element
    concentrations=[0.012e-6, 0.04e-6],  # [kg/kg] of the element in the material
    ref_time=0.0,  # [s]
    isotope_names=["U238", "Th232"])  # One label per isotope

print(model.num_isotopes)  # 2
print(model.isotope_names)  # ['U238', 'Th232'], or [] for unlabeled isotopes
print(model.heat_production)  # [9.48e-05, 2.69e-05], a list, as are the other three tables
```

In C++, `c_IsotopeRadiogenics` holds the four tables as `std::vector<double>` members and caches each isotope's specific heating at the reference time and its decay constant.

## Built-in Isotope Datasets

TidalPy provides some sets of isotopes popular in the literature. List them with `available_isotope_datasets()`, inspect one with `isotope_dataset(name)`, and build a model from one with `IsotopeRadiogenics(isotopes=name)` or `make_radiogenics("isotope", {"isotopes": name})`.

| Dataset | Isotopes | Reference time | Applicability | Source |
|---|---|---|---|---|
| `modern_day_chondritic` | U238, U235, Th232, K40 | 4600 Myr | Present-day rocky and icy bodies of broadly chondritic composition. | Hussmann and Spohn (2004); Turcotte and Schubert (2001) |
| `llri` | U238, U235, Th232, K40 | 0 Myr | Ordinary chondritic rock of a body that formed too late for the short-lived isotopes to matter. | Castillo-Rogez et al. (2007) |
| `slri` | Al26, Fe60, Mn53 | 0 Myr | The short-lived isotopes of the same rock, on their own. | Castillo-Rogez et al. (2007) |
| `llri_and_slri` | U238, U235, Th232, K40, Al26, Fe60, Mn53 | 0 Myr | Early solar system thermal evolution, where the short-lived isotopes dominate the first few million years. The union of `llri` and `slri`. | Castillo-Rogez et al. (2007) |
| `bulk_silicate_earth` | U238, U235, Th232, K40 | 4600 Myr | Present-day Earth-like silicate mantles (U 20.3 ppb, Th 79.5 ppb, K 240 ppm). | McDonough and Sun (1995) concentrations; Turcotte and Schubert (2002) rates |

> [!NOTE]
> Notice that the reference times differ. The two present-day sets quote concentrations at 4600 Myr, so evaluating them at $t = 0$ gives the heating at Solar System formation and evaluating at $t = $ `ref_time` gives today's. The three Castillo-Rogez et al. (2007) sets quote their concentrations at formation instead, so they are evaluated at formation time instead.

The Castillo-Rogez et al. (2007) sets are ordinary chondritic rock at the formation of the calcium-aluminum-rich inclusions (CAIs). Their Table 3 quotes each isotope's own concentration at that time, so heating at formation is the heat production times that concentration.

- The long-lived isotopes are entered with that concentration and a mass fraction of 1. The isotopic abundances in their Table 4 are present-day values, which do not hold at formation.
- The short-lived isotopes are entered as the initial isotopic ratio of their Table 5 times the element concentration it implies (26Al/27Al = 5 × 10$^{-5}$ of 1.2 wt\% aluminum, 60Fe/56Fe = 10$^{-6}$ of 22.5 wt\% iron, 53Mn/55Mn = 10$^{-5}$ of 0.257 wt\% manganese). The product reproduces Table 3.
- The paper explores 60Fe/56Fe from 10$^{-7}$ to 10$^{-6}$; the set uses the 10$^{-6}$ of its short-lived-isotope models. Where the paper quotes a range for a half life or a heat production, the set uses its middle.

`isotope_dataset(name)` returns a built-in dataset as the isotope model's parameters: a dict with the keys `heat_production_w_kg`, `half_lives_s`, `mass_fracs`, `concentrations`, `isotope_names`, and `ref_time_s`.

A name that is not one of the built-ins is looked up in the `[radiogenics.known_isotope_data]` section of `TidalPy_Configs.toml` (`TidalPy.config['radiogenics']['known_isotope_data']`, empty by default), and an inline dict is accepted in the same place. Those two sources store half lives and reference times in mega-years \[Myr\] and the heat production rate in \[W kg$^{-1}$\], with a `ref_time` entry and one table per isotope, keyed by its label, holding `hpr`, `half_life`, `iso_mass_fraction`, and `element_concentration`. The built-in catalog is always preferred over a same-named config entry. `isotope_dataset_parameters(dataset)` converts any of the three (a built-in name, a configured name, or an inline dict) into the same parameter dict as `isotope_dataset`, in MKS.

The isotope model takes its isotopes from the first of these that it is given:

1. The four tables. A model given any of them ignores `isotopes`.
2. `isotopes`: a built-in dataset name, a configured dataset name, or an inline dict.
3. Neither: the `[radiogenics] isotopes` dataset of `TidalPy_Configs.toml` (`modern_day_chondritic` by default). The class and `make_radiogenics` both take it.

A dataset brings its labels and its reference time. A given `isotope_names` replaces the labels, and a given `ref_time` replaces the reference time, so the dataset's abundances then apply at that time.

## Python API

```python
from TidalPy.Radiogenics import FixedRadiogenics, IsotopeRadiogenics, make_radiogenics

mass = 1.0e22  # [kg]
time = 1.0e17  # [s] from the caller's epoch

# A lumped rate with no decay
model = FixedRadiogenics(
    fixed_heat_production=1.0e-11)
heating = model.calc_heating(time, mass)  # [W]

# A literature isotope set, evaluated at its own reference time (present day)
chondritic = IsotopeRadiogenics(
    isotopes="modern_day_chondritic")
heating_today = chondritic.calc_heating(chondritic.ref_time, mass)  # [W]

# Name factory: case-insensitive, aliases accepted
model = make_radiogenics("constant", {"fixed_heat_production_w_kg": 2.0e-11})
model = make_radiogenics("isotope", {"isotopes": "bulk_silicate_earth"})
```

Constructors take the model's parameters by argument name or config key, as keywords, positionally in the order `get_parameter_info()` lists them, or as one table through `config=`. `make_radiogenics(model_name, config=None)` resolves a name or alias case-insensitively and builds the model from `config`. Absent keys take the model's defaults. An unknown name, a key the model does not read, or a value outside a parameter's bounds raises `ValueError` naming the closest accepted one. A `fixed` table that names a dataset, for example, raises `ValueError` saying that radiogenics model 'fixed' has no parameter 'isotopes'.

| Parameter | Config key | Default | Bounds | Model |
|---|---|---|---|---|
| `fixed_heat_production` | `fixed_heat_production_w_kg` | 0.0 | non-negative | Fixed: specific heat production at the reference time \[W kg$^{-1}$\]. |
| `average_half_life` | `average_half_life_s` | 0.0 | finite | Fixed: half life of the lumped rate \[s\]; 0 or below for no decay. |
| `heat_production` | `heat_production_w_kg` | `[]` | non-negative | Isotope: one per isotope (see [Isotope Tables](#isotope-tables)). |
| `half_lives` | `half_lives_s` | `[]` | positive or infinite | Isotope: one per isotope. |
| `mass_fracs` | `mass_fracs` | `[]` | unit interval | Isotope: one per isotope. |
| `concentrations` | `concentrations` | `[]` | non-negative | Isotope: one per isotope. |
| `ref_time` | `ref_time_s` | 0.0 | finite | Fixed and isotope: the time the rate or abundances apply at \[s\]. |

The isotope model also reads two keys outside its parameter table, `isotopes` and `isotope_names` (see [Built-in Isotope Datasets](#built-in-isotope-datasets)).

The parameters read as attributes (`model.ref_time`, `model.heat_production`), a table as a list, and every model has `model_name`, `parameters`, `get_parameter(name)`, `get_parameter_info()`, `with_parameters(**changes)`, `get_config_dict()`, and `save_config(path)`. `off` has no parameters. The isotope model adds `num_isotopes` and `isotope_names`, and its `with_parameters` also takes `isotope_names`. `radiogenics_model_names()` lists the canonical names and `radiogenics_config_keys(name)` the keys one model reads, `isotopes` and `isotope_names` included for the isotope model.

### Factory Internals

At the C++ level `c_radiogenics_registry()` lists each model's names (canonical first, then aliases), its binary class id, and its constructor. `c_find_radiogenics(name, params)` builds a model from a `c_ParamMap` and throws `std::invalid_argument` for an unknown name or parameter, and `c_radiogenics_canonical_name` and `c_radiogenics_model_names` resolve names. The isotope labels are not parameters, so a labeled isotope model is built with `c_make_isotope_radiogenics(params, isotope_names)`. The Python `make_radiogenics` passes the config dict to the model's class. `IsotopeRadiogenics` resolves a dataset (named, inline, or the configured default) into the four tables and the labels before it builds the C++ model, so the C++ models only ever see tables.

### Attaching a Model to a `Layer`

```python
from TidalPy.Radiogenics import IsotopeRadiogenics
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds import TerrestrialWorld

radiogenics = IsotopeRadiogenics(
    isotopes="modern_day_chondritic")
time = radiogenics.ref_time  # [s] the dataset's reference time

mantle = Layer(
    "mantle",
    0,
    0.0,
    1.8e6,
    material="simple_rock",
    temperature=1600.0,
    use_heating=True,  # Heated by its model in a thermal solve
    radiogenics=radiogenics)  # A model, a model name, or a config table
print(mantle.radiogenics_set)  # True
print(mantle.calc_radiogenic_heating(time, 8.0e22))  # [W] for the mass supplied

world = TerrestrialWorld("Io-like", 1.8e6, 8.0e22)
world.add_layer(mantle)
world.solve_eos()  # Sets each layer's mass
print(world.calc_internal_heating(time))  # [W] summed over all layers

world.mantle.radiogenics = {"model": "fixed", "fixed_heat_production_w_kg": 5.0e-12}  # Replace it
world.mantle.radiogenics = None  # Remove it
```

The `radiogenics` argument and property take a model, a model name, or a config table with a `model` key, and `None` removes the model. The layer shares the model rather than copying it (models are not changed in place), so one model can be attached to several layers and stays usable afterwards. `layer.radiogenics` returns the attached model, or `None`. Setting it on a layer of a solved world makes the world forget its solve. A layer without a model reports zero heating rather than raising, and a world sums only the layers that carry one. A layer with `use_heating` set also feeds its model to the world's thermal EOS solve, which heats the layer at the model's specific rate times the local density and reports it as `layer_heating_radiogenic`, part of `layer_heating` (see [Heat Sources](../Structures/worlds/worlds.md#heat-sources)).

The mass is an argument so the caller can choose which mass is radiogenic. This is usually the layer's own mass, but possibly one differentiated part of it. The world-level sum uses each layer's `mass` attribute, which the equation-of-state solves, so solve the world's structure first.

The declarative form is a `[layers.<name>.radiogenics]` table in a world's TOML, keyed by `model` plus any parameters (and, for the isotope model, `isotopes` or `isotope_names`). See the [TOML schema](../Structures/config/toml_schema.md) and [Layer](../Structures/layers/layer.md).

### Vectorized Evaluation

Three vectorized Python methods are defined once on the base class, so every model inherits them. They differ only in which argument they take as an array.

- `calc_heating_vectorize_time(time[], mass)`: a time sweep at constant mass.
- `calc_heating_vectorize_mass(time, mass[])`: a mass sweep at constant time.
- `calc_heating_vectorize_all(time[], mass[])`: element-wise over two equal-length arrays.

Each accepts any array-like and returns a `float64` NumPy array. All three call the one C++ method, `calc_heating_vectorize(time, mass, out)`, which broadcasts a length-one input against the other and fills the caller-supplied `std::vector<double>&`. Mismatched lengths raise `ValueError`.

```python
import numpy as np
from TidalPy.Radiogenics import IsotopeRadiogenics

ages = np.linspace(0.0, 1.45e17, 200)  # [s] 0 to 4.6 Gyr
model = IsotopeRadiogenics(
    isotopes="llri_and_slri")
curve = model.calc_heating_vectorize_time(ages, 1.0e22)  # [W] at each age
```

### Convenience Functions

For a single number without keeping a model around, each model has a lower-case module function that builds the C++ model for the one call, evaluates it, and discards it.

```python
import numpy as np
from TidalPy.Radiogenics import off, isotope, fixed

mass = 1.0e22  # [kg]
time = 1.0e17  # [s]
heating = fixed(
    time,
    mass,
    fixed_heat_production=1.0e-11)

# time and mass may each be a float or an array and are broadcast together
curve = fixed(
    np.linspace(0.0, 1.0e18, 50),
    mass,
    fixed_heat_production=1.0e-11,
    average_half_life=4.47e17)
```

Signatures:

- `off(time, mass)`
- `isotope(time, mass, heat_production=(), half_lives=(), mass_fracs=(), concentrations=(), ref_time=0.0)`
- `fixed(time, mass, fixed_heat_production=0.0, average_half_life=0.0, ref_time=0.0)`

`isotope` takes its isotopes from its four tables alone: it reads no dataset and no labels, and empty tables give zero heating. An all-scalar input returns a float, anything else a `float64` array. The model parameters themselves are always constants.

Both the class methods and the convenience functions take `time` first, then `mass`, then other model parameters.

## C++ API

The models are in `TidalPy/Radiogenics/radiogenics_.hpp` and their base in `radiogenics_base_.hpp` (namespace `tidalpy`, header only).

```cpp
#include "radiogenics_.hpp"  // pulls in radiogenics_base_.hpp

using namespace tidalpy;

c_ParamMap params;  // Config key to values; one value for a scalar
params["fixed_heat_production_w_kg"] = {1.0e-11};
params["average_half_life_s"] = {4.47e17};

std::unique_ptr<c_RadiogenicsBase> model = c_find_radiogenics("fixed", params);
const double heating = model->calc_heating(time, mass);  // [W]
```

- `c_RadiogenicsBase : c_PhysicsBase`: abstract, with `calc_heating(time, mass)` pure virtual, the virtual `get_ref_time()` (zero for a model without one), and the virtual `calc_heating_vectorize(time, mass, out_heating)`, which loops over `calc_heating`.
- `c_OffRadiogenics`, `c_IsotopeRadiogenics`, and `c_FixedRadiogenics`: each a `c_SpecModel<Model, c_RadiogenicsBase>` with its `parameter_specs()` table and `C_CLASS_ID`. `c_IsotopeRadiogenics` adds `get_num_isotopes()`, `get_isotope_names()`, and a constructor that takes the labels. It checks the table lengths and the label count in `p_validate`, and writes its labels after its parameters through `p_write_extra` and `p_read_extra`.
- `c_radiogenics_registry()`: each model's names, binary class id, and constructor.
- `c_find_radiogenics(name, params)` and `c_make_isotope_radiogenics(params, isotope_names)`: build a model as a `unique_ptr`, throwing `std::invalid_argument` for an unknown name or parameter or a value outside its bounds.
- `c_radiogenics_from_binary(stream, force)`: reconstructs from a binary record, so a model attached to a layer is restored when the layer is loaded.
- `c_radiogenics_canonical_name(name)` and `c_radiogenics_model_names()`: name lookups.
- `c_Isotope` (one row of a built-in dataset: `name`, `heat_production`, `half_life`, `mass_frac`, `concentration`), `c_IsotopeDataset` (`isotopes` and `ref_time`), `c_isotope_dataset_names()`, and `c_get_isotope_dataset(name)`: the built-in catalog, already in MKS.

## Serialization

| Call | Result |
|---|---|
| `get_config_dict()` | `model` plus the model's parameters by config key. The isotope model always writes its four tables, empty ones too, so a model with no isotopes is rebuilt with none rather than with the configured default dataset. It writes `isotope_names` when its isotopes are labeled. |
| `save_config(path)` | That dict written as TOML. |
| `save_binary(path)` / `load_binary(path)` | The model's TidalPy binary record: its parameters written by key, then, for the isotope model, its labels. |

The dict is accepted by `make_radiogenics`, so a model round-trips through it. A radiogenics model on a layer is saved and restored with the layer, in its binary record and in its config dict (as its `radiogenics` table).

## Adding a New Model

To add a radiogenics model named `Foo`:

**C++ (`TidalPy/Radiogenics/radiogenics_.hpp`)**

1. Add `c_FooRadiogenics : public c_SpecModel<c_FooRadiogenics, c_RadiogenicsBase>`: its `parameter_specs()` table (argument name, config key, member, default, bounds, and a one-line description), `C_CLASS_ID`, two constructors that call `p_initialize`, and `calc_heating`. Guard any half-life denominator with `c_guard_denominator` (`constants_.hpp`) and any growth term with `c_safe_exp` or `c_safe_pow`, so an overflow returns NaN like the existing models. Override `get_ref_time` when the model has a reference time, `p_validate` for checks across parameters, and `p_update_derived` for cached values.
2. Add one row to `c_radiogenics_registry()` with the model's names and aliases.

**C++ (`TidalPy/Utilities/binary/binary_.hpp`)**

3. Reserve a unique `BinaryClassID::FooRadiogenics` in the 50X block.

**Cython (`radiogenics.pxd` and `radiogenics.pyx`)**

4. Declare `cdef class FooRadiogenics(RadiogenicsBase)` in the `.pxd`. In the `.pyx`, add the class with a docstring and `MODEL_NAME = "foo"`, include it in the `ModelFamily` list, and add a lower-case `foo(time, mass, ...)` convenience function. Add its heating tests to `Tests/Test_Radiogenics/`. The generic tests in `Tests/Test_Utilities/Test_Classes/test_spec_models_01.py` cover its parameters, config, binary record, and errors without changes.

**Package, tests, and documentation**

5. Export `FooRadiogenics` and `foo` from `TidalPy/Radiogenics/__init__.py`.
6. Extend `Tests/Test_Radiogenics/test_radiogenics_01.py`: heating against an independent reference, the factory and aliases, vectorization, the config dict, and the binary round trip.
7. Document the model here with its formula, parameters, and references.

No build-system change is needed; `Radiogenics.radiogenics` is already registered in `cython_extensions.json`.

## References

- Turcotte, D. L., and Schubert, G. (2001, 2002). *Geodynamics*. Radiogenic heat production rates.
- Hussmann, H., and Spohn, T. (2004). Thermal-orbital evolution of Io and Europa. *Icarus*, 171(2), 391-410. Chondritic isotope data.
- Castillo-Rogez, J. C., et al. (2007). Iapetus geophysics: Rotation rate, shape, and equatorial ridge. *Icarus*, 190(1), 179-202, [doi:10.1016/j.icarus.2007.02.018](https://doi.org/10.1016/j.icarus.2007.02.018). Long- and short-lived isotope inventories.
- McDonough, W. F., and Sun, S.-s. (1995). The composition of the Earth. *Chemical Geology*, 120(3-4), 223-253, [doi:10.1016/0009-2541(94)00140-4](https://doi.org/10.1016/0009-2541(94)00140-4). Bulk silicate Earth abundances.
