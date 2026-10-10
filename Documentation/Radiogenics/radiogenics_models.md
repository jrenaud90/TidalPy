# Radiogenic Models (`Radiogenics`)

_Updated: 2026-10-07_

A radiogenics model gives the power $Q$ \[W\] released by radioactive decay inside a layer of mass $m$ at time $t$, through `calc_heating(time, mass)`.

Time is in seconds from an epoch the caller chooses, and each model carries the reference time `ref_time` at which its rates or concentrations were quoted. Only $t - t_{\text{ref}}$ enters the heating. A model built from present-day abundances with a reference time of 4600 Myr therefore gives the heating at the birth of the Solar System at $t = 0$, and today's at $t = t_{\text{ref}}$.

## Models

For a half life $t_{1/2}$ the decay constant is $\gamma = \ln(0.5) / t_{1/2}$, a negative number whose magnitude grows as the half life shortens.

| Model (aliases) | Python class | Heating $Q$ \[W\] | Parameters |
|---|---|---|---|
| `off` (`none`) | `OffRadiogenics` | $0$ | none |
| `isotope` (`isotopes`) | `IsotopeRadiogenics` | $m \sum_i q_i f_i c_i \exp[\gamma_i (t - t_{\text{ref}})]$ | `heat_production`, `half_lives`, `mass_fracs`, `concentrations`, `ref_time` |
| `fixed` (`constant`) | `FixedRadiogenics` | $m\, q \exp[\gamma (t - t_{\text{ref}})]$ | `fixed_heat_production`, `average_half_life`, `ref_time` |

The canonical names (first in each row) are the ones written to a config dict.

- `off` returns zero for any time and mass.
- `isotope` sums the decay of any number of isotopes, each with its own half life. $q_i$ is the specific heat production of pure isotope $i$, $f_i$ its mass fraction within its parent element, and $c_i$ that element's concentration in the layer material, both at the reference time. The product $q_i f_i c_i$ is the isotope's specific heating at the reference time, per kilogram of layer, and the exponential carries it forward or backward in time. A source that quotes the isotope's own concentration rather than its element's is entered with $f_i = 1$. See [Isotope Tables](#isotope-tables).
- `fixed` applies one lumped specific rate $q$ to the whole layer, optionally with one effective half life. An `average_half_life` of zero (the default) or below means no decay, and the heating is constant for all time.

## Quick Example

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

## Attaching a Model to a `Layer`

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

The layer shares the model rather than copying it (models are not changed in place), so one model can serve several layers. Setting `radiogenics` on a layer of a solved world makes the world forget its solve. A layer without a model reports zero heating rather than raising.

The mass is an argument so the caller can choose which mass is radiogenic, usually the layer's own but possibly one differentiated part of it. `world.calc_internal_heating` uses each layer's `mass`, which the equation-of-state solve sets, so solve the structure first. With `use_heating` set, the world's thermal solve heats the layer at the model's specific rate times the local density and reports it as `layer_heating_radiogenic`, part of `layer_heating` (see [Heat Sources](../Structures/worlds/worlds.md#heat-sources)).

In a world's TOML the model is a `[layers.<name>.radiogenics]` table keyed by `model` plus any parameters (and `isotopes` or `isotope_names` for the isotope model). See the [TOML schema](../Structures/config/toml_schema.md).

## Isotope Tables

The isotope model takes its isotopes as four tables of equal length, one value per isotope:

| Parameter | Config key | Meaning |
|---|---|---|
| `heat_production` | `heat_production_w_kg` | Specific heat production of the pure isotope \[W kg$^{-1}$\]. |
| `half_lives` | `half_lives_s` | Half life \[s\]; infinite for a stable isotope. |
| `mass_fracs` | `mass_fracs` | Mass fraction of the isotope within its element at the reference time \[kg kg$^{-1}$\]; 1 when `concentrations` holds the isotope's own concentration. |
| `concentrations` | `concentrations` | Concentration of the parent element in the layer material at the reference time \[kg kg$^{-1}$\]. |

`isotope_names` gives one label per isotope, or none. The labels are written to the config dict and the binary record and do not affect the heating.

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

## Built-in Isotope Datasets

TidalPy has some isotope sets popular in the literature. List them with `available_isotope_datasets()`, inspect one with `isotope_dataset(name)`, and build a model from one with `IsotopeRadiogenics(isotopes=name)` or `make_radiogenics("isotope", {"isotopes": name})`.

| Dataset | Isotopes | Reference time | Applicability | Source |
|---|---|---|---|---|
| `modern_day_chondritic` | U238, U235, Th232, K40 | 4600 Myr | Present-day rocky and icy bodies of broadly chondritic composition. | Hussmann and Spohn (2004); Turcotte and Schubert (2001) |
| `llri` | U238, U235, Th232, K40 | 0 Myr | Ordinary chondritic rock of a body that formed too late for the short-lived isotopes to matter. | Castillo-Rogez et al. (2007) |
| `slri` | Al26, Fe60, Mn53 | 0 Myr | The short-lived isotopes of the same rock, on their own. | Castillo-Rogez et al. (2007) |
| `llri_and_slri` | U238, U235, Th232, K40, Al26, Fe60, Mn53 | 0 Myr | Early solar system thermal evolution, where the short-lived isotopes dominate the first few million years. The union of `llri` and `slri`. | Castillo-Rogez et al. (2007) |
| `bulk_silicate_earth` | U238, U235, Th232, K40 | 4600 Myr | Present-day Earth-like silicate mantles (U 20.3 ppb, Th 79.5 ppb, K 240 ppm). | McDonough and Sun (1995) concentrations; Turcotte and Schubert (2002) rates |

> [!NOTE]
> The reference times differ. For the two present-day sets $t = 0$ is Solar System formation and $t = $ `ref_time` is today. The three Castillo-Rogez et al. (2007) sets quote their concentrations at formation.

The Castillo-Rogez et al. (2007) sets are ordinary chondritic rock at the formation of the calcium-aluminum-rich inclusions (CAIs). Their Table 3 quotes each isotope's own concentration at that time, so heating at formation is the heat production times that concentration.

- The long-lived isotopes are entered with that concentration and a mass fraction of 1. The isotopic abundances in their Table 4 are present-day values, which do not hold at formation.
- The short-lived isotopes are entered as the initial isotopic ratio of their Table 5 times the element concentration it implies (26Al/27Al = 5 × 10$^{-5}$ of 1.2 wt\% aluminum, 60Fe/56Fe = 10$^{-6}$ of 22.5 wt\% iron, 53Mn/55Mn = 10$^{-5}$ of 0.257 wt\% manganese). The product reproduces Table 3.
- The paper explores 60Fe/56Fe from 10$^{-7}$ to 10$^{-6}$; the set uses the 10$^{-6}$ of its short-lived-isotope models. Where the paper quotes a range for a half life or a heat production, the set uses its middle.

`isotope_dataset(name)` returns a dataset as the isotope model's parameters: a dict with the keys `heat_production_w_kg`, `half_lives_s`, `mass_fracs`, `concentrations`, `isotope_names`, and `ref_time_s`.

### Custom Datasets

A name that is not built in is looked up in the `[radiogenics.known_isotope_data]` section of `TidalPy_Configs.toml` (empty by default; a built-in of the same name wins), and an inline dict is accepted in the same place. Both hold a `ref_time` entry and one table per isotope, keyed by its label, with `hpr` \[W kg$^{-1}$\], `half_life`, `iso_mass_fraction`, and `element_concentration`; half lives and reference times are in \[Myr\]. `isotope_dataset_parameters(dataset)` converts a built-in name, a configured name, or an inline dict into the MKS parameter dict of `isotope_dataset`.

The isotope model takes its isotopes from the first of these it is given:

1. The four tables. A model given any of them ignores `isotopes`.
2. `isotopes`: a built-in dataset name, a configured dataset name, or an inline dict.
3. Neither: the `[radiogenics] isotopes` dataset of `TidalPy_Configs.toml` (`modern_day_chondritic` by default), for both the class and `make_radiogenics`.

A dataset brings its labels and its reference time. A given `isotope_names` replaces the labels, and a given `ref_time` replaces the reference time, so the dataset's abundances then apply at that time.

## Python API

Constructors take the model's parameters by argument name or config key, as keywords, positionally in `get_parameter_info()` order, or as one table through `config=`. `make_radiogenics(model_name, config=None)` resolves a name or alias case-insensitively and builds the model from `config`; absent keys take the defaults. An unknown name, a key the model does not read, or a value outside a parameter's bounds raises `ValueError` naming the closest accepted one. For example, a `fixed` table that names a dataset raises `ValueError` saying that radiogenics model 'fixed' has no parameter 'isotopes'.

| Parameter | Config key | Default | Bounds | Model |
|---|---|---|---|---|
| `fixed_heat_production` | `fixed_heat_production_w_kg` | 0.0 | non-negative | Fixed: specific heat production at the reference time \[W kg$^{-1}$\]. |
| `average_half_life` | `average_half_life_s` | 0.0 | finite | Fixed: half life of the lumped rate \[s\]; 0 or below for no decay. |
| `heat_production` | `heat_production_w_kg` | `[]` | non-negative | Isotope: one per isotope (see [Isotope Tables](#isotope-tables)). |
| `half_lives` | `half_lives_s` | `[]` | positive or infinite | Isotope: one per isotope. |
| `mass_fracs` | `mass_fracs` | `[]` | unit interval | Isotope: one per isotope. |
| `concentrations` | `concentrations` | `[]` | non-negative | Isotope: one per isotope. |
| `ref_time` | `ref_time_s` | 0.0 | finite | Fixed and isotope: the time the rate or abundances apply at \[s\]. |

The isotope model also reads `isotopes` and `isotope_names`, which are outside its parameter table (see [Built-in Isotope Datasets](#built-in-isotope-datasets)).

Parameters read as attributes (a table as a list). Every model has `model_name`, `parameters`, `get_parameter(name)`, `get_parameter_info()`, `with_parameters(**changes)`, `get_config_dict()`, `save_config(path)`, and `save_binary(path)` / `load_binary(path)`; the isotope model adds `num_isotopes` and `isotope_names`, and its `with_parameters` also takes `isotope_names`. `radiogenics_model_names()` lists the canonical names and `radiogenics_config_keys(name)` the keys one model reads (with `isotopes` and `isotope_names` for the isotope model).

`get_config_dict()` gives `model` plus the parameters by config key, and `make_radiogenics` accepts it. The isotope model always writes its four tables, even empty ones, so a model with no isotopes is rebuilt with none rather than with the default dataset, plus `isotope_names` when labeled. A model on a layer is saved and restored with the layer.

### Vectorized Evaluation

Every model has three vectorized methods that differ only in which argument is an array. Each accepts any array-like and returns a `float64` array; a length-one input is broadcast against the other, and mismatched lengths raise `ValueError`.

- `calc_heating_vectorize_time(time[], mass)`: a time sweep at constant mass.
- `calc_heating_vectorize_mass(time, mass[])`: a mass sweep at constant time.
- `calc_heating_vectorize_all(time[], mass[])`: element-wise over two equal-length arrays.

```python
import numpy as np
from TidalPy.Radiogenics import IsotopeRadiogenics

ages = np.linspace(0.0, 1.45e17, 200)  # [s] 0 to 4.6 Gyr
model = IsotopeRadiogenics(
    isotopes="llri_and_slri")
curve = model.calc_heating_vectorize_time(ages, 1.0e22)  # [W] at each age
```

### Convenience Functions

Each model has a lower-case module function for a single evaluation without keeping a model around.

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

`isotope` uses its four tables alone: it reads no dataset and no labels, and empty tables give zero heating. An all-scalar input returns a float, anything else a `float64` array. The model parameters are always constants.

## Limits and Failure Modes

- A `fixed` average half life at or below zero is treated as infinite (no decay), not zero.
- A half life that is finite but below the module's floor is clamped to that floor, so no decay constant divides by zero.
- Evaluating a decaying model far before its reference time would overflow the exponential; the model returns NaN instead, so a bad epoch shows up as NaN heating.
- Unphysical values raise `ValueError` naming the model and the parameter: `heat_production`, `concentrations`, and `fixed_heat_production` must be finite and non-negative; `half_lives` positive, or infinite for a stable isotope; `mass_fracs` from 0 to 1; `average_half_life` and `ref_time` finite. A negative half life or abundance would otherwise make the heating grow without bound or turn negative.
- The isotope model refuses tables of different lengths and labels that are not one per isotope.

## C++ API

The models are in `TidalPy/Radiogenics/radiogenics_.hpp`, and their base in `radiogenics_base_.hpp` (namespace `tidalpy`, header only).

```cpp
#include "radiogenics_.hpp"  // pulls in radiogenics_base_.hpp

using namespace tidalpy;

c_ParamMap params;  // Config key to values; one value for a scalar
params["fixed_heat_production_w_kg"] = {1.0e-11};
params["average_half_life_s"] = {4.47e17};

std::unique_ptr<c_RadiogenicsBase> model = c_find_radiogenics("fixed", params);
const double heating = model->calc_heating(time, mass);  // [W]
```

- `c_RadiogenicsBase` (abstract): `calc_heating(time, mass)`, `get_ref_time()` (zero for a model without one), and `calc_heating_vectorize(time, mass, out_heating)`, which broadcasts a length-one input into a caller-supplied `std::vector<double>&`.
- `c_OffRadiogenics`, `c_FixedRadiogenics`, and `c_IsotopeRadiogenics`, which adds `get_num_isotopes()` and `get_isotope_names()`.
- `c_find_radiogenics(name, params)` and, for a labeled isotope model, `c_make_isotope_radiogenics(params, isotope_names)` return a `unique_ptr`, throwing `std::invalid_argument` for an unknown name or parameter or an out-of-bounds value. C++ takes only the four tables; dataset names are resolved in Python.
- `c_radiogenics_from_binary(stream, force)`, `c_radiogenics_canonical_name(name)`, and `c_radiogenics_model_names()`.
- `c_Isotope` (`name`, `heat_production`, `half_life`, `mass_frac`, `concentration`), `c_IsotopeDataset` (`isotopes`, `ref_time`), `c_isotope_dataset_names()`, and `c_get_isotope_dataset(name)`: the built-in catalog, in MKS.

### Adding a New Model

A model `Foo` is a `c_SpecModel<c_FooRadiogenics, c_RadiogenicsBase>` in `radiogenics_.hpp` with a `parameter_specs()` table, `C_CLASS_ID`, and `calc_heating`. Guard half-life denominators with `c_guard_denominator` and growth terms with `c_safe_exp` or `c_safe_pow` so an overflow returns NaN, and override `get_ref_time` if it has a reference time. Then add a registry row, a `BinaryClassID` in the 50X block, a Cython class and `foo(time, mass, ...)` function exported from `TidalPy/Radiogenics/__init__.py`, tests in `Tests/Test_Radiogenics/`, and a section here.

## References

- Turcotte, D. L., and Schubert, G. (2001, 2002). *Geodynamics*. Radiogenic heat production rates.
- Hussmann, H., and Spohn, T. (2004). Thermal-orbital evolution of Io and Europa. *Icarus*, 171(2), 391-410. Chondritic isotope data.
- Castillo-Rogez, J. C., et al. (2007). Iapetus geophysics: Rotation rate, shape, and equatorial ridge. *Icarus*, 190(1), 179-202, [doi:10.1016/j.icarus.2007.02.018](https://doi.org/10.1016/j.icarus.2007.02.018). Long- and short-lived isotope inventories.
- McDonough, W. F., and Sun, S.-s. (1995). The composition of the Earth. *Chemical Geology*, 120(3-4), 223-253, [doi:10.1016/0009-2541(94)00140-4](https://doi.org/10.1016/0009-2541(94)00140-4). Bulk silicate Earth abundances.
