# Luminosity Models (`Stellar`)

_Updated: 2026-10-02_

A luminosity model maps a star's mass onto its luminosity $L$ \[W\]. That sets the effective temperature through the Stefan-Boltzmann law, and, once the star is placed in a `System`, the flux and equilibrium temperature of every world orbiting it.

All three models share the same conversions between luminosity and effective temperature, which depend only on the radius and not on the mass relation:

$$L = 4 \pi R^2 \sigma T^4, \qquad T = \left(\frac{L}{4 \pi R^2 \sigma}\right)^{1/4}$$

What separates the models is $L(M)$.

## Inheritance

```
c_TidalPyBaseClass
  └── c_PhysicsBase
        └── c_LuminosityBase  (abstract)
              ├── c_FixedLuminosity     aliases "fixed", "constant"
              ├── c_MassToLuminosity    aliases "mass_to_luminosity", "cuntz_wang", "cw"
              └── c_PowerLawLuminosity  aliases "power_law", "powerlaw"
```

The base class declares `calc_luminosity(mass)` pure virtual and supplies the Stefan-Boltzmann conversions and the vectorized mass sweep, so a new relation is a single method. Every concrete model derives from `c_SpecModel`, which gives it its parameters, config dict, and binary record from one table. The canonical names are `fixed`, `mass_to_luminosity`, and `power_law`, and those are the names written to a config dict. The Cython classes mirror the hierarchy: `LuminosityBase`, `FixedLuminosity`, `MassToLuminosity`, `PowerLawLuminosity`.

## Models

| Model (aliases) | Relation |
|---|---|
| `fixed` (`constant`) | $L$ is supplied directly and does not depend on mass. |
| `mass_to_luminosity` (`cuntz_wang`, `cw`) | The piecewise main-sequence relation below. |
| `power_law` (`powerlaw`) | $L = L_\odot\, c \left(M / M_\odot\right)^{p}$, with $c = 1$ and $p = 3.5$ by default. |

### Fixed

Reports the luminosity it was given for any mass. Used when the star's luminosity is a measured quantity or when performing a sweep with luminosity as an independent variable.

### Mass to Luminosity

A piecewise function of the star's mass, with one power law per mass range. Below, $x = M / M_\odot$.

| Range | Relation |
|---|---|
| $x < 0.2$ | $L = 0.23\, L_\odot\, x^{2.3}$ |
| $0.2 \le x < 0.85$ | $L = L_\odot\, x^{\,p(x)}$, with the Cuntz and Wang (2018) polynomial exponent $p(x) = -141.7 x^4 + 232.4 x^3 - 129.1 x^2 + 33.29 x + 0.215$ |
| $0.85 \le x < 2$ | $L = L_\odot\, x^{4}$ |
| $2 \le x < 55$ | $L = 1.4\, L_\odot\, x^{3.5}$ |
| $x \ge 55$ | $L = 3.2 \times 10^{4}\, L_\odot\, x$ |

The mass-ratio exponent in the second branch helps the fit across the whole M-dwarf range. The branches are the published fits as they stand, and not every joint is continuous:

| Joint | Below | Above | Step |
|---|---|---|---|
| $x = 0.2$ | $5.68 \times 10^{-3}\, L_\odot$ | $4.62 \times 10^{-3}\, L_\odot$ | $-18.7$ percent |
| $x = 0.85$ | | | continuous |
| $x = 2$ | $16.0\, L_\odot$ | $15.8\, L_\odot$ | $-1.0$ percent |
| $x = 55$ | $1.73 \times 10^{6}\, L_\odot$ | $1.76 \times 10^{6}\, L_\odot$ | $+1.9$ percent |

At $x = 0.2$ the luminosity falls as the mass rises across the joint, so a star evolved through it (or a sweep across it) sees a step. Keep a model's masses on one side of a joint, or blend the branches, when that matters.

### Power Law

A single power law with a caller-set prefactor and exponent. Use it to reproduce a paper that adopted one, or to isolate the effect of the mass-luminosity slope by varying $p$ directly. The defaults reproduce the textbook $L \propto M^{3.5}$ scaling for solar-type stars.

### Behavior at the Limits

A non-positive mass returns NaN rather than raising, and so does a non-positive luminosity or radius in the temperature conversions. This keeps a bad entry in a large mass sweep visible as NaN in the output rather than aborting the sweep. The Stefan-Boltzmann constant comes from the shared TidalPy config.

## Python API

```python
import numpy as np
from TidalPy.Stellar import (
    FixedLuminosity, MassToLuminosity, PowerLawLuminosity,
    make_luminosity, fixed, mass_to_luminosity, power_law)
from TidalPy.constants import mass_solar, radius_solar

model = MassToLuminosity()
luminosity = model.calc_luminosity(mass_solar)                          # [W]
temperature = model.calc_effective_temperature(mass_solar, radius_solar)  # [K]

# The two conversions that do not involve the mass relation.
model.calc_luminosity_from_temperature(5772.0, radius_solar)   # [W]
model.calc_temperature_from_luminosity(luminosity, radius_solar)  # [K]

# calc_luminosity accepts an array of masses and returns a float64 array.
masses = np.array([0.1, 0.5, 1.0, 5.0, 30.0]) * mass_solar
curve = model.calc_luminosity(masses)

# Name factory: case-insensitive, aliases accepted.
model = make_luminosity("cw")
model = make_luminosity("power_law", {"power_law_coeff": 1.0, "power_law_exponent": 4.0})
model = make_luminosity("constant", {"luminosity_w": 3.828e26})

# Build, evaluate, and discard in one call.
luminosity = mass_to_luminosity(mass_solar)
luminosity = power_law(mass_solar, coeff=1.0, exponent=4.0)
luminosity = fixed(mass_solar, luminosity=3.828e26)
```

Constructors take the model's parameters by argument name or config key, as keywords, positionally in the order `get_parameter_info()` lists them, or as one table through `config=`, so `PowerLawLuminosity(exponent=4.0)` and `PowerLawLuminosity(power_law_exponent=4.0)` build the same model. `make_luminosity(model_name, config=None)` resolves a name or alias case-insensitively and builds the model from `config`. Absent keys take the model's defaults. An unknown name, a key the model does not read, or a value outside a parameter's bounds raises `ValueError` naming the closest accepted one, so a misspelling fails loudly instead of silently building a default model. A key that another luminosity model reads is refused too: a `fixed` table that carries `power_law_exponent` raises.

| Parameter | Config key | Default | Bounds | Model |
|---|---|---|---|---|
| `luminosity` | `luminosity_w` | 0.0 | non-negative | Fixed: the luminosity reported at every mass \[W\]. |
| `coeff` | `power_law_coeff` | 1.0 | positive | Power law: the prefactor $c$, in solar luminosities. |
| `exponent` | `power_law_exponent` | 3.5 | finite | Power law: the exponent $p$ of the mass ratio. |

The parameters read as attributes (`model.exponent`), and every model has `model_name`, `parameters`, `get_parameter(name)`, `get_parameter_info()`, `with_parameters(**changes)`, `get_config_dict()`, and `save_config(path)`. `mass_to_luminosity` has no parameters. `luminosity_model_names()` lists the canonical names and `luminosity_config_keys(name)` the keys one model reads.

The convenience functions `fixed(mass, luminosity=0.0)`, `mass_to_luminosity(mass)`, and `power_law(mass, coeff=1.0, exponent=3.5)` each build the C++ model for the one call, evaluate it, and discard it. The mass may be a float or an array; the model parameters are always constants.

### Attaching a Model to a `StarWorld`

```python
from TidalPy.Structures.worlds.stellar import StarWorld
from TidalPy.Stellar import MassToLuminosity
from TidalPy.constants import mass_solar, radius_solar

star = StarWorld("sun", radius_solar, mass_solar)
star.set_luminosity_model(MassToLuminosity())   # the star shares the model
star.luminosity_model_set                       # True

star.calc_luminosity_from_mass()                # [W]
star.calc_effective_temperature_from_mass()     # [K]
star.update_luminosity_from_mass()              # writes both onto the star's own fields
star.luminosity, star.effective_temperature
```

`set_luminosity_model` shares the model with the star rather than copying it (models are not changed in place), so one model can serve several stars. The two mass-derived calculations raise `RuntimeError` when no model has been attached; `set_luminosity` and `set_effective_temperature` remain available and keep the star's two scalars consistent through the Stefan-Boltzmann relation without one.

Once the star is part of a `System`, `calc_insolation_flux(world)` and `calc_equilibrium_temperature(world)` use its luminosity to give each orbiting world its incident flux and gray-body equilibrium temperature.

## Serialization

| Call | Result |
|---|---|
| `get_config_dict()` | `model` plus the model's parameters by config key. |
| `save_config(path)` | That dict written as TOML. |
| `save_binary(path)` / `load_binary(path)` | The model's TidalPy binary record, its parameters written by key. |

The dict is accepted by `make_luminosity`, so a model round-trips through it. A star's luminosity model is saved and restored with the star's binary record.

## C++ API

The models are in `TidalPy/Stellar/luminosity_.hpp` and their base in `luminosity_base_.hpp` (namespace `tidalpy`, header only).

```cpp
#include "luminosity_.hpp"  // pulls in luminosity_base_.hpp

using namespace tidalpy;

c_ParamMap params;  // Config key to values; one value for a scalar
params["power_law_exponent"] = {4.0};

std::unique_ptr<c_LuminosityBase> model = c_find_luminosity("power_law", params);
const double luminosity = model->calc_luminosity(mass);  // [W]
const double temperature = model->calc_effective_temperature(mass, radius);  // [K]
```

- `c_LuminosityBase : c_PhysicsBase`: abstract, with `calc_luminosity(mass)` pure virtual plus the shared Stefan-Boltzmann conversions and `calc_luminosity_vectorize_mass`.
- `c_FixedLuminosity`, `c_MassToLuminosity`, and `c_PowerLawLuminosity`: each a `c_SpecModel<Model, c_LuminosityBase>` with its `parameter_specs()` table and `C_CLASS_ID`. The two relations are also free functions, `lum_from_mass(mass)` and `lum_from_power_law(mass, coeff, exponent)`.
- `c_luminosity_registry()`: each model's names (canonical first, then aliases), binary class id, and constructor.
- `c_find_luminosity(name, params)`: builds a model from a `c_ParamMap` as a `unique_ptr`, throwing `std::invalid_argument` for an unknown name or parameter or a value outside its bounds.
- `c_luminosity_from_binary(stream, force)`: reconstructs from a binary record.
- `c_luminosity_canonical_name(name)` and `c_luminosity_model_names()`: name lookups.

The solar anchors come from `TidalPyConstants::d_MASS_SOLAR` and `d_LUMINOSITY_SOLAR`, and the Stefan-Boltzmann constant from the shared config singleton (`tidalpy_config_ptr->d_SBC`).

## Adding a New Model

**C++ (`TidalPy/Stellar/luminosity_.hpp`)**

To add a luminosity model named `Foo`:

1. Add a free function implementing the relation, returning NaN for a non-positive mass.
2. Add `c_FooLuminosity : public c_SpecModel<c_FooLuminosity, c_LuminosityBase>`: its `parameter_specs()` table (argument name, config key, member, default, bounds, and a one-line description), `C_CLASS_ID`, two constructors that call `p_initialize`, and `calc_luminosity` (calling the free function). Override `p_validate` for checks across parameters and `p_update_derived` for cached values.
3. Add one row to `c_luminosity_registry()` with the model's names and aliases.

**C++ (`TidalPy/Utilities/binary/binary_.hpp`)**

4. Reserve a unique `BinaryClassID::FooLuminosity` in the 100X block.

**Cython (`luminosity.pxd` and `luminosity.pyx`)**

5. Declare `cdef class FooLuminosity(LuminosityBase)` in the `.pxd`. In the `.pyx`, add the class with a docstring and `MODEL_NAME = "foo"`, include it in the `ModelFamily` list, and add a lower-case `foo(mass, ...)` convenience function. Add its luminosity tests to `Tests/Test_Stellar/`. The generic tests in `Tests/Test_Utilities/Test_Classes/test_spec_models_01.py` cover its parameters, config, binary record, and errors without changes.

**Package, tests, and documentation**

6. Export `FooLuminosity` and `foo` from `TidalPy/Stellar/__init__.py`.
7. Extend `Tests/Test_Stellar/Test_Luminosity/test_luminosity_01.py`: luminosity against an independent reference, the factory and aliases, the temperature conversions, the config dict, and the binary round trip.
8. Document the model here with its relation, parameters, and references.

## References

- Cuntz, M., and Wang, Z. (2018). The mass-luminosity relation for a refined set of late-K/M stars. *Research Notes of the AAS*, 2(1), 19, [doi:10.3847/2515-5172/aaaa67](https://doi.org/10.3847/2515-5172/aaaa67). The low-mass polynomial exponent.
- Salaris, M., and Cassisi, S. (2005). *Evolution of Stars and Stellar Populations*. The piecewise main-sequence regimes.
