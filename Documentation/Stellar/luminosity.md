# Luminosity Models (`Stellar`)

_Updated: 2026-10-02_

A luminosity model maps a star's mass onto its luminosity $L$ \[W\]. That sets the effective temperature through the Stefan-Boltzmann law and, once the star is in a `System`, the flux and equilibrium temperature of every world orbiting it. All models share the conversions between luminosity and effective temperature, which depend only on the radius:

$$L = 4 \pi R^2 \sigma T^4, \qquad T = \left(\frac{L}{4 \pi R^2 \sigma}\right)^{1/4}$$

The models differ only in $L(M)$.

## Example

```python
import numpy as np
from TidalPy.Stellar import (
    FixedLuminosity, MassToLuminosity, PowerLawLuminosity,
    make_luminosity, fixed, mass_to_luminosity, power_law)
from TidalPy.constants import mass_solar, radius_solar

model = MassToLuminosity()
luminosity = model.calc_luminosity(mass_solar)                          # [W]
temperature = model.calc_effective_temperature(mass_solar, radius_solar)  # [K]

# The two conversions that do not involve the mass relation
model.calc_luminosity_from_temperature(5772.0, radius_solar)   # [W]
model.calc_temperature_from_luminosity(luminosity, radius_solar)  # [K]

# An array of masses returns a float64 array
masses = np.array([0.1, 0.5, 1.0, 5.0, 30.0]) * mass_solar
curve = model.calc_luminosity(masses)

# Name factory: case-insensitive, aliases accepted
model = make_luminosity("cw")
model = make_luminosity("power_law", {"power_law_coeff": 1.0, "power_law_exponent": 4.0})
model = make_luminosity("constant", {"luminosity_w": 3.828e26})

# Build, evaluate, and discard in one call
luminosity = mass_to_luminosity(mass_solar)
luminosity = power_law(mass_solar, coeff=1.0, exponent=4.0)
luminosity = fixed(mass_solar, luminosity=3.828e26)
```

The convenience functions `fixed(mass, luminosity=0.0)`, `mass_to_luminosity(mass)`, and `power_law(mass, coeff=1.0, exponent=3.5)` take a float or an array of masses; their model parameters are constants.

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

`set_luminosity_model` shares the model rather than copying it, so one model can serve several stars. The mass-derived calculations raise `RuntimeError` with no model attached; `set_luminosity` and `set_effective_temperature` work without one and keep the two scalars consistent. In a `System`, `calc_insolation_flux(world)` and `calc_equilibrium_temperature(world)` give each world its incident flux and gray-body equilibrium temperature.

## Models

| Model (aliases) | Relation |
|---|---|
| `fixed` (`constant`) | $L$ is supplied directly and does not depend on mass. |
| `mass_to_luminosity` (`cuntz_wang`, `cw`) | The piecewise main-sequence relation below. |
| `power_law` (`powerlaw`) | $L = L_\odot\, c \left(M / M_\odot\right)^{p}$, with $c = 1$ and $p = 3.5$ by default. |

The Python classes are `FixedLuminosity`, `MassToLuminosity`, and `PowerLawLuminosity`; a config dict holds the first (canonical) name.

### Fixed

Reports its given luminosity at any mass. Use it when the luminosity is measured or is the independent variable of a sweep.

### Mass to Luminosity

A piecewise function of the mass, with one power law per mass range. Below, $x = M / M_\odot$.

| Range | Relation |
|---|---|
| $x < 0.2$ | $L = 0.23\, L_\odot\, x^{2.3}$ |
| $0.2 \le x < 0.85$ | $L = L_\odot\, x^{\,p(x)}$, with the Cuntz and Wang (2018) polynomial exponent $p(x) = -141.7 x^4 + 232.4 x^3 - 129.1 x^2 + 33.29 x + 0.215$ |
| $0.85 \le x < 2$ | $L = L_\odot\, x^{4}$ |
| $2 \le x < 55$ | $L = 1.4\, L_\odot\, x^{3.5}$ |
| $x \ge 55$ | $L = 3.2 \times 10^{4}\, L_\odot\, x$ |

The polynomial exponent of the second branch fits the whole M-dwarf range. The branches are the published fits as they stand, and not every joint is continuous:

| Joint | Below | Above | Step |
|---|---|---|---|
| $x = 0.2$ | $5.68 \times 10^{-3}\, L_\odot$ | $4.62 \times 10^{-3}\, L_\odot$ | $-18.7$ percent |
| $x = 0.85$ | | | continuous |
| $x = 2$ | $16.0\, L_\odot$ | $15.8\, L_\odot$ | $-1.0$ percent |
| $x = 55$ | $1.73 \times 10^{6}\, L_\odot$ | $1.76 \times 10^{6}\, L_\odot$ | $+1.9$ percent |

At $x = 0.2$ the luminosity falls as the mass rises, so a star evolved through it (or a sweep across it) sees a step. Keep masses on one side of a joint, or blend the branches, when that matters.

### Power Law

A single power law with a caller-set prefactor and exponent, to reproduce a paper that adopted one or to vary the mass-luminosity slope $p$ directly. The defaults give the textbook $L \propto M^{3.5}$ scaling for solar-type stars.

### Parameters

Constructors take parameters by argument name or config key, as keywords, positionally (in `get_parameter_info()` order), or as a `config=` table: `PowerLawLuminosity(exponent=4.0)` equals `PowerLawLuminosity(power_law_exponent=4.0)`. In `make_luminosity(model_name, config=None)`, absent keys take the defaults. An unknown name, a key the model does not read (even one another luminosity model reads), or an out-of-bounds value raises `ValueError` naming the closest accepted one.

| Parameter | Config key | Default | Bounds | Model |
|---|---|---|---|---|
| `luminosity` | `luminosity_w` | 0.0 | non-negative | Fixed: the luminosity reported at every mass \[W\]. |
| `coeff` | `power_law_coeff` | 1.0 | positive | Power law: the prefactor $c$, in solar luminosities. |
| `exponent` | `power_law_exponent` | 3.5 | finite | Power law: the exponent $p$ of the mass ratio. |

`mass_to_luminosity` has no parameters. Parameters read as attributes (`model.exponent`), and every model has `model_name`, `parameters`, `get_parameter(name)`, `get_parameter_info()`, and `with_parameters(**changes)`. `luminosity_model_names()` lists the canonical names and `luminosity_config_keys(name)` a model's keys.

### Limits

A non-positive mass returns NaN rather than raising, as does a non-positive luminosity or radius in the temperature conversions, so a bad entry in a large sweep shows as NaN instead of aborting it. The Stefan-Boltzmann constant comes from the shared TidalPy config.

## Serialization

`get_config_dict()` returns `model` plus the parameters by config key, which `make_luminosity` accepts back; `save_config(path)` writes it as TOML. `save_binary(path)` and `load_binary(path)` write and read the binary record. A star's luminosity model is saved and restored with the star's binary record.

## C++ API

The models are in `TidalPy/Stellar/luminosity_.hpp` (namespace `tidalpy`, header only).

```cpp
#include "luminosity_.hpp"  // pulls in luminosity_base_.hpp

using namespace tidalpy;

c_ParamMap params;  // Config key to values; one value for a scalar
params["power_law_exponent"] = {4.0};

std::unique_ptr<c_LuminosityBase> model = c_find_luminosity("power_law", params);
const double luminosity = model->calc_luminosity(mass);  // [W]
const double temperature = model->calc_effective_temperature(mass, radius);  // [K]
```

- `c_LuminosityBase : c_PhysicsBase`: abstract, with `calc_luminosity(mass)` pure virtual, the Stefan-Boltzmann conversions, and `calc_luminosity_vectorize_mass`.
- `c_FixedLuminosity`, `c_MassToLuminosity`, and `c_PowerLawLuminosity`. The two relations are also free functions, `lum_from_mass(mass)` and `lum_from_power_law(mass, coeff, exponent)`.
- `c_find_luminosity(name, params)` builds a model as a `unique_ptr`, throwing `std::invalid_argument` for an unknown name or parameter or an out-of-bounds value. `c_luminosity_from_binary(stream, force)` rebuilds one from a binary record. `c_luminosity_canonical_name(name)` and `c_luminosity_model_names()` look up names.

## Adding a New Model

To add a model named `Foo`:

1. In `TidalPy/Stellar/luminosity_.hpp`, add a free function for the relation (NaN for a non-positive mass) and `c_FooLuminosity : public c_SpecModel<c_FooLuminosity, c_LuminosityBase>` with its `parameter_specs()` table, `C_CLASS_ID`, two constructors that call `p_initialize`, and `calc_luminosity` (override `p_validate` and `p_update_derived` as needed). Add a row to `c_luminosity_registry()`.
2. Reserve `BinaryClassID::FooLuminosity` in the 100X block of `TidalPy/Utilities/binary/binary_.hpp`.
3. Declare `cdef class FooLuminosity(LuminosityBase)` in `luminosity.pxd`; in `luminosity.pyx` add the class (docstring, `MODEL_NAME = "foo"`), add it to the `ModelFamily` list, and add a `foo(mass, ...)` function. Export both from `TidalPy/Stellar/__init__.py`.
4. Extend `Tests/Test_Stellar/Test_Luminosity/test_luminosity_01.py` (an independent reference, aliases, temperature conversions, config, binary round trip; the generic `test_spec_models_01.py` covers the rest) and document the model here.

## References

- Cuntz, M., and Wang, Z. (2018). The mass-luminosity relation for a refined set of late-K/M stars. *Research Notes of the AAS*, 2(1), 19, [doi:10.3847/2515-5172/aaaa67](https://doi.org/10.3847/2515-5172/aaaa67). The low-mass polynomial exponent.
- Salaris, M., and Cassisi, S. (2005). *Evolution of Stars and Stellar Populations*. The piecewise main-sequence regimes.
