# Luminosity Models (`stellar_x`)

_Updated: 2026-09-12_

A luminosity model maps a star's mass onto its luminosity $L$ [W]. That single relation, plus the star's radius, fixes everything else this module reports: the effective temperature through the Stefan-Boltzmann law, and, once the star is placed in a `System`, the flux and equilibrium temperature of every world orbiting it.

All three models share the same conversions between luminosity and effective temperature, which depend only on the radius and not on the mass relation:

$$L = 4 \pi R^2 \sigma T^4, \qquad T = \left(\frac{L}{4 \pi R^2 \sigma}\right)^{1/4}$$

What separates the models is $L(M)$.

## Inheritance

```
c_TidalPyBaseClass
  └── c_PhysicsBase
        └── c_LuminosityBase           (abstract)
              ├── c_FixedLuminosity      alias "constant"
              ├── c_MassToLuminosity     aliases "cuntz_wang", "cw"
              └── c_PowerLawLuminosity   alias "powerlaw"
```

The base class declares `calc_luminosity(mass)` pure virtual and supplies the Stefan-Boltzmann conversions and the vectorized mass sweep, so a new relation is a single method. The Cython classes mirror the hierarchy: `LuminosityBase`, `FixedLuminosity`, `MassToLuminosity`, `PowerLawLuminosity`.

## The three models

| Model (aliases) | Relation |
|---|---|
| `fixed` (`constant`) | $L$ is supplied directly and does not depend on mass. |
| `mass_to_luminosity` (`cuntz_wang`, `cw`) | The piecewise main-sequence relation below. |
| `power_law` (`powerlaw`) | $L = L_\odot\, c \left(M / M_\odot\right)^{p}$, with $c = 1$ and $p = 3.5$ by default. |

### Fixed

Reports the luminosity it was given for any mass. Use it when the star's luminosity is a measured quantity rather than something to be predicted, which is the common case for a named host star, and when a sweep needs luminosity as an independent variable.

### Mass to luminosity

The main-sequence relation is not a single power law. The exponent of $L \propto M^{p}$ falls from roughly 4 for solar-type stars toward 2.3 at the bottom of the main sequence and toward 1 at the top, because the dominant energy transport and opacity regimes change across that range. The model is therefore piecewise in $x = M / M_\odot$:

| Range | Relation |
|---|---|
| $x < 0.2$ | $L = 0.23\, L_\odot\, x^{2.3}$ |
| $0.2 \le x < 0.85$ | $L = L_\odot\, x^{\,p(x)}$, with the Cuntz and Wang (2018) polynomial exponent $p(x) = -141.7 x^4 + 232.4 x^3 - 129.1 x^2 + 33.29 x + 0.215$ |
| $0.85 \le x < 2$ | $L = L_\odot\, x^{4}$ |
| $2 \le x < 55$ | $L = 1.4\, L_\odot\, x^{3.5}$ |
| $x \ge 55$ | $L = 3.2 \times 10^{4}\, L_\odot\, x$ |

The mass-ratio exponent in the second branch is what makes the low-mass end usable: a single power law fit across the whole M-dwarf range misses the observed luminosities badly, and the polynomial exponent absorbs the curvature. The last two branches meet at $x = 55$, where both give about $1.75 \times 10^{6}\, L_\odot$, so the relation is continuous there.

This is the model to use when the star's mass is what is known, and the one to use for any population study where the host mass is being varied.

### Power law

A single power law with a caller-set prefactor and exponent. Use it to reproduce a paper that adopted one, or to isolate the effect of the mass-luminosity slope by varying $p$ directly. The defaults reproduce the classic $L \propto M^{3.5}$ scaling for solar-type stars.

### Behavior at the limits

A non-positive mass returns NaN rather than raising, and so does a non-positive luminosity or radius in the temperature conversions. The Stefan-Boltzmann constant comes from the shared TidalPy config, so a degenerate configuration also yields NaN. This keeps a bad entry in a large mass sweep visible as NaN in the output rather than aborting the sweep.

## Choosing a model

Use `mass_to_luminosity` whenever mass is the input and the star is on the main sequence. It is the default choice and needs no parameters.

Use `fixed` when the luminosity is measured, tabulated, or being scanned as a free parameter.

Use `power_law` when matching a published relation, or when the exponent itself is the subject of the study.

## Python API

```python
import numpy as np
from TidalPy.stellar_x import (
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

`make_luminosity(model_name, config=None)` maps the name or alias to the C++ enum factory and returns the matching subclass. Its config keys are the ones `get_config_dict()` emits, which are not always the constructor's argument names.

| Config key | Model | Constructor argument | Meaning |
|---|---|---|---|
| `luminosity_w` | fixed | `luminosity` | Luminosity to report [W]. |
| `power_law_coeff` | power law | `coeff` | Dimensionless prefactor. |
| `power_law_exponent` | power law | `exponent` | Dimensionless exponent. |

`mass_to_luminosity` takes no parameters. A config key that no luminosity model reads raises `ValueError` naming the closest accepted key, so a misspelling fails loudly instead of silently building a default model.

The convenience functions `fixed(mass, luminosity=0.0)`, `mass_to_luminosity(mass)`, and `power_law(mass, coeff=1.0, exponent=3.5)` each build a stack-allocated C++ model, evaluate it, and discard it. The mass may be a float or an array; the model parameters are always constants.

### On a star

```python
from TidalPy.structures_x.worlds.stellar import StarWorld
from TidalPy.stellar_x import MassToLuminosity
from TidalPy.constants import mass_solar, radius_solar

star = StarWorld("sun", radius_solar, mass_solar)
star.set_luminosity_model(MassToLuminosity())   # transfers ownership of the model
star.luminosity_model_set                       # True

star.calc_luminosity_from_mass()                # [W]
star.calc_effective_temperature_from_mass()     # [K]
star.update_luminosity_from_mass()              # writes both onto the star's own fields
star.luminosity, star.effective_temperature
```

`set_luminosity_model` moves ownership of the C++ model into the star, leaving the passed wrapper an empty shell. The two mass-derived calculations raise `RuntimeError` when no model has been attached; `set_luminosity` and `set_effective_temperature` remain available and keep the star's two scalars consistent through the Stefan-Boltzmann relation without one.

Once the star is part of a `System`, `calc_insolation_flux(world)` and `calc_equilibrium_temperature(world)` use its luminosity to give each orbiting world its incident flux and gray-body equilibrium temperature.

## Serialization

Every model supports the standard interfaces inherited from the base class.

- `get_config_dict()` returns the model name under the key `model` plus its parameters. The dict is accepted by `make_luminosity`, so a model round-trips through it.
- `save_config(path)` writes the same content as TOML.
- `save_binary(path)` and `load_binary(path, force=False)` use the TidalPy binary format. All three models serialize through the shared `c_PhysicsBase` scalar helpers, and the layer back-pointer is not written.

## C++ API

```cpp
#include "luminosity_.hpp"   // pulls in luminosity_base_.hpp

using namespace tidalpy;

c_LuminosityConfig config;
config.power_law_exponent = 4.0;

const c_LuminosityModel model_id = c_luminosity_model_from_name("power_law");
std::unique_ptr<c_LuminosityBase> model = c_find_luminosity(model_id, config);

const double luminosity  = model->calc_luminosity(mass);                    // [W]
const double temperature = model->calc_effective_temperature(mass, radius); // [K]
```

- `c_LuminosityBase : c_PhysicsBase`: abstract, with `calc_luminosity(mass)` pure virtual plus the shared Stefan-Boltzmann conversions and `calc_luminosity_vectorize_mass`.
- Concrete models `c_FixedLuminosity`, `c_MassToLuminosity`, and `c_PowerLawLuminosity`.
- `enum class c_LuminosityModel { Fixed, MassToLuminosity, PowerLaw }`.
- `c_luminosity_model_from_name(name)`: name or alias to enum, throwing `std::invalid_argument` on an unknown name.
- `c_find_luminosity(model, config)`: heap-allocates the model as a `unique_ptr`.
- `c_luminosity_from_binary(stream, force)`: reconstructs from a binary record.

The solar anchors come from `TidalPyConstants::d_MASS_SOLAR` and `d_LUMINOSITY_SOLAR`, and the Stefan-Boltzmann constant from the shared config singleton (`tidalpy_config_ptr->d_SBC`). Binary class ids 1000 through 1003 are reserved for this module.

## Adding a new model

**C++ (`TidalPy/stellar_x/luminosity_.hpp`)**

1. Add any new parameters to `c_LuminosityConfig` with sensible defaults. The single combined config is shared by all models.
2. Add a free function implementing the relation, returning NaN for a non-positive mass.
3. Add the model class deriving from `c_LuminosityBase`: a default constructor and one taking the config, `get_*` accessors, the `calc_luminosity` override, and `write_binary` / `read_binary` through the `c_PhysicsBase` helpers.
4. Add the enum value, the name and alias branch in `c_luminosity_model_from_name`, and the cases in `c_find_luminosity` and `c_luminosity_from_binary`.

**C++ (`TidalPy/Utilities_x/binary_x/binary_.hpp`)**

5. Add a unique `BinaryClassID` in the 100X block.

**Cython (`luminosity.pxd` and `luminosity.pyx`)**

6. Declare the C++ class and the new enum value in the `.pxd`.
7. Add the `cdef class` wrapper with parameter properties, the adoption branch in `make_luminosity`, and the lower-case convenience function.

**Package, tests, and docs**

8. Export the class and the function from `__init__.py`, and the C++ names from `__init__.pxd`.
9. Extend `Tests/Test_Stellar_x`: model name, luminosity against an independent reference, factory and aliases, the temperature conversions, config dict, and binary round trip.
10. Document the model here and add a changelog entry.

## References

- Cuntz, M., and Wang, Z. (2018). The mass-luminosity relation for a refined set of late-K/M stars. *Research Notes of the AAS*, 2(1), 19, [doi:10.3847/2515-5172/aaaa67](https://doi.org/10.3847/2515-5172/aaaa67). The low-mass polynomial exponent.
- Salaris, M., and Cassisi, S. (2005). *Evolution of Stars and Stellar Populations*. The piecewise main-sequence regimes.
