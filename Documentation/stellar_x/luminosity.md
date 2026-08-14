# Stellar Luminosity (`stellar_x`)

The `stellar_x` module holds the stellar physics models. The luminosity hierarchy turns a star's mass
into its luminosity and effective temperature. All quantities are MKS: mass [kg], radius [m],
temperature [K], luminosity [W].

Every model shares the Stefan-Boltzmann conversions between a star's effective surface temperature and
its luminosity (both need the stellar radius) and implements a mass-to-luminosity relation:

* `L = 4 pi R^2 sigma T^4`
* `T = (L / (4 pi R^2 sigma))^(1/4)`
* `L = f(M)` (model specific)

## Models

| Name (aliases) | Relation |
|----------------|----------|
| `fixed` (`constant`) | Luminosity supplied directly, independent of mass. |
| `mass_to_luminosity` (`cuntz_wang`, `cw`) | Piecewise main-sequence `L(M)`: `0.23 (M/Msun)^2.3` below 0.2 Msun, the Cuntz and Wang (2018) polynomial exponent up to 0.85 Msun, `(M/Msun)^4` up to 2 Msun, `1.4 (M/Msun)^3.5` up to 20 Msun, and `3.2e4 (M/Msun)` above. |
| `power_law` (`powerlaw`) | Single power law `L = Lsun * coeff * (M/Msun)^exponent` (defaults `coeff=1`, `exponent=3.5`). |

The solar anchors `Msun` and `Lsun` come from the TidalPy constants (`d_MASS_SOLAR`,
`d_LUMINOSITY_SOLAR`). A non-positive mass returns NaN.

## Python API

```python
from TidalPy.stellar_x import make_luminosity, MassToLuminosity, mass_to_luminosity

model = make_luminosity("mass_to_luminosity")     # or MassToLuminosity()
lum = model.calc_luminosity(mass_kg)              # float or ndarray (mass may be an array)
temp = model.calc_effective_temperature(mass_kg, radius_m)
lum_from_t = model.calc_luminosity_from_temperature(temperature_k, radius_m)
temp_from_l = model.calc_temperature_from_luminosity(luminosity_w, radius_m)

# Direct convenience functions (build, solve, discard):
lum = mass_to_luminosity(mass_kg)                 # float or ndarray
```

`make_luminosity(name, config)` maps the (case-insensitive) name/alias to the C++ enum factory and
returns the matching rich subclass. Config keys: `luminosity_w` (fixed), `power_law_coeff` /
`power_law_exponent` (power law). Each model round-trips through `save_binary` / `load_binary` and
reports its parameters via `get_config_dict()`.

### On a star

`StarWorld` can hold a luminosity model and derive its luminosity and effective temperature from its
own mass and radius:

```python
from TidalPy.structures_x.worlds.stellar import StarWorld
from TidalPy.stellar_x import MassToLuminosity

star = StarWorld("sun", radius_m, mass_kg)
star.set_luminosity_model(MassToLuminosity())     # transfers ownership of the model
lum = star.calc_luminosity_from_mass()            # [W]
temp = star.calc_effective_temperature_from_mass()# [K]
star.update_luminosity_from_mass()                # writes L and T onto the star's scalar fields
```

`set_luminosity_model` moves ownership of the C++ model into the star; the passed wrapper becomes an
empty shell. The mass-derived calculations raise `RuntimeError` if no model has been attached. The
model is optional: a star still keeps a consistent scalar temperature and luminosity (via
Stefan-Boltzmann) without one.

## C++ API

```cpp
#include "luminosity_.hpp"   // pulls in luminosity_base_.hpp

using namespace tidalpy;

c_LuminosityConfig config;
config.power_law_exponent = 4.0;
std::unique_ptr<c_LuminosityBase> model = c_find_luminosity("power_law", config);

const double lum  = model->calc_luminosity(mass_kg);
const double temp = model->calc_effective_temperature(mass_kg, radius_m);
```

* `c_LuminosityBase : c_PhysicsBase` — abstract; pure virtual `calc_luminosity(mass)`, plus the shared
  Stefan-Boltzmann conversions and `calc_luminosity_vectorize_mass`.
* Concrete models `c_FixedLuminosity`, `c_MassToLuminosity`, `c_PowerLawLuminosity`.
* `enum class c_LuminosityModel { Fixed, MassToLuminosity, PowerLaw }`.
* `c_luminosity_model_from_name(name)` — name/alias to enum (throws `std::invalid_argument` on unknown).
* `c_find_luminosity(model_or_name, config)` — heap-allocates the model as a `unique_ptr`.
* `c_luminosity_from_binary(stream, force)` — reconstruct from a binary record (class ids 1001-1003).

The Stefan-Boltzmann constant is read from the shared TidalPy config (`tidalpy_config_ptr->d_SBC`); a
null config pointer or a non-positive/degenerate input yields NaN.

## Adding a new model

1. Add a concrete `c_<Name>Luminosity : c_LuminosityBase` in `luminosity_.hpp` overriding
   `calc_luminosity(mass)` and the binary `write_binary`/`read_binary` (via the `c_PhysicsBase`
   helpers). Register a unique `BinaryClassID`.
2. Add an enum value to `c_LuminosityModel`, a name/alias branch to `c_luminosity_model_from_name`, and
   switch cases to `c_find_luminosity` and `c_luminosity_from_binary`.
3. Add a Cython wrapper subclass in `luminosity.pyx`/`.pxd`, wire it into `make_luminosity`, and export
   it from `__init__.py`/`__init__.pxd`.
4. Add tests and update this page and the changelog.

## References

* Cuntz and Wang (2018), doi:10.3847/2515-5172/aaaa67 — low-mass mass-luminosity polynomial exponent.
* Wikipedia mass-luminosity relation — the high/low-mass piecewise regimes.
