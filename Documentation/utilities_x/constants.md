# Constants (`TidalPy.constants`)

_Updated: 2026-09-12_

TidalPy's constants come in three kinds, and the distinction matters because they are set at different times and from different sources. Mathematical and floating-point limits are fixed at compile time. Physical constants are pulled from SciPy when the package initializes, so TidalPy always agrees with the reference values SciPy ships. Numerical floors and ceilings, the values that keep a solver from dividing by a vanishing viscosity or evaluating a mode at zero frequency, come from the configuration file and can be changed by the user.

All of them are stored in one place at the C++ level: a struct of static members in `constants_.hpp`, reachable by every compiled module through the shared pointer `tidalpy_config_ptr`. That arrangement is what lets a value set from Python reach code running inside a `nogil` inner loop without a lookup.

## Using them from Python

```python
import TidalPy.constants as constants

constants.G                 # [N m2 kg-2] Newton's constant, from SciPy
constants.pi
constants.mass_solar        # [kg]   also M_sol
constants.radius_earth      # [m]    also R_earth
constants.luminosity_solar  # [W]    also L_sol
constants.year              # [s]    Julian year
constants.min_viscosity     # [Pa s] configurable floor
```

Most constants have a short alias beside the descriptive name (`M_sol` for `mass_solar`, `Au` for `au`, `SBC` for `sbc`, `k` for `k_boltzmann`), because both spellings appear throughout the literature and the older code. They refer to the same value.

The corresponding C++ names carry a `d_` prefix and full capitals: `d_MASS_SOLAR`, `d_LUMINOSITY_SOLAR`, `d_PI`, `d_NAN`, `d_EPS`. A physics module reads them through `TidalPyConstants::d_MASS_SOLAR` for the compile-time values and through `tidalpy_config_ptr->d_G` for the runtime ones.

## What is fixed and what is not

**Compile-time, read-only.** Mathematical constants ($\pi$, infinity, NaN) and floating-point limits (the largest and smallest normal double, the machine epsilon, the mantissa digit count) are `constexpr`. So are the solar-system body properties: the masses and radii of the Sun, Earth, Jupiter, Pluto, and Io, and the solar luminosity, all set to the IAU nominal values.

**Set at initialization, from SciPy.** The gravitational constant, the astronomical unit, the Stefan-Boltzmann constant, the molar gas constant, and the Boltzmann constant are read from `scipy.constants` each time TidalPy initializes. The Julian year comes from the same place. This is deliberate: a physical constant should have one authoritative source, and duplicating CODATA values into a header is how packages end up disagreeing with each other in the sixth digit.

**Set from configuration.** The numerical guards are read from the configuration file: the minimum and maximum tidal frequency, the minimum spin-orbit frequency difference, the minimum viscosity and modulus, and the minimum layer thickness. These are not physics. They are the thresholds below which a quantity is treated as zero or a mode is dropped, and their right values depend on the problem, so they are exposed rather than hard-coded.

## Updating them

Two functions repopulate the shared C++ struct, and both run automatically during package initialization.

`update_constants()` reads the classic configuration and SciPy. It sets the physical constants and the numerical guards from the legacy config sections.

`update_constants_x()` then copies the `[numerical]` section of the new configuration, `TidalPy.config_x` loaded from `TidalPy_Configs_x.toml`, over the shared numerical fields. It runs second, so for the frequency, viscosity, modulus, and thickness floors the new configuration wins. It does not touch the physical constants, which stay as SciPy set them.

There is one process-wide C++ config singleton shared by both code paths, which is why the ordering is what decides the result. Call either function again after editing a configuration in a running session; each call reconfigures in place.

## Where the values live

The authoritative lists are short enough to read directly, and they change often enough that copying them into this page would guarantee it goes stale:

- [`constants_.hpp`](https://github.com/jrenaud90/TidalPy/blob/main/TidalPy/constants_.hpp) for the compile-time values and the runtime struct layout.
- [`constants.pyx`](https://github.com/jrenaud90/TidalPy/blob/main/TidalPy/constants.pyx) for the Python names and their aliases.
- `defaultc_x.py` for the default values of everything the configuration file controls. See the [TOML schema](../structures_x/config/toml_schema.md).
