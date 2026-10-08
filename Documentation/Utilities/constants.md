# Constants (`TidalPy.constants`)

_Updated: 2026-09-30_

TidalPy's constants come in three kinds. Mathematical and floating-point limits are fixed at compile time. Physical constants are read from SciPy (mostly) when the package initializes, so TidalPy agrees with its dependencies. Numerical floors and ceilings, which keep a solver from dividing by a vanishing modulus or evaluating a mode at zero frequency, come from the configuration file and can be changed. A value set from Python reaches every compiled module.

## Python Usage

```python
import TidalPy.constants as constants

constants.G                 # [N m2 kg-2] Newton's constant, from SciPy
constants.pi
constants.mass_solar        # [kg]   also M_sol
constants.radius_earth      # [m]    also R_earth
constants.luminosity_solar  # [W]    also L_sol
constants.year              # [s]    Julian year
constants.seconds_per_myr   # [s]    one Julian mega-year, exact
constants.min_frequency     # [rad s-1] configurable floor
```

Most constants also have a short alias (`M_sol` for `mass_solar`, `Au` for `au`, `SBC` for `sbc`, `k` for `k_boltzmann`) with the same value.

In C++ the names carry a `d_` prefix and full capitals (`d_MASS_SOLAR`, `d_PI`, `d_NAN`, `d_EPS`). Compile-time values are read as `TidalPyConstants::d_MASS_SOLAR` and runtime ones through `tidalpy_config_ptr->d_G`.

## Compile-Time and Runtime Values

### Compile-Time Values

Mathematical constants ($\pi$, infinity, NaN), floating-point limits (the largest and smallest normal double, the machine epsilon, the mantissa digit count), the masses and radii of the Sun, Earth, Jupiter, Pluto, and Io, and the solar luminosity (IAU nominal values) are fixed at compile time. So is the number of seconds in a Julian mega-year, which is exact.

### Values from SciPy

The gravitational constant, the astronomical unit, the Stefan-Boltzmann constant, the molar gas constant, and the Boltzmann constant are read from `scipy.constants` each time TidalPy initializes. The Julian year comes from the same place.

### Values from Configuration

The numerical guards (the frequency extremes, the minimum modulus and solid rigidity, the minimum layer thickness, the numerical floor, and the rest) and the `[eos_solver]` and `[radial_solver]` defaults are read from the configuration file (see the [configuration page](../Overview/2_TidalPy_Configurations.md#numerical-settings)). They are the thresholds below which a quantity counts as zero or a mode is dropped, and the right values depend on the problem.

The numerical floor is not a physical threshold. It is the smallest magnitude a denominator may take before a guard substitutes it, read by `Rheology`, `Cooling`, and `Radiogenics` where a zero forcing frequency, layer thickness, or half life would otherwise divide to infinity. Its default, `1e-100`, only replaces a true zero in practice.

## Updating

`update_constants()` copies the `[numerical]`, `[eos_solver]`, and `[radial_solver]` sections of `TidalPy.config` into the C++ code, rereads the physical constants from SciPy, and refreshes the Python names in `TidalPy.constants`. It runs when the package initializes and whenever `TidalPy.reinit` changes the configuration. Call it yourself after editing `TidalPy.config` directly in a running session.

## Files

[`constants_.hpp`](https://github.com/jrenaud90/TidalPy/blob/main/TidalPy/constants_.hpp) holds the compile-time values, [`constants.pyx`](https://github.com/jrenaud90/TidalPy/blob/main/TidalPy/constants.pyx) the Python names and aliases, and [`defaultc.py`](https://github.com/jrenaud90/TidalPy/blob/main/TidalPy/defaultc.py) the configuration defaults.
