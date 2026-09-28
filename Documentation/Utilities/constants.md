# Constants (`TidalPy.constants`)

_Updated: 2026-09-25_

TidalPy's constants come in three kinds, set at different times and from different sources. Mathematical and floating-point limits are fixed at compile time. Physical constants are pulled from third-party sources (mostly SciPy) when the package initializes, so TidalPy agrees with the reference values of those dependencies. Numerical floors and ceilings, the values that keep a solver from dividing by a vanishing modulus or evaluating a mode at zero frequency, come from the configuration file and can be changed by the user.

All of them are stored in one place at the C++ level: a struct of static members in `constants_.hpp`, reachable by every compiled module through the shared pointer `tidalpy_config_ptr`, so a value set from Python reaches code running inside a `nogil` inner loop without a lookup.

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

Most constants have a short alias beside the descriptive name (`M_sol` for `mass_solar`, `Au` for `au`, `SBC` for `sbc`, `k` for `k_boltzmann`), because both spellings appear throughout the literature. They refer to the same value.

The corresponding C++ names carry a `d_` prefix (indicating they are doubles) and full capitals: `d_MASS_SOLAR`, `d_LUMINOSITY_SOLAR`, `d_PI`, `d_NAN`, `d_EPS`. A physics module reads them through `TidalPyConstants::d_MASS_SOLAR` for the compile-time values and through `tidalpy_config_ptr->d_G` for the runtime ones.

## Compile-Time and Runtime Values

### Compile-Time Values

Mathematical constants ($\pi$, infinity, NaN) and floating-point limits (the largest and smallest normal double, the machine epsilon, the mantissa digit count) are `constexpr`. So are the solar-system body properties: the masses and radii of the Sun, Earth, Jupiter, Pluto, and Io, and the solar luminosity, all set to the IAU nominal values. The number of seconds in a Julian mega-year is compile-time as well, because the Julian year is exact by definition.

### Values from SciPy

The gravitational constant, the astronomical unit, the Stefan-Boltzmann constant, the molar gas constant, and the Boltzmann constant are read from `scipy.constants` each time TidalPy initializes. The Julian year comes from the same place.

### Values from Configuration

The numerical guards are read from the `[numerical]` section of the configuration file: the minimum and maximum tidal frequency, the minimum modulus and solid rigidity, the minimum layer thickness, the shared numerical floor, the layer-boundary continuity tolerance, and the other settings listed on the [configuration page](../Overview/2_TidalPy_Configurations.md#numerical-settings). The `[eos_solver]` and `[radial_solver]` defaults are stored in the same struct. These are the thresholds below which a quantity is treated as zero or a mode is dropped, and their right values depend on the problem, so they are exposed rather than hard-coded.

The numerical floor is not a physical threshold. It is the smallest magnitude a denominator may take before a guard substitutes it, and `Rheology`, `Cooling`, and `Radiogenics` all read it: a zero forcing frequency, a zero layer thickness, and a zero half life each reach a division that would otherwise produce infinity. Its default, `1e-100`, sits far below any physical value, so in practice it only replaces a true zero.

## Updating

`update_constants()` repopulates the shared C++ struct. It copies the `[numerical]`, `[eos_solver]`, and `[radial_solver]` sections of `TidalPy.config` (loaded from `TidalPy_Configs.toml`) into the struct, reads the physical constants from SciPy, and refreshes the Python names in `TidalPy.constants`.

It runs automatically when the package initializes and whenever `TidalPy.reinit` changes the configuration. Call it again after editing `TidalPy.config` directly in a running session; each call reconfigures in place.

## Files

Where the various variables are defined.

- [`constants_.hpp`](https://github.com/jrenaud90/TidalPy/blob/main/TidalPy/constants_.hpp) for the compile-time values and the runtime struct layout.
- [`constants.pyx`](https://github.com/jrenaud90/TidalPy/blob/main/TidalPy/constants.pyx) for the Python names and their aliases.
- [`defaultc.py`](https://github.com/jrenaud90/TidalPy/blob/main/TidalPy/defaultc.py) for the default values of everything the configuration file controls. See [TidalPy Configurations](../Overview/2_TidalPy_Configurations.md).
