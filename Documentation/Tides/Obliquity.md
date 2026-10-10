# Obliquity Functions

_Updated: 2026-10-06_

The obliquity functions $F_{l,m,p}(I)$ are the other forcing component of the tidal potential, alongside the [eccentricity functions](Eccentricity.md); see $F_{lmp}(I)$ in Eq. 1 of [Kaula (1964)](http://doi.wiley.com/10.1029/RG002i004p00661). TidalPy returns the unsquared $F_{l,m,p}(I)$, the form the potential uses; `obliquity_squared_func` returns the squares the global heating uses. Several $(l, m, p)$ modes stay non-zero at $I = 0$, so the functions cannot be skipped for an aligned body.

> [!NOTE]
> The wider literature usually calls these inclination functions. For tides the relevant angle is between the deformed body's spin axis and its orbital plane about its host, normally called the body's obliquity. Inclination, an orbital element, still sets the relative obliquity of a host planet with respect to its tidal target's orbit.

## Example

The usual route is a world's `[tides]` configuration (`obliquity_trunc_lvl` in TOML, which also accepts `"off"` and `"gen"`; see [Global Tides](global_tides.md) and the [TOML schema](../Structures/config/toml_schema.md)):

```python
world.set_tide_config(max_degree_l=2, eccentricity_truncation=6, obliquity_truncation=0)
```

Calling the functions directly:

```python
import numpy as np
from TidalPy.Tides.obliquity import obliquity_func, obliquity_squared_func

obliquity = np.radians(20.0)           # radians
modes_by_lmp, modes_by_lm = obliquity_func(obliquity, degree_l=2, truncation="gen")

print(len(modes_by_lmp))                     # 9 active modes
print(modes_by_lmp[(2, 2, 0)])               # one mode by its (l, m, p) key

for (l, m), by_p in modes_by_lm.items():     # walk the modes grouped by (l, m)
    for (p,), value in by_p:
        print(f"F({l},{m},{p}) = {value:0.3e}")

squares_by_lmp, _ = obliquity_squared_func(obliquity, degree_l=2, truncation=2)
print(squares_by_lmp[(2, 2, 0)])             # 9 - 9 I^2
```

Degrees $l = 2$ through $10$ are supported. The truncation is an integer level, `"off"`, or `"gen"` (the `obliquity_func` default). The two lookups hold the active modes: `modes_by_lmp` keyed by `(l, m, p)`, and `modes_by_lm` keyed by `(l, m)`, iterating as `((p,), value)` pairs. C++ and Cython callers have the same functions in `Tides/obliquity/obliquity_common`, `obliquity_accuracy_`, and `obliquity_driver`.

## Choosing a Truncation

These series converge at every obliquity, and the general functions cost little more than a table (at most two extra modes at degree 2). The truncated levels only save modes, and Love number solves, at small obliquity.

The accuracy columns give the largest obliquity up to which the degree-2 heating stays within the stated relative error of the general functions' heating: the worst case over the tide models and spin rates of the [eccentricity table](Eccentricity.md#choosing-a-truncation), at eccentricities from 0 to 0.2. Every limit is set by synchronous rotation at low eccentricity, where the obliquity tide carries the heating, with $-\mathrm{Im}\,k$ falling as one over the frequency (a stiff Maxwell body such as the Io of the truncation demo). Level 2 errs high there and level 4 low.

| Truncation | Aliases | Heating modes at $l=2$ | Error below $10^{-4}$ to | Error below 1% to | Warning above (max $l$ = 2 / 3 / 10) |
|---|---|---|---|---|---|
| `0` | `"off"` | 2 | $I = 0$ only | $I = 0$ only | any nonzero $I$ |
| `2` | | 4 | 0.010 rad (0.6°) | 0.11 rad (6.3°) | 0.35 / 0.24 / 0.075 rad (20° / 14° / 4.3°) |
| `4` | | 7 | 0.13 rad (7.4°) | 0.405 rad (23.2°) | 0.695 / 0.48 / 0.155 rad (40° / 27.5° / 8.9°) |
| `"gen"` | `"general"` | 9 | exact | exact | never |

At zero obliquity with truncation `0`, the two surviving degree-2 terms are $F_{2,0,1} = -1/2$ and $F_{2,2,0} = 3$.

A body spinning faster than near-synchronous (the spin bands of [By Spin Rate](Eccentricity.md#by-spin-rate)) holds both levels further: at degree 2, level 2 holds 10% to 0.58 rad (33°) between 1.5 and 5 times its mean motion and to 0.565 rad (32°) between 5 and 30 times it, and level 4 to 0.795 rad (46°) and 0.76 rad (44°). The near-synchronous band has the table's limits.

`recommend_obliquity_truncation(obliquity, tolerance=0.01, max_degree_l=2, spin_ratio=None)` returns the lowest level that holds a tolerance (0 at zero obliquity), or `"gen"`, and `obliquity_accuracy_limit(level, tolerance, max_degree_l, spin_ratio=None)` the largest obliquity a level holds a tolerance to. Tolerances, degrees, and `spin_ratio` work as in [Recommending a Level](Eccentricity.md#recommending-a-level).

```python
import numpy as np
from TidalPy.Tides.obliquity import recommend_obliquity_truncation

recommend_obliquity_truncation(np.radians(5.0))                    # 2: heating within 1% at 5 degrees
recommend_obliquity_truncation(np.radians(5.0), max_degree_l=3)    # 4: degree 3 needs more
recommend_obliquity_truncation(np.radians(8.0))                    # 4
recommend_obliquity_truncation(np.radians(8.0), max_degree_l=10)   # 'gen': degrees up to 10 need more still
recommend_obliquity_truncation(np.radians(23.4))                   # 'gen'
recommend_obliquity_truncation(np.radians(23.4), spin_ratio=10.0)  # 4: a fast rotator holds level 4 further
```

### Defaults and Warnings

- The default is `"off"` (the configured `[tides]` `obliquity_trunc_lvl`) for every world and for the standalone functions (`global_potential`, `collapse_global_tides`, `tidal_potential_3d_modes`).
- A world running `calc_tides` at nonzero obliquity with the truncation off logs a warning once, since its obliquity tides are ignored.
- At level 2 or 4, `calc_tides` warns once per world and level past the level's 10% limit for its spin band and degree range (`obliquity_accuracy_limit(level, 0.1, max_degree_l, spin_ratio)`).
- A level passed directly must be tabulated. An untabulated level read from a configuration (e.g., the earlier levels 1 and 10) is promoted with a once-per-session warning: 1 to 2, and anything past 4 to `"gen"`.

> [!WARNING]
> The obliquity truncation compounds with the eccentricity truncation: high levels of both give many tidal modes and slow calculations.

## Physics

Kaula (1966, Eq. 3.62) writes the obliquity functions as a triple sum. TidalPy evaluates the equivalent single sum in the half-angle form (Gooding and Wagner 2008; Renaud et al. 2021, Eq. C7),

$$F_{lmp}(I) = (-1)^{\lfloor (l-m+1)/2 \rfloor}\,\frac{(l+m)!}{2^{l}\,p!\,(l-p)!}\sum_{\lambda=\lambda_{1}}^{\lambda_{2}}(-1)^{\lambda}\binom{2l-2p}{\lambda}\binom{2p}{l-m-\lambda}\cos^{3l-m-2p-2\lambda}\!\left(\frac{I}{2}\right)\sin^{m-l+2p+2\lambda}\!\left(\frac{I}{2}\right),$$

with $\lambda_{1} = \max(0,\, l-m-2p)$ and $\lambda_{2} = \min(l-m,\, 2l-2p)$. The sign factor makes $F_{lmp}$ equal to Kaula's triple sum, sign included. The global (1D) heating uses $F_{lmp}^{2}$ and does not depend on it, but the [3D potential](multilayer_3d_heating.md) is linear in $F_{lmp}$ and does. The Taylor series of $F_{lmp}$ in $I$ starts at $I^{\lvert l - m - 2p \rvert}$ and steps by $I^2$, so the modes with $l - m - 2p \neq 0$, such as $m = 1$ at $l = 2$, vanish at $I = 0$.

## Truncation Rule

As for the [eccentricity functions](Eccentricity.md#truncation-rule), level $N$ keeps every product of two obliquity functions through $I^N$ and nothing past it, so the heating, quadratic in the potential, is complete through $I^N$.

- Level $N$ holds every $F_{lmp}$ whose series starts at or below $I^N$, each through $I^N$, so a potential built from them is complete through $I^N$.
- The global (1D) heating uses each $F_{lmp}^2$ cut at $I^N$; only functions starting at or below $I^{N/2}$ enter it.
- The 3D secular heating cuts each product of two obliquity functions at $I^N$ the same way, alongside the eccentricity cut at $e^N$. The instantaneous 3D fields use the unsquared functions.
- Level `0` (`"off"`) is $F$ at $I = 0$.
- `"gen"` (`"general"`, `OBLIQUITY_GENERAL`) evaluates the half-angle sum above, exact at any obliquity, and its products are plain products.

For example, at degree 2 and level 2 the heating takes $F_{2,2,0}^2 = 9 - 9 I^2$ (from $F_{2,2,0} = 3 - \tfrac{3}{2} I^2 + \dots$), $F_{2,1,0}^2 = F_{2,1,1}^2 = \tfrac{9}{4} I^2$ (from $F_{2,1,0} = \tfrac{3}{2} I + \dots$ and $F_{2,1,1} = -\tfrac{3}{2} I + \dots$), and $F_{2,0,1}^2 = \tfrac{1}{4} - \tfrac{3}{4} I^2$. The $I^2$ corrections of the aligned modes matter for a body not in synchronous rotation. In synchronous rotation only the $m = 1$ modes dissipate through obliquity, and level 2 gives their full $I^2$ heating, the traditional obliquity-tide term.

All coefficients are exact rationals, from the code generator of Renaud et al. (2021) (available on request). The half-angle sum agrees with a 40-digit evaluation to $6 \times 10^{-13}$ relative for every degree from 2 to 10.

## References

- Kaula, W. M. (1964). Tidal dissipation by solid friction and the resulting orbital evolution. *Reviews of Geophysics*, 2(4), 661-685. Eq. 1, the tidal potential.
- Kaula, W. M. (1966). *Theory of Satellite Geodesy: Applications of Satellites to Geodesy*. Blaisdell. Eq. 3.62, the triple sum.
- Gooding, R. H., and Wagner, C. A. (2008). On the inclination functions and a rapid stable procedure for their evaluation together with derivatives. *Celestial Mechanics and Dynamical Astronomy*, 101. The single-sum form.
- Renaud, J. P., et al. (2021). Tidal dissipation in dual-body, highly eccentric, and nonsynchronously rotating systems: Applications to Pluto-Charon and the exoplanet TRAPPIST-1e. *The Planetary Science Journal*, 2(1), 4. Appendix C.
