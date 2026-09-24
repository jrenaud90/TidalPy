# Obliquity Functions

_Updated: 2026-09-24_

The obliquity functions $F_{l,m,p}(I)$ are the other forcing component of the tidal potential, alongside the [eccentricity functions](eccentricity.md); see $F_{lmp}(I)$ in Eq. 1 of [Kaula (1964)](http://doi.wiley.com/10.1029/RG002i004p00661). The potential built from them gives tidal strain, heating, and spin-orbit evolution.

> [!NOTE]
> The wider literature usually calls these inclination functions. For tides the relevant angle is between the deformed body's spin axis and its orbital plane about its host, which is normally called the body's obliquity, so TidalPy calls them obliquity functions. Inclination is an orbital element, but it still sets the relative obliquity of a host planet with respect to its tidal target's orbit.

Several $(l, m, p)$ modes stay non-zero when $I = 0$, so the functions cannot be skipped for an aligned body.

The values returned are the unsquared $F_{l,m,p}(I)$.

As with the eccentricity functions, TidalPy uses analytically pre-computed terms rather than evaluating the definition at runtime. That is faster and identifies the modes that cannot contribute: if $F_{lmp}(I)$ is identically zero for some $(l, m, p)$, that mode never enters the potential.

## Physics

Kaula (1966, Eq. 3.62) writes the obliquity functions as a triple sum. TidalPy evaluates the equivalent single sum in the half-angle form (Gooding and Wagner 2008; Renaud et al. 2021, Eq. C7),

$$F_{lmp}(I) = (-1)^{\lfloor (l-m+1)/2 \rfloor}\,\frac{(l+m)!}{2^{l}\,p!\,(l-p)!}\sum_{\lambda=\lambda_{1}}^{\lambda_{2}}(-1)^{\lambda}\binom{2l-2p}{\lambda}\binom{2p}{l-m-\lambda}\cos^{3l-m-2p-2\lambda}\!\left(\frac{I}{2}\right)\sin^{m-l+2p+2\lambda}\!\left(\frac{I}{2}\right),$$

with $\lambda_{1} = \max(0,\, l-m-2p)$ and $\lambda_{2} = \min(l-m,\, 2l-2p)$. The sign factor makes $F_{lmp}$ equal to Kaula's triple sum, sign included. The global (1D) heating uses $F_{lmp}^{2}$ and does not depend on it, but the [3D potential](multilayer_3d_heating.md) is linear in $F_{lmp}$ and does. $F_{lmp}$ is even in $I$ when $l - m$ is even and odd when $l - m$ is odd, so the odd-$(l - m)$ modes, such as $m = 1$ at $l = 2$, vanish at $I = 0$.

The truncated tables are the Taylor series of this expression in $I$ with exact rational coefficients, cut after $I^{n}$. Truncation `off` evaluates it at $I = 0$, and `gen` compiles the closed form itself.

## Choosing a Truncation

Truncation $n$ keeps every term of the Taylor expansion of $F_{l,m,p}(I)$ through $I^{n}$. Heating goes as the potential squared, so it is an even function of $I$. At synchronous rotation only the $m = 1$ modes dissipate through obliquity, and truncation `1` already gives their full $I^2$ heating, the traditional obliquity-tide term. A non-synchronous body also dissipates through the aligned modes (such as the $(2, 2, 0, 0)$ tide), whose $F$ carry $I^2$ corrections ($F_{2,2,0} = 3 - 3 I^2/2 + \dots$) that only truncation `2` keeps. Pick the lowest truncation that covers your problem; the counts below are the active modes at degree $l = 2$. Unlike the eccentricity functions, these functions do converge and can be written down completely.

| Truncation | Aliases | Modes at $l=2$ | Heating complete through | When to use it |
|---|---|---|---|---|
| `0` | `"off"` | 2 | $I^{0}$ | You know $I = 0$. Far cheaper than evaluating the general form at zero, and still returns the non-zero aligned modes. |
| `1` | `"1"` | 4 | $I^{2}$ at synchronous rotation; otherwise only the $I^{0}$ terms | Terms through $I$. The traditional obliquity heating of a synchronous body, and the usual choice in the literature when obliquity is included at all. |
| `2` | `"2"` | 7 | $I^{2}$ | Terms through $I^{2}$. Small to moderate obliquity at any spin rate. |
| `10` | `"gen"`, `"general"` | 9 | Exact | No truncation: exact, and the slowest. |

At zero obliquity with truncation `0`, the two surviving degree-2 terms are $F_{2,0,1} = -1/2$ and $F_{2,2,0} = 3$.

The default everywhere is `"off"`: the `[tides]` `obliquity_trunc_lvl` of the TidalPy configuration, which a built world, a directly constructed world, and the standalone functions (`global_potential`, `collapse_global_tides`, `tidal_potential_3d_modes`) all take when no truncation is given. A world that runs `calc_tides` with a nonzero obliquity while its truncation is off logs a warning once, since its obliquity tides are then ignored.

> [!WARNING]
> The truncation level chosen for the obliquity functions is compounded by the truncation used on eccentricity. A high eccentricity truncation coupled with a high obliquity truncation will result in many tidal modes leading to slower calculations.

## Example

```python
import numpy as np
from TidalPy.Tides_x.obliquity import obliquity_func

obliquity = np.radians(20.0)           # radians
modes_by_lmp, modes_by_lm = obliquity_func(obliquity, degree_l=2, truncation=10)

print(len(modes_by_lmp))                     # 9 active modes
print(modes_by_lmp[(2, 2, 0)])               # one mode by its (l, m, p) key

for (l, m), by_p in modes_by_lm.items():     # walk the modes grouped by (l, m)
    for (p,), value in by_p:
        print(f"F({l},{m},{p}) = {value:0.3e}")
```

Degrees $l = 2$ through $10$ are supported, and the truncation accepts either the integer or its string alias. The return is a pair of lookup objects over the same numbers: `modes_by_lmp` keyed by the full `(l, m, p)` mode, and `modes_by_lm` a dict keyed by `(l, m)` whose values iterate as `((p,), value)` pairs. Only active modes appear.

The C++ and Cython entry points mirror the Python ones: `Tides_x/obliquity/obliquity_common` holds the types, `obliquity_driver` dispatches on degree and truncation, and the remaining files hold the generated terms.

## Where Obliquity Functions are Used

As with eccentricity, the usual route is a world's `[tides]` configuration rather than a direct call:

```python
world.set_tide_config(max_degree_l=2, eccentricity_truncation=6, obliquity_truncation=0)
```

The TOML spelling is `obliquity_trunc_lvl`, which also accepts `"off"`. See [Global Tides](global_tides.md) and the [TOML schema](../structures_x/config/toml_schema.md).

## References

- Kaula, W. M. (1966). *Theory of Satellite Geodesy: Applications of Satellites to Geodesy*. Blaisdell. Eq. 3.62, the triple sum.
- Gooding, R. H., and Wagner, C. A. (2008). On the inclination functions and a rapid stable procedure for their evaluation together with derivatives. *Celestial Mechanics and Dynamical Astronomy*, 101. The single-sum form.
- Renaud, J. P., et al. (2021). Tidal dissipation in dual-body, highly eccentric, and nonsynchronously rotating systems: Applications to Pluto-Charon and the exoplanet TRAPPIST-1e. *The Planetary Science Journal*, 2(1), 4. Appendix C.
