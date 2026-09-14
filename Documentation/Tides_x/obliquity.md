# Obliquity Functions

_Updated: 2026-09-13_

The obliquity functions $F_{l,m,p}(I)$ are the other forcing component of the tidal potential, alongside the [eccentricity functions](eccentricity.md); see $F_{lmp}(I)$ in Eq. 1 of [Kaula (1964)](http://doi.wiley.com/10.1029/RG002i004p00661). The potential built from them gives tidal strain, heating, and spin-orbit evolution.

> [!NOTE]
> A note on the name: the wider literature usually calls these inclination functions. For tides the relevant angle is between the deformed body's spin axis and its orbital plane about its host, which is normally called the body's obliquity. Inclination is an orbital element, but it is still what sets the relative obliquity of a host planet with respect to its tidal target's orbit. For this reason we choose to call this obliquity functions for the context of tides.

They matter even at zero obliquity. Several $(l, m, p)$ modes stay non-zero when $I = 0$, so the functions cannot simply be skipped for an aligned body. 

The values returned are the **unsquared** $F_{l,m,p}(I)$.

As with the eccentricity functions, TidalPy uses analytically pre-computed terms rather than evaluating the definition at runtime. That is much faster and it identifies the modes that cannot contribute: if $F_{lmp}(I)$ is identically zero for some $(l, m, p)$, that mode never enters the potential.

## Choosing a Truncation

The truncation sets how many terms of the Taylor expansion in $I$ are kept. Pick the lowest one that covers your problem; the counts below are the active modes at degree $l = 2$. Unlike the eccentricity functions, these functions do converge and can be written down completely.

| Truncation | Aliases | Modes at $l=2$ | When to use it |
|---|---|---|---|
| `0` | `"off"` | 2 | You know $I = 0$. Far cheaper than evaluating the general form at zero, and still returns the non-zero aligned modes. |
| `2` | `"2"` | 4 | Terms through $I$. Accurate for small obliquity, and the usual choice in the literature when obliquity is included at all. |
| `4` | `"4"` | 8 | Terms through $I^{3}$. Accurate for moderate to small obliquity. |
| `10` | `"gen"`, `"general"` | 9 | No truncation. The most accurate and the slowest. |

At zero obliquity with truncation `0`, the two surviving degree-2 terms are $F_{2,0,1} = -1/2$ and $F_{2,2,0} = 3$.

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

## Use in `World`s

As with eccentricity, the usual route is a world's `[tides]` configuration rather than a direct call:

```python
world.set_tide_config(max_degree_l=2, eccentricity_truncation=3, obliquity_truncation=0)
```

The TOML spelling is `obliquity_trunc_lvl`, which also accepts `"off"`. See [Global Tides](global_tides.md) and the [TOML schema](../structures_x/config/toml_schema.md).
