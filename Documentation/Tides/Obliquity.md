# Inclination Functions (classic)

The classic module keeps its pre-computed inclination terms in `TidalPy.tides.inclination_funcs`, organised by degree (`orderl2` through `orderl7`) with a separate zero-obliquity form for each. `get_inclination_func(tidal_order_lvl, inclination_nonzero)` returns the right one, and the function takes the obliquity in radians and returns the non-zero $F_{l,m,p}(I)$ terms as a flat dict keyed by the $(m, p)$ pair.

These are the **squared** terms $F^{2}_{l,m,p}(I)$, matching the classic heating expression. The new backend returns the unsquared $F_{l,m,p}(I)$ instead: at zero obliquity the classic degree-2 terms are $1/4$ and $9$ where the new ones are $-1/2$ and $3$.

These are the same functions the wider literature calls inclination functions. For tides the angle that matters is between the deformed body's spin axis and its orbital plane about its host, which is normally called its obliquity. Several modes stay non-zero at zero obliquity, which is why a dedicated zero-obliquity form exists rather than simply passing zero.

```python
from TidalPy.tides.inclination_funcs import get_inclination_func

inclination_func = get_inclination_func(tidal_order_lvl=2, inclination_nonzero=False)
terms = inclination_func(0.0)
```

## The maintained version

New development happens in the C++ backend, whose equivalent is documented in [Obliquity Functions](../Tides_x/obliquity.md). That page covers the truncation levels, their mode counts, and why the functions matter even for an aligned body.
