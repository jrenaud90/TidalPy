# Eccentricity Functions

_Updated: 2026-09-12_

The eccentricity functions $G_{l,p,q}(e)$ are one of the two angular ingredients of the tidal potential, alongside the [obliquity functions](obliquity.md); see $G_{lpq}(e)$ in Eq. 1 of [Kaula (1964)](http://doi.wiley.com/10.1029/RG002i004p00661). Unlike the obliquity functions they are defined by an infinite sum over $q$ and cannot be written down exactly, so a truncation level has to be chosen. As long as $e < 1$ that choice trades accuracy against the number of active tidal modes, and therefore against computation time.

The values returned are the **unsquared** $G_{l,p,q}(e)$, which is the difference to watch when comparing against the classic `TidalPy.tides.eccentricity_funcs`: that module returns the squared terms its heating expression consumes, and numbers its truncation levels differently.

TidalPy does not evaluate the series at runtime. The terms for each degree and truncation are generated analytically ahead of time and compiled in, which is both far faster and a way to discard modes that cannot contribute: a mode whose $G_{lpq}(e)$ is identically zero is never returned.

## Choosing a truncation

The count below is the number of non-zero modes at degree $l = 2$; higher degrees activate more. Because the potential is squared to obtain heating, a truncation at $e^{n}$ in the potential gives heating accurate to $e^{2n}$.

| Truncation | Modes at $l=2$ | Potential terms | Heating terms | Notes |
|---|---|---|---|---|
| 1 | 3 | $e^{1}$ | $e^{2}$ | Close to the traditional formula, though not identical: exact factors such as $(1-e^{2})^{-3/2}$ are kept rather than expanded. |
| 2 | 9 | $e^{2}$ | $e^{4}$ | Common choice. |
| 3 | 13 | $e^{3}$ | $e^{6}$ | Common choice, and TidalPy's default. |
| 4 | 19 | $e^{4}$ | $e^{8}$ | |
| 5 | 25 | $e^{5}$ | $e^{10}$ | |
| 10 | 55 | $e^{10}$ | $e^{20}$ | Used in [Renaud et al. (2021)](https://iopscience.iop.org/article/10.3847/1538-4357/abc0f2). |
| 15 | 85 | $e^{15}$ | $e^{30}$ | |
| 20 | 115 | $e^{20}$ | $e^{40}$ | |

Three things are worth knowing before reaching for a high truncation.

The very large truncations can carry numerical error. Treat their results with suspicion and test sensitivity by pushing $e$ into the range where the extra terms start to matter, then check that the answer is stable.

High truncations are not merely academic. Renaud et al. (2021) found that $e^{10}$ terms matter once eccentricity passes roughly 0.5: at $e = 0.8$, heating computed at truncation 5 is up to two orders of magnitude below the truncation 10 result.

The cost is worse than linear. Each new truncation activates new tidal modes, and new modes can introduce new unique forcing frequencies, each of which needs its own Love number solve. Doubling the mode count more than doubles the work, and the effect compounds at $l > 2$ because higher degrees activate more modes of their own.

## Using the calculator

```python
from TidalPy.Tides_x.eccentricity import eccentricity_func

eccentricity = 0.1                     # 0 <= e < 1
by_lpq, by_lp = eccentricity_func(eccentricity, degree_l=2, truncation=4)

print(len(by_lpq))                     # 19 non-zero modes
print(by_lpq[(2, 0, 0)])               # one mode by its (l, p, q) key

for (l, p), by_q in by_lp.items():     # walk the modes grouped by (l, p)
    for (q,), value in by_q:
        print(f"G({l},{p},{q}) = {value:0.3e}")
```

Degrees $l = 2$ through $10$ are supported. The function returns a pair of lookup objects holding the same numbers two ways: `by_lpq` is keyed by the full `(l, p, q)` mode, and `by_lp` is a dict keyed by `(l, p)` whose values iterate as `((q,), value)` pairs, which is the convenient form when you want every $q$ at a given $(l, p)$. Only non-zero modes appear in either.

The same functions exist in C++ and Cython for callers who need them: `Tides_x/eccentricity/eccentricity_common` carries the types, `eccentricity_driver` dispatches on degree and truncation, and the remaining files hold the generated terms.

## Where the truncation is set in practice

Most users never call `eccentricity_func` directly. A world's `[tides]` configuration sets the truncation for every tidal solve it runs:

```python
world.set_tide_config(max_degree_l=2, eccentricity_truncation=3, obliquity_truncation=0)
```

The same value can be given in a world TOML file as `eccentricity_trunc_lvl`. See [Global Tides](global_tides.md) and the [TOML schema](../structures_x/config/toml_schema.md).
