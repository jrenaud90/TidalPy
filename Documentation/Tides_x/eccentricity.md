# Eccentricity Functions

_Updated: 2026-09-18_

The eccentricity functions $G_{l,p,q}(e)$ are one of the two drivers of the tidal potential, alongside the [obliquity functions](obliquity.md); see $G_{lpq}(e)$ in Eq. 1 of [Kaula (1964)](http://doi.wiley.com/10.1029/RG002i004p00661). Unlike the obliquity functions they are defined by an infinite sum over $q$ and cannot be written down exactly, so a truncation level has to be chosen. As long as $e < 1$ that choice trades accuracy against the number of active tidal modes, and therefore against computation time.

The values returned are the unsquared $G_{l,p,q}(e)$.

TidalPy does not evaluate the series at runtime. The terms for each degree and truncation are generated analytically ahead of time and compiled in, which is faster and discards modes that cannot contribute: a mode whose $G_{lpq}(e)$ is identically zero is never returned.

> [!NOTE]
> TidalPy's eccentricity functions require $0 \le e < 1$; parabolic and hyperbolic orbits are outside their domain.

## Choosing a Truncation

The count below is the number of non-zero modes at degree $l = 2$; higher degrees activate more. Truncation $n$ keeps every term of $G_{l,p,q}(e)$ through $e^{n}$. Heating goes as the potential squared, so it is complete through at least $e^{n}$ as well.

| Truncation | Modes at $l=2$ | Notes |
|---|---|---|
| 1 | 9 | Reproduces the traditional synchronous heating formula, $(21/2)(k_2/Q) G M^2 R^5 n e^2 / a^6$, at small $e$. Exact factors such as $(1-e^{2})^{-3/2}$ are kept rather than expanded. |
| 2 | 13 | Common choice. |
| 3 | 19 | Common choice, and TidalPy's default. |
| 4 | 25 | |
| 5 | 31 | |
| 10 | 61 | Used in [Renaud et al. (2021)](https://iopscience.iop.org/article/10.3847/1538-4357/abc0f2). |
| 15 | 91 | |
| 20 | 121 | |

The very large truncations can carry numerical error. Treat their results with suspicion and test sensitivity by pushing $e$ into the range where the extra terms start to matter, then check that the answer is stable.

Renaud et al. (2021) found that $e^{10}$ terms matter once eccentricity passes roughly 0.5: at $e = 0.8$, heating computed at truncation 5 is up to two orders of magnitude below the truncation 10 result.

The cost is worse than linear. Each new truncation activates new tidal modes, and new modes can introduce new unique forcing frequencies, each of which needs its own Love number solve. Doubling the mode count more than doubles the work, and the effect compounds at $l > 2$ because higher degrees activate more modes of their own. We recommend the lowest truncation that covers your eccentricity values. In a numerical integration where eccentricity may be driven higher (_e.g._, in a mean motion resonance), a truncation adequate for the initial eccentricity may become inaccurate as the eccentricity rises.

## Example

```python
from TidalPy.Tides_x.eccentricity import eccentricity_func

eccentricity = 0.1                         # 0 <= e < 1
modes_by_lpq, modes_by_lp = eccentricity_func(eccentricity, degree_l=2, truncation=4)

print(len(modes_by_lpq))                   # 25 non-zero modes
print(modes_by_lpq[(2, 0, 0)])             # one mode by its (l, p, q) key

for (l, p), by_q in modes_by_lpq.items():  # walk the modes grouped by (l, p)
    for (q,), value in by_q:
        print(f"G({l},{p},{q}) = {value:0.3e}")
```

Degrees $l = 2$ through $10$ are supported. Higher degrees require generating and compiling more terms; open a GitHub issue if you need them.

The function returns a pair of lookup objects holding the same numbers two ways: `modes_by_lpq` is keyed by the full `(l, p, q)` mode, and `modes_by_lp` is a dict keyed by `(l, p)` whose values iterate as `((q,), value)` pairs, which is the convenient form when you want every $q$ at a given $(l, p)$. Only non-zero modes appear in either.

The same functions exist in C++ and Cython for callers who need them: `Tides_x/eccentricity/eccentricity_common` carries the types, `eccentricity_driver` dispatches on degree and truncation, and the remaining files hold the generated terms.

## Where Eccentricity Functions are Used

Most users never call `eccentricity_func` directly. A world's `[tides]` configuration sets the truncation for every tidal solve it runs:

```python
world.set_tide_config(max_degree_l=2, eccentricity_truncation=3, obliquity_truncation=0)
```

The same value can be given in a world TOML file as `eccentricity_trunc_lvl`. See [Global Tides](global_tides.md) and the [TOML schema](../structures_x/config/toml_schema.md).
