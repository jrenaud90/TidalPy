# Eccentricity Functions

_Updated: 2026-09-18_

The eccentricity functions $G_{l,p,q}(e)$ are one of the two drivers of the tidal potential, alongside the [obliquity functions](obliquity.md); see $G_{lpq}(e)$ in Eq. 1 of [Kaula (1964)](http://doi.wiley.com/10.1029/RG002i004p00661). Unlike the obliquity functions they are defined by an infinite sum over $q$ and cannot be written down exactly, so a truncation level has to be chosen. As long as $e < 1$ that choice trades accuracy against the number of active tidal modes, and therefore against computation time.

The values returned are the unsquared $G_{l,p,q}(e)$.

TidalPy does not evaluate the series at runtime. The terms for each degree and truncation are generated analytically ahead of time and compiled in, which is faster and discards modes that cannot contribute: a mode whose $G_{lpq}(e)$ is identically zero is never returned.

> [!NOTE]
> TidalPy's eccentricity functions require $0 \le e < 1$; parabolic and hyperbolic orbits are outside their domain.

## Physics

The eccentricity functions are Hansen coefficients (Kaula 1964),

$$G_{lpq}(e) = X^{-(l+1),\,l-2p}_{l-2p+q}(e), \qquad X^{n,m}_{k}(e) = \frac{1}{2\pi}\int_{0}^{2\pi}\left(\frac{r}{a}\right)^{n} e^{i(mf - k\mathcal{M})}\,d\mathcal{M},$$

where $r$ is the orbital distance, $a$ the semi-major axis, $f$ the true anomaly, and $\mathcal{M}$ the mean anomaly. $G_{lpq}$ is of order $e^{|q|}$. The functions are symmetric, $G_{lpq} = G_{l,\,l-p,\,-q}$.

For $k \equiv l - 2p + q = 0$ the coefficient has a closed form (Laskar and Boué 2010; Renaud et al. 2021, Eq. C3). With $n' = -n = l + 1$ and $0 \le m \le n' - 2$ (the coefficient is even in $m$ and zero for larger $|m|$),

$$X^{-n',\,m}_{0}(e) = \left(1-e^{2}\right)^{3/2-n'}\sum_{j=0}^{\lfloor (n'-2-m)/2 \rfloor}\frac{(n'-2)!}{j!\,(m+j)!\,(n'-2-m-2j)!}\left(\frac{e}{2}\right)^{m+2j}.$$

The prefactor is kept exact, which is why a factor such as $(1-e^{2})^{-3/2}$ appears unexpanded, and only the sum is truncated.

For $k \neq 0$ the coefficient is a double series in Bessel functions of the first kind (Veras et al. 2019; Renaud et al. 2021, Eq. C5),

$$X^{n,m}_{k}(e) = \left(1+\beta^{2}\right)^{-(n+1)}\sum_{j=0}^{\infty}(-\beta)^{j}\sum_{h=0}^{j}\binom{1+n+m}{j-h}\binom{1+n-m}{h}J_{k-m+j-2h}(ke), \qquad \beta = \frac{1-\sqrt{1-e^{2}}}{e},$$

with $J_{-\nu} = (-1)^{\nu}J_{\nu}$ and generalized binomial coefficients for negative upper arguments.

The tables are generated ahead of time with exact rational arithmetic: $\beta$, every Bessel function, and every product are expanded in $e$ and cut after $e^{n}$, so truncation $n$ holds every term of every $G_{lpq}$ through $e^{n}$, including every $|q| \le n$. The coefficients are then printed to 25 significant digits and compiled.

## Choosing a Truncation

The count below is the number of non-zero modes at degree $l = 2$; higher degrees activate more. Truncation $n$ keeps every term of $G_{l,p,q}(e)$ through $e^{n}$. Heating goes as the potential squared, so it is complete through at least $e^{n}$ as well.

| Truncation | Modes at $l=2$ | Notes |
|---|---|---|
| 1 | 9 | Reproduces the traditional synchronous heating formula, $(21/2)(k_2/Q) G M^2 R^5 n e^2 / a^6$, at small $e$. Exact factors such as $(1-e^{2})^{-3/2}$ are kept rather than expanded. |
| 2 | 13 | Common choice. |
| 3 | 19 | Common choice, and TidalPy's default. |
| 4 | 25 | |
| 5 | 31 | |
| 10 | 61 | Used in [Renaud et al. (2021)](https://doi.org/10.3847/PSJ/abc0f3). |
| 15 | 91 | |
| 20 | 121 | |

The very large truncations can carry numerical error. Treat their results with suspicion and test sensitivity by pushing $e$ into the range where the extra terms start to matter, then check that the answer is stable.

Every level is a power series in $e$ cut at a fixed order, and none is reliable at high eccentricity. Against the closed-form constant-time-lag heating of Hut (1981), exact in $e$ (see the benchmark `Benchmarks_x/Tides/Hut1981_Constant_Time_Lag.ipynb`), the error in synchronous-rotation heating is:

| Truncation | Error below $10^{-6}$ up to | Error below 1% up to |
|---|---|---|
| 2 | - | $e pprox 0.05$ |
| 5 | $e pprox 0.05$ | $e pprox 0.21$ |
| 10 | $e pprox 0.25$ | $e pprox 0.37$ |
| 15 | $e pprox 0.31$ | $e pprox 0.44$ |
| 20 | $e pprox 0.37$ | $e pprox 0.48$ |

Above $e pprox 0.5$ every level is off by tens of percent or more, and from $e pprox 0.55$ a higher level is worse, not better (truncation 20 is 8 times too high at $e = 0.6$): the series converge slowly toward the Laplace limit, $e pprox 0.663$. Renaud et al. (2021) found that the $e^{10}$ terms matter once the eccentricity passes about 0.1 to 0.3; results beyond $e pprox 0.5$ from any tabulated level, theirs included, should be treated as qualitative.

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

## References

- Kaula, W. M. (1964). Tidal dissipation by solid friction and the resulting orbital evolution. *Reviews of Geophysics*, 2(4), 661-685.
- Laskar, J., and Boué, G. (2010). Explicit expansion of the three-body disturbing function for arbitrary eccentricities and inclinations. *Astronomy and Astrophysics*, 522, A60. The closed form at $k = 0$.
- Veras, D., et al. (2019). Orbital relaxation and excitation of planets tidally interacting with white dwarfs. *Monthly Notices of the Royal Astronomical Society*, 486. The Bessel series at $k \neq 0$.
- Renaud, J. P., et al. (2021). Tidal dissipation in dual-body, highly eccentric, and nonsynchronously rotating systems: Applications to Pluto-Charon and the exoplanet TRAPPIST-1e. *The Planetary Science Journal*, 2(1), 4. Appendix C collects the forms used here.
