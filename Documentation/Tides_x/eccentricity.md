# Eccentricity Functions

_Updated: 2026-09-24_

The eccentricity functions $G_{l,p,q}(e)$ are one of the two drivers of the tidal potential, alongside the [obliquity functions](obliquity.md); see $G_{lpq}(e)$ in Eq. 1 of [Kaula (1964)](http://doi.wiley.com/10.1029/RG002i004p00661). Unlike the obliquity functions they are defined by an infinite sum over $q$ and cannot be written down exactly, so a truncation level has to be chosen. As long as $e < 1$ that choice trades accuracy against the number of active tidal modes, and therefore against computation time.

The functions TidalPy returns are unsquared, the form a tidal potential uses. Truncation levels are named by the heating, which goes as their square: level $N$ keeps every product of two eccentricity functions through $e^N$.

TidalPy does not evaluate the series at runtime. The coefficients for each degree and truncation are generated ahead of time with exact rational arithmetic and compiled in, and a mode whose $G_{lpq}(e)$ is identically zero, or whose series starts past the truncation, is never returned.

> [!NOTE]
> TidalPy's eccentricity functions require $0 \le e < 1$; parabolic and hyperbolic orbits are outside their domain.

## Physics

The eccentricity functions are Hansen coefficients (Kaula 1964),

$$G_{lpq}(e) = X^{-(l+1),\,l-2p}_{l-2p+q}(e), \qquad X^{n,m}_{k}(e) = \frac{1}{2\pi}\int_{0}^{2\pi}\left(\frac{r}{a}\right)^{n} e^{i(mf - k\mathcal{M})}\,d\mathcal{M},$$

where $r$ is the orbital distance, $a$ the semi-major axis, $f$ the true anomaly, and $\mathcal{M}$ the mean anomaly. $G_{lpq}$ is $e^{\lvert q \rvert}$ times a power series in $e^2$. The functions are symmetric, $G_{lpq} = G_{l,\,l-p,\,-q}$.

For $k \equiv l - 2p + q = 0$ the coefficient has a closed form (Laskar and Boué 2010; Renaud et al. 2021, Eq. C3). With $n' = -n = l + 1$ and $0 \le m \le n' - 2$ (the coefficient is even in $m$ and zero for larger $|m|$),

$$X^{-n',\,m}_{0}(e) = \left(1-e^{2}\right)^{3/2-n'}\sum_{j=0}^{\lfloor (n'-2-m)/2 \rfloor}\frac{(n'-2)!}{j!\,(m+j)!\,(n'-2-m-2j)!}\left(\frac{e}{2}\right)^{m+2j}.$$

TidalPy keeps these modes exact: the prefactor stays unexpanded and the sum is finite.

For $k \neq 0$ the coefficient is a double series in Bessel functions of the first kind (Veras et al. 2019; Renaud et al. 2021, Eq. C5),

$$X^{n,m}_{k}(e) = \left(1+\beta^{2}\right)^{-(n+1)}\sum_{j=0}^{\infty}(-\beta)^{j}\sum_{h=0}^{j}\binom{1+n+m}{j-h}\binom{1+n-m}{h}J_{k-m+j-2h}(ke), \qquad \beta = \frac{1-\sqrt{1-e^{2}}}{e},$$

with $J_{-\nu} = (-1)^{\nu}J_{\nu}$ and generalized binomial coefficients for negative upper arguments. $\beta$, every Bessel function, and every product are expanded in $e$ with exact rational coefficients.

## Truncation Rule

Truncation level $N$ (an even number) keeps every product of two eccentricity functions through $e^N$ and nothing past it. Since the heating and the orbital and spin rates are sums of such products, they are the Taylor series of the exact result through $e^N$.

- The unsquared functions returned by `eccentricity_func` hold every mode with $\lvert q \rvert \le N$, each through $e^N$, so a tidal potential built from them is complete through $e^N$.
- The squares used by the global (1D) heating, `eccentricity_squared_func`, hold every mode with $\lvert q \rvert \le N/2$, each square cut at $e^N$. A cut square can be negative for the highest-$\lvert q \rvert$ modes at large $e$; the sum over modes is still the heating's Taylor series.
- The 3D secular heating forms products of pairs of modes that share a frequency and cuts each at $e^N$ in the same way, so at zero obliquity its volume integral equals the 1D heating at any eccentricity (with obliquity, see the cross terms in [3D Tidal Heating](multilayer_3d_heating.md)). The instantaneous 3D fields (displacements, stress, strain, instantaneous power) are linear in the potential and use the unsquared functions.
- The $k = 0$ modes are exact, and so is a product of two of them. A product of a $k = 0$ mode with another keeps the exact factor whole and cuts the other factor's series; such products occur only in the 3D heating, in pointwise cross terms and at spin rates commensurate with the mean motion.

This is the definition of truncation used by Renaud et al. (2021).

### Rule Choice

Two other rules were tested against the exact heating before this one was adopted, with the exact functions computed by quadrature over the eccentric anomaly and checked against Hut (1981):

- Cutting each function at $e^N$ and squaring it without a further cut adds incomplete $e^{N+1}$ to $e^{2N}$ terms, mostly from the modes with $\lvert q \rvert$ near $N$, which carry only the first one or two terms of their series. Above $e \approx 0.45$ those terms overshoot and grow with $N$, so a higher level becomes worse: every function through $e^{20}$ overestimates the synchronous constant-time-lag heating by 760% at $e = 0.6$.
- Carrying each mode well past its leading power and squaring without a cut converges from below, but it needs about twice the modes of this rule for the same accuracy.

At a fixed number of modes the product cut was 10 to 100 times more accurate than either, errs low at every eccentricity, and never gets worse as the level rises.

### Tables

The tables in `Tides_x/eccentricity/eccentricity_func_l{2..10}_.hpp` come from `cpp_eccentricity_func_builder.py` in the `Tidal Derivations` folder of the [derivations repository](https://github.com/jrenaud90/derivations). They hold data only: for each (degree, level), one record per mode $(p, q)$, dense in $p$ and $q$, pointing at its coefficients $c_j$ of $e^{\lvert q \rvert + 2j}$ (converted from exact rationals to the nearest double). `eccentricity_common_.hpp` evaluates them: `c_eccentricity_mode_value` gives $G_{lpq}(e)$ and `c_eccentricity_cut_product` any product cut at $e^N$. `eccentricity_driver_.hpp` fills the `(l, p, q)` and `(l, p)` lookup maps the tide engines read.

### Validation

The generated coefficients match the exact Taylor coefficients of the Hansen integral, taken independently by a Cauchy integral in the complex $e$ plane, to $10^{-11}$ or better at $l = 2$ to 4. The synchronous constant-time-lag heating of every level stays below Hut's closed form and never gets worse as the level rises (`Tests/Test_Tides/test_b_eccentricity_funcs.py`, `Benchmarks_x/Tides/Hut1981_Constant_Time_Lag.ipynb`), and the collapsed 3D heating equals the 1D heating at $e = 0.6$ at level 50.

### High Eccentricity

The series of the highest-$\lvert q \rvert$ modes converge slowly at large $e$, and their cut squares alternate in sign. The heating is then a sum of large terms of both signs, which cancel. At level 50 and degree 2 the terms are about $10^4$ times the total at $e = 0.7$ and $5 \times 10^5$ times at $e = 0.8$, and about ten times more at degree 3. Rounding and any error in the Love numbers are multiplied by that factor, so no tabulated level is reliable past about $e \approx 0.78$ (0.75 at degree 3), even where the truncation alone would allow it (level 50 stays within 10% to $e = 0.80$). The exact functions below have no such limit.

## Exact Functions

`eccentricity_trunc_lvl = "exact"` (in Python `truncation="exact"`, or `ECCENTRICITY_EXACT`) takes the eccentricity functions from the exact Kepler orbit instead of a table. $G_{lpq}$ is the $k$-th Fourier coefficient in mean anomaly of $(r/a)^{-(l+1)} e^{imf}$, so sampling that function (Kepler's equation solved at each sample) and taking one fast Fourier transform per $(l, p)$ gives every mode at once. Nothing is truncated in $e$, and a product of two functions is the plain product, so the result holds at any $e < 1$ with no cancellation between modes.

The mode range follows from the tolerance `eccentricity_exact_tolerance` (`[tides]`, default $10^{-4}$): the modes kept are those whose $q^2$-weighted squares leave a tail below that fraction of the whole. That weight is the synchronous constant-time-lag heating's, the most demanding of the tide models measured, so the tolerance bounds its relative error; against Hut (1981) the heating is within the tolerance from $e = 0.1$ to $0.9$. The price is modes: at $10^{-4}$ it keeps about $\lvert q \rvert \le 55$ at $e = 0.7$, 105 at 0.8, and 320 at 0.9, against 25 for level 50. Each distinct forcing frequency needs its own Love number solve, so at $e = 0.9$ a rheology tide solves a few hundred of them per call. The transform itself costs about 0.1 ms per degree at small $e$ and 2 ms at $e = 0.9$ (7 ms at degree 10).

Past about $e = 0.99$ the functions would need more than 20000 modes and are refused.

## Choosing a Truncation

The accuracy columns give the largest eccentricity at which the degree-2 heating stays within the stated relative error of the exact value. They are the worst case over constant-phase-lag, constant-time-lag, and Maxwell tides at spin rates of 0.5, 1, and 2.3 times the mean motion. Every level errs low. The modes are those that enter the heating at degree 2 ($\lvert q \rvert \le N/2$); higher degrees activate more.

| Truncation | Heating modes at $l=2$ | Error below $10^{-6}$ to | Error below 1% to | Warning above ($l = 2$ / $l \ge 3$) |
|---|---|---|---|---|
| 2 | 9 | below 0.005 | 0.02 | 0.075 / 0.06 |
| 4 | 13 | 0.005 | 0.095 | 0.19 / 0.155 |
| 6 | 19 | 0.035 | 0.175 | 0.285 / 0.24 |
| 8 | 25 | 0.07 | 0.245 | 0.36 / 0.31 |
| 10 | 31 | 0.11 | 0.31 | 0.405 / 0.36 |
| 20 | 61 | 0.295 | 0.49 | 0.595 / 0.535 |
| 50 | 151 | 0.58 | 0.745 | 0.78 / 0.75 |
| `"exact"` | depends on $e$ | any $e < 0.99$, to the tolerance | | never |

Level 10 is TidalPy's default. Level 2 is the traditional $e^2$ theory: its synchronous heating is exactly $(21/2)(k_2/Q) G M^2 R^5 n e^2 / a^6$.

`recommend_eccentricity_truncation(eccentricity, tolerance=0.01, max_degree_l=2)` returns the lowest level that holds a tolerance at an eccentricity, or `"exact"` when none does, and `eccentricity_accuracy_limit(level, tolerance, max_degree_l)` the largest eccentricity a level holds a tolerance to. Both read the measurements behind the table, at tolerances of $10^{-8}$, $10^{-6}$, $10^{-4}$, $10^{-3}$, $10^{-2}$, and $10^{-1}$; a tolerance in between uses the next smaller one. The limits barely depend on the spin rate, so the helper takes none:

```python
from TidalPy.Tides_x.eccentricity import recommend_eccentricity_truncation

recommend_eccentricity_truncation(0.3)                    # 10: heating within 1% at e = 0.3
recommend_eccentricity_truncation(0.3, max_degree_l=3)    # 20: degree 3 needs more
recommend_eccentricity_truncation(0.3, tolerance=1e-6)    # 50
recommend_eccentricity_truncation(0.85)                   # 'exact'
```

A world's `calc_tides` logs a warning, once per world, when its eccentricity is past the last column (the 3D calls do not check it): the point where that level's heating can be 10% or more below the exact value. Degree 3 loses accuracy at a lower eccentricity than degree 2, so a world whose tides include degree 3 or higher uses the second value, measured at degree 3.

The cost is worse than linear. Each new truncation activates new tidal modes, and new modes can introduce new unique forcing frequencies, each of which needs its own Love number solve. Doubling the mode count more than doubles the work, and the effect compounds at $l > 2$ because higher degrees activate more modes of their own. We recommend the lowest truncation that covers your eccentricity values. In a numerical integration where eccentricity may be driven higher (_e.g._, in a mean motion resonance), a truncation adequate for the initial eccentricity may become inaccurate as the eccentricity rises; the warning above flags it.

## Example

```python
from TidalPy.Tides_x.eccentricity import ECCENTRICITY_TRUNCATIONS, eccentricity_func, eccentricity_squared_func

print(ECCENTRICITY_TRUNCATIONS)            # the tabulated levels

eccentricity = 0.1                         # 0 <= e < 1
modes_by_lpq, modes_by_lp = eccentricity_func(eccentricity, degree_l=2, truncation=4)

print(len(modes_by_lpq))                   # 25 non-zero modes, |q| <= 4
print(modes_by_lpq[(2, 0, 0)])             # one mode by its (l, p, q) key

for (l, p), by_q in modes_by_lp.items():   # walk the modes grouped by (l, p)
    for (q,), value in by_q:
        print(f"G({l},{p},{q}) = {value:0.3e}")

squares_by_lpq, _ = eccentricity_squared_func(eccentricity, degree_l=2, truncation=4)
print(len(squares_by_lpq))                 # 13 modes enter the heating, |q| <= 2
```

Degrees $l = 2$ through $10$ are supported. Higher degrees require generating and compiling more terms; open a GitHub issue if you need them.

Both functions return a pair of lookup objects holding the same numbers two ways: `modes_by_lpq` is keyed by the full `(l, p, q)` mode, and `modes_by_lp` is a dict keyed by `(l, p)` whose values iterate as `((q,), value)` pairs, which is the convenient form when you want every $q$ at a given $(l, p)$. Only non-zero modes appear in either. `truncation` defaults to the `[tides]` `eccentricity_trunc_lvl` of the TidalPy configuration; `exact_tolerance` sets the mode range of `"exact"` (default: the `[tides]` `eccentricity_exact_tolerance`). `validate_eccentricity_truncation` checks a level, and `promote_eccentricity_truncation` resolves an untabulated configured level to the next tabulated one with a warning, as the world builder does.

The same functions exist in C++ and Cython for callers who need them: `Tides_x/eccentricity/eccentricity_common` carries the table format and its evaluation, `eccentricity_driver` dispatches on degree and truncation, and the remaining files hold the generated tables.

## Where Eccentricity Functions are Used

Most users never call `eccentricity_func` directly. A world's `[tides]` configuration sets the truncation for every tidal solve it runs:

```python
world.set_tide_config(max_degree_l=2, eccentricity_truncation=20, obliquity_truncation=0)
```

The same value can be given in a world TOML file as `eccentricity_trunc_lvl`. See [Global Tides](global_tides.md) and the [TOML schema](../structures_x/config/toml_schema.md).

## References

- Hut, P. (1981). Tidal evolution in close binary systems. *Astronomy and Astrophysics*, 99, 126-140. The exact constant-time-lag heating used to check each level.
- Kaula, W. M. (1964). Tidal dissipation by solid friction and the resulting orbital evolution. *Reviews of Geophysics*, 2(4), 661-685.
- Laskar, J., and Boué, G. (2010). Explicit expansion of the three-body disturbing function for arbitrary eccentricities and inclinations. *Astronomy and Astrophysics*, 522, A60. The closed form at $k = 0$.
- Veras, D., et al. (2019). Orbital relaxation and excitation of planets tidally interacting with white dwarfs. *Monthly Notices of the Royal Astronomical Society*, 486. The Bessel series at $k \neq 0$.
- Renaud, J. P., et al. (2021). Tidal dissipation in dual-body, highly eccentric, and nonsynchronously rotating systems: Applications to Pluto-Charon and the exoplanet TRAPPIST-1e. *The Planetary Science Journal*, 2(1), 4. Appendix C collects the forms used here.
