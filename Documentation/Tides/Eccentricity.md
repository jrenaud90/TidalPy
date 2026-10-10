# Eccentricity Functions

_Updated: 2026-10-06_

The eccentricity functions $G_{l,p,q}(e)$ are one of the two drivers of the tidal potential, alongside the [obliquity functions](Obliquity.md); see $G_{lpq}(e)$ in Eq. 1 of [Kaula (1964)](http://doi.wiley.com/10.1029/RG002i004p00661). They are an infinite sum over $q$, so a truncation level has to be chosen, trading accuracy against the number of tidal modes and so against run time.

TidalPy returns the unsquared functions, the form a tidal potential uses, but by tradition names levels for the heating, which goes as their square: level $N$ keeps every product of two functions through $e^N$ ([Truncation Rule](#truncation-rule)). The coefficients are precomputed; only the [exact option](#exact-functions) evaluates at runtime. Modes that are identically zero, or start past the truncation, are never returned.

> [!NOTE]
> TidalPy's eccentricity functions require $0 \le e < 1$; parabolic and hyperbolic orbits are outside their domain.

## Example

Most users only set the truncation of a world (`eccentricity_trunc_lvl` in a TOML file; see [Global Tides](global_tides.md) and the [TOML schema](../Structures/config/toml_schema.md)):

```python
world.set_tide_config(max_degree_l=2, eccentricity_truncation=20, obliquity_truncation=0)
```

Calling the functions directly:

```python
from TidalPy.Tides.eccentricity import ECCENTRICITY_TRUNCATIONS, eccentricity_func, eccentricity_squared_func

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

Both return two lookups of the same non-zero modes: `modes_by_lpq` keyed by `(l, p, q)`, and `modes_by_lp`, keyed by `(l, p)`, whose values iterate as `((q,), value)` pairs. `truncation` defaults to the configured `[tides]` `eccentricity_trunc_lvl`, and `exact_tolerance` (default `[tides]` `eccentricity_exact_tolerance`) sets the mode range of `"exact"`. `validate_eccentricity_truncation` checks a level; `promote_eccentricity_truncation` resolves an untabulated configured level to the next tabulated one with a warning, as the world builder does. Degrees $l = 2$ through $10$ are supported (open a GitHub issue if you need more). C++ and Cython callers have the same functions in `Tides/eccentricity/eccentricity_common` and `eccentricity_driver`.

## Choosing a Truncation

The accuracy columns give the largest eccentricity up to which the degree-2 heating stays within the stated relative error of the exact value, at any spin rate. They are the worst case over constant-phase-lag, constant-time-lag, Maxwell (relaxation times $10^{-3}$ to $10^{3}$ inverse mean motions, which also covers any body whose $-\mathrm{Im}\,k$ is a positive sum of such peaks), and Andrade tides at spin rates from $-10$ to 30 times the mean motion. Errors can go either way: low for constant-time-lag tides, high for some strongly frequency-dependent ones. Mode counts are those of the degree-2 heating ($\lvert q \rvert \le N/2$).

| Truncation | Heating modes at $l=2$ | Error below $10^{-6}$ to | Error below 1% to | Warning above (max $l$ = 2 / 3 / 10) |
|---|---|---|---|---|
| 2 | 9 | 0.0002 | 0.02 | 0.075 / 0.06 / 0.025 |
| 4 | 13 | 0.005 | 0.095 | 0.18 / 0.145 / 0.065 |
| 6 | 19 | 0.035 | 0.17 | 0.265 / 0.22 / 0.1 |
| 8 | 25 | 0.07 | 0.235 | 0.335 / 0.28 / 0.135 |
| 10 | 31 | 0.11 | 0.29 | 0.395 / 0.335 / 0.16 |
| 20 | 61 | 0.27 | 0.44 | 0.5 / 0.435 / 0.24 |
| 50 | 151 | 0.445 | 0.545 | 0.57 / 0.525 / 0.345 |
| `"exact"` | depends on $e$ | any $e < 0.99$, to the tolerance | | never |

Level 10 is TidalPy's default. Level 2 is the traditional $e^2$ theory (its synchronous heating is exactly $(21/2)(k_2/Q) G M^2 R^5 n e^2 / a^6$).

Cost grows faster than linearly with the level, since new modes bring new forcing frequencies and each needs its own Love number solve, more so at $l > 2$. Use the lowest level that covers your eccentricities. If the eccentricity can rise during an integration (in a mean motion resonance, for example), the [range warning](#range-warning) flags a level that has become inaccurate.

### By Spin Rate

The limits are also measured for three bands of the spin ratio $\lvert \Omega / n \rvert$ (either sign; retrograde spins measured at -0.5, -1, -2, -5, and -10), edges included:

- Near synchronous: up to 1.5, which holds synchronous rotation and the 3:2 resonance.
- Moderate: 1.5 to 5.
- Fast: 5 to 30.

A ratio past 30, not measured, takes the limits for any spin rate. The 10% (and 1%) limits at degree 2 by band:

| Truncation | Near synchronous | Moderate | Fast | Any spin rate |
|---|---|---|---|---|
| 2 | 0.075 (0.02) | 0.165 (0.095) | 0.235 (0.125) | 0.075 (0.02) |
| 4 | 0.18 (0.095) | 0.21 (0.14) | 0.335 (0.215) | 0.18 (0.095) |
| 6 | 0.265 (0.17) | 0.27 (0.19) | 0.415 (0.285) | 0.265 (0.17) |
| 8 | 0.335 (0.235) | 0.335 (0.24) | 0.475 (0.35) | 0.335 (0.235) |
| 10 | 0.4 (0.29) | 0.395 (0.29) | 0.52 (0.4) | 0.395 (0.29) |
| 20 | 0.61 (0.505) | 0.505 (0.445) | 0.5 (0.44) | 0.5 (0.44) |
| 50 | 0.805 (0.75) | 0.8 (0.73) | 0.57 (0.545) | 0.57 (0.545) |

The low levels hold furthest for a fast rotator, and the high levels for a slow one.

### Recommending a Level

`recommend_eccentricity_truncation(eccentricity, tolerance=0.01, max_degree_l=2, spin_ratio=None)` returns the lowest level that holds a tolerance at an eccentricity, or `"exact"`. `eccentricity_accuracy_limit(level, tolerance, max_degree_l, spin_ratio=None)` returns the largest eccentricity a level holds a tolerance to. Both use the measurements for degrees 2 to 10 at tolerances $10^{-8}$, $10^{-6}$, $10^{-4}$, $10^{-3}$, $10^{-2}$, and $10^{-1}$ (a tolerance in between uses the next smaller one), taking the tightest limit of degrees 2 to `max_degree_l`. The eccentricity grid is log-spaced from $10^{-5}$ to 0.005, then steps by 0.005 (level 2 holds $10^{-4}$ to $e = 0.002$); a limit of 0.0 means the level fails already at $e = 10^{-5}$. `spin_ratio` picks the spin band; None uses the any-spin limits:

```python
from TidalPy.Tides.eccentricity import recommend_eccentricity_truncation

recommend_eccentricity_truncation(0.25)                   # 10: heating within 1% at e = 0.25
recommend_eccentricity_truncation(0.25, max_degree_l=4)   # 20: degree 4 needs more
recommend_eccentricity_truncation(0.3)                    # 20
recommend_eccentricity_truncation(0.3, max_degree_l=10)   # 50: degrees up to 10 need more still
recommend_eccentricity_truncation(0.3, tolerance=1e-6)    # 50
recommend_eccentricity_truncation(5e-4, tolerance=1e-4)   # 2
recommend_eccentricity_truncation(0.6)                    # 'exact'
recommend_eccentricity_truncation(0.6, spin_ratio=1.0)    # 50: a synchronous body holds level 50 further
```

### Range Warning

A world's `calc_tides` logs a warning, once per world and level, when its eccentricity is past the level's 10% limit for its spin band and degree range (`eccentricity_accuracy_limit(level, 0.1, max_degree_l, spin_ratio)`); the standalone functions warn once per session and level. Higher degrees lose accuracy sooner, so a world with more degrees warns sooner. Its total heating is usually closer to exact than that, since the higher degrees are weaker, but its per-degree values and 3D pattern are not.

## Physics

The eccentricity functions are Hansen coefficients (Kaula 1964),

$$G_{lpq}(e) = X^{-(l+1),\,l-2p}_{l-2p+q}(e), \qquad X^{n,m}_{k}(e) = \frac{1}{2\pi}\int_{0}^{2\pi}\left(\frac{r}{a}\right)^{n} e^{i(mf - k\mathcal{M})}\,d\mathcal{M},$$

where $r$ is the orbital distance, $a$ the semi-major axis, $f$ the true anomaly, and $\mathcal{M}$ the mean anomaly. $G_{lpq}$ is $e^{\lvert q \rvert}$ times a power series in $e^2$, and $G_{lpq} = G_{l,\,l-p,\,-q}$.

For $k \equiv l - 2p + q = 0$ the coefficient has a closed form (Laskar and Boué 2010; Renaud et al. 2021, Eq. C3). With $n' = -n = l + 1$ and $0 \le m \le n' - 2$ (the coefficient is even in $m$ and zero for larger $|m|$),

$$X^{-n',\,m}_{0}(e) = \left(1-e^{2}\right)^{3/2-n'}\sum_{j=0}^{\lfloor (n'-2-m)/2 \rfloor}\frac{(n'-2)!}{j!\,(m+j)!\,(n'-2-m-2j)!}\left(\frac{e}{2}\right)^{m+2j}.$$

TidalPy keeps these modes exact: the prefactor stays unexpanded and the sum is finite.

For $k \neq 0$ the coefficient is a double series in Bessel functions of the first kind (Veras et al. 2019; Renaud et al. 2021, Eq. C5),

$$X^{n,m}_{k}(e) = \left(1+\beta^{2}\right)^{-(n+1)}\sum_{j=0}^{\infty}(-\beta)^{j}\sum_{h=0}^{j}\binom{1+n+m}{j-h}\binom{1+n-m}{h}J_{k-m+j-2h}(ke), \qquad \beta = \frac{1-\sqrt{1-e^{2}}}{e},$$

with $J_{-\nu} = (-1)^{\nu}J_{\nu}$ and generalized binomial coefficients for negative upper arguments. $\beta$, every Bessel function, and every product are expanded in $e$ with exact rational coefficients.

## Truncation Rule

Truncation level $N$ (an even number) keeps every product of two eccentricity functions through $e^N$ and nothing past it. The heating and the orbital and spin rates are sums of such products, so they are the Taylor series of the exact result through $e^N$.

- The unsquared functions of `eccentricity_func` hold every mode with $\lvert q \rvert \le N$, each through $e^N$, so a tidal potential built from them is complete through $e^N$.
- The squares of the global (1D) heating, `eccentricity_squared_func`, hold every mode with $\lvert q \rvert \le N/2$, each square cut at $e^N$. A cut square can be negative for the highest-$\lvert q \rvert$ modes at large $e$; the sum over modes is still the heating's Taylor series.
- The 3D secular heating cuts each product of two modes that share a frequency at $e^N$ the same way, so at zero obliquity its volume integral equals the 1D heating at any eccentricity (with obliquity, see the cross terms in [3D Tidal Heating](multilayer_3d_heating.md)). The instantaneous 3D fields are linear in the potential and use the unsquared functions.
- The $k = 0$ modes are exact, and so is a product of two of them. A product of a $k = 0$ mode with another keeps the exact factor whole and cuts the other factor's series.

This is the truncation of Renaud et al. (2021), whose code generator (available on request) produced the coefficients. They match the exact Taylor coefficients of the Hansen integral to $10^{-11}$ or better at $l = 2$ to 4. At every level the synchronous constant-time-lag heating stays below the closed form of Hut (1981) and never worsens as the level rises (`Benchmarks/Tides/Hut1981_Constant_Time_Lag.ipynb`).

### Rule Choice

Cutting each function at $e^N$ and squaring without a further cut adds incomplete $e^{N+1}$ to $e^{2N}$ terms that overshoot above $e \approx 0.45$ and grow with $N$ (every function through $e^{20}$ overestimates the synchronous constant-time-lag heating by 760% at $e = 0.6$). Carrying each mode well past its leading power converges from below but needs about twice the modes. At a fixed number of modes the product cut is 10 to 100 times more accurate than either and never gets worse as the level rises.

## Exact Functions

`eccentricity_trunc_lvl = "exact"` (in Python `truncation="exact"`, or `ECCENTRICITY_EXACT`) takes the functions from the exact Kepler orbit. $G_{lpq}$ is the $k$-th Fourier coefficient in mean anomaly of $(r/a)^{-(l+1)} e^{imf}$, so sampling that function and taking one fast Fourier transform per $(l, p)$ gives every mode at once. Nothing is truncated in $e$ and products are plain products, so the result holds at any $e < 1$ with no cancellation between modes.

`eccentricity_exact_tolerance` (`[tides]`, default $10^{-4}$) sets the mode range: the kept modes leave a $q^2$-weighted tail below that fraction of the whole. That weight is the synchronous constant-time-lag heating's, the most demanding tide model measured, so the tolerance bounds its relative error (the heating is within the tolerance of Hut (1981) from $e = 0.1$ to $0.9$). The cost is modes: at $10^{-4}$, about $\lvert q \rvert \le 55$ at $e = 0.7$, 105 at 0.8, and 320 at 0.9, against 25 for level 50. Each forcing frequency needs a Love number solve, so at $e = 0.9$ a rheology tide solves a few hundred per call; the transform itself takes about 0.1 ms per degree at small $e$ and 2 ms at $e = 0.9$ (7 ms at degree 10). Past about $e = 0.99$ more than 20000 modes would be needed, and the call is refused.

> [!NOTE]
> The tolerance assumes the fast, high-$q$ modes matter less and less, which holds for fixed-Q and fixed-time-lag models and most viscoelastic bodies, but not for a body in resonance. A body solved with inertia (a dynamic layer) has natural resonant frequencies, and a tidal mode near one is amplified enormously; if it lies just outside the kept range, its heating and torque are missed. This matters only at high eccentricity, where the orbit forces very high frequencies. In the [exact-orbit benchmark](../../Benchmarks/Tides/Exact_Orbit_Tidal_Heating.ipynb), a dynamic Maxwell Io at $e = 0.8$ resonates at about 147 times its orbital frequency, and the default tolerance misses 0.13% of its heating; $10^{-10}$ includes the resonance and is right to $10^{-12}$. The tabulated levels miss such modes too. For a dynamic body at high eccentricity, tighten the tolerance until the heating converges.

### High Eccentricity

At large $e$ the highest-$\lvert q \rvert$ series converge slowly and their cut squares alternate in sign, so the heating is a sum of large canceling terms: at level 50 and degree 2 about $10^4$ times the total at $e = 0.7$ and $5 \times 10^5$ at $e = 0.8$ (ten times more at degree 3), multiplying rounding and Love-number errors. The sum is also accurate only while the modes are weighted much as a synchronous body's heating weighs them. A fast spinner with a strongly frequency-dependent rheology breaks this: a Maxwell body spinning at 11.5 times its mean motion is 7% off at $e = 0.57$ and 73% off at $e = 0.6$ at level 50. So past 5 times the mean motion no tabulated level is trusted past $e = 0.57$ (0.525 at degree 3, 0.345 at degree 10). Up to 5 times, level 50 holds 10% to about $e = 0.8$ at degree 2, where the cancellation takes over; past that only `"exact"` is safe.

## References

- Hut, P. (1981). Tidal evolution in close binary systems. *Astronomy and Astrophysics*, 99, 126-140. The exact constant-time-lag heating used to check each level.
- Kaula, W. M. (1964). Tidal dissipation by solid friction and the resulting orbital evolution. *Reviews of Geophysics*, 2(4), 661-685.
- Laskar, J., and Boué, G. (2010). Explicit expansion of the three-body disturbing function for arbitrary eccentricities and inclinations. *Astronomy and Astrophysics*, 522, A60. The closed form at $k = 0$.
- Veras, D., et al. (2019). Orbital relaxation and excitation of planets tidally interacting with white dwarfs. *Monthly Notices of the Royal Astronomical Society*, 486. The Bessel series at $k \neq 0$.
- Renaud, J. P., et al. (2021). Tidal dissipation in dual-body, highly eccentric, and nonsynchronously rotating systems: Applications to Pluto-Charon and the exoplanet TRAPPIST-1e. *The Planetary Science Journal*, 2(1), 4. Appendix C collects the forms used here.
