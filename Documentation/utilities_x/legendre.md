# Legendre Polynomials (`Utilities_x.legendre`)

_Updated: 2026-09-23_

The associated Legendre functions $P_{lm}(\cos\theta)$ are the angular part of the tidal potential: expanding a companion's gravitational potential over a body's surface gives one term per degree $l$ and order $m$. The three-dimensional stress, strain, and heating fields also need the first and second derivatives with respect to colatitude, because the tangential components of the tidal displacement are angular gradients of the potential.

This module returns all three at once. Both entry points give the triple $\left(P_{lm},\ dP_{lm}/d\theta,\ d^2P_{lm}/d\theta^2\right)$ for a colatitude in radians on $[0, \pi]$, using the unnormalized functions with the Condon-Shortley phase, the same convention as `scipy.special.assoc_legendre_p` with `branch_cut=2`. Normalization and phase conventions differ between references, and a mismatched convention produces tidal amplitudes that are wrong by a factor that looks almost right.

## Evaluation Paths

`legendre(l, m, colatitude)` uses precomputed closed forms for degrees 2 through 10. Each $(l, m)$ pair is a hard-coded polynomial in $\cos\theta$ and $\sin\theta$, generated from the Ferrers construction by `legendre/codegen/gen_legendre.py`. There are no calls to `pow` and no $1/\sin\theta$ factors that would blow up at the poles, so this path is fast enough for a mode loop evaluated at every grid point.

`legendre_generic(l, m, colatitude)` handles any degree. The value comes from the standard upward recurrence in degree at fixed order, written in $\cos\theta$ and $\sin\theta$ themselves (rebuilding $\sin\theta$ from $\cos\theta$ would keep only a few digits near the poles), and the derivatives from the order recurrence

$$\frac{dP_{l}^{m}}{d\theta} = \frac{1}{2}\left[P_{l}^{m+1} - (l+m)(l-m+1)\,P_{l}^{m-1}\right],$$

applied once for the first derivative and twice for the second, with $P_{l}^{-m} = (-1)^{m}\frac{(l-m)!}{(l+m)!}P_{l}^{m}$. Every term is smooth in $\theta$, so the triple is exact at and near the poles, where a chain rule through $x = \cos\theta$ gives NaN for odd $m$. It exists for degrees outside the precomputed range. The tidal pipeline as shipped uses only the fast path, since the mode sums are truncated well below degree 10.

Both agree with SciPy, and with each other, to about $10^{-11}$ across every supported degree and order, and with each other at the poles too. An order outside $0 \le m \le l$ returns a NaN triple rather than raising, and so does a degree outside 2 through 10 on the precomputed path.

## Python API

```python
from TidalPy.Utilities_x.legendre import legendre, legendre_generic

colatitude = 1.1   # [rad]

value, first_derivative, second_derivative = legendre(2, 1, colatitude)
value_high_degree, _, _ = legendre_generic(11, 4, colatitude)

legendre(11, 0, colatitude)   # (nan, nan, nan): outside the precomputed degrees
legendre(2, 5, colatitude)    # (nan, nan, nan): order exceeds the degree
```

| Function | Description |
|---|---|
| `legendre(l, m, colatitude)` | Precomputed triple for degrees 2 through 10. NaN outside that range or for an out-of-range order. |
| `legendre_generic(l, m, colatitude)` | The same triple for any degree, from the recurrences. NaN for an out-of-range order. |

## C++ API

Header only, namespace `tidalpy`.

```cpp
#include "legendre_driver_.hpp"    // precomputed tables, degrees 2 to 10, plus dispatch
#include "legendre_generic_.hpp"   // generic recurrence evaluator, any degree

const tidalpy::c_LegendreValue table   = tidalpy::c_legendre(2, 1, colatitude);
const tidalpy::c_LegendreValue generic = tidalpy::c_legendre_generic(11, 4, colatitude);
// table.p, table.dp_dtheta, table.d2p_dtheta2
```

| Symbol | Description |
|---|---|
| `struct c_LegendreValue { double p, dp_dtheta, d2p_dtheta2; }` | The value and its two colatitude derivatives, in `legendre_common_.hpp`. |
| `c_legendre(l, m, colatitude)` | The precomputed path, in `legendre_driver_.hpp`. |
| `c_legendre_generic(l, m, colatitude)` | The generic path, in `legendre_generic_.hpp`. |
| `c_legendre_l2` through `c_legendre_l10` | The per-degree table functions the driver dispatches to, one header each. Each takes the order and the already-computed sine and cosine of the colatitude. |

The degree-2 table reproduces the closed forms the tidal kernel has always used, exactly: $P_{20} = \tfrac{1}{2}(3\cos^2\theta - 1)$, $P_{21} = -3\cos\theta \sin\theta$, and $P_{22} = 3\sin^2\theta$.

## Where Legendre Functions are Used

The three-dimensional tidal potential engine evaluates these functions and their derivatives at the colatitude of every query point, once per mode, to build the potential and its angular derivatives. See [3D Tidal Heating](../Tides_x/multilayer_3d_heating.md).
