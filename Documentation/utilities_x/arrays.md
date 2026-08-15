# Utilities_x: Array Helpers (`arrays`)

Small, dependency-free array utilities used across the `_x` code. Currently this
holds a fast 1-D linear interpolation routine.

---

## 1-D linear interpolation

`interp` performs `numpy.interp`-style linear interpolation: given sample
coordinates `xp` (sorted ascending) and sample values `fp`, it returns the
linearly-interpolated value(s) at query coordinate(s) `x`. Queries outside
`[xp[0], xp[-1]]` clamp to the corresponding endpoint value.

It is backed by a header-only C++ implementation (`interp_.hpp`) so other C++/Cython
code can call it at full speed without the Python layer.

### Python API

```python
from TidalPy.Utilities_x.arrays import interp
import numpy as np

xp = np.array([0.0, 1.0, 2.0, 3.0])
fp = np.array([0.0, 10.0, 5.0, 7.0])

interp(1.5, xp, fp)              # 7.5  (scalar in -> float out)
interp([-1.0, 0.5, 9.0], xp, fp) # array([0. , 5. , 7. ])  (clamped at both ends)
```

#### `interp(x, xp, fp)`

| Parameter | Type | Description |
|-----------|------|-------------|
| `x` | float or array_like | Query coordinate(s) to interpolate at. |
| `xp` | array_like of float | Sample coordinates, **sorted ascending**, length >= 1. |
| `fp` | array_like of float | Sample values; same length as `xp`. |

**Returns:** a Python `float` for scalar `x`; otherwise a `float64` `numpy.ndarray`
with the shape of `x`.

**Raises:** `ValueError` if `xp` is empty or `xp` and `fp` differ in length.

**Behavior notes**

- Matches `numpy.interp` for in-range and out-of-range (clamped) queries.
- `xp` is assumed sorted ascending; results are undefined otherwise (no check, for speed).
- Equivalent to `numpy.interp` but routed through the same C++ routine the solvers use,
  so Python results and internal C++ results agree exactly.

---

## C++ API

Include `interp_.hpp` (namespace `tidalpy`). All functions are header-only and
`inline`. The domain `x_domain` must be sorted ascending with length >= 1.

```cpp
#include "interp_.hpp"

std::vector<double> r   = {0.0, 1.0e6, 2.0e6};
std::vector<double> rho = {5000.0, 4000.0, 3000.0};

double d = tidalpy::c_interp(0.5e6, r.data(), rho.data(), r.size());  // 4500.0
```

| Function | Description |
|----------|-------------|
| `double c_interp(desired_x, x_domain, dependent_values, len_x, guess = 0)` | Real linear interpolation at `desired_x`. Clamps out-of-range queries to the endpoint values; returns NaN only for an empty domain. |
| `std::complex<double> c_interp_complex(desired_x, x_domain, dependent_values, len_x, guess = 0)` | As `c_interp`, but interpolates complex values (real and imaginary parts independently). |
| `std::size_t c_binary_search_with_guess(key, array, length, guess, int& code)` | Index search underlying the interpolators. Returns `j` with `array[j] <= key < array[j+1]`; returns `length` past the right end; sets `code = -1` (returns 0) left of the array. Requires `length >= 3`. |

**The `guess` parameter** seeds the internal binary search. For a single lookup
pass `0` (the default). When interpolating a monotonic sequence of query points,
passing the previous result index accelerates the search toward `O(1)` per query.

**Edge cases handled by `c_interp` / `c_interp_complex`:**

- `len_x == 0` → NaN.
- `len_x == 1` → the single sample value.
- `len_x == 2` → direct interpolation over the one interval (the guess-search needs `len >= 3`).
- A NaN interpolation slope falls back to the value from the other bracket endpoint
  (NumPy's robustness trick).

---

## Where it is used

- `Material_x.eos.material_eos` — `c_InterpolatedEOS::calc_density` interpolates a
  layer's `density(radius)` table via `c_interp` (see
  [material EOS models](../material_x/material_eos.md)).

This is the `_x` (new class scheme) replacement for the legacy
`TidalPy.utilities.arrays` interpolation helpers.

---

## Implementation notes

- Header-only; no separate compilation unit. The search routine is adapted from
  NumPy's `compiled_interp`, so results match `numpy.interp` (including endpoint
  clamping and the NaN-slope fallback).
- The Python wrapper coerces `xp`/`fp`/`x` to contiguous `float64` arrays and loops
  over array queries in C, returning a result shaped like `x`.
