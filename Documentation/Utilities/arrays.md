# Arrays and Interpolation (`Utilities.arrays`)

_Updated: 2026-10-01_

`TidalPy.Utilities.arrays` provides linear interpolation over a sorted grid. The tabulated material laws, the layer profiles, and the radial solver's dense output use it to evaluate between samples. `interp` wraps the same C++ routine the solvers call, so Python and C++ results agree exactly.

## Python API

```python
import numpy as np
from TidalPy.Utilities.arrays import interp

sample_x = np.array([0.0, 1.0, 2.0, 3.0])
sample_y = np.array([0.0, 10.0, 5.0, 7.0])

interp(1.5, sample_x, sample_y)                # 7.5, a float for a scalar query
interp([-1.0, 0.5, 9.0], sample_x, sample_y)   # array([0., 5., 7.]), clamped at both ends
```

`interp(x, xp, fp)` takes the query coordinate or coordinates, the sample coordinates sorted ascending, and the sample values. It returns a float for a scalar query and a `float64` array shaped like `x` otherwise. Complex sample values give a complex result.

The behavior matches `numpy.interp`, including clamping out-of-range queries to the nearest endpoint value and NumPy's fallback when a slope comes out NaN. Limits:

- The sample coordinates are not checked for order (the check would cost more than the interpolation). Unsorted input gives undefined results, not an error.
- Empty sample arrays, or arrays of different lengths, raise `ValueError`.

## C++ API

The implementation is in `interp_.hpp` (namespace `tidalpy`), header only and fully inline.

```cpp
#include "interp_.hpp"

std::vector<double> radius  = {0.0, 1.0e6, 2.0e6};
std::vector<double> density = {5000.0, 4000.0, 3000.0};

const double value = tidalpy::c_interp(0.5e6, radius.data(), density.data(), radius.size());
```

| Function | Description |
|---|---|
| `c_interp(desired_x, x_domain, dependent_values, len_x, guess = 0)` | Linear interpolation of real values. Clamps out-of-range queries to the endpoint value and returns NaN only for an empty domain. |
| `c_interp_complex(...)` | The same for complex values, interpolating the real and imaginary parts independently. |
| `c_binary_search_with_guess(key, array, length, guess, int& code)` | The index search underneath both. Returns the index `j` with `array[j] <= key < array[j+1]`, returns `length` past the right end, and sets `code = -1` while returning 0 left of the array. Requires a length of at least three. |

`guess` seeds the binary search: pass zero for an isolated lookup, or the previous result index when walking a monotonic sequence of queries, which makes each search effectively constant time. Short domains skip the search: an empty domain gives NaN, a single sample gives its value, and two samples interpolate over the one interval.

## Where Interpolation is Used

The `interpolate` equation-of-state, shear-modulus, and viscosity laws (in radius) and the `interpolate` melting curve (in pressure) read their tables through it, at the same cost anywhere in a long table (see [Equation-of-State and Shear-Modulus Laws](../Material/material_eos.md)). The layer profiles, the radial solver's retained solution, and the equation-of-state solution use it to answer queries between stored slices. The search is adapted from NumPy's compiled interpolation.
