# Arrays and Interpolation (`Utilities_x.arrays`)

_Updated: 2026-09-16_

`TidalPy.Utilities_x.arrays` provides linear interpolation over a sorted grid. The tabulated equation of state uses it to get density at a radius, the layer profiles use it to get gravity and pressure between slices, and the radial solver's dense output uses it to evaluate the solution between integration steps.

`interp` is a thin wrapper over the same header-only C++ routine the solvers call, so a Python result and a C++ result agree exactly.

## Python API

```python
import numpy as np
from TidalPy.Utilities_x.arrays import interp

sample_x = np.array([0.0, 1.0, 2.0, 3.0])
sample_y = np.array([0.0, 10.0, 5.0, 7.0])

interp(1.5, sample_x, sample_y)                # 7.5, a float for a scalar query
interp([-1.0, 0.5, 9.0], sample_x, sample_y)   # array([0., 5., 7.]), clamped at both ends
```

`interp(x, xp, fp)` takes the query coordinate or coordinates, the sample coordinates sorted ascending, and the sample values. It returns a Python float for a scalar query and a `float64` array shaped like `x` otherwise, and raises `ValueError` if the sample arrays are empty or differ in length.

The behavior matches `numpy.interp`, including the clamping of out-of-range queries to the nearest endpoint value and NumPy's fallback when an interpolation slope comes out NaN. The one difference is that the sample coordinates are assumed sorted ascending, and no check is performed, because the routine sits inside inner loops where the check would cost more than the interpolation. Unsorted input gives undefined results rather than an error.

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

The `guess` argument seeds the binary search. For an isolated lookup, pass zero. When interpolating a monotonic sequence of query points, passing the previous result index turns the search from logarithmic into effectively constant time, which makes a dense-output evaluation over thousands of radial slices cheap.

Short domains are handled without the search: an empty domain gives NaN, a single sample gives that sample's value, and two samples interpolate directly over the one interval, since the guess-seeded search needs at least three points.

## Where Interpolation is Used

The tabulated equation of state interpolates its density and viscoelastic tables with `c_interp`; see [Material EOS Models](../material_x/material_eos.md). The layer equation-of-state data, the radial solver's retained solution, and the equation-of-state solution object all use it to answer queries at an arbitrary radius between stored slices.

## Implementation Notes

The search routine is adapted from NumPy's compiled interpolation, which is why the results match `numpy.interp` down to the endpoint clamping and the NaN-slope fallback. The Python wrapper coerces its three arguments to contiguous `float64` arrays, loops over the queries in C, and returns a result shaped like the query.
