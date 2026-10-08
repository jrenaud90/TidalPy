# Lookup Structures (`Utilities.lookups`)

_Updated: 2026-09-30_

Several TidalPy results are indexed by small integers rather than by position: the eccentricity and obliquity functions by the Kaula mode numbers $(l, m, p, q)$, per-degree Love numbers by $l$. The code that needs these lookups runs in C++ without the interpreter lock, where a Python dictionary is not available, so it uses `IntMap`.

## `IntMap`

An `IntMap` maps one to four small integers onto a stored value: `IntMap1` through `IntMap4`, and the complex-valued `IntMap1Complex` through `IntMap4Complex`. The C++ templates can store any type; Python exposes only `double` and `double complex`.

The key components are packed into one 64-bit integer (for `IntMap4`, 16 bits each of $l$, $m$, $p$, and $q$), and the entries sit contiguously as (packed key, value) pairs. Lookup and insertion are cheap and allocation-free once there is capacity, so the map works inside a mode loop. It is not a general-purpose map: each key integer must fit in 16 bits (-32768 to +32767), and keys outside that range silently collide.

## Python API

```python
from TidalPy.Utilities.lookups.intmap import IntMap3, IntMap3Complex

modes = IntMap3()

# Subscript access takes the key as a tuple.
modes[(1, 2, 3)] = 70.0
modes[(1, -2, 3)] = 170.0       # negative key components are fine

modes[(1, 2, 3)]                # 70.0
len(modes)                      # 2

for key, value in modes:        # iteration yields (key_tuple, value)
    print(key, value)

# The explicit methods take the same tuple key.
modes.set((1, 1, 2), 45.4)
modes.get((1, 1, 2))            # 45.4

modes.reserve(10)               # pre-allocate space to avoid repeated reallocation
modes.size()                    # same as len()
modes.clear()

# The complex variants behave identically and store complex values.
amplitudes = IntMap3Complex()
amplitudes[(1, 2, 3)] = 90.0 - 3.5j
```

Subscript get and set, `len`, and iteration work as for a dictionary. Iteration yields `(key_tuple, value)` pairs in insertion order, not sorted key order. `reserve(n)` has no dictionary analogue: call it before inserting a known number of entries to avoid repeated reallocation.

## C++ and Cython

The templates are in `intmap_.hpp` (template signatures) and the key-packing helpers in `keys_.hpp`. The stored type can be anything, a struct or a vector of complex values included. A `cdef` function that needs a map without Python objects cimports the declarations in `intmap.pxd`.

## Where Lookups are Used

The eccentricity and obliquity function results (see [Eccentricity Functions](../Tides/Eccentricity.md) and [Obliquity Functions](../Tides/Obliquity.md)) and the tidal mode sum, which at truncation 20 and degree 10 carries thousands of modes per evaluation.
