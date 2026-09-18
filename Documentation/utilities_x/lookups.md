# Lookup Structures (`Utilities_x.lookups`)

_Updated: 2026-09-16_

Several parts of TidalPy produce quantities indexed by small integers rather than by position. The eccentricity and obliquity functions are keyed by the Kaula mode numbers $(l, m, p, q)$; per-degree Love numbers are keyed by $l$. Most of the code that needs the lookup runs in C++ with the interpreter lock released, where a Python dictionary is not available; `IntMap` fills that gap.

## `IntMap`

An `IntMap` maps one to four small integers onto a stored value. There are separate classes per key count, `IntMap1` through `IntMap4`, and a complex-valued variant of each, `IntMap1Complex` through `IntMap4Complex`. The C++ templates can store any type; the Cython layer exposes only `double` and `double complex`.

These are not hash tables. Each key integer is assumed to fit in 16 bits, from -32768 to +32767, and the key components are packed into a single 64-bit integer that indexes the stored value. For `IntMap4` the layout is 16 bits of $l$, then $m$, then $p$, then $q$. The data sits contiguously as a vector of pairs, with the packed key first and the value second.

Lookup and insertion are cheap and allocation-free once the vector has capacity, so the structure is usable inside a mode loop. Keys outside the 16-bit range silently collide, so it is not a general-purpose map.

## Python API

```python
from TidalPy.Utilities_x.lookups.intmap import IntMap3, IntMap3Complex

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

modes.reserve(10)               # pre-allocate, so growth happens on your terms
modes.size()                    # same as len()
modes.clear()

# The complex variants behave identically and store complex values.
amplitudes = IntMap3Complex()
amplitudes[(1, 2, 3)] = 90.0 - 3.5j
```

The map is dictionary-like: subscript get and set, `len`, and iteration all work as expected. Iteration yields `(key_tuple, value)` pairs in storage order, which is insertion order rather than sorted key order.

`reserve(n)` is the one method with no dictionary analogue. Calling it before a loop that inserts a known number of entries avoids repeated reallocation, which matters when the loop is the inner loop of a mode sum.

## C++ and Cython

The templates live in `intmap_.hpp`, with the key-packing helpers in `keys_.hpp`. Because the C++ side is templated on the stored type, it is not restricted to the two types the Python layer exposes: a map of structs or of complex vectors is equally valid, and is how some internal mode bookkeeping is done. The Cython declarations in `intmap.pxd` are what a `cdef` function should cimport when it needs a map without going through Python objects.

See `intmap_.hpp` and `intmap.pxd` for the template signatures.

## Where Lookups are Used

The eccentricity and obliquity function results are the original consumers; see [Eccentricity Functions](../Tides_x/eccentricity.md) and [Obliquity Functions](../Tides_x/obliquity.md). The tidal mode collapse uses the same structures to carry per-mode quantities through the sum: a truncation-20 solve at degree 10 touches thousands of modes per evaluation.
