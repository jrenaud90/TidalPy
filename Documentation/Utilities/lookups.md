# Lookup Structures
Several parts of TidalPy need efficient lookup arrays (_e.g._, eccentricity and obliquity function results) in C++, where a Python dictionary is not available. The structures below fill that role, and most are wrapped in Cython so they can also be used from Cython or Python.

These structures live in `TidalPy.Utilities_x.lookups`.

## `IntMap`
`IntMap` is a lookup array, built on C++ vectors, that takes 1 to 4 integers (`l,m,p,q`; each count has its own class, `IntMap1`, `IntMap2`, and so on) and stores a double or a double complex. The C++ versions store an arbitrary data structure; only doubles and double complex numbers are exposed to Python and Cython.

These are not hash tables. The N integers, assumed to fit in 16 bits (-32768 to +32767), are packed into a single 64-bit integer, which is the key that maps to the stored value.
- The key layout for `IntMap4` is: [16 bits `l` | 16 bits `m` | 16 bits `p` | 16 bits `q`].
- Data is stored contiguously as a C++ `pair`: `pair.first` = Packed Key object, `pair.second` = Value.

Example usage in Python:
```python
from TidalPy.Utilities_x.lookups.intmap import IntMap3, IntMap3Complex

my_map = IntMap3()

# When using the python setter/getter you provide the key in a tuple
# Try setting
my_map[(1,2,3)] = 70.0

# They can be negative
my_map[(1,-2,3)] = 170.0

# Try getting
print(my_map[(1,2,3)])
print(my_map[(1,-2,3)])

# Check other methods
my_map.clear()
print(my_map.size())
my_map.reserve(10)  # Reserves memory so that memory allocation can happen on your terms.

# The set method takes each integer separately, not in a tuple.
test_map.set((1, 1, 2), 45.4)
print(test_map.get())

my_complex_map = IntMap3Complex()
my_complex_map[(1, 2, 3)] = 90 - 3.5j
print(my_complex_map[(1, 2, 3)])

```

`IntMapN` behaves much like a Python dictionary. It supports:
- Setting and getting with `my_map[...]`
- The number of items with `len(my_map)`
- Iteration (`iter`, or `for x in my_map`) with the signature `key_tuple, value in my_map`

See the .pyx, .cpp, and .hpp files for use in C++ or Cython.
