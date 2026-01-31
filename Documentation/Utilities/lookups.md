# Lookup Structures
TidalPy has several instances that require efficient lookup arrays (_e.g._, eccentricity and obliquity function
results). While Python dictionaries are great, they don't work in C++ where some of TidalPy's functionality lives.
The package provides the following lookup structures that can be used in C++. However, many of these are wrapped via
Cython so they can also be accessed in Cython or Python.

## `IntMap3`
`IntMap3` is a lookup array (based off C++ vectors) that takes in 3 integers (`l,m,p`) and stores a double. This is not
a hash table. Instead it converts the 3 integers (which are assumed to be less than 1023) into a single 32-bit integer.
This is then used as a key that is mapped to the provided double precision floating point number.
- The key layout is: [ 2 bits unused | 10 bits l | 10 bits m | 10 bits p ].
- Data is stored contiguously as a C++ `pair`: `pair.first` = Packed Key, `pair.second` = Value.

Example usage in Python:
```python
from TidalPy.utilities.lookups.intmap3 import Intmap3

my_map = IntMap3()

# When using the python setter/getter you provide the key in a tuple
# Try setting
my_map[(1,2,3)] = 70.0

# Try getting
print(my_map[(1,2,3)])

# Check other methods
my_map.clear()
print(my_map.size())
my_map.reserve(10)  # Reserves memory so that memory allocation can happen on your terms.

# The set method takes each integer separately, not in a tuple.
test_map.set(1, 1, 2, 45.4)
print(test_map.get())
```

Please take a look at the .pyx, .cpp, and .hpp to see how to use these structures in C++ or Cython.
