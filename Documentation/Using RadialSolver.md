# Radial Solver Documents
_Instructions and examples on how to use TidalPy's `RadialSolver` module._

## `radial_solver` Function

## `RadialSolverSolution` Class

### Love Numbers
In addition to the full `y` radial solution results, you can also quickly access the Love numbers using:

```python
solution.love  
# Returns all Love numbers for each solution type specified.
```

These are returned in a 2D array structured like [num_solution_types : love_number] where the first column is the
solution type. If you are only solving for tidal Love numbers then there will only be one solution in the "0" location.
If you are solving for multiple solutions then they will be stored in order (e.g., if you solve for tidal and load then
tidal will be at 0 and load at 1). The second column is used to access the three Love numbers:
- index 0: k Love number.
- index 1: h Love number.
- index 2: l Shida number.

Example:

Say you want to access the h Love number for solution 0:

```python
k = solution.love[0, 1]
```

There are shortcuts to make accessing specific Love numbers easier. The following attributes will return arrays for the
respective Love number for each solution type.

```python
solution.k = # [Solution 0's k, Solution 1's k, ...]
solution.h = # [Solution 0's h, Solution 1's h, ...]
solution.l = # [Solution 0's l, Solution 1's l, ...]
```

**Performance Note**

For all of these attributes (`.love; .k; .h; .l`) the arrays are being created each time the attribute is accessed.
If you have a code that accesses these numbers often (more than once) it is better to store them in a local variable.

```python
k = solution.k  # Creates an array for all k Love numbers.

# ... Do a bunch of stuff with the new "k" variable.
```
