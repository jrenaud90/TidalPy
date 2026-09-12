# Eccentricity Functions (classic)

The classic module keeps its pre-computed eccentricity terms in `TidalPy.tides.eccentricity_funcs`, with one submodule per degree (`orderl2` through `orderl7`). The `eccentricity_truncations` table is keyed by truncation level and then by degree, and each function takes the eccentricity and returns the non-zero terms as a dict of $p$ to a dict of $q$ to value.

These are the **squared** terms $G^{2}_{l,p,q}(e)$, which is what the classic heating expression consumes directly. The new backend returns the unsquared $G_{l,p,q}(e)$ instead, so terms from the two modules cannot be compared directly without squaring one side.

The truncation levels available here are the even numbers 2 through 22.

```python
from TidalPy.tides.eccentricity_funcs import eccentricity_truncations

eccentricity_func = eccentricity_truncations[6][2]   # truncation level 6, degree l = 2
terms = eccentricity_func(0.1)                       # terms[p][q]
```

## The maintained version

New development happens in the C++ backend, whose equivalent is documented in [Eccentricity Functions](../Tides_x/eccentricity.md). That page explains what each truncation level costs in active modes and computation, and why the choice matters more than it looks. Its levels are numbered differently, running 1 through 5 and then 10, 15, and 20, so a level number from one module should not be carried across to the other.
