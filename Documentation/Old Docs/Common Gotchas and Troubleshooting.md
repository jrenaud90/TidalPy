# TidalPy Troubleshooting Guide

## Gotchas

### Numba.njit Hiding exception Trace Backs
Many functions throughout TidalPy utilize the `numba` package's `njit` wrapper to increase performance. A downside of the increased performance is a loss of information when errors are raised or when you are just doing routine debugging. It is recommended that any changes you make to the code should first be run with the `use_numba` switch (found in the `TidalPy.configurations.py` file) set to False. Debug your changes, run some tests, and then once things seem to be working as expected you can turn Numba back on and rerun your tests. 

### Numba Array Shapes
As of `Numba v0.50` there is limited support for N-dimensional numpy arrays for `N > 1`. This means that all `numba.njit`'d TidalPy functions also do not support these arrays. However, there is a fairly easy work around:
```
import numpy as np

x = np.linspace(0, 10, 10)
y = np.linspace(0, 100, 10)
xx, yy = np.meshgrid(x, y)

# xx and yy are 2D arrays, they will not currently work with TidalPy functions.
# But we can flatten them...
shape = xx.shape
xx = xx.flatten()

# ... and then pass them to TidalPy ...
from TidalPy.rheology.complex_compliance.compliance_models import maxwell
complex_comp = maxwell(xx, 1.e-11, 1.e20)

# ... finally we can convert the results back to the original shape
#  (for example, if you want to use them in contour plotting or similar)
complex_comp = complex_comp.reshape(shape)
```

### Numpy Array Shapes
TidalPy has no way to know what shape most arrays should be (as either inputs or outputs) it must infer and adapt to whatever the user uses as input, however it can only do so much. It is, therefore, recommended that every state variable have dimensions of either 0 (scalar) or equal to the largest dimensioned variable. For example, if a planet has a temperature array with dimensions 1000x20 then every other state variable should either have a dimension of 0 (scalar) or of 1000x20. TidalPy assumes the input was made such that the state variables make physical sense at each index (i.e., in our example if you input a viscosity that is 1000x20 then each index should have a logical and physical connection to the 1000x20 indices of the temperature array).

Common state variables include:
* Temperature
* Viscosity
* Shear Modulus
* Spin Frequency
* Orbital Frequency / Semi-major Axis
* Eccentricity
* Obliquity
* Time
 
*These shapes must be the same across all layers and between a planet and its tidal host.*

**Numpy Array Errors**

If there is a mismatch in variable size then TidalPy will likely (hopefully!) throw an error. This will likely be a numpy array error.
