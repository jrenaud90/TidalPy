# Rheology Module
TidalPy calculates complex moduli from a user-provided static modulus, static viscosity, and frequency. "Static" here does not mean "constant": it refers to the purely real-valued versions of these parameters, all of which can vary with radius or time.

"Moduli" in this context refers to either the shear or the bulk modulus. Be sure to use the corresponding viscosity.

## Using TidalPy's Rheology Class
The rheology functionality is provided by compiled Cython classes, importable from either Cython or Python.

### Import and Basic Use
Example:

```python
from TidalPy.models import Andrade

# Create an instance of the class
rheology_instance = Andrade()

# Define required inputs
frequency = 1.0e-5
shear_mod = 50.0e9
viscosity = 1.0e18

# Solve for the complex shear modulus
complex_shear = rheology_instance(frequency, shear_mod, viscosity)
print(complex_shear)
```

For a more programmatic way to import a rheology, use the `find_rheology` function:

```python

from TidalPy.rheology import find_rheology

# Create an instance of the class
rheology_class = find_rheology("andrade")
rheology_instance = rheology_class()

# Define required inputs
frequency = 1.0e-5
shear_mod = 50.0e9
viscosity = 1.0e18

# Solve for the complex shear modulus
complex_shear = rheology_instance(frequency, shear_mod, viscosity)
print(complex_shear)
```

### Working with Arrays
Helper methods sweep over arrays, using multithreading where possible.

All arrays must be [C-contiguous](https://stackoverflow.com/questions/26998223/what-is-the-difference-between-contiguous-and-non-contiguous-arrays). If an array may not be, use the numpy function `arr = np.ascontiguousarray(arr)` before passing it to a rheology method.

```python
from TidalPy.models import Andrade

# Create an instance of the class
rheology_instance = Andrade()

# Create our arrays. There are two flavors: frequency arrays or radius arrays.
# Frequency Arrays. These are when the frequency input is vectorized and the modulus and viscosity are scalars:
import numpy as np
shear_mod = 50.0e9
viscosity = 1.0e18
freq_arr  = np.logspace(-6, -4, 10)

# Solve
complex_shear_arr = rheology_instance.vectorize_frequency(freq_arr, shear_mod, viscosity)

# Radius Arrays. These are when _both_ shear modulus and viscosity are vectorized (e.g, 1D slice of a planet), while
# frequency remains a scalar. Note you must provide both shear and viscosity as an array of equal size, even if one
# remains constant while the other varies.
frequency = 1.0e-5
shear_arr = np.linspace(40.0e9, 80.0e9, 10)
visco_arr = np.logspace(1.0e18, 1.0e22, 10) 1.0e18  # Viscosity and Shear must have the same shape

# Solve
complex_shear_arr = rheology_instance.vectorize_modulus_viscosity(frequency, shear_arr, visco_arr)
```

### Changing Other Rheological Parameters
Some rheologies have additional parameters, such as Andrade's $\zeta$ and $\alpha$. These were not specified above, so TidalPy used its default values. You can provide your own:

```python
from TidalPy.models import Andrade

# Andrade accepts additional args \alpha and \zeta. They must be provided in this order as a tuple.
rheology_instance = Andrade((0.2, 10.))  # Alpha = 0.2; Zeta = 10.

# You can later change these additional parameter arguments on the same instance
rheology_instance.change_args((0.5, 0.1))  # Changing Alpha to 0.5; Zeta to 0.1.
```

## Defining a New Rheological Model
Clone the TidalPy repository and add a new subclass to `\TidalPy\rheology\models.pyx`, using one of the rheologies already there as a template. Then add the model to the `find_rheology` function at the top of that file and to the `\TidalPy\rheology\models.pxd` header, so other TidalPy methods can import it. Reinstall TidalPy to compile the new model: run `pip install -v .` from a terminal pointing at the TidalPy directory you modified.

If you plan to push the new rheology to the main TidalPy GitHub, please add test cases to `Tests\Test_Functions\test_rheology.py` so TidalPy's CI system can check it in future releases.
