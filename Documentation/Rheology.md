# Rheology Functions and Classes
TidalPy provides an efficient suite of tools to calculate the complex shear modulus from a user-provided static shear, 
static viscosity, and static frequency. Note that "static" here does not mean "constant", instead it refers to the
purely real-valued versions of these parameters. All of them can change with time or radius.

## Using TidalPy's Rheology Class
TidalPy defines efficient compiled cython classes for rheology functionality. These can be imported using either cython
or python - maximizing flexibility.

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

If you are looking for a more programatic way to import a rheology, consider using the `find_rheology` function:

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
TidalPy provides helper methods to efficiently parse over arrays. These functions use multithreading when possible to 
quickly calculate results over large arrays.

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
# frequency remains a scalar:
frequency = 1.0e-5
shear_arr = np.linspace(40.0e9, 80.0e9, 10)
visco_arr = np.logspace(1.0e18, 1.0e22, 10) 1.0e18  # Viscosity and Shear must have the same shape

# Solve
complex_shear_arr = rheology_instance.vectorize_modulus_viscosity(frequency, shear_arr, visco_arr)
```

### Changing Other Rheological Parameters
Some rheologies have additional parameters, like Andrade's $\zeta$ and $\alpha$. Notice how we did not specify these
above. TidalPy provides default values that are used if no parameters are specified. You can easily provide your own
parameters however:

```python
from TidalPy.models import Andrade

# Andrade accepts additional args \alpha and \zeta. They must be provided in this order as a tuple.
rheology_instance = Andrade((0.2, 10.))  # Alpha = 0.2; Zeta = 10.

# You can later change these additional parameter arguments on the same instance
rheology_instance.change_args((0.5, 0.1))  # Changing Alpha to 0.5; Zeta to 0.1.
```

## Defining a New Rheological Model
New rheologies can be implemented by creating a new sub class by cloning the TidalPy repository and adding a new
subclass to `\TidalPy\rheology\models.pyx`. You can copy one of the rheologies already present and use them as a
template to create your new model. After you are finished you will need to add the model to the
`find_rheology` function at the top of the same file as well as the `\TidalPy\rheology\models.pxd` header file 
so they can be imported by other TidalPy methods. You will then need to reinstall TidalPy so that it can compile
your new method (use `pip install -v .` from a terminal pointing to your TidalPy directory that you made the
modifications in). 

If you plan to push this new rheology to the main TidalPy Github: please consider adding test cases to 
`Tests\Test_Functions\test_rheology.py` so that TidalPy's CI system can check that it is performing correctly in 
future releases.
