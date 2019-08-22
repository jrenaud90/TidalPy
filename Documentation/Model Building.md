# Custom Model Building
TidalPy is constructed to be very flexible in the types of models that you can construct and use seamlessly with the OOP implementation.

Some TidalPy modules (e.g., `radiogenics`, `cooling`, `compliance`, etc.) can have one or more customizable model types. You can find these by looking for the following files in a TidalPy module:

master file: "X.py"
model file: "X_models.py"

You can enter custom functions into the *model* file. The *master* file will contain the classes and functionality that allow TidalPy to access your custom functions and incorporate them in the OOP scheme.

Note: Users that are inexperienced with the TidalPy project should refrain from making changes to (or creating new) master files.

At release, all TidalPy model files will contain at least one function made by the developers. Use it, as well as the class in the master file, as a guide on how to construct your own custom functions.

## Custom Function Structure
All custom functions have the following structure:

```
def func_name(*required_inputs, *live_arguments, *const_arguments):
     """ Doc string
     
     !TPY_args live: <comma seperated argument signitures>
     !TPY_args const: <comma seperated argument signitures>
     
     --- Parameters ---
     ... 
     
     """
```

`required_inputs` represent arguments common to all functions for a given module. Your custom function must include these required inputs even if it does not end up needing them. The master file's class will **always** passes these inputs to your function so a crash will occur if you do not include them.

`!TPY_args const:`: **Constant Arguments** or *static* arguments do not change and are not updated at each call. These must be found in the `master` class' `self.config` with the same name. These could be loaded in by a planet's configuration file or added to the class' config at a later step. An example of a static argument would be a planet's `mass`.
* Look at `rheology/compliance_models.py` for some examples of functions that do and do not use const arguments.

`!TPY_args live:`: **Live Arguments** are arguments which you want to be updated *at each function call*. For example, if you are calculating convection and your material's thermal conductivity changes with temperature, then `thermal_conductivity` should be a live argument because each time the material's temperature changes then cooling will be updated with the new `thermal_conductivity`.
* Look at `thermal/cooling_models.py` for some examples of functions that use both live and const arguments.

**Notes**
* Live arguments must come before static arguments in both the function signature and the doc string.
* Live arguments work with the master class' `self`. So there must be a way for the master class to get to the desired argument via an attribute call (see example below).

## Custom Function Example

```
# While not requried, it is often a good idea to use njit on your function to increase performance
#  Do not use "from numba import njit" as there is a TidalPy configuration that won't be able to toggle njit on or off.

from TidalPy.performance import njit

@njit
def my_func(temperature, viscosity,                           # These are required inputs
            thermal_conductivity, thermal_diffusivity,        # These are live inputs
            surface_gravity, planet_mass):                    # These are static inputs
            """ my_func does X
            
            # Note that live args MUST come before static args
            !TPY_args live: self.layer.thermal_conductivity, self.layer.thermal_diffusivity
            !TPY_args const: surface_gravity, planet_mass
            
            Parameters
            ----------
            Make sure to put information about the inputs here!
            
            """
            
            code ...
            
            return ...
``` 