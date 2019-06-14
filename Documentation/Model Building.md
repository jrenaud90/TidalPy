# Custom Model Building
TidalPy is constructed to be very flexible in the types of models that you can construct and use seamlessly with OOP implementations.

A module can have one or more customizable model types. You can find these by looking for the following files in a TidalPy module:

master file: "X.py"
model file: "X_models.py"

You can enter custom functions into the model file. The master file will contain the classes and functionality that allow TidalPy to access your custom functions.

At release, all TidalPy model files will contain at least one function made by the developers. Use it, as well as the class in the master file, as a guide on how to construct your own custom functions.

## Custom Function Structure
All custom functions will have the following signature:

```
def func_name(*required_inputs, *live_arguments, *other_arguments):
     """ Doc string
     
     --- Inputs ---
     live args: <comma seperated arg signitures>
     other args: <comma seperated arg names> 
     
     """
```

Required inputs are used by all functions in for this class. They do not have to be listed in the doc string (though they should be!). The master file's class will **always** pass these inputs to your function. If you do not need one or more of them to make your calculations then simply ignore them in the function body.

**Live Arguments** are arguments which you want to be updated *at each call*. For example if you are calculating convection and your material's thermal conductivity changes with temperature, then `thermal_conductivity` should be a live argument because each time the material's temperature changes then cooling will be updated with the new `thermal_conductivity`.

**Other Arguments** or *static* arguments do not change and are not updated at each call. An example of a static argument could be a planet's `mass`.

**Notes**
* Live arguments must come before static arguments in both the function signature and the doc string.
* Live arguments work with the master class' `self`. So there must be a way for the master class to get to the desired argument via an attribute call (see example below).

## Custom Function Example

```
# While not requried, it is often a good idea to use njit on your function to increase performance
from TidalPy.performance import njit

@njit
def my_func(temperature, viscosity,                    # These are required inputs
            thermal_conductivity, thermal_diffusivity, # These are live inputs
            surface_gravity, planet_mass):                    # These are static inputs
            """ my_func does X
            
            --- Inputs ---
            # Note that live args MUST come before static args
            live args: self.layer.thermal_conductivity, self.layer.thermal_diffusivity
            other args: surface_gravity, planet_mass
            """
            
            code ...
            
            return ...
```

## Additional Reading
Confused? It may help to dig into how the custom functions are found and interpreted by TidalPy. You can find that implementation in the TidalPy/utilities/classes.py > `ModelSearcher` class. 