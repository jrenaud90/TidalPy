
# Radial Solver Helper Functions

TidalPy's radial solver function is particular on the format of its inputs. It can be easy to make a mistake which leads to solution failures, exceptions, or crashes[^1].
The helper functions described here are designed to give users an easy interface to provide data that is then translated into the inputs required by the `radial_solver` function.

These functions are not designed to be particularly efficient. So it is recommended to use them until you get comfortable with the kind of inputs `radial_solver` requires.
At that point it would be more efficient to make the inputs correctly from the beginning and forego using these functions.
That being said, benchmarking shows that they only add about a 2% -- 10% overhead depending on the layer structure.

[^1]: If a crash does occur, please report it on TidalPy's GitHub issues page. Please include the exact inputs used.

## Planet with homogeneous layers: `build_planet_constant_layers`

Import with `from TidalPy.RadialSolver import build_rs_input_homogeneous_layers`

Creates radial solver inputs based on user provided parameters for a planet with homogeneous layers (each layer has a constant density, viscosity, shear, etc.).
Checks will be performed to ensure that the inputs are valid.

Arguments and use case:
```python
from TidalPy.RadialSolver import build_rs_input_homogeneous_layers

# Create inputs that define our planet. Let's assume it is a three layer planet with layers: solid-liquid-solid.
planet_radius = 1.0e6
forcing_frequency = 2.0 * 3.14 / (86400.0 * 2.0)
density_tuple                 = (8000.0, 6000.0, 3400.0)
static_bulk_modulus_tuple     = (500.0e9, 100.0e9, 200.0e9)
static_shear_modulus_tuple    = (200.0e9, 0.0, 50.0e9)  # Note that liquid layers have 0 shear.
bulk_viscosity_tuple          = (1.0e30, 1.0e30, 1.0e30)  # Since we are using Elastic bulk rheology in every layer, these values are unused.
shear_viscosity_tuple         = (1.0e28, 1.0e3, 1.0e20)
layer_type_tuple              = ('solid', 'liquid', 'solid')
layer_is_static_tuple         = (False, True, False)
layer_is_incompressible_tuple = (False, False, False)

# The inputs to this method are pretty self explanatory except for the rheology model inputs.
# These are provided as _instantiated_ rheology classes like so:
from TidalPy.rheology import Maxwell, Elastic, Andrade
# Some rheologies have other properties that can be changed when they are created.
andrade_alpha = 0.25
andrade_zeta = 0.01
andrade_args = (andrade_alpha, andrade_zeta)
shear_rheology_model_tuple = (Maxwell(), Elastic(), Andrade(andrade_args)) 
# We assume there is no bulk dissipation so its rheology is purely elastic.
bulk_rheology_model_tuple  = (Elastic(), Elastic(), Elastic())

# There are three different ways we can define our layer structure size. Only one of these should be provided.
radius_fraction_tuple    = (0.3, 0.6, 1.0)  # By their relative upper radius / planet radius
thickness_fraction_tuple = (0.3, 0.3, 0.4)  # By their thickness / planet radius
volume_fraction_tuple    = (0.1, 0.4, 0.5)  # Or by their volume fraction (relative to total planet's volume)
# We will just use radius fraction and set the others to None
thickness_fraction_tuple = None
volume_fraction_tuple = None

# The number of sub slices in each layer can be set different for each layer by:
slices_tuple = (10, 20, 50)
# Or set the same for each layer
slices_tuple = None
slice_per_layer = 10
# Note that each layer must have at least 5 slices.

rs_inputs = build_rs_input_homogeneous_layers(
    planet_radius,                    # Radius of planet (float64) [m]
    forcing_frequency,                # Forcing frequency, used to solve for the complex shear / bulk modulus (float64) [rad s-1]
    density_tuple,                    # Tuple of floats for each layer's constant density. (Tuple[float64]; len = num_layers)
    static_bulk_modulus_tuple,        # Tuple of floats for each layer's constant static bulk modulus. (Tuple[float64]; len = num_layers)
    static_shear_modulus_tuple,       # Tuple of floats for each layer's constant static shear modulus. (Tuple[float64]; len = num_layers)
    bulk_viscosity_tuple,             # Tuple of floats for each layer's constant bulk viscosity. (Tuple[float64]; len = num_layers)
    shear_viscosity_tuple,            # Tuple of floats for each layer's constant shear viscosity. (Tuple[float64]; len = num_layers)
    layer_type_tuple,                 # Tuple of strings for each layer type. (Tuple[str]; len = num_layers)
    layer_is_static_tuple,            # Tuple of booleans for if each layer should use the static assumption. (Tuple[bool]; len = num_layers)
    layer_is_incompressible_tuple,    # Tuple of booleans for if each layer should use the incompressible assumption. (Tuple[bool]; len = num_layers)
    shear_rheology_model_tuple,       # Tuple of rheology instances for each layer's complex shear calculation. (Tuple[RheologyModelBase]; len = num_layers)
    bulk_rheology_model_tuple,        # Tuple of rheology instances for each layer's complex bulk calculation. (Tuple[RheologyModelBase]; len = num_layers)
    # One of the following tuples must be provided. They define the layer geometry
    radius_fraction_tuple,            # (optional, default=None) Tuple of floats for each layer's radius fraction (R_layer / R_Planet).  (Tuple[float64]; len = num_layers)
    thickness_fraction_tuple,         # (optional, default=None) Tuple of floats for each layer's thickness fraction ((R_layer - R_layer_below) / R_Planet).  (Tuple[float64]; len = num_layers)
    volume_fraction_tuple,            # (optional, default=None) Tuple of floats for each layer's volume fraction (V_layer / V_Planet).  (Tuple[float64]; len = num_layers)
    # Each layer is further sub-divided into slices. At least 5 slices per layer is required. 
    slices_tuple,                     # (optional, default=None) Tuple of ints for the number of subslices each layer should be built with. (Tuple[int]; len = num_layers)
    slice_per_layer,                  # (optional, default=10) Number of slices to build all layers with. Used for each layer if `slices_tuple` is not provided (int)
    perform_checks = True             # (optional, default=True) Flag to tell function to perform additional checks on inputs. There is a small performance hit, but recommended unless you are sure your input is valid. (boolean)
)

# The output is a named tuple with the following attributes:
rs_inputs.radius_array
rs_inputs.density_array
rs_inputs.complex_bulk_modulus_array
rs_inputs.complex_shear_modulus_array
rs_inputs.frequency
rs_inputs.planet_bulk_density
rs_inputs.layer_types
rs_inputs.is_static_bylayer
rs_inputs.is_incompressible_bylayer
rs_inputs.upper_radius_bylayer_array

# These are all of the required positional arguments of TidalPy's `radial_solver` function. 
# They can be easily passed to the solver
solution = radial_solver(
    *rs_inputs,
    # Any changes to keyword arguments here....
    )
```

## Planet with inhomogeneous layers: `build_rs_input_from_data`

If your planet has an interior structure already defined by data arrays (these could be from the literature or from a much more robust equation of state than TidalPy has built in) then it is usually still a good idea to parse these arrays to ensure they are properly formatted to work with `radial_solver`.
That is where the `build_rs_input_from_data` helper function comes in.

```python
from TidalPy.RadialSolver import build_rs_input_from_data
import numpy as np
# Let's assume that I have a robust interior model found by a 3rd party EOS solver. These are saved in a numpy data file.

vulcan_data = np.load('planet_vulcan.npy')
# This planet has three layers with the following upper radii [m]
layer_upper_radius_tuple = (1.0e6, 3.0e6, 6.0e6)
# The last value in the above tuple must be the radius of the planet and it must match the last value of the radius array.
radius_array               = vulcan_data[:, 0] * 1000.0  # Assume radius was saved as km so need to convert.
density_array              = vulcan_data[:, 1]
static_bulk_modulus_array  = vulcan_data[:, 2]
static_shear_modulus_array = vulcan_data[:, 3]
shear_viscosity_array      = vulcan_data[:, 4]

# The model I am using does not provide values for bulk viscosity, but I am using an elastic rheology anyways.
# I still need to provide an array though.
bulk_viscosity_array = np.zeros_like(radius_array)

# The other inputs are pretty much the same as the example above.
forcing_frequency = 2.0 * 3.14 / (86400.0 * 2.0)
layer_type_tuple              = ('solid', 'liquid', 'solid')
layer_is_static_tuple         = (False, True, False)
layer_is_incompressible_tuple = (False, False, False)

# The inputs to this method are pretty self explanatory except for the rheology model inputs.
# These are provided as _instantiated_ rheology classes like so:
from TidalPy.rheology import Maxwell, Elastic, Andrade
# Some rheologies have other properties that can be changed when they are created.
andrade_alpha = 0.25
andrade_zeta = 0.01
andrade_args = (andrade_alpha, andrade_zeta)
shear_rheology_model_tuple = (Maxwell(), Elastic(), Andrade(andrade_args)) 
# We assume there is no bulk dissipation so its rheology is purely elastic.
bulk_rheology_model_tuple  = (Elastic(), Elastic(), Elastic())

rs_inputs = build_rs_input_from_data(
    forcing_frequency,               # Forcing frequency, used to solve for the complex shear / bulk modulus (float64) [rad s-1]
    radius_array,                    # np.ndarray for radius values throughout the planet. Must end at the planet's surface (last value = planet radius) [m]. (np.ndarray[double]; len = total_slices)
    density_array,                   # np.ndarray for density values defined at each radius value [kg m-3]. (np.ndarray[double]; len = total_slices)
    static_bulk_modulus_array,       # np.ndarray for static bulk modulus values defined at each radius value [Pa]. (np.ndarray[double]; len = total_slices)
    static_shear_modulus_array,      # np.ndarray for static shear modulus values defined at each radius value [Pa]. (np.ndarray[double]; len = total_slices)
    bulk_viscosity_array,            # np.ndarray for bulk viscosity values defined at each radius value [Pa s]. (np.ndarray[double]; len = total_slices)
    shear_viscosity_array,           # np.ndarray for shear viscosity values defined at each radius value [Pa s]. (np.ndarray[double]; len = total_slices)
    layer_upper_radius_tuple,        # Tuple of floats for each layer defining its upper radius value [m]. (Tuple[double]; len = num_layers)
    layer_type_tuple,                # Tuple of strings for each layer type. (Tuple[str]; len = num_layers)
    layer_is_static_tuple,           # Tuple of booleans for if each layer should use the static assumption. (Tuple[bool]; len = num_layers)
    layer_is_incompressible_tuple,   # Tuple of booleans for if each layer should use the incompressible assumption. (Tuple[bool]; len = num_layers)
    shear_rheology_model_tuple,      # Tuple of rheology instances for each layer's complex shear calculation. (Tuple[RheologyModelBase]; len = num_layers)
    bulk_rheology_model_tuple,       # Tuple of rheology instances for each layer's complex bulk calculation. (Tuple[RheologyModelBase]; len = num_layers)
    perform_checks = True,           # (optional, default=True) Flag to tell function to perform additional checks on inputs. There is a small performance hit, but recommended unless you are sure your input is valid. (boolean)
    warnings = True                  # (optional, default=True) Flag to tell function to raise warnings if it has to make corrections to input arrays. (boolean)
)

# Let's assume in this example my radius array did not start at r=0.0.
# It also did not have layer 1's lower radius value listed twice at the interface between layer 0 and 1.
# These are both required and would break `radial_solver`. However, since we used this helper function it will automatically fix these problems for us.

# Just as with the previous helper, the output is a named tuple with the following attributes:
rs_inputs.radius_array
rs_inputs.density_array
rs_inputs.complex_bulk_modulus_array
rs_inputs.complex_shear_modulus_array
rs_inputs.frequency
rs_inputs.planet_bulk_density
rs_inputs.layer_types
rs_inputs.is_static_bylayer
rs_inputs.is_incompressible_bylayer
rs_inputs.upper_radius_bylayer_array

# These are all of the required positional arguments of TidalPy's `radial_solver` function. 
# They can be easily passed to the solver
solution = radial_solver(
    *rs_inputs,
    # Any changes to keyword arguments here....
    )
```
