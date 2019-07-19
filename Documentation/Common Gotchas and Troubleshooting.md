# TidalPy Troubleshooting Guide

## Gotchas

### Numba.njit hiding exception TraceBacks
Many functions throughout TidalPy utilize the `numba` package's `njit` wrapper to increase performance. A downside of the increased performance is a loss of information when errors are raised or when you are just doing routine debugging. It is recommended that any changes you make to the code should first be run with the `use_numba` switch (found in the TidalPy.configurations file) set to False. Debug your changes, run some tests, and then once things seem to be working as expected you can turn Numba back on and rerun your tests. 

### Numpy Array Shapes
TidalPy has no way to know what shape most arrays should be (as either inputs or outputs) it must infer and adapt to whatever the user uses as input, however it can only do so much. Therefore every state variable must have dimensions of either 0 (scalar) or diemsions == to the largest dimensioned variable. For example, if a planet has a temperature array with dimensions 1000x20 then every other state variable must either have a dimension of 0 (scalar) or of 1000x20. TidalPy assumes the input was made such that the state variables make physical sense at each index (i.e., in our example if you input a viscosity that is 1000x20 then each index should have a logical and physical connection to the 1000x20 indices of the temperature array).

Common state variables include:
* Temperature
* Viscosity
* Shear Modulus
* Spin Frequency
* Orbital Frequency / Semi-major Axis
* Eccentricity
* Inclination
* Time
 
The shape must be the same across all layers and planets.

**Numpy Array Errors**
As of 0.1.0a, TidalPy does not check if an array is the correct shape. Instead, a numpy error will fire if it encounters mismatched array shapes.
 
