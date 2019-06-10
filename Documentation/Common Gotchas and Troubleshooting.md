# TidalPy Troubleshooting Guide

## Gotchas

### Numpy Array Shapes
TidalPy has no way to know what shape most arrays should be (as either inputs or outputs) it therefore requires that the following arrays must be the same shape:
* Temperature
* Viscosity
* Shear Modulus
* Spin Frequency
* Orbital Frequency / Semi-major Axis
 
The shape must be the same across all layers and planets.

As of 0.1.0a, TidalPy does not check if an array is the correct shape. Instead, if an array is an incorrect shape, then a later numpy array error may arise.
 
