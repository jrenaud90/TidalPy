# Numerics (`Utilities.math`)

_Updated: 2026-10-01_

Three small C++ functions, in the header-only `numerics_.hpp`, that the physics modules call in place of their standard-library equivalents. There is no Python wrapper.

## Floating-point Comparison

```cpp
bool c_isclose(double value_a, double value_b, double rtol = 1e-9, double atol = 0.0);
```

Mirrors Python's `math.isclose`. Two values are close when their absolute difference is within the larger of the relative tolerance scaled by the bigger magnitude and the absolute tolerance. Exact equality short-circuits to true. Any NaN input returns false, so a NaN never compares close to anything, including itself. An infinite value is close only to the same infinity.

The default absolute tolerance is zero, which means values near zero compare close only when they are equal. Pass a non-zero `atol` when comparing against zero.

## Guarded Growth Functions

```cpp
double c_safe_pow(double base, double exponent);
double c_safe_exp(double exponent);
```

Both wrap the standard-library function and return a quiet NaN whenever the result is not finite, whether from overflow to infinity or from a domain error such as a negative base raised to a fractional power.

A model evaluated far outside its fitted range (a radiogenics model at an epoch long before its reference time, say) can overflow. An infinity can turn into a finite, wrong number several steps later, while a NaN reaches the output and makes the problem visible. The cost is that a NaN does not say whether the true answer was large or undefined; a caller that cares must check its inputs.

## Where Numerics are Used

The melt-weakening laws guard their power laws and exponentials, the equation-of-state laws their thermal scaling, and both decaying radiogenics models their decay exponentials. The viscosity models intentionally keep the plain exponential because their cold limit is an infinite (rigid) viscosity, not NaN. Each module's page notes where a NaN can appear and what it means there. See [Viscosity Models](../Viscosity/viscosity_models.md), [Melting Laws](../PartialMelt/partial_melt_models.md), and [Radiogenic Models](../Radiogenics/radiogenics_models.md).
