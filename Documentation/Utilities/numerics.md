# Numerics (`Utilities.math`)

_Updated: 2026-09-30_

Three small C++ functions, in the header-only `numerics_.hpp`, that the physics modules call in place of their standard-library equivalents.

There is no Python or Cython wrapper. This is infrastructure used from inside C++ model code.

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

Scientific formulas with exponential or power-law growth are evaluated at parameter values chosen by the user. A viscosity model given a temperature far outside the range its Arrhenius parameters were fitted to, or a radiogenics model evaluated at an epoch before its reference time, will overflow. An infinity then propagates through sums and ratios and can emerge as a finite, wrong number several steps later. A NaN propagates to everything that depends on it and shows up in the output, which makes the problem visible.

Returning NaN loses the information that the true answer was large rather than undefined, and a caller that wants to distinguish the two has to check its inputs before calling. In exchange, no silently wrong finite number leaves the function.

## Where Numerics are Used

The partial-melt models guard their power laws and both decaying radiogenics models guard their decay exponentials. The viscosity models intentionally keep the plain exponential because their cold limit is an infinite (rigid) viscosity, not NaN. Each module's page notes where a NaN can appear and what it means there. See [Viscosity Models](../Viscosity/viscosity_models.md), [Partial-Melt Models](../PartialMelt/partial_melt_models.md), and [Radiogenic Models](../Radiogenics/radiogenics_models.md).
