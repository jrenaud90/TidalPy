# Numerics (`Utilities_x.math_x`)

_Updated: 2026-09-12_

Three small C++ functions, in the header-only `numerics_.hpp`, that the physics modules call in place of their standard-library equivalents. They exist for one reason: to make a bad parameter regime visible instead of letting it propagate as a plausible-looking number.

There is no Python or Cython wrapper. This is infrastructure used from inside C++ model code.

## Floating-point comparison

```cpp
bool c_isclose(double value_a, double value_b, double rtol = 1e-9, double atol = 0.0);
```

Mirrors Python's `math.isclose`. Two values are close when their absolute difference is within the larger of the relative tolerance scaled by the bigger magnitude and the absolute tolerance. Exact equality short-circuits to true, and any NaN input returns false, so a NaN never compares close to anything including itself.

The default absolute tolerance is zero, which means values near zero compare close only when they are equal. Pass a non-zero `atol` when comparing against zero.

## Guarded growth functions

```cpp
double c_safe_pow(double base, double exponent);
double c_safe_exp(double exponent);
```

Both wrap the standard-library function and return a quiet NaN whenever the result is not finite, whether from overflow to infinity or from a domain error such as a negative base raised to a fractional power.

The motivation is specific to this kind of code. Scientific formulas with exponential or power-law growth are evaluated at parameter values a user chose, and a viscosity model handed a temperature far outside the range its Arrhenius parameters were fitted to, or a radiogenics model evaluated an epoch before its reference time, will overflow. An infinity then propagates through sums and ratios and can emerge as a finite, wrong number several steps later. A NaN cannot: it contaminates everything downstream of it and shows up in the output, which is what makes the regime visible.

The tradeoff is deliberate. Returning NaN loses the information that the true answer was large rather than undefined, and a caller that wants to distinguish the two has to check its inputs before calling. In exchange, no silently wrong finite number leaves the function.

## Where these are used

The viscosity models guard their Arrhenius exponentials, the partial-melt models guard their power laws, and the radiogenic isotope model guards its decay exponential. Each module's page notes where a NaN can appear and what it means there. See [Viscosity Models](../viscosity_x/viscosity_models.md), [Partial-Melt Models](../partial_melt_x/partial_melt_models.md), and [Radiogenic Models](../radiogenics_x/radiogenics_models.md).
