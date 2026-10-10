# Conversions and Scales (`Utilities.conversions`, `Utilities.dimensions`)

_Updated: 2026-09-30_

`conversions` converts between MKS and the units typically used in the literature, and between orbital elements related by Kepler's third law. `dimensions` builds the scale factors that turn a dimensional interior problem into a non-dimensional one, so the non-dimensionalization the solvers rely on is defined in one place.

## Unit Conversions

```python
from TidalPy.Utilities.conversions import (
    Au2m, m2Au, days2rads, rads2days, myr2sec, sec2myr
)

Au2m(1.0)          # 1.495978707e11  [m]
m2Au(1.495978707e11)   # 1.0          [AU]

days2rads(1.0)     # 7.2722e-05      [rad s-1] from an orbital or rotational period in days
rads2days(7.2722e-05)  # 1.0          [days]

myr2sec(1.0)       # 3.15576e13      [s]
sec2myr(3.15576e13) # 1.0            [Myr]
```

Periods are usually quoted in days, while every TidalPy argument named a frequency is an angular frequency in rad s$^{-1}$; `days2rads` converts between them.

The mega-year conversion uses the Julian year of 365.25 days, so one Myr is exactly $3.15576 \times 10^{13}$ s (`TidalPy.constants.seconds_per_myr`), the same constant the radiogenics datasets use.

## Orbital Elements

```python
from TidalPy.Utilities.conversions import orbital_motion2semi_a, semi_a2orbital_motion

semi_major_axis   = orbital_motion2semi_a(orbital_motion, host_mass, target_mass=0.0)
orbital_motion    = semi_a2orbital_motion(semi_major_axis, host_mass, target_mass=0.0)
```

Both are Kepler's third law, $n^2 a^3 = G (M_{\text{host}} + M_{\text{target}})$, solved for one element or the other. The orbiting body's mass defaults to zero, the test-particle limit; supply it when the mass ratio matters, as for Pluto and Charon.

The optional `G_to_use` sets the gravitational constant, to match a published result. The default, `None`, uses TidalPy's current value (from SciPy), read on every call.

A non-positive orbital motion or semi-major axis (NaN included), a non-positive host mass, a negative target mass, or a non-positive `G_to_use` raises `ValueError`.

## Non-Dimensionalization

The radial structure and deformation problems are integrated in non-dimensional variables. A radius near $10^7$, a density near $10^3$, a modulus near $10^{11}$, and a gravitational constant near $10^{-11}$ put the entries of one linear system thirty orders of magnitude apart, which costs significant digits; scaling each variable by a characteristic value brings the system to order unity. The scales come from the body's mean radius and bulk density plus the gravitational constant:

$$T^2 = \frac{1}{\pi G \bar{\rho}}, \qquad L = R, \qquad \rho_{\text{scale}} = \bar{\rho}$$

with the mass and pressure scales following as $\rho_{\text{scale}} L^3$ and $M_{\text{scale}} / (L T^2)$.

```python
from TidalPy.Utilities.dimensions import build_nondimensional_scales

scales = build_nondimensional_scales(mean_radius, bulk_density)

scales.length_conversion     # [m]        the body's mean radius
scales.length3_conversion    # [m^3]
scales.second_conversion     # [s]        sqrt(1 / (pi G rho_bulk))
scales.second2_conversion    # [s^2]
scales.density_conversion    # [kg m-3]   the bulk density
scales.mass_conversion       # [kg]
scales.pascal_conversion     # [Pa]
```

Multiply a non-dimensional value by the attribute to recover MKS, and divide to go the other way. You rarely need these: the radial solver and the equation-of-state solve convert their inputs and results themselves, so the scales matter only when writing a new solver stage.

## C++ API

The Kepler functions are `c_orbital_motion2semi_a` and `c_semi_a2orbital_motion` in `conversions_.hpp`; they skip the input checks. The scales are `c_NonDimensionalScales` in `nondimensional_.hpp`, constructed from a mean radius and a bulk density.
