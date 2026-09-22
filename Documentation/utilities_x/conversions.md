# Conversions and Scales (`Utilities_x.conversions`, `Utilities_x.dimensions`)

_Updated: 2026-09-16_

`conversions` converts between MKS and the units typically used in the literature, and between orbital elements related by Kepler's third law. `dimensions` builds the scale factors that turn a dimensional interior problem into a non-dimensional one, so the non-dimensionalization the solvers rely on is defined in one place.

## Unit Conversions

```python
from TidalPy.Utilities_x.conversions import (
    Au2m, m2Au, days2rads, rads2days, myr2sec, sec2myr)

Au2m(1.0)          # 1.495978707e11  [m]
m2Au(1.495978707e11)   # 1.0          [AU]

days2rads(1.0)     # 7.2722e-05      [rad s-1] from an orbital or rotational period in days
rads2days(7.2722e-05)  # 1.0          [days]

myr2sec(1.0)       # 3.15576e13      [s]
sec2myr(3.15576e13) # 1.0            [Myr]
```

The period conversions are the pair used most often, because orbital and rotational periods are quoted in days while every TidalPy argument named a frequency is an angular frequency in rad s$^{-1}$.

The mega-year conversion uses the Julian year of 365.25 days, so one Myr is exactly $3.15576 \times 10^{13}$ s. That is the same constant the radiogenics datasets use for their half lives, available in Python as `TidalPy.constants.seconds_per_myr`, so times converted here and times read from an isotope model agree.

## Orbital Elements

```python
from TidalPy.Utilities_x.conversions import orbital_motion2semi_a, semi_a2orbital_motion

semi_major_axis   = orbital_motion2semi_a(orbital_motion, host_mass, target_mass=0.0)
orbital_motion    = semi_a2orbital_motion(semi_major_axis, host_mass, target_mass=0.0)
```

Both are Kepler's third law, $n^2 a^3 = G (M_{\text{host}} + M_{\text{target}})$, solved for one element or the other. The orbiting body's mass defaults to zero, which is the test-particle limit and is usually what a satellite problem wants; supply it when the mass ratio is large enough to matter, as in a binary or a planet-moon pair like Pluto and Charon.

Both take an optional `G_to_use` so a comparison against a published result can use whatever value of the gravitational constant that work adopted. The default is TidalPy's own, which comes from SciPy at initialization.

A non-positive host mass or a negative target mass raises `ValueError`.

## Non-Dimensionalization

The radial structure and deformation problems are integrated in non-dimensional variables. A radius near $10^7$, a density near $10^3$, a modulus near $10^{11}$, and a gravitational constant near $10^{-11}$ put the entries of one linear system thirty orders of magnitude apart, and the solution loses significant digits to that spread. Scaling each variable by a characteristic value of its own dimension brings the system to order unity.

The scales are built from just two properties of the body, its mean radius and its bulk density, plus the gravitational constant:

$$T^2 = \frac{1}{\pi G \bar{\rho}}, \qquad L = R, \qquad \rho_{\text{scale}} = \bar{\rho}$$

with the mass and pressure scales following as $\rho_{\text{scale}} L^3$ and $M_{\text{scale}} / (L T^2)$.

```python
from TidalPy.Utilities_x.dimensions import build_nondimensional_scales

scales = build_nondimensional_scales(mean_radius, bulk_density)

scales.length_conversion     # [m]        the body's mean radius
scales.length3_conversion    # [m^3]
scales.second_conversion     # [s]        sqrt(1 / (pi G rho_bulk))
scales.second2_conversion    # [s^2]
scales.density_conversion    # [kg m-3]   the bulk density
scales.mass_conversion       # [kg]
scales.pascal_conversion     # [Pa]
```

Each attribute is the factor a non-dimensional value is multiplied by to recover MKS, and divided by to go the other way. The returned object is a thin wrapper over the C++ `c_NonDimensionalScales` struct, which is what the solvers hold internally.

Callers rarely build these by hand. The radial solver and the world equation-of-state solve non-dimensionalize their inputs, integrate, and re-dimensionalize their results before returning, so the scales are an implementation detail unless you are reading solver internals or writing a new solver stage.

## C++ API

The conversion functions have `cf_` prefixed `nogil` Cython counterparts (`cy_orbital_motion2semi_a` and the rest) for use from `cdef` code without Python overhead. The non-dimensional scales live in `nondimensional_.hpp` as `c_NonDimensionalScales`, populated by `cy_build_nondimensional_scales`, and are passed by reference into the solvers that need them.
