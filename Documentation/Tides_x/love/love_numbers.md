# Love Numbers


`TidalPy.Tides_x.love`: the Love-number storage type, the names of the Love-number solution methods, and
the closed-form homogeneous-sphere Love numbers. Full computation from a layered interior is performed by
the radial solver (`RadialSolver_x` or `LayeredWorld.solve_love_numbers`), which populates these values;
`LayeredWorld.solve_love_numbers(love_method=...)` can also use the homogeneous-sphere formulas below.

## Overview

Tidal Love numbers describe how a planetary body deforms under an external tidal potential. Three numbers characterize the elastic response:

| Symbol | Name | Description |
|--------|------|-------------|
| **k** | Potential Love number | Change in external gravitational potential due to mass redistribution |
| **h** | Radial displacement Love number | Amplitude of vertical (radial) surface deformation |
| **l** | Tangential displacement Love number | Amplitude of horizontal (tangential) surface deformation |

All three are dimensionless complex numbers. The **real part** is the elastic amplitude; the **imaginary part** represents energy dissipation at the tidal forcing frequency.

## C++ struct — `c_LoveNumbers`

Defined in `TidalPy/Tides_x/love/love_.hpp` within namespace `tidalpy`.

```cpp
struct c_LoveNumbers {
    std::complex<double> k = {0.0, 0.0};
    std::complex<double> h = {0.0, 0.0};
    std::complex<double> l = {0.0, 0.0};

    c_LoveNumbers() noexcept = default;
    c_LoveNumbers(std::complex<double> k_in,
                  std::complex<double> h_in,
                  std::complex<double> l_in) noexcept;

    bool operator==(const c_LoveNumbers& o) const noexcept;
    bool operator!=(const c_LoveNumbers& o) const noexcept;
};
```

Stack-allocated in Cython; heap-allocated as a member of `c_PhysicsLayer` and `c_SolidLiquidLayer`.

## Python class — `LoveNumbers`

```python
from TidalPy.Tides_x.love import LoveNumbers

ln = LoveNumbers(k=0.3 - 0.01j, h=0.6 - 0.02j, l=0.1 - 0.005j)

print(ln.k)   # (0.3-0.01j)
print(ln.h)   # (0.6-0.02j)
print(ln.l)   # (0.1-0.005j)

# Tuple unpacking
k, h, l = ln

# Serialization
d = ln.to_dict()
# {'love_number_k_re': 0.3, 'love_number_k_im': -0.01,
#  'love_number_h_re': 0.6, 'love_number_h_im': -0.02,
#  'love_number_l_re': 0.1, 'love_number_l_im': -0.005}
```

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `k` | `complex` | Potential Love number |
| `h` | `complex` | Radial displacement Love number |
| `l` | `complex` | Tangential displacement Love number |

### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `to_dict()` | `dict` | Six-key flat dict with `_re`/`_im` suffixes |
| `__repr__()` | `str` | Human-readable representation |
| `__eq__(other)` | `bool` | Equality check (component-wise) |
| `__iter__()` | iterator | Yields `k`, `h`, `l` for tuple unpacking |

## Integration with Layer Classes

`PhysicsLayer` and `SolidLiquidLayer` store Love numbers internally as a `c_LoveNumbers` struct:

```python
from TidalPy.structures_x.layers.physics import PhysicsLayer

pl = PhysicsLayer("mantle", 1, 3.485e6, 6.371e6, 4.043e24,
                  love_number_k=0.3 - 0.01j,
                  love_number_h=0.6 - 0.02j,
                  love_number_l=0.1 - 0.005j)

# Access as a LoveNumbers object
ln = pl.love_numbers

# Or access individual components
k = pl.love_number_k
h = pl.love_number_h
l = pl.love_number_l
```

## Binary Serialization

Love numbers are serialized in `c_PhysicsLayer::write_binary` as six consecutive `double` values (re, im for each of k, h, l), contributing `6 × 8 = 48 bytes` to the payload.

## Computing Love Numbers

Solved Love numbers come from the radial solver (`RadialSolver_x.radial_solver` for the
standalone array API, or `LayeredWorld.solve_love_numbers` on a built world). The
`c_LoveNumbers` struct is the storage target: the solver populates the three fields
after integration.

## Love-number methods

A world obtains its Love numbers by one of these methods (`LayeredWorld.solve_love_numbers(love_method=...)`
per call, `set_tide_config(love_method=...)` or the `[tides]` key `love_method` for the default that
`calc_tides` uses). `love_method_name(alias)` returns the canonical name.

| Canonical name | Aliases | Source of k, h, l |
|---|---|---|
| `radial_solver` | `shooting`, `rs` | Shooting-method radial solve (default). |
| `propagation_matrix` | `prop_matrix`, `pm`, `prop` | Propagation-matrix radial solve (single solid, static, incompressible layer). |
| `homogeneous` | `homogen` | Homogeneous-sphere formulas with the volume-averaged complex shear modulus of the tidal layers. |
| `cpl` | | Homogeneous-sphere formulas on the static modulus, then `(1 - i/Q)`. |
| `ctl` | | Homogeneous-sphere formulas on the static modulus, then `(1 - i omega dt)`. |
| `laterally_inhomogeneous` | `3d`, `lat_inhom` | Reserved for the 3D Love solver (`NotImplementedError`). |

## Homogeneous-sphere formulas

For a homogeneous incompressible sphere of radius `R`, bulk density `rho`, surface gravity `g`, and shear
modulus `mu` (Love 1911; Munk & MacDonald 1960):

```
mu_eff_l = (2 l^2 + 4 l + 3) / l * mu / (rho g R)          # 19/2 * mu / (rho g R) at l = 2
k_l = 3 / (2 (l - 1))      / (1 + mu_eff_l)
h_l = (2 l + 1) / (2 (l - 1)) / (1 + mu_eff_l)
l_l = 3 / (2 l (l - 1))    / (1 + mu_eff_l)
```

In the fluid limit (`mu -> 0`) these give `k_2 = 3/2`, `h_2 = 5/2`, `l_2 = 3/4`. With a complex, frequency
dependent shear modulus from a rheology model the Love numbers are complex; with the real static modulus they
are the static (elastic) Love numbers.

```python
from TidalPy.rheology_x import Maxwell
from TidalPy.Tides_x.love import (
    apply_fixed_dt, apply_fixed_q, calc_effective_rigidity, calc_homogeneous_love_numbers)

mu = Maxwell().calc_complex_modulus(60.0e9, 1.0e19, 1.0e-5)          # complex shear modulus [Pa]
love = calc_homogeneous_love_numbers(mu, 4000.0, 6.7, 6.0e6)          # LoveNumbers(k, h, l), degree 2
mu_eff = calc_effective_rigidity(60.0e9, 4000.0, 6.7, 6.0e6, degree_l=2)

static = calc_homogeneous_love_numbers(60.0e9, 4000.0, 6.7, 6.0e6)
cpl = apply_fixed_q(static, 50.0)            # k, h, l times (1 - i/50):     -Im[k] = Re[k] / 50
ctl = apply_fixed_dt(static, 1.0e-5, 600.0)  # k, h, l times (1 - i omega dt)
```

| Function | Returns |
|---|---|
| `calc_effective_rigidity(shear_modulus, density, gravity, radius, degree_l=2)` | `mu_eff_l` (complex for a complex modulus). |
| `calc_homogeneous_love_numbers(complex_shear_modulus, density, gravity, radius, degree_l=2)` | `LoveNumbers` k, h, l. |
| `apply_fixed_q(love_numbers, fixed_q)` | Constant phase lag applied to k, h, and l. |
| `apply_fixed_dt(love_numbers, frequency, fixed_dt)` | Constant time lag applied to k, h, and l. |
| `love_method_name(method)` | Canonical method name for an alias. |

The classic `TidalPy.tides.love1d` helpers (`complex_love`, `static_love`, `effective_rigidity`, and their
`_general` forms) map onto these: `effective_rigidity_general` in the classic code evaluated
`2 l^2 + 4 l + 3/l` (a precedence slip; correct at `l = 2` only through the dedicated degree-2 function), the
new functions use `(2 l^2 + 4 l + 3)/l` at every degree. Numerical Love numbers of a homogeneous sphere on
a radial grid are still available from `TidalPy.RadialSolver_x.homogeneous_love_numbers`.
