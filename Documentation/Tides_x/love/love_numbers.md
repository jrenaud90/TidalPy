# Love Numbers

`TidalPy.Tides_x.love`: the Love-number storage type. Full computation from material
properties is performed by the radial solver (`RadialSolver_x` or
`LayeredWorld.solve_love_numbers`), which populates these values.

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
