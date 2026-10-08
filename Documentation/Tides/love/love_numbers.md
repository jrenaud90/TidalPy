# Love Numbers (`Tides.love`)

_Updated: 2026-09-30_

`TidalPy.Tides.love` holds the Love-number storage type, the Love-number method names, and the closed-form homogeneous-sphere Love numbers. Solved Love numbers come from the [radial solver](../../RadialSolver/index.md) (`radial_solver` or `BaseWorld.solve_love_numbers`); the world method can also use the formulas below.

## Overview

Tidal Love numbers describe how a body deforms under an external potential:

| Symbol | Name | Description |
|--------|------|-------------|
| **k** | Potential Love number | Change in external gravitational potential due to mass redistribution |
| **h** | Radial displacement Love number | Amplitude of vertical (radial) surface deformation |
| **l** | Tangential displacement Love number | Amplitude of horizontal (tangential) surface deformation |

All three are dimensionless and complex: the real part is the elastic amplitude, the imaginary part the dissipation at the forcing frequency.

## Love-Number Methods

A world picks its method per call with `BaseWorld.solve_love_numbers(love_method=...)`, or sets the default `calc_tides` uses with `set_tide_config(love_method=...)` or the `[tides]` key `love_method`. `love_method_name(alias)` returns the canonical name.

| Canonical name | Aliases | Source of k, h, l |
|---|---|---|
| `radial_solver` | `shooting`, `rs` | Shooting-method radial solve (default). |
| `propagation_matrix` | `prop_matrix`, `pm`, `prop` | Propagation-matrix radial solve (single solid, static, incompressible layer). |
| `homogeneous` | `homogen` | Quasi-homogeneous: each tidal layer's homogeneous-sphere Love numbers from its averaged complex shear modulus, weighted by the layers' tidal scales (see Quasi-Homogeneous Love Numbers in [Worlds](../../Structures/worlds/worlds.md)). |
| `cpl` | | The same on the static (unrelaxed) averaged moduli, then `(1 - i/Q)`. |
| `ctl` | | The same on the static averaged moduli, then `(1 - i omega dt)`. |
| `laterally_inhomogeneous` | `3d`, `lat_inhom` | Reserved for the 3D Love solver (`NotImplementedError`). |

## Homogeneous-Sphere Formulas

For a homogeneous incompressible sphere of radius `R`, bulk density `rho`, surface gravity `g`, and shear modulus `mu` (Love 1911; Munk and MacDonald 1960):

```text
mu_eff_l = (2 l^2 + 4 l + 3) / l * mu / (rho g R)  # 19/2 * mu / (rho g R) at l = 2
k_l = 3 / (2 (l - 1))      / (1 + mu_eff_l)
h_l = (2 l + 1) / (2 (l - 1)) / (1 + mu_eff_l)
l_l = 3 / (2 l (l - 1))    / (1 + mu_eff_l)
```

In the fluid limit ($\mu \to 0$) these give $k_2 = 3/2$, $h_2 = 5/2$, and $l_2 = 3/4$. A complex shear modulus from a rheology model gives complex Love numbers; the real static modulus gives the elastic ones.

```python
from TidalPy.Rheology import Maxwell
from TidalPy.Tides.love import (
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
| `apply_fixed_q(love_numbers, fixed_q)` | k, h, l with a constant phase lag. |
| `apply_fixed_dt(love_numbers, frequency, fixed_dt)` | k, h, l with a constant time lag. |
| `love_method_name(method)` | Canonical method name for an alias. |

## Python API

### `LoveNumbers`

```python
from TidalPy.Tides.love import LoveNumbers

ln = LoveNumbers(k=0.3 - 0.01j, h=0.6 - 0.02j, l=0.1 - 0.005j)
print(ln.k)   # (0.3-0.01j); likewise ln.h, ln.l (complex)

k, h, l = ln  # Tuple unpacking

d = ln.to_dict()
# {'love_number_k_re': 0.3, 'love_number_k_im': -0.01,
#  'love_number_h_re': 0.6, 'love_number_h_im': -0.02,
#  'love_number_l_re': 0.1, 'love_number_l_im': -0.005}
```

`to_dict()` returns the flat dict above; `==` compares component-wise, and iteration yields `k`, `h`, `l`.

## C++ API

`tidalpy::c_LoveNumbers` (header `TidalPy/Tides/love/love_.hpp`):

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
