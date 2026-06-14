# Material EOS Models (`Material_x.eos.material_eos`)

A **material equation-of-state (EOS) model** returns a material's mass density
[kg/m³] as a function of the local pressure [Pa] (analytic models) or radius [m]
(interpolated model). EOS models are attached to a layer and supply the per-layer
density source consumed by the whole-planet radial EOS solve.

EOS models follow the same class pattern as the rheology / cooling / radiogenics
hierarchies: a `MaterialEOSBase` (a `PhysicsBase`) with concrete subclasses, a
name-based factory (`make_material_eos`), and shared binary serialization.

All quantities are **MKS**. The analytic models are **isothermal** (temperature is
accepted for API uniformity but currently unused; thermal expansion is deferred).

---

## Why density-from-pressure works inline

The radial structure solve integrates over **radius**, carrying pressure as one of
its ODE state variables. The density callback is evaluated at each radial step
*with the current integrated pressure available*, so an analytic `density(P)` model
is evaluated inline, no separate coupled iteration is needed beyond the solver's
existing surface-pressure convergence loop. The interpolated model instead returns
density as a function of radius.

---

## Models

| Class | Alias(es) | Density from | Parameters |
|-------|-----------|--------------|------------|
| `ConstantDensityEOS` | `constant`, `uniform` | — (incompressible) | `reference_density_kg_m3` |
| `BirchMurnaghanEOS` | `bm`, `birch_murnaghan` | pressure | `reference_density_kg_m3`, `reference_bulk_modulus_pa`, `bulk_modulus_derivative` |
| `VinetEOS` | `vinet` | pressure | same as BM |
| `InterpolatedEOS` | `interp`, `interpolate` | radius (table) | `radius_m`, `density_kg_m3` |

### Birch-Murnaghan (3rd order)

With compression `η = ρ/ρ₀ = V₀/V`,

```
P(η) = (3/2)·K₀·(η^(7/3) − η^(5/3))·[ 1 + (3/4)(K₀′−4)(η^(2/3) − 1) ]
```

`calc_density(P)` inverts this monotonic law for `η` (safeguarded Newton with a
bisection fallback) and returns `ρ₀·η`.

**Inversion settings (configurable).** The density-from-pressure inversion is a
safeguarded Newton iteration with two numerical knobs, both config fields so they
can be tuned per material:

- `invert_rtol` — relative convergence tolerance on the compression `η` (default
  `1e-13`).
- `invert_max_iters` — hard iteration cap; a termination safeguard only, since the
  iteration normally converges in well under ~10 steps (default `60`).

```python
bm = BirchMurnaghanEOS(3500.0, 1.3e11, 4.5, invert_rtol=1e-9, invert_max_iters=80)
bm.invert_rtol, bm.invert_max_iters   # (1e-9, 80)
```

Both are carried in `get_config_dict` and the binary round-trip. Their defaults
live in one place — the C++ `c_MaterialEOSConfig` (`d_EOS_INVERT_RTOL` /
`d_EOS_INVERT_MAX_ITERS`); the Python wrappers only override them when a value is
explicitly supplied.

### Vinet

With `x = (V/V₀)^(1/3) = η^(−1/3)`,

```
P(x) = 3·K₀·(1 − x)/x²·exp[ (3/2)(K₀′ − 1)(1 − x) ]
```

inverted the same way as Birch-Murnaghan.

### Interpolated

Linear interpolation of a sorted `(radius_m → density_kg_m3)` table, clamped at the
boundaries. Reproduces the legacy interpolation EOS (e.g. PREM Earth profiles).

---

## Python API

```python
from TidalPy.Material_x.eos import (
    ConstantDensityEOS, BirchMurnaghanEOS, VinetEOS, InterpolatedEOS,
    make_material_eos,
)

bm = BirchMurnaghanEOS(reference_density_kg_m3=3500.0,
                       reference_bulk_modulus_pa=1.3e11,
                       bulk_modulus_derivative=4.5)
rho = bm.calc_density(5.0e10)        # density [kg/m^3] at 50 GPa
bm.reference_bulk_modulus            # 1.3e11

# Name factory (aliases accepted):
eos = make_material_eos("vinet", {"reference_density_kg_m3": 3500.0,
                                  "reference_bulk_modulus_pa": 1.3e11,
                                  "bulk_modulus_derivative": 4.5})

# Interpolated (PREM-style):
prem = InterpolatedEOS(radius_m=[0.0, 1.0e6, 2.0e6],
                       density_kg_m3=[5000.0, 4000.0, 3000.0])
prem.calc_density(0.0, 0.0, 0.5e6)   # 4500.0 (radius interpolation)
```

`calc_density(pressure_pa, temperature_k=0.0, radius_m=0.0)` returns density.
`get_config_dict`, `save_config`, `save_binary`, and `load_binary` are inherited
from the base class. The free functions `birch_murnaghan_pressure(eta, K0, K0p)`
and `vinet_pressure(eta, K0, K0p)` expose the forward pressure laws (useful for
inversion cross-checks).

Binary class ids: `ConstantDensityEOS` 601, `BirchMurnaghanEOS` 602, `VinetEOS`
603, `InterpolatedEOS` 604.

---

## C++ API

The models live in `material_eos_.hpp` (namespace `tidalpy`, header-only). The C++
layer is the canonical one — the Cython wrappers above are thin adapters over it,
and other C++ consumers (layers attaching an EOS, the whole-planet EOS solve, binary
reconstruction) use these types directly.

```cpp
#include "material_eos_.hpp"
using namespace tidalpy;

// Build a config; defaults come from the struct itself (single source of truth).
c_MaterialEOSConfig cfg;
cfg.reference_density_kg_m3   = 3500.0;
cfg.reference_bulk_modulus_pa = 1.3e11;
cfg.bulk_modulus_derivative   = 4.5;

c_BirchMurnaghanEOS bm(cfg);
double rho = bm.calc_density(/*pressure_pa=*/5.0e10, /*temperature_k=*/0.0, /*radius_m=*/0.0);

// Or construct via the enum factory (returns a unique_ptr to the base):
std::unique_ptr<c_MaterialEOSBase> eos =
    c_find_material_eos(c_MaterialEOSModel::Vinet, cfg);
double rho2 = eos->calc_density(5.0e10, 0.0, 0.0);
```

### `c_MaterialEOSConfig`

Combined config struct shared by every model (analytic models read the modulus
fields; the interpolated model reads the table fields). Its member initializers are
the single source of default values for both C++ and Python.

| Field | Used by | Default |
|-------|---------|---------|
| `reference_density_kg_m3` | all | `3500.0` |
| `reference_bulk_modulus_pa` | BM, Vinet | `1.0e11` |
| `bulk_modulus_derivative` | BM, Vinet | `4.0` |
| `invert_rtol` | BM, Vinet | `d_EOS_INVERT_RTOL` (`1e-13`) |
| `invert_max_iters` | BM, Vinet | `d_EOS_INVERT_MAX_ITERS` (`60`) |
| `radius_m`, `density_kg_m3` | Interpolated | empty `std::vector<double>` |

### Classes

All derive from `c_MaterialEOSBase : public c_PhysicsBase` and override
`double calc_density(double pressure_pa, double temperature_k, double radius_m) const`.
Each has a default constructor and an `explicit c_XEOS(const c_MaterialEOSConfig&)`.

| Class | Key accessors |
|-------|---------------|
| `c_ConstantDensityEOS` | `get_reference_density()` |
| `c_BirchMurnaghanEOS` | `get_reference_density()`, `get_reference_bulk_modulus()`, `get_bulk_modulus_derivative()`, `get_invert_rtol()`, `get_invert_max_iters()` |
| `c_VinetEOS` | same as `c_BirchMurnaghanEOS` |
| `c_InterpolatedEOS` | `get_num_points()` |

Binary serialization (`write_binary`/`read_binary`, and `save_binary`/`load_binary`)
is inherited from `c_PhysicsBase` via the shared `write_physics_binary` /
`read_physics_binary` helpers.

### Free functions

| Function | Description |
|----------|-------------|
| `double eos_bm_pressure(eta, K0, K0_prime)` | 3rd-order Birch-Murnaghan pressure [Pa] at compression `η`. |
| `double eos_vinet_pressure(eta, K0, K0_prime)` | Vinet pressure [Pa] at `η`. |
| `double eos_invert_eta(pressure_target_pa, K0, K0_prime, pressure_fn, rtol, max_iters)` | Inverts a monotonic `P(η)` law for `η` (safeguarded Newton + bisection bracket). Reused by both analytic models; pass `eos_bm_pressure` / `eos_vinet_pressure` as `pressure_fn`. |

### Factory

| Function | Description |
|----------|-------------|
| `c_MaterialEOSModel c_material_eos_model_from_name(const std::string&)` | Maps a (case-insensitive) name/alias to the enum; throws `std::invalid_argument` on an unknown name. |
| `std::unique_ptr<c_MaterialEOSBase> c_find_material_eos(c_MaterialEOSModel, const c_MaterialEOSConfig&)` | Switch-dispatch construction (via `std::make_unique`). A name-string overload is also provided. |
| `std::unique_ptr<c_MaterialEOSBase> c_material_eos_from_binary(std::istream&, bool force=false)` | Peeks the binary class id and reconstructs the matching model (used when a layer with an attached EOS is loaded). |

---

## Adding a new EOS model

To add a model `Foo`:

1. **`material_eos_.hpp`** — if the law is analytic, add a free
   `eos_foo_pressure(eta, K0, K0p)` (monotonic in `η`) so the shared
   `eos_invert_eta` inverter can be reused; otherwise compute density directly. Add
   the class `c_FooEOS : public c_MaterialEOSBase` with constructors, getters, an
   override of `calc_density(pressure, temperature, radius)`, and
   `write_binary`/`read_binary` via the shared `c_PhysicsBase` helpers.
2. **`binary_.hpp`** — add a unique `BinaryClassID::FooEOS` (EOS models occupy 60X).
3. **`material_eos_.hpp` factory** — add `Foo` to the `c_MaterialEOSModel` enum, map
   its name/aliases in `c_material_eos_model_from_name`, add a `case` to
   `c_find_material_eos`, and add a `case` to `c_material_eos_from_binary` (so a
   `Foo` attached to a layer can be reconstructed on load).
4. **`material_eos.pxd` / `material_eos.pyx`** — declare `c_FooEOS` and the enum
   value; add the `cdef class FooEOS(MaterialEOSBase)` wrapper (params, getters,
   `get_config_dict`) and an adoption branch in `make_material_eos`.
5. **`__init__.py`** — export `FooEOS`.
6. **Tests** (`Tests/Test_Material_x/test_material_eos_01.py`) — density evaluation,
   an inversion cross-check (analytic models), factory/alias, config dict, binary
   round-trip, isinstance.
7. **Docs** — add the model here with its formula and references.

No build-system change is needed — `Material_x.eos.material_eos` is already
registered in `cython_extensions.json`.

---

## References

- Birch (1947), Phys. Rev. 71, 809 — finite-strain (Birch-Murnaghan) EOS.
- Vinet et al. (1987), J. Geophys. Res. 92, 9319 — universal (Vinet) EOS.
