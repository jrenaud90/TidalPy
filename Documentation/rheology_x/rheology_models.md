# Rheology Models (`rheology_x`)

`TidalPy.rheology_x` provides the C++ rheology model hierarchy that maps a
material's static (background) mechanical properties to a frequency-dependent
**complex (shear/bulk) modulus** `μ*` [Pa]. This is the single quantity the
models compute and expose (`calc_complex_modulus`); the imaginary part controls
tidal dissipation. Element compliances appear only as internal intermediates for
the series composites (see below).

Each `PhysicsLayer` (and its subclasses) can hold up to two rheology models, one
for the shear response and one for the bulk response. `PhysicsLayer.calc_complex_shear_modulus(frequency)`
and `calc_complex_bulk_modulus(frequency)` delegate to the attached models'
`calc_complex_modulus`.

## Inheritance

```
c_TidalPyBaseClass
  └── c_PhysicsBase
        └── c_RheologyBase        (abstract)
              ├── c_Elastic       alias "off"
              ├── c_Viscous       alias "newton"
              ├── c_Voigt         alias "voigt-kelvin"
              ├── c_Maxwell
              ├── c_Burgers
              ├── c_Andrade
              └── c_Sundberg      alias "sundberg-cooper"
```

The abstract base declares `calc_complex_modulus(modulus, viscosity, frequency)`
as pure virtual; every model returns the complex modulus `μ*` [Pa]
directly. (Argument order is **modulus, viscosity, frequency** — modulus first.) Simple models
(Elastic, Viscous, Maxwell, Voigt) are evaluated analytically.
The series composites (Burgers, Andrade, Sundberg) combine their
constituent elements *in series*, which means the element **compliances** add
and the resulting modulus is `μ* = 1 / Σ Jᵢ`. Those element compliances are kept
as internal intermediates and are never exposed to Python.

## The seven models

All inputs are in MKS units: viscosity `η` [Pa·s], unrelaxed modulus `μ` [Pa],
and forcing frequency `ω` [rad/s]. The complex modulus has a positive storage
(real) part and a non-negative loss (imaginary) part.

| Model | Complex modulus `μ*` [Pa] | Parameters |
|-------|---------------------------|------------|
| `Elastic` (`off`) | `μ` (real) | — |
| `Viscous` (`newton`) | `i η ω` | — |
| `Maxwell` | `1 / J_maxwell` | — |
| `Voigt` (`voigt-kelvin`) | `1 / J_voigt` = `μ·f_J + i ω η_v` | `voigt_modulus_frac`, `voigt_viscosity_frac` |
| `Burgers` | `1 / (J_maxwell + J_voigt)` | `voigt_modulus_frac`, `voigt_viscosity_frac` |
| `Andrade` | `1 / J_andrade` | `alpha`, `zeta` |
| `Sundberg` (`sundberg-cooper`) | `1 / (J_andrade + J_voigt)` | `alpha`, `zeta`, `voigt_modulus_frac`, `voigt_viscosity_frac` |

The internal element compliances (mirroring TidalPy's validated legacy
`rheology.complex_compliance` math) are:

```
J_maxwell = 1/μ − i / (η ω)
J_voigt   = J_v / (1 + (J_v η_v ω)²) − i J_v² η_v ω / (1 + (J_v η_v ω)²)
J_andrade = J_maxwell + (1/μ)·(J ω η ζ)^{−α}·Γ(1+α)·[cos(απ/2) − i sin(απ/2)]
```

where `f_J = voigt_modulus_frac`, `J_v = (1/μ) / f_J` is the Voigt arm's
compliance (the layer compliance divided by the modulus fraction),
`η_v = voigt_viscosity_frac · η` is its viscosity, `J = 1/μ`, and
`Γ` is the gamma function. Defaults are `voigt_modulus_frac = 5.0`,
`voigt_viscosity_frac = 0.02`, `alpha = 0.3`, `zeta = 1.0`. The Andrade and
Sundberg families assume a positive forcing frequency.

## Usage

```python
from TidalPy.rheology_x import Maxwell, Andrade, make_rheology

viscosity = 1.0e20    # Pa·s
modulus   = 5.0e10    # Pa
frequency = 1.0e-5    # rad/s

m = Maxwell()
# complex shear modulus [Pa]
G = m.calc_complex_modulus(modulus, viscosity, frequency)

# Parameterised model
a = Andrade(alpha=0.25, zeta=2.0)

# Name/alias-based factory (case-insensitive), config keys are model-dependent
s = make_rheology("Sundberg-Cooper", {"alpha": 0.4, "zeta": 2.0})
```

`make_rheology(model_name, config=None)` resolves aliases case-insensitively and
forwards the relevant keys from `config` (`alpha`, `zeta`, `voigt_modulus_frac`,
`voigt_viscosity_frac`). Unknown names raise `ValueError`.

### Factory internals

At the C++ level the factory is enum-based. `c_RheologyModel` names one value per
model, `c_rheology_model_from_name(name)` maps a (case-insensitive) name or alias
to that enum (throwing `std::invalid_argument` on an unknown name), and
`c_find_rheology(model, config)` returns a `std::unique_ptr<c_RheologyBase>` to a
freshly heap-allocated rheology model. This is the canonical C++ factory used by
all C++ consumers (layers attaching rheology, future binary reconstruction). The
Python `make_rheology` wraps it: it builds the `c_RheologyConfig`, calls
`c_rheology_model_from_name` then `c_find_rheology`, and adopts the returned
`unique_ptr` into the matching rich Python wrapper (`Maxwell`, `Voigt`, ...).

## Vectorized evaluation

Each model also provides vectorized complex-modulus methods (defined once on
`c_RheologyBase`, so all models inherit them via virtual dispatch):

- `calc_complex_modulus_vectorize_modulus(modulus[], viscosity[], frequency)` —
  element-wise over equal-length modulus/viscosity arrays at one frequency.
- `calc_complex_modulus_vectorize_frequency(modulus, viscosity, frequency[])` —
  a frequency sweep at constant modulus and viscosity.
- `calc_complex_modulus_vectorize_all(modulus[], viscosity[], frequency[])` —
  element-wise over three equal-length arrays.

At the C++ level each fills a caller-supplied `std::vector<std::complex<double>>&`
output; mismatched input lengths throw `std::invalid_argument`. The Cython
wrappers accept array-likes and return a `complex128` NumPy array.

## Direct convenience functions

For one-shot evaluation without explicitly constructing a model object, each
model has a lower-case module function:

```python
from TidalPy.rheology_x import maxwell, andrade
import numpy as np

# Scalar in -> Python complex out.
G = maxwell(omega, mu, eta)

# Arrays in -> complex128 ndarray out (frequency, modulus, viscosity may each be
# a float or an ndarray; they are broadcast together).
G_profile = maxwell(omega, np.array([1e10, 5e10]), np.array([1e19, 1e20]))
G_sweep   = andrade(np.logspace(-7, -4, 50), mu, eta, alpha=0.3, zeta=1.0)
```

Signatures: `elastic/viscous/maxwell(frequency, modulus, viscosity)`,
`voigt/burgers(frequency, modulus, viscosity, voigt_modulus_frac=5.0, voigt_viscosity_frac=0.02)`,
`andrade(frequency, modulus, viscosity, alpha=0.3, zeta=1.0)`,
`sundberg(frequency, modulus, viscosity, alpha=0.3, zeta=1.0, voigt_modulus_frac=5.0, voigt_viscosity_frac=0.02)`.

Each builds a *stack-allocated* C++ model, solves (picking the most specific
vectorized routine for the input pattern), and returns. `frequency`, `modulus`
and `viscosity` accept floats or NumPy arrays; the model parameters
(`alpha`, `zeta`, `voigt_*`) are always scalar constants.

> **Note on argument order:** these convenience functions take `frequency` first
> (`func(frequency, modulus, viscosity, ...)`), whereas the class methods
> `calc_complex_modulus*` take `modulus` first (`(modulus, viscosity, frequency)`).

## Serialization

Every model supports the standard TidalPy interfaces:

- `get_config_dict()` → dict of `model_name` plus any model parameters.
- `save_config(path)` → TOML config file.
- `save_binary(path)` / `load_binary(path, force=False)` → TidalPy binary format
  (preserves the model name and parameters; the layer observer pointer is not
  serialized).

## Adding a new rheology model

The hierarchy is designed so a new model is a small, local addition. To add a
new rheology model named `Foo`:

**C++ (`TidalPy/rheology_x/rheology_.hpp`)**

1. If `Foo` needs new parameters, add them to `c_RheologyConfig` (with sensible
   defaults). The single combined config is shared by all models.
2. If the constitutive law is a new series combination, add an internal
   `detail::element_compliance_*` helper and/or a `rheo_modulus_foo(...)` free
   function. Simple models can compute the modulus analytically inline.
3. Add the model class `c_Foo : public c_RheologyBase`:
   - constructors `c_Foo()` and `explicit c_Foo(const c_RheologyConfig&)` that
     pass a model-name string to the base and copy any params into `p_*` members;
   - `get_*` accessors for each param;
   - override `calc_complex_modulus(modulus, viscosity, frequency)` returning the
     complex modulus [Pa] (modulus-first argument order);
   - override `write_binary` / `read_binary` using the shared base helpers
     (`this->write_physics_binary(out, class_id, {params...})` and
     `this->read_physics_binary(in, force, n_params)` from `c_PhysicsBase`).

**C++ (`TidalPy/Utilities_x/binary_x/binary_.hpp`)**

4. Add a unique `BinaryClassID::Foo` value (rheology models occupy 30X).

**C++ factory (`rheology_.hpp`)**

5. Add `Foo` to the `c_RheologyModel` enum, map its name/aliases in
   `c_rheology_model_from_name`, and add a `case` to `c_find_rheology`.
6. Add a `case BinaryClassID::Foo` to `c_rheology_from_binary` (the binary-dispatch
   factory) so a `Foo` attached to a layer can be reconstructed when the layer is
   loaded from a binary file. Forgetting this makes recursive layer loads throw
   "unknown rheology class id in binary stream".

**Cython (`rheology.pxd` / `rheology.pyx`)**

7. Declare `c_Foo` (constructors + getters) in `rheology.pxd` and add the enum
   value to the `c_RheologyModel` cimport.
8. Add the `cdef class Foo(RheologyBase)` wrapper in `rheology.pyx` (with param
   properties and a `get_config_dict` override), the adoption branch in
   `make_rheology`, and the lower-case `foo(frequency, modulus, viscosity, ...)`
   convenience function.

**Package + tests + docs**

9. Export `Foo` (and `foo`) from `TidalPy/rheology_x/__init__.py` (and the C++
   names from `__init__.pxd`).
10. Add `Foo` to the parametrized test lists in
    `Tests/Test_Rheology_x/test_rheology_01.py` (model name, modulus cross-check
    against an independent reference, factory/alias, vectorization, config dict,
    binary round-trip, isinstance). When a model can be attached to a layer, also
    cover the layer recursive-binary round-trip.
11. Document the model here (formula, parameters, references).

No build-system change is needed — `rheology_x.rheology` is already registered in
`cython_extensions.json`.

## References

- Henning, O'Connell, and Sasselov (2009), ApJ, [DOI: 10.1088/0004-637X/707/2/1000](https://doi.org/10.1088/0004-637X/707/2/1000) — Maxwell, Voigt-Kelvin, Burgers.
- Efroimsky (2012), ApJ, [DOI: 10.1088/0004-637X/746/2/150](https://doi.org/10.1088/0004-637X/746/2/150) — complex compliances and Love numbers.
- Renaud and Henning (2018), ApJ, [DOI: 10.3847/1538-4357/aab784](https://doi.org/10.3847/1538-4357/aab784) — Andrade and Sundberg-Cooper.
