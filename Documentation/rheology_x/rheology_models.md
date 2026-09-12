# Rheology Models (`rheology_x`)

_Updated: 2026-09-12_

A rheology model maps a material's static (purely real) mechanical properties onto a **complex modulus** $\mu^*(\omega)$ [Pa] at a given forcing frequency. The real part is the storage modulus, the part of the stress in phase with the strain; the imaginary part is the loss modulus, the part in quadrature, and it is what converts mechanical work into heat. Their ratio $\mathrm{Im}[\mu^*]/\mathrm{Re}[\mu^*]$ is the material's loss tangent, the inverse of its quality factor $Q$.

Everything on this page applies equally to the shear and the bulk response. The models do not know which one they are computing; supply a shear modulus with a shear viscosity, or a bulk modulus with a bulk viscosity, and the same constitutive law applies. In practice the shear response dominates tidal dissipation in solid bodies, and the bulk response is usually left elastic.

## What the models compute

Each model implements one method, `calc_complex_modulus(modulus, viscosity, frequency)`, which takes the unrelaxed modulus $\mu$ [Pa], the viscosity $\eta$ [Pa s], and the forcing frequency $\omega$ [rad s$^{-1}$], and returns $\mu^*$ [Pa]. The argument order is modulus first throughout the module, including the convenience functions and the vectorized variants.

The single scale that organizes every model is the **Maxwell time** $\tau = \eta / \mu$, the time a material takes to relax an applied stress by viscous flow. Forcing much faster than $\tau$ finds the material effectively elastic; forcing much slower finds it effectively fluid. Dissipation is largest in between, and the models differ mostly in how broad that window is and in what happens on its high-frequency side.

Two conventions matter when reading results. First, the returned imaginary part is non-negative for a positive forcing frequency: the models use the $e^{+i\omega t}$ convention, so a lagging response carries a positive imaginary modulus. Second, the models are evaluated at a single frequency with no memory of previous calls, so they are safe to use across threads and in any order.

## Inheritance

```
c_TidalPyBaseClass
  └── c_PhysicsBase
        └── c_RheologyBase  (abstract)
              ├── c_Elastic       alias "off"
              ├── c_Viscous       alias "newton"
              ├── c_Voigt         aliases "voigt-kelvin", "voigt_kelvin"
              ├── c_Maxwell
              ├── c_Burgers
              ├── c_Andrade
              └── c_Sundberg      aliases "sundberg-cooper", "sundberg_cooper"
```

`c_RheologyBase` declares `calc_complex_modulus` pure virtual and supplies the vectorized loops, the configuration export, and the binary encoding that every model inherits. The Python classes (`Elastic`, `Viscous`, `Voigt`, `Maxwell`, `Burgers`, `Andrade`, `Sundberg`) are thin Cython wrappers holding a pointer to the C++ object, and `RheologyBase` is their shared Python base.

## The seven models

Simple models (Elastic, Viscous, Maxwell, Voigt) are evaluated in closed form. The composites (Burgers, Andrade, Sundberg) place elements **in series**, which means their **compliances** add and the modulus is the reciprocal of the sum, $\mu^* = 1 / \sum_i J_i$. Those element compliances are internal intermediates and are not exposed.

| Model | Complex modulus $\mu^*$ [Pa] | Parameters | Character |
|---|---|---|---|
| `Elastic` (`off`) | $\mu$ | — | No dissipation at any frequency. |
| `Viscous` (`newton`) | $i \eta \omega$ | — | No stored energy; pure loss. |
| `Maxwell` | $1 / J_\mathrm{maxwell}$ | — | One relaxation peak at $\omega\tau = 1$; loss falls as $\omega^{-1}$ above it. |
| `Voigt` (`voigt-kelvin`) | $1 / J_\mathrm{voigt} = \mu f_J + i \omega \eta_v$ | `voigt_modulus_frac`, `voigt_viscosity_frac` | Stiffens without limit at high frequency; rarely used alone. |
| `Burgers` | $1 / (J_\mathrm{maxwell} + J_\mathrm{voigt})$ | `voigt_modulus_frac`, `voigt_viscosity_frac` | Maxwell plus a secondary peak from the Voigt arm. |
| `Andrade` | $1 / J_\mathrm{andrade}$ | `alpha`, `zeta` | Maxwell plus a transient term; loss falls only as $\omega^{-\alpha}$. |
| `Sundberg` (`sundberg-cooper`) | $1 / (J_\mathrm{andrade} + J_\mathrm{voigt})$ | `alpha`, `zeta`, `voigt_modulus_frac`, `voigt_viscosity_frac` | Andrade's high-frequency tail plus Burgers' secondary peak. |

The element compliances, which mirror the math of TidalPy's validated classic `rheology.complex_compliance` module, are

$$J_\mathrm{maxwell} = \frac{1}{\mu} - \frac{i}{\eta \omega}$$

$$J_\mathrm{voigt} = \frac{J_v}{1 + (J_v \eta_v \omega)^2} - i \frac{J_v^2 \eta_v \omega}{1 + (J_v \eta_v \omega)^2}$$

$$J_\mathrm{andrade} = J_\mathrm{maxwell} + \frac{1}{\mu} \left( \frac{\eta \omega \zeta}{\mu} \right)^{-\alpha} \Gamma(1 + \alpha) \left[ \cos\frac{\pi\alpha}{2} - i \sin\frac{\pi\alpha}{2} \right]$$

where $f_J$ is `voigt_modulus_frac`, $J_v = (1/\mu) / f_J$ is the Voigt arm's compliance, $\eta_v = $ `voigt_viscosity_frac` $\times\, \eta$ is its viscosity, and $\Gamma$ is the gamma function.

| Parameter | Default | Meaning |
|---|---|---|
| `alpha` | 0.3 | Andrade exponent. Sets how slowly the transient loss decays with frequency; laboratory values for silicates cluster near 0.2 to 0.4. |
| `zeta` | 1.0 | Ratio of the Andrade timescale to the Maxwell time. Larger values push the transient term down. |
| `voigt_modulus_frac` | 5.0 | Stiffness of the Voigt arm's spring relative to the main spring. The arm's compliance is the material compliance divided by this value. |
| `voigt_viscosity_frac` | 0.02 | Viscosity of the Voigt arm's dashpot as a fraction of the material viscosity. |

### Behavior at the limits

Two edge cases are worth knowing before they surprise you.

At zero frequency Maxwell, Burgers, Andrade, and Sundberg return effectively zero: with unlimited time to flow, a viscoelastic body supports no static rigidity. Elastic returns $\mu$, Viscous returns zero, and Voigt returns $\mu f_J$.

At **negative** frequency Elastic, Viscous, Voigt, Maxwell, and Burgers simply mirror the imaginary part, but Andrade and Sundberg return `NaN`: their transient term raises a negative quantity to a fractional power. Always pass the absolute value of the forcing frequency. TidalPy's own tidal solvers do this for you; a direct call does not.

## Choosing a model

`Elastic` is the right choice when you want deformation without dissipation, for example to isolate the elastic part of a Love number or to model a layer you believe is effectively rigid on the forcing timescale. It is also what a `PhysicsLayer` behaves like when no rheology is attached.

`Maxwell` is the standard first choice and the one to reach for when comparing against published Love numbers, since most of the literature uses it. Its weakness is the high-frequency tail: dissipation falls as $\omega^{-1}$, which underestimates the response of real silicates to fast forcing.

`Andrade` and `Sundberg` are the models to use when the forcing is fast compared with the Maxwell time, which is the usual situation for a cool, stiff, or rapidly forced body. Their loss falls only as $\omega^{-\alpha}$, and for tidal problems that difference can be orders of magnitude in the heating rate.

`Burgers` and `Voigt` are mainly useful for reproducing published work that used them, or for deliberately placing a secondary relaxation peak at a chosen frequency. `Viscous` exists for completeness and for the fluid limit.

## Building a model

Instantiate a class directly, or resolve one by name.

```python
from TidalPy.rheology_x import Maxwell, Andrade, make_rheology

maxwell_model = Maxwell()
andrade_model = Andrade(alpha=0.25, zeta=2.0)

# Case-insensitive; hyphens and underscores are interchangeable in aliases.
sundberg_model = make_rheology("Sundberg-Cooper", {"alpha": 0.4, "zeta": 2.0})
```

`make_rheology(model_name, config=None)` recognizes every name and alias in the inheritance tree above and reads the keys `alpha`, `zeta`, `voigt_modulus_frac`, and `voigt_viscosity_frac` from `config`. Keys a model does not use are ignored, keys a model does use but that are absent fall back to its default, and an unrecognized model name raises `ValueError`.

Model parameters are fixed at construction and exposed as read-only properties (`andrade_model.alpha`, `sundberg_model.voigt_viscosity_frac`). To change one, build a new model.

### Factory internals

At the C++ level the factory is enum-based. `c_RheologyModel` names one value per model; `c_rheology_model_from_name(name)` maps a case-insensitive name or alias onto that enum, throwing `std::invalid_argument` for an unknown name; and `c_find_rheology(model, config)` returns a `std::unique_ptr<c_RheologyBase>` to a freshly heap-allocated model. Every C++ consumer uses this path, including layers attaching a rheology and the binary loader rebuilding one. The Python `make_rheology` wraps it: it fills a `c_RheologyConfig`, calls the two C++ functions, and adopts the returned pointer into the matching Python wrapper.

## Evaluating a model

The scalar call takes three floats and returns a Python `complex`.

```python
from TidalPy.rheology_x import Maxwell

modulus   = 50.0e9    # Pa
viscosity = 1.0e20    # Pa s
frequency = 1.0e-5    # rad s-1

model = Maxwell()
complex_shear = model.calc_complex_modulus(modulus, viscosity, frequency)
```

Three vectorized methods, defined once on the base class so every model inherits them, cover the array patterns that actually occur in a tidal calculation.

| Method | Vectorized over | Typical use |
|---|---|---|
| `calc_complex_modulus_vectorize_modulus(modulus[], viscosity[], frequency)` | Equal-length modulus and viscosity arrays at one frequency | A radial profile of a planet at one forcing frequency. |
| `calc_complex_modulus_vectorize_frequency(modulus, viscosity, frequency[])` | Frequency array at one modulus and viscosity | A frequency sweep of one material. |
| `calc_complex_modulus_vectorize_all(modulus[], viscosity[], frequency[])` | Three equal-length arrays, element by element | Paired inputs already zipped by the caller. |

```python
import numpy as np
from TidalPy.rheology_x import Andrade

model = Andrade()
radial_moduli     = np.linspace(40.0e9, 80.0e9, 100)
radial_viscosity  = np.logspace(18.0, 22.0, 100)

profile = model.calc_complex_modulus_vectorize_modulus(radial_moduli, radial_viscosity, 1.0e-5)
sweep   = model.calc_complex_modulus_vectorize_frequency(50.0e9, 1.0e20, np.logspace(-7, -4, 50))
```

Each fills a caller-supplied `std::vector<std::complex<double>>` at the C++ level; the Cython wrappers accept array-likes and return a `complex128` NumPy array. Mismatched input lengths raise `ValueError`.

## Convenience functions

Each model also has a lower-case free function that builds a stack-allocated C++ model, evaluates it, and returns, with no Python object left behind. These are the fastest way to evaluate a rheology once, and the most convenient way to explore one.

```python
import numpy as np
from TidalPy.rheology_x import maxwell, andrade

# Scalars in, Python complex out.
complex_shear = maxwell(50.0e9, 1.0e20, 1.0e-5)

# Arrays in, complex128 ndarray out. Any of the three may be an array.
profile = maxwell(np.array([1.0e10, 5.0e10]), np.array([1.0e19, 1.0e20]), 1.0e-5)
sweep   = andrade(50.0e9, 1.0e20, np.logspace(-7, -4, 50), alpha=0.3, zeta=1.0)
```

The signatures follow the classes: `elastic/viscous/maxwell(modulus, viscosity, frequency)`, `voigt/burgers(modulus, viscosity, frequency, voigt_modulus_frac=5.0, voigt_viscosity_frac=0.02)`, `andrade(modulus, viscosity, frequency, alpha=0.3, zeta=1.0)`, and `sundberg(modulus, viscosity, frequency, alpha=0.3, zeta=1.0, voigt_modulus_frac=5.0, voigt_viscosity_frac=0.02)`. The model parameters are always scalars; `modulus`, `viscosity`, and `frequency` may each be a float or an array and are broadcast together, with the most specific vectorized routine chosen for the pattern supplied.

## Attaching a rheology to a layer

A rheology becomes part of a planet when it is attached to a layer.

```python
from TidalPy.rheology_x import Maxwell, make_rheology
from TidalPy.structures_x.layers.physics import PhysicsLayer

mantle = PhysicsLayer("mantle", 0, 0.0, 1.0e6, 2.1e19,
                      shear_modulus_static=50.0e9, bulk_modulus_static=100.0e9,
                      shear_viscosity_static=1.0e20, bulk_viscosity_static=1.0e20)

mantle.set_shear_rheology(Maxwell())
mantle.set_bulk_rheology(make_rheology("andrade", {"alpha": 0.3}))

complex_shear = mantle.calc_complex_shear_modulus(1.0e-5)
```

Ownership of the C++ model transfers into the layer: the Python wrapper becomes an empty shell and cannot be attached again. Until a rheology is set, `calc_complex_shear_modulus` returns the static modulus as a purely real complex number, which is elastic behavior. The equivalent declarative form is a `[layers.<name>.shear_rheology]` table in a world's TOML, keyed by `model` plus any parameters; see the [TOML schema](../structures_x/config/toml_schema.md) and [PhysicsLayer](../structures_x/layers/physics_layer.md).

## Serialization

Every model supports the standard TidalPy interfaces.

| Call | Result |
|---|---|
| `get_config_dict()` | A dict of `model` plus the model's own parameters, in the same form the world builder reads. |
| `save_config(path)` | That dict written as a TOML file. |
| `save_binary(path)` | The TidalPy binary format, preserving the model name and parameters. |
| `load_binary(path, force=False)` | Reads a saved model back into an existing instance. `force=True` accepts a file written by a different schema version. |

A rheology attached to a layer is written as part of that layer's binary record and reconstructed recursively when the layer is loaded, so a saved planet round-trips with its rheologies intact. See [Binary serialization](../utilities_x/binary_x.md).

## Adding a new rheology model

The hierarchy is designed so that a new model is a small, local addition. To add one named `Foo`:

**C++ (`TidalPy/rheology_x/rheology_.hpp`)**

1. If `Foo` needs new parameters, add them to `c_RheologyConfig` with sensible defaults. The single combined config is shared by all models.
2. If the constitutive law is a new series combination, add an internal `detail::element_compliance_*` helper. A model with a closed form can compute its modulus inline instead.
3. Add the class `c_Foo : public c_RheologyBase` with constructors `c_Foo()` and `explicit c_Foo(const c_RheologyConfig&)` that pass a model-name string to the base and copy any parameters into `p_*` members, a `get_*` accessor per parameter, an override of `calc_complex_modulus(modulus, viscosity, frequency)` returning the complex modulus, and overrides of `write_binary` / `read_binary` built on the `c_PhysicsBase` helpers.

**C++ (`TidalPy/Utilities_x/binary_x/binary_.hpp`)**

4. Add a unique `BinaryClassID::Foo` value. Rheology models occupy the 30X range.

**C++ factory (`rheology_.hpp`)**

5. Add `Foo` to the `c_RheologyModel` enum, map its name and aliases in `c_rheology_model_from_name`, and add a `case` to `c_find_rheology`.
6. Add a `case BinaryClassID::Foo` to `c_rheology_from_binary` so a `Foo` attached to a layer can be rebuilt when the layer is loaded. Omitting this makes recursive layer loads throw "unknown rheology class id in binary stream".

**Cython (`rheology.pxd` / `rheology.pyx`)**

7. Declare `c_Foo` (constructors and getters) in `rheology.pxd` and add the enum value to the `c_RheologyModel` cimport.
8. Add the `cdef class Foo(RheologyBase)` wrapper in `rheology.pyx` with its parameter properties, the adoption branch in `make_rheology`, and the lower-case `foo(modulus, viscosity, frequency, ...)` convenience function. The config dict comes from the C++ `append_config_entries` override, so no Cython override is needed.

**Package, tests, and documentation**

9. Export `Foo` and `foo` from `TidalPy/rheology_x/__init__.py`, and the C++ names from `__init__.pxd`.
10. Add `Foo` to the parametrized lists in `Tests/Test_Rheology_x/test_rheology_01.py`, which cover the model name, the modulus against an independent reference, the factory and aliases, vectorization, the config dict, the binary round trip, and `isinstance`. If the model can be attached to a layer, also cover the layer's recursive binary round trip.
11. Document the model here with its formula, parameters, and references.

No build-system change is needed; `rheology_x.rheology` is already registered in `cython_extensions.json`.

## References

- Henning, W. G., O'Connell, R. J., and Sasselov, D. D. (2009). Tidally heated terrestrial exoplanets: Viscoelastic response models. *The Astrophysical Journal*, 707(2), 1000-1015. [DOI](https://doi.org/10.1088/0004-637X/707/2/1000) — Maxwell, Voigt-Kelvin, and Burgers.
- Efroimsky, M. (2012). Tidal dissipation compared to seismic dissipation: In small bodies, Earths, and super-Earths. *The Astrophysical Journal*, 746(2), 150. [DOI](https://doi.org/10.1088/0004-637X/746/2/150) — complex compliances and Love numbers.
- Renaud, J. P., and Henning, W. G. (2018). Increased tidal dissipation using advanced rheological models: Implications for Io and tidally active exoplanets. *The Astrophysical Journal*, 857(2), 98. [DOI](https://doi.org/10.3847/1538-4357/aab784) — Andrade and Sundberg-Cooper.
