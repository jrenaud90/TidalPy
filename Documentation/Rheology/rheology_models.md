# Rheology Models (`Rheology`)

_Updated: 2026-10-02_

A rheology model maps a material's static (purely real) mechanical properties onto a complex modulus $\mu^*(\omega)$ \[Pa\] at a given forcing frequency. The real part is the storage modulus, the part of the stress in phase with the strain; the imaginary part is the loss, and it is what converts mechanical work into frictional heat. Their ratio $\mathrm{Im}[\mu^*]/\mathrm{Re}[\mu^*]$ is the material's loss tangent, the inverse of its quality factor $Q$ (Efroimsky 2012).

Everything on this page applies equally to the shear and the bulk response. The models do not know which one they are computing; supply a shear modulus with a shear viscosity, or a bulk modulus with a bulk viscosity, and the same constitutive law applies. In practice the shear response dominates tidal dissipation in solid bodies, and the bulk response is usually left elastic, however this is a new and active area of research.

## Complex Modulus

Each model implements one method, `calc_complex_modulus(modulus, viscosity, frequency)`, which takes the unrelaxed modulus $\mu$ \[Pa\], the viscosity $\eta$ \[Pa s\], and the forcing frequency $\omega$ \[rad s$^{-1}$\], and returns $\mu^*$ \[Pa\]. The argument order is modulus first throughout the module, including the convenience functions and the vectorized variants.

The scale common to every model is the Maxwell time $\tau = \eta / \mu$, the time a material takes to relax an applied stress by viscous flow. Forcing much faster than $\tau$ finds the material effectively elastic; forcing much slower finds it effectively fluid. Dissipation is largest when the forcing timescale is approximately the Maxwell time; the models differ mostly in how broad that window is and in what happens on its high-frequency side.

Two conventions apply to the results. The returned imaginary part is non-negative for a positive forcing frequency: the models use the $e^{+i\omega t}$ convention, so a lagging response carries a positive imaginary modulus. The models are evaluated at a single frequency with no memory of previous calls, so they are safe to use across threads and in any order.

## Physics

A linear viscoelastic material forced at a frequency $\omega$, with the $e^{i\omega t}$ convention, responds as $\sigma = \mu^{*}(\omega)\,\varepsilon$ (the correspondence principle), where the complex modulus is the reciprocal of the complex compliance, $\mu^{*} = 1/J^{*}$. The complex compliance follows from the creep function $J(t)$, the strain that follows a unit step in stress (Efroimsky 2012):

$$J^{*}(\omega) = J(0) + \int_{0}^{\infty}\dot{J}(t)\,e^{-i\omega t}\,dt.$$

With the unrelaxed compliance $J = 1/\mu$ and the Maxwell time $\tau_{M} = \eta/\mu$, the elements TidalPy combines are:

- Maxwell, a spring and a dashpot in series, with $J(t) = J + t/\eta$:

  $$J_\mathrm{maxwell} = J - \frac{i}{\eta\,\omega}.$$

- Voigt-Kelvin, a spring and a dashpot in parallel, with the compliance $J_{v} = J/f_{J}$ and the viscosity $\eta_{v} = f_{\eta}\,\eta$ set by `voigt_modulus_frac` ($f_{J}$) and `voigt_viscosity_frac` ($f_{\eta}$), and $J(t) = J_{v}\left(1 - e^{-t/(J_{v}\eta_{v})}\right)$:

  $$J_\mathrm{voigt} = \frac{J_{v}}{1 + i\,J_{v}\eta_{v}\,\omega} = \frac{J_{v}}{1 + (J_{v}\eta_{v}\omega)^{2}} - i\,\frac{J_{v}^{2}\eta_{v}\,\omega}{1 + (J_{v}\eta_{v}\omega)^{2}}.$$

- Andrade, a Maxwell element plus the transient creep $\beta t^{\alpha}$, with $\beta = J\,(\zeta\tau_{M})^{-\alpha}$, where $\zeta$ is the ratio of the Andrade timescale to the Maxwell time, and $J(t) = J + t/\eta + \beta t^{\alpha}$ (Efroimsky 2012; Renaud and Henning 2018):

  $$J_\mathrm{andrade} = J_\mathrm{maxwell} + J\left(\zeta\tau_{M}\,\omega\right)^{-\alpha}\Gamma(1+\alpha)\left[\cos\frac{\pi\alpha}{2} - i\sin\frac{\pi\alpha}{2}\right],$$

  where $\Gamma$ is the gamma function.

- Zener, the standard linear solid: a spring of the relaxed modulus $r\mu$ in parallel with a Maxwell element whose spring is the rest, $(1-r)\mu$, set by `relaxed_modulus_frac` ($r$). With $\tau = \eta/((1-r)\mu)$ (Nowick and Berry 1972):

  $$\mu^{*}_\mathrm{zener} = r\mu + (1-r)\mu\,\frac{i\omega\tau}{1 + i\omega\tau}.$$

  It is $\mu$ at high frequency and relaxes to $r\mu$, not to zero, at low frequency; its loss peaks at $\omega\tau = 1$ with $\mathrm{Im}[\mu^{*}] = (1-r)\mu/2$. $r = 0$ is Maxwell and $r = 1$ is elastic.

- Seismic Q, which takes its loss from a quality factor rather than a viscosity. It reads its viscosity input as $Q_\mathrm{ref}$, the quality factor measured at a reference frequency $\omega_\mathrm{ref}$ (`reference_frequency_rad_s`), and $\mu$ is the modulus measured there too. With $s = \omega_\mathrm{ref}/|\omega|$ and the exponent $a$ (`q_frequency_exponent`, in $[0, 1)$):

  $$Q(\omega) = Q_\mathrm{ref}\,s^{-a}, \qquad \mathrm{Re}[\mu^{*}_\mathrm{sq}] = \frac{\mu}{1 + D(s)/Q_\mathrm{ref}}, \qquad \mathrm{Im}[\mu^{*}_\mathrm{sq}] = \frac{\mathrm{Re}[\mu^{*}_\mathrm{sq}]}{Q(\omega)},$$

  $$D(s) = \cot\frac{\pi a}{2}\left(s^{a} - 1\right), \qquad D(s) \to \frac{2}{\pi}\ln s \quad (a = 0).$$

  $D$ is the dispersion that causality (Kramers-Kronig) requires of a loss with that frequency dependence: even a $Q$ that never changes softens the modulus logarithmically toward low frequency, by about 5% from 1 s to a semidiurnal tide at $Q = 143$. To first order in $1/Q$ this is the form of Kanamori and Anderson (1977) for $a = 0$ and of Wahr and Bergen (1986) for $a > 0$; written as a compliance, the storage modulus stays positive. At $\omega = \omega_\mathrm{ref}$ it returns $\mu(1 + i/Q_\mathrm{ref})$ exactly.

Elements in series add their compliances: Burgers is Maxwell plus Voigt-Kelvin, and Sundberg-Cooper is Andrade plus Voigt-Kelvin (Sundberg and Cooper 2010). The single-element limits are the elastic spring, $\mu^{*} = \mu$, and the Newtonian dashpot, $\mu^{*} = i\eta\omega$. The loss tangent is $\mathrm{Im}[\mu^{*}]/\mathrm{Re}[\mu^{*}] = 1/Q$; seismic Q specifies it directly.

## Inheritance

```
c_TidalPyBaseClass
  └── c_PhysicsBase
        └── c_RheologyBase  (abstract)
              ├── c_Elastic   alias "off"
              ├── c_Viscous   alias "newton"
              ├── c_Voigt     aliases "voigt-kelvin", "voigt_kelvin"
              ├── c_Maxwell
              ├── c_Burgers
              ├── c_Andrade
              ├── c_Sundberg  aliases "sundberg-cooper", "sundberg_cooper"
              ├── c_Zener     aliases "sls", "standard_linear_solid"
              └── c_SeismicQ  aliases "constant_q", "power_law_q"; name "seismic_q"
```

`c_RheologyBase` declares `calc_complex_modulus` pure virtual and supplies the vectorized loops, the configuration export, and the binary encoding that every model inherits. The Python classes (`Elastic`, `Viscous`, `Voigt`, `Maxwell`, `Burgers`, `Andrade`, `Sundberg`, `Zener`, `SeismicQ`) are thin Cython wrappers holding a pointer to the C++ object, and `RheologyBase` is their shared Python base.

## Models

Simple models (Elastic, Viscous, Maxwell, Voigt, Zener, SeismicQ) are evaluated in closed form. The composites (Burgers, Andrade, Sundberg) place elements in series, so their compliances add and the modulus is the reciprocal of the sum, $\mu^* = 1 / \sum_i J_i$. Those element compliances are internal intermediates and are not exposed.

| Model | Complex modulus $\mu^*$ [Pa] | Parameters | Character |
|---|---|---|---|
| `Elastic` (`off`) | $\mu$ | - | No dissipation at any frequency. |
| `Viscous` (`newton`) | $i \eta \omega$ | - | No stored energy; pure loss. |
| `Maxwell` | $1 / J_\mathrm{maxwell}$ | - | One relaxation peak at $\omega\tau = 1$; loss falls as $\omega^{-1}$ above it. |
| `Voigt` (`voigt-kelvin`) | $1 / J_\mathrm{voigt} = \mu f_J + i \omega \eta_v$ | `voigt_modulus_frac`, `voigt_viscosity_frac` | Stiffens without limit at high frequency; rarely used alone. |
| `Burgers` | $1 / (J_\mathrm{maxwell} + J_\mathrm{voigt})$ | `voigt_modulus_frac`, `voigt_viscosity_frac` | Maxwell plus a secondary peak from the Voigt arm. |
| `Andrade` | $1 / J_\mathrm{andrade}$ | `alpha`, `zeta` | Maxwell plus a transient term; loss falls only as $\omega^{-\alpha}$. |
| `Sundberg` (`sundberg-cooper`) | $1 / (J_\mathrm{andrade} + J_\mathrm{voigt})$ | `alpha`, `zeta`, `voigt_modulus_frac`, `voigt_viscosity_frac` | Andrade's high-frequency tail plus Burgers' secondary peak. |
| `Zener` (`sls`) | $\mu^{*}_\mathrm{zener}$ | `relaxed_modulus_frac` | One relaxation peak like Maxwell, but relaxes to $r\mu$ instead of zero. |
| `SeismicQ` (`seismic_q`, `constant_q`) | $\mu^{*}_\mathrm{sq}$ | `reference_frequency`, `q_frequency_exponent` | Its viscosity input is a quality factor. A set $Q$ at every frequency ($a = 0$), or one falling toward low frequency as $\omega^{a}$. |

The element compliances $J_\mathrm{maxwell}$, $J_\mathrm{voigt}$, and $J_\mathrm{andrade}$ and the Zener and seismic Q moduli are defined in the Physics section above.

| Parameter | Default | Meaning |
|---|---|---|
| `alpha` | 0.3 | Andrade exponent. Sets how slowly the transient loss decays with frequency; laboratory values for silicates cluster near 0.2 to 0.4. |
| `zeta` | 1.0 | Ratio of the Andrade timescale to the Maxwell time. Larger values push the transient term down. |
| `voigt_modulus_frac` | 5.0 | Stiffness of the Voigt arm's spring relative to the main spring. The arm's compliance is the material compliance divided by this value. |
| `voigt_viscosity_frac` | 0.02 | Viscosity of the Voigt arm's dashpot as a fraction of the material viscosity. |
| `relaxed_modulus_frac` | 0.5 | Zener relaxed modulus as a fraction of the unrelaxed one, in \[0, 1\]; a value outside raises `ValueError`. |
| `reference_frequency` | $2\pi$ | Seismic Q: the frequency \[rad s⁻¹\] at which its quality factor and modulus were measured (config key `reference_frequency_rad_s`). The default is a 1 s period, PREM's. Must be positive and finite. |
| `q_frequency_exponent` | 0.0 | Seismic Q: the exponent $a$ of $Q \propto \omega^{a}$, in \[0, 1). Values near 0.1 to 0.3 describe Earth's mantle between seismic and tidal periods. |

> [!WARNING]
> Rheologies and their parameters are a very active area of research. The properties can vary greatly for different material and even for the same material that has had different histories (previous cracking, is porous, is hydrated or desiccated, etc.). TidalPy's defaults are roughly those applicable to Earth's upper mantle, but the uncertainties are large. We highly encourage users to read up on the latest research for the material under investigation or treat these as free parameters rather than stick with TidalPy's defaults.

### Behavior at the Limits

At zero frequency Maxwell, Burgers, Andrade, and Sundberg return approx. zero. In this scenario there is unlimited time to flow, a viscoelastic body supports no static rigidity. Elastic returns $\mu$, Viscous returns zero, Voigt returns $\mu f_J$, and Zener returns $r\mu$. Seismic Q treats zero frequency as no forcing and returns $\mu$ with no loss, as it does for an infinite $Q$; a $Q$ that is not positive returns `NaN`.

At negative frequency Elastic, Viscous, Voigt, Maxwell, Burgers, Zener, and SeismicQ mirror the imaginary part, but Andrade and Sundberg return `NaN`: their transient term raises a negative quantity to a fractional power. Always pass the absolute value of the forcing frequency. TidalPy's own tidal solvers do this; a direct call does not.

### Choosing a Model

`Elastic` gives deformation without dissipation. It should be used to isolate the elastic part of a Love number or for a layer that is effectively rigid on the forcing timescale.

`Maxwell` is the traditional rheology used in tidal studies. It is a good choice when comparing against published Love numbers, since most of the literature uses it. Its weakness is the high-frequency tail which has dissipation fall off quickly, as $\omega^{-1}$, which underestimates the dissipation response of real silicates during fast forcing.

`Andrade` and `Sundberg` are the models to use when the forcing is fast compared with the Maxwell time, which is the usual situation for a cool, stiff, or rapidly forced body. Their loss falls only as $\omega^{-\alpha}$, and for tidal problems that difference can be orders of magnitude in the heating rate (Renaud & Henning 2018).

`Zener` suits a response that relaxes only partway. A Maxwell bulk rheology lets a layer's bulk modulus relax to zero at long periods, which is not realistic. A Zener bulk rheology relaxes it to $r K$. Melt-driven compaction is the usual case: a partially molten rock's bulk modulus relaxes from its unrelaxed (undrained) value toward its drained one as melt moves, and a material's compaction law can supply the bulk viscosity that sets the rate (see [Bulk Mixing](../PartialMelt/partial_melt_models.md#bulk-mixing)). Pick $r$ as the drained-to-unrelaxed ratio; for melt in isolated pockets it is near 0.9 at 10% melt, and melt films lower it.

`SeismicQ` takes the loss straight from a measured quality factor, with no viscosity and no model of the relaxation behind it. It is what a seismic profile such as PREM supports, and a world built from a radial data file gives it to every solid layer when the world sets `q_provided = true` (see the [TOML schema](../Structures/config/toml_schema.md)). What it assumes is how $Q$ changes between the reference frequency and the forcing frequency: $a = 0$ keeps the seismic $Q$, while laboratory and geodetic constraints put Earth's mantle nearer $a = 0.1$ to $0.3$, so a tidal $Q$ several times lower than the seismic one. That choice, not the rest of the model, sets the tidal dissipation. Because its viscosity input is a quality factor, pair it only with a material whose viscosity slot holds one; a layer whose material carries a viscosity model would have that viscosity read as $Q$.

`Burgers` and `Voigt` are mainly useful for reproducing published work that used them, or for deliberately placing a secondary relaxation peak at a chosen frequency. `Viscous` exists for completeness and for the fluid limit.

## Example Usage

Instantiate a class directly, or resolve one by name.

```python
from TidalPy.Rheology import Maxwell, Andrade, make_rheology

maxwell_model = Maxwell()
andrade_model = Andrade(alpha=0.25, zeta=2.0)

# Case-insensitive; every registered alias is accepted.
sundberg_model = make_rheology("Sundberg-Cooper", {"alpha": 0.4, "zeta": 2.0})
```

Constructors take a model's parameters positionally, in the order `get_parameter_info()` lists them, or as keywords by argument name or config key. `make_rheology(model_name, config=None)` recognizes every name and alias in the inheritance tree above and builds that model from `config`, keyed by config key (`alpha`, `zeta`, `voigt_modulus_frac`, `voigt_viscosity_frac`, `relaxed_modulus_frac`, `reference_frequency_rad_s`, `q_frequency_exponent`, as each model reads them). Absent keys take the model's default. An unrecognized model name, a key the model does not read, or a value outside a parameter's bounds (an Andrade `alpha` outside (0, 1], a Zener `relaxed_modulus_frac` outside [0, 1], a seismic Q exponent outside [0, 1)) raises `ValueError` naming the closest accepted name or key.

Model parameters are fixed at construction and read as attributes (`andrade_model.alpha`, `sundberg_model.voigt_viscosity_frac`), and `parameters` returns them all. `with_parameters(**changes)` returns a new model with some changed, leaving this one as it was.

### Factory Internals

Each model declares its parameters once, in a `parameter_specs()` table (`c_SpecModel`, `Utilities/classes/spec_model_.hpp`), which gives its construction, validation, config entries, and binary record. `c_rheology_registry()` lists each model's names (canonical first, then aliases), binary class id, and constructor; `c_find_rheology(name, params)` returns a `std::unique_ptr<c_RheologyBase>` from it and throws `std::invalid_argument` for an unknown name or parameter, and `c_rheology_from_binary` rebuilds a saved model. Every C++ consumer uses these, including layers attaching a rheology and the binary loader. The Python classes and `make_rheology` are thin wrappers over them.

### Vectorized Evaluation

The scalar call takes three floats and returns a Python `complex`.

```python
from TidalPy.Rheology import Maxwell

modulus   = 50.0e9    # Pa
viscosity = 1.0e20    # Pa s
frequency = 1.0e-5    # rad s-1

model = Maxwell()
complex_shear = model.calc_complex_modulus(modulus, viscosity, frequency)
```

`calc_complex_modulus` also takes arrays, broadcast together, and returns a `complex128` array of their shape. Three vectorized methods, defined once on the base class so every model inherits them, name the array patterns that occur in a tidal calculation and always return a one-dimensional array.

| Method | Vectorized over | Typical use |
|---|---|---|
| `calc_complex_modulus_vectorize_modulus(modulus[], viscosity[], frequency)` | Equal-length modulus and viscosity arrays at one frequency | A radial profile of a planet at one forcing frequency. |
| `calc_complex_modulus_vectorize_frequency(modulus, viscosity, frequency[])` | Frequency array at one modulus and viscosity | A frequency sweep of one material. |
| `calc_complex_modulus_vectorize_all(modulus[], viscosity[], frequency[])` | Three equal-length arrays, element by element | Paired inputs already zipped by the caller. |

```python
import numpy as np
from TidalPy.Rheology import Andrade

model = Andrade()
radial_moduli     = np.linspace(40.0e9, 80.0e9, 100)
radial_viscosity  = np.logspace(18.0, 22.0, 100)

profile = model.calc_complex_modulus_vectorize_modulus(radial_moduli, radial_viscosity, 1.0e-5)
sweep   = model.calc_complex_modulus_vectorize_frequency(50.0e9, 1.0e20, np.logspace(-7, -4, 50))
```

The three methods differ only in which arguments they take as arrays. Each accepts array-likes, returns a `complex128` NumPy array, and calls the one C++ method, `c_RheologyBase::calc_complex_modulus_vectorize(modulus, viscosity, frequency, out)`, which broadcasts a length-one input against the others and fills the caller-supplied `std::vector<std::complex<double>>`. Mismatched input lengths raise `ValueError`.

### Convenience Functions

Each model also has a lower-case free function that builds the C++ model for the one call, evaluates it, and returns, with no Python object left behind.

```python
import numpy as np
from TidalPy.Rheology import maxwell, andrade

# Scalars in, Python complex out.
complex_shear = maxwell(50.0e9, 1.0e20, 1.0e-5)

# Arrays in, complex128 ndarray out. Any of the three may be an array.
profile = maxwell(np.array([1.0e10, 5.0e10]), np.array([1.0e19, 1.0e20]), 1.0e-5)
sweep   = andrade(50.0e9, 1.0e20, np.logspace(-7, -4, 50), alpha=0.3, zeta=1.0)
```

The signatures follow the classes: `elastic/viscous/maxwell(modulus, viscosity, frequency)`, `voigt/burgers(modulus, viscosity, frequency, voigt_modulus_frac=5.0, voigt_viscosity_frac=0.02)`, `andrade(modulus, viscosity, frequency, alpha=0.3, zeta=1.0)`, `sundberg(modulus, viscosity, frequency, alpha=0.3, zeta=1.0, voigt_modulus_frac=5.0, voigt_viscosity_frac=0.02)`, `zener(modulus, viscosity, frequency, relaxed_modulus_frac=0.5)`, and `seismic_q(modulus, quality_factor, frequency, reference_frequency_rad_s=2π, q_frequency_exponent=0.0)`. The model parameters are always scalars; `modulus`, `viscosity`, and `frequency` may each be a float or an array and are broadcast together, with the most specific vectorized routine chosen for the pattern supplied.

### Attaching a Rheology to a `Layer`

```python
from TidalPy.Material import Material, Phase
from TidalPy.Rheology import Maxwell, make_rheology
from TidalPy.Structures.layers import Layer

# The static moduli and viscosities the rheology works on belong to the layer's material.
rock = Phase(
    eos={"model": "constant", "reference_density_kg_m3": 3300.0, "bulk_modulus_pa": 1.0e11},
    shear_modulus={"model": "constant", "shear_modulus_pa": 5.0e10},
    shear_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e20},
    bulk_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e20},
    shear_rheology="andrade")                                   # The phase's default rheology
mantle = Layer(
    "mantle",
    0,
    0.0,
    1.0e6,
    material=Material(solid=rock),
    temperature=1600.0)

mantle.shear_rheology = Maxwell()                               # Overrides the phase's default
mantle.bulk_rheology = make_rheology("andrade", {"alpha": 0.3})
complex_shear = mantle.calc_complex_shear_modulus(1.0e-5)       # [Pa] at 1e-5 rad s-1

mantle.shear_rheology = None                                    # Back to the phase's Andrade
```

A layer's rheology is its own override when it has one, else its material's default (the solid phase's, or the liquid's for a liquid-only material). With neither, `calc_complex_shear_modulus` returns the static modulus as a purely real complex number, which is elastic behavior. The layer and the phase share the model rather than copying it, so one model can serve several layers. The one-argument `calc_complex_shear_modulus(frequency)` evaluates the material at zero pressure and the layer's own temperature. After a world's EOS solve, `calc_complex_shear_modulus(radius, frequency)` uses the post-melt modulus and viscosity the solve found at that radius. The equivalent declarative forms are a `[layers.<name>.shear_rheology]` table in a world's TOML, keyed by `model` plus any parameters, for the override, and a `shear_rheology` table in a phase of the material for the default; see the [TOML schema](../Structures/config/toml_schema.md) and [Layer](../Structures/layers/layer.md).

## Serialization

Every model supports the standard TidalPy interfaces.

| Call | Result |
|---|---|
| `get_config_dict()` | A dict of `model` plus the model's own parameters, in the same form the world builder reads. |
| `save_config(path)` | That dict written as a TOML file. |
| `save_binary(path)` / `load_binary(path)` | The model's TidalPy binary record, its parameters written by key. A load reads a record of the same model into the wrapper; a record saved by another model raises `IOError`. |

## Adding a New Rheology

To add one named `Foo`:

**C++ (`TidalPy/Rheology/rheology_.hpp`)**

1. If the constitutive law is a new series combination, add an internal `detail::element_compliance_*` helper, and a free `rheo_modulus_foo(modulus, viscosity, frequency, ...)` function. A model with a closed form can compute its modulus inline instead.
2. Add `c_Foo : public c_SpecModel<c_Foo, c_RheologyBase>`: its `parameter_specs()` table (argument name, config key, member, default, bounds, one-line description per parameter), `C_CLASS_ID`, two constructors that call `p_initialize`, and `calc_complex_modulus(modulus, viscosity, frequency)`. Override `p_validate` for checks across parameters and `p_update_derived` for values cached from them (the Andrade factors, say).
3. Reserve a unique `BinaryClassID::Foo` in `Utilities/binary/binary_.hpp`; rheology models occupy the 30X range.
4. Add one row to `c_rheology_registry()` with the model's names and aliases.

**Cython (`rheology.pyx`)**

5. Add `cdef class Foo(RheologyBase)` with a docstring and `MODEL_NAME = "foo"`, include it in the `ModelFamily` list, and add the lower-case `foo(modulus, viscosity, frequency, ...)` convenience function.

**Package, tests, and documentation**

6. Export `Foo` and `foo` from `TidalPy/Rheology/__init__.py`.
7. Add `Foo` to the physics tests in `Tests/Test_Rheology/test_rheology_01.py` (the modulus against an independent reference). The generic tests in `Tests/Test_Utilities/Test_Classes/test_spec_models_01.py` cover its parameters, config, binary record, and errors without changes.
8. Document the model here with its formula, parameters, and references.

## References

- Henning, W. G., O'Connell, R. J., and Sasselov, D. D. (2009). Tidally heated terrestrial exoplanets: Viscoelastic response models. *The Astrophysical Journal*, 707(2), 1000-1015. [DOI](https://doi.org/10.1088/0004-637X/707/2/1000). Maxwell, Voigt-Kelvin, and Burgers.
- Efroimsky, M. (2012). Tidal dissipation compared to seismic dissipation: In small bodies, Earths, and super-Earths. *The Astrophysical Journal*, 746(2), 150. [DOI](https://doi.org/10.1088/0004-637X/746/2/150). Complex compliances and Love numbers.
- Renaud, J. P., and Henning, W. G. (2018). Increased tidal dissipation using advanced rheological models: Implications for Io and tidally active exoplanets. *The Astrophysical Journal*, 857(2), 98. [DOI](https://doi.org/10.3847/1538-4357/aab784). Andrade and Sundberg-Cooper.
- Nowick, A. S., and Berry, B. S. (1972). *Anelastic Relaxation in Crystalline Solids*. Academic Press. The standard linear solid (Zener model).
- Kanamori, H., and Anderson, D. L. (1977). Importance of physical dispersion in surface wave and free oscillation problems: Review. *Reviews of Geophysics*, 15(1), 105-112. [DOI](https://doi.org/10.1029/RG015i001p00105). The dispersion of a constant $Q$.
- Wahr, J., and Bergen, Z. (1986). The effects of mantle anelasticity on nutations, earth tides, and tidal variations in rotation rate. *Geophysical Journal of the Royal Astronomical Society*, 87(2), 633-668. [DOI](https://doi.org/10.1111/j.1365-246X.1986.tb06642.x). A power-law $Q$ and its dispersion carried from seismic to tidal frequencies.
- Sundberg, M., and Cooper, R. F. (2010). A composite viscoelastic model for incorporating grain boundary sliding and transient diffusion creep; correlating creep and attenuation responses for materials with a fine grain size. *Philosophical Magazine*, 90. The Sundberg-Cooper composite.
