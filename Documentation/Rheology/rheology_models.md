# Rheology Models (`Rheology`)

_Updated: 2026-10-07_

A rheology model maps a material's static (purely real) modulus and viscosity onto a complex modulus $\mu^*(\omega)$ \[Pa\] at a forcing frequency. The real part is the storage modulus, the stress in phase with the strain; the imaginary part is the loss, which converts mechanical work into frictional heat. Their ratio $\mathrm{Im}[\mu^*]/\mathrm{Re}[\mu^*]$ is the loss tangent, the inverse of the quality factor $Q$ (Efroimsky 2012).

The models do not know whether they compute a shear or a bulk response: give a shear modulus with a shear viscosity, or a bulk modulus with a bulk viscosity. The shear response dominates tidal dissipation in solid bodies, and the bulk response is usually left elastic, though this is an active area of research.

```python
from TidalPy.Rheology import Maxwell, Andrade, make_rheology

maxwell_model = Maxwell()
andrade_model = Andrade(alpha=0.25, zeta=2.0)

# Case-insensitive; every registered alias is accepted.
sundberg_model = make_rheology("Sundberg-Cooper", {"alpha": 0.4, "zeta": 2.0})

modulus   = 50.0e9    # Pa
viscosity = 1.0e20    # Pa s
frequency = 1.0e-5    # rad s-1
complex_shear = maxwell_model.calc_complex_modulus(modulus, viscosity, frequency)
```

Every model has one method, `calc_complex_modulus(modulus, viscosity, frequency)`: the unrelaxed modulus $\mu$ \[Pa\], the viscosity $\eta$ \[Pa s\], and the forcing frequency $\omega$ \[rad s$^{-1}$\] in, $\mu^*$ \[Pa\] out. The argument order is modulus first throughout the module. Three floats give a Python `complex`.

## Models

Elastic, Viscous, Maxwell, Voigt, Zener, and SeismicQ are closed forms. Burgers, Andrade, and Sundberg place elements in series, so $\mu^* = 1 / \sum_i J_i$, with the element compliances of [Physics](#physics).

| Class (name, aliases) | Complex modulus $\mu^*$ [Pa] | Parameters | Character and use |
|---|---|---|---|
| `Elastic` (`elastic`, `off`) | $\mu$ | - | No dissipation at any frequency. The elastic part of a Love number, or a layer effectively rigid on the forcing timescale. |
| `Viscous` (`viscous`, `newton`) | $i \eta \omega$ | - | No stored energy; pure loss. The fluid limit. |
| `Maxwell` (`maxwell`) | $1 / J_\mathrm{maxwell}$ | - | One relaxation peak at $\omega\tau = 1$; loss falls as $\omega^{-1}$ above it, which underestimates real silicates under fast forcing. The traditional tidal rheology, for comparison with published Love numbers. |
| `Voigt` (`voigt`, `voigt-kelvin`, `voigt_kelvin`) | $1 / J_\mathrm{voigt} = \mu f_J + i \omega \eta_v$ | `voigt_modulus_frac`, `voigt_viscosity_frac` | Stiffens without limit at high frequency; rarely used alone. |
| `Burgers` (`burgers`) | $1 / (J_\mathrm{maxwell} + J_\mathrm{voigt})$ | `voigt_modulus_frac`, `voigt_viscosity_frac` | Maxwell plus a secondary peak from the Voigt arm. For reproducing published work, or a secondary peak at a chosen frequency. |
| `Andrade` (`andrade`) | $1 / J_\mathrm{andrade}$ | `alpha`, `zeta` | Maxwell plus a transient term; loss falls only as $\omega^{-\alpha}$. For forcing fast compared with the Maxwell time (a cool, stiff, or rapidly forced body), where it can raise the heating by orders of magnitude (Renaud and Henning 2018). |
| `Sundberg` (`sundberg`, `sundberg-cooper`, `sundberg_cooper`) | $1 / (J_\mathrm{andrade} + J_\mathrm{voigt})$ | `alpha`, `zeta`, `voigt_modulus_frac`, `voigt_viscosity_frac` | Andrade's high-frequency tail plus Burgers' secondary peak; used as Andrade is. |
| `Zener` (`zener`, `sls`, `standard_linear_solid`) | $\mu^{*}_\mathrm{zener}$ | `relaxed_modulus_frac` | One relaxation peak like Maxwell, but relaxes to $r\mu$ instead of zero. |
| `SeismicQ` (`seismic_q`, `constant_q`, `power_law_q`) | $\mu^{*}_\mathrm{sq}$ | `reference_frequency`, `q_frequency_exponent` | Its viscosity input is a quality factor: a set $Q$ at every frequency ($a = 0$), or one falling toward low frequency as $\omega^{a}$. |

| Parameter | Default | Meaning |
|---|---|---|
| `alpha` | 0.3 | Andrade exponent: how slowly the transient loss decays with frequency. Laboratory values for silicates cluster near 0.2 to 0.4. |
| `zeta` | 1.0 | Ratio of the Andrade timescale to the Maxwell time. Larger values push the transient term down. |
| `voigt_modulus_frac` | 5.0 | Stiffness of the Voigt arm's spring relative to the main spring; the arm's compliance is the material's divided by this value. |
| `voigt_viscosity_frac` | 0.02 | Viscosity of the Voigt arm's dashpot as a fraction of the material viscosity. |
| `relaxed_modulus_frac` | 0.5 | Zener relaxed modulus as a fraction of the unrelaxed one, in \[0, 1\]. |
| `reference_frequency` | $2\pi$ | Seismic Q: the frequency \[rad s⁻¹\] at which its quality factor and modulus were measured (config key `reference_frequency_rad_s`). The default is PREM's 1 s period. Must be positive and finite. |
| `q_frequency_exponent` | 0.0 | Seismic Q: the exponent $a$ of $Q \propto \omega^{a}$, in \[0, 1). Values near 0.1 to 0.3 describe Earth's mantle between seismic and tidal periods. |

### Zener and Seismic Q

A Maxwell bulk rheology relaxes a layer's bulk modulus to zero at long periods, which is not realistic; a Zener bulk rheology relaxes it to $r K$. The usual case is melt-driven compaction: a partially molten rock's bulk modulus relaxes from its unrelaxed (undrained) value toward its drained one as melt moves, at a rate set by a bulk viscosity that a material's compaction law can supply (see [Bulk Mixing](../PartialMelt/partial_melt_models.md#bulk-mixing)). Set $r$ to the drained-to-unrelaxed ratio: near 0.9 at 10% melt for melt in isolated pockets, lower with melt films.

`SeismicQ` takes the loss from a measured quality factor, with no model of the relaxation behind it. It is what a seismic profile such as PREM supports, and a world built from a radial data file gives it to every solid layer when the world sets `q_provided = true` (see the [TOML schema](../Structures/config/toml_schema.md)). Its one assumption, which sets the tidal dissipation, is how $Q$ changes between the reference and the forcing frequency: $a = 0$ keeps the seismic $Q$, while laboratory and geodetic constraints put Earth's mantle nearer $a = 0.1$ to $0.3$, a tidal $Q$ several times lower. Pair it only with a material whose viscosity slot holds a quality factor; a viscosity model's output would be read as $Q$.

> [!WARNING]
> Rheologies and their parameters are an active area of research. They vary greatly between materials, and for one material with its history (cracking, porosity, hydration). TidalPy's defaults roughly fit Earth's upper mantle, with large uncertainties. Read the latest research for the material under study, or treat these as free parameters.

## Common Tasks

Constructors take a model's parameters positionally, in the order `get_parameter_info()` lists them, or as keywords by argument name or config key. `make_rheology(model_name, config=None)` accepts every name and alias in the [Models](#models) table and builds the model from `config`, keyed by config key; absent keys take the defaults. An unknown name, a key the model does not read, or a value outside a parameter's bounds (an Andrade `alpha` outside (0, 1], a Zener `relaxed_modulus_frac` outside [0, 1], a seismic Q exponent outside [0, 1)) raises `ValueError` naming the closest accepted name or key.

Parameters are fixed at construction and read as attributes (`andrade_model.alpha`); `parameters`, `with_parameters(**changes)`, and the other members shared by every physics model are described with [`PhysicsBase`](../Utilities/classes.md#physicsbase).

### Vectorized Evaluation

`calc_complex_modulus` also takes arrays, broadcast together, and returns a `complex128` array of their shape. Three methods on the base class name the array patterns of a tidal calculation and always return a one-dimensional `complex128` array; a length-one input is broadcast against the others, and mismatched lengths raise `ValueError`.

| Method | Vectorized over | Typical use |
|---|---|---|
| `calc_complex_modulus_vectorize_modulus(modulus[], viscosity[], frequency)` | Equal-length modulus and viscosity arrays at one frequency | A radial profile at one forcing frequency. |
| `calc_complex_modulus_vectorize_frequency(modulus, viscosity, frequency[])` | Frequency array at one modulus and viscosity | A frequency sweep of one material. |
| `calc_complex_modulus_vectorize_all(modulus[], viscosity[], frequency[])` | Three equal-length arrays, element by element | Paired inputs. |

```python
import numpy as np
from TidalPy.Rheology import Andrade

model = Andrade()
radial_moduli     = np.linspace(40.0e9, 80.0e9, 100)
radial_viscosity  = np.logspace(18.0, 22.0, 100)

profile = model.calc_complex_modulus_vectorize_modulus(radial_moduli, radial_viscosity, 1.0e-5)
sweep   = model.calc_complex_modulus_vectorize_frequency(50.0e9, 1.0e20, np.logspace(-7, -4, 50))
```

### Convenience Functions

Each model has a lower-case free function that evaluates it once, with no Python object left behind.

```python
import numpy as np
from TidalPy.Rheology import maxwell, andrade

# Scalars in, Python complex out.
complex_shear = maxwell(50.0e9, 1.0e20, 1.0e-5)

# Arrays in, complex128 ndarray out. Any of the three may be an array.
profile = maxwell(np.array([1.0e10, 5.0e10]), np.array([1.0e19, 1.0e20]), 1.0e-5)
sweep   = andrade(50.0e9, 1.0e20, np.logspace(-7, -4, 50), alpha=0.3, zeta=1.0)
```

The signatures are `elastic/viscous/maxwell(modulus, viscosity, frequency)`, `voigt/burgers(modulus, viscosity, frequency, voigt_modulus_frac=5.0, voigt_viscosity_frac=0.02)`, `andrade(modulus, viscosity, frequency, alpha=0.3, zeta=1.0)`, `sundberg(modulus, viscosity, frequency, alpha=0.3, zeta=1.0, voigt_modulus_frac=5.0, voigt_viscosity_frac=0.02)`, `zener(modulus, viscosity, frequency, relaxed_modulus_frac=0.5)`, and `seismic_q(modulus, quality_factor, frequency, reference_frequency_rad_s=2π, q_frequency_exponent=0.0)`. Model parameters are scalars; `modulus`, `viscosity`, and `frequency` may each be a float or an array and are broadcast together.

### Attaching a Rheology to a `Layer`

A phase can carry default rheologies (`Phase(shear_rheology="andrade")`, or a `shear_rheology` table in a phase of a TOML material), and a layer's `shear_rheology` and `bulk_rheology` (in TOML, a `[layers.<name>.shear_rheology]` table) override them with a model, a model name, or a config table; setting `None` restores the default. A layer reads the solid phase's default, or the liquid's for a liquid-only material. With no rheology in effect, `calc_complex_shear_modulus` returns the static modulus as a purely real complex number (elastic). Models are shared, not copied, so one can serve several layers.

```python
from TidalPy.Rheology import Maxwell, make_rheology

# mantle: a Layer whose material's solid phase carries a default rheology
mantle.shear_rheology = Maxwell()                               # Overrides the phase's default
mantle.bulk_rheology = make_rheology("andrade", {"alpha": 0.3})
complex_shear = mantle.calc_complex_shear_modulus(1.0e-5)       # [Pa] at 1e-5 rad s-1
mantle.shear_rheology = None                                    # Back to the phase's default
```

`calc_complex_shear_modulus(frequency)` evaluates the material at zero pressure and the layer's temperature; after a world's EOS solve, `calc_complex_shear_modulus(radius, frequency)` uses the post-melt modulus and viscosity at that radius. See [Rheology Overrides](../Structures/layers/layer.md#rheology-overrides) and the [TOML schema](../Structures/config/toml_schema.md).

## Physics

The scale common to every model is the Maxwell time $\tau = \eta / \mu$, the time a material takes to relax a stress by viscous flow. Forcing much faster than $\tau$ finds the material effectively elastic, and much slower effectively fluid. Dissipation peaks when the forcing timescale is near the Maxwell time; the models differ mostly in how broad that window is and in what happens on its high-frequency side.

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

  $D$ is the dispersion that causality (Kramers-Kronig) requires of a loss with that frequency dependence: even a constant $Q$ softens the modulus logarithmically toward low frequency, by about 5% from 1 s to a semidiurnal tide at $Q = 143$. To first order in $1/Q$ this is the form of Kanamori and Anderson (1977) for $a = 0$ and of Wahr and Bergen (1986) for $a > 0$; written as a compliance, the storage modulus stays positive. At $\omega = \omega_\mathrm{ref}$ it returns $\mu(1 + i/Q_\mathrm{ref})$ exactly.

Elements in series add their compliances: Burgers is Maxwell plus Voigt-Kelvin, and Sundberg-Cooper is Andrade plus Voigt-Kelvin (Sundberg and Cooper 2010). The single-element limits are the elastic spring, $\mu^{*} = \mu$, and the Newtonian dashpot, $\mu^{*} = i\eta\omega$.

With the $e^{+i\omega t}$ convention a lagging response has a positive imaginary modulus, so the returned imaginary part is non-negative for a positive frequency. The models keep no memory between calls, so they are safe across threads and in any order.

## Behavior at the Limits

- Zero frequency: Maxwell, Burgers, Andrade, and Sundberg return approximately zero (with unlimited time to flow, a viscoelastic body has no static rigidity). Elastic returns $\mu$, Viscous zero, Voigt $\mu f_J$, and Zener $r\mu$. Seismic Q treats it as no forcing and returns $\mu$ with no loss.
- Seismic Q with an infinite $Q$ returns $\mu$ with no loss; a $Q$ that is not positive returns `NaN`.
- Negative frequency: Elastic, Viscous, Voigt, Maxwell, Burgers, Zener, and SeismicQ mirror the imaginary part, but Andrade and Sundberg return `NaN` (their transient term raises a negative number to a fractional power). Pass the absolute frequency; TidalPy's tidal solvers do, a direct call does not.
- A NaN viscosity (a phase without a viscosity law) gives a NaN modulus for a viscoelastic model.

## Serialization

`get_config_dict()` returns `model` plus the model's parameters, the form the world builder reads, and `save_config(path)` writes it as TOML. `save_binary(path)` / `load_binary(path)` write and read the model's TidalPy binary record, parameters by key; loading a record saved by another model raises `IOError`.

## C++ API

The models are `c_Elastic`, `c_Viscous`, `c_Voigt`, `c_Maxwell`, `c_Burgers`, `c_Andrade`, `c_Sundberg`, `c_Zener`, and `c_SeismicQ`, on the abstract `c_RheologyBase` in `TidalPy/Rheology/rheology_.hpp`; the Python classes wrap them. `c_RheologyBase` declares `calc_complex_modulus(modulus, viscosity, frequency)` and provides `calc_complex_modulus_vectorize(modulus, viscosity, frequency, out)`, which broadcasts a length-one input and fills a caller-supplied `std::vector<std::complex<double>>`. `c_find_rheology(name, params)` builds a model by name or alias as a `std::unique_ptr<c_RheologyBase>` and throws `std::invalid_argument` for an unknown name or parameter; `c_rheology_from_binary` rebuilds a saved model.

## Adding a New Rheology

To add one named `Foo`:

1. In `rheology_.hpp`, for a new series combination add an internal `detail::element_compliance_*` helper and a free `rheo_modulus_foo(modulus, viscosity, frequency, ...)`; a closed form can be computed inline instead.
2. Add `c_Foo : public c_SpecModel<c_Foo, c_RheologyBase>` with its `parameter_specs()` table (argument name, config key, member, default, bounds, one-line description), `C_CLASS_ID`, two constructors that call `p_initialize`, and `calc_complex_modulus`. Override `p_validate` for cross-parameter checks and `p_update_derived` for cached values.
3. Reserve a `BinaryClassID::Foo` in `Utilities/binary/binary_.hpp` (rheologies use the 30X range), and add a row to `c_rheology_registry()`.
4. In `rheology.pyx`, add `cdef class Foo(RheologyBase)` with a docstring and `MODEL_NAME = "foo"`, include it in the `ModelFamily` list, and add the `foo(modulus, viscosity, frequency, ...)` convenience function.
5. Export `Foo` and `foo` from `TidalPy/Rheology/__init__.py`.
6. Add `Foo` to the physics tests in `Tests/Test_Rheology/test_rheology_01.py` (the modulus against an independent reference). The generic tests in `Tests/Test_Utilities/Test_Classes/test_spec_models_01.py` cover parameters, config, binary record, and errors.
7. Document the model here with its formula, parameters, and references.

## References

- Henning, W. G., O'Connell, R. J., and Sasselov, D. D. (2009). Tidally heated terrestrial exoplanets: Viscoelastic response models. *The Astrophysical Journal*, 707(2), 1000-1015. [DOI](https://doi.org/10.1088/0004-637X/707/2/1000). Maxwell, Voigt-Kelvin, and Burgers.
- Efroimsky, M. (2012). Tidal dissipation compared to seismic dissipation: In small bodies, Earths, and super-Earths. *The Astrophysical Journal*, 746(2), 150. [DOI](https://doi.org/10.1088/0004-637X/746/2/150). Complex compliances and Love numbers.
- Renaud, J. P., and Henning, W. G. (2018). Increased tidal dissipation using advanced rheological models: Implications for Io and tidally active exoplanets. *The Astrophysical Journal*, 857(2), 98. [DOI](https://doi.org/10.3847/1538-4357/aab784). Andrade and Sundberg-Cooper.
- Nowick, A. S., and Berry, B. S. (1972). *Anelastic Relaxation in Crystalline Solids*. Academic Press. The standard linear solid (Zener model).
- Kanamori, H., and Anderson, D. L. (1977). Importance of physical dispersion in surface wave and free oscillation problems: Review. *Reviews of Geophysics*, 15(1), 105-112. [DOI](https://doi.org/10.1029/RG015i001p00105). The dispersion of a constant $Q$.
- Wahr, J., and Bergen, Z. (1986). The effects of mantle anelasticity on nutations, earth tides, and tidal variations in rotation rate. *Geophysical Journal of the Royal Astronomical Society*, 87(2), 633-668. [DOI](https://doi.org/10.1111/j.1365-246X.1986.tb06642.x). A power-law $Q$ and its dispersion carried from seismic to tidal frequencies.
- Sundberg, M., and Cooper, R. F. (2010). A composite viscoelastic model for incorporating grain boundary sliding and transient diffusion creep; correlating creep and attenuation responses for materials with a fine grain size. *Philosophical Magazine*, 90. The Sundberg-Cooper composite.
