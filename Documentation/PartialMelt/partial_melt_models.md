# Partial-Melt Models (`PartialMelt`)

_Updated: 2026-09-29_

A partial-melt model maps a material's pre-melt (solid) viscosity and shear modulus, together with its temperature, onto the post-melt viscosity and shear modulus, and reports the volumetric melt fraction it used to get there. Every model can also mix the melt into the density, weaken the bulk modulus with its own, much weaker law, and let melt set a compaction bulk viscosity; each is off by default (see Density and Bulk Response).

These quantities depend only on the temperature and pressure state fixed by the equation-of-state solve, not on the forcing frequency, so they are computed once per solve and cached. Only the downstream [rheology](../Rheology/index.md) step, which produces the complex modulus (frequency-dependent), is recomputed for each tidal mode.

## Melt Fraction

The volumetric melt fraction is model-independent.

$$\phi = \mathrm{clip}\left( \frac{T - T_\mathrm{solidus}}{T_\mathrm{liquidus} - T_\mathrm{solidus}},\; 0,\; 1 \right)$$

Below the solidus $\phi = 0$, above the liquidus $\phi = 1$, and a degenerate envelope with the solidus at or above the liquidus returns $\phi = 0$, which is the fully solid answer. A non-finite temperature has no melt state: the melt fraction is NaN, and so are the Spohn and Henning strengths (Off passes its pre-melt strengths through).

## Models

| Model | Aliases | Behavior |
|---|---|---|
| `OffPartialMelt` | `off`, `none` | No weakening. Post-melt strength equals pre-melt strength, but the melt fraction is still reported. |
| `SpohnPartialMelt` | `spohn`, `fischer`, `fischer_spohn` | Fischer and Spohn (1990). Above the solidus, post-melt viscosity and shear modulus depend on temperature alone, not on the pre-melt values. |
| `HenningPartialMelt` | `henning` | Henning et al. (2009) and Renaud and Henning (2018). Three regimes separated by a critical melt fraction. |

### Off

The melt fraction is computed and returned, and the strengths pass through untouched. Use it to turn off melt weakening, or for a layer that stays below its solidus.

### Spohn (Fischer and Spohn 1990)

$$\eta_\mathrm{post} = 10^{\,L_\eta + s_\eta (1/T - 1/T_\mathrm{solidus})}, \qquad \mu_\mathrm{post} = 10^{\,L_\mu + s_\mu (1/T - 1/T_\mathrm{solidus})}$$

with the slopes $s$ given by `fs_visc_power_slope` and `fs_shear_power_slope`, and the base-10 logarithms of the strengths at the solidus, $L$, by `fs_visc_log10_at_solidus` and `fs_shear_log10_at_solidus`. Both results are floored at the liquid limits, `liquid_viscosity` and `liquid_shear`.

Fischer and Spohn (1990) fit silicates with absolute laws, $10^{27000/T - 1}$ Pa s and $10^{82000/T - 40.6}$ Pa. Those are the form above at a 1600 K solidus ($L_\eta$ = 15.875, $L_\mu$ = 10.65), so the defaults reproduce them there. Anchoring the law at the model's own solidus keeps it usable for other materials: the absolute fit at an icy solidus of 273 K would give a shear modulus of $10^{260}$ Pa.

The law applies above the solidus ($\phi > 0$). At or below it the pre-melt viscosity and shear modulus are returned unchanged; the law itself grows without bound as the temperature falls (10$^{41}$ Pa at 1000 K with the defaults). Above the solidus the pre-melt values do not appear on the right-hand side, so the layer's viscosity model has no effect wherever the material is partially molten.

### Henning (2009, 2018)

Three regimes in the melt fraction, with a transition band running from `crit_melt_frac` to `crit_melt_frac + crit_melt_frac_width`. Write $\phi_c$ for the critical melt fraction and $T_\mathrm{break} = T_\mathrm{solidus} + \phi_c (T_\mathrm{liquidus} - T_\mathrm{solidus})$ for the temperature at which it is reached.

| Regime | Post-melt viscosity | Post-melt shear modulus |
|---|---|---|
| $\phi \le 0$ | $\eta_\mathrm{pre}$ | $\mu_\mathrm{pre}$ |
| $0 < \phi < \phi_c$ | $\eta_\mathrm{pre} \exp(-a_\eta \phi)$ | $\mu_\mathrm{pre} \exp[b_1 (1/T - 1/T_\mathrm{solidus})]$ |
| $\phi_c \le \phi \le \phi_c + w$ | $\eta_\mathrm{pre} \exp(-a_\eta \phi_c) \exp(-f_\eta (\phi - \phi_c))$ | $\mu_\mathrm{pre} \exp[b_1 (1/T_\mathrm{break} - 1/T_\mathrm{solidus})] \exp(-f_\mu (\phi - \phi_c))$ |
| $\phi > \phi_c + w$ | $\eta_\mathrm{liquid}$ | `liquid_shear` |

Here $a_\eta$ is `hn_visc_slope_1`, $b_1$ is `hn_shear_param_1`, and $f_\eta$ and $f_\mu$ are `hn_visc_falloff_slope` and `hn_shear_falloff_slope`. Every branch is floored at the liquid limits.

The shear law is Henning et al. (2009) Eq. 20, $\exp(40000/T - 25)$, anchored at the solidus: their constant 25 is $40000/1600$, the silicate solidus they calibrated to, so the default 1600 K solidus reproduces it exactly. Written this way the shear modulus is continuous at any solidus, where the published constant would stiffen an icy layer by a factor of $e^{100}$ as it began to melt. Configurations from a 0.8.0 pre-release that still set `hn_shear_param_2` get a warning and the key is ignored.

Below the critical melt fraction, melt sits in isolated pockets and weakens the solid framework gradually. Above it the framework loses contact and the material behaves as a crystal-laden liquid, a drop of many orders of magnitude. The breakdown band is a steep but finite bridge between the two regimes; its width is a numerical convenience that keeps the transition differentiable, not a measured quantity.

## Parameters

Each parameter carries two names: the constructor keyword, which is also the read-only property, and the config key used in a TOML table, a `make_partial_melt` config dictionary, and `get_config_dict()`. A dimensional config key ends in its unit; the code name does not.

| Parameter | Config key | Default | Units | Used by |
|---|---|---|---|---|
| `solidus` | `solidus_k` | 1600.0 | K | All |
| `liquidus` | `liquidus_k` | 2000.0 | K | All |
| `liquid_shear` | `liquid_shear_pa` | 1.0e-5 | Pa | All |
| `liquid_viscosity` | `liquid_viscosity_pas` | 0.2 | Pa s | All |
| `fs_visc_power_slope` | `fs_visc_power_slope_k` | 27000.0 | K | Spohn |
| `fs_visc_log10_at_solidus` | `fs_visc_log10_at_solidus` | 15.875 | log10 Pa s | Spohn |
| `fs_shear_power_slope` | `fs_shear_power_slope_k` | 82000.0 | K | Spohn |
| `fs_shear_log10_at_solidus` | `fs_shear_log10_at_solidus` | 10.65 | log10 Pa | Spohn |
| `crit_melt_frac` | `crit_melt_frac` | 0.5 | - | Henning |
| `crit_melt_frac_width` | `crit_melt_frac_width` | 0.05 | - | Henning |
| `hn_visc_slope_1` | `hn_visc_slope_1` | 13.5 | - | Henning |
| `hn_visc_falloff_slope` | `hn_visc_falloff_slope` | 370.0 | - | Henning |
| `hn_shear_param_1` | `hn_shear_param_1_k` | 40000.0 | K | Henning |
| `hn_shear_falloff_slope` | `hn_shear_falloff_slope` | 700.0 | - | Henning |

`liquid_shear` and `liquid_viscosity` are the shear modulus and viscosity assigned to material treated as pure liquid, and they are also the floors every model applies to its post-melt shear modulus and viscosity. The shear default is small but not exactly zero. The viscosity default is a molten silicate's; the packaged configuration gives each material type its own (rock 0.2, iron 1.3e-2, ice and high-pressure ice 8.9e-4 Pa s).

## Density and Bulk Response

The shear laws above are empirical fits in temperature. Melt's effect on the density and on the bulk response is of another kind: it follows from mixing two phases with their own equations of state. The models describe the melt phase once and use it three ways, each behind its own switch, all off by default:

| Switch | Effect |
|---|---|
| `density_melt_mixing` | The melt enters the density, so the structure solve sees it. |
| `bulk_melt_weakening` | The unrelaxed bulk modulus is the two-phase bound. |
| `bulk_viscosity_melt_weakening` | Melt sets a compaction bulk viscosity, which lets a bulk rheology dissipate. |

### Order of Evaluation

At every point of a world's equation-of-state solve the material evaluates, in this order:

1. The density law gives the solid density and bulk modulus at the local pressure (and temperature, for a thermal law).
2. The melt fraction follows from the temperature, and with `density_melt_mixing` the density becomes the two-phase mixture. The structure iteration and the dense readout call the same function, so the mass the solve integrated and the density `get_density` reports agree.
3. The shear modulus and viscosities come from their laws and models; the melt model then weakens the shear pair.
4. With `bulk_melt_weakening` the bulk modulus becomes the two-phase bound, using the post-melt shear modulus as the framework's. With `bulk_viscosity_melt_weakening` the bulk viscosity takes the compaction term, using the post-melt shear viscosity.

So the melt fraction is known before the density, and the bulk quantities are weakened after the shear ones they depend on.

### The Melt Phase

The melt follows a Murnaghan (1944) law, with its zero-pressure density $\rho_{l0}$ (`liquid_density`), bulk modulus $K_{l0}$ (`liquid_bulk_modulus`), and pressure derivative $K_l'$ (`liquid_bulk_modulus_derivative`):

$$K_l(P) = K_{l0} + K_l' P, \qquad \rho_l(P) = \rho_{l0} \left(1 + \frac{K_l' P}{K_{l0}}\right)^{1/K_l'}.$$

$K_l' = 0$ gives $\rho_l = \rho_{l0} e^{P/K_{l0}}$. Under tension, which only the structure solve's trial central pressures reach, the law continues at constant $K_{l0}$ so it stays finite. The law is isothermal. Because $\rho_l / (d\rho_l/dP) = K_l$, a layer that is wholly melt, with both switches on, is neutrally stratified under the bulk modulus the tidal equations see (see [dynamic liquids](../RadialSolver/dense_radial_solution.md#dynamic-liquid-layers-at-long-forcing-periods)).

### Density

With `density_melt_mixing = true` the two phases mix by volume at the same pressure:

$$\rho = (1 - \phi)\,\rho_s + \phi\,\rho_l(P),$$

where $\rho_s$ is the density law's value. Silicate melt is lighter than its source at low pressure, so melt lowers a rocky layer's density; water is denser than ice I, so it raises an ice shell's. A world whose mass was fitted without mixing will hold a different mass with it, and a partially molten layer couples its density to its temperature, so a thermal solve takes more passes to settle. The phase change's own compressibility (the melt fraction changing with pressure) and latent heat are not included.

### Bulk Modulus

Melt lowers the bulk modulus far less than the shear modulus: a melt keeps a bulk modulus of the same order as the solid's (about 20 GPa for a silicate melt at low pressure against about 130 GPa for the rock), while its shear modulus vanishes (Mavko 1980; Takei 2002). With `bulk_melt_weakening = true` the post-melt bulk modulus is the Hashin and Shtrikman (1963) bound for melt of bulk modulus $K_l(P)$ in a solid framework of bulk modulus $K_s$ (the pre-melt value), evaluated with the framework's post-melt shear modulus $\mu$:

$$K = K_s + \frac{\phi}{\dfrac{1}{K_l - K_s} + \dfrac{1 - \phi}{K_s + \tfrac{4}{3}\mu}}$$

While the framework holds, $\mu$ is close to the solid's and this is the upper bound for isolated melt pockets, a weak reduction (about 16 percent at $\phi = 0.1$ for $K_s$ = 130, $\mu$ = 60, $K_l$ = 20 GPa). Once the melt model has collapsed the framework's shear modulus (Henning past the critical melt fraction), the same expression becomes the Reuss (Wood 1955) average of a crystal suspension, and it reaches $K_l(P)$ at $\phi = 1$. Since $K_l$ rises with pressure faster than a rock's, the bound weakens less at depth, and not at all where the melt is as stiff as the solid.

This is the unrelaxed (undrained) modulus: the melt is held in place over a forcing cycle. Isolated pockets are the stiffest geometry, and melt that wets grain edges or forms films weakens the framework more (Mavko 1980; Takei 2002), so the bound is an upper limit.

### Bulk Viscosity

A melt-free rock has no viscous compaction, so its bulk response is elastic; melt adds one as it moves through the matrix. With `bulk_viscosity_melt_weakening = true` melt contributes a compaction bulk viscosity in series with the pre-melt one,

$$\frac{1}{\zeta} = \frac{1}{\zeta_\mathrm{pre}} + \frac{\phi^n}{c\,\eta},$$

with $\eta$ the post-melt shear viscosity, $c$ `melt_bulk_viscosity_coefficient`, and $n$ `melt_bulk_viscosity_exponent`. $n = 1$ with $c$ of order one is the classic compaction viscosity $\eta/\phi$ (McKenzie 1984); micromechanical models give a bulk viscosity of the same order as the shear viscosity instead (Takei and Holtzman 2009), which is $n = 0$. The series form is continuous at the solidus for $n > 0$. A non-finite pre-melt bulk viscosity counts as no pre-melt dashpot, so melt alone sets $\zeta$.

The bulk viscosity reaches the tides only through the layer's bulk rheology, which is `elastic` by default. A Maxwell bulk rheology would let the bulk modulus relax to zero at long periods, which no rock does. The [Zener](../Rheology/rheology_models.md) (standard linear solid) rheology relaxes it to a set fraction $r$ of the unrelaxed modulus instead, the drained-to-undrained ratio; for isolated pockets at 10% melt that ratio is near 0.9, and melt films lower it. Bulk dissipation in a partially molten layer can rival the shear dissipation (Kervazo et al. 2021).

| Parameter | Config key | Default | Units | Used by |
|---|---|---|---|---|
| `liquid_density` | `liquid_density_kg_m3` | 2750.0 | kg/m$^3$ | All, when `density_melt_mixing` is true |
| `liquid_bulk_modulus` | `liquid_bulk_modulus_pa` | 2.0e10 | Pa | All, when `density_melt_mixing` or `bulk_melt_weakening` is true |
| `liquid_bulk_modulus_derivative` | `liquid_bulk_modulus_derivative` | 5.0 | - | As `liquid_bulk_modulus` |
| `density_melt_mixing` | `density_melt_mixing` | false | - | All |
| `bulk_melt_weakening` | `bulk_melt_weakening` | false | - | All |
| `bulk_viscosity_melt_weakening` | `bulk_viscosity_melt_weakening` | false | - | All |
| `melt_bulk_viscosity_coefficient` | `melt_bulk_viscosity_coefficient` | 1.0 | - | All, when `bulk_viscosity_melt_weakening` is true |
| `melt_bulk_viscosity_exponent` | `melt_bulk_viscosity_exponent` | 1.0 | - | As `melt_bulk_viscosity_coefficient` |

The rock defaults are roughly an ultramafic silicate melt's. The packaged configuration gives liquid iron 7019 kg/m$^3$, 1.1e11 Pa, and $K'$ = 4.66, and liquid water 999.84 kg/m$^3$, 2.2e9 Pa, and $K'$ = 6.8. The melt and bulk dissipation demo (`Demos/Physics/19_melt_and_bulk_dissipation.ipynb`) works through each switch.

## Python API

```python
from TidalPy.PartialMelt import (
    HenningPartialMelt, OffPartialMelt, SpohnPartialMelt, make_partial_melt)

melt_model = HenningPartialMelt(solidus=1600.0, liquidus=2000.0, liquid_shear=1.0e-5)

phi = melt_model.calc_melt_fraction(1800.0)                 # 0.5
phi, post_viscosity, post_shear = melt_model.calc_partial_melt(
    temperature=1700.0,
    premelt_viscosity=1.0e22,   # Pa s
    premelt_shear=6.0e10,       # Pa
)

# Name factory: case-insensitive, aliases accepted.
spohn_model = make_partial_melt("fischer", {"solidus_k": 1500.0})
```

Constructors take the melt envelope plus their own parameters, all with the defaults from the table: 

`OffPartialMelt(solidus=1600.0, liquidus=2000.0, liquid_shear=1.0e-5, liquid_viscosity=0.2, bulk_melt_weakening=False, liquid_bulk_modulus=2.0e10, <melt phase>)`

`SpohnPartialMelt(solidus, liquidus, liquid_shear, fs_visc_power_slope=27000.0, fs_visc_log10_at_solidus=15.875, fs_shear_power_slope=82000.0, fs_shear_log10_at_solidus=10.65, liquid_viscosity, bulk_melt_weakening, liquid_bulk_modulus, <melt phase>)`

`HenningPartialMelt(solidus, liquidus, liquid_shear, crit_melt_frac=0.5, crit_melt_frac_width=0.05, hn_visc_slope_1=13.5, hn_visc_falloff_slope=370.0, hn_shear_param_1=40000.0, hn_shear_falloff_slope=700.0, liquid_viscosity, bulk_melt_weakening, liquid_bulk_modulus, <melt phase>)`

where `<melt phase>` is `liquid_bulk_modulus_derivative=5.0, liquid_density=2750.0, density_melt_mixing=False, bulk_viscosity_melt_weakening=False, melt_bulk_viscosity_coefficient=1.0, melt_bulk_viscosity_exponent=1.0`.

| Member | Returns | Description |
|---|---|---|
| `calc_melt_fraction(temperature)` | `float` | Melt fraction in [0, 1]. |
| `calc_partial_melt(temperature, premelt_viscosity, premelt_shear)` | `(phi, viscosity, shear_modulus)` | Melt fraction, post-melt viscosity [Pa s], post-melt shear modulus [Pa]. |
| `calc_liquid_density(pressure)`, `calc_liquid_bulk_modulus(pressure)` | `float` | The melt phase's density [kg/m$^3$] and bulk modulus [Pa]. |
| `calc_mixture_density(temperature, pressure, solid_density)` | `float` | Density [kg/m$^3$]; `solid_density` unless `density_melt_mixing` is on. |
| `calc_bulk_modulus_melt(temperature, pressure, premelt_bulk_modulus, framework_shear_modulus)` | `float` | Post-melt bulk modulus [Pa]; the pre-melt value unless `bulk_melt_weakening` is on. |
| `calc_bulk_viscosity_melt(temperature, premelt_bulk_viscosity, postmelt_shear_viscosity)` | `float` | Post-melt bulk viscosity [Pa s]; the pre-melt value unless `bulk_viscosity_melt_weakening` is on. |
| `solidus`, `liquidus`, `liquid_shear`, `liquid_viscosity`, `bulk_melt_weakening`, `liquid_bulk_modulus`, `liquid_bulk_modulus_derivative`, `liquid_density`, `density_melt_mixing`, `bulk_viscosity_melt_weakening`, `melt_bulk_viscosity_coefficient`, `melt_bulk_viscosity_exponent` | `float`, `bool` | The melt envelope, liquid limits, melt phase, and switches, read-only. |
| `fs_visc_power_slope`, `fs_visc_log10_at_solidus`, `fs_shear_power_slope`, `fs_shear_log10_at_solidus` | `float` | The Spohn model's parameters, read-only. |
| `crit_melt_frac`, `crit_melt_frac_width`, `hn_visc_slope_1`, `hn_visc_falloff_slope`, `hn_shear_param_1`, `hn_shear_falloff_slope` | `float` | The Henning model's parameters, read-only. |
| `model_name` | `str` | The resolved model name (`off`, `spohn`, `henning`). |
| `get_config_dict()` | `dict` | `model` plus every parameter the model carries, under the config keys from the table above. |
| `save_config(path)` | - | That dict written as TOML. |

`make_partial_melt(model_name, config=None)` resolves a name or alias case-insensitively; absent keys fall back to the model defaults, and both an unrecognized name and a key that no partial-melt model reads raise `ValueError`. Note that the configuration keys for the melt envelope carry their units (`solidus_k`, `liquidus_k`, `liquid_shear_pa`, `liquid_viscosity_pas`, `liquid_bulk_modulus_pa`, `liquid_density_kg_m3`), matching the TOML the world builder reads, while the constructor keywords do not. The model-specific parameters use one name everywhere: constructor keyword, configuration key, and property.

### Attaching a Melt Model to a `Layer`

```python
from TidalPy.Material.eos import ConstantDensityEOS
from TidalPy.PartialMelt import make_partial_melt
from TidalPy.Structures.layers.physics import PhysicsLayer

mantle = PhysicsLayer("mantle", 0, 0.0, 1.0e6, 2.1e19)
mantle.set_eos(ConstantDensityEOS(shear_modulus_static=50.0e9, bulk_modulus_static=100.0e9))

mantle.set_partial_melt(make_partial_melt("henning", {"solidus_k": 1500.0}))
```

A partial-melt model belongs to the layer's material, which is its EOS model: the layer's `set_partial_melt` is a helper that hands the model to the attached EOS (so attach the EOS first), and the same method is on the EOS model itself. The world's equation-of-state solve applies the model as it integrates, to the shear pair and (when switched on) to the density, the bulk modulus, and the bulk viscosity, and `get_melt_fraction(radius)` reads the result back. The declarative form is a `[layers.<name>.material.partial_melt]` table in the world's TOML; see the [TOML schema](../Structures/config/toml_schema.md).

## C++ API

`c_PartialMeltConfig` (in `partial_melt_base_.hpp`) carries every parameter for every model in one struct with the defaults listed above. The call passes `c_PartialMeltInputs { temperature, premelt_viscosity, premelt_shear }` and returns `c_PartialMeltResult { melt_fraction, postmelt_viscosity, postmelt_shear_modulus }`.

`c_PartialMeltBase : c_PhysicsBase` (in `partial_melt_base_.hpp`) holds the envelope and liquid limits, implements `calc_melt_fraction(temperature)` and `calc_bulk_modulus_melt(temperature, premelt_bulk, framework_shear)` for every model, declares `calc_partial_melt(const c_PartialMeltInputs&) const` pure virtual, and adds `calc_partial_melt_vectorize(temperature, premelt_viscosity, premelt_shear, out_results)`, a radial sweep. Accessors `get_solidus`, `get_liquidus`, `get_liquid_shear`, `get_liquid_viscosity`, `get_bulk_melt_weakening`, and `get_liquid_bulk_modulus` are on the base.

The concrete models `c_OffPartialMelt`, `c_SpohnPartialMelt`, and `c_HenningPartialMelt` live in `partial_melt_.hpp` with a getter per parameter. The factory follows the same shape as the other physics modules: `c_partial_melt_model_from_name(name)` maps onto the `c_PartialMeltModel` enum and throws `std::invalid_argument` for an unknown name, `c_find_partial_melt(model, config)` returns a `std::unique_ptr<c_PartialMeltBase>` (a name overload does both), and `c_partial_melt_from_binary(stream, force=false)` peeks the class id, builds, and reads.

## Adding a New Model

1. Add the model's parameters to `c_PartialMeltConfig` in `partial_melt_base_.hpp`, with defaults.
2. Add `c_<Name>PartialMelt : c_PartialMeltBase` implementing `calc_partial_melt`, `append_config_entries`, `get_binary_class_id`, and `get_binary_params` / `set_binary_params` appending the model's parameters to the base's ([Binary Serialization](../Utilities/binary.md)).
3. Reserve the next free `BinaryClassID` in the 700 block in `Utilities/binary/binary_.hpp`.
4. Register a `c_PartialMeltModel::<Name>` enum value and wire it into `c_partial_melt_model_from_name`, `c_find_partial_melt`, and `c_partial_melt_from_binary`.
5. Add the Cython `cdef class` with its parameter properties, the adoption branch in `make_partial_melt`, tests in `Tests/Test_PartialMelt/`, and an entry on this page.

## References

- Fischer, H.-J., and Spohn, T. (1990). Thermal-orbital histories of viscoelastic models of Io. *Icarus*, 83(1), 39-65.
- Hashin, Z., and Shtrikman, S. (1963). A variational approach to the theory of the elastic behaviour of multiphase materials. *Journal of the Mechanics and Physics of Solids*, 11(2), 127-140.
- Henning, W. G., O'Connell, R. J., and Sasselov, D. D. (2009). Tidally heated terrestrial exoplanets: Viscoelastic response models. *The Astrophysical Journal*, 707(2), 1000-1015.
- Mavko, G. M. (1980). Velocity and attenuation in partially molten rocks. *Journal of Geophysical Research*, 85(B10), 5173-5189.
- Renaud, J. P., and Henning, W. G. (2018). Increased tidal dissipation using advanced rheological models: Implications for Io and tidally active exoplanets. *The Astrophysical Journal*, 857(2), 98.
- Takei, Y. (2002). Effect of pore geometry on VP/VS: From equilibrium geometry to crack. *Journal of Geophysical Research*, 107(B2), 2043.
- Wood, A. B. (1955). *A Textbook of Sound*. G. Bell and Sons.
