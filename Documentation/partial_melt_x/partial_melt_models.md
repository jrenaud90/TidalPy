# Partial-Melt Models (`partial_melt_x`)

_Updated: 2026-09-23_

A partial-melt model maps a material's pre-melt (solid) viscosity and shear modulus, together with its temperature, onto the post-melt viscosity and shear modulus, and reports the volumetric melt fraction it used to get there. Every model can also weaken the bulk modulus with its own, much weaker law; that is off by default (see Bulk Modulus).

These quantities depend only on the temperature and pressure state fixed by the equation-of-state solve, not on the forcing frequency, so they are computed once per solve and cached. Only the downstream [rheology](../rheology_x/index.md) step, which produces the complex modulus (frequency-dependent), is recomputed for each tidal mode.

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

$$\eta_\mathrm{post} = 10^{\,(s_\eta / T) - p_\eta}, \qquad \mu_\mathrm{post} = 10^{\,(s_\mu / T) - p_\mu}$$

with the slopes $s$ and phases $p$ given by `fs_visc_power_slope`, `fs_visc_power_phase`, `fs_shear_power_slope`, and `fs_shear_power_phase`. Both results are floored at the liquid limits, `liquid_viscosity` and `liquid_shear`.

The law applies above the solidus ($\phi > 0$). At or below it the pre-melt viscosity and shear modulus are returned unchanged; the law itself grows without bound as the temperature falls (10$^{41}$ Pa at 1000 K with the defaults). Above the solidus the pre-melt values do not appear on the right-hand side, so the layer's viscosity model has no effect wherever the material is partially molten.

### Henning (2009, 2018)

Three regimes in the melt fraction, with a transition band running from `crit_melt_frac` to `crit_melt_frac + crit_melt_frac_width`. Write $\phi_c$ for the critical melt fraction and $T_\mathrm{break} = T_\mathrm{solidus} + \phi_c (T_\mathrm{liquidus} - T_\mathrm{solidus})$ for the temperature at which it is reached.

| Regime | Post-melt viscosity | Post-melt shear modulus |
|---|---|---|
| $\phi \le 0$ | $\eta_\mathrm{pre}$ | $\mu_\mathrm{pre}$ |
| $0 < \phi < \phi_c$ | $\eta_\mathrm{pre} \exp(-a_\eta \phi)$ | $\mu_\mathrm{pre} \exp(b_1/T - b_2)$ |
| $\phi_c \le \phi \le \phi_c + w$ | $\eta_\mathrm{pre} \exp(-a_\eta \phi_c) \exp(-f_\eta (\phi - \phi_c))$ | $\mu_\mathrm{pre} \exp(b_1/T_\mathrm{break} - b_2) \exp(-f_\mu (\phi - \phi_c))$ |
| $\phi > \phi_c + w$ | $\eta_\mathrm{liquid}$ | `liquid_shear` |

Here $a_\eta$ is `hn_visc_slope_1`, $b_1$ and $b_2$ are `hn_shear_param_1` and `hn_shear_param_2`, and $f_\eta$ and $f_\mu$ are `hn_visc_falloff_slope` and `hn_shear_falloff_slope`. Every branch is floored at the liquid limits.

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
| `fs_visc_power_phase` | `fs_visc_power_phase` | 1.0 | - | Spohn |
| `fs_shear_power_slope` | `fs_shear_power_slope_k` | 82000.0 | K | Spohn |
| `fs_shear_power_phase` | `fs_shear_power_phase` | 40.6 | - | Spohn |
| `crit_melt_frac` | `crit_melt_frac` | 0.5 | - | Henning |
| `crit_melt_frac_width` | `crit_melt_frac_width` | 0.05 | - | Henning |
| `hn_visc_slope_1` | `hn_visc_slope_1` | 13.5 | - | Henning |
| `hn_visc_falloff_slope` | `hn_visc_falloff_slope` | 370.0 | - | Henning |
| `hn_shear_param_1` | `hn_shear_param_1_k` | 40000.0 | K | Henning |
| `hn_shear_param_2` | `hn_shear_param_2` | 25.0 | - | Henning |
| `hn_shear_falloff_slope` | `hn_shear_falloff_slope` | 700.0 | - | Henning |

`liquid_shear` and `liquid_viscosity` are the shear modulus and viscosity assigned to material treated as pure liquid, and they are also the floors every model applies to its post-melt shear modulus and viscosity. The shear default is small but not exactly zero. The viscosity default is a molten silicate's; the packaged configuration gives each material type its own (rock 0.2, iron 1.3e-2, ice and high-pressure ice 8.9e-4 Pa s).

## Bulk Modulus

Melt lowers the bulk modulus far less than the shear modulus: a melt keeps a bulk modulus of the same order as the solid's (about 20 GPa for a silicate melt at low pressure against about 130 GPa for the rock), while its shear modulus vanishes (Mavko 1980; Takei 2002). The partial-melt models therefore leave the bulk modulus unchanged by default. With `bulk_melt_weakening = true` the post-melt bulk modulus is the Hashin and Shtrikman (1963) bound for melt of bulk modulus $K_l$ (`liquid_bulk_modulus`) in a solid framework of bulk modulus $K_s$ (the pre-melt value), evaluated with the framework's post-melt shear modulus $\mu$:

$$K = K_s + \frac{\phi}{\dfrac{1}{K_l - K_s} + \dfrac{1 - \phi}{K_s + \tfrac{4}{3}\mu}}$$

While the framework holds, $\mu$ is close to the solid's and this is the upper bound for isolated melt pockets, a weak reduction (about 16 percent at $\phi = 0.1$ for $K_s$ = 130, $\mu$ = 60, $K_l$ = 20 GPa). Once the melt model has collapsed the framework's shear modulus (Henning past the critical melt fraction), the same expression becomes the Reuss (Wood 1955) average of a crystal suspension, and it reaches $K_l$ at $\phi = 1$. The bulk viscosity is not changed by melt.

| Parameter | Config key | Default | Units | Used by |
|---|---|---|---|---|
| `bulk_melt_weakening` | `bulk_melt_weakening` | false | - | All |
| `liquid_bulk_modulus` | `liquid_bulk_modulus_pa` | 2.0e10 | Pa | All, when `bulk_melt_weakening` is true |

The packaged configuration gives liquid iron 1.1e11 Pa and liquid water 2.2e9 Pa.

## Python API

```python
from TidalPy.partial_melt_x import (
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

`OffPartialMelt(solidus=1600.0, liquidus=2000.0, liquid_shear=1.0e-5, liquid_viscosity=0.2, bulk_melt_weakening=False, liquid_bulk_modulus=2.0e10)`

`SpohnPartialMelt(solidus, liquidus, liquid_shear, fs_visc_power_slope=27000.0, fs_visc_power_phase=1.0, fs_shear_power_slope=82000.0, fs_shear_power_phase=40.6, liquid_viscosity, bulk_melt_weakening, liquid_bulk_modulus)`

`HenningPartialMelt(solidus, liquidus, liquid_shear, crit_melt_frac=0.5, crit_melt_frac_width=0.05, hn_visc_slope_1=13.5, hn_visc_falloff_slope=370.0, hn_shear_param_1=40000.0, hn_shear_param_2=25.0, hn_shear_falloff_slope=700.0, liquid_viscosity, bulk_melt_weakening, liquid_bulk_modulus)`

| Member | Returns | Description |
|---|---|---|
| `calc_melt_fraction(temperature)` | `float` | Melt fraction in [0, 1]. |
| `calc_partial_melt(temperature, premelt_viscosity, premelt_shear)` | `(phi, viscosity, shear_modulus)` | Melt fraction, post-melt viscosity [Pa s], post-melt shear modulus [Pa]. |
| `calc_bulk_modulus_melt(temperature, premelt_bulk_modulus, framework_shear_modulus)` | `float` | Post-melt bulk modulus [Pa]; the pre-melt value unless `bulk_melt_weakening` is on. |
| `solidus`, `liquidus`, `liquid_shear`, `liquid_viscosity`, `bulk_melt_weakening`, `liquid_bulk_modulus` | `float`, `bool` | The melt envelope and liquid limits, read-only. |
| `fs_visc_power_slope`, `fs_visc_power_phase`, `fs_shear_power_slope`, `fs_shear_power_phase` | `float` | The Spohn model's parameters, read-only. |
| `crit_melt_frac`, `crit_melt_frac_width`, `hn_visc_slope_1`, `hn_visc_falloff_slope`, `hn_shear_param_1`, `hn_shear_param_2`, `hn_shear_falloff_slope` | `float` | The Henning model's parameters, read-only. |
| `model_name` | `str` | The resolved model name (`off`, `spohn`, `henning`). |
| `get_config_dict()` | `dict` | `model` plus every parameter the model carries, under the config keys from the table above. |
| `save_config(path)` | - | That dict written as TOML. |
| `save_binary(path)` / `load_binary(path, force=False)` | - | TidalPy binary format; see [Binary serialization](../utilities_x/binary_x.md). |

`make_partial_melt(model_name, config=None)` resolves a name or alias case-insensitively; absent keys fall back to the model defaults, and both an unrecognized name and a key that no partial-melt model reads raise `ValueError`. Note that the configuration keys for the melt envelope carry their units (`solidus_k`, `liquidus_k`, `liquid_shear_pa`, `liquid_viscosity_pas`, `liquid_bulk_modulus_pa`), matching the TOML the world builder reads, while the constructor keywords do not. The model-specific parameters use one name everywhere: constructor keyword, configuration key, and property.

### Attaching a Melt Model to a `Layer`

```python
from TidalPy.Material_x.eos import ConstantDensityEOS
from TidalPy.partial_melt_x import make_partial_melt
from TidalPy.structures_x.layers.physics import PhysicsLayer

mantle = PhysicsLayer("mantle", 0, 0.0, 1.0e6, 2.1e19)
mantle.set_eos(ConstantDensityEOS(shear_modulus_static=50.0e9, bulk_modulus_static=100.0e9))

mantle.set_partial_melt(make_partial_melt("henning", {"solidus_k": 1500.0}))
```

A partial-melt model belongs to the layer's material, which is its EOS model: the layer's `set_partial_melt` is a helper that hands the model to the attached EOS (so attach the EOS first), and the same method is on the EOS model itself. The world's equation-of-state solve applies the model as it integrates, to the shear pair and (when switched on) to the bulk modulus, and `get_melt_fraction(radius)` reads the result back. The declarative form is a `[layers.<name>.material.partial_melt]` table in the world's TOML; see the [TOML schema](../structures_x/config/toml_schema.md).

## C++ API

`c_PartialMeltConfig` (in `partial_melt_base_.hpp`) carries every parameter for every model in one struct with the defaults listed above. The call passes `c_PartialMeltInputs { temperature, premelt_viscosity, premelt_shear }` and returns `c_PartialMeltResult { melt_fraction, postmelt_viscosity, postmelt_shear_modulus }`.

`c_PartialMeltBase : c_PhysicsBase` (in `partial_melt_base_.hpp`) holds the envelope and liquid limits, implements `calc_melt_fraction(temperature)` and `calc_bulk_modulus_melt(temperature, premelt_bulk, framework_shear)` for every model, declares `calc_partial_melt(const c_PartialMeltInputs&) const` pure virtual, and adds `calc_partial_melt_vectorize(temperature, premelt_viscosity, premelt_shear, out_results)`, a radial sweep. Accessors `get_solidus`, `get_liquidus`, `get_liquid_shear`, `get_liquid_viscosity`, `get_bulk_melt_weakening`, and `get_liquid_bulk_modulus` are on the base, and every model's binary payload starts with those six values.

The concrete models `c_OffPartialMelt`, `c_SpohnPartialMelt`, and `c_HenningPartialMelt` live in `partial_melt_.hpp` with a getter per parameter, and occupy binary class ids 701 through 703. The factory follows the same shape as the other physics modules: `c_partial_melt_model_from_name(name)` maps onto the `c_PartialMeltModel` enum and throws `std::invalid_argument` for an unknown name, `c_find_partial_melt(model, config)` returns a `std::unique_ptr<c_PartialMeltBase>` (a name overload does both), and `c_partial_melt_from_binary(stream, force=false)` peeks the class id, builds, and reads.

## Adding a New Model

1. Add the model's parameters to `c_PartialMeltConfig` in `partial_melt_base_.hpp`, with defaults.
2. Add `c_<Name>PartialMelt : c_PartialMeltBase` implementing `calc_partial_melt`, `append_config_entries`, `write_binary`, and `read_binary`.
3. Reserve the next free `BinaryClassID` in the 700 block in `Utilities_x/binary_x/binary_.hpp`.
4. Register a `c_PartialMeltModel::<Name>` enum value and wire it into `c_partial_melt_model_from_name`, `c_find_partial_melt`, and `c_partial_melt_from_binary`.
5. Add the Cython `cdef class` with its parameter properties, the adoption branch in `make_partial_melt`, tests in `Tests/Test_PartialMelt_x/`, and an entry on this page.

## References

- Fischer, H.-J., and Spohn, T. (1990). Thermal-orbital histories of viscoelastic models of Io. *Icarus*, 83(1), 39-65.
- Hashin, Z., and Shtrikman, S. (1963). A variational approach to the theory of the elastic behaviour of multiphase materials. *Journal of the Mechanics and Physics of Solids*, 11(2), 127-140.
- Henning, W. G., O'Connell, R. J., and Sasselov, D. D. (2009). Tidally heated terrestrial exoplanets: Viscoelastic response models. *The Astrophysical Journal*, 707(2), 1000-1015.
- Mavko, G. M. (1980). Velocity and attenuation in partially molten rocks. *Journal of Geophysical Research*, 85(B10), 5173-5189.
- Renaud, J. P., and Henning, W. G. (2018). Increased tidal dissipation using advanced rheological models: Implications for Io and tidally active exoplanets. *The Astrophysical Journal*, 857(2), 98.
- Takei, Y. (2002). Effect of pore geometry on VP/VS: From equilibrium geometry to crack. *Journal of Geophysical Research*, 107(B2), 2043.
- Wood, A. B. (1955). *A Textbook of Sound*. G. Bell and Sons.
