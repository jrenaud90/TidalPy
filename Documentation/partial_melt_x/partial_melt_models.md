# Partial-Melt Models (`partial_melt_x`)

_Updated: 2026-09-12_

A partial-melt model maps a material's pre-melt (solid) viscosity and shear modulus, together with its temperature, onto the post-melt viscosity and shear modulus, and reports the volumetric melt fraction it used to get there. In short, it applies melt weakening.

These quantities depend only on the temperature and pressure state fixed by the equation-of-state solve, not on the forcing frequency, so they are computed once per solve and cached. Only the downstream [rheology](../rheology_x/index.md) step, which produces the complex modulus, is recomputed for each tidal mode.

The math mirrors the validated classic implementation in `TidalPy/rheology/partial_melt/melting_models.py`.

## Melt fraction

The volumetric melt fraction is model-independent. It is the position of the temperature within the material's melting envelope, clipped to a physical range:

$$\phi = \mathrm{clip}\left( \frac{T - T_\mathrm{solidus}}{T_\mathrm{liquidus} - T_\mathrm{solidus}},\; 0,\; 1 \right)$$

Below the solidus $\phi = 0$, above the liquidus $\phi = 1$, and a degenerate envelope with the solidus at or above the liquidus returns $\phi = 0$, which is the fully solid answer. This linear form is a convenience, not a claim about melting thermodynamics; real silicate systems melt non-linearly across the envelope, and the models below are what carry the physics.

## The three models

| Model | Aliases | Behavior |
|---|---|---|
| `OffPartialMelt` | `off`, `none` | No weakening. Post-melt strength equals pre-melt strength, but the melt fraction is still reported. |
| `SpohnPartialMelt` | `spohn`, `fischer`, `fischer_spohn` | Fischer and Spohn (1990). Post-melt viscosity and shear modulus depend on temperature alone, not on the pre-melt values. |
| `HenningPartialMelt` | `henning` | Henning et al. (2009) and Renaud and Henning (2018). Three regimes separated by a critical melt fraction. |

### Off

The melt fraction is computed and returned, and the strengths pass through untouched. Use it to isolate the effect of melt weakening by turning it off, or for a layer you have reason to believe stays below its solidus.

### Spohn (Fischer and Spohn 1990)

$$\eta_\mathrm{post} = 10^{\,(s_\eta / T) - p_\eta}, \qquad \mu_\mathrm{post} = 10^{\,(s_\mu / T) - p_\mu}$$

with the slopes $s$ and phases $p$ given by `fs_visc_power_slope`, `fs_visc_power_phase`, `fs_shear_power_slope`, and `fs_shear_power_phase`. Both results are floored at the liquid limits, the supplied liquid viscosity and the model's `liquid_shear`.

This model overwrites rather than weakens: the pre-melt viscosity and shear modulus do not appear on the right-hand side. That makes it self-contained and cheap, and it makes the layer's viscosity model irrelevant wherever this model is active.

### Henning (2009, 2018)

Three regimes in the melt fraction, with a transition band running from `crit_melt_frac` to `crit_melt_frac + crit_melt_frac_width`. Write $\phi_c$ for the critical melt fraction and $T_\mathrm{break} = T_\mathrm{solidus} + \phi_c (T_\mathrm{liquidus} - T_\mathrm{solidus})$ for the temperature at which it is reached.

| Regime | Post-melt viscosity | Post-melt shear modulus |
|---|---|---|
| $\phi \le 0$ | $\eta_\mathrm{pre}$ | $\mu_\mathrm{pre}$ |
| $0 < \phi < \phi_c$ | $\eta_\mathrm{pre} \exp(-a_\eta \phi)$ | $\mu_\mathrm{pre} \exp(b_1/T - b_2)$ |
| $\phi_c \le \phi \le \phi_c + w$ | $\eta_\mathrm{pre} \exp(-a_\eta \phi_c) \exp(-f_\eta (\phi - \phi_c))$ | $\mu_\mathrm{pre} \exp(b_1/T_\mathrm{break} - b_2) \exp(-f_\mu (\phi - \phi_c))$ |
| $\phi > \phi_c + w$ | $\eta_\mathrm{liquid}$ | `liquid_shear` |

Here $a_\eta$ is `hn_visc_slope_1`, $b_1$ and $b_2$ are `hn_shear_param_1` and `hn_shear_param_2`, and $f_\eta$ and $f_\mu$ are `hn_visc_falloff_slope` and `hn_shear_falloff_slope`. Every branch is floored at the liquid limits.

The structure encodes the disaggregation transition. Below the critical melt fraction, melt sits in isolated pockets and weakens the solid framework gradually. Above it the framework loses contact and the material behaves as a crystal-laden liquid, which is a drop of many orders of magnitude. The breakdown band is a deliberately steep but finite bridge between those two pictures; its width is a numerical convenience that keeps the transition differentiable, not a measured quantity.

## Parameters

| Parameter | Default | Units | Used by |
|---|---|---|---|
| `solidus` | 1600.0 | K | All |
| `liquidus` | 2000.0 | K | All |
| `liquid_shear` | 1.0e-5 | Pa | All |
| `fs_visc_power_slope` | 27000.0 | K | Spohn |
| `fs_visc_power_phase` | 1.0 | — | Spohn |
| `fs_shear_power_slope` | 82000.0 | K | Spohn |
| `fs_shear_power_phase` | 40.6 | — | Spohn |
| `crit_melt_frac` | 0.5 | — | Henning |
| `crit_melt_frac_width` | 0.05 | — | Henning |
| `hn_visc_slope_1` | 13.5 | — | Henning |
| `hn_visc_falloff_slope` | 370.0 | — | Henning |
| `hn_shear_param_1` | 40000.0 | K | Henning |
| `hn_shear_param_2` | 25.0 | — | Henning |
| `hn_shear_falloff_slope` | 700.0 | — | Henning |

`liquid_shear` is the shear modulus assigned to material treated as pure liquid. Its default is small but not exactly zero, and it doubles as the floor every model applies to its post-melt shear modulus. The matching viscosity floor is not a model parameter: it is the `liquid_viscosity` passed in with each call, which the equation-of-state solve currently supplies as the pre-melt viscosity.

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
    liquid_viscosity=0.2,       # Pa s
)

# Name factory: case-insensitive, aliases accepted.
spohn_model = make_partial_melt("fischer", {"solidus_k": 1500.0})
```

Constructors take the melt envelope plus their own parameters, all with the defaults from the table: `OffPartialMelt(solidus=1600.0, liquidus=2000.0, liquid_shear=1.0e-5)`, `SpohnPartialMelt(..., fs_visc_power_slope=27000.0, fs_visc_power_phase=1.0, fs_shear_power_slope=82000.0, fs_shear_power_phase=40.6)`, and `HenningPartialMelt(..., crit_melt_frac=0.5, crit_melt_frac_width=0.05, hn_visc_slope_1=13.5, hn_visc_falloff_slope=370.0, hn_shear_param_1=40000.0, hn_shear_param_2=25.0, hn_shear_falloff_slope=700.0)`.

| Member | Returns | Description |
|---|---|---|
| `calc_melt_fraction(temperature)` | `float` | Melt fraction in [0, 1]. |
| `calc_partial_melt(temperature, premelt_viscosity, premelt_shear, liquid_viscosity)` | `(phi, viscosity, shear_modulus)` | Melt fraction, post-melt viscosity [Pa s], post-melt shear modulus [Pa]. |
| `solidus`, `liquidus`, `liquid_shear` | `float` | The melt envelope, read-only. |
| `model_name` | `str` | The resolved model name (`off`, `spohn`, `henning`). |
| `get_config_dict()` | `dict` | `model` plus every parameter the model carries. |
| `save_config(path)` | — | That dict written as TOML. |
| `save_binary(path)` / `load_binary(path, force=False)` | — | TidalPy binary format; see [Binary serialization](../utilities_x/binary_x.md). |

`make_partial_melt(model_name, config=None)` resolves a name or alias case-insensitively; absent keys fall back to the model defaults and an unrecognized name raises `ValueError`. Note that the configuration keys for the melt envelope carry their units (`solidus_k`, `liquidus_k`, `liquid_shear_pa`), matching the TOML the world builder reads, while the constructor keywords do not. The model-specific parameters use the same name in both places.

## Attaching a melt model to a layer

```python
from TidalPy.partial_melt_x import make_partial_melt
from TidalPy.structures_x.layers.physics import PhysicsLayer

mantle = PhysicsLayer("mantle", 0, 0.0, 1.0e6, 2.1e19,
                      shear_modulus_static=50.0e9, bulk_modulus_static=100.0e9)

mantle.set_partial_melt(make_partial_melt("henning", {"solidus_k": 1500.0}))
```

Ownership of the C++ model transfers into the layer. During the world's equation-of-state solve the model is applied at every radial slice, first to the shear pair and then to the bulk pair. The declarative form is a `[layers.<name>.partial_melt]` table in the world's TOML; see the [TOML schema](../structures_x/config/toml_schema.md).

## C++ API

`c_PartialMeltConfig` carries every parameter for every model in one struct with the defaults listed above. The call passes `c_PartialMeltInputs { temperature, premelt_viscosity, premelt_shear, liquid_viscosity }` and returns `c_PartialMeltResult { melt_fraction, postmelt_viscosity, postmelt_shear_modulus }`.

`c_PartialMeltBase : c_PhysicsBase` (in `partial_melt_base_.hpp`) implements `calc_melt_fraction(temperature)` for every model, declares `calc_partial_melt(const c_PartialMeltInputs&) const` pure virtual, and adds `calc_partial_melt_vectorize(temperature, premelt_viscosity, premelt_shear, liquid_viscosity, out_results)`, the radial sweep the equation-of-state solve uses. Accessors `get_solidus`, `get_liquidus`, and `get_liquid_shear` are on the base.

The concrete models `c_OffPartialMelt`, `c_SpohnPartialMelt`, and `c_HenningPartialMelt` live in `partial_melt_.hpp` with a getter per parameter, and occupy binary class ids 701 through 703. The factory follows the same shape as the other physics modules: `c_partial_melt_model_from_name(name)` maps onto the `c_PartialMeltModel` enum and throws `std::invalid_argument` for an unknown name, `c_find_partial_melt(model, config)` returns a `std::unique_ptr<c_PartialMeltBase>` (a name overload does both), and `c_partial_melt_from_binary(stream, force=false)` peeks the class id, builds, and reads.

## Adding a new model

1. Add the model's parameters to `c_PartialMeltConfig` in `partial_melt_.hpp`, with defaults.
2. Add `c_<Name>PartialMelt : c_PartialMeltBase` implementing `calc_partial_melt`, `append_config_entries`, `write_binary`, and `read_binary`.
3. Reserve the next free `BinaryClassID` in the 700 block in `Utilities_x/binary_x/binary_.hpp`.
4. Register a `c_PartialMeltModel::<Name>` enum value and wire it into `c_partial_melt_model_from_name`, `c_find_partial_melt`, and `c_partial_melt_from_binary`.
5. Add the Cython `cdef class` with its parameter properties, the adoption branch in `make_partial_melt`, tests in `Tests/Test_PartialMelt_x/`, and an entry on this page.

## References

- Fischer, H.-J., and Spohn, T. (1990). Thermal-orbital histories of viscoelastic models of Io. *Icarus*, 83(1), 39-65.
- Henning, W. G., O'Connell, R. J., and Sasselov, D. D. (2009). Tidally heated terrestrial exoplanets: Viscoelastic response models. *The Astrophysical Journal*, 707(2), 1000-1015.
- Renaud, J. P., and Henning, W. G. (2018). Increased tidal dissipation using advanced rheological models: Implications for Io and tidally active exoplanets. *The Astrophysical Journal*, 857(2), 98.
