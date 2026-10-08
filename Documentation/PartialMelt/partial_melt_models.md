# Melting Laws (`PartialMelt`)

_Updated: 2026-10-07_

The melting laws describe how a [material](../Material/materials.md) melts. Two melting curves give the solidus and liquidus \[K\] at a pressure \[Pa\], and the melt fraction $\phi$ runs linearly between them. A melt-weakening law gives the partially molten aggregate's shear modulus \[Pa\] and viscosity \[Pa s\] from the solid's and the liquid's. Optional bulk-mixing laws give its bulk modulus \[Pa\] and bulk viscosity \[Pa s\]. The material combines the laws, which know nothing of the phases' equations of state.

| Family | Model | Aliases | Python class |
|---|---|---|---|
| Melting curve | `constant` | `const` | `ConstantMeltingCurve` |
| | `simon_glatzel` | `simon-glatzel` | `SimonGlatzelCurve` |
| | `simon_glatzel_2` | `simon-glatzel-2` | `SimonGlatzel2Curve` |
| | `interpolate` | `interp`, `interpolated` | `InterpolatedMeltingCurve` |
| Melt weakening | `none` | `off` | `NoMeltWeakening` |
| | `spohn` | `fischer`, `fischer_spohn` | `SpohnMeltWeakening` |
| | `henning` | | `HenningMeltWeakening` |
| Bulk-modulus mixing | `hashin_shtrikman` | `hs`, `hashin-shtrikman` | `HashinShtrikmanMixing` |
| Bulk-viscosity mixing | `compaction` | `mckenzie` | `CompactionViscosity` |

## Quick Example

```python
import numpy as np

from TidalPy.PartialMelt import (
    HenningMeltWeakening,
    SimonGlatzelCurve,
    make_bulk_modulus_mixing,
    make_melt_weakening,
    make_melting_curve)

# Peridotite's solidus (Monteux et al. 2016) by name, with its parameters by config key
solidus = make_melting_curve(
    "simon_glatzel_2",
    {"temperature_k": 1661.2,
     "simon_a_pa": 1.336e9,
     "simon_c": 7.437,
     "transition_pressure_pa": 20.0e9,
     "high_temperature_k": 2081.8,
     "high_simon_a_pa": 1.0169e11,
     "high_simon_c": 1.226})
print(solidus.calc_melting_temperature(135.0e9))   # [K] about 4150 at Earth's core-mantle boundary
print(solidus.calc_melting_slope(5.0e9))           # [K Pa-1] dT_m/dP at 5 GPa

# Ice Ih's melting curve falls with pressure, held past its triple point with ice III
ice_melting = SimonGlatzelCurve(
    temperature=273.16,
    simon_a=-4.15e8,
    simon_c=8.25,
    maximum_pressure=2.08566e8)
print(ice_melting.calc_melting_temperature(np.array([0.0, 1.0e8, 2.0e8, 5.0e8])))   # [K]

# The aggregate's shear modulus and viscosity a quarter of the way through the melting range
weakening = HenningMeltWeakening()
shear_modulus, viscosity = weakening.calc_weakening(
    temperature=1700.0,
    solidus=1600.0,
    liquidus=2000.0,
    solid_shear=6.0e10,        # [Pa]
    solid_viscosity=1.0e20,    # [Pa s]
    liquid_shear=0.0,
    liquid_viscosity=0.1)

# Fischer and Spohn by an alias
spohn = make_melt_weakening("fischer")

# The bulk modulus of rock with 10 percent melt
mixing = make_bulk_modulus_mixing("hs")
print(mixing.calc_bulk_modulus(1.3e11, 2.0e10, 6.0e10, 0.1))   # [Pa], about 84 percent of the rock's
```

## Attaching Melting Laws to a `Material`

```python
from TidalPy.Material import Material, Phase
from TidalPy.PartialMelt import make_melt_weakening, make_melting_curve

rock = Phase(
    eos={"model": "constant", "reference_density_kg_m3": 3300.0},
    shear_modulus={"model": "constant", "shear_modulus_pa": 6.0e10},
    shear_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e20})
melt = Phase(
    eos={"model": "murnaghan", "reference_density_kg_m3": 2750.0},
    shear_viscosity={"model": "constant", "reference_viscosity_pas": 0.1})

melting_rock = Material(
    solid=rock,
    liquid=melt,
    solidus=make_melting_curve("constant", {"temperature_k": 1600.0}),
    liquidus=make_melting_curve("constant", {"temperature_k": 2000.0}),
    weakening=make_melt_weakening("henning"),
    bulk_modulus_mixing="hashin_shtrikman",   # A name builds the law at its defaults
    bulk_viscosity_mixing={"model": "compaction", "exponent": 1.0})

state = melting_rock.calc_state(
    1.0e9,
    1700.0,
    use_melting=True)   # A layer passes its own switches
print(state["melt_fraction"], state["shear_modulus"], state["bulk_viscosity"])
```

A material shares its laws rather than copying them. In a MatPack file or a world TOML they go in the material's `melting` table:

```toml
# Henning weakening with the peridotite melting curves of Monteux et al. (2016)
[layers.mantle.material.melting.solidus]
model = "simon_glatzel_2"
temperature_k = 1661.2
simon_a_pa = 1.336e9
simon_c = 7.437
transition_pressure_pa = 20.0e9
high_temperature_k = 2081.8
high_simon_a_pa = 1.0169e11
high_simon_c = 1.226

[layers.mantle.material.melting.liquidus]
model = "simon_glatzel_2"
temperature_k = 1982.1
simon_a_pa = 6.594e9
simon_c = 5.374
transition_pressure_pa = 20.0e9
high_temperature_k = 78.74
high_simon_a_pa = 4.054e6
high_simon_c = 2.44

[layers.mantle.material.melting.weakening]
model = "henning"
```

The layer's `use_melting` and `use_pressure_melting` decide whether it uses them. See [Phases and Materials](../Material/materials.md) and the [TOML schema](../Structures/config/toml_schema.md).

## Melting Curves

The solidus and liquidus are independent curves. A layer reads them at the local pressure when it sets `use_pressure_melting`, and at zero pressure otherwise; the melt fraction, the weakening laws (anchored at the local solidus), the density mixing, and the bulk effects all use that reading.

`constant` gives $T_m(P) = T_0$ \[K\] at every pressure. `interpolate` is linear in pressure between the points of a `pressure` and `temperature` table, held at the end values beyond it.

### Simon and Glatzel

The Simon and Glatzel (1929) law fits most planetary melting data:

$$T_m(P) = T_0 \left(1 + \frac{P - P_\mathrm{ref}}{a}\right)^{1/c},$$

with $T_0$ the melting temperature at the reference pressure $P_\mathrm{ref}$ (usually zero), $a$ \[Pa\] the pressure at which the curve begins to rise, and $c$ how quickly the rise flattens. Below $P_\mathrm{ref}$ the curve holds $T_0$; above `maximum_pressure_pa` (default: none), the end of the fitted range, it holds its value there.

A negative $a$ gives a curve that falls with pressure, as ice Ih's does (MatPack's `ice_ih`: $T_0$ = 273.16 K, $a$ = -415 MPa, $c$ = 8.25). A falling curve reaches 0 K at $P_\mathrm{ref} - a$, so it needs a `maximum_pressure_pa` below that; `ice_ih` ends at its triple point with ice III and liquid water, 208.566 MPa, holding about 251 K.

`simon_glatzel_2` joins two branches at a transition pressure $P_t$, both in absolute pressure as the published fits are, so a fit can change slope across a phase transition (such as the one near 20 GPa at the top of Earth's lower mantle):

$$T_m(P) = \begin{cases} T_0 \left(1 + (P - P_\mathrm{ref}) / a\right)^{1/c} & P \le P_t \\ T_{0,\mathrm{high}} \left(1 + (P - P_\mathrm{ref,high}) / a_\mathrm{high}\right)^{1/c_\mathrm{high}} & P > P_t \end{cases}$$

### Peridotite Fits

Monteux et al. (2016) fit two-branch laws to the peridotite and chondritic-mantle experiments of Fiquet et al. (2010) and Andrault et al. (2011). MatPack's `peridotite` and `lower_mantle` use them, and the solidus fit is the `simon_glatzel_2` default.

| Curve | $T_0$ \[K\] | $a$ \[Pa\] | $c$ | $P_t$ \[Pa\] | $T_{0,\mathrm{high}}$ \[K\] | $a_\mathrm{high}$ \[Pa\] | $c_\mathrm{high}$ |
|---|---|---|---|---|---|---|---|
| Solidus | 1661.2 | 1.336e9 | 7.437 | 20.0e9 | 2081.8 | 1.0169e11 | 1.226 |
| Liquidus | 1982.1 | 6.594e9 | 5.374 | 20.0e9 | 78.74 | 4.054e6 | 2.44 |

| Pressure \[GPa\] | Solidus \[K\] | Liquidus \[K\] | Where |
|---|---|---|---|
| 0 | 1661 | 1982 | Surface |
| 5 | 2048 | 2202 | Base of Io's mantle |
| 20 | 2411 | 2569 | Branch transition |
| 60 | 3039 | 4030 | Mid lower mantle of the Earth |
| 135 | 4147 | 5619 | Earth's core-mantle boundary |

The experiments reach about 140 GPa, so the curves are extrapolated in the deeper mantles of planets larger than the Earth.

### Melting Slope

`calc_melting_slope(pressure)` returns $dT_m/dP$ \[K Pa$^{-1}$\]: $T_0 / (a c) \, (1 + (P - P_\mathrm{ref}) / a)^{1/c - 1}$ for a Simon and Glatzel branch, and the interval's slope for a table. It is zero for a constant curve and wherever a curve is held flat. A material uses the slopes for the latent heat's share of an adiabat's expansivity inside a melting range (see [Latent Heat](../Material/materials.md#latent-heat)).

### Choosing Melting Curves

Use pressure-dependent curves, with `use_pressure_melting` on, for a rocky layer whose base is deeper than a few GPa: at the base of Io's mantle (about 5 GPa) the peridotite solidus is already about 390 K above its zero-pressure value. The bundled worlds with warm silicate mantles do so. Pair them with the [compressed expansivity](../Material/material_eos.md#expansivity-under-compression): with a constant expansivity an Earth-like mantle's adiabat climbs above the peridotite solidus at depth. A constant curve suits a thin layer, a comparison with a published model that used one, or a material without a measured curve.

## Melt Weakening

`none` keeps the solid's shear modulus and viscosity until fully molten. `spohn` (Fischer and Spohn 1990) lowers both with temperature from their values at the solidus. `henning` (Henning et al. 2009; Renaud and Henning 2018) has three regimes split by a critical melt fraction.

Every law returns the solid pair at $\phi \le 0$ and the liquid pair at $\phi \ge 1$, floored at the liquid's values in between. A material without a weakening law behaves as `none`. The material decides from the returned shear modulus where the radial solver treats the aggregate as a liquid. Both temperature laws are anchored at the solidus, so they carry over to materials whose solidus is not the 1600 K of the published silicate fits.

`henning` is the usual choice for a silicate mantle: it weakens the framework gradually until the critical melt fraction, then collapses it. `spohn` reproduces the earlier Io models built on Fischer and Spohn's fits. `none` suits ices and other materials whose partially molten strength is not constrained, and a material that melts as a step.

### Breakdown Band

The Spohn and Henning laws describe the partially molten framework up to a rheological transition, the breakdown band from $\phi_c$ (`crit_melt_frac`) to $\phi_c + w$ (`crit_melt_frac + crit_melt_frac_width`). Across the band the framework's values $\eta_f$, $\mu_f$ blend into the liquid's $\eta_l$, $\mu_l$, reached at the band's end and held past it. With $s = (\phi - \phi_c) / w$:

$$\eta = \eta_f^{\,1 - s} \, \eta_l^{\,s}, \qquad \mu = (1 - s)\,\mu_f + s\,\mu_l$$

The viscosity blends log-linearly (it spans many orders of magnitude) and the shear modulus linearly (a liquid's is usually zero). Both laws are therefore continuous in temperature through the melting range, which a time integration needs: a viscosity step of about $10^{10}$ makes an implicit integrator crawl where the mantle sits on it. A zero `crit_melt_frac_width` makes the transition a step into the liquid at $\phi_c$.

### None

The solid's values hold until full melt, and the melt fraction is still reported. Use it when the strength change is not worth modeling, or to isolate the density and thermal terms.

### Spohn (Fischer and Spohn 1990)

Below the [breakdown band](#breakdown-band),

$$\eta = 10^{\,L_\eta + s_\eta (1/T - 1/T_\mathrm{sol})}, \qquad \mu = 10^{\,L_\mu + s_\mu (1/T - 1/T_\mathrm{sol})}$$

with slopes $s$ \[K\] (`visc_power_slope`, `shear_power_slope`) and $L$ the base-10 logarithms of the strengths at the solidus (`visc_log10_at_solidus`, `shear_log10_at_solidus`). Unset (the default), each $L$ is the solid phase's own value at the solidus and local pressure, so the law continues the solid's values into the melting range without a step.

Fischer and Spohn (1990) fit silicates with $10^{27000/T - 1}$ Pa s and $10^{82000/T - 40.6}$ Pa: the form above at a 1600 K solidus with $L_\eta$ = 15.875 and $L_\mu$ = 10.65. Setting those two reproduces the published fits there, with a step at the solidus from the solid's values to the fit's. The absolute fit would fail elsewhere: at an icy 273 K solidus it gives a shear modulus of $10^{260}$ Pa.

### Henning (2009, 2018)

Three regimes in the melt fraction, split by the [breakdown band](#breakdown-band). With $T_\mathrm{break} = T_\mathrm{sol} + \phi_c (T_\mathrm{liq} - T_\mathrm{sol})$ the temperature where $\phi_c$ is reached, the framework's pair is

| Regime | Viscosity $\eta_f$ | Shear modulus $\mu_f$ |
|---|---|---|
| $\phi \le 0$ | $\eta_s$ | $\mu_s$ |
| $0 < \phi < \phi_c$ | $\eta_s \exp(-a_\eta \phi)$ | $\mu_s \exp[b_1 (1/T - 1/T_\mathrm{sol})]$ |
| $\phi_c \le \phi < \phi_c + w$ | $\eta_s \exp(-a_\eta \phi_c) \exp(-f_\eta (\phi - \phi_c))$ | $\mu_s \exp[b_1 (1/T_\mathrm{break} - 1/T_\mathrm{sol})] \exp(-f_\mu (\phi - \phi_c))$ |

which the band blends into the liquid's pair, held from $\phi_c + w$ on. $\eta_s$ and $\mu_s$ are the solid's values, $a_\eta$ is `visc_slope_1`, $b_1$ is `shear_param_1` \[K\], and $f_\eta$, $f_\mu$ are `visc_falloff_slope`, `shear_falloff_slope`. Every branch is floored at the liquid's values. The falloff alone reaches only about $10^{-11}$ of the solid's viscosity by the band's end (defaults), far above a melt's; the blend closes that gap.

The shear law is Henning et al. (2009) Eq. 20, $\exp(40000/T - 25)$, anchored at the solidus: their 25 is $40000/1600$, their silicate solidus, so a 1600 K solidus reproduces it exactly. The published constant would stiffen an icy layer with a 273 K solidus by about $e^{121}$ ($e^{40000/273 - 25}$) as it began to melt.

Below $\phi_c$ melt sits in isolated pockets; above it the framework loses contact and the material behaves as a crystal-laden liquid (see [Partial Melting](index.md)). The band's width is a numerical convenience that keeps this transition continuous, not a measured quantity.

## Bulk Mixing

Melt lowers the bulk modulus far less than the shear modulus: a silicate melt's bulk modulus (about 20 GPa at low pressure) is of the same order as the rock's (about 130 GPa), while its shear modulus vanishes (Mavko 1980; Takei 2002).

Without a bulk-modulus mixing law, the bulk modulus blends linearly from the solid's to the liquid's across the weakening law's breakdown band, so the mush the radial solver treats as a liquid has the liquid's bulk modulus too. With no weakening law it steps at full melt, with the shear modulus. Without a bulk-viscosity mixing law, the bulk viscosity is the solid's until fully molten, then the liquid's.

### Hashin-Shtrikman

The Hashin and Shtrikman (1963) bound for melt of bulk modulus $K_l$ in a framework of bulk modulus $K_s$ and post-melt shear modulus $\mu$:

$$K = K_s + \frac{\phi}{\dfrac{1}{K_l - K_s} + \dfrac{1 - \phi}{K_s + \tfrac{4}{3}\mu}}$$

The material applies it to both the isothermal and the adiabatic moduli. While the framework holds, $\mu$ is near the solid's and this is the upper bound for isolated melt pockets, a weak reduction (about 16 percent at $\phi = 0.1$ for $K_s$ = 130, $\mu$ = 60, $K_l$ = 20 GPa). Once the weakening law collapses $\mu$ (Henning past $\phi_c$), it becomes the Reuss (Wood 1955) average of a crystal suspension and reaches $K_l$ at $\phi = 1$.

This is the unrelaxed (undrained) modulus, with the melt held in place over a forcing cycle. Melt that wets grain edges or forms films weakens the framework more than isolated pockets do (Mavko 1980; Takei 2002), so it is an upper limit. Its relaxation as melt moves is a bulk rheology's job, at the rate the bulk viscosity sets.

### Compaction Viscosity

Melt-free rock has no viscous compaction; melt adds one as it moves through the matrix. The law adds a matrix bulk viscosity in series with the solid's:

$$\frac{1}{\zeta} = \frac{1}{\zeta_s} + \frac{\phi^n}{c\,\eta},$$

with $\eta$ the post-melt shear viscosity, $c$ `coefficient`, and $n$ `exponent`. $n = 1$ with $c$ of order 1 is the classic compaction viscosity $\eta/\phi$ (McKenzie 1984); $n = 0$ gives a bulk viscosity of the order of the shear viscosity, as micromechanical models do (Takei and Holtzman 2009). The form is continuous at the solidus for $n > 0$. A non-finite or non-positive solid bulk viscosity counts as no solid dashpot, so melt alone sets $\zeta$.

The bulk viscosity reaches the tides only through the layer's bulk rheology, which is elastic unless the layer or its material names one. A Maxwell bulk rheology would relax the bulk modulus to zero at long periods, which is unphysical. A [Zener](../Rheology/rheology_models.md#models) bulk rheology relaxes it to a set fraction $r$ of the unrelaxed modulus, the drained-to-undrained ratio: near 0.9 for isolated pockets at 10% melt, lower with melt films. Bulk dissipation in a partially molten layer, while understudied, has been shown to rival the shear dissipation (Kervazo et al. 2021).

## Parameters

Each parameter has a constructor keyword (also an attribute) and a config key (TOML, factory config dict, `get_config_dict()`), which ends in its unit when dimensional. Either name is accepted.

**Melting curves**

| Parameter | Config key | Default | Units | Used by |
|---|---|---|---|---|
| `temperature` | `temperature_k` | 1600.0 (`simon_glatzel_2` 1661.2; `[1600.0]` for interpolated) | K | All |
| `simon_a` | `simon_a_pa` | 1.0e9 (`simon_glatzel_2` 1.336e9) | Pa | Simon and Glatzel |
| `simon_c` | `simon_c` | 5.0 (`simon_glatzel_2` 7.437) | - | Simon and Glatzel |
| `reference_pressure` | `reference_pressure_pa` | 0.0 | Pa | Simon and Glatzel |
| `maximum_pressure` | `maximum_pressure_pa` | unset (no limit) | Pa | `simon_glatzel` |
| `transition_pressure` | `transition_pressure_pa` | 20.0e9 | Pa | `simon_glatzel_2` |
| `high_temperature` | `high_temperature_k` | 2081.8 | K | `simon_glatzel_2` |
| `high_simon_a` | `high_simon_a_pa` | 1.0169e11 | Pa | `simon_glatzel_2` |
| `high_simon_c` | `high_simon_c` | 1.226 | - | `simon_glatzel_2` |
| `high_reference_pressure` | `high_reference_pressure_pa` | 0.0 | Pa | `simon_glatzel_2` |
| `pressure` | `pressure_pa` | `[0.0]` | Pa | Interpolated |

**Melt weakening**

| Parameter | Config key | Default | Units | Used by |
|---|---|---|---|---|
| `visc_power_slope` | `fs_visc_power_slope_k` | 27000.0 | K | Spohn |
| `visc_log10_at_solidus` | `fs_visc_log10_at_solidus` | unset (the solid's own) | log10 Pa s | Spohn |
| `shear_power_slope` | `fs_shear_power_slope_k` | 82000.0 | K | Spohn |
| `shear_log10_at_solidus` | `fs_shear_log10_at_solidus` | unset (the solid's own) | log10 Pa | Spohn |
| `crit_melt_frac` | `crit_melt_frac` | 0.5 | - | Spohn, Henning |
| `crit_melt_frac_width` | `crit_melt_frac_width` | 0.05 | - | Spohn, Henning |
| `visc_slope_1` | `hn_visc_slope_1` | 13.5 | - | Henning |
| `visc_falloff_slope` | `hn_visc_falloff_slope` | 370.0 | - | Henning |
| `shear_param_1` | `hn_shear_param_1_k` | 40000.0 | K | Henning |
| `shear_falloff_slope` | `hn_shear_falloff_slope` | 700.0 | - | Henning |

**Bulk mixing**

| Parameter | Config key | Default | Units | Used by |
|---|---|---|---|---|
| `coefficient` | `coefficient` | 1.0 | - | Compaction |
| `exponent` | `exponent` | 1.0 | - | Compaction |

Hashin-Shtrikman has no parameters. The liquid's shear modulus and viscosity come from the material's liquid phase (a phase with no shear-modulus law has a shear modulus of 0).

## Python API

Constructors take parameters by name or config key, as keywords or positionally in `get_parameter_info()` order. `make_melting_curve`, `make_melt_weakening`, `make_bulk_modulus_mixing`, and `make_bulk_viscosity_mixing` take `(model_name, config=None)` and resolve a name or alias case-insensitively; absent keys take the defaults. An unknown name, an unread key, or an out-of-bounds value raises `ValueError` naming the closest accepted name or key.

| Member | Returns | Description |
|---|---|---|
| `calc_melting_temperature(pressure)` | `float` or `np.ndarray` \[K\] | Melting temperature; same shape as `pressure`. |
| `calc_melting_slope(pressure)` | `float` or `np.ndarray` \[K Pa$^{-1}$\] | $dT_m/dP$; same shape as `pressure`. |
| `calc_weakening(temperature, solidus, liquidus, solid_shear, solid_viscosity, liquid_shear, liquid_viscosity, melt_fraction=None)` | `(shear_modulus, viscosity)` \[Pa, Pa s\] | Temperatures in K, shear moduli in Pa, viscosities in Pa s, $\phi$ in m$^3$ m$^{-3}$. |
| `calc_bulk_modulus(solid_bulk_modulus, liquid_bulk_modulus, framework_shear_modulus, melt_fraction)` | `float` \[Pa\] | `framework_shear_modulus` is the post-melt one. |
| `calc_bulk_viscosity(solid_bulk_viscosity, postmelt_shear_viscosity, melt_fraction)` | `float` \[Pa s\] | |
| `model_name`, `parameters`, `get_parameter(name)`, `get_parameter_info()`, `with_parameters(**changes)`, `get_config_dict()`, `save_config(path)`, `save_binary(path)`, `load_binary(path)` | | As for every law (see [Equation-of-State and Shear-Modulus Laws](../Material/material_eos.md#python-api)). |

Without `melt_fraction`, `calc_weakening` uses the linear $\phi$ from the temperature and the curves, which is NaN at the melting temperature when the solidus equals the liquidus. A material treats that step itself (solid at the melting temperature, liquid above), so pass `melt_fraction` to reproduce it. The weakening and mixing laws take floats; to sweep a melting range, use a material's vectorized `calc_state` (see [Phases and Materials](../Material/materials.md#vectorized-evaluation)).

`get_config_dict()` gives `model` plus every parameter by config key, ready for the factory. A material saves its laws with itself, in its `melting` table and binary record.

`melting_curve_model_names()` and `melt_weakening_model_names()` list the canonical names; `TidalPy.PartialMelt.melting` also has `canonical_melting_curve_name`, `canonical_melt_weakening_name`, `bulk_modulus_mixing_model_names`, and `bulk_viscosity_mixing_model_names`.

## Limits and Failure Modes

- Tension ($P < P_\mathrm{ref}$), which trial central pressures of a structure solve can reach, holds $T_0$.
- A falling Simon and Glatzel curve without a `maximum_pressure_pa`, or with one at or past $P_\mathrm{ref} - a$, raises `ValueError`, as does an `a` of 0. Past `maximum_pressure_pa` a curve holds its value there, but a material is still meant for the pressure range its curves were fitted over.
- A non-finite pressure gives a NaN melting temperature and slope, and a material then stays solid.
- A material refuses a liquidus below its solidus at zero pressure, where every layer without `use_pressure_melting` reads them (`ValueError`). Deeper, where the liquidus falls to or below the solidus, it melts as a step at the solidus.
- A NaN solid value (a solid phase with no viscosity law) stays NaN through a weakening law rather than taking the liquid's.

## C++ API

Header-only, namespace `tidalpy`: `TidalPy/PartialMelt/melting_curve_.hpp`, `melt_weakening_.hpp`, and `melt_mixing_.hpp`. Classes carry a `c_` prefix (`c_SimonGlatzelCurve`, `c_HenningMeltWeakening`, ...).

- `c_MeltingCurveBase`: `calc_melting_temperature(pressure)`, `calc_melting_slope(pressure)`, and their `_vectorize` forms. `c_simon_glatzel(pressure, temperature, simon_a, simon_c, reference_pressure)` and `c_simon_glatzel_slope(...)` evaluate one branch.
- `c_MeltWeakeningBase::calc_weakening(const c_MeltWeakeningInputs&)`: the Python arguments as one struct, returning a `c_MeltWeakeningResult` (`shear_modulus`, `viscosity`).
- `c_BulkModulusMixingBase::calc_bulk_modulus(...)` and `c_BulkViscosityMixingBase::calc_bulk_viscosity(...)`, with the Python arguments.
- Per family (`melting_curve`, `melt_weakening`, `bulk_modulus_mixing`, `bulk_viscosity_mixing`): `c_find_<family>(name, params)` (throws `std::invalid_argument` for an unknown name or parameter), `c_<family>_from_binary(stream, force)`, `c_<family>_canonical_name(name)`, and, for the first two, `c_<family>_model_names()`.

### Adding a New Model

A new law is a `c_SpecModel<c_FooMeltWeakening, c_MeltWeakeningBase>` (or the family's base) with a `parameter_specs()` table, `C_CLASS_ID`, and its law (`p_calc_partial` for a weakening law; the base handles the range ends and the liquid floor). Then add a `BinaryClassID` (melting curves 71X, weakening 72X, bulk-modulus mixing 73X, bulk-viscosity mixing 74X), a registry row, a Cython subclass exported from `TidalPy/PartialMelt/__init__.py`, tests in `Tests/Test_PartialMelt/`, and a section here.

## References

- Andrault, D., Bolfan-Casanova, N., Lo Nigro, G., Bouhifd, M. A., Garbarino, G., and Mezouar, M. (2011). Solidus and liquidus profiles of chondritic mantle: Implication for melting of the Earth across its history. *Earth and Planetary Science Letters*, 304(1-2), 251-259.
- Fiquet, G., Auzende, A. L., Siebert, J., Corgne, A., Bureau, H., Ozawa, H., and Garbarino, G. (2010). Melting of peridotite to 140 gigapascals. *Science*, 329(5998), 1516-1518.
- Fischer, H.-J., and Spohn, T. (1990). Thermal-orbital histories of viscoelastic models of Io. *Icarus*, 83(1), 39-65.
- Hashin, Z., and Shtrikman, S. (1963). A variational approach to the theory of the elastic behaviour of multiphase materials. *Journal of the Mechanics and Physics of Solids*, 11(2), 127-140.
- Henning, W. G., O'Connell, R. J., and Sasselov, D. D. (2009). Tidally heated terrestrial exoplanets: Viscoelastic response models. *The Astrophysical Journal*, 707(2), 1000-1015.
- Kervazo, M., Tobie, G., Choblet, G., Dumoulin, C., and Běhounková, M. (2021). Solid tides in Io's partially molten interior. *Astronomy and Astrophysics*, 650, A72.
- Mavko, G. M. (1980). Velocity and attenuation in partially molten rocks. *Journal of Geophysical Research*, 85(B10), 5173-5189.
- McKenzie, D. (1984). The generation and compaction of partially molten rock. *Journal of Petrology*, 25(3), 713-765.
- Monteux, J., Andrault, D., and Samuel, H. (2016). On the cooling of a deep terrestrial magma ocean. *Earth and Planetary Science Letters*, 448, 140-149.
- Renaud, J. P., and Henning, W. G. (2018). Increased tidal dissipation using advanced rheological models: Implications for Io and tidally active exoplanets. *The Astrophysical Journal*, 857(2), 98.
- Simon, F., and Glatzel, G. (1929). Bemerkungen zur Schmelzdruckkurve. *Zeitschrift für anorganische und allgemeine Chemie*, 178(1), 309-316.
- Takei, Y. (2002). Effect of pore geometry on VP/VS: From equilibrium geometry to crack. *Journal of Geophysical Research*, 107(B2), 2043.
- Takei, Y., and Holtzman, B. K. (2009). Viscous constitutive relations of solid-liquid composites in terms of grain boundary contiguity. *Journal of Geophysical Research*, 114, B06205.
- Wood, A. B. (1955). *A Textbook of Sound*. G. Bell and Sons.
