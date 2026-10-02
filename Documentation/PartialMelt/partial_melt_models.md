# Melting Laws (`PartialMelt`)

_Updated: 2026-10-02_

The melting laws describe how a [material](../Material/materials.md) melts. Its two melting curves give the solidus and liquidus temperatures \[K\] at a pressure \[Pa\], and the melt fraction $\phi$ runs linearly between them. Its melt-weakening law gives the partially molten aggregate's shear modulus \[Pa\] and viscosity \[Pa s\] from the solid's and the liquid's values, and its optional bulk-mixing laws give the aggregate's bulk modulus \[Pa\] and bulk viscosity \[Pa s\]. Each law is a model of its own family, and the material combines them. The laws themselves know nothing of the phases' equations of state.

## Inheritance

```
c_TidalPyBaseClass
  └── c_PhysicsBase
        ├── c_MeltingCurveBase  (abstract)
        │     ├── c_ConstantMeltingCurve       aliases "constant", "const"
        │     ├── c_SimonGlatzelCurve          aliases "simon_glatzel", "simon-glatzel"
        │     ├── c_SimonGlatzel2Curve         aliases "simon_glatzel_2", "simon-glatzel-2"
        │     └── c_InterpolatedMeltingCurve   aliases "interpolate", "interp", "interpolated"
        ├── c_MeltWeakeningBase  (abstract)
        │     ├── c_NoMeltWeakening            aliases "none", "off"
        │     ├── c_SpohnMeltWeakening         aliases "spohn", "fischer", "fischer_spohn"
        │     └── c_HenningMeltWeakening       alias "henning"
        ├── c_BulkModulusMixingBase  (abstract)
        │     └── c_HashinShtrikmanMixing      aliases "hashin_shtrikman", "hs", "hashin-shtrikman"
        └── c_BulkViscosityMixingBase  (abstract)
              └── c_CompactionViscosity        aliases "compaction", "mckenzie"
```

Every concrete law derives from `c_SpecModel`, which gives it its parameters, config dict, and binary record from one table (see [C++ API](#c-api)). The Cython classes mirror the hierarchy: `MeltingCurveBase`, `ConstantMeltingCurve`, `SimonGlatzelCurve`, `SimonGlatzel2Curve`, `InterpolatedMeltingCurve`; `MeltWeakeningBase`, `NoMeltWeakening`, `SpohnMeltWeakening`, `HenningMeltWeakening`; `BulkModulusMixingBase`, `HashinShtrikmanMixing`; and `BulkViscosityMixingBase`, `CompactionViscosity`.

## Inputs and Result

A melting curve takes a pressure and returns a temperature, `calc_melting_temperature(pressure)`, and its slope $dT_m/dP$ \[K Pa$^{-1}$\], `calc_melting_slope(pressure)`.

A weakening law sees everything the material knows at a point, bundled in C++ into a `c_MeltWeakeningInputs` struct:

| Input | Units | Meaning |
|---|---|---|
| `temperature` | K | Local temperature. |
| `solidus`, `liquidus` | K | The melting curves at the local pressure. |
| `melt_fraction` | m$^3$ m$^{-3}$ | $\phi$. |
| `solid_shear`, `solid_viscosity` | Pa, Pa s | The solid phase's shear modulus and shear viscosity. |
| `liquid_shear`, `liquid_viscosity` | Pa, Pa s | The liquid phase's, which floor the result. |

It returns a `c_MeltWeakeningResult` of the aggregate's `shear_modulus` and `viscosity`. The bulk-modulus mixing law takes the solid's and the liquid's bulk moduli, the framework's (post-melt) shear modulus, and the melt fraction. The bulk-viscosity mixing law takes the solid's bulk viscosity, the post-melt shear viscosity, and the melt fraction.

## Melting Curves

The solidus and liquidus are two independent curves. Rock and iron melt at higher temperatures under pressure: peridotite's solidus, for example, rises from about 1660 K at the surface to about 4150 K at Earth's core-mantle boundary. A layer reads the curves at the local pressure when it sets `use_pressure_melting`, and at zero pressure otherwise. The melt fraction, the weakening laws (each anchored at the local solidus), the density mixing, and the bulk effects all use the curves at that pressure.

| Model | Melting temperature $T_m(P)$ \[K\] |
|---|---|
| `constant` | $T_0$ at every pressure |
| `simon_glatzel` | $T_0 \left(1 + (P - P_\mathrm{ref}) / a\right)^{1/c}$ |
| `simon_glatzel_2` | one Simon and Glatzel branch below a transition pressure $P_t$ and another above it |
| `interpolate` | linear in pressure between the points of a `pressure` and `temperature` table, held at the end values beyond it |

### Simon and Glatzel

The Simon and Glatzel (1929) law,

$$T_m(P) = T_0 \left(1 + \frac{P - P_\mathrm{ref}}{a}\right)^{1/c},$$

fits most planetary melting data. $T_0$ is the melting temperature at the reference pressure $P_\mathrm{ref}$ (usually zero), $a$ \[Pa\] sets the pressure at which the curve begins to rise, and $c$ sets how quickly the rise flattens. A negative $a$ gives a curve that falls with pressure, as ice Ih's does: MatPack's `ice_ih` uses $T_0$ = 273.16 K, $a$ = -415 MPa, and $c$ = 8.25. Below $P_\mathrm{ref}$ the curve holds $T_0$, and above `maximum_pressure_pa` (default: none), the end of the range the curve was fitted over, it holds its value there. A falling curve would reach 0 K at $P_\mathrm{ref} - a$, so it needs a `maximum_pressure_pa` below that; `ice_ih` ends at its triple point with ice III and liquid water, 208.566 MPa, where it holds about 251 K.

`simon_glatzel_2` joins two branches at a transition pressure $P_t$:

$$T_m(P) = \begin{cases} T_0 \left(1 + (P - P_\mathrm{ref}) / a\right)^{1/c} & P \le P_t \\ T_{0,\mathrm{high}} \left(1 + (P - P_\mathrm{ref,high}) / a_\mathrm{high}\right)^{1/c_\mathrm{high}} & P > P_t \end{cases}$$

Both branches are written in the absolute pressure, as the published fits are. The second branch lets a fit change slope across a phase transition, such as the one near 20 GPa at the top of Earth's lower mantle.

### Peridotite Fits

Monteux et al. (2016) fit two-branch Simon and Glatzel laws to the peridotite and chondritic-mantle melting experiments of Fiquet et al. (2010) and Andrault et al. (2011). MatPack's `peridotite` and `lower_mantle` use them, and they are the `simon_glatzel_2` defaults for the solidus.

| Curve | $T_0$ \[K\] | $a$ \[Pa\] | $c$ | $P_t$ \[Pa\] | $T_{0,\mathrm{high}}$ \[K\] | $a_\mathrm{high}$ \[Pa\] | $c_\mathrm{high}$ |
|---|---|---|---|---|---|---|---|
| Solidus | 1661.2 | 1.336e9 | 7.437 | 20.0e9 | 2081.8 | 1.0169e11 | 1.226 |
| Liquidus | 1982.1 | 6.594e9 | 5.374 | 20.0e9 | 78.74 | 4.054e6 | 2.44 |

These fits give the following melting temperatures.

| Pressure \[GPa\] | Solidus \[K\] | Liquidus \[K\] | Where |
|---|---|---|---|
| 0 | 1661 | 1982 | Surface |
| 5 | 2048 | 2202 | Base of Io's mantle |
| 20 | 2411 | 2569 | Branch transition |
| 60 | 3039 | 4030 | Mid lower mantle of the Earth |
| 135 | 4147 | 5619 | Earth's core-mantle boundary |

The experiments reach about 140 GPa, so the curves are extrapolated in the deeper mantles of planets larger than the Earth.

### Melting Slope

`calc_melting_slope(pressure)` returns $dT_m/dP$ \[K Pa$^{-1}$\]: for a Simon and Glatzel branch $T_0 / (a c) \, (1 + (P - P_\mathrm{ref}) / a)^{1/c - 1}$, and for a table the slope of the interval holding the pressure. It is zero for a constant curve and wherever a curve is held flat (below its reference pressure, past the end of a falling curve, or beyond a table's ends). A material uses the slopes for the latent heat's share of an adiabat's expansivity inside a melting range (see [Latent Heat](../Material/materials.md#latent-heat)).

### Behavior at the Limits

- Tension ($P < P_\mathrm{ref}$), which the structure solve's trial central pressures can reach, holds the reference temperature $T_0$.
- Past `maximum_pressure_pa` a Simon and Glatzel curve holds its value at that pressure, and a falling curve without one, or with one at or past $P_\mathrm{ref} - a$ where it would reach 0 K, raises `ValueError`. A material is still meant to be used inside the pressure range its curves were fitted over.
- A non-finite pressure gives a NaN melting temperature and slope, and a material then stays solid.
- The two curves are independent, but a material refuses a liquidus below its solidus at zero pressure, where every layer reads them without `use_pressure_melting` (`ValueError`). Deeper, where the liquidus falls to or below the solidus, the material melts as a step at the solidus.
- An `a` of 0 is refused when the curve is built.

### Choosing Melting Curves

We recommend pressure-dependent curves, and `use_pressure_melting` on, for any rocky layer whose base is more than a few GPa deep: the peridotite solidus at the base of Io's mantle (about 5 GPa) is already about 390 K above its zero-pressure value. The curves matter most together with the [compressed expansivity](../Material/material_eos.md#expansivity-under-compression), which sets how much a convecting mantle warms with depth. With a constant expansivity an Earth-like mantle's adiabat climbs above the peridotite solidus at depth. The bundled worlds with warm silicate mantles use the peridotite curves with pressure melting on. A constant curve is the usual choice for a thin layer, for a comparison with a published model that used one, and for a material without a measured curve.

## Melt Weakening

| Model | Behavior |
|---|---|
| `NoMeltWeakening` | The solid's shear modulus and viscosity until fully molten, then the liquid's. |
| `SpohnMeltWeakening` | Fischer and Spohn (1990). Above the solidus, both strengths fall with temperature from their values at the solidus. |
| `HenningMeltWeakening` | Henning et al. (2009) and Renaud and Henning (2018). Three regimes separated by a critical melt fraction. |

Every law returns the solid pair with no melt ($\phi \le 0$) and the liquid pair when fully molten ($\phi \ge 1$), and floors its result at the liquid's values in between. A NaN solid value (a solid phase with no viscosity law) stays NaN rather than taking the liquid's. A material without a weakening law behaves as `none`. The material, not the law, decides from the returned shear modulus where the radial solver treats the aggregate as a liquid. Both temperature laws are anchored at the solidus, so they carry over to materials whose solidus is not the 1600 K silicate value the published fits assume.

### Breakdown Band

The Spohn and Henning laws describe the partially molten framework up to a rheological transition, the breakdown band of melt fraction from `crit_melt_frac` ($\phi_c$) to `crit_melt_frac + crit_melt_frac_width` ($\phi_c + w$). Across the band the framework's pair blends into the liquid's, so the aggregate reaches the liquid's values at the band's end and holds them past it. With $s = (\phi - \phi_c) / w$, the framework's values $\eta_f$ and $\mu_f$, and the liquid's $\eta_l$ and $\mu_l$:

$$\eta = \eta_f^{\,1 - s} \, \eta_l^{\,s}, \qquad \mu = (1 - s)\,\mu_f + s\,\mu_l$$

The viscosity blends log-linearly, since it spans many orders of magnitude, and the shear modulus linearly, since a liquid's is usually zero. Both laws are therefore continuous in temperature from below the solidus to above the liquidus, which an integration through the melting range (a thermal evolution, say) needs: a step in viscosity at the band's end of about $10^{10}$ makes an implicit integrator crawl where the mantle sits on it. A zero `crit_melt_frac_width` makes the transition a step into the liquid at $\phi_c$.

### None

The solid's values pass through untouched inside the melting range, and the melt fraction is still reported. Use it for a material that melts without a strength change worth modeling, or to isolate the effect of the density and thermal terms.

### Spohn (Fischer and Spohn 1990)

$$\eta = 10^{\,L_\eta + s_\eta (1/T - 1/T_\mathrm{sol})}, \qquad \mu = 10^{\,L_\mu + s_\mu (1/T - 1/T_\mathrm{sol})}$$

below the [breakdown band](#breakdown-band), with the slopes $s$ given by `visc_power_slope` and `shear_power_slope`, and the base-10 logarithms of the strengths at the solidus, $L$, by `visc_log10_at_solidus` and `shear_log10_at_solidus`. Left unset (the default), each $L$ is the solid phase's own value at the solidus and the local pressure, so the law continues the solid's viscosity and shear modulus into the melting range without a step.

Fischer and Spohn (1990) fit silicates with absolute laws, $10^{27000/T - 1}$ Pa s and $10^{82000/T - 40.6}$ Pa. Those are the form above at a 1600 K solidus with $L_\eta$ = 15.875 and $L_\mu$ = 10.65, so setting those two reproduces the published fits there; the strengths then step at the solidus from the solid's values to the fit's. Anchoring the law at the material's own solidus keeps it usable for other materials: the absolute fit at an icy solidus of 273 K would give a shear modulus of $10^{260}$ Pa.

### Henning (2009, 2018)

Three regimes in the melt fraction, separated by the [breakdown band](#breakdown-band) from $\phi_c$ to $\phi_c + w$. Write $T_\mathrm{break} = T_\mathrm{sol} + \phi_c (T_\mathrm{liq} - T_\mathrm{sol})$ for the temperature at which $\phi_c$ is reached. The framework's pair is

| Regime | Viscosity $\eta_f$ | Shear modulus $\mu_f$ |
|---|---|---|
| $\phi \le 0$ | $\eta_s$ | $\mu_s$ |
| $0 < \phi < \phi_c$ | $\eta_s \exp(-a_\eta \phi)$ | $\mu_s \exp[b_1 (1/T - 1/T_\mathrm{sol})]$ |
| $\phi_c \le \phi < \phi_c + w$ | $\eta_s \exp(-a_\eta \phi_c) \exp(-f_\eta (\phi - \phi_c))$ | $\mu_s \exp[b_1 (1/T_\mathrm{break} - 1/T_\mathrm{sol})] \exp(-f_\mu (\phi - \phi_c))$ |

which the band then blends into the liquid's $\eta_l$ and $\mu_l$, held from $\phi_c + w$ on. Here $\eta_s$ and $\mu_s$ are the solid's values, $a_\eta$ is `visc_slope_1`, $b_1$ is `shear_param_1`, and $f_\eta$ and $f_\mu$ are `visc_falloff_slope` and `shear_falloff_slope`. Every branch is floored at the liquid's values. On its own the falloff reaches only about $10^{-11}$ of the solid's viscosity by the band's end (for the defaults), far above a melt's; the blend closes that gap.

The shear law is Henning et al. (2009) Eq. 20, $\exp(40000/T - 25)$, anchored at the solidus: their constant 25 is $40000/1600$, the silicate solidus they calibrated to, so a 1600 K solidus reproduces it exactly. Written this way the shear modulus is continuous at any solidus, where the published constant would stiffen an icy layer with a 273 K solidus by a factor of about $e^{121}$ ($e^{40000/273 - 25}$) as it began to melt.

Below the critical melt fraction, melt sits in isolated pockets and weakens the solid framework gradually. Above it the framework loses contact and the material behaves as a crystal-laden liquid, a drop of many orders of magnitude. The breakdown band is a steep but finite bridge between the two regimes. Its width is a numerical convenience that keeps the transition continuous, not a measured quantity.

### Choosing a Weakening Law

`henning` is the usual choice for a silicate mantle: it weakens the framework gradually until the critical melt fraction and then collapses it. `spohn` reproduces the earlier Io models built on Fischer and Spohn's fits. `none` suits ices and other materials whose partially molten strength is not constrained, and a material that melts as a step, where no partially molten state exists.

## Bulk Mixing

Melt lowers the bulk modulus far less than the shear modulus. A silicate melt's bulk modulus is of the same order as the rock's (about 20 GPa at low pressure against about 130 GPa) while its shear modulus vanishes (Mavko 1980; Takei 2002). Without a bulk modulus mixing law, a material's bulk modulus blends linearly from the solid's into the liquid's across its weakening law's breakdown band (`calc_band_blend`), the same band over which the shear modulus reaches the liquid's, so the mush the radial solver treats as a liquid has the liquid's bulk modulus too; with no weakening law it steps at full melt, with the shear modulus. Without a bulk viscosity mixing law, the bulk viscosity is the solid's until fully molten, then the liquid's.

### Hashin-Shtrikman

The post-melt bulk modulus is the Hashin and Shtrikman (1963) bound for melt of bulk modulus $K_l$ in a solid framework of bulk modulus $K_s$, evaluated with the framework's post-melt shear modulus $\mu$:

$$K = K_s + \frac{\phi}{\dfrac{1}{K_l - K_s} + \dfrac{1 - \phi}{K_s + \tfrac{4}{3}\mu}}$$

The material applies it to both the isothermal and the adiabatic moduli, each phase's own. While the framework holds, $\mu$ is close to the solid's and this is the upper bound for isolated melt pockets, a weak reduction (about 16 percent at $\phi = 0.1$ for $K_s$ = 130, $\mu$ = 60, $K_l$ = 20 GPa). Once the weakening law has collapsed the framework's shear modulus (Henning past the critical melt fraction), the same expression becomes the Reuss (Wood 1955) average of a crystal suspension, and it reaches $K_l$ at $\phi = 1$.

This is the unrelaxed (undrained) modulus: the melt is held in place over a forcing cycle. Isolated pockets are the stiffest geometry, and melt that wets grain edges or forms films weakens the framework more (Mavko 1980; Takei 2002), so the bound is an upper limit. Its relaxation as melt moves is a bulk rheology's job, at the rate the bulk viscosity sets.

### Compaction Viscosity

A melt-free rock has no viscous compaction, so its bulk response is elastic; melt adds one as it moves through the matrix. The compaction law adds a matrix bulk viscosity in series with the solid's,

$$\frac{1}{\zeta} = \frac{1}{\zeta_s} + \frac{\phi^n}{c\,\eta},$$

with $\eta$ the post-melt shear viscosity, $c$ `coefficient`, and $n$ `exponent`. $n = 1$ with $c$ of order 1 is the classic compaction viscosity $\eta/\phi$ (McKenzie 1984). Micromechanical models give a bulk viscosity of the same order as the shear viscosity instead (Takei and Holtzman 2009), which is $n = 0$. The series form is continuous at the solidus for $n > 0$. A non-finite or non-positive solid bulk viscosity counts as no solid dashpot, so melt alone sets $\zeta$.

The bulk viscosity reaches the tides only through the layer's bulk rheology, which is elastic unless the layer or its material names one. A Maxwell bulk rheology would let the bulk modulus relax to zero at long periods, which is unphysical. The [Zener](../Rheology/rheology_models.md) (standard linear solid) rheology relaxes it to a set fraction $r$ of the unrelaxed modulus instead, the drained-to-undrained ratio. For isolated pockets at 10% melt that ratio is near 0.9, and melt films lower it further. Bulk dissipation in a partially molten layer, while understudied, has been shown to rival the shear dissipation (Kervazo et al. 2021).

## Parameters

Each parameter carries two names: the constructor keyword, which also reads as an attribute, and the config key used in a TOML table, a factory config dict, and `get_config_dict()`. A dimensional config key ends in its unit while the code name does not. Either name is accepted wherever a parameter is given.

**Melting curves**

| Parameter | Config key | Default | Units | Used by |
|---|---|---|---|---|
| `temperature` | `temperature_k` | 1600.0 (`simon_glatzel_2` 1661.2; `[1600.0]` for interpolated) | K | All |
| `simon_a` | `simon_a_pa` | 1.0e9 (`simon_glatzel_2` 1.336e9) | Pa | Simon and Glatzel |
| `simon_c` | `simon_c` | 5.0 (`simon_glatzel_2` 7.437) | - | Simon and Glatzel |
| `reference_pressure` | `reference_pressure_pa` | 0.0 | Pa | Simon and Glatzel |
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

The Hashin-Shtrikman law has no parameters. The liquid's shear modulus and viscosity, which every weakening law falls to, belong to the material's liquid phase (a phase with no shear-modulus law has a shear modulus of 0).

## Python API

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

Constructors take every parameter their law uses, by its argument name or config key, as keywords or positionally in the order `get_parameter_info()` lists them, with the defaults from the tables above. `make_melting_curve`, `make_melt_weakening`, `make_bulk_modulus_mixing`, and `make_bulk_viscosity_mixing` each take `(model_name, config=None)`, resolve a name or alias case-insensitively, and build the law from `config`. Absent keys take the law's defaults. An unknown name, a key the law does not read, or a value outside a parameter's bounds raises `ValueError` naming the closest accepted name or key.

| Member | Returns | Description |
|---|---|---|
| `calc_melting_temperature(pressure)` | `float` or `np.ndarray` \[K\] | The melting temperature (melting curves). |
| `calc_melting_slope(pressure)` | `float` or `np.ndarray` \[K Pa$^{-1}$\] | $dT_m/dP$ (melting curves). |
| `calc_weakening(temperature, solidus, liquidus, solid_shear, solid_viscosity, liquid_shear, liquid_viscosity, melt_fraction=None)` | `(shear_modulus, viscosity)` | The aggregate's strengths \[Pa, Pa s\] (weakening laws). Without `melt_fraction` it is the linear $\phi$ from the temperature and the curves, which is NaN at the melting temperature when the solidus equals the liquidus. A material treats that step itself (solid at the melting temperature, liquid above it), so pass `melt_fraction` to reproduce it. |
| `calc_bulk_modulus(solid_bulk_modulus, liquid_bulk_modulus, framework_shear_modulus, melt_fraction)` | `float` \[Pa\] | The aggregate's bulk modulus (bulk-modulus mixing). |
| `calc_bulk_viscosity(solid_bulk_viscosity, postmelt_shear_viscosity, melt_fraction)` | `float` \[Pa s\] | The aggregate's bulk viscosity (bulk-viscosity mixing). |
| `model_name`, `parameters`, `get_parameter(name)`, `get_parameter_info()`, `with_parameters(**changes)`, `get_config_dict()`, `save_config(path)` | | As for every law (see [Equation-of-State and Shear-Modulus Laws](../Material/material_eos.md#python-api)). |

`melting_curve_model_names()` and `melt_weakening_model_names()` list the canonical names, and `TidalPy.PartialMelt.melting` also holds `canonical_melting_curve_name`, `canonical_melt_weakening_name`, `bulk_modulus_mixing_model_names`, and `bulk_viscosity_mixing_model_names`.

### Factory Internals

At the C++ level each family has one registry (`c_melting_curve_registry()`, `c_melt_weakening_registry()`, `c_bulk_modulus_mixing_registry()`, `c_bulk_viscosity_mixing_registry()`) listing each law's names (canonical first, then aliases), its binary class id, and its constructor. `c_find_melting_curve(name, params)` and its siblings build a law from a `c_ParamMap` and throw `std::invalid_argument` for an unknown name or parameter. The Python factories pass the config dict through to them and wrap the result in the matching class.

### Vectorized Evaluation

`calc_melting_temperature` and `calc_melting_slope` take a float or an array and return the same shape. The weakening and mixing laws take floats; a material evaluates them point by point inside its own vectorized `calc_state` (see [Phases and Materials](../Material/materials.md#vectorized-evaluation)), which is the usual way to sweep a melting range.

### Attaching Melting Laws to a `Material`

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

A material holding laws shares them rather than copying them. The declarative form is the `melting` table of a material, which a MatPack file or a world TOML holds:

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

## C++ API

The laws are in `TidalPy/PartialMelt/melting_curve_.hpp`, `melt_weakening_.hpp`, and `melt_mixing_.hpp` (namespace `tidalpy`, header only).

- `c_MeltingCurveBase` declares `calc_melting_temperature(pressure)` and `calc_melting_slope(pressure)` pure virtual and provides `calc_melting_temperature_vectorize` and `calc_melting_slope_vectorize`. The free functions `c_simon_glatzel(pressure, temperature, simon_a, simon_c, reference_pressure)` and `c_simon_glatzel_slope(...)` evaluate one branch.
- `c_MeltWeakeningBase::calc_weakening(const c_MeltWeakeningInputs&)` handles the ends of the range and the liquid floor and calls the protected `p_calc_partial(inputs, out)` inside it, the one method a weakening law implements.
- `c_BulkModulusMixingBase::calc_bulk_modulus(solid_bulk_modulus, liquid_bulk_modulus, framework_shear_modulus, melt_fraction)` and `c_BulkViscosityMixingBase::calc_bulk_viscosity(solid_bulk_viscosity, postmelt_shear_viscosity, melt_fraction)` are pure virtual.

Each family has `c_find_<family>(name, params)`, `c_<family>_from_binary(stream, force)` (used when a material is loaded), and `c_<family>_canonical_name(name)`, where `<family>` is `melting_curve`, `melt_weakening`, `bulk_modulus_mixing`, or `bulk_viscosity_mixing`; the first two also have `c_<family>_model_names()`. `c_Material` (`Material/material_.hpp`) calls the laws from its `calc_state`.

## Serialization

| Call | Result |
|---|---|
| `get_config_dict()` | `model` plus every parameter under its config key, ready for the factory. |
| `save_config(path)` | That dict written as TOML. |
| `save_binary(path)` / `load_binary(path)` | The law's TidalPy binary record, its parameters written by key. |

A material saves and restores its melting laws with itself, under its `melting` table in the config dict and as optional records in its binary record.

## Adding a New Model

To add a weakening law named `Foo` (a melting curve or a mixing law is the same in its own header):

1. Add `c_FooMeltWeakening : public c_SpecModel<c_FooMeltWeakening, c_MeltWeakeningBase>` to `melt_weakening_.hpp`: its `parameter_specs()` table, `C_CLASS_ID`, two constructors that call `p_initialize`, and `p_calc_partial`. Override `p_validate` for checks across parameters and `p_update_derived` for cached values.
2. Reserve a unique `BinaryClassID` in `Utilities/binary/binary_.hpp`: melting curves occupy the 71X range, weakening laws 72X, bulk-modulus mixing 73X, and bulk-viscosity mixing 74X.
3. Add one row to `c_melt_weakening_registry()`.
4. Add a two-line Cython subclass (a docstring and `MODEL_NAME`) to `melting.pyx`, include it in the family's `ModelFamily` list, export it from `TidalPy/PartialMelt/__init__.py`, add its physics tests to `Tests/Test_PartialMelt/`, and document it here. The generic tests in `Tests/Test_Utilities/Test_Classes/test_spec_models_01.py` cover its parameters, config, binary record, and errors without changes.

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
