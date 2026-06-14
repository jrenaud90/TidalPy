# Partial-Melt Models (`partial_melt_x`)

A **partial-melt** model maps a material's pre-melt (solid) viscosity and shear
modulus, plus its temperature, to the **post-melt** viscosity and shear modulus
(i.e. it applies melt weakening), and reports the volumetric **melt fraction**.

These quantities are **frequency-independent**, they depend only on the
temperature/pressure state fixed by the EOS solve, so in the whole-planet
love-number pipeline they are computed once and cached; only the downstream
[rheology](../rheology_x/rheology_models.md) (complex modulus) step is recomputed
per forcing frequency.

All quantities are **MKS**. The math mirrors the validated legacy implementation
in `TidalPy/rheology/partial_melt/melting_models.py`.

---

## Melt fraction

The volumetric melt fraction is a model-independent function of temperature and
the material's solidus/liquidus:

```
φ = clip( (T − T_solidus) / (T_liquidus − T_solidus), 0, 1 )
```

A degenerate envelope (`solidus >= liquidus`) yields `φ = 0` (fully solid).

---

## Models

| Model | Alias | Behaviour |
|-------|-------|-----------|
| `OffPartialMelt` | `none` | No weakening: post-melt strength equals pre-melt. |
| `SpohnPartialMelt` | `fischer`, `fischer_spohn` | Fischer & Spohn (1990): post-melt viscosity/shear depend only on temperature. |
| `HenningPartialMelt` | — | Henning (2009/2010): three-regime melt weakening. |

### Spohn (Fischer & Spohn 1990)

```
η_post = 10^( fs_visc_power_slope / T − fs_visc_power_phase )
G_post = 10^( fs_shear_power_slope / T − fs_shear_power_phase )
```
floored at the liquid viscosity / `liquid_shear`.

### Henning (2009/2010)

Three regimes in the melt fraction `φ`, with the transition band
`[crit_melt_frac, crit_melt_frac + crit_melt_frac_width]`:

- `φ ≤ 0` — pre-melt strengths returned unchanged.
- `0 < φ < crit` — sub-critical exponential weakening:
  `η = η_pre·exp(−hn_visc_slope_1·φ)`,
  `G = G_pre·exp(hn_shear_param_1/T − hn_shear_param_2)`.
- `crit ≤ φ ≤ crit+width` — breakdown band: the maximum sub-critical effect at the
  break-down temperature `T_break = solidus + crit·(liquidus − solidus)`, times a
  steep falloff `exp(−hn_*_falloff_slope·(φ − crit))`.
- `φ > crit+width` — liquid-like: `η = η_liquid`, `G = liquid_shear`.

All branches are floored at the liquid limits.

**References:** Fischer and Spohn (1990), *Icarus* 83, 39; Henning, O'Connell &
Sasselov (2009); Renaud & Henning (2018), *ApJ* 857, 98.

---

## Python API

```python
from TidalPy.partial_melt_x import (
    OffPartialMelt, SpohnPartialMelt, HenningPartialMelt, make_partial_melt)

melt = HenningPartialMelt(solidus_k=1600.0, liquidus_k=2000.0, liquid_shear_pa=1.0e-5)

phi = melt.calc_melt_fraction(1800.0)               # 0.5
phi, post_visc, post_shear = melt.calc_partial_melt(
    temperature_k=1700.0,
    premelt_viscosity=1.0e22,   # [Pa·s]
    premelt_shear=6.0e10,       # [Pa]
    liquid_viscosity=0.2,       # [Pa·s]
)

# Name factory (aliases, case-insensitive):
melt = make_partial_melt("fischer", {"solidus_k": 1500.0})
```

**Constructors / parameters**

- `OffPartialMelt(solidus_k=1600, liquidus_k=2000, liquid_shear_pa=1e-5)`
- `SpohnPartialMelt(..., fs_visc_power_slope=27000, fs_visc_power_phase=1.0,
  fs_shear_power_slope=82000, fs_shear_power_phase=40.6)`
- `HenningPartialMelt(..., crit_melt_frac=0.5, crit_melt_frac_width=0.05,
  hn_visc_slope_1=13.5, hn_visc_falloff_slope=370, hn_shear_param_1=40000,
  hn_shear_param_2=25, hn_shear_falloff_slope=700)`

**Methods / properties**

| Member | Returns | Description |
|--------|---------|-------------|
| `calc_melt_fraction(T)` | float | Melt fraction φ ∈ [0, 1]. |
| `calc_partial_melt(T, premelt_visc, premelt_shear, liquid_visc)` | `(φ, η_post, G_post)` | Post-melt viscosity [Pa·s] + shear modulus [Pa]. |
| `solidus`, `liquidus`, `liquid_shear` | float | Melt envelope. |
| `get_config_dict()` | dict | Model name + all parameters. |
| `save_config(path)` / `save_binary(path)` / `load_binary(path)` | — | Inherited from `TidalPyBaseClass`. |

`make_partial_melt(name, config=None)` returns the concrete subclass; absent
config keys fall back to the C++ defaults. Unknown names raise `ValueError`.

---

## C++ API

The C++ layer is canonical; the Cython classes are thin adapters.

**Config struct** `tidalpy::c_PartialMeltConfig` (defaults in parentheses):
`solidus_k` (1600), `liquidus_k` (2000), `liquid_shear_pa` (1e-5);
Spohn: `fs_visc_power_slope` (27000), `fs_visc_power_phase` (1.0),
`fs_shear_power_slope` (82000), `fs_shear_power_phase` (40.6);
Henning: `crit_melt_frac` (0.5), `crit_melt_frac_width` (0.05),
`hn_visc_slope_1` (13.5), `hn_visc_falloff_slope` (370),
`hn_shear_param_1` (40000), `hn_shear_param_2` (25),
`hn_shear_falloff_slope` (700).

**Input / result** `c_PartialMeltInputs { temperature_k, premelt_viscosity,
premelt_shear, liquid_viscosity }` → `c_PartialMeltResult { melt_fraction,
postmelt_viscosity, postmelt_shear_modulus }`.

**Base class** `c_PartialMeltBase : c_PhysicsBase` (in `partial_melt_base_.hpp`):
`calc_melt_fraction(double temperature_k) const`,
pure virtual `calc_partial_melt(const c_PartialMeltInputs&) const`,
`calc_partial_melt_vectorize(temperature, premelt_visc, premelt_shear,
liquid_viscosity, out_results)` (the primary radial sweep — one entry per slice),
and accessors `get_solidus/get_liquidus/get_liquid_shear`.

**Concrete models** (`partial_melt_.hpp`): `c_OffPartialMelt`,
`c_SpohnPartialMelt` (+ `get_visc_power_slope` etc.), `c_HenningPartialMelt`
(+ `get_crit_melt_frac` etc.). Binary class ids **701–703**.

**Enum factory:**
`c_partial_melt_model_from_name(const std::string&) -> c_PartialMeltModel`
(`Off`/`Spohn`/`Henning`, throws `std::invalid_argument` on unknown names);
`c_find_partial_melt(c_PartialMeltModel, const c_PartialMeltConfig&) ->
std::unique_ptr<c_PartialMeltBase>` (and a name overload);
`c_partial_melt_from_binary(std::istream&, bool force=false)` (peek class id →
build → `read_binary`).

---

## Adding a new model

1. Add a `c_<Name>PartialMelt : c_PartialMeltBase` in `partial_melt_.hpp`
   implementing `calc_partial_melt`, `write_binary`, `read_binary`.
2. Add its parameters to `c_PartialMeltConfig`.
3. Register a `c_PartialMeltModel::<Name>` enum value and wire it into
   `c_partial_melt_model_from_name`, `c_find_partial_melt`, and
   `c_partial_melt_from_binary`.
4. Reserve a unique `BinaryClassID` (next free in the 700-block) in
   `Utilities_x/binary_x/binary_.hpp`.
5. Add the Cython `cdef class`, factory branch, config-dict, and tests + this doc.
