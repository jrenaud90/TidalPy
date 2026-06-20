# Global (1D) Tidal Dissipation (`Tides_x.classes`)

The **global** (or "1D potential") approach computes a body's total tidal heating and the
three orbital potential derivatives (`dU/dM`, `dU/dw`, `dU/dO`) by summing over the active
tidal **modes** `(l, m, p, q)`. Each mode carries an orbital/spin frequency
`omega_lmpq = (l − 2p + q)·n − m·spin` and a precomputed potential weight; a **tide model**
supplies the per-mode dissipation multiplier `−Im[k_l(omega)]` that the
**collapse** multiplies in and sums.

The model-independent per-mode weights (eccentricity functions `G_lpq`, obliquity functions
`F_lmp`, and the common coefficient `G_lpq²·F_lmp²·(l−m)!/(l+m)!·(R/a)^(2l+1)·G·M_host/a`)
come from [`c_global_potential`](../Tides/index.md). This page documents the tide models and
the collapse that turn those weights into heating and torque.

The full complex Love-number suite **(k, h, l)** is always the transport type
([`c_LoveNumbers`](love/love_numbers.md)), even though only `k` drives heating and orbital
dynamics — so the displacement Love numbers from the radial solver are never discarded. The
analytic models cannot produce `h`, `l` (no radial solution) and return them as `NaN`.

All quantities are **MKS**; frequencies in rad s⁻¹. Supported degrees are `l = 2..10`.

**References:** Renaud et al. (2021, PSJ) — global dual-body dissipation (collapse form);
Efroimsky & Makarov (2013) — CPL/CTL frequency dependence.

---

## Models

| Model | Alias | Complex Love number `k_l(omega)` | `−Im[k_l]` |
|-------|-------|----------------------------------|------------|
| `RheologyTide` | `rheology` | supplied by the radial solver | `−Im[k_l]` from the solver |
| `FixedQTide` | `cpl`, `fixed_q` | `k_l·(1 − i/Q_l)` | `k_l/Q_l` (frequency independent) |
| `FixedLagTide` | `ctl`, `fixed_dt` | `k_l·(1 − i·omega·dt_l)` | `k_l·omega·dt_l` |
| `CTLQTide` | `ctl_q`, `fixed_dt_q` | `k_l·(1 − i·omega·dt_l/Q_l)` | `k_l·omega·dt_l/Q_l` |

Fixed per-degree parameters `k_l` (static Love number), `Q_l` (quality factor), and `dt_l`
(time lag [s]) are supplied as lists indexed from degree `l = 2` (index 0: `l = 2`). The
`rheology` model needs the radial solver and is driven by the world's `calc_tides` method,
not the standalone collapse below.

A zero/absent `Q_l` is treated as purely elastic (no dissipation) rather than a divide by
zero.

---

## Collapse

For a synchronously rotating, low-eccentricity body the `fixed_q` collapse reproduces the
classic CPL tidal-heating rate exactly:

```
E_dot = (21/2) · (k2/Q) · G · M_host² · R⁵ · n · e² / a⁶
```

The collapse sums over every active (nonzero-frequency) mode:

```
E_dot += E_dot_term_lmpq · (−Im[k_l])           # heating  [W]
dU/dX += dU_dX_term_lmpq · (−Im[k_l])           # X = M, w, O  [J kg⁻¹ rad⁻¹]
```

Layer-level heat partitioning (`tidal_scale`) is applied by the world afterward; the
whole-body collapse uses the unscaled `−Im[k]`.

---

## Python API

```python
from TidalPy.Tides_x.classes import (
    RheologyTide, FixedQTide, FixedLagTide, CTLQTide,
    make_tide, collapse_global_tides)
from TidalPy.constants import G

# Build a model directly or by name (aliases, case-insensitive):
tide = make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [50.0]})

love = tide.calc_love_numbers(degree_l=2, frequency=4.1e-5)   # LoveNumbers(k=0.3-0.006j, h=nan, l=nan)
neg_imk = tide.calc_neg_imk(degree_l=2, frequency=4.1e-5)     # 0.006
# The rheology model returns the supplied radial-solver suite unchanged:
#   make_tide("rheology").calc_love_numbers(2, w, solver_love) -> solver_love (k, h, l)

# Standalone global collapse for an analytic model:
result = collapse_global_tides(
    planet_radius=1.82e6,
    semi_major_axis=4.22e8,
    orbital_frequency=4.11e-5,
    spin_frequency=4.11e-5,   # synchronous
    obliquity=0.0,
    eccentricity=0.0041,
    host_mass=1.898e27,
    G_to_use=G,
    tide_model="cpl",
    tide_config={"fixed_k": [0.3], "fixed_q": [50.0]},
    max_degree_l=2,
    eccentricity_truncation=2)
# {"tidal_heating": W, "dUdM": ..., "dUdw": ..., "dUdO": ..., "num_modes": int}
```

**Constructors / parameters**

- `RheologyTide()`
- `FixedQTide(fixed_k=None, fixed_q=None)`
- `FixedLagTide(fixed_k=None, fixed_dt=None)`
- `CTLQTide(fixed_k=None, fixed_dt=None, fixed_q=None)`

where each `fixed_*` is a list indexed from `l = 2`.

**Methods / properties**

| Member | Returns | Description |
|--------|---------|-------------|
| `calc_love_numbers(degree_l, frequency, solver_love=None)` | `LoveNumbers` | `(k, h, l)`; analytic models set `h, l = NaN`. `solver_love` (a `LoveNumbers`) is returned unchanged by `RheologyTide`. |
| `calc_neg_imk(degree_l, frequency, solver_love=None)` | float | Dissipation multiplier `−Im[k_l]`. |
| `needs_radial_solve` | bool | `True` only for `RheologyTide`. |
| `get_fixed_k/q/dt(degree_l)` | float | Per-degree parameter (analytic models). |
| `get_config_dict()` | dict | Model name + per-degree parameters. |
| `save_binary(path)` / `load_binary(path)` | — | Inherited from `TidalPyBaseClass`. |

`make_tide(name, config=None)` returns the concrete subclass; unknown names raise
`ValueError`. `collapse_global_tides(...)` supports the analytic models only — the
`rheology` model raises `NotImplementedError` (use the world's `calc_tides`).

---

## C++ API

The C++ layer is canonical; the Cython classes are thin adapters.

**Config struct** `tidalpy::c_TideModelConfig` (in `tide_.hpp`): per-degree `std::vector<double>`
`fixed_k`, `fixed_q`, `fixed_dt` (indexed from `l = 2`; entries beyond `l = 10` ignored).

**Base class** `c_TideBase : c_PhysicsBase` (in `tide_base_.hpp`): pure virtual
`calc_love_numbers(int degree_l, double frequency, const c_LoveNumbers& solver_love) const`
(returns the full `c_LoveNumbers` suite) and `needs_radial_solve() const`; non-virtual
`calc_neg_imk(...)` = `−Im[calc_love_numbers(...).k]`.

**Models** (`tide_.hpp`): `c_RheologyTide`, `c_FixedQTide`, `c_FixedLagTide`,
`c_CTLQTide` (each with per-degree getters). Binary class ids **901–904**.

**Enum factory:** `c_tide_model_from_name(const std::string&) -> c_TideModel`
(`Rheology`/`FixedQ`/`FixedLag`/`CTLQ`, throws `std::invalid_argument` on unknown names);
`c_find_tide(c_TideModel, const c_TideModelConfig&) -> std::unique_ptr<c_TideBase>` (and a
name overload); `c_tide_from_binary(std::istream&, bool force=false)` (peek class id -> build
→ `read_binary`).

**Collapse** (`tide_collapse_.hpp`):
`c_collapse_global_tides(const c_GlobalPotentialStorage&, const c_TideBase&,
const c_IntMap<c_Key4, c_LoveNumbers>* solver_love_by_lmpq = nullptr) -> c_GlobalTideResult`.
The `solver_love_by_lmpq` map supplies the radial-solver Love numbers (k, h, l) per mode for
the `rheology` model; pass `nullptr` for the analytic models. `c_GlobalTideResult` holds
`tidal_heating`, `dU_dM`, `dU_dw`, `dU_dO`, `num_modes`, `error_code`.

---

## Adding a new model

1. Add a `c_<Name>Tide : c_TideBase` in `tide_.hpp` implementing `calc_love_numbers`
   (return a `c_LoveNumbers`; set `h, l = NaN` if no radial solution), `needs_radial_solve`,
   `write_binary`, `read_binary`.
2. Add any new parameters to `c_TideModelConfig`.
3. Register a `c_TideModel::<Name>` enum value and wire it into `c_tide_model_from_name`,
   `c_find_tide`, and `c_tide_from_binary`.
4. Reserve a unique `BinaryClassID` (next free in the 900-block) in
   `Utilities_x/binary_x/binary_.hpp`.
5. Add the Cython `cdef class`, factory branch, config-dict, and tests + this doc.
