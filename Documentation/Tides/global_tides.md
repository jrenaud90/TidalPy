# Global (1D) Tidal Dissipation (`Tides.classes`)

_Updated: 2026-10-02_

The global (or "1D potential") method computes a body's total tidal heating and the three orbital potential derivatives (`dU/dM`, `dU/dw`, `dU/dO`) by summing over the active tidal modes `(l, m, p, q)`, degrees `l = 2..10`. Each mode has a model-independent weight, and a tide model supplies its dissipation, $-\mathrm{Im}[k_{l}(\omega)]$. For depth-resolved heating, see [3D Tidal Heating](multilayer_3d_heating.md). Every model returns the full complex Love-number suite ($k$, $h$, $l$; see [Love numbers](love/love_numbers.md)), though only $k$ drives heating and dynamics; the analytic models return $h$ and $l$ as `NaN`.

## Example

```python
from TidalPy.Tides.classes import (
    RheologyTide, FixedQTide, FixedLagTide, CTLQTide,
    make_tide, collapse_global_tides)

# Build a model directly or by name (aliases, case-insensitive):
tide = make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [50.0]})
same_tide = make_tide(
    "cpl",
    fixed_k=[0.3],
    fixed_q=[50.0])                                              # Parameters as keywords
print(tide)                                                      # FixedQTide('fixed_q', fixed_k=[0.3], fixed_q=[50])

love    = tide.calc_love_numbers(degree_l=2, frequency=4.1e-5)   # LoveNumbers(k=0.3-0.006j, h=nan, l=nan)
neg_imk = tide.calc_neg_imk(degree_l=2, frequency=4.1e-5)        # 0.006
# The rheology model returns the supplied radial-solver suite unchanged:
#   make_tide("rheology").calc_love_numbers(2, w, solver_love) -> solver_love (k, h, l)

# Standalone global collapse for an analytic model
# (orbital-state arguments in the order of the world's calc_tides)
result = collapse_global_tides(
    planet_radius=1.82e6,
    orbital_frequency=4.11e-5,
    spin_frequency=4.11e-5,   # synchronous
    eccentricity=0.0041,
    obliquity=0.0,
    semi_major_axis=4.22e8,
    host_mass=1.898e27,
    G_to_use=None,            # None takes the configured G (SciPy's)
    tide_model="cpl",
    tide_config={"fixed_k": [0.3], "fixed_q": [50.0]},
    max_degree_l=2,
    eccentricity_truncation=2)
# {"tidal_heating": W, "dUdM": ..., "dUdw": ..., "dUdO": ..., "num_modes": int}
```

`collapse_global_tides` supports the analytic models only; the `rheology` model raises `NotImplementedError` (use the world's `calc_tides`). On a world, attach a model with `set_tide_model` and run `calc_tides`; see [Worlds](../Structures/worlds/worlds.md).

## Models

| Model | Name (aliases) | Complex Love number $k_{l}(\omega)$ | $-\mathrm{Im}[k_{l}]$ |
|-------|-------|----------------------------------|------------|
| `RheologyTide` | `rheology` | supplied by the radial solver | $-\mathrm{Im}[k_{l}]$ from the solver |
| `FixedQTide` | `fixed_q` (`cpl`, `constant_phase_lag`) | $k_{l}\,(1 - i/Q_{l})$ | $k_{l}/Q_{l}$ (frequency independent) |
| `FixedLagTide` | `fixed_dt` (`ctl`, `constant_time_lag`) | $k_{l}\,(1 - i\,\lvert\omega\rvert\,\Delta t_{l})$ | $k_{l}\,\lvert\omega\rvert\,\Delta t_{l}$ |
| `CTLQTide` | `fixed_dt_q` (`ctl_q`, `constant_time_lag_and_q`) | $k_{l}\,(1 - i\,\lvert\omega\rvert\,\Delta t_{l}/Q_{l})$ | $k_{l}\,\lvert\omega\rvert\,\Delta t_{l}/Q_{l}$ |

The first name is canonical (`model_name`, config dicts). Each model depends on the frequency's magnitude only; the mode's sign is in the collapse coefficients. The analytic models take per-degree lists, index 0 being `l = 2`, of at most nine values; a degree past the end of a list, or a zero $Q_{l}$, has no dissipation. The `rheology` model takes no parameters and is driven by the world's `calc_tides`.

### Parameters

Constructors take parameters by argument name or config key, as keywords, positionally (in `get_parameter_info()` order), or as a `config=` table. `make_tide(model_name, config=None, **parameters)` merges keyword parameters over `config`. A list left out (or `None`) takes its configured `[tides]` value, so `FixedQTide([0.3])` has the configured Q.

| Parameter | Config key | Default | Bounds | Models |
|---|---|---|---|---|
| `fixed_k` | `fixed_k` | `[tides] fixed_k` | non-negative, at most 9 values | `fixed_q`, `fixed_dt`, `fixed_dt_q`: static potential Love number $k_l$ of each degree. |
| `fixed_q` | `fixed_q` | `[tides] fixed_q` | non-negative, at most 9 values | `fixed_q`, `fixed_dt_q`: quality factor $Q_l$ of each degree; 0 is no dissipation. |
| `fixed_dt` | `fixed_dt_s` | `[tides] fixed_dt_s` | non-negative, at most 9 values | `fixed_dt`, `fixed_dt_q`: time lag $\Delta t_l$ of each degree \[s\]. |

In `CTLQTide` the positional order is `fixed_k`, `fixed_dt`, `fixed_q`. A key the model does not read raises `ValueError` naming the closest one it does (`make_tide("fixed_q", {"fixed_dt_s": [600.0]})` raises), as do an unknown name and an out-of-bounds value; `collapse_global_tides` behaves the same. `tide_model_names()` lists the canonical names and `tide_config_keys(name)` a model's keys. A world built from TOML passes its tide model only the `[tides]` lists it reads.

## Physics

### Tidal Potential

A host of mass $M_h$ on an orbit of semi-major axis $a$, eccentricity $e$, and obliquity $I$ relative to the body's equator raises, at the surface of a body of radius $R$, the tidal potential (Kaula 1964; Efroimsky and Williams 2009, Eq. 18)

$$U(\theta, \phi, t) = \frac{G M_h}{a}\sum_{l=2}^{\infty}\left(\frac{R}{a}\right)^{l}\sum_{m=0}^{l}\frac{(l-m)!}{(l+m)!}\left(2-\delta_{0m}\right)P_{lm}(\cos\theta)\sum_{p=0}^{l}F_{lmp}(I)\sum_{q=-\infty}^{\infty}G_{lpq}(e)\,\mathcal{T}_{lm}\!\left(\omega_{lmpq}t - m\phi\right),$$

where $\theta$ is the colatitude, $\phi$ the east longitude, $P_{lm}$ the associated Legendre functions without the Condon-Shortley phase, $\mathcal{T}_{lm}$ is $\cos$ for even $l - m$ and $\sin$ for odd $l - m$, and $F_{lmp}$ and $G_{lpq}$ are the [obliquity](Obliquity.md) and [eccentricity](Eccentricity.md) functions. Each $(l, m, p, q)$ is a tidal mode with the forcing frequency

$$\omega_{lmpq} = (l - 2p + q)\,n - m\,\dot{\theta},$$

where $n$ is the orbital mean motion and $\dot{\theta}$ the spin rate \[rad s$^{-1}$\]. Periapse and node precession are neglected. The degree range and the two truncation levels decide which modes are kept.

### Dissipation and Orbital Derivatives

The body couples each mode with the complex Love number $k_l$ at the forcing frequency $\chi_{lmpq} = |\omega_{lmpq}|$, and the tide model supplies $K_{l} = -\mathrm{Im}[k_{l}(\chi_{lmpq})]$. Averaged over the orbit and over apsidal precession, the tidal heating $\dot{E}$ \[W\] and the derivatives of the tidal potential with respect to the mean anomaly $\mathcal{M}$, the argument of periapse $\varpi$, and the node $\Omega$ \[J kg$^{-1}$ rad$^{-1}$\] are (Renaud et al. 2021, Eq. 7)

$$\begin{bmatrix} \partial U / \partial \mathcal{M} \\ \partial U / \partial \varpi \\ \partial U / \partial \Omega \\ \dot{E} \end{bmatrix} = \frac{G M_h}{a}\sum_{l=2}^{l_{\max}}\left(\frac{R}{a}\right)^{2l+1}\sum_{m=0}^{l}\frac{(l-m)!}{(l+m)!}\left(2-\delta_{0m}\right)\sum_{p=0}^{l}F_{lmp}^{2}(I)\sum_{q}G_{lpq}^{2}(e)\begin{bmatrix} (l-2p+q)\,\mathrm{sgn}(\omega_{lmpq})\,K_{l} \\ (l-2p)\,\mathrm{sgn}(\omega_{lmpq})\,K_{l} \\ m\,\mathrm{sgn}(\omega_{lmpq})\,K_{l} \\ \chi_{lmpq}\,M_h\,K_{l} \end{bmatrix}.$$

Zero-frequency modes do not dissipate and are skipped. The derivatives set the orbital and spin rates ([Dynamics](../Dynamics/dynamics.md)). `collapse_global_tides` and a world's `calc_tides` return this sum; a layered world then splits the heating among its layers by their `tidal_scale`.

For a synchronously rotating body at zero obliquity, only the $q = \pm 1$ modes of degree 2 dissipate to leading order in $e$, and with $K_2 = k_2/Q$ the sum reduces to the constant-phase-lag heating (Peale and Cassen 1978)

$$\dot{E} = \frac{21}{2}\,\frac{k_{2}}{Q}\,\frac{G M_h^{2} R^{5} n\, e^{2}}{a^{6}},$$

which eccentricity truncation level 2 at degree 2 reproduces.

## Methods and Properties

| Member | Returns | Description |
|--------|---------|-------------|
| `calc_love_numbers(degree_l, frequency, solver_love=None)` | `LoveNumbers` | `(k, h, l)`; analytic models set `h, l = NaN`. `solver_love` (a `LoveNumbers`) is returned unchanged by `RheologyTide`. |
| `calc_neg_imk(degree_l, frequency, solver_love=None)` | float | Dissipation multiplier `−Im[k_l]`. |
| `needs_radial_solve` | bool | `True` only for `RheologyTide`. |
| `get_fixed_k(degree_l)`, `get_fixed_q(degree_l)`, `get_fixed_dt(degree_l)` | float | That degree's value, 0 past the end of the list; NaN for a model without the parameter, which tells a world's `cpl` or `ctl` [Love method](love/love_numbers.md) whether it can fall back to the tide model. |
| `model_name` | str | The canonical name, for example `fixed_q`. |
| `get_config_dict()` | dict | `model` plus the per-degree lists as given (not padded). |
| `save_config(path)`, `save_binary(path)`, `load_binary(path)`, `get_schema_version_str()` | - | Shared by every physics model; see [Base Classes](../Utilities/classes.md). |

Parameters also read as attributes (`tide.fixed_k`), with `parameters`, `get_parameter(name)`, `get_parameter_info()`, and `with_parameters(**changes)`.

## C++ API

```cpp
#include "tide_.hpp"

using namespace tidalpy;

c_ParamMap params;  // Config key to values; a per-degree list from l = 2
params["fixed_k"] = {0.3};
params["fixed_q"] = {50.0};

std::unique_ptr<c_TideBase> tide = c_find_tide("cpl", params);
const double neg_imk = tide->calc_neg_imk(2, 4.1e-5, c_LoveNumbers());  // 0.006
```

- `c_TideBase` (`tide_base_.hpp`): `calc_love_numbers(int degree_l, double frequency, const c_LoveNumbers& solver_love) const`, `needs_radial_solve() const`, `calc_neg_imk(...)`, and `get_fixed_k`, `get_fixed_q`, `get_fixed_dt` (NaN unless overridden).
- Models (`tide_.hpp`): `c_RheologyTide`, `c_FixedQTide`, `c_FixedLagTide`, and `c_CTLQTide`.
- `c_find_tide(const std::string& name, const c_ParamMap& params)` returns a `std::unique_ptr<c_TideBase>` and throws `std::invalid_argument` for an unknown name or parameter or an out-of-bounds value. It does not read the TidalPy configuration: a list it is not given is empty. `c_tide_from_binary(std::istream&, bool force=false)` rebuilds a model from its binary record; `c_tide_canonical_name(name)` and `c_tide_model_names()` resolve names.
- `c_collapse_global_tides(const c_GlobalPotentialStorage&, const c_TideBase&, const c_IntMap<c_Key4, c_LoveNumbers>* solver_love_by_lmpq = nullptr)` (`tide_collapse_.hpp`) returns a `c_GlobalTideResult` (`tidal_heating`, `dU_dM`, `dU_dw`, `dU_dO`, `num_modes`, `error_code`). Pass the radial-solver Love numbers per mode for `rheology`, `nullptr` otherwise.

## Adding a New Model

1. In `tide_.hpp`, add `c_<Name>Tide` deriving from `c_AnalyticTide<c_<Name>Tide>` (analytic) or `c_SpecModel<c_<Name>Tide, c_TideBase>`, with its `parameter_specs()` table, `C_CLASS_ID`, two constructors that call `p_initialize`, and `calc_love_numbers` (`h, l = NaN` with no radial solution). A non-analytic model also implements `needs_radial_solve`; override `get_fixed_k`, `get_fixed_q`, or `get_fixed_dt` for values it carries.
2. Add a row to `c_tide_registry()` with its names and aliases.
3. Reserve a `BinaryClassID` (next free in the 900-block) in `Utilities/binary/binary_.hpp`.
4. Add a Cython `cdef class` (docstring and `MODEL_NAME`) to `tide.pyx` and `tide.pxd`, add it to the `ModelFamily` list, export it from `Tides/classes/__init__.py`, add tests to `Tests/Test_Tides/Test_Classes/` (the generic `test_spec_models_01.py` covers parameters, config, binary, and errors), add any `[tides]` default keys to the TidalPy configuration, and document it here.

## References

- Renaud, J. P., et al. (2021). Tidal dissipation in dual-body, highly eccentric, and nonsynchronously rotating systems: Applications to Pluto-Charon and the exoplanet TRAPPIST-1e. *The Planetary Science Journal*, 2(1), 4. The collapse form of global dual-body dissipation.
- Efroimsky, M., and Makarov, V. V. (2013). Tidal friction and tidal lagging. Applicability limitations of a popular formula for the tidal torque. *The Astrophysical Journal*, 764(1), 26. The CPL and CTL frequency dependence.
- Kaula, W. M. (1964). Tidal dissipation by solid friction and the resulting orbital evolution. *Reviews of Geophysics*, 2(4), 661-685. The tidal potential expansion.
- Efroimsky, M., and Williams, J. G. (2009). Tidal torques: A critical review of some techniques. *Celestial Mechanics and Dynamical Astronomy*, 104, 257-289. Eq. 18, the tidal potential with mode-dependent Love numbers.
- Peale, S. J., and Cassen, P. (1978). Contribution of tidal dissipation to lunar thermal history. *Icarus*, 36(2), 245-269. The synchronous heating formula.
