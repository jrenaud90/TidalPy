# Global (1D) Tidal Dissipation (`Tides_x.classes`)

_Updated: 2026-09-19_

The global (or "1D potential") approach computes a body's total tidal heating and the three orbital potential derivatives (`dU/dM`, `dU/dw`, `dU/dO`) by summing over the active tidal modes `(l, m, p, q)`. Each mode carries a forcing frequency $\omega_{lmpq} = (l - 2p + q)\,n - m\,\dot{\theta}$ and a precomputed potential weight; a tide model supplies the per-mode dissipation multiplier $-\mathrm{Im}[k_{l}(\omega)]$ that the collapse multiplies in and sums. Harmonic degrees `l = 2..10` are supported.

The model-independent per-mode weights (the eccentricity functions $G_{lpq}$, the obliquity functions $F_{lmp}$, and the common coefficient $G_{lpq}^{2}F_{lmp}^{2}\frac{(l-m)!}{(l+m)!}\left(2-\delta_{0m}\right)\left(\frac{R}{a}\right)^{2l+1}\frac{G M_{h}}{a}$) come from `c_global_potential`; the Physics section below derives them. This page documents the tide models and the collapse that turn those weights into heating and torque.

The full complex Love-number suite (k, h, l) is always the transport type ([`c_LoveNumbers`](love/love_numbers.md)), even though only `k` drives heating and orbital dynamics, so the displacement Love numbers from the radial solver are never discarded. The analytic models cannot produce `h` and `l` (no radial solution) and return them as `NaN`.

## Physics

### Tidal Potential

A host of mass $M_h$ on an orbit of semi-major axis $a$, eccentricity $e$, and obliquity $I$ relative to the body's equator raises, at the surface of a body of radius $R$, the tidal potential (Kaula 1964; Efroimsky and Williams 2009, Eq. 18)

$$U(\theta, \phi, t) = \frac{G M_h}{a}\sum_{l=2}^{\infty}\left(\frac{R}{a}\right)^{l}\sum_{m=0}^{l}\frac{(l-m)!}{(l+m)!}\left(2-\delta_{0m}\right)P_{lm}(\cos\theta)\sum_{p=0}^{l}F_{lmp}(I)\sum_{q=-\infty}^{\infty}G_{lpq}(e)\,\mathcal{T}_{lm}\!\left(\omega_{lmpq}t - m\phi\right),$$

where $\theta$ is the colatitude, $\phi$ the east longitude, $P_{lm}$ the associated Legendre functions without the Condon-Shortley phase, $\mathcal{T}_{lm}$ is $\cos$ for even $l - m$ and $\sin$ for odd $l - m$, and $F_{lmp}$ and $G_{lpq}$ are the [obliquity](obliquity.md) and [eccentricity](eccentricity.md) functions. Each $(l, m, p, q)$ is a tidal mode with the forcing frequency

$$\omega_{lmpq} = (l - 2p + q)\,n - m\,\dot{\theta},$$

where $n$ is the orbital mean motion and $\dot{\theta}$ the spin rate \[rad s$^{-1}$\]. Periapse and node precession are neglected. The degree range and the two truncation levels decide which modes are kept.

### Dissipation and Orbital Derivatives

The body answers each mode with the complex Love number $k_l$ at the forcing frequency $\chi_{lmpq} = |\omega_{lmpq}|$, and the tide model supplies its dissipative part $K_{l} = -\mathrm{Im}[k_{l}(\chi_{lmpq})]$ (see Models below). Averaged over the orbit and over apsidal precession, the tidal heating $\dot{E}$ \[W\] and the derivatives of the tidal potential with respect to the mean anomaly $\mathcal{M}$, the argument of periapse $\varpi$, and the node $\Omega$ \[J kg$^{-1}$ rad$^{-1}$\] are (Renaud et al. 2021, Eq. 7)

$$\begin{bmatrix} \partial U / \partial \mathcal{M} \\ \partial U / \partial \varpi \\ \partial U / \partial \Omega \\ \dot{E} \end{bmatrix} = \frac{G M_h}{a}\sum_{l=2}^{l_{\max}}\left(\frac{R}{a}\right)^{2l+1}\sum_{m=0}^{l}\frac{(l-m)!}{(l+m)!}\left(2-\delta_{0m}\right)\sum_{p=0}^{l}F_{lmp}^{2}(I)\sum_{q}G_{lpq}^{2}(e)\begin{bmatrix} (l-2p+q)\,\mathrm{sgn}(\omega_{lmpq})\,K_{l} \\ (l-2p)\,\mathrm{sgn}(\omega_{lmpq})\,K_{l} \\ m\,\mathrm{sgn}(\omega_{lmpq})\,K_{l} \\ \chi_{lmpq}\,M_h\,K_{l} \end{bmatrix}.$$

Modes at zero frequency do not dissipate and are skipped. The derivatives set the orbital and spin rates (see [Dynamics](../dynamics_x/dynamics.md)). `collapse_global_tides` and a world's `calc_tides` return this sum; a layered world then splits the heating among its layers by their `tidal_scale`, so the whole-body sum itself uses the unscaled $K_l$.

For a synchronously rotating body at zero obliquity, only the $q = \pm 1$ modes of degree 2 dissipate to leading order in $e$, and with $K_2 = k_2/Q$ the sum reduces to the constant-phase-lag heating (Peale and Cassen 1978)

$$\dot{E} = \frac{21}{2}\,\frac{k_{2}}{Q}\,\frac{G M_h^{2} R^{5} n\, e^{2}}{a^{6}},$$

Using TidalPy with an eccentricity truncation of 1 reproduces this classic formula.

## Models

| Model | Alias | Complex Love number $k_{l}(\omega)$ | $-\mathrm{Im}[k_{l}]$ |
|-------|-------|----------------------------------|------------|
| `RheologyTide` | `rheology` | supplied by the radial solver | $-\mathrm{Im}[k_{l}]$ from the solver |
| `FixedQTide` | `cpl`, `fixed_q` | $k_{l}\,(1 - i/Q_{l})$ | $k_{l}/Q_{l}$ (frequency independent) |
| `FixedLagTide` | `ctl`, `fixed_dt` | $k_{l}\,(1 - i\,\omega\,\Delta t_{l})$ | $k_{l}\,\omega\,\Delta t_{l}$ |
| `CTLQTide` | `ctl_q`, `fixed_dt_q` | $k_{l}\,(1 - i\,\omega\,\Delta t_{l}/Q_{l})$ | $k_{l}\,\omega\,\Delta t_{l}/Q_{l}$ |

Fixed per-degree parameters $k_{l}$ (static Love number, `fixed_k`), $Q_{l}$ (quality factor, `fixed_q`), and $\Delta t_{l}$ (time lag \[s\], `fixed_dt`) are supplied as lists indexed from degree `l = 2` (index 0 is `l = 2`). The `rheology` model needs the radial solver and is driven by the world's `calc_tides` method, not the standalone collapse below.

A zero or absent $Q_{l}$ is treated as purely elastic (no dissipation) rather than a divide by zero.

## Python API

```python
from TidalPy.Tides_x.classes import (
    RheologyTide, FixedQTide, FixedLagTide, CTLQTide,
    make_tide, collapse_global_tides)
from TidalPy.constants import G

# Build a model directly or by name (aliases, case-insensitive):
tide = make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [50.0]})

love    = tide.calc_love_numbers(degree_l=2, frequency=4.1e-5)   # LoveNumbers(k=0.3-0.006j, h=nan, l=nan)
neg_imk = tide.calc_neg_imk(degree_l=2, frequency=4.1e-5)        # 0.006
# The rheology model returns the supplied radial-solver suite unchanged:
#   make_tide("rheology").calc_love_numbers(2, w, solver_love) -> solver_love (k, h, l)

# Standalone global collapse for an analytic model:
# The orbital-state arguments follow the world's calc_tides order.
result = collapse_global_tides(
    planet_radius=1.82e6,
    orbital_frequency=4.11e-5,
    spin_frequency=4.11e-5,   # synchronous
    eccentricity=0.0041,
    obliquity=0.0,
    semi_major_axis=4.22e8,
    host_mass=1.898e27,
    G_to_use=G,
    tide_model="cpl",
    tide_config={"fixed_k": [0.3], "fixed_q": [50.0]},
    max_degree_l=2,
    eccentricity_truncation=2)
# {"tidal_heating": W, "dUdM": ..., "dUdw": ..., "dUdO": ..., "num_modes": int}
```

**Constructors**

- `RheologyTide()`
- `FixedQTide(fixed_k=None, fixed_q=None)`
- `FixedLagTide(fixed_k=None, fixed_dt=None)`
- `CTLQTide(fixed_k=None, fixed_dt=None, fixed_q=None)`

where each `fixed_*` is a list indexed from `l = 2`.

**Methods and properties**

| Member | Returns | Description |
|--------|---------|-------------|
| `calc_love_numbers(degree_l, frequency, solver_love=None)` | `LoveNumbers` | `(k, h, l)`; analytic models set `h, l = NaN`. `solver_love` (a `LoveNumbers`) is returned unchanged by `RheologyTide`. |
| `calc_neg_imk(degree_l, frequency, solver_love=None)` | float | Dissipation multiplier `−Im[k_l]`. |
| `needs_radial_solve` | bool | `True` only for `RheologyTide`. |
| `get_fixed_k(degree_l)` | float | Static Love number for that degree (analytic models). |
| `get_fixed_q(degree_l)`, `get_fixed_dt(degree_l)` | float | Quality factor and time lag [s] for that degree. Defined on every model: a model that does not carry the parameter returns NaN, which is how a world's `cpl` or `ctl` Love method decides whether it can fall back to the attached tide model. See [Love numbers](love/love_numbers.md). |
| `model_name` | str | The model's registered name, for example `cpl`. |
| `get_config_dict()` | dict | Model name plus per-degree parameters. |
| `save_config(path)`, `get_schema_version_str()` | — | Configuration output and schema version, shared by every physics model; see [Base Classes](../utilities_x/classes_x.md). |
| `save_binary(path)` / `load_binary(path)` | — | Inherited from `TidalPyBaseClass`. |

`make_tide(name, config=None)` returns the concrete subclass; unknown names, and config keys other than `fixed_k`, `fixed_q`, and `fixed_dt`, raise `ValueError`. `collapse_global_tides(...)` supports the analytic models only: the `rheology` model raises `NotImplementedError` (use the world's `calc_tides`).

## C++ API

The C++ layer is canonical; the Cython classes are thin adapters.

- Config struct `tidalpy::c_TideModelConfig` (`tide_.hpp`): per-degree `std::vector<double>` `fixed_k`, `fixed_q`, `fixed_dt` (indexed from `l = 2`; entries beyond `l = 10` are ignored).
- Base class `c_TideBase : c_PhysicsBase` (`tide_base_.hpp`): pure virtual `calc_love_numbers(int degree_l, double frequency, const c_LoveNumbers& solver_love) const` (returns the full `c_LoveNumbers` suite) and `needs_radial_solve() const`; non-virtual `calc_neg_imk(...)`, which is `−Im[calc_love_numbers(...).k]`.
- Models (`tide_.hpp`): `c_RheologyTide`, `c_FixedQTide`, `c_FixedLagTide`, `c_CTLQTide`, each with per-degree getters. Binary class ids 901–904.
- Enum factory: `c_tide_model_from_name(const std::string&)` returns a `c_TideModel` (`Rheology`, `FixedQ`, `FixedLag`, `CTLQ`) and throws `std::invalid_argument` on an unknown name; `c_find_tide(c_TideModel, const c_TideModelConfig&)` returns a `std::unique_ptr<c_TideBase>` (a name overload also exists); `c_tide_from_binary(std::istream&, bool force=false)` peeks the class id, builds the model, and calls `read_binary`.
- Collapse (`tide_collapse_.hpp`): `c_collapse_global_tides(const c_GlobalPotentialStorage&, const c_TideBase&, const c_IntMap<c_Key4, c_LoveNumbers>* solver_love_by_lmpq = nullptr)` returns a `c_GlobalTideResult`. The `solver_love_by_lmpq` map supplies the radial-solver Love numbers (k, h, l) per mode for the `rheology` model; pass `nullptr` for the analytic models. `c_GlobalTideResult` holds `tidal_heating`, `dU_dM`, `dU_dw`, `dU_dO`, `num_modes`, and `error_code`.

## Adding a New Model

1. Add a `c_<Name>Tide : c_TideBase` in `tide_.hpp` implementing `calc_love_numbers` (return a `c_LoveNumbers`; set `h, l = NaN` if no radial solution), `needs_radial_solve`, `write_binary`, `read_binary`.
2. Add any new parameters to `c_TideModelConfig`.
3. Register a `c_TideModel::<Name>` enum value and wire it into `c_tide_model_from_name`, `c_find_tide`, and `c_tide_from_binary`.
4. Reserve a unique `BinaryClassID` (next free in the 900-block) in `Utilities_x/binary_x/binary_.hpp`.
5. Add the Cython `cdef class`, factory branch, config-dict, and tests, and document the model here.

## References

- Renaud, J. P., et al. (2021). Tidal dissipation in dual-body, highly eccentric, and nonsynchronously rotating systems: Applications to Pluto-Charon and the exoplanet TRAPPIST-1e. *The Planetary Science Journal*, 2(1), 4. The collapse form of global dual-body dissipation.
- Efroimsky, M., and Makarov, V. V. (2013). Tidal friction and tidal lagging. Applicability limitations of a popular formula for the tidal torque. *The Astrophysical Journal*, 764(1), 26. The CPL and CTL frequency dependence.
- Kaula, W. M. (1964). Tidal dissipation by solid friction and the resulting orbital evolution. *Reviews of Geophysics*, 2(4), 661-685. The tidal potential expansion.
- Efroimsky, M., and Williams, J. G. (2009). Tidal torques: A critical review of some techniques. *Celestial Mechanics and Dynamical Astronomy*, 104, 257-289. Eq. 18, the tidal potential with mode-dependent Love numbers.
- Peale, S. J., and Cassen, P. (1978). Contribution of tidal dissipation to lunar thermal history. *Icarus*, 36(2), 245-269. The synchronous heating formula.
