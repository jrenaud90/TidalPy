# Tides (`Tides_x`)

_Updated: 2026-09-12_

`TidalPy.Tides_x` turns an orbital state into tidal dissipation. It carries the angular ingredients of the tidal potential (the eccentricity and obliquity functions), the global one-dimensional dissipation models that collapse the potential into heating and orbital derivatives, the depth-resolved three-dimensional stress, strain, and heating kernel, and the Love-number container shared with the radial solver.

| Page | Covers |
|---|---|
| [Global Tides](global_tides.md) | The tide models (`rheology`, `cpl`, `ctl`, `ctl_q`), the mode collapse, and the world's `calc_tides`. |
| [3D Tidal Heating](multilayer_3d_heating.md) | The depth-resolved kernel, the secular and instantaneous heating paths, and the displacement grid. |
| [Love Numbers](love/love_numbers.md) | The `LoveNumbers` container, the Love-number solution methods, and the closed-form homogeneous-sphere formulas. |
| [Eccentricity Functions](eccentricity.md) | $G_{l,p,q}(e)$, the truncation levels, and what each costs. |
| [Obliquity Functions](obliquity.md) | $F_{l,m,p}(I)$, the truncation levels, and why they matter even at zero obliquity. |

```{toctree}
:maxdepth: 1

Global Tides <global_tides.md>
3D Tidal Heating <multilayer_3d_heating.md>
Love Numbers <love/love_numbers.md>
Eccentricity Functions <eccentricity.md>
Obliquity Functions <obliquity.md>
```

## How the pieces fit

A tidal solve starts from the potential. The Kaula expansion writes it as a sum over modes indexed by $(l, m, p, q)$, each carrying an amplitude built from an obliquity function $F_{lmp}(I)$, an eccentricity function $G_{lpq}(e)$, and a forcing frequency $\omega_{lmpq}$ set by the orbital and spin rates. The truncation levels decide how many of those modes are kept, and therefore both the accuracy and the cost of everything downstream.

Each active mode then needs a response. That is the Love number $k_l(\omega)$, which comes either from a radial solve of the interior or from one of the analytic models. Summing $-\mathrm{Im}[k_l]$ against the per-mode potential terms gives the global heating and the orbital derivatives; carrying the full radial functions instead gives the depth-resolved strain, stress, heating, and displacement fields.

## Reference

- Kaula, W. M. (1964). Tidal dissipation by solid friction and the resulting orbital evolution. *Reviews of Geophysics*, 2(4), 661-685.
- Renaud, J. P., et al. (2021). Tidal dissipation in dual-body, highly eccentric, and non-synchronously rotating systems. *The Astrophysical Journal*, 902(2), 122.
