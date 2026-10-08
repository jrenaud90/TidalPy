# Tides (`Tides`)

_Updated: 2026-10-06_

`TidalPy.Tides` calculates tidal dissipation from an orbital state: the eccentricity and obliquity functions that drive the tidal potential, the global (1D) models that collapse it into heating and orbital derivatives, the depth-resolved 3D stress, strain, and heating kernel, and the Love-number container shared with the radial solver.

| Page | Covers |
|---|---|
| [Global Tides](global_tides.md) | The tide models (`rheology`, `cpl`, `ctl`, `ctl_q`), the mode collapse, and the world's `calc_tides`. |
| [3D Tidal Heating](multilayer_3d_heating.md) | The depth-resolved kernel, the secular and instantaneous heating paths, the displacement grid, and the stress and strain grids. |
| [Love Numbers](love/love_numbers.md) | The `LoveNumbers` container, the Love-number solution methods, and the closed-form homogeneous-sphere formulas. |
| [Eccentricity Functions](Eccentricity.md) | $G_{l,p,q}(e)$, the truncation levels, and what each costs. |
| [Obliquity Functions](Obliquity.md) | $F_{l,m,p}(I)$, the truncation levels, and why they matter even at zero obliquity. |

```{toctree}
:maxdepth: 1

Global Tides <global_tides.md>
3D Tidal Heating <multilayer_3d_heating.md>
Love Numbers <love/love_numbers.md>
Eccentricity Functions <Eccentricity.md>
Obliquity Functions <Obliquity.md>
```

## Imports

The package re-exports its public entry points:

```python
from TidalPy.Tides import (
    make_tide, collapse_global_tides,                   # global tide models and the mode collapse
    LoveNumbers, calc_homogeneous_love_numbers,         # Love-number container and closed forms
    global_potential, tidal_potential_3d_modes,         # tidal potentials
    strain_stress_heating_point, volumetric_heating,    # point-wise 3D kernels
    eccentricity_func, obliquity_func)                  # the functions behind the potential
```

Each also lives in a subpackage (`classes`, `love`, `potential`, `multilayer`, `eccentricity`, `obliquity`), which the pages below import from.

## Structure

A tidal solve starts from the potential. The Kaula expansion writes it as a sum over modes $(l, m, p, q)$, each with an amplitude built from an obliquity function $F_{lmp}(I)$ and an eccentricity function $G_{lpq}(e)$, and a forcing frequency $\omega_{lmpq}$ set by the orbital and spin rates. The truncation levels decide how many modes are kept, and so the accuracy and cost of everything downstream.

Each mode then needs a response, the Love number $k_l(\omega)$, from a radial solve of the interior or an analytic model. Summing $-\mathrm{Im}[k_l]$ against the per-mode potential terms gives the global heating and orbital derivatives; carrying the full radial functions instead gives the depth-resolved strain, stress, heating, and displacement fields.

## Where Tides are Used

A world's tide model is attached with `set_tide_model` and configured with `set_tide_config` or the `[tides]` TOML table; `calc_tides` runs the global solve and the `calc_3d_*` methods the depth-resolved one. See [Worlds](../Structures/worlds/worlds.md). A `System` supplies the orbital state for each call and turns the potential derivatives into orbital and spin rates; see [System](../Structures/system/system.md) and [Dynamics](../Dynamics/index.md).

## Examples

`Demos/Physics/P02_tidal_basics.ipynb` attaches a fixed-Q tide model and computes heating, `P03_rheology_io.ipynb` and `P04_gasgiant_fixedQ_dt.ipynb` compare tide models, and `P06_tidal_heating_3d.ipynb` and `P08_tidal_maps_3d.ipynb` map the three-dimensional heating, stress, and displacement.

## References

- Kaula, W. M. (1964). Tidal dissipation by solid friction and the resulting orbital evolution. *Reviews of Geophysics*, 2(4), 661-685.
- Renaud, J. P., et al. (2021). Tidal dissipation in dual-body, highly eccentric, and nonsynchronously rotating systems: Applications to Pluto-Charon and the exoplanet TRAPPIST-1e. *The Planetary Science Journal*, 2(1), 4.
