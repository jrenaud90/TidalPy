# TidalPy.tides Documentation

**TidalPy's classic Tides module**

[Auto Generated API](https://tidalpy.readthedocs.io/en/latest/API/generated/TidalPy.tides.html)

`TidalPy.tides` is the tidal machinery of the classic backend: the tidal potential and its mode decomposition, the analytic Love-number formulas, the dissipation and heating calculators, and the world-attached tide methods.

| Piece | What it holds |
|---|---|
| `tides.love1d` | Analytic homogeneous Love numbers: `complex_love`, `static_love`, `effective_rigidity`, and their `_general` degree-l forms. |
| `tides.methods` | `TidesBase` and the two implementations, `GlobalApproxTides` (fixed Q or fixed time lag) and `LayeredTides`. |
| `tides.dissipation`, `tides.heating` | `calc_tidal_susceptibility` and `calculate_volumetric_heating`. |
| `tides.modes`, `tides.potential` | Mode bookkeeping and the tidal potential forms. |
| `tides.eccentricity_funcs`, `tides.inclination_funcs` | The pre-computed eccentricity and inclination (obliquity) terms, selected by degree and truncation. |
| `tides.multilayer` | Grid-based stress, strain, displacement, and heating from a radial solution. |
| `tides.ctl_funcs` | Frequency laws for the constant-time-lag models. |

```{toctree}
:maxdepth: 2
:caption: Contents

Eccentricity Functions <Eccentricity.md>
Inclination Functions <Obliquity.md>
```

## Moving to the new backend

The replacement is [`TidalPy.Tides_x`](../Tides_x/index.md), which is where new development happens. It covers the same ground with a C++ core: global dissipation models, the depth-resolved three-dimensional kernel, the Love-number container and solution methods, and its own [eccentricity](../Tides_x/eccentricity.md) and [obliquity](../Tides_x/obliquity.md) functions. The [porting guide](../future_structure.md) maps the old names onto the new ones.
