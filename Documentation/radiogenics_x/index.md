# Radiogenics (`radiogenics_x`)

_Updated: 2026-09-12_

`TidalPy.radiogenics_x` adds functionality to calculate internal heating due to the decay of radioactive isotopes (both long- and short-duration isotopes). Each model in this module uses a layer's mass and the elapsed time to find the radiogenic heating $Q$ [W] released inside that layer.

| Page | Covers |
|---|---|
| [Radiogenic Models](radiogenics_models.md) | The three models, the isotope value type, the built-in literature datasets, the factory, vectorized and one-shot evaluation, serialization, and how to add a model. |

```{toctree}
:maxdepth: 1

Radiogenic Models <radiogenics_models.md>
```

## Where Radiogenics is Used

A radiogenics model is attached to a `SolidLiquidLayer` with `set_radiogenics`, alongside the layer's cooling model. See [SolidLiquidLayer](../structures_x/layers/solidliquid_layer.md). The layer then answers `calc_radiogenic_heating(time, mass)`, and `LayeredWorld.calc_internal_heating(time)` sums the contributions of every layer that carries a model. Layers without one contribute zero rather than raising.

When a world is built from a TOML file or a config dict, the `[layers.<name>.radiogenics]` table names the model and its parameters, and anything the user omits falls back to the material defaults in `TidalPy_Configs_x.toml`. The shipped defaults give a rock mantle the chondritic isotope set and turn radiogenics off in iron cores and ice shells. See the [TOML schema](../structures_x/config/toml_schema.md).

Radiogenic heating is deliberately kept separate from the tidal solve. The two sources are computed independently and summed by whatever drives the thermal state, so a study can hold one fixed while varying the other.

## References

- Turcotte, D. L., and Schubert, G. (2002). *Geodynamics*, second edition. Radiogenic heat production rates and the crustal heat budget.
- Hussmann, H., and Spohn, T. (2004). Thermal-orbital evolution of Io and Europa. *Icarus*, 171(2), 391-410. Chondritic isotope abundances for icy satellites.
- Castillo-Rogez, J. C., et al. (2007). Iapetus' geophysics: Rotation rate, shape, and equatorial ridge. *Icarus*, 190(1), 179-202. Long- and short-lived isotope inventories for early thermal evolution.
- McDonough, W. F., and Sun, S.-s. (1995). The composition of the Earth. *Chemical Geology*, 120(3-4), 223-253. Bulk silicate Earth elemental abundances.
