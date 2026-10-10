# Radiogenics (`Radiogenics`)

_Updated: 2026-10-07_

`TidalPy.Radiogenics` calculates internal heating from the decay of radioactive isotopes, long- and short-lived. Each model uses a layer's mass and the elapsed time to find the radiogenic heating $Q$ \[W\] released inside that layer.

| Page | Covers |
|---|---|
| [Radiogenic Models](radiogenics_models.md) | The three models and their parameters, the isotope tables, the built-in literature datasets, vectorized and one-shot evaluation, and the C++ API. |

```{toctree}
:maxdepth: 1

Radiogenic Models <radiogenics_models.md>
```

## Where Radiogenics is Used

A layer holds at most one radiogenics model, alongside its cooling model: `Layer(..., radiogenics=...)` or `layer.radiogenics = ...` takes a model, a model name, or a config table, and `None` removes it (see [Layer](../Structures/layers/layer.md)). The layer then provides `calc_radiogenic_heating(time, mass)`, and `BaseWorld.calc_internal_heating(time)` sums every layer that carries a model; layers without one contribute zero. In a world's thermal solve the model heats its layer when the layer sets `use_heating` (see [Heat Sources](../Structures/worlds/worlds.md#heat-sources)).

In a TOML file or config dict, the `[layers.<name>.radiogenics]` table names the model and its parameters; a layer without one has no radiogenics. An `isotope` model given neither a dataset nor isotope tables takes the `[radiogenics] isotopes` dataset of `TidalPy_Configs.toml` (`modern_day_chondritic` by default). See the [TOML schema](../Structures/config/toml_schema.md).

## Examples

`Demos/Physics/P10_thermal_interior.ipynb` compares the radiogenic models over the age of the Solar System and attaches an isotope dataset to a layer of the bundled Io.

## References

- Turcotte, D. L., and Schubert, G. (2002). *Geodynamics*, second edition. Radiogenic heat production rates and the crustal heat budget.
- Hussmann, H., and Spohn, T. (2004). Thermal-orbital evolution of Io and Europa. *Icarus*, 171(2), 391-410. Chondritic isotope abundances for icy satellites.
- Castillo-Rogez, J. C., et al. (2007). Iapetus' geophysics: Rotation rate, shape, and equatorial ridge. *Icarus*, 190(1), 179-202. Long- and short-lived isotope inventories for early thermal evolution.
- McDonough, W. F., and Sun, S.-s. (1995). The composition of the Earth. *Chemical Geology*, 120(3-4), 223-253. Bulk silicate Earth elemental abundances.
