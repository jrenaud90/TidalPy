# Stellar (`stellar_x`)

_Updated: 2026-09-16_

`TidalPy.stellar_x` contains the stellar physics used to calculate stellar heating on a planet in long-term thermal-orbital evolution models. Its luminosity models map a star's mass \[kg\] onto its luminosity \[W\] using published relationships, which is helpful for exoplanet hosts whose stellar properties are not fully published.

| Page | Covers |
|---|---|
| [Luminosity Models](luminosity.md) | The three models, the Stefan-Boltzmann conversions they share, the factory, attaching a model to a star, serialization, the C++ surface, and how to add a model. |

```{toctree}
:maxdepth: 1

Luminosity Models <luminosity.md>
```

## Where Luminosity is Used

A luminosity model is attached to a `StarWorld` with `set_luminosity_model`. The star then derives its own luminosity and effective temperature from its mass and radius. See [Worlds](../structures_x/worlds/worlds.md).

The model is optional. A star given a luminosity or an effective temperature directly keeps the two consistent through the Stefan-Boltzmann relation without a model. A model lets the star's mass set those numbers instead, as a population study or an evolving system needs.

`System` uses the star's luminosity to compute the orbit-averaged flux at each of its worlds and the gray-body equilibrium temperature that follows from it. See [System](../structures_x/system/system.md).

## Scope

This module covers luminosity and the temperature conversions tied to it. Stellar structure, evolution along the main sequence, spectra, and the wavelength dependence of the radiation field are not included. The mass-to-luminosity relations are empirical fits to main-sequence stars.

## Examples

`Demos_x/Basics/02_world_building.ipynb` builds a star and reads its luminosity and effective temperature, and `Demos_x/Physics/04_orbits_insolation.ipynb` uses the luminosity for insolation and equilibrium temperatures across stellar types.

## References

- Cuntz, M., and Wang, Z. (2018). The mass-luminosity relation for a refined set of late-K/M stars. *Research Notes of the AAS*, 2(1), 19, [doi:10.3847/2515-5172/aaaa67](https://doi.org/10.3847/2515-5172/aaaa67).
- Salaris, M., and Cassisi, S. (2005). *Evolution of Stars and Stellar Populations*. Main-sequence mass-luminosity scaling across the mass range.
