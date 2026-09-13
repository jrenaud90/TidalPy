# Stellar (`stellar_x`)

_Updated: 2026-09-13_

`TidalPy.stellar_x` holds the stellar physics. Its job is to calculate stellar heating on a planet for use in long term thermal-orbital evolution models. It also contains helper functionality to turn a star's mass into its luminosity using published relationships. This is helpful when working with exoplanets where some stellar properties are not published.

| Page | Covers |
|---|---|
| [Luminosity Models](luminosity.md) | The three models, the Stefan-Boltzmann conversions they share, the factory, attaching a model to a star, serialization, the C++ surface, and how to add a model. |

```{toctree}
:maxdepth: 1

Luminosity Models <luminosity.md>
```

## Usage

A luminosity model is attached to a `StarWorld` with `set_luminosity_model`. The star then derives its own luminosity and effective temperature from its mass and radius. See [Worlds](../structures_x/worlds/worlds.md).

The model is optional. A star that is given a luminosity directly, or an effective temperature, stays internally consistent through the Stefan-Boltzmann relation without any model attached. A model is what lets the star's mass drive those numbers instead, which is what a population study or an evolving system needs.

Downstream, `System` uses the star's luminosity to compute the orbit-averaged flux at each of its worlds and the gray-body equilibrium temperature that follows from it. See [System](../structures_x/system/system.md). The star is where the energy budget starts; everything the tidal machinery computes is added on top of it.

## Scope

This module covers luminosity and the temperature conversions tied to it. Nothing related to stellar structure, evolution along the main sequence, spectra, and the wavelength dependence of the radiation field are currently included in TidalPy. The mass-to-luminosity relations here are empirical fits to main-sequence stars and say nothing about how a star got there or where it goes next.

## References

- Cuntz, M., and Wang, Z. (2018). The mass-luminosity relation for a refined set of late-K/M stars. *Research Notes of the AAS*, 2(1), 19, [doi:10.3847/2515-5172/aaaa67](https://doi.org/10.3847/2515-5172/aaaa67).
- Salaris, M., and Cassisi, S. (2005). *Evolution of Stars and Stellar Populations*. Main-sequence mass-luminosity scaling across the mass range.
