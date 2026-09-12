# Viscosity (`viscosity_x`)

_Updated: 2026-09-12_

`TidalPy.viscosity_x` turns local conditions into a viscosity. Each model maps a temperature [K] and a pressure [Pa] onto a dynamic viscosity [Pa s], which is the single material property that decides how readily a planet's interior flows and therefore how much tidal energy it converts into heat.

Viscosity deserves its own module because it is the most uncertain and the most strongly varying quantity in the whole calculation. A silicate mantle's viscosity changes by ten orders of magnitude across the temperature range a tidally heated body can occupy, while its shear modulus changes by less than one. The choice of viscosity law, and the activation energy inside it, usually matters more to a predicted heating rate than any other input.

| Page | Covers |
|---|---|
| [Viscosity Models](viscosity_models.md) | The three models, their formulas and parameters, the factory, the Python and C++ surfaces, and how to add a model. |

```{toctree}
:maxdepth: 1

Viscosity Models <viscosity_models.md>
```

## Where viscosity fits

The viscosity models sit one step before the rheology models. A layer holds a viscosity model for its shear response and one for its bulk response, attached with `set_shear_viscosity` and `set_bulk_viscosity`. During a whole-planet equation-of-state solve each radial slice arrives with a temperature and a pressure, the viscosity model converts them into that slice's pre-melt viscosity, and the [partial-melt](../partial_melt_x/partial_melt_models.md) model then weakens both the viscosity and the shear modulus wherever melt is present. The resulting post-melt values are what [`rheology_x`](../rheology_x/index.md) consumes to produce a complex modulus.

Viscosity is frequency-independent, which is why it can be resolved once per equation-of-state solve and reused across every tidal forcing frequency. That separation is what makes a many-mode tidal calculation affordable.

A layer with a rheology but no viscosity model falls back to its static viscosity, which is NaN unless one was supplied at construction. An equation of state that supplies its own viscosity profile as extra output overrides the model's value slice by slice.

## References

- Moore, W. B. (2006). Thermal equilibrium in Europa's ice shell. *Icarus*, 180(1), 141-146.
- Henning, W. G., O'Connell, R. J., and Sasselov, D. D. (2009). Tidally heated terrestrial exoplanets: Viscoelastic response models. *The Astrophysical Journal*, 707(2), 1000-1015.
