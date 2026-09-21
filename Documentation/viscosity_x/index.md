# Viscosity (`viscosity_x`)

_Updated: 2026-09-16_

`TidalPy.viscosity_x` contains viscosity models which map a temperature \[K\] and a pressure \[Pa\] onto a dynamic viscosity \[Pa s\].

A silicate mantle's viscosity changes by about ten orders of magnitude across the temperature range a tidally heated body can occupy, while its shear modulus changes by less than one, so the viscosity law and its activation energy usually matter more to a predicted heating rate than any other input.

| Page | Covers |
|---|---|
| [Viscosity Models](viscosity_models.md) | The three models, their formulas and parameters, the factory, the Python and C++ surfaces, and how to add a model. |

```{toctree}
:maxdepth: 1

Viscosity Models <viscosity_models.md>
```

## Where Viscosity is Used

A layer's material (its EOS model) holds a viscosity model for its shear response and one for its bulk response, attached with `set_shear_viscosity` and `set_bulk_viscosity` on the EOS model or through the layer's helpers of the same name. During a whole-planet equation-of-state solve the viscosity model converts each radial slice's temperature and pressure into that slice's pre-melt viscosity, and the [partial-melt](../partial_melt_x/partial_melt_models.md) model then weakens both the viscosity and the shear modulus wherever melt is present. The post-melt values are what [`rheology_x`](../rheology_x/index.md) consumes to produce a complex modulus.

Viscosity is frequency-independent, so it is resolved once per equation-of-state solve and reused across every tidal forcing frequency.

A layer with a rheology but no viscosity model falls back to its static viscosity, which is NaN unless one was supplied at construction. An equation of state that supplies its own viscosity profile as extra output overrides the model's value slice by slice.

## Examples

`Demos_x/Systems/12_thermal_orbital_evolution.ipynb` builds a temperature-dependent viscosity with `make_viscosity` and feeds it through a Maxwell rheology and the radial solver.

## References

- Moore, W. B. (2006). Thermal equilibrium in Europa's ice shell. *Icarus*, 180(1), 141-146.
- Henning, W. G., O'Connell, R. J., and Sasselov, D. D. (2009). Tidally heated terrestrial exoplanets: Viscoelastic response models. *The Astrophysical Journal*, 707(2), 1000-1015.
