# Viscosity (`Viscosity`)

_Updated: 2026-10-01_

`TidalPy.Viscosity` contains viscosity models which map a temperature \[K\] and a pressure \[Pa\] onto a dynamic viscosity \[Pa s\].

A silicate mantle's viscosity changes by about ten orders of magnitude across the temperature range a tidally heated body can occupy, while its shear modulus changes by less than one, so the viscosity law and its activation energy usually matter more to a predicted heating rate than any other input.

| Page | Covers |
|---|---|
| [Viscosity Models](viscosity_models.md) | The models, their formulas and parameters, the factory, the Python and C++ surfaces, and how to add a model. |

```{toctree}
:maxdepth: 1

Viscosity Models <viscosity_models.md>
```

## Where Viscosity is Used

Each phase of a layer's material holds a viscosity model for its shear response and one for its bulk response, its `shear_viscosity` and `bulk_viscosity` slots (see [Phases and Materials](../Material/materials.md)). The liquid phase's shear viscosity is the melt's. During a whole-planet equation-of-state solve the material evaluates them at the local temperature and pressure, and its [melt weakening](../PartialMelt/partial_melt_models.md#melt-weakening) then lowers the viscosity and the shear modulus wherever melt is present. The post-melt values are then used by [`Rheology`](../Rheology/index.md) to calculate a complex modulus.

Viscosity is frequency-independent, so it is resolved once per equation-of-state solve and reused across every tidal forcing frequency.

A phase with no viscosity model reports a NaN viscosity, so a layer that tides through it with a viscoelastic rheology needs one. A viscosity that varies with radius inside a layer, such as a seismic profile's, is an `interpolate` model.

## Examples

`Demos/Systems/12_thermal_orbital_evolution.ipynb` builds a temperature-dependent viscosity with `make_viscosity` and feeds it through a Maxwell rheology and the radial solver.

## References

- Moore, W. B. (2006). Thermal equilibrium in Europa's ice shell. *Icarus*, 180(1), 141-146.
- Henning, W. G., O'Connell, R. J., and Sasselov, D. D. (2009). Tidally heated terrestrial exoplanets: Viscoelastic response models. *The Astrophysical Journal*, 707(2), 1000-1015.
