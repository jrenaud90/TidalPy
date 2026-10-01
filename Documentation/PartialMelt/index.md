# Partial Melting (`PartialMelt`)

_Updated: 2026-10-01_

`TidalPy.PartialMelt` has functionality to modify planetary material's strength once it begins to experience partial melt. Partial melt occurs in minerals where one element begins to melt while others remain solid. It is a fraction between 0 and 1, and for many materials has a critical melt fraction where a material switches from behaving like a solid with pockets of melt to a liquid with chunks of solid ($f_{crit}$ is often between 0.4 to 0.6 for rocks and ices). This transition causes a dramatic change in shear modulus and viscosity. 

Each model takes the pre-melt (solid) viscosity and shear modulus at a point, together with the temperature and pressure, and returns the post-melt values plus the volumetric melt fraction. The melt fraction is set by the solidus and liquidus, which are constant by default or can rise with pressure following a Simon and Glatzel law.

Partial melting is an important feedback that makes solid-body tidal heating self-limiting. Heating raises the temperature, the temperature raises the melt fraction, the melt fraction drops the viscosity and shear modulus by orders of magnitude, and a weaker body deforms more but has weaker dissipation (less friction). Whether a body runs away to a magma ocean or settles into a warm steady state is largely decided by the shape of the weakening curve near the critical melt fraction, which is why the models differ most sharply right there.

| Page | Covers |
|---|---|
| [Partial-Melt Models](partial_melt_models.md) | The melt-fraction definition, the pressure-dependent melting curves, the three models, their parameters, the Python and C++ surfaces, and how to add a model. |

```{toctree}
:maxdepth: 1

Partial-Melt Models <partial_melt_models.md>
```

## Where Partial Melt is Used

A partial-melt model is attached to a layer's material (its EOS model) with `set_partial_melt`, on the EOS model or through the layer's helper of the same name, and applied during the whole-planet equation-of-state solve. At each radial slice the [viscosity model](../Viscosity/index.md) supplies the pre-melt viscosity and the equation of state supplies the pre-melt moduli; the melt model then rewrites the shear modulus and viscosity and, behind switches that are off by default, mixes the melt into the density, weakens the bulk modulus (by a separate, much weaker law), and sets a compaction bulk viscosity. The post-melt values are what [`Rheology`](../Rheology/index.md) turns into a complex modulus, and a layer keeps both sets so you can compare them (`get_premelt_shear_viscosity` against `get_shear_viscosity`).

Like viscosity, melt weakening is frequency-independent and therefore resolved once per equation-of-state solve rather than once per tidal mode.

## Examples

`Demos/Physics/10_thermal_eos.ipynb` attaches a Henning melt model to a mantle and sweeps its temperature through the melting range.

## References

- Fischer, H.-J., and Spohn, T. (1990). Thermal-orbital histories of viscoelastic models of Io. *Icarus*, 83(1), 39-65.
- Henning, W. G., O'Connell, R. J., and Sasselov, D. D. (2009). Tidally heated terrestrial exoplanets: Viscoelastic response models. *The Astrophysical Journal*, 707(2), 1000-1015.
- Monteux, J., Andrault, D., and Samuel, H. (2016). On the cooling of a deep terrestrial magma ocean. *Earth and Planetary Science Letters*, 448, 140-149.
- Renaud, J. P., and Henning, W. G. (2018). Increased tidal dissipation using advanced rheological models: Implications for Io and tidally active exoplanets. *The Astrophysical Journal*, 857(2), 98.
- Simon, F., and Glatzel, G. (1929). Bemerkungen zur Schmelzdruckkurve. *Zeitschrift für anorganische und allgemeine Chemie*, 178(1), 309-316.
