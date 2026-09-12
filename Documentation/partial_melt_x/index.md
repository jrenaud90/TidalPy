# Partial Melting (`partial_melt_x`)

_Updated: 2026-09-12_

`TidalPy.partial_melt_x` handles what happens to a material's strength once part of it has melted. Each model takes the pre-melt (solid) viscosity and shear modulus at a point, together with the temperature, and returns the post-melt values plus the volumetric melt fraction.

Melt weakening is the feedback that makes tidal heating self-limiting. Heating raises the temperature, the temperature raises the melt fraction, the melt fraction drops the viscosity and shear modulus by orders of magnitude, and a weaker body deforms more but dissipates differently. Whether a body runs away to a magma ocean or settles into a warm steady state is largely decided by the shape of the weakening curve near the critical melt fraction, which is why the models differ most sharply right there.

| Page | Covers |
|---|---|
| [Partial-Melt Models](partial_melt_models.md) | The melt-fraction definition, the three models, their parameters, the Python and C++ surfaces, and how to add a model. |

```{toctree}
:maxdepth: 1

Partial-Melt Models <partial_melt_models.md>
```

## Where melt weakening fits

A partial-melt model is attached to a layer with `set_partial_melt` and applied during the whole-planet equation-of-state solve. At each radial slice the [viscosity model](../viscosity_x/index.md) supplies the pre-melt viscosity and the equation of state supplies the pre-melt moduli; the melt model then rewrites both, once for the shear pair and once for the bulk pair. The post-melt values are what [`rheology_x`](../rheology_x/index.md) turns into a complex modulus, and a layer keeps both sets so you can compare them (`get_premelt_shear_viscosity` against `get_shear_viscosity`).

Like viscosity, melt weakening is frequency-independent and therefore resolved once per equation-of-state solve rather than once per tidal mode.

## References

- Fischer, H.-J., and Spohn, T. (1990). Thermal-orbital histories of viscoelastic models of Io. *Icarus*, 83(1), 39-65.
- Henning, W. G., O'Connell, R. J., and Sasselov, D. D. (2009). Tidally heated terrestrial exoplanets: Viscoelastic response models. *The Astrophysical Journal*, 707(2), 1000-1015.
- Renaud, J. P., and Henning, W. G. (2018). Increased tidal dissipation using advanced rheological models: Implications for Io and tidally active exoplanets. *The Astrophysical Journal*, 857(2), 98.
