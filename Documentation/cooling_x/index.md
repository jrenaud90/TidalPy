# Cooling (`cooling_x`)

_Updated: 2026-09-12_

`TidalPy.cooling_x` answers the other half of a thermal history. Tidal dissipation and radioactive decay put energy into a layer; a cooling model says how fast that layer can get rid of it. Each model maps the layer's physical state onto a surface heat flux [W m$^{-2}$], the thickness of the thermal boundary layer that carries it, and the Rayleigh and Nusselt numbers that describe the transport regime.

This is what closes the thermal-orbital feedback loop. A layer that heats faster than it cools warms up, which drops its viscosity, which changes both its tidal dissipation and its ability to convect. Whether that feedback settles into an equilibrium or runs away depends on the relative slopes of the heating and cooling curves, so a cooling model is not a detail bolted onto the end of a tidal calculation. It is half the physics.

| Page | Covers |
|---|---|
| [Cooling Models](cooling_models.md) | The three models, their inputs and result structure, vectorized and one-shot evaluation, serialization, and how to add a model. |

```{toctree}
:maxdepth: 1

Cooling Models <cooling_models.md>
```

## Where cooling fits

A cooling model is attached to a `SolidLiquidLayer` with `set_cooling`, alongside the layer's radiogenic-heating model. See [SolidLiquidLayer](../structures_x/layers/solidliquid_layer.md). Unlike the viscosity and melt models, cooling is not evaluated during the equation-of-state solve: it depends on a temperature drop across the layer, which is a property of the thermal state being evolved rather than of the static structure.

The models are parameterized, not resolved. They reduce the whole of mantle convection to a boundary-layer scaling with a Rayleigh number and a handful of fitted constants, which is the standard approach for the timescales planetary evolution deals in. When their assumptions fail, usually in a thin layer, at a vanishing temperature drop, or in a regime where the scaling constants were never calibrated, the models degrade to a defined edge case rather than to a numerical error.

## References

- Turcotte, D. L., and Schubert, G. (2002). *Geodynamics*, second edition. Rayleigh and Nusselt convection scaling, and conduction.
- Solomatov, V. S. (1995). Scaling of temperature- and stress-dependent viscosity convection. *Physics of Fluids*, 7(2), 266-274.
- Schubert, G., Turcotte, D. L., and Olson, P. (2001). *Mantle Convection in the Earth and Planets*. Boundary-layer theory.
