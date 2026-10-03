# Cooling (`Cooling`)

_Updated: 2026-10-01_

`TidalPy.Cooling` contains cooling models which quantify how heat moves through a layer. Each model maps a layer's physical state onto a surface heat flux [W m$^{-2}$], the thickness of the thermal boundary layer that carries it, and the Rayleigh and Nusselt numbers that describe the transport regime. In a world's thermal solve each model also builds its layer's temperature profile: one temperature throughout (`off`), two conducting halves (`conduction`), or conducting boundary layers around an adiabatic interior whose top is at the layer's temperature (`convection`), with a magma-ocean scaling when that interior is liquid.

| Page | Covers |
|---|---|
| [Cooling Models](cooling_models.md) | The three models, their flux laws and the profiles they build, the convective reference state and the magma-ocean branch, vectorized and one-shot evaluation, serialization, and how to add a model. |

```{toctree}
:maxdepth: 1

Cooling Models <cooling_models.md>
```

## Where Cooling is Used

A layer holds at most one cooling model: `Layer(..., cooling=...)` or `layer.cooling = ...` takes a model, a model name, or a config table, and `layer.cooling = None` removes it. A world built from a TOML file reads it from the layer's `[layers.<name>.cooling]` table. See [Layer](../Structures/layers/layer.md). A layer without one holds one temperature, as with `off`.

The models act in a world's thermal solve, `solve_eos(solve_temperature=True, surface_temperature=...)`. On each pass the thermal network asks every layer's model for its profile against the structure the pass solved: which stretches conduct and which are adiabatic, the resistance of each conducting stretch, and the temperature at the base of a convecting interior. The model reads its layer's material through the network, so it chooses where to evaluate it: the convection model takes its viscosity at the top of its adiabatic interior, where the layer's temperature applies. The network then joins the layers into a chain of thermal resistances, integrates the temperature and heat flow through the planet, and reports each layer's rate of temperature change. See [Temperature and Heat Flow](../Structures/worlds/worlds.md#temperature-and-heat-flow).

The models are parameterized. They reduce the whole of mantle convection to a boundary-layer scaling with a Rayleigh number and a handful of fitted constants, which is the standard approach for the timescales planetary evolution deals in. When their assumptions fail, usually in a thin layer, at a vanishing temperature drop, or with a viscosity the material cannot give, the models degrade to a defined edge case rather than to a numerical error, and the solve reports it (`layer_boundary_fallback`).

## Examples

`Demos/Physics/15_thermal_interior.ipynb` compares the models' flux laws and solves the temperature and heat flow profile of a layered world, and `Demos/Systems/12_thermal_orbital_evolution.ipynb` evolves a convecting, melting mantle.

## References

- Turcotte, D. L., and Schubert, G. (2002). *Geodynamics*, second edition. Rayleigh and Nusselt convection scaling, and conduction.
- Solomatov, V. S. (1995). Scaling of temperature- and stress-dependent viscosity convection. *Physics of Fluids*, 7(2), 266-274.
- Schubert, G., Turcotte, D. L., and Olson, P. (2001). *Mantle Convection in the Earth and Planets*. Boundary-layer theory.
- Stevenson, D. J., Spohn, T., and Schubert, G. (1983). Magnetism and thermal evolution of the terrestrial planets. *Icarus*, 54(3), 466-489. The upper-mantle temperature of parameterized convection.
- Solomatov, V. S. (2000). Fluid dynamics of a terrestrial magma ocean. In *Origin of the Earth and Moon*, University of Arizona Press, 323-338. Heat transport in a liquid interior.
