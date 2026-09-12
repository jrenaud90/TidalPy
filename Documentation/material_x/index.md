# Material and Equation of State (`Material_x`)

_Updated: 2026-09-12_

`TidalPy.Material_x` describes what a planet is made of, in the one sense a tidal calculation needs: how dense the material is at the pressure it finds itself under. Each model maps the local state onto a mass density [kg m$^{-3}$], and the whole-planet solve integrates those densities from the center outward to produce the body's radial structure.

That structure is not a detail of the setup. It is the coefficient set of every later calculation. The density profile fixes the gravity and pressure profiles, which fix the moment of inertia, which is one of the few interior quantities a spacecraft can actually measure. It also fixes the coefficients of the radial functions the Love-number solver integrates, so two bodies with the same mass and radius but different internal density distributions have measurably different $k_2$. Getting the equation of state wrong does not shift a tidal answer slightly; it changes the interior the answer describes.

| Page | Covers |
|---|---|
| [Material EOS Models](material_eos.md) | The four models, the pressure inversion they share, the interpolated tables, the factory, serialization, the C++ surface, and how to add a model. |

```{toctree}
:maxdepth: 1

Material EOS Models <material_eos.md>
```

## Where the equation of state fits

An EOS model is attached to a layer with `BaseLayer.set_eos`. Once every layer has one, `LayeredWorld.solve_eos()` integrates the planet's radial structure from the center to the surface and populates each layer's density, gravity, pressure, mass, and moment-of-inertia profiles. That solve, its convergence loop, and its results are documented with the world class; see [Worlds](../structures_x/worlds/worlds.md).

The integration runs over radius and carries pressure as one of its state variables, so an analytic density law that depends on pressure is simply evaluated at each step with the pressure the integrator has already reached. No separate coupled iteration is needed beyond the solver's outer loop, which adjusts the central pressure until the integrated surface pressure matches the requested boundary value.

A world built from a TOML file gets its EOS models from the `[layers.<name>.eos]` table, with material defaults filling in anything the user omits. See the [TOML schema](../structures_x/config/toml_schema.md).

## Scope

The analytic models are isothermal. They accept a temperature for interface uniformity and currently ignore it, so thermal expansion is not part of the density they return. Phase transitions, composition gradients within a layer, and self-consistent thermal structure are all outside the module. The interpolated model is the escape hatch: a density profile computed by any external tool, including a full mineral-physics package, can be loaded as a table and used exactly like an analytic law.

## References

- Birch, F. (1947). Finite elastic strain of cubic crystals. *Physical Review*, 71(11), 809-824. The finite-strain equation of state.
- Vinet, P., Ferrante, J., Rose, J. H., and Smith, J. R. (1987). Compressibility of solids. *Journal of Geophysical Research*, 92(B9), 9319-9325. The universal equation of state.
- Dziewonski, A. M., and Anderson, D. L. (1981). Preliminary reference Earth model. *Physics of the Earth and Planetary Interiors*, 25(4), 297-356. The tabulated profile shipped with TidalPy.
