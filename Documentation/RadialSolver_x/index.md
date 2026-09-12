# RadialSolver_x: Love Numbers and Radial Functions

_Updated: 2026-09-12_

`TidalPy.RadialSolver_x` solves the viscoelastic-gravitational problem for a layered, spherically symmetric planet. It returns the radial functions y1 through y6 throughout the interior and the Love numbers k, h, and l at the surface. Those numbers set the magnitude of tidal dissipation, the speed of orbital and rotational evolution, and the predicted gravity and displacement signals of a body.

There are three ways in, in increasing order of control:

| Entry point | Use it when |
|---|---|
| `LayeredWorld.solve_love_numbers(...)` | You have a built world. The layer rheologies supply the complex moduli and the equation of state is already solved. See [Worlds](../structures_x/worlds/worlds.md). |
| `radial_solver(...)` | You have arrays rather than a world, or you want to drive the solver directly. Documented in [Calculating Love Numbers](calculating_love_numbers.md). |
| `homogeneous_love_numbers(...)` | You want a quick estimate for a uniform sphere without building anything. |

Every path returns a [`RadialSolverSolution`](solution_class.md), which carries the radial functions, the Love numbers, the equation-of-state profiles, and the diagnostics that tell you whether to trust them.

```{toctree}
:maxdepth: 2
:caption: Contents

Calculating Love Numbers <calculating_love_numbers.md>
Solution Class <solution_class.md>
Helper Functions <build_inputs.md>
Dense Radial Solutions <dense_radial_solution.md>
```

## What the solver does

The default method integrates a set of viscoelastic-gravitational ordinary differential equations from a starting radius near the center out to the surface, one layer at a time. Each layer contributes a fixed number of independent solutions depending on its assumptions, and the physical solution in the layer is a linear combination of them. The combination coefficients are fixed by the surface boundary condition and then propagated downward through every interface, so each layer inherits a consistent set of constants. The propagation-matrix alternative is quasi-analytic and restricted to a single solid, static, incompressible layer; the analytic `homogeneous`, `cpl`, and `ctl` methods skip the interior solve altogether. [Calculating Love Numbers](calculating_love_numbers.md) explains how to choose between them, and [Dense Radial Solutions](dense_radial_solution.md) covers the numerics in more depth.

## References

The methods implemented here come from the following work.

**Numerical shooting method**
- Takeuchi, H., and Saito, M. (1972). Seismic Surface Waves. In *Methods in Computational Physics: Advances in Research and Applications*, 11, 217-295.
- Tobie, G., Mocquet, A., and Sotin, C. (2005). Tidal dissipation within large icy satellites: Applications to Europa and Titan. *Icarus*, 177(2), 534-549.
- Kervazo, M., Tobie, G., Choblet, G., Dumoulin, C., and Běhounková, M. (2021). Solid tides in Io's partially molten interior. *Astronomy and Astrophysics*, 650, A72.

**Interfaces, constants, and assumptions**
- Saito, M. (1974). Some problems of static deformation of the earth. *Journal of Physics of the Earth*, 22(1), 123-140.
- Beuthe, M. (2015). Tidal Love numbers of membrane worlds: Europa, Titan, and Co. *Icarus*, 258, 239-266.

**Starting conditions**
- Kamata, S., Matsuyama, I., and Nimmo, F. (2015). Tidal resonance in icy satellites with subsurface oceans. *Journal of Geophysical Research: Planets*, 120(9), 1528-1542.
- Martens, H. R. (2016). *Using Earth deformation caused by surface mass loading to constrain the elastic structure of the crust and mantle*. PhD thesis, California Institute of Technology.

**Propagation matrix method**
- Sabadini, R., and Vermeersen, B. (2004). *Global Dynamics of the Earth: Applications of Normal Mode Relaxation Theory to Solid-Earth Geophysics*.
- Roberts, J. H., and Nimmo, F. (2008). Tidal heating and the long-term stability of a subsurface ocean on Enceladus. *Icarus*, 194(2), 675-689.
- Henning, W. G., and Hurford, T. (2014). Tidal heating in multilayered terrestrial exoplanets. *The Astrophysical Journal*, 789(1), 30.
- Sabadini, R., Vermeersen, B., and Cambiotti, G. (2016). *Global Dynamics of the Earth*, second edition.

**Homogeneous sphere**
- Love, A. E. H. (1911). *Some Problems of Geodynamics*.
- Munk, W. H., and MacDonald, G. J. F. (1960). *The Rotation of the Earth: A Geophysical Discussion*.

## Learning by example

The `Demos (_x)` notebooks work through the solver from both ends: `Physics/08_love_numbers_1d.ipynb` drives it directly, and the world notebooks reach it through `solve_love_numbers`. The `Benchmarks (_x)` pages validate the results against published Earth and Enceladus models.
