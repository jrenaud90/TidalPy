# Love Numbers and Radial Functions (`RadialSolver`)

_Updated: 2026-10-07_

`TidalPy.RadialSolver` solves the viscoelastic-gravitational problem for a layered, spherically symmetric planet. It returns the radial functions $y_1$ through $y_6$ throughout the interior and the Love numbers $k$, $h$, and $l$ at the surface. These set the magnitude of tidal dissipation, the speed of orbital and rotational evolution, and a body's predicted gravity and displacement signals.

| Entry point | Use it when |
|---|---|
| `BaseWorld.solve_love_numbers(...)` | You have a built world. Its equation of state is solved, its solid and liquid zones are the solver's layers, and each layer's rheology supplies the complex moduli. See [Worlds](../Structures/worlds/worlds.md). |
| `radial_solver(...)` | You have arrays rather than a world, or want to drive the solver directly. See [Calculating Love Numbers](calculating_love_numbers.md); the [helper functions](build_inputs.md) build its inputs. |
| `homogeneous_love_numbers(...)` | You want a quick estimate for a uniform sphere. Not realistic for most worlds. |

Every path returns a [`RadialSolverSolution`](solution_class.md) with the radial functions, the Love numbers, the equation-of-state profiles, and diagnostics.

```{toctree}
:maxdepth: 2
:caption: Contents

Calculating Love Numbers <calculating_love_numbers.md>
Starting Conditions <starting_conditions.md>
Solution Class <solution_class.md>
Helper Functions <build_inputs.md>
Dense Radial Solutions <dense_radial_solution.md>
Parallel Love Solves <parallel.md>
```

## Solution Methods

The default "shooting method" integrates the viscoelastic-gravitational equations from a starting radius near the center out to the surface, one layer at a time, then fixes the combination of each layer's independent solutions from the surface boundary condition (see [Dense Radial Solutions](dense_radial_solution.md)). [Starting Conditions](starting_conditions.md) covers how it starts. The quasi-analytic propagation matrix is restricted to a single solid, static, incompressible layer. Three analytic methods, `homogeneous`, `cpl`, and `ctl`, skip the interior solve. [Calculating Love Numbers](calculating_love_numbers.md#choosing-a-method) compares them.

## References

A starting point, not a comprehensive list. The starting-condition references are on [Starting Conditions](starting_conditions.md#references).

**Numerical shooting method**
- Takeuchi, H., and Saito, M. (1972). Seismic Surface Waves. In *Methods in Computational Physics: Advances in Research and Applications*, 11, 217-295.
- Tobie, G., Mocquet, A., and Sotin, C. (2005). Tidal dissipation within large icy satellites: Applications to Europa and Titan. *Icarus*, 177(2), 534-549.
- Kervazo, M., Tobie, G., Choblet, G., Dumoulin, C., and Běhounková, M. (2021). Solid tides in Io's partially molten interior. *Astronomy and Astrophysics*, 650, A72.

**Interfaces, constants, and assumptions**
- Saito, M. (1974). Some problems of static deformation of the earth. *Journal of Physics of the Earth*, 22(1), 123-140.
- Beuthe, M. (2015). Tidal Love numbers of membrane worlds: Europa, Titan, and Co. *Icarus*, 258, 239-266.

**Propagation matrix method**
- Sabadini, R., and Vermeersen, B. (2004). *Global Dynamics of the Earth: Applications of Normal Mode Relaxation Theory to Solid-Earth Geophysics*.
- Roberts, J. H., and Nimmo, F. (2008). Tidal heating and the long-term stability of a subsurface ocean on Enceladus. *Icarus*, 194(2), 675-689.
- Henning, W. G., and Hurford, T. (2014). Tidal heating in multilayered terrestrial exoplanets. *The Astrophysical Journal*, 789(1), 30.
- Sabadini, R., Vermeersen, B., and Cambiotti, G. (2016). *Global Dynamics of the Earth*, second edition.

**Homogeneous sphere**
- Love, A. E. H. (1911). *Some Problems of Geodynamics*.
- Munk, W. H., and MacDonald, G. J. F. (1960). *The Rotation of the Earth: A Geophysical Discussion*.

## Examples

`Demos/Physics/P05_love_numbers_1d.ipynb` calls the solver directly; the world notebooks reach it through `solve_love_numbers`. The notebooks in `Benchmarks/RadialSolver/` validate the results against published Earth and Enceladus models and the closed-form homogeneous sphere.
