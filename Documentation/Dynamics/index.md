# Dynamics (`Dynamics`)

_Updated: 2026-10-05_

`TidalPy.Dynamics` contains functionality to couple tidal dissipation with orbit-spin evolution. It takes the tidal-potential derivatives produced by a tidal solve and turns them into the instantaneous rates of change of a body's spin and orbit.

Tides raised on a viscoelastic body lag behind the raising potential, and that lag exerts a torque. The torque exchanges angular momentum between the body's rotation and the orbit, while the associated frictional heat removes energy from both bodies.

This module just computes derivative rates, it does not perform integrations directly. To evolve a system, use these rates in an integrator of your choice. `System` supplies the state at each step and collects the rates, which is the usual starting point for a thermal-orbital evolution model.

| Page | Covers |
|---|---|
| [Spin and Orbital Rates](dynamics.md) | The `Spin` and `OrbitSolver` classes, the rate equations they implement, the moment of inertia they use, energy conservation, and how the world and system classes drive them. |

```{toctree}
:maxdepth: 1

Spin and Orbital Rates <dynamics.md>
```

## Where Dynamics is Used

The solve runs in one direction. A world's [EOS](../Material/index.md) provides static moduli, interior structure, etc. [Rheology](../Rheology/index.md) converts these to complex moduli. The [radial solver](../RadialSolver/index.md) turns those into Love numbers. The [tidal solve](../Tides/index.md) collapses the Kaula mode sum into a heating rate and the three potential derivatives $\partial U / \partial M$, $\partial U / \partial \omega$, and $\partial U / \partial \Omega$ with respect to the mean anomaly, the argument of pericenter, and the longitude of the node. This module takes all of those to find the change in spin and orbital elements.

Every world holds a `Spin` model and drives it with its own moment of inertia from the equation-of-state solve, so the spin rate uses the structure-resolved value rather than a uniform-density estimate. See [Worlds](../Structures/worlds/worlds.md). A [System](../Structures/system/system.md) attaches an `OrbitSolver`, pulls the orbital state and potential derivatives from its worlds, and reports the full set of rates with an energy-balance diagnostic.

## Examples

`Demos/Physics/07_gasgiant_fixedQ_dt.ipynb` reads the orbital rates of a gas giant from `System.calc_world_evolution`, `Demos/Systems/11_multi_world.ipynb` reads them for every world in a system, and `Demos/Systems/12_thermal_orbital_evolution.ipynb` integrates them in time with `System.evolve`.

## References

- Boué, G., and Efroimsky, M. (2019). Tidal evolution of the Keplerian elements. *Celestial Mechanics and Dynamical Astronomy*, 131(7), 30. The disturbing-function form of the orbital rate equations used here.
- Ferraz-Mello, S., Rodríguez, A., and Hussmann, H. (2008). Tidal friction in close-in satellites and exoplanets. *Celestial Mechanics and Dynamical Astronomy*, 101(1-2), 171-201. The spin-rate equation.
- Murray, C. D., and Dermott, S. F. (1999). *Solar System Dynamics*. The disturbing function and Lagrange's planetary equations.
