# Dynamics (`dynamics_x`)

_Updated: 2026-09-12_

`TidalPy.dynamics_x` is where tidal dissipation stops being a heating rate and becomes an orbit. It takes the tidal-potential derivatives produced by a tidal solve and turns them into the instantaneous rates of change of a body's spin and of the orbit it sits in.

The physical statement behind this is simple to write and easy to lose track of. Tides raised on a body lag behind the raising potential, and that lag exerts a torque. The torque exchanges angular momentum between the body's rotation and the orbit, while the associated friction removes energy from the pair and deposits it as heat. Everything in this module follows from that exchange: a satellite spins down to synchronous rotation, its orbit expands or decays, and its eccentricity is damped or pumped.

| Page | Covers |
|---|---|
| [Spin and Orbital Rates](dynamics.md) | The `Spin` and `OrbitSolver` classes, the rate equations they implement, the moment of inertia they use, energy conservation, and how the world and system classes drive them. |

```{toctree}
:maxdepth: 1

Spin and Orbital Rates <dynamics.md>
```

## Rates, not trajectories

This module computes rates only. It never integrates in time. That separation is deliberate: a tidal evolution over billions of years is an integration problem whose stiffness, step control, and stopping conditions depend entirely on the study, and hard-coding one choice would make the rest impossible. Getting $\dot{a}$, $\dot{e}$, $\dot{n}$, and $\dot{\Omega}_{\text{spin}}$ at a given state is the part that is genuinely TidalPy's, and it is what these classes provide.

To evolve a system, hand those rates to an integrator of your choice. `System` supplies the state at each step and collects the rates.

## Where dynamics fits

The chain runs in one direction. A world's [rheology](../rheology_x/index.md) gives its complex moduli; the [radial solver](../RadialSolver_x/index.md) turns those into Love numbers; the [tidal solve](../Tides_x/index.md) collapses the Kaula mode sum into a heating rate and the three potential derivatives $\partial U / \partial M$, $\partial U / \partial \omega$, and $\partial U / \partial \Omega$ with respect to the mean anomaly, the argument of pericenter, and the longitude of the node. This module consumes those three numbers.

A `LayeredWorld` holds a `Spin` model and drives it with its own moment of inertia from the equation-of-state solve, so the spin rate uses the structure-resolved value rather than a uniform-density estimate. See [Worlds](../structures_x/worlds/worlds.md). A [System](../structures_x/system/system.md) attaches an `OrbitSolver`, pulls the orbital state and potential derivatives from its worlds, and reports the full set of rates with an energy-balance diagnostic.

## References

- Boué, G., and Efroimsky, M. (2019). Tidal evolution of the Keplerian elements. *Celestial Mechanics and Dynamical Astronomy*, 131(7), 30. The disturbing-function form of the orbital rate equations used here.
- Ferraz-Mello, S., Rodríguez, A., and Hussmann, H. (2008). Tidal friction in close-in satellites and exoplanets. *Celestial Mechanics and Dynamical Astronomy*, 101(1-2), 171-201. The spin-rate equation.
- Murray, C. D., and Dermott, S. F. (1999). *Solar System Dynamics*. The disturbing function and Lagrange's planetary equations.
