# Spin and Orbital Rates (`Dynamics`)

_Updated: 2026-10-09_

`Spin` gives the rate of change of a body's rotation and `OrbitSolver` that of its orbit. Both take the derivatives of the tidal potential with respect to the mean anomaly, the argument of pericenter, and the longitude of the node, $\partial U / \partial M$, $\partial U / \partial \omega$, and $\partial U / \partial \Omega$ \[J kg$^{-1}$ rad$^{-1}$\], from a completed tidal solve: `world.calc_tides(...)` returns them as `dU_dM`, `dU_dw`, and `dU_dO`, and `world.get_tidal_potential_derivatives()` reads them back.

## Spin

A tidal torque changes a body's rotation rate at $\dot{\Omega}_{\text{spin}} = M_{\text{host}} \, (\partial U / \partial \Omega) / C$, where $C$ is the body's polar moment of inertia. _Usually_ a body spinning faster than its mean motion is despun, and a slower one spun up, until synchronous rotation, where the torque vanishes (`calc_synchronous_spin` returns the mean motion). With many tidal modes, depending on the thermal state and the spin-orbit ratio, the spin may not follow these trends (see Renaud et al. 2021).

```python
from TidalPy.Dynamics import Spin

model = Spin()                                               # factor 0.4, a uniform sphere
moment = model.calc_moment_of_inertia(mass, radius)          # [kg m2]
spin_rate = model.calc_dspin_dt(host_mass, dU_dO, moment)    # [rad s-2]
synchronous = model.calc_synchronous_spin(orbital_frequency) # [rad s-1]
```

### Moment of Inertia

`Spin` carries a simple estimate of the moment of inertia,

$$C = f M R^2$$

where $f$ is the constructor's `moment_of_inertia_factor`, $C / (M R^2)$. A uniform sphere has $f = 0.4$, the default. A centrally condensed body has less (the Earth about 0.331), and no body with non-negative density can exceed $2/3$, the value for a thin surface shell. A factor that is not finite or lies outside $(0, 2/3]$ raises `ValueError`.

This is a fallback: after an equation-of-state solve, `world.get_moment_of_inertia()` returns the structure-resolved value.

```python
world.set_spin_model(Spin())  # Has a fallback MOI
world.solve_eos()             # Now has a more accurate MOI
world.calc_tides(orbital_frequency, spin_frequency, eccentricity, obliquity,
                 semi_major_axis, host_mass)

world.get_moment_of_inertia()                    # [kg m2] EOS value once solved
world.calc_spin_derivative(host_mass)            # [rad s-2] uses the stored dU/dO
world.calc_synchronous_spin(orbital_frequency)   # [rad s-1]
world.planet_moi_eos                             # [kg m2] the raw EOS value, NaN if unsolved
```

`set_spin_model` copies the model into the world, so the Python object stays usable. `calc_spin_derivative` requires a completed tidal solve, since it reads the stored $\partial U / \partial \Omega$.

## OrbitSolver

The orbital rates follow from Lagrange's planetary equations applied to the tidal part of the disturbing function (Boué and Efroimsky 2019, Eqs. 116 and 117). With

$$\frac{\partial R}{\partial X} = -\frac{M_{\text{target}} + M_{\text{host}}}{M_{\text{target}}} \frac{\partial U}{\partial X}$$

the three rates are

$$\dot{a} = \frac{2}{n a} \frac{\partial R}{\partial M}, \qquad \dot{e} = \frac{\sqrt{1 - e^2}}{n a^2 e} \left( \sqrt{1 - e^2} \frac{\partial R}{\partial M} - \frac{\partial R}{\partial \omega} \right), \qquad \dot{n} = -\frac{3}{2} \frac{n}{a} \dot{a}$$

The last is Kepler's third law differentiated: the mean motion follows from the semi-major axis.

```python
from TidalPy.Dynamics import OrbitSolver

solver = OrbitSolver()
da_dt = solver.calc_da_dt(orbital_frequency, semi_major_axis, eccentricity, target_mass, host_mass, dU_dM)   # [m s-1]
de_dt = solver.calc_de_dt(orbital_frequency, semi_major_axis, eccentricity, target_mass, host_mass, dU_dM, dU_dw)   # [s-1]
dn_dt = solver.calc_dn_dt(orbital_frequency, semi_major_axis, da_dt)   # [rad s-2]
rates = solver.calc_derivatives(orbital_frequency, semi_major_axis, eccentricity, target_mass, host_mass, dU_dM, dU_dw)

print(rates["da_dt"], rates["de_dt"], rates["dn_dt"])
```

`calc_derivatives` computes all three in one call, as the system class does. The $1/e$ factor of the eccentricity rate is indeterminate at $e = 0$, so a circular or degenerate orbit returns exactly zero rather than raise a divide-by-zero error.

At small eccentricity the two terms of $\dot{e}$ nearly cancel, and forming their difference from the separate sums loses precision. `calc_de_dt` and `calc_derivatives` take the per-mode sum of $\partial U / \partial M - \partial U / \partial \omega$ as an optional last argument, `dU_dM_minus_dw`, which keeps the rate exact. A world returns it from `world.get_tidal_dU_dM_minus_dw()` after `calc_tides`, and `System` passes it for you.

When both bodies dissipate, their contributions add in the disturbing-function derivatives; `System.calc_pair_evolution` combines the tide raised on the host with the one raised on the target world.

## Energy Conservation

Together the spin and orbital rates must account for all the energy dissipated by tides,

$$Q_{\text{tidal}} = -\left( \dot{E}_{\text{orbit}} + \dot{E}_{\text{spin}} \right), \qquad E_{\text{orbit}} = -\frac{G M_{\text{target}} M_{\text{host}}}{2 a}, \qquad E_{\text{spin}} = \frac{1}{2} C \Omega_{\text{spin}}^2$$

Every evolution dict from `System` reports `dE_orbit_dt`, `dE_spin_dt`, and the `energy_residual` between them and the tidal heating, so a failure anywhere upstream shows as a non-zero residual.

## Driving the Rates from a `System`

For more than one world, `System` holds the orbital state and does the bookkeeping.

```python
tide = system.calc_dissipation(world)        # one world's tide on its orbit about its tidal host, no rates
rates = system.calc_world_evolution(world)   # one dissipating world, rigid host
rates = system.calc_pair_evolution(world)    # both bodies (the `world` and its tidal host, or a given partner) dissipate on a shared orbit
all_rates = system.calc_system_evolution()   # every world, in index order
```

`calc_world_evolution` returns one dict with the orbital state used, the tidal heating, the potential derivatives, the four rates, the moment of inertia, and the energy-balance diagnostics. `calc_pair_evolution` adds each body's own contribution under `worlds`, keyed by world name. An entry that cannot be evolved (the host's own row, a world with no usable orbit) has `evolved` set to `False` rather than raising. See [System](../Structures/system/system.md). To integrate the rates of a pair in time, through capture into and release from spin-orbit locks, use `System.evolve` (see [Evolving a World About Its Host](../Structures/system/system.md#evolving-a-world-about-its-host)).

## C++ API

```cpp
#include "spin_.hpp"
#include "orbit_solver_.hpp"

using namespace tidalpy;

c_SpinConfig spin_config;
spin_config.moment_of_inertia_factor = 0.3307;
const c_Spin spin(spin_config);
const double dspin_dt = spin.calc_dspin_dt(host_mass, dU_dO, moment_of_inertia);

c_OrbitState state;   // orbital_frequency, semi_major_axis, eccentricity, target_mass, host_mass
const c_OrbitSolver solver;
const c_OrbitDerivatives rates = solver.calc_derivatives(state, dU_dM, dU_dw);
```

- `c_Spin` with `c_SpinConfig`: `calc_moment_of_inertia`, `calc_dspin_dt`, `calc_synchronous_spin`. The config constructor throws `std::invalid_argument` for a factor outside $(0, 2/3]$.
- `c_OrbitSolver`, with the `c_OrbitState` input and `c_OrbitDerivatives` result structs: `calc_da_dt`, `calc_de_dt`, `calc_dn_dt`, `calc_derivatives`.

Both are stateless value types.

## References

- Boué, G., and Efroimsky, M. (2019). Tidal evolution of the Keplerian elements. *Celestial Mechanics and Dynamical Astronomy*, 131(7), 30. Equations 116 and 117, the rate equations implemented here.
- Ferraz-Mello, S., Rodríguez, A., and Hussmann, H. (2008). Tidal friction in close-in satellites and exoplanets. *Celestial Mechanics and Dynamical Astronomy*, 101(1-2), 171-201. The spin-rate equation.
- Murray, C. D., and Dermott, S. F. (1999). *Solar System Dynamics*. Lagrange's planetary equations and the disturbing function.
