# Spin and Orbital Rates (`dynamics_x`)

_Updated: 2026-09-12_

The module holds two calculators. `Spin` answers how fast a body's rotation is changing; `OrbitSolver` answers how fast its orbit is changing. Both take the tidal-potential derivatives from a completed tidal solve and return rates in MKS units. Neither integrates in time.

The inputs both classes depend on are $\partial U / \partial M$, $\partial U / \partial \omega$, and $\partial U / \partial \Omega$, the derivatives of the tidal potential with respect to the mean anomaly, the argument of pericenter, and the longitude of the node, each in J kg$^{-1}$ rad$^{-1}$. They come from `world.calc_tides(...)` and are read back with `world.get_tidal_potential_derivatives()`. Everything on this page assumes that solve has already run.

## Spin

A tidal torque changes a body's rotation rate at $\dot{\Omega}_{\text{spin}} = M_{\text{host}} \, (\partial U / \partial \Omega) / I$, where $I$ is the body's polar moment of inertia. The sign follows from the potential derivative: a body spinning faster than its orbital mean motion is despun, and one spinning slower is spun up, until the two match and the torque vanishes. That fixed point is synchronous rotation, and `calc_synchronous_spin` simply returns the orbital mean motion, since a tidally locked body rotates once per orbit.

```python
from TidalPy.dynamics_x import Spin

model = Spin()                                               # factor 0.4, a uniform sphere
moment = model.calc_moment_of_inertia(mass, radius)          # [kg m2]
spin_rate = model.calc_dspin_dt(host_mass, dU_dO, moment)    # [rad s-2]
synchronous = model.calc_synchronous_spin(orbital_frequency) # [rad s-1]
```

### The moment of inertia

Everything about the spin rate except the moment of inertia comes from the tidal solve, and the moment of inertia is where a body's internal structure enters. `Spin` carries only a simple estimate of it:

$$I = f M R^2$$

where $f$ is the constructor's `moment_of_inertia_factor`, the conventional dimensionless factor $C / (M R^2)$. A uniform sphere has $f = 0.4$, which is the default. A centrally condensed body has less, with the Earth at 0.3307, and no body with non-negative density can exceed $2/3$, the value for all of its mass in a thin surface shell. A factor that is not finite or lies outside $(0, 2/3]$ raises `ValueError`. That bound deliberately rejects 1.0, so code written against the older ratio-to-a-uniform-sphere convention fails loudly instead of producing a moment of inertia 2.5 times too large.

This estimate exists as a fallback. A `LayeredWorld` that has solved its equation of state has the real structure-resolved moment of inertia, and `world.get_moment_of_inertia()` returns that instead, falling back to the model's formula only when no solve has run.

```python
world.set_spin_model(Spin())
world.solve_eos()
world.calc_tides(orbital_frequency, spin_frequency, eccentricity, obliquity,
                 semi_major_axis, host_mass)

world.get_moment_of_inertia()                    # [kg m2] EOS value once solved
world.calc_spin_derivative(host_mass)            # [rad s-2] uses the stored dU/dO
world.calc_synchronous_spin(orbital_frequency)   # [rad s-1]
world.planet_moi_eos                             # [kg m2] the raw EOS value, NaN if unsolved
```

`set_spin_model` copies the model into the world, so the Python object stays usable afterward. `calc_spin_derivative` requires a completed tidal solve, because it reads the world's stored $\partial U / \partial \Omega$.

## OrbitSolver

The orbital rates follow from Lagrange's planetary equations applied to the tidal part of the disturbing function. Following Boué and Efroimsky (2019), equations 116 and 117, with

$$\frac{\partial R}{\partial X} = -\frac{M_{\text{target}} + M_{\text{host}}}{M_{\text{target}}} \frac{\partial U}{\partial X}$$

the three rates are

$$\dot{a} = \frac{2}{n a} \frac{\partial R}{\partial M}, \qquad \dot{e} = \frac{\sqrt{1 - e^2}}{n a^2 e} \left( \sqrt{1 - e^2} \frac{\partial R}{\partial M} - \frac{\partial R}{\partial \omega} \right), \qquad \dot{n} = -\frac{3}{2} \frac{n}{a} \dot{a}$$

The last is Kepler's third law differentiated, so the mean motion is not an independent quantity: it follows from the semi-major axis.

```python
from TidalPy.dynamics_x import OrbitSolver

solver = OrbitSolver()
da_dt = solver.calc_da_dt(orbital_frequency, semi_major_axis, eccentricity,
                          target_mass, host_mass, dU_dM)                     # [m s-1]
de_dt = solver.calc_de_dt(orbital_frequency, semi_major_axis, eccentricity,
                          target_mass, host_mass, dU_dM, dU_dw)              # [s-1]
dn_dt = solver.calc_dn_dt(orbital_frequency, semi_major_axis, da_dt)         # [rad s-2]

rates = solver.calc_derivatives(orbital_frequency, semi_major_axis, eccentricity,
                                target_mass, host_mass, dU_dM, dU_dw)
rates["da_dt"], rates["de_dt"], rates["dn_dt"]
```

`calc_derivatives` computes all three in one call and is the entry point the system class uses. The eccentricity rate carries a $1/e$ factor that is indeterminate at $e = 0$, so a circular or degenerate orbit returns exactly zero rather than a division by zero.

For a system where both bodies dissipate, the two contributions are additive in the disturbing-function derivatives: solve each body's tides with the other as the raiser and sum the rates. `System.calc_pair_evolution` does exactly that.

## Energy conservation

The spin and orbital rates are not independent. Together they must account for all the energy the tidal solve says is being dissipated:

$$Q_{\text{tidal}} = -\left( \dot{E}_{\text{orbit}} + \dot{E}_{\text{spin}} \right), \qquad E_{\text{orbit}} = -\frac{G M_{\text{target}} M_{\text{host}}}{2 a}, \qquad E_{\text{spin}} = \frac{1}{2} I \Omega_{\text{spin}}^2$$

This is the strongest available check on the whole chain, because it ties the rate equations here back to the heating computed by an entirely separate code path. Every evolution dict returned by `System` reports `dE_orbit_dt`, `dE_spin_dt`, and the `energy_residual` between them and the tidal heating, so a failure anywhere upstream shows up as a non-zero residual. The balance is verified to machine precision in `Tests/Test_Structures_x/Test_Worlds/test_world_spin_01.py`.

## Driving the rates from a system

Calling the two classes by hand is useful for a single body in isolation. For anything with more than one world, `System` holds the orbital state and does the bookkeeping.

```python
rates = system.calc_world_evolution(world)   # one dissipating world, rigid host
rates = system.calc_pair_evolution(world)    # both bodies dissipate on a shared orbit
all_rates = system.calc_system_evolution()   # every world, in index order
```

`calc_world_evolution` returns the orbital state it used, the tidal heating, the three potential derivatives, the four rates, the moment of inertia, and the energy-balance diagnostics, all in one dict. `calc_pair_evolution` returns the combined orbital rates plus each body's full single-body contribution under the keys `world` and `host`. Entries that could not be evolved, such as the host's own row or a world with no usable orbit, come back with `evolved` set to `False` rather than raising. See [System](../structures_x/system/system.md).

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
- `c_OrbitSolver` with the `c_OrbitState` input struct and the `c_OrbitDerivatives` result struct: `calc_da_dt`, `calc_de_dt`, `calc_dn_dt`, `calc_derivatives`.

Both classes are small, stateless value types with no heap allocation and no base class. The orbital state is passed as a struct rather than as five loose arguments, and the derivatives come back as a struct rather than through output parameters.

## References

- Boué, G., and Efroimsky, M. (2019). Tidal evolution of the Keplerian elements. *Celestial Mechanics and Dynamical Astronomy*, 131(7), 30. Equations 116 and 117, the rate equations implemented here.
- Ferraz-Mello, S., Rodríguez, A., and Hussmann, H. (2008). Tidal friction in close-in satellites and exoplanets. *Celestial Mechanics and Dynamical Astronomy*, 101(1-2), 171-201. The spin-rate equation.
- Murray, C. D., and Dermott, S. F. (1999). *Solar System Dynamics*. Lagrange's planetary equations and the disturbing function.
