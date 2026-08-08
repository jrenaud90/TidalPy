# Dynamics (`dynamics_x`)

The `dynamics_x` module holds the spin and orbital rate equations that turn a tidal solve into the
instantaneous evolution of a body's rotation and orbit. It computes **rates only** (no time
integration; the System class or manual scripts can perform integrations). All quantities are MKS.

The tidal solve (`world.calc_tides`) produces, per the global mode collapse, the tidal heating and the
three tidal-potential derivatives `dU/dM`, `dU/dw`, `dU/dO` (with respect to the mean anomaly, argument
of pericenter, and longitude of the node) [J kg-1 rad-1]. These feed the dynamics rates below.

## Spin

`Spin` is the spin-dynamics model. A world holds one and drives it with _its own EOS-based moment of
inertia_, so the spin-rate change uses the structure-resolved value rather than the uniform-density
formula:

```python
from TidalPy.dynamics_x import Spin

world.set_spin_model(Spin(moment_of_inertia_factor=1.0))   # factor is the pre-EOS fallback
world.solve_eos()
world.calc_tides(orbital_frequency, spin_frequency, eccentricity, obliquity, semi_major_axis, host_mass)

moi = world.get_moment_of_inertia()                       # [kg m2] EOS value once solved, else uniform fallback
dspin_dt = world.calc_spin_derivative(host_mass)          # [rad s-2]
n_sync = world.calc_synchronous_spin(orbital_frequency)   # [rad s-1]
```

* **Moment of inertia** — `world.get_moment_of_inertia()` returns the EOS-integrated moment of inertia
  (`planet_moi_eos`) once the EOS has been solved, otherwise the model's uniform-density fallback
  `factor * (2/5) M (R_o^5 - R_i^5)/(R_o^3 - R_i^3)`. The dimensionless `moment_of_inertia_factor`
  (`C/(M R^2)`, `1.0` for a uniform sphere) is only used for that fallback.
* **Tidal spin-rate change** — `world.calc_spin_derivative(host_mass) = M_host * dU/dO / I`
  (Ferraz-Mello et al. 2008), using the world's stored `dU/dO` (from the last `calc_tides`) and its
  moment of inertia. Requires a completed tidal solve.
* **Synchronous spin** — a tidally locked body spins at the orbital mean motion, so the synchronous spin
  equals the orbital frequency.

The `Spin` model can also be used standalone (`Spin().calc_moment_of_inertia(...)`,
`.calc_dspin_dt(...)`, `.calc_synchronous_spin(...)`), which is what the world delegates to.

## OrbitSolver

`OrbitSolver` computes the orbital rates from a dissipating body's tidal-potential derivatives. It is
the low-level engine; the **System** class attaches it and pulls the orbital state and
`dU/dM`, `dU/dw` from its worlds. Following Boué & Efroimsky (2019, CMDA) Eqs. 116-117, with
`dR/dX = -((M_target + M_host)/M_target) dU/dX`:

```python
from TidalPy.dynamics_x import OrbitSolver

orbit = OrbitSolver()
da_dt = orbit.calc_da_dt(n, a, e, target_mass, host_mass, dU_dM)               # [m s-1]
de_dt = orbit.calc_de_dt(n, a, e, target_mass, host_mass, dU_dM, dU_dw)        # [s-1]
dn_dt = orbit.calc_dn_dt(n, a, da_dt)                                          # [rad s-2]
rates = orbit.calc_derivatives(n, a, e, target_mass, host_mass, dU_dM, dU_dw)  # dict
```

* `da/dt = (2 / (n a)) dR/dM`
* `de/dt = (sqrt(1-e^2) / (n a^2 e)) ( sqrt(1-e^2) dR/dM - dR/dw )` (0 for a circular orbit)
* `dn/dt = -(3/2)(n / a) da/dt` (Kepler's third law)

For a dual-dissipation system the two bodies' rates are additive in the disturbing-function
derivatives; the System sums the per-body contributions.

## Energy conservation

The world spin-rate and the orbital rates together conserve energy:
`tidal_heating = -(dE_orbit/dt + dE_spin/dt)` with `E_orbit = -G M_target M_host / (2a)` and
`E_spin = (1/2) I spin^2`. This is validated to machine precision in
`Tests/Test_Structures_x/Test_Worlds/test_world_spin_01.py`.
