# System (`structures_x.system`)

A `System` links two or more worlds (a star, planets, moons) into a gravitationally bound group. It
tracks two roles independently:

* the **tidal host** — the body that raises a tide on a target world (_e.g._, the Moon raising a tide on Earth or
vice versa).
* the **star** — star provides insolation heating on each of the worlds in the system.

The star need not be the tidal host. For the Earth-Moon system the tidal host can be either the Earth
or the Moon, depending on which world is the focus of your study. For an exoplanet orbiting its star
the star is acts as the tidal host, and all planets in the system orbit about the host and about the
star coincide. Each world therefore carries **two** two-body orbits: one about the tidal host and one
about the star. Worlds interact only with the host / the star, not with one another. The System is the
container on which orbital evolution and insolation are computed. All quantities are MKS.

## Building a system

```python
from TidalPy.structures_x.system import System
from TidalPy.structures_x.worlds.stellar import StarWorld
from TidalPy.structures_x.worlds.layered import LayeredWorld

AU = 1.495978707e11

system = System("sol")
star = StarWorld("sun", 6.957e8, 1.988e30)
earth = LayeredWorld("earth", 6.371e6, 5.972e24)

system.add_world(star, is_host=True, is_star=True)   # exoplanet: the star is also the tidal host
system.add_world(earth, semi_major_axis=AU, eccentricity=0.0167)
```

`add_world(world, is_host=False, is_star=False, semi_major_axis=None, eccentricity=0.0)` returns the
world's index. `is_host` and `is_star` are independent roles — a world can be both (the last world
flagged with each wins). The `semi_major_axis`/`eccentricity` describe the orbit about the **tidal
host**; the orbit about the star is set separately (see below). The system **co-owns** each world with
its Python wrapper, so the wrapper you passed stays fully usable and is the same object the system hands
back.

For a system where the star is a separate body from the tidal host (e.g. Earth-Moon-Sun):

```python
system = System()
system.add_world(sun, is_star=True)                                   # star, not the tidal host
system.add_world(earth, is_host=True)                                 # tidal host, not the star
system.add_world(moon, semi_major_axis=3.844e8, eccentricity=0.0549)  # moon about Earth (tidal)
system.set_stellar_semi_major_axis("moon", AU)                        # moon about the Sun (insolation)
system.set_stellar_eccentricity("moon", 0.0167)
```

## Worlds, Host, and Identification

Most methods accept a world by **index** (int, negatives allowed), **name** (str), or the **world
object** itself:

```python
system.host                       # the host world wrapper (or None)
system.host_index                 # index of the host, or -1
system.set_host("sun")            # designate the host by index / name / object
system.num_worlds                 # 2
system.worlds                     # [star, earth]
```

The system is a sequence over its worlds:

```python
for world in system: ...          # iterate members
len(system)                       # 2
system[0]                         # the star
system[1:]                        # [earth]
system.earth                      # attribute access by world name
system["earth"]                   # or by name via indexing
```

## Orbital elements and Kepler helpers

```python
system.get_semi_major_axis("earth")      # a [m]
system.get_eccentricity("earth")         # e
system.set_semi_major_axis("earth", 1.5 * AU)
system.set_eccentricity("earth", 0.1)

mu = system.calc_gravitational_parameter("earth")            # G (M_host + M_world) [m^3 s-2]
n  = system.calc_orbital_frequency("earth")                  # sqrt(mu / a^3) [rad s-1]
a  = system.calc_semi_major_axis_from_frequency("earth", n)  # inverse [m]
```

The mean motion follows Kepler's third law using the combined host + world mass.
`calc_gravitational_parameter` raises `RuntimeError` if no host is set and returns NaN for the host's
own entry (the host has no orbit within the system); `calc_orbital_frequency` returns NaN for an unset
or non-positive semi-major axis.

## Star and insolation

The star is designated with `is_star` (or `set_star`) and drives insolation. Each world carries its own
orbit about the star (independent of the tidal-host orbit):

```python
system.star                              # the star world (or None)
system.star_index
system.set_star("sun")                   # by index / name / object
system.get_star_luminosity()             # L_star [W]

system.set_stellar_semi_major_axis("earth", AU)
system.set_stellar_eccentricity("earth", 0.0167)
system.get_stellar_semi_major_axis("earth")
system.calc_stellar_orbital_frequency("earth")   # n about the star [rad s-1]

flux = system.calc_insolation_flux("earth")            # W/m^2 (orbit-averaged)
temp = system.calc_equilibrium_temperature("earth")    # K
```

`calc_insolation_flux` is the orbit-averaged incident stellar flux
`F = L_star / (4 π a² √(1-e²))` (the `√(1-e²)` is the time-average of `1/r²` over the eccentric orbit,
Mendez & Rivera-Valentin 2017), using the world's orbital elements **about the star**.
`calc_equilibrium_temperature` applies the world's own grey-body radiative balance
`T = ((1-A) F / (4 ε σ))^(1/4)` to that flux (using the world's albedo `A` and emissivity `ε`). Both
raise `RuntimeError` if no star is set and return NaN for the star's own entry, an unset stellar
semi-major axis, or a star with no luminosity.

## Orbital and spin evolution

`calc_world_evolution(world)` evolves a single orbiting world. It solves the world's global tides in the
current system state — mean motion from Kepler's third law, spin and obliquity from the world,
eccentricity and semi-major axis from the orbit about the host, host mass from the host world — then
turns the tidal-potential derivatives into the orbital rates and the world's spin rate. Only the
orbiting world raises tides; the host is treated as a point mass.

```python
ev = system.calc_world_evolution("moon")
ev["da_dt"], ev["de_dt"], ev["dn_dt"]     # orbital rates [m/s], [1/s], [rad/s^2]
ev["dspin_dt"]                            # spin rate [rad/s^2] (0 without a spin model)
ev["tidal_heating"]                       # [W]
ev["energy_residual"]                     # heating + dE_orbit/dt + dE_spin/dt (~0 under conservation)
```

The returned dict also carries the state used (`orbital_frequency`, `semi_major_axis`, `eccentricity`,
`spin_frequency`, `host_mass`, `target_mass`), the raw tidal outputs (`dU_dM`, `dU_dw`, `dU_dO`), the
`moment_of_inertia` and `has_spin` flag, and the energy terms (`dE_orbit_dt`, `dE_spin_dt`). `evolved`
is `False` for the host's own entry or a world with no usable orbit about the host; its rates are then
zero. `calc_system_evolution()` returns one such dict per world, in index order.

The rates follow the orbital rate engine (`dynamics_x`): `da/dt = (2/(n a)) dR/dM`,
`de/dt = (√(1-e²)/(n a² e))(√(1-e²) dR/dM − dR/dw)` (zero at `e = 0`), and `dn/dt = −(3/2)(n/a) da/dt`,
with the disturbing-function derivative `dR/dX = −((M_target + M_host)/M_target) dU/dX`. The spin rate
comes from the world's attached spin model (`dspin/dt = M_host dU/dO / I`). The heating and the orbit +
spin energy loss balance to `heating = −(dE_orbit/dt + dE_spin/dt)` with
`E_orbit = −G M_host M_world / (2 a)` and `E_spin = ½ I spin²`.

Each orbiting world evolves on its own two-body orbit about the host and dissipates independently.

### Dual-body dissipation

`calc_pair_evolution(world)` evolves an orbiting world **together with its host**, with both bodies
raising a tide on their shared orbit. Each body's tides are solved with the other body as the tide
raiser (masses swapped), so their orbital-rate contributions add and each body evolves its own spin:

```python
pair = system.calc_pair_evolution("moon")
pair["da_dt"], pair["de_dt"], pair["dn_dt"]     # combined shared-orbit rates
pair["tidal_heating_total"]                     # heating in both bodies
pair["energy_residual"]                         # ~0: heating_total + dE_orbit/dt + dE_spin_total/dt
pair["world"]                                   # the orbiting world's single-body contribution (a dict)
pair["host"]                                    # the host's single-body contribution (a dict)
```

The `world` and `host` entries are each a full `calc_world_evolution`-style dict (their own solve, spin,
heating, and share of the orbital rates + energy). The combined balance is the sum of the two
single-body balances, `heating_world + heating_host = −(dE_orbit/dt + dE_spin_world/dt + dE_spin_host/dt)`.
A body with no tide model attached is **rigid** and contributes nothing, so a rigid host reduces
`calc_pair_evolution` to the single-body `calc_world_evolution` result.

## C++ API

`c_System : c_TidalPyBaseClass` (`structures_x/system/system_.hpp`):

* `add_world(shared_ptr<c_BaseWorld>, is_host, is_star, a, e)` — owns worlds through `shared_ptr` so the
  C++ system and the Python wrappers co-own the same world.
* `get_num_worlds`, `get_world(i)`, `find_world_index(name)`.
* `set_host(i)`, `has_host`, `get_host_index`, `get_host`, `get_host_mass`.
* `set_star(i)`, `has_star`, `get_star_index`, `get_star`, `get_star_mass`, `get_star_luminosity`.
* `set/get_semi_major_axis(i)`, `set/get_eccentricity(i)` (orbit about the tidal host).
* `set/get_stellar_semi_major_axis(i)`, `set/get_stellar_eccentricity(i)` (orbit about the star).
* `calc_gravitational_parameter(i)`, `calc_orbital_frequency(i)`,
  `calc_semi_major_axis_from_frequency(i, n)` (host orbit), and the
  `calc_stellar_gravitational_parameter(i)` / `calc_stellar_orbital_frequency(i)` counterparts.
* `calc_insolation_flux(i)`, `calc_equilibrium_temperature(i)`.
* `calc_world_evolution(i)` / `calc_system_evolution()` — run the world's tidal solve for the current
  system state and return the orbital rates, spin rate, and energy terms as a `c_WorldEvolution` struct
  (a layered world is resolved via `dynamic_cast` so the rheology solve runs and the spin model is
  reached; a layerless world uses the analytic solve and contributes no spin). `calc_orbital_energy_derivative`,
  `calc_spin_energy_derivative`, `calc_energy_residual` compute the energy-balance terms.
* `calc_pair_evolution(i)` — dual-body evolution returning a `c_PairEvolution` (both bodies' `c_WorldEvolution`
  contributions + the combined shared-orbit rates and energy balance). Built on the shared
  `calc_dissipation(dissipator_i, companion_mass, n, a, e)` primitive, which computes one body's tidal
  solve + rate + spin contribution (a body with no tide model is rigid and contributes zero).

## References

* Boué & Efroimsky (2019), *Celestial Mechanics and Dynamical Astronomy* — orbital rate equations used
  by the orbital-evolution wiring built on this container.
