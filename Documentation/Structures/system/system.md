# System (`Structures.system`)

_Updated: 2026-10-08_

A `System` groups worlds (a star, planets, moons) and computes their orbital evolution and insolation. Each world names its own tidal host (the body raising its tides), or none, and one world may be the star (the source of insolation). In Earth-Moon-Sun the Moon's host is the Earth, the Earth's can be the Moon, and the star is the Sun; for an exoplanet the star is also the host. Each world has two two-body orbits, about its host and about the star, and interacts with nothing else. TidalPy has no N-body dynamics: three or more mutually interacting bodies need external scripting.

## Building a `System`

```python
from TidalPy.constants import au
from TidalPy.Structures import StarWorld, System, TerrestrialWorld

system = System("sol")
star = StarWorld("sun", 6.957e8, 1.988e30)
earth = TerrestrialWorld("earth", 6.371e6, 5.972e24)

system.add_world(
    star,
    is_star=True)          # Nothing raises the star's tides here
system.add_world(
    earth,
    tidal_host=star,       # The star is also the host
    semi_major_axis=au,
    eccentricity=0.0167)
```

`add_world(world, tidal_host=None, is_star=False, semi_major_axis=None, eccentricity=None, synchronous=False, stellar_semi_major_axis=None, stellar_eccentricity=None)` returns the world's index.

`tidal_host` is a world already in the system (index, name, or object), changed with `set_tidal_host`; a world without one is not tidally forced and evolution skips it. `semi_major_axis` and `eccentricity` give the orbit about the host, and the `stellar_` pair the orbit about the [star](#star-and-insolation). `None` leaves an element unset (an unset eccentricity reads as 0). `synchronous=True` sets the spin to the mean motion ([Synchronous Rotation](#synchronous-rotation)).

When the star is not a world's tidal host (_e.g._, Earth-Moon-Sun):

```python
from TidalPy.Structures import build_world

sun = build_world("sol")
earth = build_world("earth_simple")
moon = build_world("luna")

earth_moon_sun = System()
earth_moon_sun.add_world(
    sun,
    is_star=True)
earth_moon_sun.add_world(earth)
earth_moon_sun.add_world(
    moon,
    tidal_host=earth,
    semi_major_axis=3.84748e8,                   # About the Earth (3.844e8 m is the mean distance)
    eccentricity=0.0549,
    synchronous=True,
    stellar_semi_major_axis=au,                  # About the Sun
    stellar_eccentricity=0.0167)
earth_moon_sun.set_tidal_host(earth, moon)       # A mutual pair
```

### Mutual Pairs

Two worlds that host each other share one orbit: each element can come from either member (both must agree), and setting it on one sets it on both. `is_mutual_pair(world)` reports the relation. `calc_system_evolution` gives each member its own entry, and their contributions to the shared orbit add, as in `calc_pair_evolution`.

### Synchronous Rotation

`set_synchronous_rotation(world)` (or `add_world(..., synchronous=True)`, or `synchronous = true` in a system file) sets the world's spin to `calc_orbital_frequency(world)`, as for most large moons, and returns it \[rad s-1\]. It is set once, so call it again after changing the orbit or the masses. A world with no tidal host or semi-major axis raises `ValueError`.

```python
earth_moon_sun.set_synchronous_rotation(moon)    # [rad s-1]
```

## Building a `System` from TOML (`build_system`)

A system file names each member in a `[worlds.<name>]` table: its `world` (a bundled name, a world file, or an inline table), `tidal_host`, star role, orbital elements, and whether it rotates synchronously ([System Schema](../config/toml_schema.md#system-schema)). The table key becomes the world's name, so one world template can be reused under several names.

```python
from TidalPy.Structures import build_system

system = build_system("sol_system")      # a bundled system name, a .toml or binary path, or a dict
# equivalently: System.build("sol_system")
system.calc_insolation_flux("earth")     # ~1361 W/m^2

prem_system = build_system(
    "sol_system",
    overrides={"worlds": {"earth": {"world": "earth_prem"}}})   # The same system with the PREM Earth
```

`build_system(source, overrides=None, force=False)` mirrors `build_world`: it validates the source and builds each member with `build_world`; `overrides` needs only the keys it changes. A binary file is loaded with `load_system` ([Serialization](#serialization)). An unknown name or key raises an error naming the closest accepted one.

`system.save_to_toml(path)` writes the system as it is now, changed orbits included. An unchanged member built from a bundled name or file keeps that reference (relative paths rewritten for the destination); any other member is written inline. A world whose host is the star is written with its tidal elements only. `system.get_config_dict()` inlines every world's live `get_config_dict()`, and `build_system` rebuilds from it. `source_config` keeps the normalized build configuration. An inline world's relative `data_file` is found beside the system file.

## Worlds, Host, and Identification

Most methods take a world by index (any integer type, negatives allowed; not a `bool`), by name, or as the object.

```python
system.get_tidal_host("earth")        # the world that raises Earth's tides (or None)
system.get_tidal_host_index("earth")  # its index, or -1
system.has_tidal_host("sun")          # False: nothing forces the star here
system.set_tidal_host("earth", "sun") # name a host by index / name / object (None removes it)
system.num_worlds                 # 3
system.worlds                     # [sun, earth, jupiter]

for world in system: ...          # a sequence over its worlds
len(system)                       # 3
system[0]                         # the star
system[1:]                        # [earth, jupiter]
system.earth                      # by name, as an attribute
system["earth"]                   # or by indexing
```

## Orbital Elements and Kepler Helpers

```python
system.get_semi_major_axis("earth")      # a [m]
system.get_eccentricity("earth")         # e
system.set_semi_major_axis("earth", 1.5 * au)
system.set_eccentricity("earth", 0.1)

mu = system.calc_gravitational_parameter("earth")            # G (M_host + M_world) [m^3 s-2]
n  = system.calc_orbital_frequency("earth")                  # sqrt(mu / a^3) [rad s-1]
a  = system.calc_semi_major_axis_from_frequency("earth", n)  # inverse [m]
```

The mean motion uses Kepler's third law with the combined mass. `calc_gravitational_parameter` and `calc_orbital_frequency` return NaN for a world with no tidal host, the latter also for an unset or non-positive semi-major axis.

## Star and Insolation

The star is set with `is_star` or `set_star` and drives insolation. Each world has its own orbit about the star, separate from its tidal-host orbit unless its host is the star.

```python
system.star                              # the star world (or None)
system.star_index
system.set_star("sun")                   # by index / name / object
system.get_star_luminosity()             # L_star [W]

system.set_stellar_semi_major_axis("earth", au)
system.set_stellar_eccentricity("earth", 0.0167)
system.get_stellar_semi_major_axis("earth")
system.calc_stellar_orbital_frequency("earth")   # n about the star [rad s-1]

flux = system.calc_insolation_flux("earth")            # W/m^2 (orbit-averaged)
temp = system.calc_equilibrium_temperature("earth")    # K
```

`calc_insolation_flux` is the orbit-averaged incident stellar flux, from the world's orbit about the star,

$$F = \frac{L_{\star}}{4\pi a^{2}\sqrt{1-e^{2}}},$$

where $\sqrt{1-e^{2}}$ comes from the time average $\langle 1/r^{2} \rangle = 1/(a^{2}\sqrt{1-e^{2}})$ over the eccentric orbit (Méndez and Rivera-Valentín 2017). A star gives the same flux without a system: `star.calc_insolation_flux(distance, eccentricity=0.0)`. `calc_equilibrium_temperature` applies the world's gray-body radiative balance,

$$T = \left(\frac{(1-A)\,F}{4\,\varepsilon\,\sigma}\right)^{1/4},$$

with the world's albedo $A$, emissivity $\varepsilon$, and the Stefan-Boltzmann constant $\sigma$. Both raise `RuntimeError` if no star is set, and return NaN for the star itself, an unset stellar semi-major axis, or a star with no luminosity.

## Orbital and Spin Evolution

Run `world.solve_eos()` on each world on the `rheology` tide model (the builder's default for terrestrial worlds) before any evolution method, which otherwise raises `RuntimeError`. A later `solve_eos` discards the tidal result, and the next evolution call solves again. A world constructed in Python has no tide model until one is set.

`calc_world_evolution(world)` solves the world's tides in the current system state (Kepler mean motion, the world's spin and obliquity, its orbit about the host, the host's mass) and returns the orbital-element and spin rates. The host is a point mass. A world with no host or usable orbit returns `evolved = False`. Threads sharing a world take turns on it.

```python
earth.solve_eos()                         # Needed for the Love numbers
moon.solve_eos()
ev = earth_moon_sun.calc_world_evolution(moon)  # The Moon about the Earth
ev["da_dt"], ev["de_dt"], ev["dn_dt"]     # orbital rates [m/s], [1/s], [rad/s^2]
ev["dspin_dt"]                            # spin rate [rad/s^2]
ev["tidal_heating"]                       # [W]
ev["energy_residual"]                     # heating + dE_orbit/dt + dE_spin/dt (~0 under conservation)
```

The dict also holds `world_index`, `world_name`, the state used (`orbital_frequency`, `semi_major_axis`, `eccentricity`, `spin_frequency`, `host_mass`, `target_mass`), `dU_dM`, `dU_dw`, `dU_dO`, `moment_of_inertia` (from `solve_eos`, else the spin model's `moment_of_inertia_factor`; every world, stars included, has a spin model), `has_spin`, `has_tide_model`, `dE_orbit_dt`, and `dE_spin_dt`. `calc_system_evolution()` returns one dict per world, in index order. A rigid world (no tide model) returns `evolved = True` with zero rates and `has_tide_model = False`, with a one-time warning. `world.get_tide_state()` gives a member's current state as a dict in the argument order of `calc_tides` (`None` without a system, host, or usable orbit).

The rates follow the orbital rate engine ([Dynamics](../../Dynamics/dynamics.md)). The tidal-potential derivatives $\partial U/\partial X$ of the [global tides](../../Tides/global_tides.md) become disturbing-function derivatives $\partial\mathcal{R}/\partial X = -\frac{M_{w} + M_{h}}{M_{w}}\,\partial U/\partial X$, for world mass $M_{w}$ and host mass $M_{h}$, and

$$\frac{da}{dt} = \frac{2}{na}\,\frac{\partial\mathcal{R}}{\partial\mathcal{M}}, \qquad \frac{de}{dt} = \frac{\sqrt{1-e^{2}}}{na^{2}e}\left(\sqrt{1-e^{2}}\,\frac{\partial\mathcal{R}}{\partial\mathcal{M}} - \frac{\partial\mathcal{R}}{\partial\varpi}\right), \qquad \frac{dn}{dt} = -\frac{3}{2}\,\frac{n}{a}\,\frac{da}{dt},$$

with $de/dt = 0$ at $e = 0$. At small $e$ the bracket is of order $e^{2}$ while each term is of order one (non-synchronous spin), so it is evaluated as $-\frac{e^{2}}{1 + \sqrt{1-e^{2}}}\,\partial\mathcal{R}/\partial\mathcal{M} + \partial\mathcal{R}/\partial(\mathcal{M} - \varpi)$, the last term summed mode by mode, which keeps $de/dt / e$ exact at any small $e$. The spin rate is $\ddot{\theta} = (M_{h}/C)\,\partial U/\partial\Omega$, with $C$ the polar moment of inertia. Heating balances the orbit and spin energy loss,

$$\dot{E} = -\left(\frac{dE_\mathrm{orbit}}{dt} + \frac{dE_\mathrm{spin}}{dt}\right), \qquad E_\mathrm{orbit} = -\frac{G M_{h} M_{w}}{2a}, \qquad E_\mathrm{spin} = \frac{1}{2}\,C\,\dot{\theta}^{2}.$$

> [!NOTE]
> A system changes a world's `spin_frequency` only when asked (`synchronous`). A bundled world's spin is the synchronous rate of a rounded period, about 2e-4 off the Kepler mean motion for Io and Europa. In a strongly dissipative world that adds a slow tide at 2(n - spin) that can dominate the heating, so call `system.set_synchronous_rotation(world)` first. A spin within 1e-3 of the mean motion, but not equal, warns once per world.

### Dual-Body Dissipation

`calc_pair_evolution(world)` lets a world and its host raise tides on each other: each is solved with the other as raiser, their orbital rates add, and each evolves its own spin:

```python
pair = earth_moon_sun.calc_pair_evolution(moon)
pair["da_dt"], pair["de_dt"], pair["dn_dt"]     # combined shared-orbit rates
pair["tidal_heating_total"]                     # heating in both bodies
pair["energy_residual"]                         # ~0: heating_total + dE_orbit/dt + dE_spin_total/dt
pair["world"], pair["host"]                     # each body's own contribution (dicts)
```

The result also has `world_name` and `host_name`. The balance is $\dot{E}_{w} + \dot{E}_{h} = -\left(dE_\mathrm{orbit}/dt + dE_{\mathrm{spin},w}/dt + dE_{\mathrm{spin},h}/dt\right)$. A rigid body contributes nothing, so a rigid host (a star as a point mass) reduces the result to `calc_world_evolution`'s, without a warning. The top-level `has_tide_model` is `True` when either body has a tide model; when neither does, every rate is zero and a warning is logged once per world.

### Evolving a World About Its Host

`System.evolve(world, (t_start, t_end))` integrates the `calc_pair_evolution` rates with CyRK's LSODA: $a$, $e$, both spins, and, with `evolve_thermal=True` (the default), each layer's temperature. It returns an `EvolutionResult` of every step and leaves the system at the final state.

```python
import numpy as np
from TidalPy.constants import au, year

planet = build_world(
    "earth_thermal",
    {"radial_solver": {"rtol": 1.0e-10, "atol": 1.0e-14}})       # A tight radial solve (see below)
sun = build_world("sol")
sun.set_spin_frequency(2.0 * np.pi / (25.38 * 86400.0))          # Sidereal rotation [rad s-1]
system = System("Exoplanet")
system.add_world(
    sun,
    is_star=True)
system.add_world(
    planet,
    tidal_host=sun,
    semi_major_axis=0.15 * au,
    eccentricity=0.2)
system.set_tidal_host(sun, planet)                                 # The Sun raises tides too
planet.set_spin_frequency(10.0 * system.calc_orbital_frequency(planet))

result = system.evolve(planet, (0.0, 5.0e9 * year))               # 5 Gyr
result.time, result.semi_major_axis, result.eccentricity          # [s], [m], -
result.spin_ratio, result.tracked                                  # spin / n; True where on an equilibrium
result.temperature, result.tidal_heating                           # [K] (layers by steps), [W]
result.segments                                                    # Captures, recenterings, releases
```

Near a commensurability $s = k/2$ the resonant mode's torque changes sign, and the spin can have a stable equilibrium: a zero of its balance $b(s) = ds/dt = (\ddot{\theta} - s\,\dot{n})/n$ with a restoring torque on both sides. The orbit-averaged spin equation is first order, so a despinning world that reaches one stays. The spin relaxes in years to kyr while the orbit and interior change over Myr to Gyr, and in a cold mantle the equilibrium is about $10^{-6}$ wide in $s$, narrower than an integrator's usual tolerance.

`evolve` therefore runs in segments. Between equilibria the spin is *free*, a state variable. On one it is *tracked*:

* Only the slow state is integrated.
* Each rate evaluation finds that state's equilibrium $s^*$ by a bracketed root search (Illinois regula falsi to `root_tolerance`, warm-started from the last root).
* The rates are the combination of the two sides' rates that zeroes the balance (the Filippov combination), so the spin and orbit take one torque.

This neglects the spin's drift with its equilibrium, of relative size its relaxation time over the evolution time (as do Walterova and Behounkova 2020). Transitions are integration events on the state alone:

* Capture: a free spin within `capture_band` of the next $k/2$ probes its balance, in the direction it points, on a grid refined toward $k/2$; the first sign change brackets the stable equilibrium it will reach. Within `capture_margin` of it the spin is set on it; otherwise a free segment runs until it comes within `capture_margin` or enters another band.
* Recentering and release: a tracked spin carries a window about its equilibrium in which samples spaced by factors of 4 show no other sign change. When the balance at an edge loses its restoring sign, the equilibrium is solved again and the window rebuilt, or, when no equilibrium is left nearby (it met its unstable neighbor, as in a warming mantle), the spin is released.

The net torque near an equilibrium is a small difference of large mode torques, so the balance carries the radial solve's error times about the quality factor: for `earth_thermal`, about $10^{-3}$ per Myr at a radial `rtol` of $10^{-8}$, and $6 \times 10^{-6}$ at $10^{-10}$. Windows are tested no closer than `resolution` to the equilibrium. A radial `rtol` of $10^{-10}$ costs about 1.3 times one of $10^{-8}$ and keeps a warm equilibrium from being released on noise. The orbit, which changes by fractions of a percent over Gyr, uses `orbit_rtol` ($10^{-9}$); temperatures and a free spin use `thermal_rtol` and `spin_rtol` ($10^{-4}$).

With `evolve_thermal`, each new state re-solves the EOS with its temperature profile, a surface at the insolation temperature (or `surface_temperature`), and its heat sources at that time, and each layer warms at `calc_layer_temperature_rate`; without it the solved structure is kept. The host's spin is always free. A tracked step costs about 20 to 30 pair solves; most run time goes to recentering near a vanishing equilibrium and to steps at a lock whose resonant mode is nearly static (see `minimum_complex_rigidity` in the [configuration](../../Overview/2_TidalPy_Configurations.md)). `Demos/Systems/S02_thermal_orbital_evolution.ipynb` runs the example above. A solve failing at an integrator trial state only restarts the segment from its last step; a run that cannot continue returns `success` False and a `message` ([Limits and Failure Modes](#limits-and-failure-modes)).

## Serialization

A system's binary file ([Binary Serialization](../../Utilities/binary.md)) holds the name, star, hosts, orbits, and each world's full record, so each world returns as its own type with its models and settings. Solved state is not saved: re-run `solve_eos` on each layered world before evolving.

```python
import copy
import pickle
from TidalPy.Structures import load_system

system.save_binary("sol_system.tpyb")            # A str or an os.PathLike path
loaded = load_system("sol_system.tpyb")          # A new System, each world its own class
same = build_system("sol_system.tpyb")           # build_system loads a binary file too
twin = system.copy()                             # Also used by copy.copy, copy.deepcopy, and pickle
pickled = pickle.loads(pickle.dumps(system))     # e.g., for a process pool
print(system)                                    # System('Sol System', worlds=['sun', 'earth', 'jupiter'], star='sun')
```

`System.load_binary(path)` loads into an existing system. `copy()` uses the binary record, so the copy's worlds are unsolved, but it keeps the build configurations.

## Limits and Failure Modes

* `add_world` refuses, adding nothing: an element out of range, `synchronous=True` without a `tidal_host` and `semi_major_axis`, and stellar elements for a world whose host is the star.
* `ValueError`: a semi-major axis that is not positive or an eccentricity outside $[0, 1)$ (from `add_world`, a setter, or a file); a duplicate world name or object; a mutual pair whose members disagree on an element (from any method reading the orbit); a system file stating `semi_major_axis_m`, `eccentricity`, or `synchronous` with no `tidal_host`, or a world as its own host.
* `RuntimeError`: an evolution method on a `rheology` world whose EOS is not solved; insolation with no star set.
* `evolve` raises `ValueError` for a world with no tidal host or no prograde spin, a time span that does not increase, with `evolve_thermal` a world with no layers or a layer temperature that is not finite and positive, or capture settings out of order (`root_tolerance` < `resolution` < `capture_margin` <= `capture_band` < 0.25).
* `evolve` stops with `success` False on a tide or EOS solve failing at a reached state, on `max_wall_time`, or after five segments in a row end where they began or ten end on failed evaluations.
* `load_system` raises `IOError` naming what a non-system file holds. `load_binary` raises `IOError` for a corrupt record (a bad host or star index, two worlds with one name, an unbound orbit, trailing data) and leaves the system unchanged.

## C++ API

`c_System : c_TidalPyBaseClass, c_TideStateProvider`:

* `add_world(shared_ptr<c_BaseWorld>, is_star, a, e)`: co-owns the world and becomes its tide-state provider.
* `get_num_worlds`, `get_world(i)`, `find_world_index(name)`.
* `set_tidal_host(i, host_i)`, `clear_tidal_host(i)`, `has_tidal_host(i)`, `get_tidal_host_index(i)`, `get_tidal_host(i)`, `get_tidal_host_mass(i)`, `is_mutual_pair(i)`, and `get_host_orbit(i)` (the elements about the host, merged element by element with a mutual partner's).
* `get_tide_state(i, state_out)` and `get_equilibrium_temperature(i)`: the `c_TideStateProvider` interface, which a world reaches through `c_BaseWorld::get_tide_state(state_out)`.
* `set_star(i)`, `has_star`, `get_star_index`, `get_star`, `get_star_mass`, `get_star_luminosity`.
* `set/get_semi_major_axis(i)`, `set/get_eccentricity(i)` (orbit about the tidal host); `set/get_stellar_semi_major_axis(i)`, `set/get_stellar_eccentricity(i)` (orbit about the star).
* `calc_gravitational_parameter(i)`, `calc_orbital_frequency(i)`, `calc_semi_major_axis_from_frequency(i, n)`, and `calc_stellar_gravitational_parameter(i)`, `calc_stellar_orbital_frequency(i)`.
* `calc_insolation_flux(i)`, `calc_equilibrium_temperature(i)`.
* `calc_world_evolution(i)` and `calc_system_evolution()`: return the orbital rates, spin rate, and energy terms as a `c_WorldEvolution` (which also carries `dU_dM_minus_dw`; `has_tide_model` is false for a rigid world). `calc_orbital_energy_derivative`, `calc_spin_energy_derivative`, and `calc_energy_residual` give the energy-balance terms.
* `calc_pair_evolution(i)`: a `c_PairEvolution` (both bodies' `c_WorldEvolution` plus the combined rates and balance), built on `calc_dissipation(dissipator_i, companion_mass, n, a, e)`, one body's solve and contribution at its own spin (zero, unwarned, for a rigid body).
* Loading a corrupt binary record throws `std::runtime_error`. `System.evolve` is Python only.

## References

* Boué, G., and Efroimsky, M. (2019). Tidal evolution of the Keplerian elements. *Celestial Mechanics and Dynamical Astronomy*, 131(7), 30. The orbital rate equations used by the evolution methods on this container.
* Filippov, A. F. (1988). *Differential Equations with Discontinuous Righthand Sides*. Kluwer Academic Publishers. The combination of the two sides of a tracked equilibrium.
* Makarov, V. V., and Efroimsky, M. (2013). No pseudosynchronous rotation for terrestrial planets and moons. *The Astrophysical Journal*, 764(1), 27. Spin-orbit equilibria of viscoelastic planets.
* Méndez, A., and Rivera-Valentín, E. G. (2017). The equilibrium temperature of planets in elliptical orbits. *The Astrophysical Journal Letters*, 837(1), L1. The orbit-averaged insolation.
* Walterová, M., and Běhounková, M. (2020). Thermal and orbital evolution of low-mass exoplanets. *The Astrophysical Journal*, 900(1), 24. A quasi-static spin in coupled thermal-orbital evolution.
