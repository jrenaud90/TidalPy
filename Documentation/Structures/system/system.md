# System (`Structures.system`)

_Updated: 2026-10-09_

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

Three methods share one calculation of a world's tides on its orbit about its tidal host (in a mutual pair, the orbit both share):

* `calc_dissipation(world)` solves the world's tides in the current system state (Kepler mean motion, the world's spin and obliquity, its orbit about the host, the host's mass) and returns the tidal outputs without rates: `tidal_heating`, `dU_dM`, `dU_dw`, `dU_dO`, `dU_dM_minus_dw`, `moment_of_inertia`, the state used, and `companion_name`. `solved` is `False` for a world with no host or usable orbit. The per-layer heating stays on the world (`get_layer_tidal_heating`).
* `calc_world_evolution(world)` turns that into the orbital-element and spin rates with only this world dissipating: the host is a point mass whose state stays as it is, apart from the orbit they share. A world with no host or usable orbit returns `evolved = False`.
* `calc_pair_evolution(world, partner=None)` lets both dissipate ([Dual-Body Dissipation](#dual-body-dissipation)).

Threads sharing a world take turns on it.

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

`calc_pair_evolution(world, partner)` lets two worlds raise tides on each other: each is solved with the other as raiser (`calc_dissipation`), their orbital rates add, and each evolves its own spin. The two are treated alike. One must be the other's tidal host, usually each the other's; their orbit is the hosted world's (either one, in a mutual pair). `partner` defaults to `world`'s tidal host. Two worlds neither of which hosts the other, or one world twice, raise `ValueError`:

```python
pair = earth_moon_sun.calc_pair_evolution(moon, earth)
pair["da_dt"], pair["de_dt"], pair["dn_dt"]     # combined shared-orbit rates
pair["tidal_heating_total"]                     # heating in both bodies
pair["energy_residual"]                         # ~0: heating_total + dE_orbit/dt + dE_spin_total/dt
pair["world_names"]                             # ("moon", "earth"): world, then partner
pair["worlds"]["earth"]                         # each body's own contribution, keyed by world name
```

Each entry of `worlds` is a `calc_world_evolution`-style dict. A world with no tidal host and no `partner` returns `evolved = False`, `world_names` ending in `None`, and an empty `worlds`. The balance is $\dot{E}_{w} + \dot{E}_{h} = -\left(dE_\mathrm{orbit}/dt + dE_{\mathrm{spin},w}/dt + dE_{\mathrm{spin},h}/dt\right)$. A rigid body contributes nothing, so a rigid host (a star as a point mass) reduces the result to `calc_world_evolution`'s, without a warning. The top-level `has_tide_model` is `True` when either body has a tide model; when neither does, every rate is zero and a warning is logged once per world.

### Evolving a World About Its Host

`System.evolve(world, (t_start, t_end))` integrates the `calc_pair_evolution` rates of `world` and its tidal host in C++ with CyRK's implicit solvers. The state is $a$, $e$, the spin of each body with a tide model, and, with `evolve_thermal` (the default), the temperature of every layer of each body that has layers. A body without a tide model is rigid and keeps its spin. It returns a `PairedEvolutionResult`, a mapping of each world's name to its `EvolutionResult`, and leaves the system at the final state.

```python
import numpy as np
from TidalPy.constants import au, year

planet = build_world("earth_thermal", {"name": "Planet"})
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

pair = system.evolve(planet, (0.0, 5.0e9 * year))                 # 5 Gyr
pair.time, pair.semi_major_axis, pair.eccentricity                # [s], [m], -
pair["Planet"].spin_ratio, pair["Planet"].tidal_heating           # spin / n, [W]
pair["Planet"].temperature                                        # [K], layers by steps
pair.num_rhs_calls, pair.num_jacobians, pair.segments             # The integration's cost and restarts
```

A degree-$l$ mode of order $m$ is resonant where $m s$ is an integer, with $s$ the spin frequency over the mean motion, so a body's spin-orbit commensurabilities are the ratios $j/m$ with $m$ up to its largest tidal degree. Near one the resonant mode's torque changes sign, and the spin can have a stable equilibrium (a lock) where the torque restores it from both sides. A mode slower than the world's continuation frequency has its dissipation fall linearly to zero (see [Numerical Settings](../../Overview/2_TidalPy_Configurations.md#numerical-settings)), so the torque is a smooth, odd function through each commensurability and a lock is an ordinary stable equilibrium. The spin relaxes onto it within years while the orbit and interior change over Myr to Gyr, which makes the spin a stiff variable that an implicit integrator holds with long steps. Capture, passage, and release follow from the equations themselves.

Each spin is carried as its offset $\delta = s - j/m$ from the nearest commensurability. The offset enters the mode frequencies exactly, so a spin within rounding of a lock (Charon's synchronous offset is near $10^{-16}$) keeps a smooth torque. The orbit is carried as its change, $a/a_0 - 1$ and $e - e_\mathrm{ref}$, with $e_\mathrm{ref}$ the initial $e$. The integration restarts, beginning a new entry of `segments`, when

* a spin crosses its commensurability, so a lock narrower than a step is not stepped over,
* a spin moves past the midpoint to the neighboring commensurability, which becomes its reference,
* $e$ falls to a tenth of $e_\mathrm{ref}$, which then becomes the new $e_\mathrm{ref}$, so a circularizing orbit stays accurate relative to $e$.

Settings left as `None` come from the system file's `[evolution]` table, then from `TidalPy.config["evolution"]` (see [Evolution](../../Overview/2_TidalPy_Configurations.md#evolution)). `method` is `"LSODA"` by default, which was the fastest and never failed over a benchmark of captures, held locks, thermal pairs, and a 5 Gyr thermal Earth. `"Radau"` takes the fewest steps on a long cold lock (cold Pluto-Charon over 1 Gyr in 1.7 s against 6 s), and `"BDF"` is also accepted. `semi_major_axis_rtol` and `eccentricity_rtol` apply to the changes in $a$ and $e$, and `spin_rtol` to a spin's offset from its commensurability. During a run each world's Love solves use `radial_rtol` and `radial_atol`, since the rates must be smooth at the scale of the integrator's difference steps, and the world's own settings return afterward.

With `evolve_thermal`, a body's EOS is solved again whenever its temperatures, the time, or its surface temperature change. Its surface sits at the insolation temperature of the system's star, and its layers warm at `calc_layer_temperature_rate`. Without a star no heat leaves the surface, and the run warns. Without `evolve_thermal` each world keeps the structure it has. Radiogenic heating counts time from the body's formation, so a present-day run starts near $4.5 \times 10^{9}$ years. `Demos/Systems/S02_thermal_orbital_evolution.ipynb` runs a thermal example. A run that cannot continue returns `success` False and a `message` ([Limits and Failure Modes](#limits-and-failure-modes)).

> [!NOTE]
> The bundled worlds' cooling models are simple. With them, a thermal Pluto has a conducting lid about 1 km thick and loses heat about 700 times faster than its radiogenic heating replaces it. Over 1 Gyr its hydrosphere cools from about 268 K to about 151 K, freezing the ocean, while its core, which has no melting law, warms from 600 K to about 1700 K. Charon behaves similarly (see [WorldPack](../config/worldpack.md)).

`System.evolve` assumes that

* the worlds have no permanent (triaxial) figure and no rotational flattening, so the spins follow the tidal torques alone,
* layer boundaries are fixed (a two-phase layer melts and freezes inside its own radii),
* a thermal body's structure is solved without its tidal heat, which enters only its temperature rates,
* each world's moment of inertia is constant in the spin equation, $\ddot{\theta} = (M_{h}/C)\,\partial U/\partial\Omega$ (no $\dot{C}$ term).

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
* `evolve` raises `ValueError` for a world with no tidal host or no prograde spin, a time span that does not increase, with `evolve_thermal` a layer temperature that is not finite and positive, a tolerance that is not positive, or a method that is not implicit.
* `evolve` stops with `success` False on a tide or EOS solve failing at a reached state, on `max_wall_time`, on a keyboard interrupt while `progress=True` shows its bar (`tqdm`, over the simulated time), or after five segments in a row end on failed evaluations or ten end where they began. It returns what was integrated.
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
* `calc_dissipation(i)`: one world's tide on its orbit about its tidal host, as a `c_TidalDissipation` (the state used and the tidal outputs, no rates; `solved` is false with no host or usable orbit). Built on the protected `p_dissipation(dissipator_i, companion_i, orbit_i)`, which `p_evolution` turns into a `c_WorldEvolution`.
* `calc_world_evolution(i)` and `calc_system_evolution()`: return the orbital rates, spin rate, and energy terms as a `c_WorldEvolution` (which also carries `dU_dM_minus_dw`; `has_tide_model` is false for a rigid world). `calc_orbital_energy_derivative`, `calc_spin_energy_derivative`, and `calc_energy_residual` give the energy-balance terms.
* `calc_pair_evolution(i, j)`: a `c_PairEvolution` (`first` and `second`, each body's `c_WorldEvolution`, plus the combined rates and balance), each body solved by `p_dissipation` with the other raising its tide (zero, unwarned, for a rigid body); `std::invalid_argument` when neither hosts the other or `i == j`. `calc_pair_evolution(i)` pairs a world with its tidal host and is unevolved without one.
* `c_evolve_pair(system_ptr, i, t_start, t_end, c_PairEvolveSettings)` (`evolution_.hpp`): `System.evolve`'s driver, returning a shared `c_PairEvolutionRecord`. `p_dissipation` also takes a spin as an exact commensurability plus offset (`c_mode_frequency`).
* Loading a corrupt binary record throws `std::runtime_error`.

## References

* Boué, G., and Efroimsky, M. (2019). Tidal evolution of the Keplerian elements. *Celestial Mechanics and Dynamical Astronomy*, 131(7), 30. The orbital rate equations used by the evolution methods on this container.
* Hairer, E., and Wanner, G. (1996). *Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems* (2nd ed.). Springer. The implicit Radau method (`method="Radau"`).
* Makarov, V. V., and Efroimsky, M. (2013). No pseudosynchronous rotation for terrestrial planets and moons. *The Astrophysical Journal*, 764(1), 27. Spin-orbit equilibria of viscoelastic planets.
* Méndez, A., and Rivera-Valentín, E. G. (2017). The equilibrium temperature of planets in elliptical orbits. *The Astrophysical Journal Letters*, 837(1), L1. The orbit-averaged insolation.
* Walterová, M., and Běhounková, M. (2020). Thermal and orbital evolution of low-mass exoplanets. *The Astrophysical Journal*, 900(1), 24. Coupled thermal-orbital evolution with a viscoelastic spin.
