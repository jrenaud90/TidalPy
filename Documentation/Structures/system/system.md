# System (`Structures.system`)

_Updated: 2026-10-06_

A `System` links two or more worlds (a star, planets, moons) into a gravitationally bound group. It tracks two roles independently:

- A tidal host per world: the body that raises a tide on that world. Every world names its own host, or none.
- The star: the body that provides insolation to each world in the system.

In the Earth-Moon-Sun system the Moon's tidal host is the Earth, the Earth's can be the Moon, and the star is the Sun. For an exoplanet orbiting its star, the star is also its tidal host, and the orbits about the host and about the star coincide. Each world therefore carries two two-body orbits: one about its tidal host and one about the star. A world interacts only with its tidal host and the star. The `System` is the container on which orbital evolution and insolation are computed.

TidalPy has no N-body dynamics and instead utilizes semi-analytic two-body evolution (orbit/spin) rates. A problem with three or more mutually interacting bodies can be built from TidalPy worlds but needs external scripting to account for the additional forcing.

## Building a `System`

```python
from TidalPy.constants import au
from TidalPy.Structures import StarWorld, System, TerrestrialWorld

system = System("sol")
star = StarWorld("sun", 6.957e8, 1.988e30)
earth = TerrestrialWorld("earth", 6.371e6, 5.972e24)

system.add_world(
    star,
    is_star=True)          # The star is tidally forced by nothing here, so it names no host
system.add_world(
    earth,
    tidal_host=star,       # Exoplanet case: the star is also the world's tidal host
    semi_major_axis=au,
    eccentricity=0.0167)
```

`add_world(world, tidal_host=None, is_star=False, semi_major_axis=None, eccentricity=None, synchronous=False, stellar_semi_major_axis=None, stellar_eccentricity=None)` returns the world's index. `tidal_host` is a world already in the system, given by index, name, or object; `set_tidal_host(world, tidal_host)` names or changes it afterwards. There is no system-wide host. A world with no tidal host is not tidally forced and is skipped by the evolution methods. `is_star` marks the insolation source. `semi_major_axis` and `eccentricity` describe the orbit about the tidal host, and `stellar_semi_major_axis` and `stellar_eccentricity` the orbit about the star (see [Star and Insolation](#star-and-insolation)). `None` leaves an element unset, and an eccentricity that no element set reads as 0. `synchronous=True` sets the world's spin to its mean motion about its tidal host (see [Synchronous Rotation](#synchronous-rotation)). The system co-owns each world with its Python wrapper, so the wrapper you passed stays usable and is the same object the system hands back.

A world that `add_world` refuses is not added: an orbital element out of range, `synchronous=True` without a `tidal_host` and a `semi_major_axis`, or stellar elements for a world whose tidal host is the star (its orbit about the star is its tidal orbit, given by `semi_major_axis` and `eccentricity`).

For a system where the star is a separate body from a world's tidal host (_e.g._, Earth-Moon-Sun):

```python
from TidalPy.Structures import build_world

sun = build_world("sol")                         # Bundled Sun, named "Sol"
earth = build_world("earth_simple")              # Bundled layered Earth, named "Earth-Simple"
moon = build_world("luna")                       # Bundled layered Moon, named "Luna"

earth_moon_sun = System()
earth_moon_sun.add_world(
    sun,
    is_star=True)                                # The star, and nobody's tidal host
earth_moon_sun.add_world(earth)                  # The Moon's tidal host, named below
earth_moon_sun.add_world(
    moon,
    tidal_host=earth,                            # The Earth raises the Moon's tides
    semi_major_axis=3.84748e8,                   # Moon about the Earth (tidal); 3.844e8 m is the mean distance
    eccentricity=0.0549,
    synchronous=True,                            # The Moon's spin is its mean motion about the Earth
    stellar_semi_major_axis=au,                  # Moon about the Sun (insolation)
    stellar_eccentricity=0.0167)
earth_moon_sun.set_tidal_host(earth, moon)       # The Moon raises the Earth's tides in turn
```

### Mutual Pairs

Two worlds that host each other share one orbit, so each of its elements only needs to be specified once, so the semi-major axis can come from one member and the eccentricity from the other. An element both members give must agree, and a disagreement raises `ValueError` from any method that reads the orbit. Setting an element on either member (`set_semi_major_axis`, `set_eccentricity`) sets it on both, so a pair built with the elements on both members, as a saved system file has them, can still be updated from one side. `is_mutual_pair(world)` reports the relation. `calc_system_evolution` gives each member its own entry, and their contributions to the shared orbit add, which is what `calc_pair_evolution` returns for either member.

### Synchronous Rotation

A synchronously rotating world (most large moons) spins at its orbital mean motion. `set_synchronous_rotation(world)` sets the world's spin frequency to `calc_orbital_frequency(world)` and returns it [rad s-1]; `add_world(..., synchronous=True)` does the same as the world is added. The spin is set once, from the orbit at the time of the call, so call it again after changing the orbit or the masses. It raises `ValueError` for a world with no tidal host or no semi-major axis about it.

```python
earth_moon_sun.set_synchronous_rotation(moon)    # The Moon's spin [rad s-1], now its mean motion
```

## Building a `System` from TOML (`build_system`)

A whole system can be described in a single TOML, including all worlds (and their layers), and built in one call, mirroring `build_world`. Each `[worlds.<name>]` table names a `world` (a bundled world name, a path to a world TOML, or an inline world config) plus its `tidal_host` (the table key of the world that raises its tides), its star role, and its orbital elements. `system.get_config_dict()` returns the same schema for the system as it stands now (each world inlined with its own live `get_config_dict()`), and `build_system(config)` rebuilds the system from it. A host may be declared after the worlds it hosts. A world that states `semi_major_axis_m` or `eccentricity` must name a `tidal_host` for them to be about. The table key becomes the world's name within the system, so a bundled world template can be reused under different names.

```toml
schema_version = "0.2.0"
name = "Sol System"

[worlds.sun]
world = "sol"  # a bundled world name (also accepts a path or an inline [worlds.sun.world] table)
is_star = true

[worlds.earth]
world = "earth_simple"
tidal_host = "sun"                  # the world that raises this one's tides
semi_major_axis_m = 1.495978707e11  # orbit about the tidal host
eccentricity      = 0.0167  # the star is the host, so this is also the orbit about the star (insolation)
```

```python
from TidalPy.Structures import build_system

system = build_system("sol_system")      # a bundled system name, a .toml or binary path, or a dict
# equivalently: System.build("sol_system")
system.calc_insolation_flux("earth")     # ~1361 W/m^2 (the solar constant)

prem_system = build_system(
    "sol_system",
    overrides={"worlds": {"earth": {"world": "earth_prem"}}})   # The same system with the PREM Earth
```

`build_system(source, overrides=None, force=False)`, a thin wrapper over `System.build` that mirrors `build_world` and `BaseWorld.build`, resolves the source, validates it (schema version and structure), and builds each member world with `build_world`. `overrides` is a nested dict merged over the configuration before the build, table by table, so it only needs the keys it changes. A dict built in Python without a `schema_version` targets the current schema; a file without one is warned about. A path to a binary file is loaded with `load_system` (see [Serialization](#serialization)). An unknown bundled name, or an unknown key in a system or member table, raises an error naming the closest accepted one. To make the star and the tidal host different bodies, give a world a `tidal_host` other than the world marked `is_star` (_e.g._, a moon whose tidal host is its planet but whose insolation comes from the system star), and give each world both a tidal-host orbit (`semi_major_axis_m` and `eccentricity`) and a stellar orbit (`stellar_semi_major_axis_m` and `stellar_eccentricity`). A member's `world` may be a bundled name or a path to a world file. A relative path resolves against the system file's folder first, then the working directory, as a world's `data_file` resolves against the world file's folder first.

A system refuses what would give it no bound orbit or an ambiguous member, raising `ValueError`: a semi-major axis that is not positive, an eccentricity outside $[0, 1)$ (from `add_world`, the orbit setters, or a file), a second world with a name already in the system, and the same world object added twice. A world is named by its index (any integer type, numpy's included), its name, or the object itself; a `bool` is refused rather than read as index 0 or 1.

A built system retains its normalized configuration on `source_config` and can be written back out:

```python
system.save_to_toml("my_system.toml")
system_cfg = system.get_config_dict()
```

`save_to_toml` writes the system as it is now (`get_save_config(destination_dir)`): every member's current tidal host, star role, and orbital elements, so an orbit changed after the build is saved. A member built from a world reference (a bundled name or a file) and unchanged since its build keeps that reference, and a relative file path is rewritten to find the same file from the folder saved into. Any other member, a system assembled directly in Python included, is written inline as its world's `get_save_config()`. `get_config_dict` is the self-contained expansion that inlines every world's live `get_config_dict()` with its roles and orbital elements; it rebuilds through `build_system` too. A world whose tidal host is the star is written with its tidal elements only (`semi_major_axis_m` and `eccentricity`), since its orbit about the star is the same orbit. A world given inline in a system file finds a relative `data_file` beside the system file.

## Worlds, Host, and Identification

Most methods accept a world by index (int, negatives allowed), by name (str), or as the world object itself:

```python
system.get_tidal_host("earth")        # the world that raises Earth's tides (or None)
system.get_tidal_host_index("earth")  # its index, or -1
system.has_tidal_host("sun")          # False: nothing forces the star here
system.set_tidal_host("earth", "sun") # name a host by index / name / object (None removes it)
system.num_worlds                 # 3
system.worlds                     # [sun, earth, jupiter]
```

The system is a sequence over its worlds:

```python
for world in system: ...          # iterate members
len(system)                       # 3
system[0]                         # the star
system[1:]                        # [earth, jupiter]
system.earth                      # attribute access by world name
system["earth"]                   # or by name via indexing
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

The mean motion follows Kepler's third law using the combined host and world mass. `calc_gravitational_parameter` returns NaN for a world with no tidal host (it has no orbit within the system). `calc_orbital_frequency` also returns NaN for an unset or non-positive semi-major axis.

## Star and Insolation

The star is designated with `is_star` (or `set_star`) and drives insolation. Each world carries its own orbit about the star, independent of the tidal-host orbit, except for a world whose tidal host is the star.

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

`calc_insolation_flux` is the orbit-averaged incident stellar flux, using the world's orbital elements about the star,

$$F = \frac{L_{\star}}{4\pi a^{2}\sqrt{1-e^{2}}},$$

where the factor $\sqrt{1-e^{2}}$ comes from the time average of $1/r^{2}$ over the eccentric orbit, $\langle 1/r^{2} \rangle = 1/(a^{2}\sqrt{1-e^{2}})$ (Méndez and Rivera-Valentín 2017). A star gives the same flux at any orbit without a system: `star.calc_insolation_flux(distance, eccentricity=0.0)`. `calc_equilibrium_temperature` applies the world's own gray-body radiative balance to that flux,

$$T = \left(\frac{(1-A)\,F}{4\,\varepsilon\,\sigma}\right)^{1/4},$$

with the world's albedo $A$, its emissivity $\varepsilon$, and the Stefan-Boltzmann constant $\sigma$. Both raise `RuntimeError` if no star is set and return NaN for the star's own entry, an unset stellar semi-major axis, or a star with no luminosity.

## Orbital and Spin Evolution

A world whose tide model is `rheology` (the builder's default for a terrestrial world; a world constructed directly in Python has no tide model until one is set) takes its Love numbers from its interior, so run `world.solve_eos()` on each of these before any of the evolution methods below, they raise `RuntimeError` otherwise. A later `solve_eos` retires the world's tidal result, and the next evolution call solves it again.

`calc_world_evolution(world)` evolves a single world. It solves the world's global tides in the current system state (mean motion from Kepler's third law, spin and obliquity from the world, eccentricity and semi-major axis from the orbit about its tidal host, host mass from that host), then turns the tidal-potential derivatives into the orbital element time derivatives and the world's spin rate derivative. Its host is treated as a point mass with no tidal derivatives calculated (its dissipation does not affect the orbit). A world with no tidal host, or no usable orbit about it, comes back with `evolved = False`. The evolution methods hold each world's call lock from its tidal solve through the read of its result, so evolution calls on threads that share a world take turns on it and each reads its own solve.

A world that belongs to a system can be asked for the same state directly: `world.get_tide_state()` returns it as a dict in the argument order of `calc_tides`, or `None` for a world outside a system, with no tidal host, or with no usable orbit. Orbital state is never stored on a world; the system supplies it on request and stops doing so when it is deleted.

```python
earth.solve_eos()                         # Interior solve for the Earth's Love numbers
moon.solve_eos()                          # Interior solve for the Moon's Love numbers
ev = earth_moon_sun.calc_world_evolution(moon)  # The Moon about the Earth, from the Earth-Moon-Sun system above
ev["da_dt"], ev["de_dt"], ev["dn_dt"]     # orbital rates [m/s], [1/s], [rad/s^2]
ev["dspin_dt"]                            # spin rate [rad/s^2]
ev["tidal_heating"]                       # [W]
ev["energy_residual"]                     # heating + dE_orbit/dt + dE_spin/dt (~0 under conservation)
```

The returned dict also carries the world (`world_index`, `world_name`), the state used (`orbital_frequency`, `semi_major_axis`, `eccentricity`, `spin_frequency`, `host_mass`, `target_mass`), the raw tidal outputs (`dU_dM`, `dU_dw`, `dU_dO`), the `moment_of_inertia`, the `has_spin` and `has_tide_model` flags, and the energy terms (`dE_orbit_dt`, `dE_spin_dt`). Every world, stars included, carries a spin model. Either the `solve_eos`-derived MOI or the spin model's `moment_of_inertia_factor` (which comes from `[worlds]` in `TidalPy_Configs.toml`), is used for the spin evolution. `calc_system_evolution()` returns one such dict per world, in index order.

A rigid world (no tide model attached) raises no tide. Its entry comes back with `evolved = True`, zero rates and energy terms, and `has_tide_model = False`, and the system logs a warning the first time each such world is evolved.

The rates follow the orbital rate engine (`Dynamics`). With the tidal-potential derivatives $\partial U/\partial X$ of the [global tides](../../Tides/global_tides.md) converted to disturbing-function derivatives $\partial\mathcal{R}/\partial X = -\frac{M_{w} + M_{h}}{M_{w}}\,\partial U/\partial X$, for the world mass $M_{w}$ and the host mass $M_{h}$,

$$\frac{da}{dt} = \frac{2}{na}\,\frac{\partial\mathcal{R}}{\partial\mathcal{M}}, \qquad \frac{de}{dt} = \frac{\sqrt{1-e^{2}}}{na^{2}e}\left(\sqrt{1-e^{2}}\,\frac{\partial\mathcal{R}}{\partial\mathcal{M}} - \frac{\partial\mathcal{R}}{\partial\varpi}\right), \qquad \frac{dn}{dt} = -\frac{3}{2}\,\frac{n}{a}\,\frac{da}{dt},$$

with $de/dt = 0$ at $e = 0$. At small eccentricity the bracket is of order $e^{2}$ while each of its terms is of order one (for a spin that is not synchronous), so it is evaluated as $-\frac{e^{2}}{1 + \sqrt{1-e^{2}}}\,\partial\mathcal{R}/\partial\mathcal{M} + \partial\mathcal{R}/\partial(\mathcal{M} - \varpi)$, the last term summed mode by mode in the tidal collapse, where it holds no cancellation. $de/dt / e$ then stays exact however small $e$ gets. The spin rate comes from the world's attached spin model, $\ddot{\theta} = (M_{h}/C)\,\partial U/\partial\Omega$, with $C$ the polar moment of inertia. The heating and the orbit and spin energy loss balance,

$$\dot{E} = -\left(\frac{dE_\mathrm{orbit}}{dt} + \frac{dE_\mathrm{spin}}{dt}\right), \qquad E_\mathrm{orbit} = -\frac{G M_{h} M_{w}}{2a}, \qquad E_\mathrm{spin} = \frac{1}{2}\,C\,\dot{\theta}^{2}.$$

Each world evolves on its own two-body orbit about its tidal host and dissipates independently.

> [!NOTE]
> A world's spin is its own `spin_frequency`, which a system does not change. A bundled world's spin is the synchronous rate of a rounded rotation period, which differs slightly from the mean motion the system computes by Kepler's third law (by about 2e-4 for Io and Europa). For a strongly dissipative world that small difference adds a slow tide at 2(n - spin) that can dominate the heating, so set the spin from the system first to keep the world synchronous: `system.set_synchronous_rotation(world)`. A tide solve whose spin is within 1e-3 of the mean motion, but not equal to it, logs a warning once per world.

### Dual-Body Dissipation

`calc_pair_evolution(world)` evolves a world together with its own tidal host, with both bodies raising a tide on their shared orbit. Each body's tides are solved with the other body as the tide raiser (masses swapped), so their orbital-rate contributions add and each body evolves its own spin:

```python
pair = earth_moon_sun.calc_pair_evolution(moon)  # The Moon and the Earth raising tides on each other
pair["da_dt"], pair["de_dt"], pair["dn_dt"]     # combined shared-orbit rates
pair["tidal_heating_total"]                     # heating in both bodies
pair["energy_residual"]                         # ~0: heating_total + dE_orbit/dt + dE_spin_total/dt
pair["world"]                                   # the orbiting world's single-body contribution (a dict)
pair["host"]                                    # the host's single-body contribution (a dict)
```

The result names both bodies (`world_name`, and `host_name`, `None` for a world with no tidal host). The `world` and `host` entries are each a full `calc_world_evolution`-style dict (their own solve, spin, heating, and share of the orbital rates and energy). The combined balance is the sum of the two single-body balances, $\dot{E}_{w} + \dot{E}_{h} = -\left(dE_\mathrm{orbit}/dt + dE_{\mathrm{spin},w}/dt + dE_{\mathrm{spin},h}/dt\right)$. A body with no tide model attached is rigid and contributes nothing, so a rigid host reduces `calc_pair_evolution` to the single-body `calc_world_evolution` result. Each body's `has_tide_model` flag is in its `world` or `host` entry. The top-level `has_tide_model` is `True` when at least one body carries a tide model; when neither does, every rate is zero, the flag is `False`, and a warning is logged once per world. A rigid host beside a dissipating world (a star treated as a point mass) is a normal setup and is not warned about.

### Evolving a World About Its Host

`System.evolve(world, (t_start, t_end))` integrates the rates of `calc_pair_evolution` through time. It evolves the shared orbit ($a$, $e$), the host's spin, the world's spin, and, with `evolve_thermal=True` (the default), the temperature of each of the world's layers, with CyRK's LSODA, and returns an `EvolutionResult` holding every integration step. The system starts from its current state and is left at the final one.

```python
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
result.spin_ratio, result.tracked                                  # spin / n, and True where it sat on an equilibrium
result.temperature, result.tidal_heating                           # [K] (layers by steps), [W]
result.segments                                                    # The captures, recenterings, and releases
```

Near a spin-ratio commensurability $s = k/2$ the resonant tidal mode's torque changes sign, and the spin can have a stable equilibrium: a zero of its balance $b(s) = ds/dt = (\ddot{\theta} - s\,\dot{n})/n$ with a restoring torque on both sides. The orbit-averaged spin equation is first order, so a despinning world that reaches one stays. The spin relaxes onto it within years to kyr while the orbit and the interior change over Myr to Gyr, and in a cold mantle the equilibrium is about $10^{-6}$ wide in $s$, narrower than an implicit integrator's usual tolerance on $s$.

`evolve` therefore runs in segments. Between equilibria the spin is *free*, a state variable. On one it is *tracked*:

* Only the slow state is integrated.
* Each rate evaluation finds the equilibrium $s^*$ of that state by a bracketed root search warm-started from the last one (Illinois regula falsi to `root_tolerance`).
* The rates are the combination of the rates on the root's two sides that zeroes the balance (the Filippov combination), so the spin and the orbit take one torque.

The neglected term is the spin's drift with its equilibrium, of relative size its relaxation time over the evolution time (Walterova and Behounkova 2020 make the same approximation). The transitions are integration events, each a function of the state alone:

* A free spin entering the band $|s - k/2| <$ `capture_band` of the next commensurability runs a capture test: from the spin, in the direction its balance points, the balance is probed on a grid refined toward $k/2$, and the first sign change brackets the stable equilibrium the spin reaches. Within `capture_margin` of it the spin is set on it; farther, a free segment runs until it comes within `capture_margin` or enters another band.
* A tracked spin carries a window: an interval about its equilibrium in which sample points spaced by factors of 4 show no other sign change. When the balance at an edge loses its restoring sign, the equilibrium is solved again and the window rebuilt about it (a recentering), or, when no equilibrium is left nearby (it vanished with its unstable neighbor, as a warming mantle moves them together), the spin is released.

The net torque near an equilibrium is a small difference of large mode torques, so the balance carries the radial solve's error amplified by about the quality factor (for `earth_thermal` about $10^{-3}$ per Myr at a radial-solver `rtol` of $10^{-8}$, and $6 \times 10^{-6}$ at $10^{-10}$). Windows measure that noise and are tested no closer than `resolution` to the equilibrium. A radial `rtol` of $10^{-10}$, as above, costs about 1.3 times a solve at $10^{-8}$ and keeps a warm equilibrium from being released on noise. The orbit changes by fractions of a percent over Gyr, so it takes its own tight tolerance (`orbit_rtol`, $10^{-9}$); the temperatures and a free spin take `thermal_rtol` and `spin_rtol` ($10^{-4}$).

With `evolve_thermal`, every new state solves the world's EOS with its temperature profile, its surface at the system's insolation temperature (or `surface_temperature`), and its heat sources at the integration time, and each layer warms at `calc_layer_temperature_rate`; without it the world keeps its solved structure. The host's spin is always free. A tracked step costs a root search, about 20 to 30 pair solves, and most of a run's time goes to recentering near a vanishing equilibrium and to steps at a lock whose resonant mode is nearly static (see `minimum_complex_rigidity` in the [configuration](../../Overview/2_TidalPy_Configurations.md)). `Demos/Systems/S02_thermal_orbital_evolution.ipynb` evolves the example above for 5 Gyr.

`evolve` raises `ValueError` for a world with no tidal host or no prograde spin, a time span that does not increase, with `evolve_thermal` a world with no layers or a layer temperature that is not finite and positive, or capture settings out of order (`root_tolerance` < `resolution` < `capture_margin` <= `capture_band` < 0.25). A run that cannot continue returns with `success` False and its `message`:

* a tide or EOS solve failing at a reached state;
* `max_wall_time` reached;
* five segments in a row ending where they began, or ten ending on failed evaluations.

A solve failing at a trial state of the integrator only ends the segment, which restarts from its last step.

## Serialization

A system's binary file ([Binary Serialization](../../Utilities/binary.md)) carries the container state (name, the star role, and each world's tidal host and orbital elements about both that host and the star) and each world's complete record, so every world comes back as its own type with its models and settings. A loaded system evolves as the saved one did once each layered world has re-run `solve_eos`: solved state, the EOS profiles included, is not saved.

```python
import copy
import pickle
from TidalPy.Structures import load_system

system.save_binary("sol_system.tpyb")            # A str or an os.PathLike path
loaded = load_system("sol_system.tpyb")          # A new System, each world its own class
same = build_system("sol_system.tpyb")           # build_system loads a binary file too
twin = system.copy()                             # An independent copy (copy.copy, copy.deepcopy, and pickle use it)
pickled = pickle.loads(pickle.dumps(system))     # So a system can be sent to a process pool
print(system)                                    # System('Sol System', worlds=['sun', 'earth', 'jupiter'], star='sun')
```

`load_system(path, force=False)` reads a system file without a placeholder object and raises `IOError` naming what a file holds when it is not a system's (a world file, or a file that is not a TidalPy binary file). `System.load_binary(path)` loads a file into an existing system instead. `copy()` goes through the same binary record held in memory, so the copy's worlds are unsolved, and it keeps each world's and the system's own copies of the configurations they were built from.

`load_binary` raises `IOError` for a file whose system record is corrupt: a tidal host index that names no other world, an out-of-range star index, two worlds with one name, or an orbit that is not bound (a semi-major axis that is not positive or an eccentricity outside $[0, 1)$). Any failed load, trailing data after the record included, leaves the system, its worlds, and their wrappers unchanged: the file is read into a new system first.

## C++ API

`c_System : c_TidalPyBaseClass, c_TideStateProvider` (`Structures/system/system_.hpp`):

* `add_world(shared_ptr<c_BaseWorld>, is_star, a, e)` (the elements checked by `c_check_orbit`): owns worlds through `shared_ptr` so the C++ system and the Python wrappers co-own the same world, and registers itself as the world's tide-state provider.
* `get_num_worlds`, `get_world(i)`, `find_world_index(name)`.
* `set_tidal_host(i, host_i)`, `clear_tidal_host(i)`, `has_tidal_host(i)`, `get_tidal_host_index(i)`, `get_tidal_host(i)`, `get_tidal_host_mass(i)`, `is_mutual_pair(i)`, and `get_host_orbit(i)` (the elements about the host, merged element by element with a mutual partner's; `c_merge_orbit_elements`).
* `get_tide_state(i, state_out)` and `get_equilibrium_temperature(i)`: the `c_TideStateProvider` interface (`Tides/classes/tide_result_.hpp`). A world reaches it through `c_BaseWorld::get_tide_state(state_out)`; the system clears the world's pointer in its destructor.
* `set_star(i)`, `has_star`, `get_star_index`, `get_star`, `get_star_mass`, `get_star_luminosity`.
* `set/get_semi_major_axis(i)`, `set/get_eccentricity(i)` (orbit about the tidal host).
* `set/get_stellar_semi_major_axis(i)`, `set/get_stellar_eccentricity(i)` (orbit about the star).
* `calc_gravitational_parameter(i)`, `calc_orbital_frequency(i)`, `calc_semi_major_axis_from_frequency(i, n)` (host orbit), and the `calc_stellar_gravitational_parameter(i)` / `calc_stellar_orbital_frequency(i)` counterparts.
* `calc_insolation_flux(i)`, `calc_equilibrium_temperature(i)`.
* `calc_world_evolution(i)` and `calc_system_evolution()`: run the world's tidal solve for the current system state and return the orbital rates, spin rate, and energy terms as a `c_WorldEvolution` struct (every world runs the same `c_BaseWorld::calc_tides` and carries a spin model). `calc_orbital_energy_derivative`, `calc_spin_energy_derivative`, and `calc_energy_residual` compute the energy-balance terms. `c_WorldEvolution::has_tide_model` is false for a rigid world, which is evolved with zero rates and warned about once per world through `TIDALPY_LOG_WARN`.
* `calc_pair_evolution(i)`: dual-body evolution returning a `c_PairEvolution` (both bodies' `c_WorldEvolution` contributions plus the combined shared-orbit rates and energy balance). Built on the shared `calc_dissipation(dissipator_i, companion_mass, n, a, e)` primitive, which computes one body's tidal solve, rate, and spin contribution at its own spin (a body with no tide model is rigid and contributes zero; the primitive itself does not warn). `c_PairEvolution::has_tide_model` is true when either body carries a tide model.
* `c_WorldEvolution` also carries `dU_dM_minus_dw`. `System.evolve` is Python (`TidalPy/Structures/system/evolution.py`) on CyRK's `pysolve_ivp`.
* The binary record rebuilds the world list through `c_world_from_binary` (`Structures/worlds/factory_.hpp`), which peeks each record's `BinaryClassID` and constructs the matching world type; `c_world_kind` gives a loaded world's concrete type so the Cython layer can pick the matching wrapper. It validates the host and star indices, the world names, and every orbit (`c_check_orbit`) before committing, and throws `std::runtime_error` on corrupt data.

## References

* Boué, G., and Efroimsky, M. (2019). Tidal evolution of the Keplerian elements. *Celestial Mechanics and Dynamical Astronomy*, 131(7), 30. The orbital rate equations used by the evolution methods on this container.
* Filippov, A. F. (1988). *Differential Equations with Discontinuous Righthand Sides*. Kluwer Academic Publishers. The combination of the two sides of a tracked equilibrium.
* Makarov, V. V., and Efroimsky, M. (2013). No pseudosynchronous rotation for terrestrial planets and moons. *The Astrophysical Journal*, 764(1), 27. Spin-orbit equilibria of viscoelastic planets.
* Méndez, A., and Rivera-Valentín, E. G. (2017). The equilibrium temperature of planets in elliptical orbits. *The Astrophysical Journal Letters*, 837(1), L1. The orbit-averaged insolation.
* Walterová, M., and Běhounková, M. (2020). Thermal and orbital evolution of low-mass exoplanets. *The Astrophysical Journal*, 900(1), 24. A quasi-static spin in coupled thermal-orbital evolution.
