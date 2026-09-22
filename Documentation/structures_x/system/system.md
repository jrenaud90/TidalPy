# System (`structures_x.system`)

_Updated: 2026-09-21_

A `System` links two or more worlds (a star, planets, moons) into a gravitationally bound group. It tracks two roles independently:

* a tidal host per world: the body that raises a tide on that world. Every world names its own host, or none.
* the star: the body that provides insolation to each world in the system.

The star need not be a world's tidal host. In the Earth-Moon-Sun system the Moon's tidal host is the Earth, the Earth's can be the Moon, and the star is the Sun. For an exoplanet orbiting its star, the star is also its tidal host, and the orbits about the host and about the star coincide. Each world therefore carries two two-body orbits: one about its tidal host and one about the star. A world interacts only with its tidal host and the star. The `System` is the container on which orbital evolution and insolation are computed.

TidalPy has no N-body dynamics: its dynamics are semi-analytic two-body rates. A problem with three or more mutually interacting bodies can be built from TidalPy worlds but needs external scripting to account for the additional forcing.

## Building a `System`

```python
from TidalPy.constants import au
from TidalPy.structures_x.system import System
from TidalPy.structures_x.worlds.stellar import StarWorld
from TidalPy.structures_x.worlds.layered import LayeredWorld

system = System("sol")
star = StarWorld("sun", 6.957e8, 1.988e30)
earth = LayeredWorld("earth", 6.371e6, 5.972e24)

system.add_world(
    star,
    is_star=True)                 # The star is tidally forced by nothing here, so it names no host
system.add_world(
    earth,
    tidal_host=star,              # Exoplanet case: the star is also the world's tidal host
    semi_major_axis=au,
    eccentricity=0.0167)
```

`add_world(world, tidal_host=None, is_star=False, semi_major_axis=None, eccentricity=0.0)` returns the world's index. `tidal_host` is a world already in the system, given by index, name, or object; `set_tidal_host(world, tidal_host)` names or changes it afterwards, which is how a host added after the worlds it hosts is named. There is no system-wide host: a world with no tidal host is not tidally forced and is skipped by the evolution methods. `is_star` marks the insolation source, and the last world flagged wins. `semi_major_axis` and `eccentricity` describe the orbit about the tidal host; the orbit about the star is set separately (see below). The system co-owns each world with its Python wrapper, so the wrapper you passed stays usable and is the same object the system hands back.

For a system where the star is a separate body from a world's tidal host (_e.g._, Earth-Moon-Sun):

```python
system = System()
system.add_world(
    sun,
    is_star=True)                                # The star, and nobody's tidal host
system.add_world(earth)                          # The Moon's tidal host, named below
system.add_world(
    moon,
    tidal_host=earth,                            # The Earth raises the Moon's tides
    semi_major_axis=3.844e8,                     # Moon about the Earth (tidal)
    eccentricity=0.0549)
system.set_tidal_host(earth, moon)               # The Moon raises the Earth's tides in turn
system.set_stellar_semi_major_axis("moon", au)   # Moon about the Sun (insolation)
system.set_stellar_eccentricity("moon", 0.0167)
```

### Mutual Pairs

Two worlds that host each other share one orbit, so its elements need to be given on only one of them: the member that carries no semi-major axis of its own takes its partner's elements. When both carry elements they must agree, and a disagreement raises `ValueError` from any method that reads the orbit. `is_mutual_pair(world)` reports the relation. `calc_system_evolution` gives each member its own entry, and their contributions to the shared orbit add, which is what `calc_pair_evolution` returns for either member.

## Building a `System` from TOML (`build_system`)

A whole system can be described in TOML and built in one call, mirroring `build_world`. Each `[worlds.<name>]` table names a `world` (a bundled world name, a path to a world TOML, or an inline world config) plus its `tidal_host` (the table key of the world that raises its tides), its star role, and its orbital elements. `system.get_config_dict()` returns the same schema for the system as it stands now (each world inlined with its own live `get_config_dict()`), and `build_system_from_dict(config)` rebuilds the system from it. A host may be declared after the worlds it hosts. A world that states `semi_major_axis_m` or `eccentricity` must name a `tidal_host` for them to be about. The table key becomes the world's name within the system, so a bundled world template can be reused under different names.

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
eccentricity      = 0.0167
stellar_semi_major_axis_m = 1.495978707e11  # orbit about the star (for insolation)
stellar_eccentricity      = 0.0167
```

```python
from TidalPy.structures_x.configs import build_system

system = build_system("sol_system")      # a bundled system name, a .toml path, or a dict
# equivalently: System.build("sol_system")
system.calc_insolation_flux("earth")     # ~1361 W/m^2 (the solar constant)
```

`build_system(source, force=False)`, a thin wrapper over `System.build` that mirrors `build_world` and `BaseWorld.build`, resolves the source, validates it (schema version and structure), and builds each member world with `build_world`. `construct_system(config)` does the same from an already-parsed `dict`. To make the star and the tidal host different bodies, set `is_host` and `is_star` on different worlds (_e.g._, a moon whose tidal host is its planet but whose insolation comes from the system star), and give each world both a tidal-host orbit (`semi_major_axis_m` and `eccentricity`) and a stellar orbit (`stellar_semi_major_axis_m` and `stellar_eccentricity`).

A built system retains its normalized configuration on `source_config` and can be written back out:

```python
system.save_to_toml("my_system.toml")
system.get_config_dict()
```

`save_to_toml` writes the retained `source_config` when present (the original world references). For a system assembled directly in Python it falls back to `get_config_dict`, the self-contained expansion that inlines each world's full config together with its roles and orbital elements. Each inlined world is its live `get_config_dict()`, which is builder-valid, so the expansion rebuilds through `build_system`.

## Worlds, Host, and Identification

Most methods accept a world by index (int, negatives allowed), by name (str), or as the world object itself:

```python
system.get_tidal_host("earth")        # the world that raises Earth's tides (or None)
system.get_tidal_host_index("earth")  # its index, or -1
system.has_tidal_host("sun")          # False: nothing forces the star here
system.set_tidal_host("earth", "sun") # name a host by index / name / object (None removes it)
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

The mean motion follows Kepler's third law using the combined host and world mass. `calc_gravitational_parameter` returns NaN for a world with no tidal host (it has no orbit within the system); `calc_orbital_frequency` also returns NaN for an unset or non-positive semi-major axis.

## Star and Insolation

The star is designated with `is_star` (or `set_star`) and drives insolation. Each world carries its own orbit about the star, independent of the tidal-host orbit:

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

where the factor $\sqrt{1-e^{2}}$ comes from the time average of $1/r^{2}$ over the eccentric orbit, $\langle 1/r^{2} \rangle = 1/(a^{2}\sqrt{1-e^{2}})$ (Méndez and Rivera-Valentín 2017). `calc_equilibrium_temperature` applies the world's own gray-body radiative balance to that flux,

$$T = \left(\frac{(1-A)\,F}{4\,\varepsilon\,\sigma}\right)^{1/4},$$

with the world's albedo $A$, its emissivity $\varepsilon$, and the Stefan-Boltzmann constant $\sigma$. Both raise `RuntimeError` if no star is set and return NaN for the star's own entry, an unset stellar semi-major axis, or a star with no luminosity.

## Orbital and Spin Evolution

A layered world whose tide model is `rheology` (the default for a terrestrial world) takes its Love numbers from its interior, so run `world.solve_eos()` on every such member before any of the evolution methods below; they raise `RuntimeError` otherwise. A later `solve_eos` retires the world's tidal result, and the next evolution call solves it again.

`calc_world_evolution(world)` evolves a single world. It solves the world's global tides in the current system state (mean motion from Kepler's third law, spin and obliquity from the world, eccentricity and semi-major axis from the orbit about its tidal host, host mass from that host), then turns the tidal-potential derivatives into the orbital rates and the world's spin rate. Only this world raises tides; its host is treated as a point mass. A world with no tidal host, or no usable orbit about it, comes back with `evolved = False`.

A world that belongs to a system can be asked for the same state directly: `world.get_tide_state()` returns it as a dict in the argument order of `calc_tides`, or `None` for a world outside a system, with no tidal host, or with no usable orbit. Orbital state is never stored on a world; the system supplies it on request and stops doing so when it is deleted.

```python
ev = system.calc_world_evolution("moon")
ev["da_dt"], ev["de_dt"], ev["dn_dt"]     # orbital rates [m/s], [1/s], [rad/s^2]
ev["dspin_dt"]                            # spin rate [rad/s^2] (0 without a spin model)
ev["tidal_heating"]                       # [W]
ev["energy_residual"]                     # heating + dE_orbit/dt + dE_spin/dt (~0 under conservation)
```

The returned dict also carries the state used (`orbital_frequency`, `semi_major_axis`, `eccentricity`, `spin_frequency`, `host_mass`, `target_mass`), the raw tidal outputs (`dU_dM`, `dU_dw`, `dU_dO`), the `moment_of_inertia` and `has_spin` flag, and the energy terms (`dE_orbit_dt`, `dE_spin_dt`). `evolved` is `False` for the host's own entry or a world with no usable orbit about the host; its rates are then zero. `calc_system_evolution()` returns one such dict per world, in index order.

The rates follow the orbital rate engine (`dynamics_x`). With the tidal-potential derivatives $\partial U/\partial X$ of the [global tides](../../Tides_x/global_tides.md) converted to disturbing-function derivatives $\partial\mathcal{R}/\partial X = -\frac{M_{w} + M_{h}}{M_{w}}\,\partial U/\partial X$, for the world mass $M_{w}$ and the host mass $M_{h}$,

$$\frac{da}{dt} = \frac{2}{na}\,\frac{\partial\mathcal{R}}{\partial\mathcal{M}}, \qquad \frac{de}{dt} = \frac{\sqrt{1-e^{2}}}{na^{2}e}\left(\sqrt{1-e^{2}}\,\frac{\partial\mathcal{R}}{\partial\mathcal{M}} - \frac{\partial\mathcal{R}}{\partial\varpi}\right), \qquad \frac{dn}{dt} = -\frac{3}{2}\,\frac{n}{a}\,\frac{da}{dt},$$

with $de/dt = 0$ at $e = 0$. The spin rate comes from the world's attached spin model, $\ddot{\theta} = (M_{h}/C)\,\partial U/\partial\Omega$, with $C$ the polar moment of inertia. The heating and the orbit and spin energy loss balance,

$$\dot{E} = -\left(\frac{dE_\mathrm{orbit}}{dt} + \frac{dE_\mathrm{spin}}{dt}\right), \qquad E_\mathrm{orbit} = -\frac{G M_{h} M_{w}}{2a}, \qquad E_\mathrm{spin} = \frac{1}{2}\,C\,\dot{\theta}^{2}.$$

Each world evolves on its own two-body orbit about its tidal host and dissipates independently.

### Dual-Body Dissipation

`calc_pair_evolution(world)` evolves a world together with its own tidal host, with both bodies raising a tide on their shared orbit. Each body's tides are solved with the other body as the tide raiser (masses swapped), so their orbital-rate contributions add and each body evolves its own spin:

```python
pair = system.calc_pair_evolution("moon")
pair["da_dt"], pair["de_dt"], pair["dn_dt"]     # combined shared-orbit rates
pair["tidal_heating_total"]                     # heating in both bodies
pair["energy_residual"]                         # ~0: heating_total + dE_orbit/dt + dE_spin_total/dt
pair["world"]                                   # the orbiting world's single-body contribution (a dict)
pair["host"]                                    # the host's single-body contribution (a dict)
```

The `world` and `host` entries are each a full `calc_world_evolution`-style dict (their own solve, spin, heating, and share of the orbital rates and energy). The combined balance is the sum of the two single-body balances, $\dot{E}_{w} + \dot{E}_{h} = -\left(dE_\mathrm{orbit}/dt + dE_{\mathrm{spin},w}/dt + dE_{\mathrm{spin},h}/dt\right)$. A body with no tide model attached is rigid and contributes nothing, so a rigid host reduces `calc_pair_evolution` to the single-body `calc_world_evolution` result.

## Serialization

A system serializes to and loads from TidalPy's binary format, reconstructing its heterogeneous world list: each world's concrete type (star, layered, gas giant) is recovered from the stream.

```python
system.save_binary("system.tpyb")

loaded = System()
loaded.load_binary("system.tpyb")
loaded["earth"]        # comes back as a LayeredWorld (with its layers), the Sun as a StarWorld, ...
```

`System` inherits the binary machinery (`save_binary`, `load_binary`, `get_schema_version_str`, `save_config`) from the shared `TidalPyBaseClass`. As for a directly loaded world, the sub-models a world does not serialize (the layer EOS profile data and the tide, spin, and luminosity models) are not carried in the binary and are reattached after load; the container state (name, the star role, and each world's tidal host and orbital elements about both that host and the star) and each world's own fields (including, for a star, its effective temperature, so insolation survives) are.

## C++ API

`c_System : c_TidalPyBaseClass, c_TideStateProvider` (`structures_x/system/system_.hpp`):

* `add_world(shared_ptr<c_BaseWorld>, is_star, a, e)`: owns worlds through `shared_ptr` so the C++ system and the Python wrappers co-own the same world, and registers itself as the world's tide-state provider.
* `get_num_worlds`, `get_world(i)`, `find_world_index(name)`.
* `set_tidal_host(i, host_i)`, `clear_tidal_host(i)`, `has_tidal_host(i)`, `get_tidal_host_index(i)`, `get_tidal_host(i)`, `get_tidal_host_mass(i)`, `is_mutual_pair(i)`, and `get_host_orbit(i)` (the elements about the host, a mutual partner's when the world carries none).
* `get_tide_state(i, state_out)` and `get_equilibrium_temperature(i)`: the `c_TideStateProvider` interface (`Tides_x/classes/tide_result_.hpp`). A world reaches it through `c_BaseWorld::get_tide_state(state_out)`; the system clears the world's pointer in its destructor.
* `set_star(i)`, `has_star`, `get_star_index`, `get_star`, `get_star_mass`, `get_star_luminosity`.
* `set/get_semi_major_axis(i)`, `set/get_eccentricity(i)` (orbit about the tidal host).
* `set/get_stellar_semi_major_axis(i)`, `set/get_stellar_eccentricity(i)` (orbit about the star).
* `calc_gravitational_parameter(i)`, `calc_orbital_frequency(i)`, `calc_semi_major_axis_from_frequency(i, n)` (host orbit), and the `calc_stellar_gravitational_parameter(i)` / `calc_stellar_orbital_frequency(i)` counterparts.
* `calc_insolation_flux(i)`, `calc_equilibrium_temperature(i)`.
* `calc_world_evolution(i)` and `calc_system_evolution()`: run the world's tidal solve for the current system state and return the orbital rates, spin rate, and energy terms as a `c_WorldEvolution` struct (a layered world is resolved through `dynamic_cast` so the rheology solve runs and the spin model is reached; a layerless world uses the analytic solve and contributes no spin). `calc_orbital_energy_derivative`, `calc_spin_energy_derivative`, and `calc_energy_residual` compute the energy-balance terms.
* `calc_pair_evolution(i)`: dual-body evolution returning a `c_PairEvolution` (both bodies' `c_WorldEvolution` contributions plus the combined shared-orbit rates and energy balance). Built on the shared `calc_dissipation(dissipator_i, companion_mass, n, a, e)` primitive, which computes one body's tidal solve, rate, and spin contribution (a body with no tide model is rigid and contributes zero).
* `write_binary` and `read_binary` (the inherited `save_binary` and `load_binary` open the file). `read_binary` rebuilds the world list through `c_world_from_binary` (`structures_x/worlds/factory_.hpp`), which peeks each record's `BinaryClassID` and constructs the matching world type; `c_world_kind` gives a loaded world's concrete type so the Cython layer can pick the matching wrapper.

## References

* Boué, G., and Efroimsky, M. (2019). Tidal evolution of the Keplerian elements. *Celestial Mechanics and Dynamical Astronomy*, 131(7), 30. The orbital rate equations used by the evolution methods on this container.
* Méndez, A., and Rivera-Valentín, E. G. (2017). The equilibrium temperature of planets in elliptical orbits. *The Astrophysical Journal Letters*, 837(1), L1. The orbit-averaged insolation.
