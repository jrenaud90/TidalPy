# TidalPy Demos

Tutorial notebooks for TidalPy's API. We recommend starting at the top of `Basics/` and working your way down through `Physics/` and `Systems/`, with each notebook building on the ones before it. Each notebook covers one to three concepts.

These demos use the Python API. TidalPy's C++ API is exposed through headers that other packages can include, and it has no demonstrations here.

## Order

### Basics
- `B01_config`: the TidalPy configuration system: where the world, system, and material TOML files live and how their installed copies behave, their format, loading and editing a config in a live session and building a world from a name, a path, or a dictionary, the stated and solved mass, validation errors, and the global configuration (`TidalPy.config`, `reinit`, `save_config`).
- `B02_world_building`: building each world type (terrestrial, gas giant, star) and their shared bulk calculations; solving and plotting a terrestrial world's interior structure and moment of inertia; and a star's insolation and a planet's equilibrium temperature.
- `B03_save_load`: saving and loading worlds and systems as TOML config files and as binary files (`load_world`, `load_system`), copying a world in memory, and what round-trips through each.

### Physics
- `P01_orbits_insolation`: linking a planet to its star in a system; the orbital frequency, the orbit-averaged insolation, and the equilibrium temperature; and the conservative habitable zone of Kopparapu et al. (2013) across stellar types.
- `P02_tidal_basics`: fixed-Q tidal heating of Io in a system; heating against eccentricity with two truncation levels and the exact eccentricity functions; the tides of non-synchronous rotation and obliquity; and heating against $Q$ next to Io's observed output.
- `P03_rheology_io`: a uniform incompressible Io, with tidal heating across shear modulus and viscosity for the Maxwell, Burgers, Andrade, and Sundberg-Cooper rheologies, and each rheology's attenuation across forcing frequency.
- `P04_gasgiant_fixedQ_dt`: a short-period Saturn-mass planet around a Sun-like star with constant-phase-lag (fixed-Q) and constant-time-lag (fixed-dt) models; heating and circularization timescale versus orbital period, and their power-law slopes.
- `P05_love_numbers_1d`: the 1D radial solver for global Love numbers of multi-layer terrestrial worlds, and how static vs dynamic and compressible vs incompressible assumptions change them; Love numbers solved on a built world with each Love method; tidal, loading, and free Love numbers; the radial functions, the solution diagnostics, and $k_2$ and its quality factor across forcing periods from minutes (the normal modes) to past the Maxwell peak.
- `P06_tidal_heating_3d`: the secular 3D tidal heating of a homogeneous Io, with profiles in colatitude and depth, a longitude-resolved surface map, the power per unit radius, and the volume integral checked against the 1D total.
- `P07_thermal_eos`: changing a world's equation of state and temperature profile, how each reshapes the interior structure, and how melting changes $k_2$ and the dissipation; a solved mass that drifts from the stated one, and a density inversion from melt in the density.
- `P08_tidal_maps_3d`: maps of the 3D tidal response of a homogeneous Io, covering secular heating across eccentricity, spin, and obliquity; stress through an orbit; the largest tension over an orbit; stress-strain loops and where dissipation appears; the surface displacement; and the same grids timed on one thread and on every core.
- `P09_bundled_worlds`: a tour of the bundled worlds, checking each solved interior's mass, moment of inertia, and Love number against the observations it was fitted to or predicts; the measured moments of inertia with their uncertainties; Jupiter's layered interior and Juno's $k_2$, Pluto's ocean, Io's heat budget and where Jupiter's share comes from, and the TRAPPIST-1 planets at the eccentricities of Agol et al. (2021).
- `P10_thermal_interior`: the cooling and radiogenic model families, convection and conduction across the critical Rayleigh number, radiogenic heating including the short-lived isotopes, the temperature and heat-flow profile carried through the equation-of-state solve, tidal heating by layer, and a simple heat budget.
- `P11_tidal_truncations`: the eccentricity and obliquity truncation levels; heating and $\dot{e}$ against eccentricity at each level and with the exact eccentricity functions, and the error of each; how many modes the exact functions need; the cost of each level with closed-form and radial-solver Love numbers; and the helpers that pick a level, with the measured limit of each.
- `P12_dynamic_liquids`: static and dynamic liquid layers across forcing periods; why a constant-density liquid's dynamic solve breaks down at long periods (its stratification) and the warning TidalPy logs; incompressible liquids and the `_dynamic` bundled worlds, whose liquids follow a pressure-dependent law; the tolerance a large dynamic core needs at long periods; and capturing TidalPy's warnings with `TidalPy.capture_log`.
- `P13_melt_and_bulk_dissipation`: the Zener and Maxwell rheologies; a material point through its melting range with each of the partial-melt switches (melt in the density, the bulk modulus, and a compaction bulk viscosity); Io with a partially molten asthenosphere, against the $k_2$ and $Q$ Juno measured; why a Maxwell bulk rheology misleads and when a bulk response matters; and what mixing the melt into the density does to the structure.
- `P14_seismic_q`: the seismic Q rheology, which builds complex moduli from quality factors instead of a viscosity, and the dispersion it carries from seismic to tidal periods; PREM's quality factors; PREM with its own attenuation (`earth_prem_q`, `q_provided`) against the elastic PREM at 1 s and at tidal periods; the Q frequency exponent that reproduces Earth's measured M2 body-tide lag; and per-layer settings.

### Systems
- `S01_multi_world`: constructing systems of multiple worlds with synchronous rotation, the instantaneous spin, semi-major axis, and eccentricity derivatives and their timescales, a moon added to a planet as its tidal host, saving and loading systems as TOML and binary files, and dropping a world by rebuilding the system.
- `S02_thermal_orbital_evolution`: the coupled thermal and orbital evolution of the bundled `earth_thermal` about a star over 5 Gyr with `System.evolve`: the orbit, the spin captured into and released from spin-orbit resonances, and the layer temperatures, melting, and tidal heating.
- `S03_earth_moon_sun`: a star that is nobody's tidal host, a mutual Earth-Moon pair with single-body and dual-body evolution, the orbit about the star used for insolation, obliquity tides, and how this solid-body model's rates compare with lunar laser ranging.

## Running the Notebooks

Run the notebooks in the environment TidalPy is installed into. Figures render inline with `%matplotlib inline`, and each notebook is committed with its outputs so it previews without a running kernel. To pan and zoom interactively, switch the first cell to `%matplotlib widget` (requires `ipympl`) and rerun.
