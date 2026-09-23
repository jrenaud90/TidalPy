# TidalPy Demos

Tutorial notebooks for TidalPy's API, numbered and grouped: start at the top of `Basics/` and work down through `Physics/` and `Systems/`, with each notebook building on the ones before it. Each notebook covers one to three concepts.

These demos target the Python API. C++ API demos may be added later.

## Order

### Basics
1. `01_config`: the TidalPy configuration system, where the TOML config files live, their format, and how to load and edit configuration in a live session with dictionaries.
2. `02_world_building`: building the different world types, swapping layer physics models and their parameters, and running basic calculations.
3. `03_save_load`: saving and loading worlds and systems as TOML config files and as binary files, and what round-trips through each.

### Physics
4. `04_orbits_insolation`: linking planets, tidal hosts, and stars; equilibrium and insolation temperatures across stellar types and orbital distances.
5. `05_tidal_basics`: tidal heating for simple worlds and how it changes with rheology.
6. `06_rheology_io`: a simplified Io, with tidal heating across shear modulus and viscosity for the Maxwell, Burgers, Andrade, and Sundberg-Cooper rheologies, and across eccentricity.
7. `07_gasgiant_fixedQ_dt`: a short-period gas giant with constant-phase-lag (fixed-Q) and constant-time-lag (fixed-dt) models; heating and circularization timescale versus orbital period.
8. `08_love_numbers_1d`: the 1D radial solver for global Love numbers of multi-layer terrestrial worlds, and how static vs dynamic and compressible vs incompressible assumptions change them; Love numbers solved on a built world with each Love method; tidal, loading, and free Love numbers; the radial functions, the solution diagnostics, and the frequency dependence of $k_2$.
9. `09_tidal_heating_3d`: the secular 3D tidal heating of a homogeneous Io, with profiles in colatitude and depth, a longitude-resolved surface map, and the volume integral checked against the 1D total.
10. `10_thermal_eos`: changing a world's equation of state and temperature profile, and how each reshapes the interior structure.

14. `14_bundled_worlds`: a tour of the bundled worlds, checking each solved interior's mass, moment of inertia, and Love number against the observations it was fitted to or predicts; Jupiter's layered interior, Pluto's ocean, Io's heat budget, and the TRAPPIST-1 planets.
15. `15_thermal_interior`: the cooling and radiogenic model families, the temperature and heat-flow profile carried through the equation-of-state solve, tidal heating by layer, and a simple heat budget.

### Systems
11. `11_multi_world`: constructing systems of multiple worlds, changing tidal hosts, adding and removing planets, saving and loading systems, and the instantaneous spin, semi-major axis, and eccentricity derivatives.
12. `12_thermal_orbital_evolution`: a simple ODE for a system's coupled thermal and orbital evolution (spin, eccentricity, semi-major axis, mantle temperature, tidal heating) integrated over roughly a billion years with CyRK.
16. `16_earth_moon_sun`: a star that is nobody's tidal host, a mutual Earth-Moon pair with single-body and dual-body evolution, the orbit about the star used for insolation, and obliquity tides.

### Three-Dimensional Maps
13. `13_tidal_maps_3d` (in `Physics/`): maps of the 3D tidal response of a homogeneous Io, covering secular heating across eccentricity, spin, and obliquity; stress through an orbit; the largest tension over an orbit; stress-strain loops and where dissipation appears; the surface displacement; and the same grids timed on one thread and on every core.

## Running the Notebooks

Run the notebooks in the environment TidalPy is installed into. Figures render inline with `%matplotlib inline`, and each notebook is committed with its outputs so it previews without a running kernel. To pan and zoom interactively, switch the first cell to `%matplotlib widget` (requires `ipympl`) and rerun.
