# TidalPy Demos

A guided, hands-on tour of TidalPy's API. The notebooks are numbered and grouped so a newcomer can
start at the top of `Basics/` and work down through `Physics/` and `Systems/`, with each notebook
building on the ones before it. Each notebook is kept short and focused on one to three concepts.

These demos target the **Python** API. C++ API demos may be added later.

## Order

### Basics
1. `01_config` - the TidalPy configuration system: where the TOML config files live, their format,
   and how to load and edit configuration in a live session with dictionaries.
2. `02_world_building` - the ways to build the different world types, swap layer physics models and
   their parameters, and run basic calculations.
3. `03_save_load` - saving and loading worlds and systems as TOML config files and as binary files,
   and what round-trips through each.

### Physics
4. `04_orbits_insolation` - linking planets, tidal hosts, and stars; equilibrium and insolation
   temperatures across stellar types and orbital distances.
5. `05_tidal_basics` - tidal heating for simple worlds and how it responds to changes in rheology.
6. `06_rheology_io` - a simplified Io: tidal heating across shear modulus and viscosity for the
   Maxwell, Burgers, Andrade, and Sundberg-Cooper rheologies, and across eccentricity.
7. `07_gasgiant_fixedQ_dt` - a short-period gas giant with constant-phase-lag (fixed-Q) and
   constant-time-lag (fixed-dt) models; heating and circularization timescale versus orbital period.
8. `08_love_numbers_1d` - the 1D radial solver for global Love numbers of multi-layer terrestrial
   worlds, and how static vs dynamic and compressible vs incompressible assumptions change them.
9. `09_tidal_heating_3d` - tidal heating and stress from a tidal potential using the 1D radial
   solver for Love numbers.
10. `10_thermal_eos` - changing a world's equation of state and temperature profile, and how each
    reshapes the interior structure.

### Systems
11. `11_multi_world` - constructing systems of multiple worlds, changing tidal hosts, adding and
    removing planets, saving and loading systems, and the instantaneous spin, semi-major axis, and
    eccentricity derivatives.
12. `12_thermal_orbital_evolution` - a simple ODE for a system's coupled thermal and orbital
    evolution (spin, eccentricity, semi-major axis, mantle temperature, tidal heating) integrated
    over roughly a billion years with CyRK.

## Running the notebooks

Use the environment TidalPy is installed into. Figures render inline with `%matplotlib inline`; each
notebook is committed with its outputs so it previews without a running kernel. To pan and zoom
interactively, switch the first cell to `%matplotlib widget` (needs `ipympl`) and rerun.
