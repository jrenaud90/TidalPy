# TidalPy Documentation

```{include} Overview/Readme.md
:start-after: TidalPy
:end-before: Overview
```

The pages below cover the package. The source is on [GitHub](https://github.com/jrenaud90/TidalPy).

TidalPy 0.8.0 replaced the Python, Cython, and numba code of TidalPy 0.7.X and earlier with a C++ backend. Users coming from 0.7.X should read the [migration guide](future_structure.md).

The <a href="code_map.html">interactive code map</a> shows the main components of TidalPy, which functions call which, and why.

```{toctree}
:maxdepth: 2
:caption: TidalPy

Overview <Overview/index.md>
Style Guide <Overview/style.md>
Migrating from 0.7.X <future_structure.md>
```

```{toctree}
:maxdepth: 2
:caption: Modules

Structures <Structures/index.md>
Materials <Material/index.md>
Rheology <Rheology/index.md>
Viscosity <Viscosity/index.md>
Partial Melting <PartialMelt/index.md>
Cooling <Cooling/index.md>
Radiogenics <Radiogenics/index.md>
RadialSolver (Love Numbers) <RadialSolver/index.md>
Tides <Tides/index.md>
Dynamics <Dynamics/index.md>
Stellar <Stellar/index.md>
Utilities <Utilities/index.md>
```

```{toctree}
:maxdepth: 2
:caption: Demos

Demos/Basics/B01_config.ipynb
Demos/Basics/B02_world_building.ipynb
Demos/Basics/B03_save_load.ipynb
Demos/Physics/P01_orbits_insolation.ipynb
Demos/Physics/P02_tidal_basics.ipynb
Demos/Physics/P03_rheology_io.ipynb
Demos/Physics/P04_gasgiant_fixedQ_dt.ipynb
Demos/Physics/P05_love_numbers_1d.ipynb
Demos/Physics/P06_tidal_heating_3d.ipynb
Demos/Physics/P07_thermal_eos.ipynb
Demos/Physics/P08_tidal_maps_3d.ipynb
Demos/Physics/P09_bundled_worlds.ipynb
Demos/Physics/P10_thermal_interior.ipynb
Demos/Physics/P11_tidal_truncations.ipynb
Demos/Physics/P12_dynamic_liquids.ipynb
Demos/Physics/P13_melt_and_bulk_dissipation.ipynb
Demos/Physics/P14_seismic_q.ipynb
Demos/Systems/S01_multi_world.ipynb
Demos/Systems/S02_thermal_orbital_evolution.ipynb
Demos/Systems/S03_earth_moon_sun.ipynb
```

```{toctree}
:maxdepth: 2
:caption: Benchmarks

Benchmarks/RadialSolver/Earth_Love_Numbers.ipynb
Benchmarks/RadialSolver/Enceladus_Tobie_Roberts.ipynb
Benchmarks/RadialSolver/Homogeneous_Viscoelastic_Love_Numbers.ipynb
Benchmarks/EOS/EOS_vs_BurnMan.ipynb
Benchmarks/Tides/Renaud2021_Dual_Body_Eccentric.ipynb
Benchmarks/Tides/Hut1981_Constant_Time_Lag.ipynb
Benchmarks/Tides/Exact_Orbit_Tidal_Heating.ipynb
Benchmarks/Performance/Perf_Trends.ipynb
```

```{toctree}
:maxdepth: 2
:caption: Additional Info
:hidden:

Change Log <Changes.md>
```
