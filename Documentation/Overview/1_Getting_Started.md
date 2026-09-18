# Getting Started with TidalPy

## Installation

```{include} Readme_raw.md
:start-after: How to Install
:end-before: Using TidalPy
```

## After Installation
A comprehensive "getting started guide" is deferred until `TidalPy` is closer to a 1.0 release. Until then we recommend the rest of this documentation and the example scripts in the `Demos` [folder](https://github.com/jrenaud90/TidalPy/tree/main/Demos).

If you find an issue, have a question, or want to share an idea for a new feature, please open a GitHub issue [here](https://github.com/jrenaud90/TidalPy/issues). TidalPy also has a slack channel for developers and users. Contact us at [TidalPy@gmail.com](mailto:TidalPy@gmail.com) if you would like to be invited.

### Package Structure
TidalPy is divided into several modules, some of which rely on each other.

- `TidalPy.Extending`: Provides support for 3rd party packages.
- `TidalPy.cooling`: Functions related to a planet/layer's cooling (convective, conductive, etc.).
- `TidalPy.dynamic`: Functions related to a planet's orbital and spin evolution.
    - Read more about the dynamics module [here](https://tidalpy.readthedocs.io/en/latest/Dynamics/index.html).
- `TidalPy.radiogenics`: Functions related to a planet/layer's radiogenic heating.
- `TidalPy.RadialSolver` : Functions related to solving for a planet's Love numbers.
    - Read more about the RadialSolver module [here](https://tidalpy.readthedocs.io/en/latest/RadialSolver/index.html)
- `TidalPy.rheology`: Functions related to a planet/layer's rheological properties (complex shear, viscosity, etc.).
    - Read more about the rheology module [here](https://tidalpy.readthedocs.io/en/latest/Rheology/index.html)
- `TidalPy.stellar`: Functions related to calculating insolation and habitable zones.
- `TidalPy.structures`: TidalPy's object-oriented implementation: classes for layers and planets.
- `TidalPy.tides`: Functions related to calculating tidal dissipation (using a global approx or a multilayer approach).
- `TidalPy.toolbox`: Helper functions to quickly access various calculations with just a few function calls.
- `TidalPy.utilities`: Tools used internally by TidalPy. You should not need to interact with these unless you are developing new TidalPy functionality.
- `TidalPy.WorldPack`: Not a module: a location for the planetary configuration files used when `TidalPy.build_world` is called.
