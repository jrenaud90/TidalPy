# Getting Started with TidalPy

## Installation
Up to date installation instructions can be found
[here](https://tidalpy.readthedocs.io/en/adding-docs/Readme.html#how-to-install) section.

## After Installation
We are deferring the development of a comprehensive "getting started guide" for `TidalPy` until it is closer to a 1.0
version. Until then we recommend looking through the rest of the documentation and the example scripts found in the
`Demos` [folder](https://github.com/jrenaud90/TidalPy/tree/main/Demos).

If you find any issues, have a question, or want to share an idea about a new feature then feel free to leave a new
github issue [here](https://github.com/jrenaud90/TidalPy/issues). TidalPy also has a slack channel for developers and
users. Please contact us at [TidalPy@gmail.com](mailto:TidalPy@gmail.com) if you would like to be invited.

### Package Structure
The TidalPy package is divided into several modules some of which rely on each other.
Below is a basic breakdown of the current modules.

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
- `TidalPy.structures`: The heart of TidalPy's OOP implementation --- various classes for layers and planets.
- `TidalPy.tides`: Functions related to calculating tidal dissipation (using a global approx or a multilayer approach).
- `TidalPy.toolbox`: Helper functions to quickly access various calculations with just a few function calls.
- `TidalPy.utilities`: Various tools used internally inside TidalPy. Generally the user should not need to interact
  with these unless they are developing new TidalPy functionality.
- `TidalPy.WorldPack`: Not a real module, just a location to store planetary configuration files which are used
  when `TidalPy.build_world` is called.
