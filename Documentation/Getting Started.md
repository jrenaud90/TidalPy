# Getting Started with TidalPy

## Installation
Please see the package's `readme.md` file for the most up to date installation information.
You can find this file, along with the project's source code, at its [GitHub page](https://github.com/jrenaud90/TidalPy).

## After Installation
At the current stage of `TidalPy` development things are changing too rapidly to make a comprehensive getting started guide, sorry!

For now look to the `Demos` examples on how to use `TidalPy`. If you find any issues or have an idea about a new feature then feel free to leave a new github issue [here](https://github.com/jrenaud90/TidalPy/issues). Stuck on a problem or have a question? Feel free to contact joe.p.renaud@gmail.com try to include a minimal working snippet of code (if relevant to your question).

## TidalPy OOP structure
In addition to basic function calls, TidalPy contains an object-oriented system to quickly build planets, calculate physics, and return results. This system is fairly complicated under the hood, making changes to one object may end up breaking functionality with other objects. We have created a flow chart to help explain the OOP structure and hopefully mitigate issues. This chart can be found at `Documentation\OOP Structure (<version #>).pdf`.

This flow chart should be updated (via changing the `.svg`) whenever there is a change to TidalPy's OOP structure.

### Basic Modules
The TidalPy package is divided into several modules some of which rely on each other. Below is a basic breakdown of the current modules.

- `TidalPy.Extending`: Provides support for 3rd party packages.
- `TidalPy.cooling`: Functions related to a planet/layer's cooling (convective, conductive, etc.).
- `TidalPy.dynamic`: Functions related to a planet's orbital and spin evolution.
- `TidalPy.radiogenics`: Functions related to a planet/layer's radiogenic heating.
- `TidalPy.rheology`: Functions related to a planet/layer's rheological properties (complex shear, viscosity, etc.).
- `TidalPy.stellar`: Functions related to calculating insolation and habitable zones.
- `TidalPy.structures`: The heart of TidalPy's OOP implementation --- various classes for layers and planets.
- `TidalPy.tides`: Functions related to calculating tidal dissipation (using a global approx or a multilayer approach).
- `TidalPy.toolbox`: Helper functions to quickly access various calculations with just a few function calls.
- `TidalPy.utilities`: Various tools used internally inside TidalPy. Generally the user should not need to interact with these unless they are developing new TidalPy functionality.
- `TidalPy.WorldPack`: Not a real module, just a location to store planetary configuration files which are used when `TidalPy.build_world` is called.

## Functionality still under development
In this early iteration of `TidalPy` you will occasionally come across functions or operations that are not fully implemented or are in a very early alpha state. These will generally be denoted with a `_dev` suffix. Be very careful when using any functionality that is not fully implemented and tested.
