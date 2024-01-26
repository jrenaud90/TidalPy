# Getting Started with TidalPy

## Installation
Please see the package's `readme.md` file for the most up to date installation information.
You can find this file, along with the project's source code, at its [GitHub page](https://github.com/jrenaud90/TidalPy).

## After Installation
At the current stage of `TidalPy` development things are changing too rapidly to make a comprehensive getting started guide, sorry!

For now look to the `Cookbooks` examples on how to use `TidalPy`. If you find any issues or have an idea about a new feature then feel free to leave a new github issue [here](https://github.com/jrenaud90/TidalPy/issues). Stuck on a problem or have a question? Feel free to contact joe.p.renaud@gmail.com try to include a minimal working snippet of code (if relevant to your question).

## TidalPy OOP structure
In addition to basic function calls, TidalPy contains an object-oriented system to quickly build planets, calculate physics, and return results. This system is fairly complicated under the hood, making changes to one object may end up breaking functionality with other objects. This is particularly problemsome in the current version of TidalPy as the `Tests` suite does not fully capture all OOP operations, so you can not count on it catching potential issues. We have created a flow chart to help explain the OOP structure and hopefully mitigate issues. This chart can be found at `Documentation\OOP Structure (<version #>).pdf`.

This flow chart should be updated (via changing the `.svg`) whenever there is a change to TidalPy's OOP structure.

### Basic Modules
The TidalPy package is divided into several modules some of which rely on each other. Below is a basic breakdown of the current modules
- `TidalPy.burnman_interface`: Provides functionality to interface `Burnman` functions into TidalPy.
- `TidalPy.cooling`: Functions related to a planet/layer's cooling (convective, conductive, etc.).
- `TidalPy.dynamic`: Functions related to a planet's orbital and spin evolution.
- `TidalPy.integration_dev`: This is very in-development module that is mostly broken! It will eventually become `TidalPy.simulations` which will handle time integrations.
- `TidalPy.radiogenics`: Functions related to a planet/layer's radiogenic heating.
- `TidalPy.rheology`: Functions related to a planet/layer's rheological properties (complex compliance, viscosity, etc.).
- `TidalPy.stellar`: Functions related to calculating insolation and habitable zones.
- `TidalPy.structures`: The heart of TidalPy's OOP implementation --- various classes for layers and planets.
- `TidalPy.tides`: Functions related to calculating tidal dissipation (using a global approx or a multilayer approach).
- `TidalPy.toolbox`: Helper functions to quickly access various calculations with just a few function calls.
- `TidalPy.utilities`: Various tools used internally inside TidalPy. Generally the user should not need to interact with these unless they are developing new TidalPy functionality.
- `TidalPy.WorldConfigs`: Not a real module, just a location to store planetary configuration files which are used when `TidalPy.build_world` is called.
- `TidalPy`: There are a few top level files that are not contained within a separate module.
    - `configurations.py`: These are important top level configurations that determine how TidalPy runs. The configurations are imported once when TidalPy is initially imported. Generally these configurations should be left alone unless the user is testing out new functionality or wants to change logging behaviour.
    - `constants.py`: Scientific and mathematical constants, largely imported from `SciPy`.
    - `exceptions.py`: Python exceptions used throughout TidalPy. These are all built off of the `TidalPyException` class, so you can easily catch TidalPy exceptions by using `except TidalPyException`.
    - `initialize.py`: Various steps performed when TidalPy is first imported.
    - `io_helper.py`: Functions to help TidalPy interact with the system's directories. 
    - `logger.py`: Functions to help TidalPy record runtime information.
    - `version.py`: Functions to help TidalPy load `version.txt` so it knows what version of TidalPy is being run.

## Functionality still under development
In this early iteration of `TidalPy` you will occasionally come across functions or operations that are not fully implemented or are in a very early alpha state. These will generally be denoted with a `_dev` suffix. Be very careful when using any functionality that is not fully implemented and tested.