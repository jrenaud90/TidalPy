# TidalPy.Tides_x Documentation

**TidalPy's Tides Module**

_Note: as of TidalPy v0.7.x this module is named "Tides_x" to not cause overlap with the current TidalPy.tides module.
A future release of TidalPy will replace the old tides module in favor of this one which contains the new C++ / cython
functions. We will refactor it to `TidalPy.Tides` to follow the same capitalization scheme as `TidalPy.RadialSolver`._

[Auto Generated API](https://tidalpy.readthedocs.io/en/latest/API/generated/TidalPy.Tides_x.html)

TidalPy's Tides module contains all functions related to calculating the tidal potential and tidal heating. It includes
obliquity and eccentricity functions which are critical components of the tidal potential. It also provides
functionality to determine which tidal "modes" (forcing frequencies that carry sign) are activate and important for
a given problem.

**As of 0.7.x much of the above functionality is better accessed via the old `TidalPy.tides` module. This module is 
still under development.**

```{toctree}
:maxdepth: 2
:caption: Contents

Obliquity Functions <Obliquity.md>
```