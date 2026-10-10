# Utilities (`Utilities`)

_Updated: 2026-09-30_

`TidalPy.Utilities` holds the shared infrastructure of the other modules: the base classes that give every object logging, configuration export (`get_config_dict`, `save_config`), and binary serialization (`save_binary`, `load_binary`); the numerical primitives the inner loops call; and the conversions, constants, and plotting helpers.

| Page | Covers |
|---|---|
| [Base Classes](classes.md) | The three-level class hierarchy every C++ object inherits, and what each level adds. |
| [Constants](constants.md) | Physical, mathematical, and numerical constants, where they come from, and how they are updated. |
| [Arrays and Interpolation](arrays.md) | The linear interpolation routines used by the tabulated equation of state and the solvers. |
| [Numerics](numerics.md) | Floating-point comparison and the guarded growth functions that turn overflow into visible NaN. |
| [Conversions and Scales](conversions.md) | Unit and orbital-element conversions, and the non-dimensionalization scales the solvers integrate in. |
| [Legendre Polynomials](legendre.md) | Associated Legendre functions and their derivatives for the tidal potential. |
| [Lookup Structures](lookups.md) | Integer-keyed maps for mode-indexed results, usable from C++, Cython, and Python. |
| [Binary Serialization](binary.md) | The on-disk binary format, its header, schema versioning, and class ids. |
| [Logging](logging.md) | The C++ logging system, its configuration, and how Python and C++ share one set of sinks. |
| [Graphics](graphics.md) | Plotting helpers for radial functions, interior profiles, and surface maps. |

```{toctree}
:maxdepth: 1

Base Classes <classes.md>
Constants <constants.md>
Arrays and Interpolation <arrays.md>
Numerics <numerics.md>
Conversions and Scales <conversions.md>
Legendre Polynomials <legendre.md>
Lookup Structures <lookups.md>
Binary Serialization <binary.md>
Logging <logging.md>
Graphics <graphics.md>
```

## Where Utilities are Used

The base classes define the contract a new model class has to satisfy (see [Base Classes](classes.md) and the "adding a new model" section of any physics module page). The radial structure and deformation problems are integrated in non-dimensional variables and their results returned in MKS (see [Conversions and Scales](conversions.md)).
