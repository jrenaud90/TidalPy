# Utilities (`Utilities_x`)

_Updated: 2026-09-16_

`TidalPy.Utilities_x` contains the shared infrastructure the other modules are built on: the base classes that give every object logging, configuration export, and binary serialization; the numerical primitives the inner loops call; and the conversions, constants, and plotting helpers.

Every model class in TidalPy inherits `get_config_dict`, `save_config`, `save_binary`, and `load_binary` from the base classes defined here, along with its schema-version check and its logging.

| Page | Covers |
|---|---|
| [Base Classes](classes_x.md) | The three-level class hierarchy every C++ object inherits, and what each level adds. |
| [Constants](constants.md) | Physical, mathematical, and numerical constants, where they come from, and how they are updated. |
| [Arrays and Interpolation](arrays.md) | The linear interpolation routines used by the tabulated equation of state and the solvers. |
| [Numerics](numerics.md) | Floating-point comparison and the guarded growth functions that turn overflow into visible NaN. |
| [Conversions and Scales](conversions.md) | Unit and orbital-element conversions, and the non-dimensionalization scales the solvers integrate in. |
| [Legendre Polynomials](legendre.md) | Associated Legendre functions and their derivatives for the tidal potential. |
| [Lookup Structures](lookups.md) | Integer-keyed maps for mode-indexed results, usable from C++, Cython, and Python. |
| [Binary Serialization](binary_x.md) | The on-disk binary format, its header, schema versioning, and class ids. |
| [Logging](logging_x.md) | The C++ logging system, its configuration, and how Python and C++ share one set of sinks. |
| [Graphics](graphics_x.md) | Plotting helpers for radial functions, interior profiles, and surface maps. |

```{toctree}
:maxdepth: 1

Base Classes <classes_x.md>
Constants <constants.md>
Arrays and Interpolation <arrays.md>
Numerics <numerics.md>
Conversions and Scales <conversions.md>
Legendre Polynomials <legendre.md>
Lookup Structures <lookups.md>
Binary Serialization <binary_x.md>
Logging <logging_x.md>
Graphics <graphics_x.md>
```

## Where Utilities are Used

The base classes define the contract a new model class has to satisfy; see [Base Classes](classes_x.md) and the "adding a new model" section on any physics module page.

The radial structure and deformation problems are integrated in non-dimensional variables, because the dimensional ones span thirty orders of magnitude and destroy the conditioning of the linear algebra; results are converted back to MKS before they reach the caller. See [Conversions and Scales](conversions.md).
