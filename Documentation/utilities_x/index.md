# Utilities (`Utilities_x`)

_Updated: 2026-09-13_

`TidalPy.Utilities_x` is the shared infrastructure every other module is built on. It provides is the machinery the physics modules would otherwise each reinvent. The base classes that give every object logging, configuration export, and binary serialization; the numerical primitives the inner loops call; and the conversions, constants, and plotting helpers that sit at the boundary between a calculation and a person reading its result.

It is worth knowing what lives here even if you never import it directly, because its conventions show up everywhere else. Every model class in TidalPy has `get_config_dict`, `save_config`, `save_binary`, and `load_binary`. They are inherited from the base classes defined in this module. Every version check that refuses to load a stale file, every log line, and every non-dimensionalized radius traces back here.

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
| [Graphics](graphics_x.md) | Plotting helpers for radial functions and interior profiles. |

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

## Usage

The base classes are the ones to read first if you plan to add a model of any kind, because they define the contract a new class has to satisfy. See [Base Classes](classes_x.md) and the "adding a new model" section on any physics module page.

The conversion and non-dimensionalization helpers matter when reading solver internals. TidalPy integrates the radial structure and deformation problems in non-dimensional variables, because the dimensional ones span thirty orders of magnitude and destroy the conditioning of the linear algebra. Results are converted back to MKS before they reach the caller.
