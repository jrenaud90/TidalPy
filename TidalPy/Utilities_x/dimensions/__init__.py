"""TidalPy.Utilities_x.dimensions - non-dimensionalization utilities (C++/Cython).

Exposes the ``c_NonDimensionalScales`` conversion-scale struct (C++ level) and its Python wrapper
``NonDimensionalScalesClass`` with the ``build_nondimensional_scales`` builder.
"""

from TidalPy.Utilities_x.dimensions.nondimensional import (
    NonDimensionalScalesClass,
    build_nondimensional_scales,
)

__all__ = ["NonDimensionalScalesClass", "build_nondimensional_scales"]
