"""TidalPy.Utilities_x.arrays — array helper utilities (C++/Cython).

Exposes ``interp``, a ``numpy.interp``-style 1-D linear interpolation backed by a
header-only C++ implementation (``interp_.hpp``).
"""

from TidalPy.Utilities_x.arrays.interp import interp

__all__ = ["interp"]
