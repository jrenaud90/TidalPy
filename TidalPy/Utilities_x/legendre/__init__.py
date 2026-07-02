"""TidalPy.Utilities_x.legendre — associated-Legendre utilities (C++/Cython).

Provides the unnormalized associated Legendre functions ``P_lm(cos theta)`` with their first and
second colatitude derivatives (Condon-Shortley phase). ``legendre`` uses fast precomputed closed-form
tables (degrees l = 2..10); ``legendre_generic`` uses the vendored ``xsf`` library for any degree.
"""

from TidalPy.Utilities_x.legendre.legendre import legendre, legendre_generic

__all__ = ["legendre", "legendre_generic"]
