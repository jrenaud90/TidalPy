"""Unnormalized associated Legendre functions with their first and second colatitude derivatives.

Provides the unnormalized associated Legendre functions ``P_lm(cos theta)`` with their first and second
colatitude derivatives (Condon-Shortley phase). ``legendre`` uses the precomputed tables (l = 2..10) and
``legendre_generic`` the degree and order recurrences for any degree.
"""

from TidalPy.Utilities_x.legendre.legendre import legendre, legendre_generic

__all__ = ["legendre", "legendre_generic"]
