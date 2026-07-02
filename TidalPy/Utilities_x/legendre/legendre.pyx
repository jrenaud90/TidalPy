# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
legendre.pyx
Python/Cython wrapper for TidalPy's associated-Legendre utilities (Utilities_x/legendre).

Both entry points return the triple ``(P_lm(cos theta), dP_lm/dtheta, d2P_lm/dtheta2)`` using the
standard unnormalized associated Legendre functions with the Condon-Shortley phase (matching
``scipy.special.assoc_legendre_p`` with ``branch_cut=2``):

- ``legendre(l, m, colatitude)`` uses the fast precomputed closed-form tables (degrees l = 2..10).
- ``legendre_generic(l, m, colatitude)`` uses the vendored ``xsf`` library and works for any degree
  (a fallback for degrees outside the precomputed range; kept beside the tables so the convention is
  single-sourced).

All angles in radians; colatitude ``theta`` in ``[0, pi]``. Out-of-range orders (m < 0 or m > l)
return NaNs.
"""

from TidalPy.Utilities_x.legendre.legendre cimport (
    c_LegendreValue,
    c_legendre,
    c_legendre_generic,
)


def legendre(int degree_l, int order_m, double colatitude):
    """Associated Legendre triple from the precomputed tables (degrees l = 2..10).

    Parameters
    ----------
    degree_l : int
        Spherical-harmonic degree l (supported: 2..10 for the table path).
    order_m : int
        Order m (0 <= m <= l).
    colatitude : float
        Colatitude theta [radians], in [0, pi].

    Returns
    -------
    tuple of float
        ``(P_lm(cos theta), dP_lm/dtheta, d2P_lm/dtheta2)``. NaNs if m is out of range or l is
        outside the supported table range.
    """
    cdef c_LegendreValue value = c_legendre(degree_l, order_m, colatitude)
    return (value.p, value.dp_dtheta, value.d2p_dtheta2)


def legendre_generic(int degree_l, int order_m, double colatitude):
    """Associated Legendre triple from the generic xsf evaluator (any degree l).

    Parameters
    ----------
    degree_l : int
        Spherical-harmonic degree l (>= 0).
    order_m : int
        Order m (0 <= m <= l).
    colatitude : float
        Colatitude theta [radians], in [0, pi].

    Returns
    -------
    tuple of float
        ``(P_lm(cos theta), dP_lm/dtheta, d2P_lm/dtheta2)``. NaNs if m is out of range.
    """
    cdef c_LegendreValue value = c_legendre_generic(degree_l, order_m, colatitude)
    return (value.p, value.dp_dtheta, value.d2p_dtheta2)
