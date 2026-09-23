# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Python wrappers for TidalPy's associated-Legendre utilities.

Both entry points return ``(P_lm(cos theta), dP_lm/dtheta, d2P_lm/dtheta2)`` unnormalized with the
Condon-Shortley phase, matching ``scipy.special.assoc_legendre_p`` with ``branch_cut=2``.
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
        NaNs if m is out of range or l is outside the supported table range.
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
        NaNs if m is out of range.
    """
    cdef c_LegendreValue value = c_legendre_generic(degree_l, order_m, colatitude)
    return (value.p, value.dp_dtheta, value.d2p_dtheta2)
