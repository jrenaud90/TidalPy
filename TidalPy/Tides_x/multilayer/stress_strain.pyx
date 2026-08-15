# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""On-demand 3D tidal stress/strain/heating - point evaluation.

Thin Cython layer over the C++ strain/stress/heating kernel (kernel_.hpp + strain_radial_.hpp). These
are point evaluations: no grids are materialized. The radial y-functions and complex moduli at the point
come from the (dense) radial solver / EOS, and the tidal potential from a tidal-potential model
(TidalPy.Tides_x.potential.tidal_potential); the caller supplies both.
"""
import numpy as np
from libcpp cimport bool as cpp_bool
cimport numpy as cnp
cnp.import_array()


def volumetric_heating(double complex[::1] stress not None, double complex[::1] strain not None):
    """Volumetric tidal heating [W m-3] from the 6 complex stress and 6 complex strain components.

    Heats once from the (mode-summed) tensors: ``h = |sum_k Im(sigma_k) Re(eps_k) - Re(sigma_k) Im(eps_k)|``
    with a factor 2 on the three off-diagonal components (Europa-book Eq. 42).
    """
    cdef double[12] stress12
    cdef double[12] strain12
    cdef Py_ssize_t k
    for k in range(6):
        stress12[2 * k] = stress[k].real
        stress12[2 * k + 1] = stress[k].imag
        strain12[2 * k] = strain[k].real
        strain12[2 * k + 1] = strain[k].imag
    return c_volumetric_heating_flat(&stress12[0], &strain12[0])


def angular_gram(int degree_l, int order_m):
    """The symmetric 6x6 angular Gram matrix ``G_ij(l, m) = int_0^pi f_i f_j sin(theta) dtheta``.

    The bounded 6-function angular basis is ``f1=P_lm, f2=dP/dtheta, f3=d2P/dtheta2, f4=P/sin,
    f5=-m^2 P/sin^2 + cot dP, f6=(dP - cot P)/sin``. This is the precomputed table backing the analytic
    colatitude collapse (:meth:`LayeredWorld.calc_3d_tides`). Returns a ``(6, 6)`` float64 array; raises
    ``ValueError`` if ``(l, m)`` is outside the tabulated range (``l = 2..10``, ``m = 0..l``).
    """
    cdef double[36] gram36
    if c_angular_gram_flat(degree_l, order_m, &gram36[0]) == 0:
        raise ValueError(f"angular Gram table has no entry for l={degree_l}, m={order_m} "
                         "(supported l=2..10, m=0..l)")
    cdef cnp.ndarray[cnp.float64_t, ndim=2] gram = np.empty((6, 6), dtype=np.float64)
    cdef Py_ssize_t i, j
    for i in range(6):
        for j in range(6):
            gram[i, j] = gram36[i * 6 + j]
    return gram


def strain_stress_heating_point(
        double complex[::1] y not None,
        double complex shear,
        double complex bulk,
        double radius,
        double degree_l,
        cpp_bool is_solid,
        cpp_bool is_incompressible,
        tuple potential6,
        double colatitude):
    """6 complex strains, 6 complex stresses, and the volumetric heating [W m-3] at one point.

    y holds the radial-solver y-functions (only y1..y4 are used). shear/bulk are the COMPLEX
    viscoelastic moduli at this radius and the mode frequency. potential6 is one mode's potential row
    (U + 5 derivatives) from a tidal-potential model's calc_modes. Solid layers only (liquids return NaN).
    """
    cdef double[12] y_ri
    cdef Py_ssize_t k
    for k in range(6):
        if k < y.shape[0]:
            y_ri[2 * k]     = y[k].real
            y_ri[2 * k + 1] = y[k].imag
        else:
            y_ri[2 * k] = 0.0
            y_ri[2 * k + 1] = 0.0
    cdef double[6] pot6
    for k in range(6):
        pot6[k] = potential6[k]
    cdef double[12] strain12
    cdef double[12] stress12
    cdef double heating = 0.0
    c_strain_stress_heating(
        &y_ri[0], shear.real, shear.imag, bulk.real, bulk.imag,
        radius, degree_l, 1 if is_solid else 0, 1 if is_incompressible else 0,
        &pot6[0], colatitude, &strain12[0], &stress12[0], &heating)

    cdef cnp.ndarray[cnp.complex128_t, ndim=1] strain = np.empty(6, dtype=np.complex128)
    cdef cnp.ndarray[cnp.complex128_t, ndim=1] stress = np.empty(6, dtype=np.complex128)
    for k in range(6):
        strain[k] = complex(strain12[2 * k], strain12[2 * k + 1])
        stress[k] = complex(stress12[2 * k], stress12[2 * k + 1])
    return strain, stress, heating
