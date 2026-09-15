# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Point-wise 3D tidal strain, stress, heating, and displacement kernels.

A thin Cython layer over the C++ kernel (``kernel_.hpp`` and ``strain_radial_.hpp``) that the world's 3D methods also
use. These are point evaluations and materialize no grids. The caller supplies the radial functions and complex
moduli at the point (``LayeredWorld.get_love_radial_y``, ``calc_complex_shear_modulus``, and
``calc_complex_bulk_modulus`` after a radial-solver Love solve at the mode's degree and ``|frequency|``) and one tidal
mode's potential row from ``TidalPy.Tides_x.potential.tidal_potential_3d_modes``.

Potential rows are complex phasor amplitudes ``(U, dU/dtheta, dU/dphi, d2U/dtheta2, d2U/dphi2, d2U/dtheta_dphi)``
with ``U(t) = Re[U_c e^{i omega t}]``; a real row is a phasor with zero phase. The strains, stresses, and
displacements returned are complex amplitudes in the same convention. To assemble several modes:

- Evaluate the moduli and radial functions at ``|frequency|`` and conjugate the row of any mode whose frequency is
  negative, so every amplitude is at ``+|frequency|``.
- Sum the strain and stress amplitudes of all modes that share ``|frequency|`` and pass the sums to
  :func:`volumetric_heating`. ``|frequency| / 2`` times its result is that frequency's cycle-averaged heating
  [W m-3], and the heating of different frequencies adds.
- A field at time ``t`` is ``Re[amplitude e^{i |frequency| t}]`` summed over the modes.
"""
import numpy as np
from libcpp cimport bool as cpp_bool
cimport numpy as cnp
cnp.import_array()


cdef int cy_potential_row_to_flat(object potential6, double* potential12) except -1:
    """Copy one mode's potential row into 12 doubles (real, imaginary per entry), keeping the imaginary parts.

    Parameters
    ----------
    potential6 : array-like of complex or float
        ``(U, dU/dtheta, dU/dphi, d2U/dtheta2, d2U/dphi2, d2U/dtheta_dphi)``.
    potential12 : double*
        Output buffer of 12 doubles.

    Raises
    ------
    ValueError
        If the row does not hold exactly 6 values.
    """
    row = np.asarray(potential6, dtype=np.complex128)
    if row.ndim != 1 or row.shape[0] != 6:
        raise ValueError(
            f"potential6 must hold 6 values (U and its five angular derivatives); got shape {row.shape}.")
    cdef double complex[::1] row_view = np.ascontiguousarray(row)
    cdef Py_ssize_t k
    for k in range(6):
        potential12[2 * k] = row_view[k].real
        potential12[2 * k + 1] = row_view[k].imag
    return 0


def volumetric_heating(double complex[::1] stress not None, double complex[::1] strain not None):
    """Magnitude of the weighted bilinear form of 6 complex stress and 6 complex strain components.

    Parameters
    ----------
    stress : numpy.ndarray of complex128, shape (6,)
        Stress amplitudes [Pa] ordered rr, theta-theta, phi-phi, r-theta, r-phi, theta-phi.
    strain : numpy.ndarray of complex128, shape (6,)
        Strain amplitudes in the same order.

    Returns
    -------
    float
        ``|sum_k w_k Im(stress_k conj(strain_k))|`` [Pa] with ``w_k = 2`` on the three off-diagonal components (Europa
        book Eq. 42). For the summed amplitudes of every mode at one forcing frequency, ``|frequency| / 2`` times this
        is that frequency's cycle-averaged volumetric heating [W m-3].

    Raises
    ------
    ValueError
        If either array does not hold exactly 6 components.
    """
    if stress.shape[0] != 6 or strain.shape[0] != 6:
        raise ValueError(
            f"stress and strain must each hold 6 components; got {stress.shape[0]} and {strain.shape[0]}.")
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
        potential6,
        double colatitude):
    """Complex strain and stress amplitudes, and their heating form, at one point for one tidal mode.

    Parameters
    ----------
    y : numpy.ndarray of complex128
        Radial functions y1 to y6 [SI] at ``radius`` from a radial-solver Love solve at the mode's degree and
        ``|frequency|``. Only y1 to y4 are used; a shorter array is zero-padded.
    shear, bulk : complex
        Complex shear and bulk moduli [Pa] at ``radius`` and ``|frequency|``.
    radius : float
        Radius [m].
    degree_l : float
        Harmonic degree of the mode.
    is_solid, is_incompressible : bool
        Assumptions of the layer containing ``radius``; they select dy1/dr. A liquid point returns NaN.
    potential6 : array-like of complex or float
        One mode's potential row ``(U, dU/dtheta, dU/dphi, d2U/dtheta2, d2U/dphi2, d2U/dtheta_dphi)`` as returned
        by ``tidal_potential_3d_modes`` [m2 s-2, per radian for each angular derivative]. Pass the conjugated row
        for a mode with a negative frequency so the amplitudes are at ``+|frequency|``.
    colatitude : float
        Colatitude [rad].

    Returns
    -------
    strain : numpy.ndarray of complex128, shape (6,)
        Strain amplitudes ordered rr, theta-theta, phi-phi, r-theta, r-phi, theta-phi.
    stress : numpy.ndarray of complex128, shape (6,)
        Stress amplitudes [Pa] in the same order.
    heating : float
        ``|sum_k w_k Im(stress_k conj(strain_k))|`` [Pa], as :func:`volumetric_heating`. ``|frequency| / 2`` times
        this is the cycle-averaged volumetric heating [W m-3] of this mode alone; modes sharing a frequency must be
        summed first.

    Raises
    ------
    ValueError
        If ``potential6`` does not hold exactly 6 values.

    Assumptions
    -----------
    - Isotropic linear viscoelasticity, ``sigma = 2 mu eps + lambda tr(eps) delta``, with the Kervazo et al. (2021)
      correction to the Tobie et al. (2005) theta-phi and phi-phi strain forms. Solid layers only.
    - The potential's r^2 factor is taken at the surface radius; the depth dependence is carried by y.
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
    cdef double[12] potential12
    cy_potential_row_to_flat(potential6, &potential12[0])
    cdef double[12] strain12
    cdef double[12] stress12
    cdef double heating = 0.0
    c_strain_stress_heating(
        &y_ri[0],
        shear.real,
        shear.imag,
        bulk.real,
        bulk.imag,
        radius,
        degree_l,
        1 if is_solid else 0,
        1 if is_incompressible else 0,
        &potential12[0],
        colatitude,
        &strain12[0],
        &stress12[0],
        &heating)

    cdef cnp.ndarray[cnp.complex128_t, ndim=1] strain = np.empty(6, dtype=np.complex128)
    cdef cnp.ndarray[cnp.complex128_t, ndim=1] stress = np.empty(6, dtype=np.complex128)
    for k in range(6):
        strain[k] = complex(strain12[2 * k], strain12[2 * k + 1])
        stress[k] = complex(stress12[2 * k], stress12[2 * k + 1])
    return strain, stress, heating


def displacement_point(
        double complex[::1] y not None,
        potential6,
        double colatitude):
    """Complex tidal displacement amplitudes (radial, polar, azimuthal) [m] at one point for one tidal mode.

    Parameters
    ----------
    y : numpy.ndarray of complex128
        Radial functions [SI] at the point's radius from a radial-solver Love solve at the mode's degree and
        ``|frequency|``. Only y1 and y3 are used; a shorter array is zero-padded.
    potential6 : array-like of complex or float
        One mode's potential row as for :func:`strain_stress_heating_point`, conjugated for a negative frequency.
    colatitude : float
        Colatitude [rad].

    Returns
    -------
    numpy.ndarray of complex128, shape (3,)
        ``u_r = y1 U``, ``u_theta = y3 dU/dtheta``, ``u_phi = y3 dU/dphi / sin(colatitude)`` [m]. The displacement at
        time ``t`` is ``Re[u e^{i |frequency| t}]``. At a pole ``u_phi`` is 0 when ``dU/dphi`` is 0 and NaN otherwise.

    Raises
    ------
    ValueError
        If ``potential6`` does not hold exactly 6 values.

    Assumptions
    -----------
    - The radial functions and the potential are in SI (y1 and y3 in s2 m-1, U in m2 s-2).
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
    cdef double[12] potential12
    cy_potential_row_to_flat(potential6, &potential12[0])
    cdef double[6] disp6
    c_displacements_flat(&y_ri[0], &potential12[0], colatitude, &disp6[0])
    cdef cnp.ndarray[cnp.complex128_t, ndim=1] displacement = np.empty(3, dtype=np.complex128)
    for k in range(3):
        displacement[k] = complex(disp6[2 * k], disp6[2 * k + 1])
    return displacement
