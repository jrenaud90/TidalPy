# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
love.pyx
Cython/Python wrapper for TidalPy's Love numbers container.

LoveNumbers: an immutable Python view of the c_LoveNumbers C++ struct,
holding the three complex tidal Love numbers k, h, and l.
"""

from libcpp.complex cimport complex as cpp_complex

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


# =====================================================================================================================
# LoveNumbers
# =====================================================================================================================
cdef class LoveNumbers:
    """Container for the three complex tidal Love numbers k, h, and l.

    The Love numbers describe the elastic response of a planetary body to a
    tidal potential:

    - **k** (potential Love number): describes how the external gravitational
      potential changes due to the body's redistribution of mass.
    - **h** (radial displacement Love number): describes radial (vertical)
      surface deformation.
    - **l** (tangential displacement Love number): describes horizontal
      (tangential) surface deformation.

    All three are dimensionless complex numbers. The real part is the elastic
    amplitude; the imaginary part represents energy dissipation at the tidal
    forcing frequency.

    Parameters
    ----------
    k : complex, optional
        Potential Love number [dimensionless]. Default ``0+0j``.
    h : complex, optional
        Radial displacement Love number [dimensionless]. Default ``0+0j``.
    l : complex, optional
        Tangential displacement Love number [dimensionless]. Default ``0+0j``.
    """

    def __init__(self, complex k=0+0j, complex h=0+0j, complex l=0+0j):
        self._love.k = cpp_complex[double](k.real, k.imag)
        self._love.h = cpp_complex[double](h.real, h.imag)
        self._love.l = cpp_complex[double](l.real, l.imag)

    # ------------------------------------------------------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def k(self) -> complex:
        """Potential Love number [dimensionless]."""
        return complex(self._love.k.real(), self._love.k.imag())

    @property
    def h(self) -> complex:
        """Radial displacement Love number [dimensionless]."""
        return complex(self._love.h.real(), self._love.h.imag())

    @property
    def l(self) -> complex:
        """Tangential displacement Love number [dimensionless]."""
        return complex(self._love.l.real(), self._love.l.imag())

    # ------------------------------------------------------------------------------------------------------------------
    # Python protocol
    # ------------------------------------------------------------------------------------------------------------------
    def __repr__(self) -> str:
        return (f"LoveNumbers("
                f"k={complex(self._love.k.real(), self._love.k.imag())}, "
                f"h={complex(self._love.h.real(), self._love.h.imag())}, "
                f"l={complex(self._love.l.real(), self._love.l.imag())})")

    def __eq__(self, other) -> bool:
        if not isinstance(other, LoveNumbers):
            return NotImplemented
        return (self._love.k == (<LoveNumbers>other)._love.k and
                self._love.h == (<LoveNumbers>other)._love.h and
                self._love.l == (<LoveNumbers>other)._love.l)

    def __iter__(self):
        """Iterate as (k, h, l) — enables tuple unpacking."""
        yield complex(self._love.k.real(), self._love.k.imag())
        yield complex(self._love.h.real(), self._love.h.imag())
        yield complex(self._love.l.real(), self._love.l.imag())

    # ------------------------------------------------------------------------------------------------------------------
    # Serialisation helpers
    # ------------------------------------------------------------------------------------------------------------------
    cpdef dict to_dict(self):
        """Return all components as a flat dict with re/im suffixes.

        Returns
        -------
        dict
            Keys: ``love_number_k_re``, ``love_number_k_im``,
            ``love_number_h_re``, ``love_number_h_im``,
            ``love_number_l_re``, ``love_number_l_im``.
        """
        return {
            "love_number_k_re": self._love.k.real(),
            "love_number_k_im": self._love.k.imag(),
            "love_number_h_re": self._love.h.real(),
            "love_number_h_im": self._love.h.imag(),
            "love_number_l_re": self._love.l.real(),
            "love_number_l_im": self._love.l.imag(),
        }


# =====================================================================================================================
# Love-number methods and the homogeneous-sphere formulas
# =====================================================================================================================

cdef LoveNumbers _wrap_love(c_LoveNumbers love):
    """Wrap a C++ Love-number struct in the Python container."""
    cdef LoveNumbers out = LoveNumbers()
    out._love = love
    return out


def love_method_name(str method) -> str:
    """Canonical name of a Love-number method given any accepted alias.

    Parameters
    ----------
    method : str
        ``'radial_solver'`` (``'shooting'``, ``'rs'``), ``'propagation_matrix'`` (``'prop_matrix'``,
        ``'pm'``, ``'prop'``), ``'homogeneous'`` (``'homogen'``), ``'cpl'``, ``'ctl'``, or
        ``'laterally_inhomogeneous'`` (``'3d'``, ``'lat_inhom'``); case-insensitive.

    Returns
    -------
    str
        The canonical method name. Unknown names raise ``ValueError``.
    """
    return c_love_method_name_int(c_parse_love_method_int(method.encode('utf-8'))).decode('utf-8')


def calc_effective_rigidity(
        shear_modulus,
        double density,
        double gravity,
        double radius,
        int degree_l = 2):
    """Degree-l effective rigidity of a homogeneous sphere, ``(2 l^2 + 4 l + 3) / l * mu / (rho g R)``.

    Parameters
    ----------
    shear_modulus : float or complex
        Shear modulus [Pa]; a complex (viscoelastic) modulus gives a complex effective rigidity.
    density : float
        Bulk density [kg m-3].
    gravity : float
        Surface gravity [m s-2].
    radius : float
        Radius [m].
    degree_l : int, default 2
        Harmonic degree (>= 2). At degree 2 the prefactor is 19/2.

    Returns
    -------
    float or complex
        Effective rigidity [dimensionless].

    Assumptions
    -----------
    - Homogeneous incompressible sphere (Love 1911).
    """
    cdef complex mu_c
    cdef cpp_complex[double] out_c
    if isinstance(shear_modulus, complex):
        mu_c = shear_modulus
        out_c = c_calc_effective_rigidity_complex(
            cpp_complex[double](mu_c.real, mu_c.imag), density, gravity, radius, degree_l)
        return complex(out_c.real(), out_c.imag())
    return c_calc_effective_rigidity_real(<double>shear_modulus, density, gravity, radius, degree_l)


def calc_homogeneous_love_numbers(
        complex complex_shear_modulus,
        double density,
        double gravity,
        double radius,
        int degree_l = 2) -> LoveNumbers:
    """Love numbers k_l, h_l, l_l of a homogeneous incompressible sphere.

    ``k_l = 3 / (2 (l - 1)) / (1 + mu_eff)``, ``h_l = (2 l + 1) / (2 (l - 1)) / (1 + mu_eff)``,
    ``l_l = 3 / (2 l (l - 1)) / (1 + mu_eff)`` with the effective rigidity of :func:`calc_effective_rigidity`.
    A complex shear modulus (from a rheology model at the forcing frequency) gives complex Love numbers.

    Parameters
    ----------
    complex_shear_modulus : complex
        Shear modulus [Pa]; pass a real value for the static (elastic) Love numbers.
    density : float
        Bulk density [kg m-3].
    gravity : float
        Surface gravity [m s-2].
    radius : float
        Radius [m].
    degree_l : int, default 2
        Harmonic degree (>= 2).

    Returns
    -------
    LoveNumbers

    Assumptions
    -----------
    - Homogeneous incompressible sphere; the fluid limit (mu -> 0) gives k_2 = 3/2, h_2 = 5/2, l_2 = 3/4.
    """
    return _wrap_love(c_calc_homogeneous_love_numbers(
        cpp_complex[double](complex_shear_modulus.real, complex_shear_modulus.imag),
        density, gravity, radius, degree_l))


def apply_fixed_q(LoveNumbers love_numbers not None, double fixed_q) -> LoveNumbers:
    """Constant phase lag: multiply k, h, and l by ``(1 - i/Q)`` so that ``-Im[k] = Re[k] / Q``.

    Parameters
    ----------
    love_numbers : LoveNumbers
        Static (real) Love numbers.
    fixed_q : float
        Tidal quality factor (> 0).

    Returns
    -------
    LoveNumbers
    """
    return _wrap_love(c_apply_fixed_q(love_numbers._love, fixed_q))


def apply_fixed_dt(LoveNumbers love_numbers not None, double frequency, double fixed_dt) -> LoveNumbers:
    """Constant time lag: multiply k, h, and l by ``(1 - i |omega| dt)`` so that ``-Im[k] = Re[k] |omega| dt``.

    Parameters
    ----------
    love_numbers : LoveNumbers
        Static (real) Love numbers.
    frequency : float
        Tidal forcing frequency [rad s-1] (its magnitude is used).
    fixed_dt : float
        Time lag [s] (>= 0).

    Returns
    -------
    LoveNumbers
    """
    return _wrap_love(c_apply_fixed_dt(love_numbers._love, frequency, fixed_dt))
