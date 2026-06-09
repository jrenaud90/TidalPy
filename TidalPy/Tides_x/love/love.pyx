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
