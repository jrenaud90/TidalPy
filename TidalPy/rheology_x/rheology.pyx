# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
rheology.pyx
Cython/Python wrappers for TidalPy's rheology model hierarchy.

Exposes the seven complex-compliance rheology models:

- ``Elastic``  (alias ``"off"``)              — purely elastic, no dissipation.
- ``Viscous``  (alias ``"newton"``)           — purely viscous (Newtonian fluid).
- ``Voigt``    (alias ``"voigt-kelvin"``)     — Voigt-Kelvin element.
- ``Maxwell``                                 — standard Maxwell body.
- ``Burgers``                                 — Maxwell + Voigt in series.
- ``Andrade``                                 — Maxwell + Andrade transient term.
- ``Sundberg`` (alias ``"sundberg-cooper"``)  — Andrade + Voigt.

Each model computes the complex (shear/bulk) modulus mu* [Pa] directly via
``calc_complex_modulus``. Simple models are evaluated analytically; the series
composites (Burgers, Andrade, Sundberg) invert the sum of their element
compliances internally.

References
----------
- Henning, O'Connell, and Sasselov (2009), ApJ, DOI: 10.1088/0004-637X/707/2/1000
- Efroimsky (2012), ApJ, DOI: 10.1088/0004-637X/746/2/150
- Renaud and Henning (2018), ApJ, DOI: 10.3847/1538-4357/aab784
"""

from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex
from libcpp.memory cimport unique_ptr
from libcpp.utility cimport move
from libcpp.vector cimport vector

import numpy as np

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_TidalPyBaseClass

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


# =====================================================================================================================
# Internal helpers for vectorized solving
# =====================================================================================================================
cdef void _fill_vector(double[::1] src, vector[double]& dst) noexcept:
    """Copy a contiguous 1-D float64 memoryview into a std::vector[double]."""
    cdef Py_ssize_t n = src.shape[0]
    cdef Py_ssize_t i
    dst.resize(n)
    for i in range(n):
        dst[i] = src[i]


cdef object _complex_vector_to_ndarray(vector[cpp_complex[double]]& src, tuple shape):
    """Build a complex128 ndarray (of the given shape) from a std::vector."""
    cdef Py_ssize_t n = <Py_ssize_t>src.size()
    cdef Py_ssize_t i
    out = np.empty(n, dtype=np.complex128)
    cdef double complex[::1] mv = out
    for i in range(n):
        mv[i] = src[i].real() + 1j * src[i].imag()
    return out.reshape(shape)


cdef object _solve_complex_modulus(
        c_RheologyBase* model, object frequency, object modulus, object viscosity):
    """Solve the complex modulus for float and/or np.ndarray inputs.

    Picks the most specific C++ vectorized routine for the input pattern:

    - all scalar                       -> scalar calc_complex_modulus (returns complex)
    - frequency array, modulus/visc scalar -> vectorize_frequency
    - modulus/viscosity array(s), frequency scalar -> vectorize_modulus
    - any other mix of arrays          -> broadcast all three -> vectorize_all

    Returns a Python ``complex`` for all-scalar input, else an ``np.ndarray``
    (complex128) broadcast to the common shape.
    """
    cdef cpp_bool f_arr = isinstance(frequency, np.ndarray)
    cdef cpp_bool m_arr = isinstance(modulus, np.ndarray)
    cdef cpp_bool v_arr = isinstance(viscosity, np.ndarray)

    cdef cpp_complex[double] scalar_result
    cdef vector[double] vmod, vvisc, vfreq
    cdef vector[cpp_complex[double]] vout
    cdef double[::1] mv

    # All scalar -> scalar result.
    if not (f_arr or m_arr or v_arr):
        scalar_result = model.calc_complex_modulus(
            <double>modulus, <double>viscosity, <double>frequency)
        return complex(scalar_result.real(), scalar_result.imag())

    # Frequency varies; modulus and viscosity constant.
    if f_arr and not m_arr and not v_arr:
        freq_arr = np.ascontiguousarray(frequency, dtype=np.float64)
        mv = freq_arr.ravel()
        _fill_vector(mv, vfreq)
        model.calc_complex_modulus_vectorize_frequency(
            <double>modulus, <double>viscosity, vfreq, vout)
        return _complex_vector_to_ndarray(vout, freq_arr.shape)

    # Modulus and/or viscosity vary; frequency constant.
    if (m_arr or v_arr) and not f_arr:
        mod_b, visc_b = np.broadcast_arrays(
            np.asarray(modulus, dtype=np.float64),
            np.asarray(viscosity, dtype=np.float64))
        mod_c  = np.ascontiguousarray(mod_b)
        visc_c = np.ascontiguousarray(visc_b)
        mv = mod_c.ravel();  _fill_vector(mv, vmod)
        mv = visc_c.ravel(); _fill_vector(mv, vvisc)
        model.calc_complex_modulus_vectorize_modulus(
            vmod, vvisc, <double>frequency, vout)
        return _complex_vector_to_ndarray(vout, mod_c.shape)

    # General case: broadcast all three and vary everything.
    f_b, m_b, v_b = np.broadcast_arrays(
        np.asarray(frequency, dtype=np.float64),
        np.asarray(modulus, dtype=np.float64),
        np.asarray(viscosity, dtype=np.float64))
    f_c = np.ascontiguousarray(f_b)
    m_c = np.ascontiguousarray(m_b)
    v_c = np.ascontiguousarray(v_b)
    mv = f_c.ravel(); _fill_vector(mv, vfreq)
    mv = m_c.ravel(); _fill_vector(mv, vmod)
    mv = v_c.ravel(); _fill_vector(mv, vvisc)
    model.calc_complex_modulus_vectorize_all(vmod, vvisc, vfreq, vout)
    return _complex_vector_to_ndarray(vout, f_c.shape)


# =====================================================================================================================
# RheologyBase
# =====================================================================================================================
cdef class RheologyBase(PhysicsBase):
    """Abstract base for all rheology models.

    Holds the owning ``unique_ptr`` to the most-derived C++ model object and
    exposes the shared complex-compliance / complex-modulus calculations. Not
    directly instantiable — construct a concrete model (e.g. ``Maxwell``).
    """

    def __cinit__(self, *args, **kwargs):
        pass  # unique_ptr<c_RheologyBase> auto-inits to nullptr; concrete models set it

    def __init__(self, *args, **kwargs):
        raise TypeError(
            "RheologyBase is abstract; instantiate a concrete model "
            "(Elastic, Viscous, Voigt, Maxwell, Burgers, Andrade, Sundberg)."
        )

    def __dealloc__(self):
        # unique_ptr frees the most-derived C++ object; null the base observer ptr.
        self._rheology_ptr.reset()
        self._ptr = NULL

    # ------------------------------------------------------------------------------------------------------------------
    # Calculations
    # ------------------------------------------------------------------------------------------------------------------
    def calc_complex_modulus(self, double modulus_pa,
                             double viscosity_pas,
                             double frequency_rad_s) -> complex:
        """Complex (shear/bulk) modulus mu* [Pa] at the given forcing frequency.

        Parameters
        ----------
        modulus_pa : float
            Unrelaxed (static) modulus [Pa].
        viscosity_pas : float
            Reference dynamic viscosity [Pa·s].
        frequency_rad_s : float
            Tidal forcing frequency [rad/s].

        Returns
        -------
        complex
            Complex modulus [Pa]; real = storage (in-phase), imag = loss
            (out-of-phase, positive = energy loss).

        Assumptions
        -----------
        - Linear viscoelastic regime, single forcing frequency.
        - Andrade-family models assume a positive forcing frequency.
        """
        self._check_ptr()
        cdef cpp_complex[double] result = self._rheology_ptr.get().calc_complex_modulus(
            modulus_pa, viscosity_pas, frequency_rad_s)
        return complex(result.real(), result.imag())

    def calc_complex_modulus_vectorize_modulus(self, modulus_pa, viscosity_pas,
                                               double frequency_rad_s):
        """Complex modulus over (modulus, viscosity) pairs at one frequency.

        Parameters
        ----------
        modulus_pa, viscosity_pas : array_like
            Equal-length 1-D sequences of unrelaxed modulus [Pa] and reference
            viscosity [Pa·s].
        frequency_rad_s : float
            Constant forcing frequency [rad/s].

        Returns
        -------
        numpy.ndarray
            Complex modulus [Pa] for each pair (complex128, same length).
        """
        self._check_ptr()
        cdef vector[double] vmod, vvisc
        cdef vector[cpp_complex[double]] vout
        cdef double[::1] mv
        mod_c  = np.ascontiguousarray(modulus_pa,    dtype=np.float64).ravel()
        visc_c = np.ascontiguousarray(viscosity_pas, dtype=np.float64).ravel()
        mv = mod_c;  _fill_vector(mv, vmod)
        mv = visc_c; _fill_vector(mv, vvisc)
        self._rheology_ptr.get().calc_complex_modulus_vectorize_modulus(
            vmod, vvisc, frequency_rad_s, vout)
        return _complex_vector_to_ndarray(vout, mod_c.shape)

    def calc_complex_modulus_vectorize_frequency(self, double modulus_pa,
                                                 double viscosity_pas, frequency_rad_s):
        """Complex modulus over a frequency sweep at constant modulus/viscosity.

        Parameters
        ----------
        modulus_pa : float
            Constant unrelaxed modulus [Pa].
        viscosity_pas : float
            Constant reference viscosity [Pa·s].
        frequency_rad_s : array_like
            1-D sequence of forcing frequencies [rad/s].

        Returns
        -------
        numpy.ndarray
            Complex modulus [Pa] for each frequency (complex128, same length).
        """
        self._check_ptr()
        cdef vector[double] vfreq
        cdef vector[cpp_complex[double]] vout
        cdef double[::1] mv
        freq_c = np.ascontiguousarray(frequency_rad_s, dtype=np.float64).ravel()
        mv = freq_c; _fill_vector(mv, vfreq)
        self._rheology_ptr.get().calc_complex_modulus_vectorize_frequency(
            modulus_pa, viscosity_pas, vfreq, vout)
        return _complex_vector_to_ndarray(vout, freq_c.shape)

    def calc_complex_modulus_vectorize_all(self, modulus_pa, viscosity_pas,
                                           frequency_rad_s):
        """Complex modulus over element-wise (modulus, viscosity, frequency).

        Parameters
        ----------
        modulus_pa, viscosity_pas, frequency_rad_s : array_like
            Equal-length 1-D sequences (MKS units).

        Returns
        -------
        numpy.ndarray
            Complex modulus [Pa] for each triple (complex128, same length).
        """
        self._check_ptr()
        cdef vector[double] vmod, vvisc, vfreq
        cdef vector[cpp_complex[double]] vout
        cdef double[::1] mv
        mod_c  = np.ascontiguousarray(modulus_pa,      dtype=np.float64).ravel()
        visc_c = np.ascontiguousarray(viscosity_pas,   dtype=np.float64).ravel()
        freq_c = np.ascontiguousarray(frequency_rad_s, dtype=np.float64).ravel()
        mv = mod_c;  _fill_vector(mv, vmod)
        mv = visc_c; _fill_vector(mv, vvisc)
        mv = freq_c; _fill_vector(mv, vfreq)
        self._rheology_ptr.get().calc_complex_modulus_vectorize_all(
            vmod, vvisc, vfreq, vout)
        return _complex_vector_to_ndarray(vout, mod_c.shape)

    # ------------------------------------------------------------------------------------------------------------------
    # Config
    # ------------------------------------------------------------------------------------------------------------------
    cpdef dict get_config_dict(self):
        """Return configuration dict with the model name."""
        return PhysicsBase.get_config_dict(self)


# =====================================================================================================================
# Elastic (alias "off")
# =====================================================================================================================
cdef class Elastic(RheologyBase):
    """Purely elastic rheology — modulus is real, no dissipation.

    Complex modulus: ``mu* = modulus + 0j`` (frequency-independent).
    """

    def __init__(self):
        cdef c_RheologyConfig config
        cdef c_Elastic* raw = new c_Elastic(config)
        self._rheology_ptr.reset(<c_RheologyBase*>raw)
        self._ptr = <c_TidalPyBaseClass*>raw


# =====================================================================================================================
# Viscous (alias "newton")
# =====================================================================================================================
cdef class Viscous(RheologyBase):
    """Purely viscous (Newtonian) rheology.

    Complex modulus: ``mu* = i * viscosity * frequency`` (static modulus unused).
    """

    def __init__(self):
        cdef c_RheologyConfig config
        cdef c_Viscous* raw = new c_Viscous(config)
        self._rheology_ptr.reset(<c_RheologyBase*>raw)
        self._ptr = <c_TidalPyBaseClass*>raw


# =====================================================================================================================
# Maxwell
# =====================================================================================================================
cdef class Maxwell(RheologyBase):
    """Standard Maxwell body.

    Complex modulus: ``mu* = 1 / J``, with Maxwell compliance
    ``J = 1/modulus - i / (viscosity * frequency)``.
    """

    def __init__(self):
        cdef c_RheologyConfig config
        cdef c_Maxwell* raw = new c_Maxwell(config)
        self._rheology_ptr.reset(<c_RheologyBase*>raw)
        self._ptr = <c_TidalPyBaseClass*>raw


# =====================================================================================================================
# Voigt (alias "voigt-kelvin")
# =====================================================================================================================
cdef class Voigt(RheologyBase):
    """Voigt-Kelvin element.

    Parameters
    ----------
    voigt_modulus_frac : float, optional
        Voigt compliance fraction relative to the layer compliance.
        Default ``0.2``.
    voigt_viscosity_frac : float, optional
        Voigt viscosity fraction relative to the layer viscosity.
        Default ``0.02``.
    """

    def __cinit__(self, *args, **kwargs):
        self._voigt_ptr = NULL

    def __init__(self, double voigt_modulus_frac=5.0,
                 double voigt_viscosity_frac=0.02):
        cdef c_RheologyConfig config
        config.voigt_modulus_frac   = voigt_modulus_frac
        config.voigt_viscosity_frac = voigt_viscosity_frac
        cdef c_Voigt* raw = new c_Voigt(config)
        self._rheology_ptr.reset(<c_RheologyBase*>raw)
        self._voigt_ptr = raw
        self._ptr       = <c_TidalPyBaseClass*>raw

    def __dealloc__(self):
        self._voigt_ptr = NULL  # base unique_ptr owns the C++ object

    @property
    def voigt_modulus_frac(self) -> float:
        """Voigt compliance fraction [dimensionless]."""
        return self._voigt_ptr.get_voigt_modulus_frac()

    @property
    def voigt_viscosity_frac(self) -> float:
        """Voigt viscosity fraction [dimensionless]."""
        return self._voigt_ptr.get_voigt_viscosity_frac()

    cpdef dict get_config_dict(self):
        """Return config dict with model name and Voigt parameters."""
        d = RheologyBase.get_config_dict(self)
        d["voigt_modulus_frac"]   = self._voigt_ptr.get_voigt_modulus_frac()
        d["voigt_viscosity_frac"] = self._voigt_ptr.get_voigt_viscosity_frac()
        return d


# =====================================================================================================================
# Burgers
# =====================================================================================================================
cdef class Burgers(RheologyBase):
    """Burgers rheology — Maxwell and Voigt elements in series.

    Parameters
    ----------
    voigt_modulus_frac : float, optional
        Voigt compliance fraction relative to the layer compliance.
        Default ``0.2``.
    voigt_viscosity_frac : float, optional
        Voigt viscosity fraction relative to the layer viscosity.
        Default ``0.02``.
    """

    def __cinit__(self, *args, **kwargs):
        self._burgers_ptr = NULL

    def __init__(self, double voigt_modulus_frac=5.0,
                 double voigt_viscosity_frac=0.02):
        cdef c_RheologyConfig config
        config.voigt_modulus_frac   = voigt_modulus_frac
        config.voigt_viscosity_frac = voigt_viscosity_frac
        cdef c_Burgers* raw = new c_Burgers(config)
        self._rheology_ptr.reset(<c_RheologyBase*>raw)
        self._burgers_ptr = raw
        self._ptr         = <c_TidalPyBaseClass*>raw

    def __dealloc__(self):
        self._burgers_ptr = NULL

    @property
    def voigt_modulus_frac(self) -> float:
        """Voigt compliance fraction [dimensionless]."""
        return self._burgers_ptr.get_voigt_modulus_frac()

    @property
    def voigt_viscosity_frac(self) -> float:
        """Voigt viscosity fraction [dimensionless]."""
        return self._burgers_ptr.get_voigt_viscosity_frac()

    cpdef dict get_config_dict(self):
        """Return config dict with model name and Voigt parameters."""
        d = RheologyBase.get_config_dict(self)
        d["voigt_modulus_frac"]   = self._burgers_ptr.get_voigt_modulus_frac()
        d["voigt_viscosity_frac"] = self._burgers_ptr.get_voigt_viscosity_frac()
        return d


# =====================================================================================================================
# Andrade
# =====================================================================================================================
cdef class Andrade(RheologyBase):
    """Andrade rheology — Maxwell body plus a transient term ~ omega^{-alpha}.

    Parameters
    ----------
    alpha : float, optional
        Andrade exponent [dimensionless]. Default ``0.3``.
    zeta : float, optional
        Andrade timescale ratio [dimensionless]. Default ``1.0``.
    """

    def __cinit__(self, *args, **kwargs):
        self._andrade_ptr = NULL

    def __init__(self, double alpha=0.3, double zeta=1.0):
        cdef c_RheologyConfig config
        config.alpha = alpha
        config.zeta  = zeta
        cdef c_Andrade* raw = new c_Andrade(config)
        self._rheology_ptr.reset(<c_RheologyBase*>raw)
        self._andrade_ptr = raw
        self._ptr         = <c_TidalPyBaseClass*>raw

    def __dealloc__(self):
        self._andrade_ptr = NULL

    @property
    def alpha(self) -> float:
        """Andrade exponent [dimensionless]."""
        return self._andrade_ptr.get_alpha()

    @property
    def zeta(self) -> float:
        """Andrade timescale ratio [dimensionless]."""
        return self._andrade_ptr.get_zeta()

    cpdef dict get_config_dict(self):
        """Return config dict with model name and Andrade parameters."""
        d = RheologyBase.get_config_dict(self)
        d["alpha"] = self._andrade_ptr.get_alpha()
        d["zeta"]  = self._andrade_ptr.get_zeta()
        return d


# =====================================================================================================================
# Sundberg (alias "sundberg-cooper")
# =====================================================================================================================
cdef class Sundberg(RheologyBase):
    """Sundberg-Cooper rheology — Andrade and Voigt elements summed.

    Parameters
    ----------
    alpha : float, optional
        Andrade exponent [dimensionless]. Default ``0.3``.
    zeta : float, optional
        Andrade timescale ratio [dimensionless]. Default ``1.0``.
    voigt_modulus_frac : float, optional
        Voigt compliance fraction relative to the layer compliance.
        Default ``0.2``.
    voigt_viscosity_frac : float, optional
        Voigt viscosity fraction relative to the layer viscosity.
        Default ``0.02``.
    """

    def __cinit__(self, *args, **kwargs):
        self._sundberg_ptr = NULL

    def __init__(self, double alpha=0.3, double zeta=1.0,
                 double voigt_modulus_frac=5.0,
                 double voigt_viscosity_frac=0.02):
        cdef c_RheologyConfig config
        config.alpha                = alpha
        config.zeta                 = zeta
        config.voigt_modulus_frac   = voigt_modulus_frac
        config.voigt_viscosity_frac = voigt_viscosity_frac
        cdef c_Sundberg* raw = new c_Sundberg(config)
        self._rheology_ptr.reset(<c_RheologyBase*>raw)
        self._sundberg_ptr = raw
        self._ptr          = <c_TidalPyBaseClass*>raw

    def __dealloc__(self):
        self._sundberg_ptr = NULL

    @property
    def alpha(self) -> float:
        """Andrade exponent [dimensionless]."""
        return self._sundberg_ptr.get_alpha()

    @property
    def zeta(self) -> float:
        """Andrade timescale ratio [dimensionless]."""
        return self._sundberg_ptr.get_zeta()

    @property
    def voigt_modulus_frac(self) -> float:
        """Voigt compliance fraction [dimensionless]."""
        return self._sundberg_ptr.get_voigt_modulus_frac()

    @property
    def voigt_viscosity_frac(self) -> float:
        """Voigt viscosity fraction [dimensionless]."""
        return self._sundberg_ptr.get_voigt_viscosity_frac()

    cpdef dict get_config_dict(self):
        """Return config dict with model name and all Sundberg parameters."""
        d = RheologyBase.get_config_dict(self)
        d["alpha"]                = self._sundberg_ptr.get_alpha()
        d["zeta"]                 = self._sundberg_ptr.get_zeta()
        d["voigt_modulus_frac"]   = self._sundberg_ptr.get_voigt_modulus_frac()
        d["voigt_viscosity_frac"] = self._sundberg_ptr.get_voigt_viscosity_frac()
        return d


# =====================================================================================================================
# Factory
# =====================================================================================================================

def make_rheology(str model_name, dict config=None):
    """Build a rheology model from a (case-insensitive) name and config dict.

    This wraps the C++ enum factory: the name (and any alias) is mapped to a
    ``c_RheologyModel`` enum by ``c_rheology_model_from_name``, the matching
    concrete model is heap-allocated by ``c_find_rheology`` (returning a
    ``unique_ptr``), and the owning pointer is adopted into the correct rich
    Python wrapper.

    Parameters
    ----------
    model_name : str
        Model name or alias. Recognized names: ``elastic`` (``off``),
        ``viscous`` (``newton``), ``voigt`` (``voigt-kelvin``), ``maxwell``,
        ``burgers``, ``andrade``, ``sundberg`` (``sundberg-cooper``).
    config : dict, optional
        Model parameters. Recognized keys (model-dependent):
        ``alpha``, ``zeta``, ``voigt_modulus_frac``, ``voigt_viscosity_frac``.
        Unused keys are ignored; missing keys fall back to model defaults.

    Returns
    -------
    RheologyBase
        A concrete rheology model instance.

    Raises
    ------
    ValueError
        If the model name is not recognized.
    """
    if config is None:
        config = {}

    # Build the config struct from the dict (missing keys use legacy defaults).
    cdef c_RheologyConfig cfg
    cfg.alpha                = config.get("alpha", 0.3)
    cfg.zeta                 = config.get("zeta", 1.0)
    cfg.voigt_modulus_frac   = config.get("voigt_modulus_frac", 5.0)
    cfg.voigt_viscosity_frac = config.get("voigt_viscosity_frac", 0.02)

    # Map name/alias -> enum (raises ValueError on unknown name via except +),
    # then build the model through the canonical C++ enum factory.
    cdef c_RheologyModel model = c_rheology_model_from_name(model_name.encode("utf-8"))
    cdef unique_ptr[c_RheologyBase] ptr = c_find_rheology(model, cfg)

    # Adopt the owning unique_ptr into the matching rich Python wrapper.
    cdef Elastic  e
    cdef Viscous  vi
    cdef Maxwell  m
    cdef Voigt    vo
    cdef Burgers  b
    cdef Andrade  a
    cdef Sundberg s

    if model == c_RheologyModel.Elastic:
        e = Elastic.__new__(Elastic)
        e._rheology_ptr = move(ptr)
        e._ptr = <c_TidalPyBaseClass*>e._rheology_ptr.get()
        return e
    elif model == c_RheologyModel.Viscous:
        vi = Viscous.__new__(Viscous)
        vi._rheology_ptr = move(ptr)
        vi._ptr = <c_TidalPyBaseClass*>vi._rheology_ptr.get()
        return vi
    elif model == c_RheologyModel.Maxwell:
        m = Maxwell.__new__(Maxwell)
        m._rheology_ptr = move(ptr)
        m._ptr = <c_TidalPyBaseClass*>m._rheology_ptr.get()
        return m
    elif model == c_RheologyModel.Voigt:
        vo = Voigt.__new__(Voigt)
        vo._voigt_ptr    = <c_Voigt*>ptr.get()
        vo._rheology_ptr = move(ptr)
        vo._ptr          = <c_TidalPyBaseClass*>vo._voigt_ptr
        return vo
    elif model == c_RheologyModel.Burgers:
        b = Burgers.__new__(Burgers)
        b._burgers_ptr  = <c_Burgers*>ptr.get()
        b._rheology_ptr = move(ptr)
        b._ptr          = <c_TidalPyBaseClass*>b._burgers_ptr
        return b
    elif model == c_RheologyModel.Andrade:
        a = Andrade.__new__(Andrade)
        a._andrade_ptr  = <c_Andrade*>ptr.get()
        a._rheology_ptr = move(ptr)
        a._ptr          = <c_TidalPyBaseClass*>a._andrade_ptr
        return a
    else:  # c_RheologyModel.Sundberg
        s = Sundberg.__new__(Sundberg)
        s._sundberg_ptr = <c_Sundberg*>ptr.get()
        s._rheology_ptr = move(ptr)
        s._ptr          = <c_TidalPyBaseClass*>s._sundberg_ptr
        return s


# =====================================================================================================================
# Direct complex-modulus convenience functions
#
# Each builds a stack-allocated C++ model from its parameters, solves for the
# complex (shear/bulk) modulus, and returns the result (the C++ model is
# destroyed when the function returns). ``frequency``, ``modulus`` and
# ``viscosity`` accept Python floats or NumPy arrays (broadcast together); model
# parameters such as ``alpha``/``zeta`` are always scalar constants. A float
# result is returned for all-scalar inputs, otherwise a complex128 ``ndarray``.
# =====================================================================================================================

def elastic(frequency, modulus, viscosity):
    """Complex shear/bulk modulus for the Elastic model [Pa]. See module notes."""
    cdef c_RheologyConfig cfg
    cdef c_Elastic model = c_Elastic(cfg)
    return _solve_complex_modulus(<c_RheologyBase*>&model, frequency, modulus, viscosity)


def viscous(frequency, modulus, viscosity):
    """Complex shear/bulk modulus for the Viscous (Newton) model [Pa]."""
    cdef c_RheologyConfig cfg
    cdef c_Viscous model = c_Viscous(cfg)
    return _solve_complex_modulus(<c_RheologyBase*>&model, frequency, modulus, viscosity)


def maxwell(frequency, modulus, viscosity):
    """Complex shear/bulk modulus for the Maxwell model [Pa]."""
    cdef c_RheologyConfig cfg
    cdef c_Maxwell model = c_Maxwell(cfg)
    return _solve_complex_modulus(<c_RheologyBase*>&model, frequency, modulus, viscosity)


def voigt(frequency, modulus, viscosity,
          double voigt_modulus_frac=5.0, double voigt_viscosity_frac=0.02):
    """Complex shear/bulk modulus for the Voigt-Kelvin model [Pa]."""
    cdef c_RheologyConfig cfg
    cfg.voigt_modulus_frac   = voigt_modulus_frac
    cfg.voigt_viscosity_frac = voigt_viscosity_frac
    cdef c_Voigt model = c_Voigt(cfg)
    return _solve_complex_modulus(<c_RheologyBase*>&model, frequency, modulus, viscosity)


def burgers(frequency, modulus, viscosity,
            double voigt_modulus_frac=5.0, double voigt_viscosity_frac=0.02):
    """Complex shear/bulk modulus for the Burgers model [Pa]."""
    cdef c_RheologyConfig cfg
    cfg.voigt_modulus_frac   = voigt_modulus_frac
    cfg.voigt_viscosity_frac = voigt_viscosity_frac
    cdef c_Burgers model = c_Burgers(cfg)
    return _solve_complex_modulus(<c_RheologyBase*>&model, frequency, modulus, viscosity)


def andrade(frequency, modulus, viscosity, double alpha=0.3, double zeta=1.0):
    """Complex shear/bulk modulus for the Andrade model [Pa]."""
    cdef c_RheologyConfig cfg
    cfg.alpha = alpha
    cfg.zeta  = zeta
    cdef c_Andrade model = c_Andrade(cfg)
    return _solve_complex_modulus(<c_RheologyBase*>&model, frequency, modulus, viscosity)


def sundberg(frequency, modulus, viscosity, double alpha=0.3, double zeta=1.0,
             double voigt_modulus_frac=5.0, double voigt_viscosity_frac=0.02):
    """Complex shear/bulk modulus for the Sundberg-Cooper model [Pa]."""
    cdef c_RheologyConfig cfg
    cfg.alpha                = alpha
    cfg.zeta                 = zeta
    cfg.voigt_modulus_frac   = voigt_modulus_frac
    cfg.voigt_viscosity_frac = voigt_viscosity_frac
    cdef c_Sundberg model = c_Sundberg(cfg)
    return _solve_complex_modulus(<c_RheologyBase*>&model, frequency, modulus, viscosity)
