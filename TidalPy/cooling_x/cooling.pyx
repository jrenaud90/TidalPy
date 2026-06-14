# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
cooling.pyx
Cython/Python wrappers for TidalPy's cooling model hierarchy.

Exposes the three heat-transport models:

- ``OffCooling``         (alias ``"none"``)  — cooling disabled (zero flux).
- ``ConvectiveCooling``  (alias ``"convective"``) — parameterized boundary-layer convection.
- ``ConductiveCooling``  (alias ``"conductive"``) — conduction across the layer.

Each model maps the layer's physical state (temperature drop, geometry, and
material/transport properties) to a ``CoolingResult`` carrying the surface heat
flux [W/m^2], boundary-layer thickness [m], and the Rayleigh and Nusselt numbers.
All quantities are MKS.

References
----------
- Turcotte and Schubert (2002), Geodynamics — Rayleigh/Nusselt convection scaling.
- Solomatov (1995); Schubert, Turcotte, and Olson (2001) — boundary-layer theory.
"""

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
# Internal helpers
# =====================================================================================================================
cdef void _fill_vector(double[::1] src, vector[double]& dst) noexcept:
    """Copy a contiguous 1-D float64 memoryview into a std::vector[double]."""
    cdef Py_ssize_t n = src.shape[0]
    cdef Py_ssize_t i
    dst.resize(n)
    for i in range(n):
        dst[i] = src[i]


cdef c_CoolingInputs _build_inputs(
        double delta_temp_k, double thickness_m, double gravity_m_s2, double density_kg_m3,
        double viscosity_pas, double thermal_conductivity_w_mk, double thermal_diffusivity_m2_s,
        double thermal_expansion_1_k):
    """Pack the eight cooling inputs into a c_CoolingInputs struct."""
    cdef c_CoolingInputs inp
    inp.delta_temp_k              = delta_temp_k
    inp.thickness_m               = thickness_m
    inp.gravity_m_s2              = gravity_m_s2
    inp.density_kg_m3             = density_kg_m3
    inp.viscosity_pas             = viscosity_pas
    inp.thermal_conductivity_w_mk = thermal_conductivity_w_mk
    inp.thermal_diffusivity_m2_s  = thermal_diffusivity_m2_s
    inp.thermal_expansion_1_k     = thermal_expansion_1_k
    return inp


cdef CoolingResult _result_to_py(c_CoolingResult res):
    """Wrap a scalar c_CoolingResult in a Python CoolingResult."""
    return CoolingResult(res.cooling_flux_w_m2, res.blt_m, res.rayleigh_number, res.nusselt_number)


cdef CoolingResult _results_to_py(vector[c_CoolingResult]& src, tuple shape):
    """Build a CoolingResult of float64 ndarrays from a std::vector of results."""
    cdef Py_ssize_t n = <Py_ssize_t>src.size()
    cdef Py_ssize_t i
    flux = np.empty(n, dtype=np.float64)
    blt  = np.empty(n, dtype=np.float64)
    ray  = np.empty(n, dtype=np.float64)
    nu   = np.empty(n, dtype=np.float64)
    cdef double[::1] m_flux = flux
    cdef double[::1] m_blt  = blt
    cdef double[::1] m_ray  = ray
    cdef double[::1] m_nu   = nu
    for i in range(n):
        m_flux[i] = src[i].cooling_flux_w_m2
        m_blt[i]  = src[i].blt_m
        m_ray[i]  = src[i].rayleigh_number
        m_nu[i]   = src[i].nusselt_number
    return CoolingResult(flux.reshape(shape), blt.reshape(shape), ray.reshape(shape), nu.reshape(shape))


cdef object _solve_cooling(c_CoolingBase* model, c_CoolingInputs base,
                           object delta_temp, object viscosity):
    """Solve cooling for float and/or np.ndarray ``delta_temp`` and ``viscosity``.

    Picks the most specific C++ vectorized routine for the input pattern:

    - both scalar               -> scalar calc_cooling
    - delta_temp array, visc scalar -> vectorize_temperature
    - visc array, delta_temp scalar -> vectorize_viscosity
    - any other mix of arrays   -> broadcast both -> vectorize_all

    Returns a ``CoolingResult`` whose fields are floats for all-scalar input, else
    float64 ndarrays broadcast to the common shape.
    """
    cdef bint d_arr = isinstance(delta_temp, np.ndarray)
    cdef bint v_arr = isinstance(viscosity, np.ndarray)

    cdef vector[double] vtemp, vvisc
    cdef vector[c_CoolingResult] vout
    cdef double[::1] mv

    # All scalar -> scalar result.
    if not (d_arr or v_arr):
        base.delta_temp_k  = <double>delta_temp
        base.viscosity_pas = <double>viscosity
        return _result_to_py(model.calc_cooling(base))

    # Temperature varies; viscosity constant.
    if d_arr and not v_arr:
        base.viscosity_pas = <double>viscosity
        temp_arr = np.ascontiguousarray(delta_temp, dtype=np.float64)
        mv = temp_arr.ravel(); _fill_vector(mv, vtemp)
        model.calc_cooling_vectorize_temperature(vtemp, base, vout)
        return _results_to_py(vout, temp_arr.shape)

    # Viscosity varies; temperature constant.
    if v_arr and not d_arr:
        base.delta_temp_k = <double>delta_temp
        visc_arr = np.ascontiguousarray(viscosity, dtype=np.float64)
        mv = visc_arr.ravel(); _fill_vector(mv, vvisc)
        model.calc_cooling_vectorize_viscosity(vvisc, base, vout)
        return _results_to_py(vout, visc_arr.shape)

    # General case: broadcast both and vary everything.
    d_b, v_b = np.broadcast_arrays(
        np.asarray(delta_temp, dtype=np.float64),
        np.asarray(viscosity, dtype=np.float64))
    d_c = np.ascontiguousarray(d_b)
    v_c = np.ascontiguousarray(v_b)
    mv = d_c.ravel(); _fill_vector(mv, vtemp)
    mv = v_c.ravel(); _fill_vector(mv, vvisc)
    model.calc_cooling_vectorize_all(vtemp, vvisc, base, vout)
    return _results_to_py(vout, d_c.shape)


# =====================================================================================================================
# CoolingResult
# =====================================================================================================================
cdef class CoolingResult:
    """Container for a cooling model's outputs.

    Each field is a Python ``float`` for scalar evaluations or a ``float64``
    ``numpy.ndarray`` for vectorized evaluations.

    Attributes
    ----------
    cooling_flux : float or numpy.ndarray
        Heat flux leaving the layer [W/m^2].
    boundary_layer_thickness : float or numpy.ndarray
        Thermal boundary-layer thickness [m].
    rayleigh : float or numpy.ndarray
        Rayleigh number [dimensionless] (0 for off/conduction).
    nusselt : float or numpy.ndarray
        Nusselt number [dimensionless] (1 for off/conduction).
    """

    def __init__(self, cooling_flux, boundary_layer_thickness, rayleigh, nusselt):
        self.cooling_flux = cooling_flux
        self.boundary_layer_thickness = boundary_layer_thickness
        self.rayleigh = rayleigh
        self.nusselt = nusselt

    def to_dict(self) -> dict:
        """Return the four outputs as a dict."""
        return {
            "cooling_flux": self.cooling_flux,
            "boundary_layer_thickness": self.boundary_layer_thickness,
            "rayleigh": self.rayleigh,
            "nusselt": self.nusselt,
        }

    def __iter__(self):
        yield self.cooling_flux
        yield self.boundary_layer_thickness
        yield self.rayleigh
        yield self.nusselt

    def __repr__(self):
        return (f"CoolingResult(cooling_flux={self.cooling_flux!r}, "
                f"boundary_layer_thickness={self.boundary_layer_thickness!r}, "
                f"rayleigh={self.rayleigh!r}, nusselt={self.nusselt!r})")


# =====================================================================================================================
# CoolingBase
# =====================================================================================================================
cdef class CoolingBase(PhysicsBase):
    """Abstract base for all cooling models.

    Holds the owning ``unique_ptr`` to the most-derived C++ model object and
    exposes the shared cooling calculations. Not directly instantiable —
    construct a concrete model (e.g. ``ConvectiveCooling``).
    """

    def __cinit__(self, *args, **kwargs):
        pass  # unique_ptr<c_CoolingBase> auto-inits to nullptr; concrete models set it

    def __init__(self, *args, **kwargs):
        raise TypeError(
            "CoolingBase is abstract; instantiate a concrete model "
            "(OffCooling, ConvectiveCooling, ConductiveCooling)."
        )

    def __dealloc__(self):
        # unique_ptr frees the most-derived C++ object; null the base observer ptr.
        self._cooling_ptr.reset()
        self._ptr = NULL

    # ------------------------------------------------------------------------------------------------------------------
    # Calculations
    # ------------------------------------------------------------------------------------------------------------------
    def calc_cooling(self, double delta_temp_k, double thickness_m, double gravity_m_s2,
                     double density_kg_m3, double viscosity_pas,
                     double thermal_conductivity_w_mk, double thermal_diffusivity_m2_s,
                     double thermal_expansion_1_k) -> CoolingResult:
        """Cooling result for the given layer state (all MKS).

        Parameters
        ----------
        delta_temp_k : float
            Temperature drop across the layer [K].
        thickness_m : float
            Layer (or sub-layer) thickness [m].
        gravity_m_s2 : float
            Gravitational acceleration [m/s^2].
        density_kg_m3 : float
            Bulk density [kg/m^3].
        viscosity_pas : float
            Dynamic viscosity [Pa·s].
        thermal_conductivity_w_mk : float
            Thermal conductivity [W/m/K].
        thermal_diffusivity_m2_s : float
            Thermal diffusivity [m^2/s].
        thermal_expansion_1_k : float
            Thermal expansivity [1/K].

        Returns
        -------
        CoolingResult
        """
        self._check_ptr()
        cdef c_CoolingInputs inp = _build_inputs(
            delta_temp_k, thickness_m, gravity_m_s2, density_kg_m3, viscosity_pas,
            thermal_conductivity_w_mk, thermal_diffusivity_m2_s, thermal_expansion_1_k)
        return _result_to_py(self._cooling_ptr.get().calc_cooling(inp))

    def calc_cooling_vectorize_temperature(self, delta_temp_k, double thickness_m,
                                           double gravity_m_s2, double density_kg_m3,
                                           double viscosity_pas, double thermal_conductivity_w_mk,
                                           double thermal_diffusivity_m2_s,
                                           double thermal_expansion_1_k) -> CoolingResult:
        """Cooling over a temperature-drop sweep at otherwise fixed state.

        ``delta_temp_k`` is an array; the remaining inputs are scalar constants.
        Returns a ``CoolingResult`` of float64 ndarrays.
        """
        self._check_ptr()
        cdef c_CoolingInputs base = _build_inputs(
            0.0, thickness_m, gravity_m_s2, density_kg_m3, viscosity_pas,
            thermal_conductivity_w_mk, thermal_diffusivity_m2_s, thermal_expansion_1_k)
        cdef vector[double] vtemp
        cdef vector[c_CoolingResult] vout
        cdef double[::1] mv
        temp_c = np.ascontiguousarray(delta_temp_k, dtype=np.float64).ravel()
        mv = temp_c; _fill_vector(mv, vtemp)
        self._cooling_ptr.get().calc_cooling_vectorize_temperature(vtemp, base, vout)
        return _results_to_py(vout, temp_c.shape)

    def calc_cooling_vectorize_viscosity(self, double delta_temp_k, double thickness_m,
                                         double gravity_m_s2, double density_kg_m3,
                                         viscosity_pas, double thermal_conductivity_w_mk,
                                         double thermal_diffusivity_m2_s,
                                         double thermal_expansion_1_k) -> CoolingResult:
        """Cooling over a viscosity sweep at otherwise fixed state.

        ``viscosity_pas`` is an array; the remaining inputs are scalar constants.
        Returns a ``CoolingResult`` of float64 ndarrays.
        """
        self._check_ptr()
        cdef c_CoolingInputs base = _build_inputs(
            delta_temp_k, thickness_m, gravity_m_s2, density_kg_m3, 0.0,
            thermal_conductivity_w_mk, thermal_diffusivity_m2_s, thermal_expansion_1_k)
        cdef vector[double] vvisc
        cdef vector[c_CoolingResult] vout
        cdef double[::1] mv
        visc_c = np.ascontiguousarray(viscosity_pas, dtype=np.float64).ravel()
        mv = visc_c; _fill_vector(mv, vvisc)
        self._cooling_ptr.get().calc_cooling_vectorize_viscosity(vvisc, base, vout)
        return _results_to_py(vout, visc_c.shape)

    def calc_cooling_vectorize_all(self, delta_temp_k, double thickness_m, double gravity_m_s2,
                                   double density_kg_m3, viscosity_pas,
                                   double thermal_conductivity_w_mk, double thermal_diffusivity_m2_s,
                                   double thermal_expansion_1_k) -> CoolingResult:
        """Cooling over element-wise (delta_temp, viscosity) pairs.

        ``delta_temp_k`` and ``viscosity_pas`` are equal-length arrays; the
        remaining inputs are scalar constants. Returns a ``CoolingResult`` of
        float64 ndarrays.
        """
        self._check_ptr()
        cdef c_CoolingInputs base = _build_inputs(
            0.0,
            thickness_m,
            gravity_m_s2,
            density_kg_m3,
            0.0,
            thermal_conductivity_w_mk,
            thermal_diffusivity_m2_s,
            thermal_expansion_1_k)
        cdef vector[double] vtemp, vvisc
        cdef vector[c_CoolingResult] vout
        cdef double[::1] mv
        temp_c = np.ascontiguousarray(delta_temp_k, dtype=np.float64).ravel()
        visc_c = np.ascontiguousarray(viscosity_pas, dtype=np.float64).ravel()
        mv = temp_c; _fill_vector(mv, vtemp)
        mv = visc_c; _fill_vector(mv, vvisc)
        self._cooling_ptr.get().calc_cooling_vectorize_all(vtemp, vvisc, base, vout)
        return _results_to_py(vout, temp_c.shape)

    # ------------------------------------------------------------------------------------------------------------------
    # Config
    # ------------------------------------------------------------------------------------------------------------------
    cpdef dict get_config_dict(self):
        """Return configuration dict with the model name."""
        return PhysicsBase.get_config_dict(self)


# =====================================================================================================================
# OffCooling (alias "none")
# =====================================================================================================================
cdef class OffCooling(CoolingBase):
    """Cooling disabled — zero heat flux; boundary layer is half the thickness."""

    def __init__(self):
        cdef c_CoolingConfig config
        cdef c_OffCooling* raw = new c_OffCooling(config)
        self._cooling_ptr.reset(<c_CoolingBase*>raw)
        self._ptr = <c_TidalPyBaseClass*>raw


# =====================================================================================================================
# ConductiveCooling (alias "conductive")
# =====================================================================================================================
cdef class ConductiveCooling(CoolingBase):
    """Conduction across the layer: flux = k · ΔT / thickness."""

    def __init__(self):
        cdef c_CoolingConfig config
        cdef c_ConductiveCooling* raw = new c_ConductiveCooling(config)
        self._cooling_ptr.reset(<c_CoolingBase*>raw)
        self._ptr = <c_TidalPyBaseClass*>raw


# =====================================================================================================================
# ConvectiveCooling (alias "convective")
# =====================================================================================================================
cdef class ConvectiveCooling(CoolingBase):
    """Parameterized boundary-layer convection via the Rayleigh number.

    Parameters
    ----------
    convection_alpha : float, optional
        Nusselt scaling prefactor (``Nu = alpha · (Ra / Ra_crit)^beta``). Default ``1.0``.
    convection_beta : float, optional
        Convection exponent. Default ``1/3``.
    critical_rayleigh : float, optional
        Critical Rayleigh number. Default ``1100.0``.
    """

    def __cinit__(self, *args, **kwargs):
        self._convective_ptr = NULL

    def __init__(self,
            double convection_alpha=1.0,
            double convection_beta=0.3333333333333333,
            double critical_rayleigh=1100.0):
        cdef c_CoolingConfig config
        config.convection_alpha  = convection_alpha
        config.convection_beta   = convection_beta
        config.critical_rayleigh = critical_rayleigh
        cdef c_ConvectiveCooling* raw = new c_ConvectiveCooling(config)
        self._cooling_ptr.reset(<c_CoolingBase*>raw)
        self._convective_ptr = raw
        self._ptr            = <c_TidalPyBaseClass*>raw

    def __dealloc__(self):
        self._convective_ptr = NULL  # base unique_ptr owns the C++ object

    @property
    def convection_alpha(self) -> float:
        """Nusselt scaling prefactor [dimensionless]."""
        return self._convective_ptr.get_convection_alpha()

    @property
    def convection_beta(self) -> float:
        """Convection exponent [dimensionless]."""
        return self._convective_ptr.get_convection_beta()

    @property
    def critical_rayleigh(self) -> float:
        """Critical Rayleigh number [dimensionless]."""
        return self._convective_ptr.get_critical_rayleigh()

    cpdef dict get_config_dict(self):
        """Return config dict with model name and convection parameters."""
        d = CoolingBase.get_config_dict(self)
        d["convection_alpha"]  = self._convective_ptr.get_convection_alpha()
        d["convection_beta"]   = self._convective_ptr.get_convection_beta()
        d["critical_rayleigh"] = self._convective_ptr.get_critical_rayleigh()
        return d


# =====================================================================================================================
# Factory
# =====================================================================================================================
def make_cooling(str model_name, dict config=None):
    """Build a cooling model from a (case-insensitive) name and config dict.

    This wraps the C++ enum factory: the name (and any alias) is mapped to a
    ``c_CoolingModel`` enum by ``c_cooling_model_from_name``, the matching
    concrete model is heap-allocated by ``c_find_cooling`` (returning a
    ``unique_ptr``), and the owning pointer is adopted into the correct rich
    Python wrapper.

    Parameters
    ----------
    model_name : str
        Model name or alias. Recognised names: ``off`` (``none``), ``convection``
        (``convective``), ``conduction`` (``conductive``).
    config : dict, optional
        Convection parameters: ``convection_alpha``, ``convection_beta``,
        ``critical_rayleigh``. Ignored by the off and conduction models.

    Returns
    -------
    CoolingBase
        A concrete cooling model instance.

    Raises
    ------
    ValueError
        If the model name is not recognised.
    """
    if config is None:
        config = {}

    cdef c_CoolingConfig cfg
    cfg.convection_alpha  = config.get("convection_alpha", 1.0)
    cfg.convection_beta   = config.get("convection_beta", 0.3333333333333333)
    cfg.critical_rayleigh = config.get("critical_rayleigh", 1100.0)

    # Map name/alias -> enum (raises ValueError on unknown name via except +),
    # then build the model through the canonical C++ enum factory.
    cdef c_CoolingModel model = c_cooling_model_from_name(model_name.encode("utf-8"))
    cdef unique_ptr[c_CoolingBase] ptr = c_find_cooling(model, cfg)

    # Adopt the owning unique_ptr into the matching rich Python wrapper.
    cdef OffCooling        o
    cdef ConvectiveCooling cv
    cdef ConductiveCooling cd

    if model == c_CoolingModel.Off:
        o = OffCooling.__new__(OffCooling)
        o._cooling_ptr = move(ptr)
        o._ptr = <c_TidalPyBaseClass*>o._cooling_ptr.get()
        return o
    elif model == c_CoolingModel.Convection:
        cv = ConvectiveCooling.__new__(ConvectiveCooling)
        cv._convective_ptr = <c_ConvectiveCooling*>ptr.get()
        cv._cooling_ptr    = move(ptr)
        cv._ptr            = <c_TidalPyBaseClass*>cv._convective_ptr
        return cv
    else:  # c_CoolingModel.Conduction
        cd = ConductiveCooling.__new__(ConductiveCooling)
        cd._cooling_ptr = move(ptr)
        cd._ptr = <c_TidalPyBaseClass*>cd._cooling_ptr.get()
        return cd


# =====================================================================================================================
# Direct cooling convenience functions
#
# Each builds a stack-allocated C++ model from its parameters, solves for the
# cooling result, and returns a CoolingResult (the C++ model is destroyed when the
# function returns). ``delta_temp_k`` (and, for convection, ``viscosity_pas``)
# accept Python floats or NumPy arrays (broadcast together); the remaining physical
# inputs and model parameters are scalar constants. A CoolingResult of floats is
# returned for all-scalar inputs, otherwise of float64 ndarrays.
# =====================================================================================================================

def cooling_off(delta_temp_k, double thickness_m):
    """Cooling result for the Off model (zero flux). See module notes."""
    cdef c_CoolingConfig cfg
    cdef c_OffCooling model = c_OffCooling(cfg)
    cdef c_CoolingInputs base = _build_inputs(0.0, thickness_m, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    return _solve_cooling(<c_CoolingBase*>&model, base, delta_temp_k, 0.0)


def conductive(
        delta_temp_k,
        double thickness_m,
        double thermal_conductivity_w_mk):
    """Cooling result for the Conduction model: flux = k · ΔT / thickness."""
    cdef c_CoolingConfig cfg
    cdef c_ConductiveCooling model = c_ConductiveCooling(cfg)
    cdef c_CoolingInputs base = _build_inputs(
        0.0, thickness_m, 0.0, 0.0, 0.0, thermal_conductivity_w_mk, 0.0, 0.0)
    return _solve_cooling(<c_CoolingBase*>&model, base, delta_temp_k, 0.0)


def convective(
        delta_temp_k,
        viscosity_pas,
        double thickness_m,
        double gravity_m_s2,
        double density_kg_m3,
        double thermal_conductivity_w_mk,
        double thermal_diffusivity_m2_s,
        double thermal_expansion_1_k,
        double convection_alpha=1.0,
        double convection_beta=0.3333333333333333,
        double critical_rayleigh=1100.0):
    """Cooling result for the parameterized Convection model.

    ``delta_temp_k`` and ``viscosity_pas`` may be floats or arrays (broadcast
    together); the remaining inputs are scalar constants.
    """
    cdef c_CoolingConfig cfg
    cfg.convection_alpha  = convection_alpha
    cfg.convection_beta   = convection_beta
    cfg.critical_rayleigh = critical_rayleigh
    cdef c_ConvectiveCooling model = c_ConvectiveCooling(cfg)
    cdef c_CoolingInputs base = _build_inputs(
        0.0, thickness_m, gravity_m_s2, density_kg_m3, 0.0,
        thermal_conductivity_w_mk, thermal_diffusivity_m2_s, thermal_expansion_1_k)
    return _solve_cooling(<c_CoolingBase*>&model, base, delta_temp_k, viscosity_pas)
