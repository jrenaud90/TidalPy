# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython and Python wrappers for TidalPy's cooling models. All quantities MKS.

References
----------
- Turcotte and Schubert (2002), Geodynamics: Rayleigh and Nusselt convection scaling.
- Solomatov (1995); Schubert, Turcotte, and Olson (2001): boundary-layer theory.
"""

from libcpp cimport bool as cpp_bool
from libcpp.memory cimport unique_ptr
from libcpp.utility cimport move
from libcpp.vector cimport vector

cimport numpy as cnp

import numpy as np

cnp.import_array()

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_TidalPyBaseClass
from TidalPy.Utilities_x.classes_x.classes import check_config_keys, factory_defaults

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


cdef void cy_fill_vector(const double[::1] src, vector[double]& dst) noexcept nogil:
    cdef Py_ssize_t n = src.shape[0]
    cdef Py_ssize_t i
    dst.resize(n)
    for i in range(n):
        dst[i] = src[i]


cdef c_CoolingInputs cy_build_inputs(
        double delta_temp,
        double thickness,
        double gravity,
        double density,
        double viscosity,
        double thermal_conductivity,
        double thermal_diffusivity,
        double thermal_expansion) noexcept nogil:
    cdef c_CoolingInputs inp
    inp.delta_temp = delta_temp
    inp.thickness = thickness
    inp.gravity   = gravity
    inp.density   = density
    inp.viscosity = viscosity
    inp.thermal_conductivity = thermal_conductivity
    inp.thermal_diffusivity = thermal_diffusivity
    inp.thermal_expansion   = thermal_expansion
    return inp


cdef CoolingResult cy_result_to_py(c_CoolingResult res):
    """Wrap a scalar c_CoolingResult in a Python CoolingResult."""
    return CoolingResult(res.cooling_flux, res.blt, res.rayleigh_number, res.nusselt_number)


cdef CoolingResult cy_results_to_py(vector[c_CoolingResult]& src, tuple shape):
    """Build a CoolingResult of float64 ndarrays from a std::vector of results."""
    cdef Py_ssize_t n = <Py_ssize_t>src.size()
    cdef Py_ssize_t i
    cdef cnp.ndarray flux = np.empty(n, dtype=np.float64)
    cdef cnp.ndarray blt  = np.empty(n, dtype=np.float64)
    cdef cnp.ndarray ray  = np.empty(n, dtype=np.float64)
    cdef cnp.ndarray nu   = np.empty(n, dtype=np.float64)
    cdef double[::1] m_flux = flux
    cdef double[::1] m_blt  = blt
    cdef double[::1] m_ray  = ray
    cdef double[::1] m_nu   = nu
    cdef c_CoolingResult* cooling_result_ptr = NULL

    with nogil:
        for i in range(n):
            cooling_result_ptr = &src[i]
            m_flux[i] = cooling_result_ptr.cooling_flux
            m_blt[i]  = cooling_result_ptr.blt
            m_ray[i]  = cooling_result_ptr.rayleigh_number
            m_nu[i]   = cooling_result_ptr.nusselt_number
    return CoolingResult(flux.reshape(shape), blt.reshape(shape), ray.reshape(shape), nu.reshape(shape))


cdef object cy_solve_cooling(c_CoolingBase* model, c_CoolingInputs base,
                           object delta_temp, object viscosity):
    """Dispatch to the most specific C++ vectorized routine for the given input pattern."""
    cdef cpp_bool d_arr = isinstance(delta_temp, np.ndarray)
    cdef cpp_bool v_arr = isinstance(viscosity, np.ndarray)

    cdef vector[double] vtemp, vvisc
    cdef vector[c_CoolingResult] vout
    cdef const double[::1] mv
    cdef cnp.ndarray temp_arr, visc_arr, d_b, v_b, d_c, v_c
    cdef tuple out_shape

    if not (d_arr or v_arr):
        base.delta_temp  = <double>delta_temp
        base.viscosity = <double>viscosity
        return cy_result_to_py(model.calc_cooling(base))

    if d_arr and not v_arr:
        base.viscosity = <double>viscosity
        temp_arr = np.ascontiguousarray(delta_temp, dtype=np.float64)
        out_shape = np.shape(temp_arr)
        mv = temp_arr.ravel()
        with nogil:
            cy_fill_vector(mv, vtemp)
            model.calc_cooling_vectorize_temperature(vtemp, base, vout)
        return cy_results_to_py(vout, out_shape)

    if v_arr and not d_arr:
        base.delta_temp = <double>delta_temp
        visc_arr = np.ascontiguousarray(viscosity, dtype=np.float64)
        out_shape = np.shape(visc_arr)
        mv = visc_arr.ravel()
        with nogil:
            cy_fill_vector(mv, vvisc)
            model.calc_cooling_vectorize_viscosity(vvisc, base, vout)
        return cy_results_to_py(vout, out_shape)

    d_b, v_b = np.broadcast_arrays(
        np.asarray(delta_temp, dtype=np.float64),
        np.asarray(viscosity, dtype=np.float64))
    d_c = np.ascontiguousarray(d_b)
    v_c = np.ascontiguousarray(v_b)
    out_shape = np.shape(d_c)
    mv = d_c.ravel()
    with nogil:
        cy_fill_vector(mv, vtemp)
    mv = v_c.ravel()
    with nogil:
        cy_fill_vector(mv, vvisc)
        model.calc_cooling_vectorize_all(vtemp, vvisc, base, vout)
    return cy_results_to_py(vout, out_shape)


cdef class CoolingResult:
    """Container for a cooling model's outputs.

    Each field is a float for scalar evaluations, a float64 ndarray for vectorized ones.

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


cdef class CoolingBase(PhysicsBase):
    """Abstract base for cooling models; holds the owning pointer to the C++ model object."""

    def __cinit__(self, *args, **kwargs):
        pass  # unique_ptr auto-inits to nullptr; concrete models set it

    def __init__(self, *args, **kwargs):
        raise TypeError(
            "CoolingBase is abstract; instantiate a concrete model "
            "(OffCooling, ConvectiveCooling, ConductiveCooling)."
        )

    def __dealloc__(self):
        # unique_ptr frees the most-derived C++ object; _ptr is only an observer.
        self._cooling_ptr.reset()
        self._ptr = NULL

    def calc_cooling(
            self,
            double delta_temp,
            double thickness,
            double gravity,
            double density,
            double viscosity,
            double thermal_conductivity,
            double thermal_diffusivity,
            double thermal_expansion) -> CoolingResult:
        """Cooling result for the given layer state (all MKS).

        Parameters
        ----------
        delta_temp : float
            Temperature drop across the layer [K].
        thickness : float
            Layer (or sub-layer) thickness [m].
        gravity : float
            Gravitational acceleration [m/s^2].
        density : float
            Bulk density [kg/m^3].
        viscosity : float
            Dynamic viscosity [Pa·s].
        thermal_conductivity : float
            Thermal conductivity [W/m/K].
        thermal_diffusivity : float
            Thermal diffusivity [m^2/s].
        thermal_expansion : float
            Thermal expansivity [1/K].

        Returns
        -------
        CoolingResult
        """
        self._check_ptr()
        cdef c_CoolingInputs inp = cy_build_inputs(
            delta_temp,
            thickness,
            gravity,
            density,
            viscosity,
            thermal_conductivity,
            thermal_diffusivity,
            thermal_expansion)
        return cy_result_to_py(self._cooling_ptr.get().calc_cooling(inp))

    def calc_cooling_vectorize_temperature(
            self,
            delta_temp,
            double thickness,
            double gravity,
            double density,
            double viscosity,
            double thermal_conductivity,
            double thermal_diffusivity,
            double thermal_expansion) -> CoolingResult:
        """Cooling over a temperature-drop sweep; the remaining inputs are scalar constants."""
        self._check_ptr()
        cdef c_CoolingInputs base = cy_build_inputs(
            0.0,
            thickness,
            gravity,
            density,
            viscosity,
            thermal_conductivity,
            thermal_diffusivity,
            thermal_expansion)
        cdef vector[double] vtemp
        cdef vector[c_CoolingResult] vout
        cdef const double[::1] mv
        cdef cnp.ndarray temp_c = np.ascontiguousarray(delta_temp, dtype=np.float64).ravel()
        mv = temp_c
        with nogil:
            cy_fill_vector(mv, vtemp)
            self._cooling_ptr.get().calc_cooling_vectorize_temperature(vtemp, base, vout)
        return cy_results_to_py(vout, (temp_c.shape[0],))

    def calc_cooling_vectorize_viscosity(
            self,
            double delta_temp,
            double thickness,
            double gravity,
            double density,
            viscosity,
            double thermal_conductivity,
            double thermal_diffusivity,
            double thermal_expansion) -> CoolingResult:
        """Cooling over a viscosity sweep; the remaining inputs are scalar constants."""
        self._check_ptr()
        cdef c_CoolingInputs base = cy_build_inputs(
            delta_temp,
            thickness,
            gravity,
            density,
            0.0,
            thermal_conductivity,
            thermal_diffusivity,
            thermal_expansion)
        cdef vector[double] vvisc
        cdef vector[c_CoolingResult] vout
        cdef const double[::1] mv
        cdef cnp.ndarray visc_c = np.ascontiguousarray(viscosity, dtype=np.float64).ravel()
        mv = visc_c
        with nogil:
            cy_fill_vector(mv, vvisc)
            self._cooling_ptr.get().calc_cooling_vectorize_viscosity(vvisc, base, vout)
        return cy_results_to_py(vout, (visc_c.shape[0],))

    def calc_cooling_vectorize_all(
            self,
            delta_temp,
            double thickness,
            double gravity,
            double density,
            viscosity,
            double thermal_conductivity,
            double thermal_diffusivity,
            double thermal_expansion) -> CoolingResult:
        """Cooling over element-wise (delta_temp, viscosity) pairs of equal length."""
        self._check_ptr()
        cdef c_CoolingInputs base = cy_build_inputs(
            0.0,
            thickness,
            gravity,
            density,
            0.0,
            thermal_conductivity,
            thermal_diffusivity,
            thermal_expansion)
        cdef vector[double] vtemp, vvisc
        cdef vector[c_CoolingResult] vout
        cdef const double[::1] temp_c_view
        cdef const double[::1] visc_c_view
        cdef cnp.ndarray temp_c = np.ascontiguousarray(delta_temp, dtype=np.float64).ravel()
        cdef cnp.ndarray visc_c = np.ascontiguousarray(viscosity, dtype=np.float64).ravel()
        temp_c_view = temp_c
        visc_c_view = visc_c
        with nogil:
            cy_fill_vector(temp_c_view, vtemp)
            cy_fill_vector(visc_c_view, vvisc)
            self._cooling_ptr.get().calc_cooling_vectorize_all(vtemp, vvisc, base, vout)
        return cy_results_to_py(vout, (temp_c.shape[0],))



cdef class OffCooling(CoolingBase):
    """Cooling disabled: zero heat flux; the boundary layer is half the layer thickness."""

    def __init__(self):
        cdef c_CoolingConfig config
        cdef unique_ptr[c_CoolingBase] ptr = c_find_cooling(c_CoolingModel.Off, config)
        self._cooling_ptr = move(ptr)
        self._ptr = <c_TidalPyBaseClass*>self._cooling_ptr.get()


cdef class ConductiveCooling(CoolingBase):
    """Conduction across the layer: flux = conductivity * delta_temp / thickness."""

    def __init__(self):
        cdef c_CoolingConfig config
        cdef unique_ptr[c_CoolingBase] ptr = c_find_cooling(c_CoolingModel.Conduction, config)
        self._cooling_ptr = move(ptr)
        self._ptr = <c_TidalPyBaseClass*>self._cooling_ptr.get()


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
        cdef unique_ptr[c_CoolingBase] ptr = c_find_cooling(c_CoolingModel.Convection, config)
        self._convective_ptr = <c_ConvectiveCooling*>ptr.get()
        self._cooling_ptr = move(ptr)
        self._ptr = <c_TidalPyBaseClass*>self._cooling_ptr.get()

    def __dealloc__(self):
        self._convective_ptr = NULL  # CoolingBase._cooling_ptr owns the object

    @property
    def convection_alpha(self) -> float:
        """Nusselt scaling prefactor [dimensionless]."""
        self._check_ptr()
        return self._convective_ptr.get_convection_alpha()

    @property
    def convection_beta(self) -> float:
        """Convection exponent [dimensionless]."""
        self._check_ptr()
        return self._convective_ptr.get_convection_beta()

    @property
    def critical_rayleigh(self) -> float:
        """Critical Rayleigh number [dimensionless]."""
        self._check_ptr()
        return self._convective_ptr.get_critical_rayleigh()


# Every config key any cooling model reads; make_cooling rejects anything else.
COOLING_CONFIG_KEYS = frozenset({"convection_alpha", "convection_beta", "critical_rayleigh"})


def _same_model(str table_name, str model_name) -> bool:
    """Whether two names (aliases included) resolve to the same model."""
    return c_cooling_model_from_name(table_name.lower().encode("utf-8")) == c_cooling_model_from_name(model_name.lower().encode("utf-8"))


def make_cooling(str model_name, dict config=None):
    """Build a cooling model from a (case-insensitive) name and config dict.

    Parameters
    ----------
    model_name : str
        Model name or alias: ``off`` (``none``), ``convection`` (``convective``), ``conduction``
        (``conductive``).
    config : dict, optional
        Convection parameters: ``convection_alpha``, ``convection_beta``, ``critical_rayleigh``.
        Ignored by the off and conduction models.

    Returns
    -------
    CoolingBase
        A concrete cooling model instance.

    Raises
    ------
    ValueError
        If the model name is not recognized, or if ``config`` holds a key that no cooling model reads.
    """
    if config is None:
        # Fall back to the same defaults the world-attached path uses.
        config = factory_defaults("cooling", COOLING_CONFIG_KEYS, model_name, _same_model)
    check_config_keys(config, COOLING_CONFIG_KEYS, "cooling")
    if config is None:
        config = {}

    # The default-constructed config carries the C++ defaults, so only override what the caller gave.
    cdef c_CoolingConfig cfg
    if "convection_alpha" in config:
        cfg.convection_alpha = config["convection_alpha"]
    if "convection_beta" in config:
        cfg.convection_beta = config["convection_beta"]
    if "critical_rayleigh" in config:
        cfg.critical_rayleigh = config["critical_rayleigh"]

    cdef c_CoolingModel model = c_cooling_model_from_name(model_name.encode("utf-8"))
    cdef unique_ptr[c_CoolingBase] ptr = c_find_cooling(model, cfg)

    # Adopt the owning unique_ptr into the matching Python wrapper.
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


# Convenience functions. Each builds a stack-allocated C++ model that dies with the call. ``delta_temp``
# (and, for convection, ``viscosity``) accept floats or ndarrays broadcast together; the rest are scalars.

def cooling_off(delta_temp, double thickness):
    """Cooling result for the Off model (zero flux)."""
    cdef c_CoolingConfig cfg
    cdef c_OffCooling model = c_OffCooling(cfg)
    cdef c_CoolingInputs base = cy_build_inputs(0.0, thickness, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    return cy_solve_cooling(<c_CoolingBase*>&model, base, delta_temp, 0.0)


def conductive(
        delta_temp,
        double thickness,
        double thermal_conductivity):
    """Cooling result for the Conduction model: flux = conductivity * delta_temp / thickness."""
    cdef c_CoolingConfig cfg
    cdef c_ConductiveCooling model = c_ConductiveCooling(cfg)
    cdef c_CoolingInputs base = cy_build_inputs(
        0.0,
        thickness,
        0.0,
        0.0,
        0.0,
        thermal_conductivity,
        0.0,
        0.0)
    return cy_solve_cooling(<c_CoolingBase*>&model, base, delta_temp, 0.0)


def convective(
        delta_temp,
        double thickness,
        double gravity,
        double density,
        viscosity,
        double thermal_conductivity,
        double thermal_diffusivity,
        double thermal_expansion,
        double convection_alpha=1.0,
        double convection_beta=0.3333333333333333,
        double critical_rayleigh=1100.0):
    """Cooling result for the parameterized Convection model. Argument order matches ``calc_cooling``."""
    cdef c_CoolingConfig cfg
    cfg.convection_alpha  = convection_alpha
    cfg.convection_beta   = convection_beta
    cfg.critical_rayleigh = critical_rayleigh
    cdef c_ConvectiveCooling model = c_ConvectiveCooling(cfg)
    cdef c_CoolingInputs base = cy_build_inputs(
        0.0,
        thickness,
        gravity,
        density,
        0.0,
        thermal_conductivity,
        thermal_diffusivity,
        thermal_expansion)
    return cy_solve_cooling(<c_CoolingBase*>&model, base, delta_temp, viscosity)
