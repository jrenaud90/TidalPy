# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython and Python wrappers for TidalPy's cooling (heat-transport) models. All quantities MKS.

Each model's parameters, defaults, bounds, and descriptions come from its C++ parameter table, so the wrappers here
only name the models and expose the family's flux law. How a model shapes its layer's profile in a thermal EOS solve
(``build_profile``) lives in C++ and runs inside the solve.

References
----------
- Turcotte and Schubert (2002), Geodynamics: Rayleigh and Nusselt convection scaling.
- Solomatov (1995); Schubert, Turcotte, and Olson (2001): boundary-layer theory.
- Solomatov (2000); Lebrun et al. (2013): the scaling of a liquid (magma ocean) interior.
"""

from libcpp cimport bool as cpp_bool
from libcpp.memory cimport unique_ptr
from libcpp.utility cimport move
from libcpp.vector cimport vector

cimport numpy as cnp

# The result arrays are built through the NumPy C API.
cnp.import_array()

from TidalPy.Utilities.logging.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities.arrays.vectors cimport cy_broadcast_inputs
from TidalPy.Utilities.classes.classes cimport (
    PhysicsBase,
    c_ParamMap,
    c_share_physics,
    cy_collect_parameters,
    cy_param_map,
)
from TidalPy.Utilities.classes.families import ModelFamily

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


cdef c_CoolingInputs cy_build_inputs(
        double delta_temp,
        double thickness,
        double gravity,
        double density,
        double viscosity,
        double thermal_conductivity,
        double thermal_diffusivity,
        double thermal_expansion,
        cpp_bool liquid) noexcept nogil:
    cdef c_CoolingInputs inp
    inp.delta_temp           = delta_temp
    inp.thickness            = thickness
    inp.gravity              = gravity
    inp.density              = density
    inp.viscosity            = viscosity
    inp.thermal_conductivity = thermal_conductivity
    inp.thermal_diffusivity  = thermal_diffusivity
    inp.thermal_expansion    = thermal_expansion
    inp.liquid               = liquid
    return inp


cdef CoolingResult cy_result_to_py(c_CoolingResult res):
    """Wrap a scalar c_CoolingResult in a Python CoolingResult."""
    return CoolingResult(res.cooling_flux, res.blt, res.rayleigh_number, res.nusselt_number)


cdef object cy_solve_cooling(
        c_CoolingBase* model, c_CoolingInputs base, object delta_temp, object viscosity, cpp_bool flatten):
    """Cooling for float or ndarray ``delta_temp`` and ``viscosity`` broadcast together (cy_broadcast_inputs) at the
    otherwise fixed state ``base``; a CoolingResult of floats when both are floats and ``flatten`` is off. A sweep
    writes straight into the four result arrays."""
    cdef vector[vector[double]] inputs
    cdef object shape = cy_broadcast_inputs((delta_temp, viscosity), inputs, flatten)
    if shape is None:
        base.delta_temp = <double>delta_temp
        base.viscosity  = <double>viscosity
        return cy_result_to_py(model.calc_cooling(base))
    # Each input holds one value or one per point.
    cdef cnp.npy_intp num_points = <cnp.npy_intp>max(inputs[0].size(), inputs[1].size())
    cdef cnp.ndarray flux     = cnp.PyArray_EMPTY(1, &num_points, cnp.NPY_FLOAT64, 0)
    cdef cnp.ndarray blt      = cnp.PyArray_EMPTY(1, &num_points, cnp.NPY_FLOAT64, 0)
    cdef cnp.ndarray rayleigh = cnp.PyArray_EMPTY(1, &num_points, cnp.NPY_FLOAT64, 0)
    cdef cnp.ndarray nusselt  = cnp.PyArray_EMPTY(1, &num_points, cnp.NPY_FLOAT64, 0)
    cdef double* flux_ptr     = <double*>cnp.PyArray_DATA(flux)
    cdef double* blt_ptr      = <double*>cnp.PyArray_DATA(blt)
    cdef double* rayleigh_ptr = <double*>cnp.PyArray_DATA(rayleigh)
    cdef double* nusselt_ptr  = <double*>cnp.PyArray_DATA(nusselt)
    with nogil:
        model.calc_cooling_vectorize(
            inputs[0], inputs[1], base, <size_t>num_points, flux_ptr, blt_ptr, rayleigh_ptr, nusselt_ptr)
    return CoolingResult(flux.reshape(shape), blt.reshape(shape), rayleigh.reshape(shape), nusselt.reshape(shape))


cdef class CoolingResult:
    """Container for a cooling model's outputs.

    Each field is a float for scalar evaluations, a float64 ndarray for vectorized ones.

    Attributes
    ----------
    cooling_flux : float or numpy.ndarray
        Heat flux leaving the layer [W/m^2].
    boundary_layer_thickness : float or numpy.ndarray
        Thermal boundary-layer thickness [m]: the conducting thickness that carries the flux across the whole
        temperature drop (thickness / Nu for convection). A layer with a boundary layer at its base and its top
        splits the drop between them, so each is half this thick.
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
    """Base for cooling models, which say how heat moves through a layer: its flux law (``calc_cooling``) and, in a
    thermal EOS solve, the profile of its layer.

    Instantiate a concrete model (``OffCooling``, ``ConductiveCooling``, ``ConvectiveCooling``) with its parameters
    positionally (in the order ``get_parameter_info()`` lists them) or as keywords (argument names or config keys), or
    build one by name with ``make_cooling``. A layer shares the model it is given (``Layer.cooling``).
    """

    # The canonical name of the model a concrete subclass builds; None on this abstract base.
    MODEL_NAME = None

    def __init__(self, *args, dict config=None, **parameters):
        cdef object model_name = type(self).MODEL_NAME
        if model_name is None:
            raise TypeError(
                "CoolingBase is abstract; instantiate a concrete model (OffCooling, ConductiveCooling, "
                "ConvectiveCooling) or call make_cooling.")
        cdef c_ParamMap param_map = cy_param_map(cy_collect_parameters(type(self), args, config, parameters))
        cdef unique_ptr[c_CoolingBase] model = c_find_cooling((<str>model_name).encode("utf-8"), param_map)
        self._set_model(c_share_physics[c_CoolingBase](move(model)))

    cdef c_CoolingBase* _cooling(self) except NULL:
        self._check_ptr()
        return <c_CoolingBase*>self._ptr

    def calc_cooling(
            self,
            double delta_temp,
            double thickness,
            double gravity,
            double density,
            double viscosity,
            double thermal_conductivity,
            double thermal_diffusivity,
            double thermal_expansion,
            cpp_bool liquid=False) -> CoolingResult:
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
        liquid : bool, optional
            The interior is liquid (a magma ocean), so the convection model takes its liquid scaling. Default
            ``False``.

        Returns
        -------
        CoolingResult
        """
        cdef c_CoolingInputs inp = cy_build_inputs(
            delta_temp,
            thickness,
            gravity,
            density,
            viscosity,
            thermal_conductivity,
            thermal_diffusivity,
            thermal_expansion,
            liquid)
        return cy_result_to_py(self._cooling().calc_cooling(inp))

    def calc_cooling_vectorize_temperature(
            self,
            delta_temp,
            double thickness,
            double gravity,
            double density,
            double viscosity,
            double thermal_conductivity,
            double thermal_diffusivity,
            double thermal_expansion,
            cpp_bool liquid=False) -> CoolingResult:
        """Cooling over a temperature-drop sweep; the remaining inputs are scalar constants."""
        cdef c_CoolingInputs base = cy_build_inputs(
            0.0,
            thickness,
            gravity,
            density,
            0.0,
            thermal_conductivity,
            thermal_diffusivity,
            thermal_expansion,
            liquid)
        return cy_solve_cooling(self._cooling(), base, delta_temp, viscosity, True)

    def calc_cooling_vectorize_viscosity(
            self,
            double delta_temp,
            double thickness,
            double gravity,
            double density,
            viscosity,
            double thermal_conductivity,
            double thermal_diffusivity,
            double thermal_expansion,
            cpp_bool liquid=False) -> CoolingResult:
        """Cooling over a viscosity sweep; the remaining inputs are scalar constants."""
        cdef c_CoolingInputs base = cy_build_inputs(
            0.0,
            thickness,
            gravity,
            density,
            0.0,
            thermal_conductivity,
            thermal_diffusivity,
            thermal_expansion,
            liquid)
        return cy_solve_cooling(self._cooling(), base, delta_temp, viscosity, True)

    def calc_cooling_vectorize_all(
            self,
            delta_temp,
            double thickness,
            double gravity,
            double density,
            viscosity,
            double thermal_conductivity,
            double thermal_diffusivity,
            double thermal_expansion,
            cpp_bool liquid=False) -> CoolingResult:
        """Cooling over element-wise (delta_temp, viscosity) pairs of equal length."""
        cdef c_CoolingInputs base = cy_build_inputs(
            0.0,
            thickness,
            gravity,
            density,
            0.0,
            thermal_conductivity,
            thermal_diffusivity,
            thermal_expansion,
            liquid)
        return cy_solve_cooling(self._cooling(), base, delta_temp, viscosity, True)


cdef class OffCooling(CoolingBase):
    """No heat transport inside the layer (alias ``"none"``): it holds one temperature. Its flux law gives zero flux
    and a boundary layer of half the layer thickness."""
    MODEL_NAME = "off"


cdef class ConductiveCooling(CoolingBase):
    """Conduction through the layer (alias ``"conductive"``): flux = conductivity * delta_temp / thickness. In a
    thermal solve the layer is two conducting halves meeting at its mid-radius, where its temperature applies."""
    MODEL_NAME = "conduction"


cdef class ConvectiveCooling(CoolingBase):
    """Parameterized boundary-layer convection (alias ``"convective"``): Nu = alpha (Ra / Ra_crit)^beta for a solid
    interior and Nu = alpha_liquid Ra^beta_liquid for a liquid one (a magma ocean; Solomatov 2000). In a thermal
    solve the layer is a conducting boundary layer at its base and its top around an adiabatic interior whose top is
    at the layer's temperature; the Rayleigh number takes the material's viscosity there, at the top of the
    interior, and the interior is liquid where the Love solve takes it as liquid (a liquid layer, or a layer that can
    melt, molten past the radial solver's solid threshold ``[numerical] minimum_solid_rigidity``)."""
    MODEL_NAME = "convection"


def _canonical_name(str model_name) -> str:
    return c_cooling_canonical_name(model_name.encode("utf-8")).decode("utf-8")


_FAMILY = ModelFamily("cooling", (OffCooling, ConvectiveCooling, ConductiveCooling), _canonical_name)

# Every config key any cooling model reads.
COOLING_CONFIG_KEYS = _FAMILY.config_keys


def _same_model(str table_name, str model_name) -> bool:
    """Whether two names (aliases included) resolve to the same model."""
    return _FAMILY.same_model(table_name, model_name)


def cooling_model_names() -> tuple:
    """The canonical names of the cooling models."""
    return _FAMILY.model_names()


def cooling_config_keys(str model_name) -> frozenset:
    """The config keys a cooling model reads, by any of its names.

    Raises
    ------
    ValueError
        Unknown model name.
    """
    return _FAMILY.config_keys_of(model_name)


def make_cooling(str model_name, dict config=None):
    """Build a cooling model from a (case-insensitive) name and config dict.

    Parameters
    ----------
    model_name : str
        Model name or alias: ``off`` (``none``), ``convection`` (``convective``), ``conduction`` (``conductive``).
    config : dict, optional
        Model parameters by config key (see each model's ``get_parameter_info()``); absent keys (all of them for
        ``None``) take the model's defaults.

    Returns
    -------
    CoolingBase

    Raises
    ------
    ValueError
        Unknown model name, or a parameter the model does not read; each message names the closest accepted one.
    """
    return _FAMILY.make(model_name, config)


# Convenience functions. Each builds the model for the one call. ``delta_temp`` (and, for convection, ``viscosity``)
# accept floats or ndarrays broadcast together; the rest are scalars.

cdef object cy_direct_cooling(
        str model_name,
        dict parameters,
        c_CoolingInputs base,
        object delta_temp,
        object viscosity):
    """Cooling from a model built for the one call (by name and parameters) at the state ``base``, swept over
    ``delta_temp`` and ``viscosity`` (see cy_solve_cooling)."""
    cdef unique_ptr[c_CoolingBase] model = c_find_cooling(model_name.encode("utf-8"), cy_param_map(parameters))
    return cy_solve_cooling(model.get(), base, delta_temp, viscosity, False)


def cooling_off(delta_temp, double thickness):
    """Cooling result for the Off model (zero flux)."""
    return cy_direct_cooling(
        "off", {}, cy_build_inputs(0.0, thickness, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, False), delta_temp, 0.0)


def conductive(
        delta_temp,
        double thickness,
        double thermal_conductivity):
    """Cooling result for the Conduction model: flux = conductivity * delta_temp / thickness."""
    return cy_direct_cooling(
        "conduction",
        {},
        cy_build_inputs(0.0, thickness, 0.0, 0.0, 0.0, thermal_conductivity, 0.0, 0.0, False),
        delta_temp,
        0.0)


def convective(
        delta_temp,
        double thickness,
        double gravity,
        double density,
        viscosity,
        double thermal_conductivity,
        double thermal_diffusivity,
        double thermal_expansion,
        convection_alpha=None,
        convection_beta=None,
        critical_rayleigh=None,
        cpp_bool liquid=False,
        liquid_convection_alpha=None,
        liquid_convection_beta=None):
    """Cooling result for the parameterized Convection model. Argument order matches ``calc_cooling``, then the model
    parameters (``None`` takes the model's default, see ``ConvectiveCooling().get_parameter_info()``); ``liquid`` takes
    the liquid (magma ocean) scaling."""
    cdef dict parameters = {
        key: value for key, value in (
            ("convection_alpha", convection_alpha),
            ("convection_beta", convection_beta),
            ("critical_rayleigh", critical_rayleigh),
            ("liquid_convection_alpha", liquid_convection_alpha),
            ("liquid_convection_beta", liquid_convection_beta))
        if value is not None}
    cdef c_CoolingInputs base = cy_build_inputs(
        0.0,
        thickness,
        gravity,
        density,
        0.0,
        thermal_conductivity,
        thermal_diffusivity,
        thermal_expansion,
        liquid)
    return cy_direct_cooling("convection", parameters, base, delta_temp, viscosity)
