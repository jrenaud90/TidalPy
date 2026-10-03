# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython and Python wrappers for TidalPy's material property laws: equations of state and shear-modulus laws.

A material's phase combines one law of each family with viscosity laws and thermal constants. Each model's
parameters, defaults, bounds, and descriptions come from its C++ parameter table; ``get_parameter_info()`` lists
them.
"""

from libcpp cimport bool as cpp_bool
from libcpp.memory cimport unique_ptr
from libcpp.utility cimport move
from libcpp.vector cimport vector

cimport numpy as cnp
cnp.import_array()

from TidalPy.Utilities.logging.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport d_NAN, set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities.arrays.vectors cimport cy_broadcast_inputs, cy_vector_to_ndarray
from TidalPy.Utilities.classes.classes cimport (
    PhysicsBase,
    c_ParamMap,
    c_ThermoPoint,
    c_share_physics,
    cy_collect_parameters,
    cy_param_map,
)
from TidalPy.Utilities.classes.families import ModelFamily

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


# =====================================================================================================================
# Equations of State
# =====================================================================================================================
cdef class EOSBase(PhysicsBase):
    """Base for equation-of-state laws, which give a phase's density [kg m-3], isothermal and adiabatic bulk moduli
    [Pa], and thermal expansivity [1/K] at a pressure, temperature, and radius.

    Instantiate a concrete law with its parameters positionally (in the order ``get_parameter_info()`` lists them) or
    as keywords (argument names or config keys), or build one by name with ``make_eos``.
    """

    # The canonical name of the model a concrete subclass builds; None on this abstract base.
    MODEL_NAME = None

    def __init__(self, *args, dict config=None, **parameters):
        cdef object model_name = type(self).MODEL_NAME
        if model_name is None:
            raise TypeError("EOSBase is abstract; instantiate a concrete law or call make_eos.")
        cdef c_ParamMap param_map = cy_param_map(cy_collect_parameters(type(self), args, config, parameters))
        cdef unique_ptr[c_EOSBase] model = c_find_eos((<str>model_name).encode("utf-8"), param_map)
        self._set_model(c_share_physics[c_EOSBase](move(model)))

    cdef c_EOSBase* _eos(self) except NULL:
        self._check_ptr()
        return <c_EOSBase*>self._ptr

    def calc_eos(self, pressure, temperature=d_NAN, radius=d_NAN, thermal=False) -> dict:
        """The law at a point.

        Parameters
        ----------
        pressure : float or np.ndarray
            Pressure [Pa].
        temperature : float or np.ndarray, optional
            Temperature [K]; NaN (the default) gives the athermal law.
        radius : float or np.ndarray, optional
            Radius [m], read only by ``InterpolatedEOS``.
        thermal : bool, optional
            Whether the density sees the temperature. The expansivity is reported either way.

        Returns
        -------
        dict
            ``density`` [kg m-3], ``bulk_modulus`` (isothermal) [Pa], ``adiabatic_bulk_modulus`` [Pa], and
            ``thermal_expansion`` [1/K]; floats when every input is a float, otherwise arrays of their broadcast shape.
        """
        cdef c_EOSBase* model_ptr = self._eos()
        cdef cpp_bool use_thermal = True if thermal else False
        cdef vector[vector[double]] inputs
        cdef vector[double] density, bulk_modulus, adiabatic_bulk_modulus, thermal_expansion
        cdef c_ThermoPoint point
        cdef c_EOSPoint law_point
        cdef object shape = cy_broadcast_inputs((pressure, temperature, radius), inputs, False)
        if shape is None:
            point.pressure    = <double>pressure
            point.temperature = <double>temperature
            point.radius      = <double>radius
            model_ptr.calc_eos(point, use_thermal, law_point)
            return {
                "density": law_point.density,
                "bulk_modulus": law_point.bulk_modulus,
                "adiabatic_bulk_modulus": law_point.adiabatic_bulk_modulus,
                "thermal_expansion": law_point.thermal_expansion,
            }
        with nogil:
            model_ptr.calc_eos_vectorize(
                inputs[0], inputs[1], inputs[2], use_thermal, density, bulk_modulus, adiabatic_bulk_modulus,
                thermal_expansion)
        return {
            "density": cy_vector_to_ndarray(density, shape),
            "bulk_modulus": cy_vector_to_ndarray(bulk_modulus, shape),
            "adiabatic_bulk_modulus": cy_vector_to_ndarray(adiabatic_bulk_modulus, shape),
            "thermal_expansion": cy_vector_to_ndarray(thermal_expansion, shape),
        }

    def calc_density(self, pressure, temperature=d_NAN, radius=d_NAN, thermal=False):
        """Density [kg m-3]; see ``calc_eos``."""
        return self.calc_eos(pressure, temperature, radius, thermal)["density"]


cdef class ConstantEOS(EOSBase):
    """Uniform density (aliases ``"uniform"``, ``"constant_density"``) with a constant bulk modulus."""
    MODEL_NAME = "constant"


cdef class BirchMurnaghanEOS(EOSBase):
    """3rd-order Birch-Murnaghan law (aliases ``"bm"``, ``"birch-murnaghan"``), inverted for the density, with a
    thermal pressure alpha0 K0 (T - T_ref)."""
    MODEL_NAME = "birch_murnaghan"


cdef class VinetEOS(EOSBase):
    """Vinet (universal) law, inverted for the density, with a thermal pressure alpha0 K0 (T - T_ref)."""
    MODEL_NAME = "vinet"


cdef class MurnaghanEOS(EOSBase):
    """Murnaghan law, K = K0 + K0' P, closed form both ways; suited to melts and liquids over modest pressures."""
    MODEL_NAME = "murnaghan"


cdef class PolytropeEOS(EOSBase):
    """Polytrope P = K rho^(1 + 1/n), for gas envelopes; zero density where the pressure is not positive."""
    MODEL_NAME = "polytrope"


cdef class ModifiedPolytropeEOS(EOSBase):
    """Modified polytrope rho = rho0 + c P^n (Seager et al. 2007; alias ``"seager"``), for cold compression to TPa."""
    MODEL_NAME = "modified_polytrope"


cdef class InterpolatedEOS(EOSBase):
    """A density profile tabulated in radius (aliases ``"interp"``, ``"interpolated"``), with an optional bulk-modulus
    table: ``radius_m``, ``density_kg_m3``, and ``bulk_modulus_pa``."""
    MODEL_NAME = "interpolate"


def _eos_canonical_name(str model_name) -> str:
    return c_eos_canonical_name(model_name.encode("utf-8")).decode("utf-8")


_EOS_FAMILY = ModelFamily(
    "equation of state",
    (ConstantEOS, BirchMurnaghanEOS, VinetEOS, MurnaghanEOS, PolytropeEOS, ModifiedPolytropeEOS, InterpolatedEOS),
    _eos_canonical_name)

# Every config key any equation-of-state law reads.
EOS_CONFIG_KEYS = _EOS_FAMILY.config_keys


def eos_model_names() -> tuple:
    """The canonical names of the equation-of-state laws."""
    return _EOS_FAMILY.model_names()


def canonical_eos_name(str model_name) -> str:
    """An equation-of-state law's canonical name from any of its names or aliases (case-insensitive).

    Raises
    ------
    ValueError
        Unknown model name; the message names the closest one.
    """
    return _EOS_FAMILY.canonical_name(model_name)


def eos_config_keys(str model_name) -> frozenset:
    """The config keys an equation-of-state law reads, by any of its names."""
    return _EOS_FAMILY.config_keys_of(model_name)


def make_eos(str model_name, dict config=None) -> EOSBase:
    """Build an equation-of-state law by name.

    Parameters
    ----------
    model_name : str
        ``"constant"`` (``"uniform"``), ``"birch_murnaghan"`` (``"bm"``), ``"vinet"``, ``"murnaghan"``,
        ``"polytrope"``, ``"modified_polytrope"`` (``"seager"``), or ``"interpolate"`` (``"interp"``).
    config : dict, optional
        Law parameters by config key; absent keys take the law's defaults.

    Raises
    ------
    ValueError
        Unknown model name, or a parameter the law does not read; each message names the closest accepted one.
    """
    return _EOS_FAMILY.make(model_name, config)


# =====================================================================================================================
# Shear Modulus Laws
# =====================================================================================================================
cdef class ShearModulusBase(PhysicsBase):
    """Base for shear-modulus laws, which give a phase's static (unrelaxed) shear modulus [Pa] at a pressure,
    temperature, and radius.

    Instantiate a concrete law with its parameters positionally (in the order ``get_parameter_info()`` lists them) or
    as keywords (argument names or config keys), or build one by name with ``make_shear_modulus``.
    """

    # The canonical name of the model a concrete subclass builds; None on this abstract base.
    MODEL_NAME = None

    def __init__(self, *args, dict config=None, **parameters):
        cdef object model_name = type(self).MODEL_NAME
        if model_name is None:
            raise TypeError("ShearModulusBase is abstract; instantiate a concrete law or call make_shear_modulus.")
        cdef c_ParamMap param_map = cy_param_map(cy_collect_parameters(type(self), args, config, parameters))
        cdef unique_ptr[c_ShearModulusBase] model = c_find_shear_modulus(
            (<str>model_name).encode("utf-8"), param_map)
        self._set_model(c_share_physics[c_ShearModulusBase](move(model)))

    cdef c_ShearModulusBase* _shear_modulus(self) except NULL:
        self._check_ptr()
        return <c_ShearModulusBase*>self._ptr

    def calc_shear_modulus(self, pressure, temperature=d_NAN, radius=d_NAN):
        """Static shear modulus [Pa]; a float when every input is a float, otherwise an array of their broadcast
        shape. A NaN temperature (the default) leaves out any temperature dependence."""
        cdef c_ShearModulusBase* model_ptr = self._shear_modulus()
        cdef vector[vector[double]] inputs
        cdef vector[double] shear_modulus
        cdef c_ThermoPoint point
        cdef object shape = cy_broadcast_inputs((pressure, temperature, radius), inputs, False)
        if shape is None:
            point.pressure    = <double>pressure
            point.temperature = <double>temperature
            point.radius      = <double>radius
            return model_ptr.calc_shear_modulus(point)
        with nogil:
            model_ptr.calc_shear_modulus_vectorize(inputs[0], inputs[1], inputs[2], shear_modulus)
        return cy_vector_to_ndarray(shear_modulus, shape)


cdef class ConstantShearModulus(ShearModulusBase):
    """A constant shear modulus (alias ``"const"``)."""
    MODEL_NAME = "constant"


cdef class LinearShearModulus(ShearModulusBase):
    """mu = mu0 + mu'_P (P - P_ref) + mu'_T (T - T_ref)."""
    MODEL_NAME = "linear"


cdef class InterpolatedShearModulus(ShearModulusBase):
    """A shear-modulus profile tabulated in radius (aliases ``"interp"``, ``"interpolated"``): ``radius_m`` and
    ``shear_modulus_pa`` tables."""
    MODEL_NAME = "interpolate"


def _shear_canonical_name(str model_name) -> str:
    return c_shear_modulus_canonical_name(model_name.encode("utf-8")).decode("utf-8")


_SHEAR_FAMILY = ModelFamily(
    "shear modulus",
    (ConstantShearModulus, LinearShearModulus, InterpolatedShearModulus),
    _shear_canonical_name)

# Every config key any shear-modulus law reads.
SHEAR_MODULUS_CONFIG_KEYS = _SHEAR_FAMILY.config_keys


def shear_modulus_model_names() -> tuple:
    """The canonical names of the shear-modulus laws."""
    return _SHEAR_FAMILY.model_names()


def canonical_shear_modulus_name(str model_name) -> str:
    """A shear-modulus law's canonical name from any of its names or aliases (case-insensitive).

    Raises
    ------
    ValueError
        Unknown model name; the message names the closest one.
    """
    return _SHEAR_FAMILY.canonical_name(model_name)


def shear_modulus_config_keys(str model_name) -> frozenset:
    """The config keys a shear-modulus law reads, by any of its names."""
    return _SHEAR_FAMILY.config_keys_of(model_name)


def make_shear_modulus(str model_name, dict config=None) -> ShearModulusBase:
    """Build a shear-modulus law by name: ``"constant"`` (``"const"``), ``"linear"``, or ``"interpolate"``
    (``"interp"``).

    Raises
    ------
    ValueError
        Unknown model name, or a parameter the law does not read; each message names the closest accepted one.
    """
    return _SHEAR_FAMILY.make(model_name, config)
