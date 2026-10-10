# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython and Python wrappers for TidalPy's melting laws: melting curves (solidus and liquidus), melt weakening of the
shear modulus and viscosity, and the optional bulk-modulus and bulk-viscosity mixing of melt.

A material with a liquid phase holds two melting curves and one weakening law, and optionally the mixing laws. Each
model's parameters, defaults, bounds, and descriptions come from its C++ parameter table; ``get_parameter_info()``
lists them.
"""

from libcpp.memory cimport unique_ptr
from libcpp.utility cimport move
from libcpp.vector cimport vector

cimport numpy as cnp
cnp.import_array()

from TidalPy.Utilities.logging.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities.arrays.vectors cimport cy_broadcast_inputs, cy_vector_to_ndarray
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


cdef c_ParamMap cy_model_parameters(object model_class, tuple args, dict config, dict parameters,
                                    str family) except *:
    """A concrete model's constructor arguments as a c_ParamMap; raises for the abstract family base."""
    if model_class.MODEL_NAME is None:
        raise TypeError(f"{model_class.__name__} is abstract; instantiate a concrete {family} model.")
    return cy_param_map(cy_collect_parameters(model_class, args, config, parameters))


# =====================================================================================================================
# Melting Curves
# =====================================================================================================================
cdef class MeltingCurveBase(PhysicsBase):
    """Base for melting curves, which give a solidus or liquidus temperature [K] at a pressure [Pa].

    Instantiate a concrete curve with its parameters positionally (in the order ``get_parameter_info()`` lists them)
    or as keywords (argument names or config keys), or build one by name with ``make_melting_curve``.
    """

    # The canonical name of the model a concrete subclass builds; None on this abstract base.
    MODEL_NAME = None

    def __init__(self, *args, dict config=None, **parameters):
        cdef c_ParamMap param_map = cy_model_parameters(type(self), args, config, parameters, "melting curve")
        cdef unique_ptr[c_MeltingCurveBase] model = c_find_melting_curve(
            (<str>type(self).MODEL_NAME).encode("utf-8"), param_map)
        self._set_model(c_share_physics[c_MeltingCurveBase](move(model)))

    cdef c_MeltingCurveBase* _curve(self) except NULL:
        self._check_ptr()
        return <c_MeltingCurveBase*>self._ptr

    def calc_melting_temperature(self, pressure):
        """Melting temperature [K] at a pressure [Pa]; a float for a float, an array of the same shape otherwise."""
        cdef c_MeltingCurveBase* model_ptr = self._curve()
        cdef vector[vector[double]] inputs
        cdef vector[double] temperature
        cdef object shape = cy_broadcast_inputs((pressure,), inputs, False)
        if shape is None:
            return model_ptr.calc_melting_temperature(<double>pressure)
        with nogil:
            model_ptr.calc_melting_temperature_vectorize(inputs[0], temperature)
        return cy_vector_to_ndarray(temperature, shape)

    def calc_melting_slope(self, pressure):
        """Slope of the curve, dT_m/dP [K Pa-1], at a pressure [Pa]; zero where the curve is held flat. A float for a
        float, an array of the same shape otherwise."""
        cdef c_MeltingCurveBase* model_ptr = self._curve()
        cdef vector[vector[double]] inputs
        cdef vector[double] slope
        cdef object shape = cy_broadcast_inputs((pressure,), inputs, False)
        if shape is None:
            return model_ptr.calc_melting_slope(<double>pressure)
        with nogil:
            model_ptr.calc_melting_slope_vectorize(inputs[0], slope)
        return cy_vector_to_ndarray(slope, shape)


cdef class ConstantMeltingCurve(MeltingCurveBase):
    """A melting temperature independent of pressure (alias ``"const"``)."""
    MODEL_NAME = "constant"


cdef class SimonGlatzelCurve(MeltingCurveBase):
    """One Simon and Glatzel (1929) branch, T(P) = T0 (1 + (P - P_ref) / a)^(1 / c); a negative ``simon_a`` gives a
    curve that falls with pressure."""
    MODEL_NAME = "simon_glatzel"


cdef class SimonGlatzel2Curve(MeltingCurveBase):
    """Two Simon and Glatzel branches joined at ``transition_pressure``, as in the Monteux et al. (2016) peridotite
    fits; the defaults are their solidus."""
    MODEL_NAME = "simon_glatzel_2"


cdef class InterpolatedMeltingCurve(MeltingCurveBase):
    """A melting curve tabulated in pressure (aliases ``"interp"``, ``"interpolated"``): ``pressure_pa`` and
    ``temperature_k`` tables."""
    MODEL_NAME = "interpolate"


# The family's name lookup: the C++ registry's alias-aware, case-insensitive canonical name.
_CURVE_FAMILY = ModelFamily(
    "melting curve",
    (ConstantMeltingCurve, SimonGlatzelCurve, SimonGlatzel2Curve, InterpolatedMeltingCurve),
    lambda model_name: c_melting_curve_canonical_name(model_name.encode("utf-8")).decode("utf-8"))


def melting_curve_model_names() -> tuple:
    """The canonical names of the melting curves."""
    return _CURVE_FAMILY.model_names()


def canonical_melting_curve_name(str model_name) -> str:
    """A melting curve's canonical name from any of its names or aliases (case-insensitive)."""
    return _CURVE_FAMILY.canonical_name(model_name)


def make_melting_curve(str model_name, dict config=None) -> MeltingCurveBase:
    """Build a melting curve by name: ``"constant"``, ``"simon_glatzel"``, ``"simon_glatzel_2"``, or
    ``"interpolate"``.

    Raises
    ------
    ValueError
        Unknown model name, or a parameter the curve does not read; each message names the closest accepted one.
    """
    return _CURVE_FAMILY.make(model_name, config)


# =====================================================================================================================
# Melt Weakening
# =====================================================================================================================
cdef class MeltWeakeningBase(PhysicsBase):
    """Base for melt-weakening laws, which give a partially molten aggregate's shear modulus [Pa] and viscosity
    [Pa s] from the solid's and the liquid's values and the melt fraction.

    Instantiate a concrete law with its parameters positionally or as keywords, or build one by name with
    ``make_melt_weakening``.
    """

    # The canonical name of the model a concrete subclass builds; None on this abstract base.
    MODEL_NAME = None

    def __init__(self, *args, dict config=None, **parameters):
        cdef c_ParamMap param_map = cy_model_parameters(type(self), args, config, parameters, "melt weakening")
        cdef unique_ptr[c_MeltWeakeningBase] model = c_find_melt_weakening(
            (<str>type(self).MODEL_NAME).encode("utf-8"), param_map)
        self._set_model(c_share_physics[c_MeltWeakeningBase](move(model)))

    cdef c_MeltWeakeningBase* _weakening(self) except NULL:
        self._check_ptr()
        return <c_MeltWeakeningBase*>self._ptr

    def calc_weakening(
            self,
            double temperature,
            double solidus,
            double liquidus,
            double solid_shear,
            double solid_viscosity,
            double liquid_shear,
            double liquid_viscosity,
            melt_fraction=None,
            solid_shear_at_solidus=None,
            solid_viscosity_at_solidus=None) -> tuple:
        """The aggregate's (shear modulus [Pa], viscosity [Pa s]).

        Parameters
        ----------
        temperature : float
            Temperature [K].
        solidus, liquidus : float
            Solidus and liquidus at the local pressure [K].
        solid_shear, solid_viscosity : float
            The solid phase's shear modulus [Pa] and viscosity [Pa s].
        liquid_shear, liquid_viscosity : float
            The liquid phase's shear modulus [Pa] and viscosity [Pa s], the floors of the result.
        melt_fraction : float, optional
            Melt fraction [m3 m-3]; by default the linear (T - T_sol) / (T_liq - T_sol), clipped to [0, 1].
        solid_shear_at_solidus, solid_viscosity_at_solidus : float, optional
            The solid phase's shear modulus [Pa] and viscosity [Pa s] at the solidus, which a law anchored there
            (Spohn without absolute anchors) continues; by default ``solid_shear`` and ``solid_viscosity``.
        """
        cdef c_MeltWeakeningInputs inputs
        inputs.temperature      = temperature
        inputs.solidus          = solidus
        inputs.liquidus         = liquidus
        inputs.solid_shear      = solid_shear
        inputs.solid_viscosity  = solid_viscosity
        inputs.liquid_shear     = liquid_shear
        inputs.liquid_viscosity = liquid_viscosity
        inputs.solid_shear_at_solidus = solid_shear if solid_shear_at_solidus is None else solid_shear_at_solidus
        inputs.solid_viscosity_at_solidus = (
            solid_viscosity if solid_viscosity_at_solidus is None else solid_viscosity_at_solidus)
        if melt_fraction is None:
            inputs.melt_fraction = min(max((temperature - solidus) / (liquidus - solidus), 0.0), 1.0)
        else:
            inputs.melt_fraction = <double>melt_fraction
        cdef c_MeltWeakeningResult result = self._weakening().calc_weakening(inputs)
        return result.shear_modulus, result.viscosity


cdef class NoMeltWeakening(MeltWeakeningBase):
    """No weakening (alias ``"off"``): the solid's values until fully molten, then the liquid's."""
    MODEL_NAME = "none"


cdef class SpohnMeltWeakening(MeltWeakeningBase):
    """Fischer and Spohn (1990) temperature law (aliases ``"fischer"``, ``"fischer_spohn"``), anchored at the
    solidus on the solid's own values unless given absolute anchors, then blended into the liquid across the breakdown
    band."""
    MODEL_NAME = "spohn"


cdef class HenningMeltWeakening(MeltWeakeningBase):
    """Henning (2009, 2010) three-regime weakening, anchored at the solidus, whose breakdown band blends into the
    liquid."""
    MODEL_NAME = "henning"


# The family's name lookup: the C++ registry's alias-aware, case-insensitive canonical name.
_WEAKENING_FAMILY = ModelFamily(
    "melt weakening",
    (NoMeltWeakening, SpohnMeltWeakening, HenningMeltWeakening),
    lambda model_name: c_melt_weakening_canonical_name(model_name.encode("utf-8")).decode("utf-8"))


def melt_weakening_model_names() -> tuple:
    """The canonical names of the melt-weakening laws."""
    return _WEAKENING_FAMILY.model_names()


def canonical_melt_weakening_name(str model_name) -> str:
    """A melt-weakening law's canonical name from any of its names or aliases (case-insensitive)."""
    return _WEAKENING_FAMILY.canonical_name(model_name)


def make_melt_weakening(str model_name, dict config=None) -> MeltWeakeningBase:
    """Build a melt-weakening law by name: ``"none"`` (``"off"``), ``"spohn"``, or ``"henning"``.

    Raises
    ------
    ValueError
        Unknown model name, or a parameter the law does not read; each message names the closest accepted one.
    """
    return _WEAKENING_FAMILY.make(model_name, config)


# =====================================================================================================================
# Bulk Mixing
# =====================================================================================================================
cdef class BulkModulusMixingBase(PhysicsBase):
    """Base for bulk-modulus mixing laws, which give a partially molten aggregate's bulk modulus [Pa]."""

    # The canonical name of the model a concrete subclass builds; None on this abstract base.
    MODEL_NAME = None

    def __init__(self, *args, dict config=None, **parameters):
        cdef c_ParamMap param_map = cy_model_parameters(type(self), args, config, parameters, "bulk-modulus mixing")
        cdef unique_ptr[c_BulkModulusMixingBase] model = c_find_bulk_modulus_mixing(
            (<str>type(self).MODEL_NAME).encode("utf-8"), param_map)
        self._set_model(c_share_physics[c_BulkModulusMixingBase](move(model)))

    cdef c_BulkModulusMixingBase* _mixing(self) except NULL:
        self._check_ptr()
        return <c_BulkModulusMixingBase*>self._ptr

    def calc_bulk_modulus(
            self,
            double solid_bulk_modulus,
            double liquid_bulk_modulus,
            double framework_shear_modulus,
            double melt_fraction) -> float:
        """The aggregate's bulk modulus [Pa] from the solid's and the liquid's [Pa], the framework's (post-melt)
        shear modulus [Pa], and the melt fraction [m3 m-3]."""
        return self._mixing().calc_bulk_modulus(
            solid_bulk_modulus, liquid_bulk_modulus, framework_shear_modulus, melt_fraction)


cdef class HashinShtrikmanMixing(BulkModulusMixingBase):
    """The Hashin-Shtrikman (1963) bound (aliases ``"hs"``, ``"hashin-shtrikman"``): weak while the framework holds,
    the Reuss average of a suspension once it has collapsed."""
    MODEL_NAME = "hashin_shtrikman"


cdef class BulkViscosityMixingBase(PhysicsBase):
    """Base for bulk-viscosity mixing laws, which give a partially molten aggregate's bulk viscosity [Pa s]."""

    # The canonical name of the model a concrete subclass builds; None on this abstract base.
    MODEL_NAME = None

    def __init__(self, *args, dict config=None, **parameters):
        cdef c_ParamMap param_map = cy_model_parameters(
            type(self), args, config, parameters, "bulk-viscosity mixing")
        cdef unique_ptr[c_BulkViscosityMixingBase] model = c_find_bulk_viscosity_mixing(
            (<str>type(self).MODEL_NAME).encode("utf-8"), param_map)
        self._set_model(c_share_physics[c_BulkViscosityMixingBase](move(model)))

    cdef c_BulkViscosityMixingBase* _mixing(self) except NULL:
        self._check_ptr()
        return <c_BulkViscosityMixingBase*>self._ptr

    def calc_bulk_viscosity(
            self,
            double solid_bulk_viscosity,
            double postmelt_shear_viscosity,
            double melt_fraction) -> float:
        """The aggregate's bulk viscosity [Pa s] from the solid's [Pa s], the post-melt shear viscosity [Pa s], and
        the melt fraction [m3 m-3]."""
        return self._mixing().calc_bulk_viscosity(solid_bulk_viscosity, postmelt_shear_viscosity, melt_fraction)


cdef class CompactionViscosity(BulkViscosityMixingBase):
    """Compaction viscosity (alias ``"mckenzie"``): 1 / zeta = 1 / zeta_solid + phi^n / (c eta)."""
    MODEL_NAME = "compaction"


# The family's name lookup: the C++ registry's alias-aware, case-insensitive canonical name.
_BULK_MODULUS_FAMILY = ModelFamily(
    "bulk-modulus mixing",
    (HashinShtrikmanMixing,),
    lambda model_name: c_bulk_modulus_mixing_canonical_name(model_name.encode("utf-8")).decode("utf-8"))
_BULK_VISCOSITY_FAMILY = ModelFamily(
    "bulk-viscosity mixing",
    (CompactionViscosity,),
    lambda model_name: c_bulk_viscosity_mixing_canonical_name(model_name.encode("utf-8")).decode("utf-8"))


def bulk_modulus_mixing_model_names() -> tuple:
    """The canonical names of the bulk-modulus mixing laws."""
    return _BULK_MODULUS_FAMILY.model_names()


def bulk_viscosity_mixing_model_names() -> tuple:
    """The canonical names of the bulk-viscosity mixing laws."""
    return _BULK_VISCOSITY_FAMILY.model_names()


def make_bulk_modulus_mixing(str model_name, dict config=None) -> BulkModulusMixingBase:
    """Build a bulk-modulus mixing law by name: ``"hashin_shtrikman"`` (``"hs"``)."""
    return _BULK_MODULUS_FAMILY.make(model_name, config)


def make_bulk_viscosity_mixing(str model_name, dict config=None) -> BulkViscosityMixingBase:
    """Build a bulk-viscosity mixing law by name: ``"compaction"`` (``"mckenzie"``)."""
    return _BULK_VISCOSITY_FAMILY.make(model_name, config)
