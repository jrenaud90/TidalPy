# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Phases and materials: the composites that turn TidalPy's property laws into a material's state at a point.

A ``Phase`` is one phase of a material: an equation of state, an optional shear-modulus law, optional shear and bulk
viscosity laws, optional default rheologies, and its thermal conductivity and heat capacity. A ``Material`` is a solid
phase, a liquid phase, or both; with both it melts between its solidus and liquidus, through a melt-weakening law and
optional bulk-mixing laws, and with one it is that phase everywhere. ``Material.calc_state`` gives every property at a
pressure, temperature, and radius, as a layer with the given physics switches sees it.

Each component may be given as a model (``make_eos("vinet", {...})``), a config table with a ``model`` key, or a
model name for its defaults; a whole phase or material may be given as one nested config table, the form
``get_config_dict()`` returns and a TOML material table takes.
"""

import difflib

from libcpp cimport bool as cpp_bool
from libcpp.memory cimport shared_ptr, unique_ptr
from libcpp.string cimport string
from libcpp.utility cimport move
from libcpp.vector cimport vector

cimport numpy as cnp
cnp.import_array()

import numpy as np

from TidalPy.Utilities.logging.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport d_NAN, set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities.arrays.vectors cimport cy_broadcast_inputs
from TidalPy.Utilities.classes.classes cimport (
    PhysicsBase,
    c_ParamMap,
    c_ThermoPoint,
    c_share_physics,
    cy_param_map,
    cy_wrap_model,
)
from TidalPy.Utilities.classes.families import ModelFamily
from TidalPy.Material.laws import make_eos, make_shear_modulus
from TidalPy.Viscosity import make_viscosity
from TidalPy.Rheology import make_rheology
from TidalPy.PartialMelt.melting import (
    make_melting_curve,
    make_melt_weakening,
    make_bulk_modulus_mixing,
    make_bulk_viscosity_mixing,
)

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


# Slot -> the factory that builds its model from a name and a config table.
_PHASE_SLOTS = {
    "eos":             make_eos,
    "shear_modulus":   make_shear_modulus,
    "shear_viscosity": make_viscosity,
    "bulk_viscosity":  make_viscosity,
    "shear_rheology":  make_rheology,
    "bulk_rheology":   make_rheology,
}
_MELTING_SLOTS = {
    "solidus":               make_melting_curve,
    "liquidus":              make_melting_curve,
    "weakening":             make_melt_weakening,
    "bulk_modulus_mixing":   make_bulk_modulus_mixing,
    "bulk_viscosity_mixing": make_bulk_viscosity_mixing,
}

# The fields Material.calc_state reports, in c_MaterialState order (phase first, as a string).
_STATE_FIELDS = (
    "density", "bulk_modulus", "adiabatic_bulk_modulus", "thermal_expansion", "heat_capacity", "thermal_conductivity",
    "shear_modulus", "shear_viscosity", "bulk_viscosity", "melt_fraction", "solidus", "liquidus")
# c_MaterialPhase values, in order.
_PHASE_NAMES = ("solid", "partial", "liquid")

# Composite class -> the parameter names and config keys it reads, filled on first use from a default instance.
_PARAMETER_KEYS = {}


def _parameter_keys(object composite_class) -> frozenset:
    if composite_class not in _PARAMETER_KEYS:
        _PARAMETER_KEYS[composite_class] = frozenset(
            key for entry in composite_class().get_parameter_info() for key in (entry["name"], entry["key"]))
    return _PARAMETER_KEYS[composite_class]


def _as_model(object value, object factory, str slot):
    """A component as a model: a model passes through, a name builds that model's defaults, and a table builds the
    model its ``model`` key names."""
    if value is None or isinstance(value, PhysicsBase):
        return value
    if isinstance(value, str):
        return factory(value, {})
    if isinstance(value, dict):
        table = dict(value)
        if "model" not in table:
            raise ValueError(f"TidalPy: the '{slot}' table needs a 'model' key naming its model.")
        return factory(str(table.pop("model")), table)
    raise TypeError(f"TidalPy: the '{slot}' slot takes a model, a model name, or a config table, "
                    f"not {type(value).__name__}.")


cdef dict cy_split_config(dict config, object slots, object parameter_keys, str family):
    """Split a composite's config table into its slot tables (returned) and its parameters (left in ``config``).

    Raises ValueError for a key that is neither, naming the closest slot or parameter.
    """
    cdef dict slot_tables = {}
    cdef str key
    for key in list(config):
        if key in slots:
            slot_tables[key] = config.pop(key)
        elif key not in parameter_keys:
            accepted = sorted(set(slots) | set(parameter_keys))
            close = difflib.get_close_matches(key, accepted, n=1)
            hint = f" (did you mean '{close[0]}'?)" if close else ""
            raise ValueError(f"TidalPy: a {family} has no slot or parameter '{key}'{hint}. "
                             f"Accepted: {', '.join(accepted)}.")
    return slot_tables


cdef c_MaterialSwitches cy_switches(object use_thermal_expansion, object use_melting, object use_pressure_melting,
                                    object use_melt_density):
    cdef c_MaterialSwitches switches
    switches.use_thermal_expansion = True if use_thermal_expansion else False
    switches.use_melting           = True if use_melting else False
    switches.use_pressure_melting  = True if use_pressure_melting else False
    switches.use_melt_density      = True if use_melt_density else False
    return switches


# =====================================================================================================================
# Phase
# =====================================================================================================================
cdef class Phase(PhysicsBase):
    """One phase of a material: an equation of state, optional shear-modulus and viscosity laws, optional default
    rheologies, and its thermal conductivity and heat capacity.

    Parameters
    ----------
    eos : EOSBase, str, or dict, optional
        Equation of state; a constant-density law with its defaults when absent.
    shear_modulus : ShearModulusBase, str, or dict, optional
        Static shear modulus; none (a fluid) when absent.
    shear_viscosity, bulk_viscosity : ViscosityBase, str, or dict, optional
        Viscosity laws; none (NaN) when absent.
    shear_rheology, bulk_rheology : RheologyBase, str, or dict, optional
        The rheologies a layer of this phase uses when it sets none of its own.
    config : dict, optional
        The whole phase as one table: slot tables plus the thermal parameters by config key.
    **parameters
        Thermal parameters by argument name or config key (``thermal_conductivity``, ``heat_capacity``, and their
        temperature exponents and reference temperature; see ``get_parameter_info()``).
    """

    MODEL_NAME = "phase"

    def __init__(self, eos=None, shear_modulus=None, shear_viscosity=None, bulk_viscosity=None,
                 shear_rheology=None, bulk_rheology=None, *, dict config=None, **parameters):
        cdef dict given = {"eos": eos, "shear_modulus": shear_modulus, "shear_viscosity": shear_viscosity,
                           "bulk_viscosity": bulk_viscosity, "shear_rheology": shear_rheology,
                           "bulk_rheology": bulk_rheology}
        cdef dict params = dict(config) if config else {}
        cdef dict slot_tables = cy_split_config(params, _PHASE_SLOTS, _parameter_keys(Phase), "phase") if params else {}
        params.update(parameters)
        cdef c_PhaseComponents components
        cdef PhysicsBase model
        cdef str slot
        for slot, factory in _PHASE_SLOTS.items():
            value = given[slot] if given[slot] is not None else slot_tables.get(slot)
            model = _as_model(value, factory, slot)
            if model is not None:
                components.set(slot.encode("utf-8"), model._model_sptr)
        cdef unique_ptr[c_Phase] phase = c_make_phase(cy_param_map(params), components)
        self._set_model(c_share_physics[c_Phase](move(phase)))

    cdef c_Phase* _phase(self) except NULL:
        self._check_ptr()
        return <c_Phase*>self._ptr

    def _component(self, str slot):
        return cy_wrap_model(self._phase().get_components().get(slot.encode("utf-8")))

    @property
    def eos(self):
        """The equation of state."""
        return self._component("eos")

    @property
    def shear_modulus(self):
        """The shear-modulus law, or None."""
        return self._component("shear_modulus")

    @property
    def shear_viscosity(self):
        """The shear-viscosity law, or None."""
        return self._component("shear_viscosity")

    @property
    def bulk_viscosity(self):
        """The bulk-viscosity law, or None."""
        return self._component("bulk_viscosity")

    @property
    def shear_rheology(self):
        """The default shear rheology, or None."""
        return self._component("shear_rheology")

    @property
    def bulk_rheology(self):
        """The default bulk rheology, or None."""
        return self._component("bulk_rheology")

    def calc_state(self, double pressure, double temperature=d_NAN, double radius=d_NAN,
                   use_thermal_expansion=False) -> dict:
        """The phase's own properties at a point (no melting): density [kg m-3], bulk moduli [Pa], expansivity
        [1/K], heat capacity [J kg-1 K-1], conductivity [W m-1 K-1], shear modulus [Pa], viscosities [Pa s]."""
        cdef c_ThermoPoint point
        cdef c_PhaseState state
        point.pressure    = pressure
        point.temperature = temperature
        point.radius      = radius
        self._phase().calc_phase_state(point, True if use_thermal_expansion else False, state)
        return {
            "density": state.density, "bulk_modulus": state.bulk_modulus,
            "adiabatic_bulk_modulus": state.adiabatic_bulk_modulus, "thermal_expansion": state.thermal_expansion,
            "heat_capacity": state.heat_capacity, "thermal_conductivity": state.thermal_conductivity,
            "shear_modulus": state.shear_modulus, "shear_viscosity": state.shear_viscosity,
            "bulk_viscosity": state.bulk_viscosity}


# =====================================================================================================================
# Material
# =====================================================================================================================
cdef class Material(PhysicsBase):
    """A material: a solid phase, a liquid phase, or both. With both it melts between its solidus and liquidus,
    through a melt-weakening law and optional bulk-mixing laws; with one it is that phase everywhere (a water ocean is
    a liquid-only material).

    Parameters
    ----------
    solid : Phase or dict, optional
        The solid phase; a default phase when neither phase is given.
    liquid : Phase or dict, optional
        The liquid phase. With a solid phase it is the melt (its ``shear_viscosity`` is the melt's), and without one
        the material is liquid everywhere.
    solidus, liquidus : MeltingCurveBase, str, or dict, optional
        Required with both phases; equal curves give a single melting temperature.
    weakening : MeltWeakeningBase, str, or dict, optional
        How the shear modulus and viscosity fall with melt; none (the solid's until fully molten) when absent.
    bulk_modulus_mixing, bulk_viscosity_mixing : str or dict or model, optional
        How melt changes the bulk modulus and bulk viscosity; the solid's until fully molten when absent.
    config : dict, optional
        The whole material as one table: ``solid`` and ``liquid`` tables, a ``melting`` table holding the melting
        slots, and the material's parameters by config key.
    **parameters
        ``latent_heat`` [J kg-1] (config key ``latent_heat_j_kg``).
    """

    MODEL_NAME = "material"

    def __init__(self, solid=None, liquid=None, solidus=None, liquidus=None, weakening=None,
                 bulk_modulus_mixing=None, bulk_viscosity_mixing=None, *, dict config=None, **parameters):
        cdef dict given = {"solidus": solidus, "liquidus": liquidus, "weakening": weakening,
                           "bulk_modulus_mixing": bulk_modulus_mixing, "bulk_viscosity_mixing": bulk_viscosity_mixing}
        cdef dict params = dict(config) if config else {}
        cdef dict tables = cy_split_config(
            params, ("solid", "liquid", "melting"), _parameter_keys(Material), "material") if params else {}
        cdef dict melting = dict(tables.get("melting") or {})
        cy_split_config(dict(melting), _MELTING_SLOTS, (), "material's melting table")
        params.update(parameters)
        cdef c_MaterialComponents components
        cdef PhysicsBase model
        cdef str slot
        for slot, value in (("solid", solid), ("liquid", liquid)):
            if value is None:
                value = tables.get(slot)
            if isinstance(value, dict):
                value = Phase(config=value)
            if value is not None and not isinstance(value, Phase):
                raise TypeError(f"TidalPy: the material's '{slot}' slot takes a Phase or a phase config table, "
                                f"not {type(value).__name__}.")
            if value is not None:
                model = <PhysicsBase>value
                components.set(slot.encode("utf-8"), model._model_sptr)
        for slot, factory in _MELTING_SLOTS.items():
            value = given[slot] if given[slot] is not None else melting.get(slot)
            model = _as_model(value, factory, slot)
            if model is not None:
                components.set(slot.encode("utf-8"), model._model_sptr)
        cdef unique_ptr[c_Material] material = c_make_material(cy_param_map(params), components)
        self._set_model(c_share_physics[c_Material](move(material)))

    cdef c_Material* _material(self) except NULL:
        self._check_ptr()
        return <c_Material*>self._ptr

    def _component(self, str slot):
        return cy_wrap_model(self._material().get_components().get(slot.encode("utf-8")))

    @property
    def solid(self):
        """The solid phase, or None for a liquid-only material."""
        return self._component("solid")

    @property
    def liquid(self):
        """The liquid phase, or None."""
        return self._component("liquid")

    @property
    def solidus(self):
        """The solidus curve, or None."""
        return self._component("solidus")

    @property
    def liquidus(self):
        """The liquidus curve, or None."""
        return self._component("liquidus")

    @property
    def weakening(self):
        """The melt-weakening law, or None."""
        return self._component("weakening")

    @property
    def bulk_modulus_mixing(self):
        """The bulk-modulus mixing law, or None."""
        return self._component("bulk_modulus_mixing")

    @property
    def bulk_viscosity_mixing(self):
        """The bulk-viscosity mixing law, or None."""
        return self._component("bulk_viscosity_mixing")

    @property
    def can_melt(self) -> bool:
        """Whether the material has a solid and a liquid phase to melt between."""
        return True if self._material().get_can_melt() else False

    @property
    def is_liquid_only(self) -> bool:
        """Whether the material has a liquid phase and no solid one, so it is liquid everywhere."""
        return True if self._material().get_is_liquid_only() else False

    def replace(self, **changes):
        """A new material with some components replaced (``None`` removes one); this one is unchanged.

        ``changes`` take the slot names (``solid``, ``liquid``, ``solidus``, ...) and the material's parameters. A
        change that would leave neither a solid nor a liquid phase raises ValueError.
        """
        cdef dict components = {slot: self._component(slot) for slot in ("solid", "liquid", *_MELTING_SLOTS)}
        cdef dict parameters = dict(self.parameters)
        for key, value in changes.items():
            if key in components:
                components[key] = value
            else:
                parameters[key] = value
        if components["solid"] is None and components["liquid"] is None:
            raise ValueError("TidalPy: the material would have neither a 'solid' nor a 'liquid' phase.")
        return Material(**{slot: value for slot, value in components.items() if value is not None}, **parameters)

    def calc_melting_range(self, double pressure, use_pressure_melting=True) -> tuple:
        """The (solidus, liquidus) [K] at a pressure [Pa]; read at zero pressure when ``use_pressure_melting`` is
        off. NaN for a material that cannot melt."""
        cdef double solidus = d_NAN
        cdef double liquidus = d_NAN
        cdef c_MaterialSwitches switches = cy_switches(False, True, use_pressure_melting, False)
        self._material().calc_melting_range(pressure, switches, solidus, liquidus)
        return solidus, liquidus

    def calc_state(self, pressure, temperature=d_NAN, radius=d_NAN, *, use_thermal_expansion=False,
                   use_melting=False, use_pressure_melting=False, use_melt_density=False) -> dict:
        """Every property of the material at a point, as a layer with these switches sees it.

        Parameters
        ----------
        pressure, temperature, radius : float or np.ndarray
            Pressure [Pa], temperature [K] (NaN: no melt state), and radius [m] (read only by tabulated laws),
            broadcast together.
        use_thermal_expansion, use_melting, use_pressure_melting, use_melt_density : bool, optional
            The layer physics switches; each off, the default, is the simpler case.

        Returns
        -------
        dict
            ``phase`` ("solid", "partial", "liquid"), ``density`` [kg m-3], ``bulk_modulus`` and
            ``adiabatic_bulk_modulus`` [Pa], ``thermal_expansion`` [1/K], ``heat_capacity`` (latent heat included)
            [J kg-1 K-1], ``thermal_conductivity`` [W m-1 K-1], ``shear_modulus`` [Pa], ``shear_viscosity`` and
            ``bulk_viscosity`` [Pa s], ``melt_fraction``, ``solidus`` and ``liquidus`` [K]; floats (and a string) for
            float inputs, otherwise arrays of the broadcast shape.
        """
        cdef c_Material* material_ptr = self._material()
        cdef c_MaterialSwitches switches = cy_switches(
            use_thermal_expansion, use_melting, use_pressure_melting, use_melt_density)
        cdef vector[vector[double]] inputs
        cdef vector[c_MaterialState] states
        cdef object shape = cy_broadcast_inputs((pressure, temperature, radius), inputs, False)
        cdef c_ThermoPoint point
        cdef c_MaterialState scalar_state
        if shape is None:
            point.pressure    = <double>pressure
            point.temperature = <double>temperature
            point.radius      = <double>radius
            material_ptr.calc_state(point, switches, scalar_state)
            return {
                "phase": _PHASE_NAMES[<int>scalar_state.phase],
                "density": scalar_state.density,
                "bulk_modulus": scalar_state.bulk_modulus,
                "adiabatic_bulk_modulus": scalar_state.adiabatic_bulk_modulus,
                "thermal_expansion": scalar_state.thermal_expansion,
                "heat_capacity": scalar_state.heat_capacity,
                "thermal_conductivity": scalar_state.thermal_conductivity,
                "shear_modulus": scalar_state.shear_modulus,
                "shear_viscosity": scalar_state.shear_viscosity,
                "bulk_viscosity": scalar_state.bulk_viscosity,
                "melt_fraction": scalar_state.melt_fraction,
                "solidus": scalar_state.solidus,
                "liquidus": scalar_state.liquidus,
            }
        with nogil:
            material_ptr.calc_state_vectorize(inputs[0], inputs[1], inputs[2], switches, states)
        cdef size_t num_points = states.size()
        cdef size_t point_i
        cdef cnp.ndarray[cnp.float64_t, ndim=2] values = np.empty((len(_STATE_FIELDS), num_points), dtype=np.float64)
        phases = []
        for point_i in range(num_points):
            values[0, point_i]  = states[point_i].density
            values[1, point_i]  = states[point_i].bulk_modulus
            values[2, point_i]  = states[point_i].adiabatic_bulk_modulus
            values[3, point_i]  = states[point_i].thermal_expansion
            values[4, point_i]  = states[point_i].heat_capacity
            values[5, point_i]  = states[point_i].thermal_conductivity
            values[6, point_i]  = states[point_i].shear_modulus
            values[7, point_i]  = states[point_i].shear_viscosity
            values[8, point_i]  = states[point_i].bulk_viscosity
            values[9, point_i]  = states[point_i].melt_fraction
            values[10, point_i] = states[point_i].solidus
            values[11, point_i] = states[point_i].liquidus
            phases.append(_PHASE_NAMES[<int>states[point_i].phase])
        cdef dict result = {field: values[field_i].reshape(shape) for field_i, field in enumerate(_STATE_FIELDS)}
        result["phase"] = np.array(phases, dtype=object).reshape(shape)
        return result


def _single_model_name(str model_name, str family) -> str:
    """The one model name a composite family has, or a ValueError naming it."""
    if model_name.lower() != family:
        close = difflib.get_close_matches(model_name.lower(), [family], n=1)
        hint = f" (did you mean '{family}'?)" if close else ""
        raise ValueError(f"TidalPy: unknown {family} model name '{model_name}'{hint}. Accepted: {family}.")
    return family


def _phase_canonical_name(str model_name) -> str:
    return _single_model_name(model_name, "phase")


def _material_canonical_name(str model_name) -> str:
    return _single_model_name(model_name, "material")


_PHASE_FAMILY = ModelFamily("phase", (Phase,), _phase_canonical_name)
_MATERIAL_FAMILY = ModelFamily("material", (Material,), _material_canonical_name)


def make_phase(dict config) -> Phase:
    """A phase from its config table (slot tables plus thermal parameters), the form ``get_config_dict()`` returns."""
    return Phase(config=config)


def make_material(dict config) -> Material:
    """A material from its config table (``solid``, ``liquid``, ``melting``, and parameters), the form
    ``get_config_dict()`` returns and a TOML material table takes."""
    return Material(config=config)
