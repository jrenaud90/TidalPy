# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrapper for TidalPy's base world class.

BaseWorld owns an ordered stack of layers (inner to outer, possibly empty) and holds the world-level identity, the
orbital and thermal scalars (albedo, emissivity, obliquity, spin frequency), bulk geometry, the equilibrium
temperature, the spin model, and the tide model. The whole-planet EOS, radial (Love number), and tidal solves are
methods on this class. TerrestrialWorld, GasGiantWorld, and StarWorld subclass it.
"""

cimport numpy as cnp
cnp.import_array()

import contextlib
import copy
import math
import os
import textwrap
import weakref

import numpy as np

from libc.stdlib cimport malloc, free
from libcpp cimport bool as cpp_bool
from libcpp.utility cimport move
from libcpp.memory cimport make_shared
from libcpp.complex cimport complex as cpp_complex
from libcpp.vector cimport vector
from cpython.complex cimport PyComplex_FromDoubles
from cython.operator cimport dereference as deref

from TidalPy.Utilities.logging.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address, d_PI
from TidalPy.Utilities.classes.classes cimport (
    StructureBase,
    c_TidalPyBaseClass,
    c_PhysicsBase,
    c_binary_class_name,
    cy_physics_model_config,
    cy_wrap_model,
)
from TidalPy.Utilities.classes.classes import did_you_mean
from TidalPy.Tides.classes.tide cimport TideBase
from TidalPy.Tides.classes.tide import make_tide, tide_config_keys, TIDE_CONFIG_KEYS
from TidalPy.Structures.layers.layer cimport (
    Layer, c_Layer, cy_eos_field, cy_eos_fields, C_EOS_DENSITY_INDEX,
    C_EOS_GRAVITY_INDEX, C_EOS_PRESSURE_INDEX, C_EOS_SHEAR_MODULUS_INDEX, C_EOS_SHEAR_VISCOSITY_INDEX,
    C_EOS_BULK_MODULUS_INDEX, C_EOS_BULK_VISCOSITY_INDEX, C_EOS_TEMPERATURE_INDEX, C_EOS_HEAT_FLOW_INDEX,
    C_EOS_MELT_FRACTION_INDEX)
from TidalPy.Structures.layers.layer import LAYER_STANDALONE_CONFIG_KEYS
from TidalPy.RadialSolver.rs_constants cimport C_MAX_NUM_YTYPES
from TidalPy.RadialSolver.rs_solution cimport RadialSolverSolution
from TidalPy.RadialSolver.rs_solution cimport cy_check_surface_solve_conditioning, cy_check_degree1_frame_residual
from TidalPy.Tides.love.love cimport (
    c_love_method_name_int, c_love_method_uses_radial_solver_int, cy_parse_love_method)
from TidalPy.Tides.eccentricity.eccentricity_driver import (
    eccentricity_truncation_name, validate_eccentricity_exact_tolerance, validate_eccentricity_truncation)
from TidalPy.Tides.obliquity.obliquity_driver import obliquity_truncation_name, validate_obliquity_truncation
from TidalPy.Utilities.logging.logger import log_warning
from TidalPy.constants import ODE_METHOD_NAMES, ode_method_from_name, STARTING_METHOD_NAMES, starting_method_from_name
from TidalPy.constants import DEGREE1_FRAME_NAMES, degree1_frame_from_name
from TidalPy.exceptions import SolutionFailedError

# Pull in the out-of-line definition of c_BaseWorld::calc_tides and the 3D tidal paths, with the heavy
# global-potential engine they use, so they compile into this extension.
cdef extern from "world_tides_.hpp" nogil:
    pass

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())

# World ``type`` values the TOML builder accepts (kept in step with configs.toml_loader.WORLD_TYPES, which
# cannot be imported here at module load without a circular import).
BUILDER_WORLD_TYPES = ("star", "gasgiant", "terrestrial", "layered")

# A layer wrapper as a non-owning view onto a layer the world owns. The view keeps the world alive (see Layer._view).
cdef Layer cy_wrap_layer_view(c_Layer* ptr, object world):
    return Layer._view(ptr, world)

# Component order of the stress and strain grids returned by BaseWorld.calc_3d_stress_strain.
STRESS_STRAIN_COMPONENTS = ("rr", "theta_theta", "phi_phi", "r_theta", "r_phi", "theta_phi")


# Copy a C++ vector[double] into a new 1D float64 ndarray.
cdef cnp.ndarray cy_vec_to_ndarray(const vector[double]& v):
    cdef Py_ssize_t n = <Py_ssize_t>v.size()
    cdef cnp.ndarray out = np.empty(n, dtype=np.float64)
    cdef double[::1] mv
    cdef Py_ssize_t i
    if n > 0:
        mv = out
        for i in range(n):
            mv[i] = v[i]
    return out

# The names of a world's layers, inner to outer.
cdef list cy_layer_names(c_BaseWorld* world_ptr):
    cdef size_t layer_i
    return [world_ptr.get_layer(layer_i).get_name().decode("utf-8") for layer_i in range(world_ptr.get_num_layers())]

# A layer of the world by name or integer index (a negative one counts from the outermost layer), as its index. Every
# method that takes a layer resolves it here: KeyError, listing the layer names, for a name the world does not have;
# IndexError for an index out of range; TypeError for anything else.
cdef size_t cy_layer_index(c_BaseWorld* world_ptr, object layer) except? 0:
    cdef list names = cy_layer_names(world_ptr)
    cdef Py_ssize_t num_layers = len(names)
    cdef Py_ssize_t index
    if isinstance(layer, str):
        if layer not in names:
            raise KeyError(
                f"TidalPy: world '{world_ptr.get_name().decode('utf-8')}' has no layer named '{layer}'"
                f"{did_you_mean(layer, names)}; its layers are {names}.")
        return <size_t>names.index(layer)
    if isinstance(layer, bool) or not isinstance(layer, (int, np.integer)):
        raise TypeError(f"TidalPy: a layer is given by its name or an integer index, not {layer!r}.")
    index = int(layer)
    if index < 0:
        index += num_layers
    if not (0 <= index < num_layers):
        raise IndexError(
            f"TidalPy: layer index {int(layer)} is out of range for world "
            f"'{world_ptr.get_name().decode('utf-8')}', which has {num_layers} layers.")
    return <size_t>index


# The orbital-state arguments of calc_tides and the 3D tide methods, in their order (also get_tide_state's keys).
TIDE_STATE_ARGUMENTS = (
    "orbital_frequency",
    "spin_frequency",
    "eccentricity",
    "obliquity",
    "semi_major_axis",
    "host_mass",
)

# The keys set_tide_config takes: get_tide_config's, and the older spellings of three of them.
TIDE_SETTING_ALIASES = {
    "eccentricity_truncation": "eccentricity_trunc_lvl",
    "obliquity_truncation": "obliquity_trunc_lvl",
    "love_fixed_dt": "love_fixed_dt_s",
}
TIDE_SETTING_KEYS = frozenset((
    "min_degree_l",
    "max_degree_l",
    "eccentricity_trunc_lvl",
    "eccentricity_exact_tolerance",
    "obliquity_trunc_lvl",
    "layer_tidal_heating",
    "love_method",
    "love_fixed_q",
    "love_fixed_dt_s",
)) | frozenset(TIDE_SETTING_ALIASES)

# The keys of a tide model table that name its model (a TOML [tides] table uses the second).
TIDE_MODEL_NAME_KEYS = ("model", "global_tidal_model")

# Points per layer at which BaseWorld.summary reads the solved density for its range.
SUMMARY_DENSITY_SAMPLES = 9


# The orbital state of a tide solve from the six orbital-state arguments in their order: a value given is used, and
# each None comes from the world's get_tide_state (the system the world belongs to). ValueError names the arguments
# still missing, for a world outside a system or without a tidal host; `method` names the caller in it.
cdef c_TideSolveConfig cy_resolve_tide_state(
        BaseWorld world,
        str method,
        object orbital_frequency,
        object spin_frequency,
        object eccentricity,
        object obliquity,
        object semi_major_axis,
        object host_mass) except *:
    cdef list values = [orbital_frequency, spin_frequency, eccentricity, obliquity, semi_major_axis, host_mass]
    cdef object system_state = None
    cdef list missing = []
    cdef Py_ssize_t argument_i
    cdef cpp_bool any_missing = False
    for argument_i in range(len(values)):
        if values[argument_i] is None:
            any_missing = True
    if any_missing:
        system_state = world.get_tide_state()
        for argument_i in range(len(values)):
            if values[argument_i] is not None:
                continue
            if system_state is None:
                missing.append(TIDE_STATE_ARGUMENTS[argument_i])
            else:
                values[argument_i] = system_state[TIDE_STATE_ARGUMENTS[argument_i]]
        if missing:
            raise ValueError(
                f"TidalPy: world '{world.name}' {method} needs {', '.join(missing)}: pass them, or add the world to a "
                f"System with a tidal host so they come from its orbit (see get_tide_state).")
    return cy_tide_state(
        <double>values[0],
        <double>values[1],
        <double>values[2],
        <double>values[3],
        <double>values[4],
        <double>values[5])

# The solid and liquid zones of a solve as dicts: the layer's name, the radii [m] and enclosed masses [kg] at the zone's
# two ends, and its state ("solid" or "liquid").
cdef list cy_zones_to_list(const vector[c_EOSZone]& zones, list layer_names):
    cdef size_t zone_i
    cdef list out = []
    for zone_i in range(zones.size()):
        out.append({
            'layer':        layer_names[zones[zone_i].layer_index],
            'radius_inner': zones[zone_i].radius_inner,
            'radius_outer': zones[zone_i].radius_outer,
            'mass_inner':   zones[zone_i].mass_inner,
            'mass_outer':   zones[zone_i].mass_outer,
            'state':        'liquid' if zones[zone_i].liquid else 'solid'})
    return out

# Width of the key list in an EOSResult summary [characters].
EOS_RESULT_REPR_WIDTH = 120


class EOSResult(dict):
    """The result of :meth:`BaseWorld.solve_eos`: a ``dict`` in every respect (item access, equality with a plain
    ``dict``, copying, pickling, ``isinstance(result, dict)``) whose ``repr`` is a short summary instead of every
    profile array, so a notebook cell ending in a solve prints a few lines.
    """

    def __repr__(self):
        radius = self.get("radius")
        # Not radius[-1]: this module compiles with wraparound=False, under which a negative index is undefined.
        surface_radius = float(radius[len(radius) - 1]) if (radius is not None) and (len(radius) > 0) else math.nan
        summary = (
            f"EOSResult(success={self.get('success')}, iterations={self.get('iterations')}, "
            f"message={self.get('message')!r})\n"
            f"    planet_mass = {self.get('planet_mass', math.nan):.6g} kg, radius = {surface_radius:.6g} m, "
            f"central_pressure = {self.get('central_pressure', math.nan):.6g} Pa\n")
        return summary + textwrap.fill(
            ", ".join(str(key) for key in self), width=EOS_RESULT_REPR_WIDTH, initial_indent="    keys: ",
            subsequent_indent="          ")


# The solve_eos result from a report copied under the world's call lock (c_BaseWorld::get_eos_report).
cdef object cy_eos_report_to_dict(const c_WorldEOSReport& report, list layer_names):
    cdef size_t j
    cdef size_t num_layers = report.layer_thermal.size()
    cdef list layer_temperature        = []
    cdef list layer_heat_flow_in       = []
    cdef list layer_heat_flow_out      = []
    cdef list layer_heating            = []
    cdef list layer_node_temperature   = []
    cdef list layer_top_temperature    = []
    cdef list layer_base_temperature   = []
    cdef list layer_boundary_thickness = []
    cdef list layer_rayleigh_number    = []
    cdef list layer_nusselt_number     = []
    cdef list layer_boundary_fallback  = []
    cdef list layer_magma_ocean        = []
    cdef list layer_ref_pressure       = []
    cdef list layer_ref_viscosity      = []
    cdef list layer_ref_melt_fraction  = []
    cdef list layer_in_thermal_network = []
    cdef list layer_heating_radiogenic = []
    cdef list layer_heating_tidal = []
    cdef list layer_heating_prescribed = []
    for j in range(num_layers):
        layer_temperature.append(report.layer_thermal[j].temperature)
        layer_heat_flow_in.append(report.layer_thermal[j].heat_flow_in)
        layer_heat_flow_out.append(report.layer_thermal[j].heat_flow_out)
        layer_heating.append(report.layer_thermal[j].heating)
        layer_node_temperature.append(report.layer_thermal[j].node_temperature)
        layer_top_temperature.append(report.layer_thermal[j].top_temperature)
        layer_base_temperature.append(report.layer_thermal[j].base_temperature)
        layer_boundary_thickness.append(report.layer_thermal[j].boundary_thickness)
        layer_rayleigh_number.append(report.layer_thermal[j].rayleigh_number)
        layer_nusselt_number.append(report.layer_thermal[j].nusselt_number)
        layer_boundary_fallback.append(bool(report.layer_thermal[j].boundary_fallback))
        layer_magma_ocean.append(bool(report.layer_thermal[j].magma_ocean))
        layer_ref_pressure.append(report.layer_thermal[j].reference_pressure)
        layer_ref_viscosity.append(report.layer_thermal[j].reference_viscosity)
        layer_ref_melt_fraction.append(report.layer_thermal[j].reference_melt_fraction)
        layer_in_thermal_network.append(bool(report.layer_thermal[j].in_network))
        layer_heating_radiogenic.append(
            report.layer_thermal[j].heating_by_source[<size_t>c_HeatSourceKind.Radiogenic])
        layer_heating_tidal.append(report.layer_thermal[j].heating_by_source[<size_t>c_HeatSourceKind.Tidal])
        layer_heating_prescribed.append(
            report.layer_thermal[j].heating_by_source[<size_t>c_HeatSourceKind.Prescribed])

    return EOSResult({
        'success':          bool(report.success),
        'message':          report.message.decode('utf-8'),
        'iterations':       report.iterations,
        'structure_integrations': report.structure_integrations,
        'max_iters_hit':    bool(report.max_iters_hit),
        'pressure_error':   report.pressure_error,
        'radius':           cy_vec_to_ndarray(report.radius),
        'gravity':          cy_vec_to_ndarray(report.gravity),
        'pressure':         cy_vec_to_ndarray(report.pressure),
        'mass':             cy_vec_to_ndarray(report.mass),
        'moi':              cy_vec_to_ndarray(report.moi),
        'density':          cy_vec_to_ndarray(report.density),
        'surface_gravity':  report.surface_gravity,
        'surface_pressure': report.surface_pressure,
        'central_pressure': report.central_pressure,
        'planet_mass':      report.planet_mass,
        'planet_moi':       report.planet_moi,
        'temperature':      cy_vec_to_ndarray(report.temperature),
        'heat_flow':        cy_vec_to_ndarray(report.heat_flow),
        'thermal_passes':   report.thermal_passes,
        'thermal_converged':      bool(report.thermal_converged),
        'zones':                  cy_zones_to_list(report.zones, layer_names),
        'layer_radius_outer':     list(report.layer_radius_outer),
        'layer_temperature':      layer_temperature,
        'layer_heat_flow_in':     layer_heat_flow_in,
        'layer_heat_flow_out':    layer_heat_flow_out,
        'layer_heating':          layer_heating,
        'layer_heating_radiogenic': layer_heating_radiogenic,
        'layer_heating_tidal':      layer_heating_tidal,
        'layer_heating_prescribed': layer_heating_prescribed,
        'layer_node_temperature':   layer_node_temperature,
        'layer_top_temperature':    layer_top_temperature,
        'layer_base_temperature':   layer_base_temperature,
        'layer_boundary_thickness': layer_boundary_thickness,
        'layer_rayleigh_number':    layer_rayleigh_number,
        'layer_nusselt_number':     layer_nusselt_number,
        'layer_boundary_fallback':  layer_boundary_fallback,
        'layer_magma_ocean':        layer_magma_ocean,
        'layer_reference_pressure':      layer_ref_pressure,
        'layer_reference_viscosity':     layer_ref_viscosity,
        'layer_reference_melt_fraction': layer_ref_melt_fraction,
        'layer_in_thermal_network': layer_in_thermal_network,
    })


def build_layered_world_from_profile(
        double[::1] radius not None,
        double[::1] density not None,
        double[::1] shear_modulus not None,
        double[::1] bulk_modulus not None,
        double[::1] upper_radius_bylayer not None,
        layer_is_solid,
        layer_is_static,
        layer_is_incompressible,
        double planet_bulk_density,
        str name = 'radial_solver_profile'):
    """Build a :class:`BaseWorld` whose layers interpolate their own slice of a radial profile.

    The layers are built in C++ (``c_build_world_from_layered_profile``), which is also what the standalone
    ``RadialSolver.radial_solver`` calls, so a world built here and the temporary one that API solves are
    the same world. Interface radii appear twice in the profile, once as the top of the lower layer and once
    as the base of the upper one, and each copy belongs to its own layer.

    The world carries its layers and their materials and nothing else: no tide model, no
    ``[worlds]`` defaults, and no retained source configuration. It is the world a supplied profile
    describes, not one a configuration file asked for; :func:`build_world` is the route that adds those.
    ``save_to_toml`` still works, rebuilding the configuration from the live world.

    Parameters
    ----------
    radius, density, shear_modulus, bulk_modulus : np.ndarray[float64]
        The profile [m, kg m-3, Pa, Pa], ascending in radius. The moduli are the static (unrelaxed) ones; a
        viscoelastic response is supplied separately to the Love solve.
    upper_radius_bylayer : np.ndarray[float64]
        Upper radius of each layer [m], inner to outer.
    layer_is_solid, layer_is_static, layer_is_incompressible : sequence of bool
        Per-layer radial-solver assumptions.
    planet_bulk_density : float
        Bulk density [kg m-3], which fixes the world's mass and so its non-dimensional scales.
    name : str, optional
        Name for the constructed world.

    Returns
    -------
    BaseWorld

    Raises
    ------
    ValueError
        If the arrays disagree in length, or a layer would hold fewer than two profile points.
    """
    cdef size_t num_slices = <size_t>radius.shape[0]
    cdef size_t num_layers = <size_t>upper_radius_bylayer.shape[0]
    if num_slices == 0:
        raise ValueError("radius must not be empty.")
    if num_layers == 0:
        raise ValueError("upper_radius_bylayer must name at least one layer.")
    if (<size_t>density.shape[0] != num_slices or <size_t>shear_modulus.shape[0] != num_slices
            or <size_t>bulk_modulus.shape[0] != num_slices):
        raise ValueError(
            f"density, shear_modulus, and bulk_modulus must all match the radius length ({num_slices}); "
            f"got {density.shape[0]}, {shear_modulus.shape[0]}, {bulk_modulus.shape[0]}.")
    if (len(layer_is_solid) != num_layers or len(layer_is_static) != num_layers
            or len(layer_is_incompressible) != num_layers):
        raise ValueError(
            f"layer_is_solid, layer_is_static, and layer_is_incompressible must each have one entry per "
            f"layer ({num_layers}).")

    # The C++ builder takes the layer-type encoding the radial-solver input check produces (0 for solid), and
    # malloc'd flags because std::vector<bool> is bit-packed and so has no bool* to hand it.
    cdef vector[int] layer_type_vec = vector[int](num_layers)
    cdef cpp_bool* is_static_ptr = <cpp_bool*>malloc(num_layers * sizeof(cpp_bool))
    cdef cpp_bool* is_incomp_ptr = <cpp_bool*>malloc(num_layers * sizeof(cpp_bool))
    if not is_static_ptr or not is_incomp_ptr:
        free(is_static_ptr)
        free(is_incomp_ptr)
        raise MemoryError("Failed to allocate the per-layer assumption arrays.")

    cdef shared_ptr[c_BaseWorld] world_sptr
    cdef size_t layer_i
    try:
        for layer_i in range(num_layers):
            layer_type_vec[layer_i] = 0 if layer_is_solid[layer_i] else 1
            is_static_ptr[layer_i]  = <cpp_bool>bool(layer_is_static[layer_i])
            is_incomp_ptr[layer_i]  = <cpp_bool>bool(layer_is_incompressible[layer_i])
        world_sptr = c_build_world_from_layered_profile(
            &radius[0],
            &density[0],
            &shear_modulus[0],
            &bulk_modulus[0],
            num_slices,
            &upper_radius_bylayer[0],
            layer_type_vec.data(),
            is_static_ptr,
            is_incomp_ptr,
            num_layers,
            planet_bulk_density,
            name.encode('utf-8'))
    finally:
        free(is_static_ptr)
        free(is_incomp_ptr)

    return BaseWorld._wrap(world_sptr)


# The world classes a binary world record is written by (c_binary_class_name spells them so).
WORLD_CLASS_NAMES = ("BaseWorld", "TerrestrialWorld", "GasGiantWorld", "StarWorld")

# The Python-side configurations a world carries beside its C++ state, which a copy or a pickle keeps.
WORLD_CONFIG_ATTRIBUTES = ("source_config", "portable_config", "built_config")


def world_from_bytes(object world_class, bytes record, dict configs=None, cpp_bool force=False):
    """A world rebuilt from one binary world record held in memory, as :meth:`BaseWorld.copy`, unpickling, and
    :func:`~TidalPy.Structures.load_world` make one.

    Parameters
    ----------
    world_class : type
        ``BaseWorld`` or a subclass: the class of the world that wrote the record.
    record : bytes
        The world's binary record, the bytes ``save_binary`` writes to a file.
    configs : dict, optional
        ``source_config``, ``portable_config``, and ``built_config`` for the new world, each copied; a key left out
        leaves that configuration ``None``.
    force : bool, optional
        Read the record even on a schema version mismatch. Default False.

    Returns
    -------
    BaseWorld
        A world of ``world_class`` with the record's state (layers, materials, models, tide model and settings,
        pinned solver settings, prescribed heating), unsolved and outside any system.

    Raises
    ------
    TypeError
        ``world_class`` is not a world class, or the record was written by another world class.
    IOError
        The record is not a world's, or it is corrupt.
    """
    if not (isinstance(world_class, type) and issubclass(world_class, BaseWorld)):
        raise TypeError(f"TidalPy: world_from_bytes builds a BaseWorld or a subclass of it, not {world_class!r}.")
    cdef shared_ptr[c_BaseWorld] world_sptr
    try:
        world_sptr = c_world_from_binary_bytes(record, force)
    except RuntimeError as exc:
        raise IOError(str(exc)) from exc
    cdef str record_class = c_binary_class_name(world_sptr.get().get_binary_class_id()).decode("utf-8")
    cdef str expected_class = next(cls.__name__ for cls in world_class.__mro__ if cls.__name__ in WORLD_CLASS_NAMES)
    if record_class != expected_class:
        raise TypeError(
            f"TidalPy: the binary record is a {record_class}, not a {expected_class}; rebuild it as a {record_class}.")
    cdef BaseWorld world = world_class.__new__(world_class)
    world._bind(world_sptr)
    cy_set_world_configs(world, configs)
    return world


# The configurations a copy or a pickle of a world carries over (WORLD_CONFIG_ATTRIBUTES), by attribute name.
cdef dict cy_world_configs(BaseWorld world):
    return {attribute: getattr(world, attribute) for attribute in WORLD_CONFIG_ATTRIBUTES}


# Gives a world copies of the configurations cy_world_configs took from another (None for any left out).
cdef void cy_set_world_configs(BaseWorld world, dict configs) except *:
    cdef str attribute
    for attribute in WORLD_CONFIG_ATTRIBUTES:
        setattr(world, attribute, copy.deepcopy((configs or {}).get(attribute)))


# A world's binary record, the bytes a copy or a pickle rebuilds it from; RuntimeError, naming the class, for a wrapper
# that holds no world.
cdef bytes cy_world_record(BaseWorld world):
    if world._world_ptr.get() == NULL:
        raise RuntimeError(f"TidalPy: this {type(world).__name__} holds no world to copy.")
    return world._world_ptr.get().write_binary_bytes()


# One 3D grid axis as a contiguous 1-D float64 array (a scalar becomes one value), with its data pointer and
# length written out; the caller keeps the array alive for as long as the pointer is used.
cdef cnp.ndarray cy_grid_axis(object values, const double** data_out, size_t* num_out):
    cdef cnp.ndarray axis = np.ascontiguousarray(np.atleast_1d(values), dtype=np.float64).ravel()
    num_out[0] = <size_t>axis.shape[0]
    data_out[0] = <const double*>cnp.PyArray_DATA(axis) if num_out[0] > 0 else NULL
    return axis


# The four axes of an instantaneous 3D grid, each holding at least one value, as arrays and as `axes`.
cdef tuple cy_grid_axes(object radii, object colatitudes, object longitudes, object times, c_Grid3DAxes* axes):
    cdef cnp.ndarray radii_arr = cy_grid_axis(radii, &axes.radii, &axes.num_radii)
    cdef cnp.ndarray colat_arr = cy_grid_axis(colatitudes, &axes.colatitudes, &axes.num_colatitudes)
    cdef cnp.ndarray lon_arr   = cy_grid_axis(longitudes, &axes.longitudes, &axes.num_longitudes)
    cdef cnp.ndarray time_arr  = cy_grid_axis(times, &axes.times, &axes.num_times)
    if axes.num_radii == 0 or axes.num_colatitudes == 0 or axes.num_longitudes == 0 or axes.num_times == 0:
        raise ValueError("radii, colatitudes, longitudes, and times must each hold at least one value")
    return radii_arr, colat_arr, lon_arr, time_arr


cdef int cy_require_arguments(str method, tuple names, tuple values) except -1:
    """Raise TypeError naming the arguments of ``names`` whose ``values`` are None: the point or grid arguments of the
    3D tide methods, which follow the optional orbital state and so take a None default."""
    cdef list missing = [names[argument_i] for argument_i in range(len(names)) if values[argument_i] is None]
    if missing:
        raise TypeError(f"TidalPy: {method} needs {', '.join(missing)}.")
    return 0


cdef int cy_check_num_threads(int num_threads) except -1:
    """Raise ValueError unless ``num_threads`` is 0 (automatic) or a positive thread count."""
    if num_threads < 0:
        raise ValueError(f"num_threads must be 0 (automatic) or a positive thread count; got {num_threads}")
    return 0

# Integration-method names are resolved at the Cython boundary through the shared tables in TidalPy.constants.
cdef str cy_integration_method_name(ODEMethod method):
    """The configuration spelling of a CyRK integration method."""
    try:
        return ODE_METHOD_NAMES[<int>method]
    except KeyError:
        raise ValueError(f"Unsupported integration method code: {<int>method}.") from None


cdef ODEMethod cy_resolve_integration_method(str integration_method) except *:
    return <ODEMethod><int>ode_method_from_name(integration_method)


cdef void cy_apply_love_solve_overrides(
        c_LoveSolveConfig* cfg,
        object starting_method,
        object degree1_frame,
        object nondimensionalize,
        object start_radius_tol,
        object integration_method,
        object rtol,
        object atol,
        object scale_rtols,
        object max_num_steps,
        object expected_size,
        object max_ram_MB) except *:
    """Overwrite the solver settings of a Love-solve config with the arguments that are not ``None``.

    The config already carries the ``[radial_solver]`` defaults of the TidalPy configuration, so a ``None``
    leaves that default in place.
    """
    if starting_method is not None:
        cfg.starting_method = <int>starting_method_from_name(starting_method)
    if degree1_frame is not None:
        cfg.degree1_frame = <int>degree1_frame_from_name(degree1_frame)
    if nondimensionalize is not None:
        cfg.nondimensionalize = <cpp_bool>bool(nondimensionalize)
    if start_radius_tol is not None:
        cfg.start_radius_tol = <double>start_radius_tol
    if integration_method is not None:
        cfg.integration_method = cy_resolve_integration_method(integration_method)
    if rtol is not None:
        cfg.rtol = <double>rtol
    if atol is not None:
        cfg.atol = <double>atol
    if scale_rtols is not None:
        cfg.scale_rtols = <cpp_bool>bool(scale_rtols)
    if max_num_steps is not None:
        cfg.max_num_steps = <size_t>int(max_num_steps)
    if expected_size is not None:
        cfg.expected_size = <size_t>int(expected_size)
    if max_ram_MB is not None:
        cfg.max_ram_MB = <size_t>int(max_ram_MB)


cdef int cy_resolve_solve_for(str solve_for) except? -999:
    # Map a surface-boundary-condition name to the radial solver's bc_model integer
    # (same names as the standalone radial_solver's solve_for entries).
    cdef str name_lower = solve_for.lower()
    if name_lower == 'tidal':
        return 1
    elif name_lower == 'loading':
        return 2
    elif name_lower == 'free':
        return 0
    raise ValueError(
        f"Unsupported solve_for: {solve_for}. Supported: 'tidal', 'loading', 'free'.")


cdef void cy_set_solve_for(c_LoveSolveConfig* cfg, solve_for) except *:
    # Accept one boundary-condition name or a sequence of them, and write them onto the config in order. A
    # sequence produces one block of radial functions each, from a single integration.
    cdef list names
    cdef vector[int] models
    if isinstance(solve_for, str):
        names = [solve_for]
    else:
        names = list(solve_for)
    if len(names) == 0:
        raise ValueError("solve_for must name at least one surface boundary condition.")
    if len(names) > C_MAX_NUM_YTYPES:
        raise ValueError(
            f"solve_for accepts at most {C_MAX_NUM_YTYPES} boundary conditions; {len(names)} were given.")
    if len(set(name.lower() for name in names)) != len(names):
        raise ValueError(f"solve_for entries must be distinct; got {tuple(names)}.")
    for name in names:
        models.push_back(cy_resolve_solve_for(name))
    cfg.set_bc_models(models.data(), models.size())


# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


# The world's profile read for cy_eos_field and cy_eos_fields (layers/base.pyx): one C++ call per input that holds
# the world's call lock throughout.
cdef void cy_world_eos_fields(
        const void* owner,
        const size_t* field_indices,
        size_t num_fields,
        const double* radii,
        size_t num_radii,
        double* values_out) noexcept nogil:
    (<const c_BaseWorld*>owner).get_eos_fields(field_indices, num_fields, radii, num_radii, values_out)


# A pinned solver setting's kinds (c_SolverSettingKind): the order of the C++ enum.
cdef int C_SOLVER_SETTING_METHOD = 0
cdef int C_SOLVER_SETTING_REAL   = 1
cdef int C_SOLVER_SETTING_COUNT  = 2
cdef int C_SOLVER_SETTING_FLAG   = 3
cdef int C_SOLVER_SETTING_STARTING_METHOD = 4
cdef int C_SOLVER_SETTING_DEGREE1_FRAME   = 5


cdef double cy_solver_setting_to_double(int kind, object value) except? -1.0:
    """A validated [eos_solver] / [radial_solver] value as the double the world's override table keeps."""
    if kind == C_SOLVER_SETTING_METHOD:
        return <double><int>cy_resolve_integration_method(str(value))
    if kind == C_SOLVER_SETTING_COUNT:
        return <double>int(value)
    if kind == C_SOLVER_SETTING_FLAG:
        return 1.0 if value else 0.0
    if kind == C_SOLVER_SETTING_STARTING_METHOD:
        return <double>starting_method_from_name(value)
    if kind == C_SOLVER_SETTING_DEGREE1_FRAME:
        return <double>degree1_frame_from_name(value)
    return <double>value


cdef object cy_solver_setting_from_double(int kind, double value):
    """A pinned setting back under its configuration type: a method name, an int, a bool, or a float."""
    if kind == C_SOLVER_SETTING_METHOD:
        return cy_integration_method_name(<ODEMethod><int>value)
    if kind == C_SOLVER_SETTING_COUNT:
        return int(value)
    if kind == C_SOLVER_SETTING_FLAG:
        return value != 0.0
    if kind == C_SOLVER_SETTING_STARTING_METHOD:
        return STARTING_METHOD_NAMES[<int>value]
    if kind == C_SOLVER_SETTING_DEGREE1_FRAME:
        return DEGREE1_FRAME_NAMES[<int>value]
    return value


cdef class BaseWorld(StructureBase):
    """A world: an ordered (inner-to-outer) stack of layers, which may be empty, with its identity, orbital and
    thermal scalars, bulk geometry, spin model, and tide model.

    Parameters
    ----------
    name : str
        World name.
    radius : float
        World radius [m].
    mass : float
        World mass [kg].
    world_type : str, optional
        Free-form type label (e.g. ``"terrestrial"``). Default ``"world"``.
    albedo : float, optional
        Bond albedo [dimensionless].
    emissivity : float, optional
        Surface emissivity [dimensionless].
    obliquity : float, optional
        Axial obliquity [rad].
    spin_frequency : float, optional
        Rotation rate [rad/s].

    Each property not given takes the ``[worlds]`` default of the TidalPy configuration for the world type
    (``[worlds.<type>]`` winning), as :func:`~TidalPy.Structures.build_world` does, and so does the spin model's
    moment-of-inertia factor. A world built here has no tide model until :meth:`set_tide_model`.

    Notes
    -----
    Add layers inner-to-outer with :meth:`add_layer`; each layer's inner radius must match the previous layer's
    outer radius (the innermost starts at 0). A world with no layers still dissipates tidally through the analytic
    tide models (``cpl``, ``ctl``, ``ctl_q``).

    Assumptions
    -----------
    - Spherically symmetric world.
    """

    def __cinit__(self, *args, **kwargs):
        pass  # the shared_ptr starts empty; __init__ or _wrap binds it

    def __init__(
            self,
            str    name,
            double radius,
            double mass,
            str    world_type = "world",
            albedo            = None,
            emissivity        = None,
            obliquity         = None,
            spin_frequency    = None):
        cdef c_WorldConfig config
        cdef dict defaults = cy_fill_world_config(
            &config,
            name,
            radius,
            mass,
            world_type,
            albedo,
            emissivity,
            obliquity,
            spin_frequency)
        # _world_ptr is a shared_ptr (a System co-owns the world), so build it with make_shared.
        self._bind(make_shared[c_BaseWorld](config))
        cy_set_default_spin(self, defaults)

    def __dealloc__(self):
        self._world_ptr.reset()
        self._ptr = NULL

    cdef void _bind(self, shared_ptr[c_BaseWorld] ptr):
        self._world_ptr = ptr
        self._ptr = <c_TidalPyBaseClass*>ptr.get()
        self._layer_views = None
        self._layer_view_by_name = None

    @staticmethod
    cdef BaseWorld _wrap(shared_ptr[c_BaseWorld] ptr):
        """Wrap an already-constructed C++ world (no new C++ object is built)."""
        cdef BaseWorld world = BaseWorld.__new__(BaseWorld)
        world._bind(ptr)
        return world

    def add_layer(self, Layer layer not None):
        """Add a layer to the world (inner to outer).

        Ownership of the C++ layer and its attached physics models moves to the world. ``layer`` stays usable: it
        becomes a non-owning view of the layer the world now holds, the same kind ``world.<layer name>`` returns.

        Parameters
        ----------
        layer : Layer
            A layer constructed without a ``layer_index`` takes its place in the stack (the number of layers already
            added), and one without a ``radius_inner`` starts at the current outermost radius (0 for the first layer).
            A given index or inner radius must match those.

        Raises
        ------
        ValueError
            If ``layer`` has already been added or moved, is a view of another world's layer, has another index than
            its place in the stack, does not continue the existing stack, reaches past the world radius, or shares
            its name with a layer already added.
        """
        if layer._layer_ptr.get() == NULL:
            raise ValueError(
                "This layer holds no C++ object (already added to a world or moved).")
        if layer._is_view:
            raise ValueError(
                "This layer is a non-owning view into a world-owned layer and cannot be "
                "added to another world. Construct a new layer instead.")
        cdef c_Layer* added_layer_ptr = layer._layer_ptr.get()
        cdef size_t position = self._world_ptr.get().get_num_layers()
        cdef double stack_top = 0.0
        if position > 0:
            stack_top = self._world_ptr.get().get_layer(position - 1).get_radius_outer()
        cdef int index_before = added_layer_ptr.get_layer_index()
        cdef double inner_before = added_layer_ptr.get_radius_inner()
        if (not layer.p_inner_given) and (added_layer_ptr.get_radius_outer() < stack_top):
            raise ValueError(
                f"TidalPy: layer '{layer.name}' ends at {added_layer_ptr.get_radius_outer()} m, below the top of the "
                f"stack of world '{self.name}' at {stack_top} m.")
        # What the constructor was not given comes from the layer's place in the stack; what it was given is checked.
        if not layer.p_index_given:
            added_layer_ptr.set_layer_index(<int>position)
        if not layer.p_inner_given:
            added_layer_ptr.set_radii(stack_top, added_layer_ptr.get_radius_outer())
        # Validate before transferring ownership so a rejected layer stays usable, as it was constructed (the C++
        # add_layer would otherwise consume it on throw).
        cdef string rejection = self._world_ptr.get().layer_rejection_reason(deref(added_layer_ptr))
        if rejection.size() > 0:
            added_layer_ptr.set_layer_index(index_before)
            added_layer_ptr.set_radii(inner_before, added_layer_ptr.get_radius_outer())
            raise ValueError(rejection.decode('utf-8'))
        self._world_ptr.get().add_layer(move(layer._layer_ptr))
        layer._init_view(added_layer_ptr, self)
        self._track_view(layer)
        # The layer set changed; drop the cached views so they rebuild on next access.
        self._layer_views = None
        self._layer_view_by_name = None

    def copy(self):
        """An independent world of the same class with the same state.

        The copy goes through the binary record (:meth:`save_binary` without the file), so it has this world's
        layers, materials, physics models, tide model and tide settings, pinned solver settings, spin model, and
        prescribed heating, and its own copies of :attr:`source_config`, :attr:`portable_config`, and
        :attr:`built_config`. Like a loaded world it is unsolved (run ``solve_eos`` and ``calc_tides`` again), holds
        no tidal heat source, and belongs to no system. ``copy.copy``, ``copy.deepcopy``, and ``pickle`` use it, so a
        world can be handed to a process pool.

        Returns
        -------
        BaseWorld
            A world of this world's class.
        """
        return world_from_bytes(type(self), cy_world_record(self), cy_world_configs(self))

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def __reduce__(self):
        # Pickled as its binary record and its configurations, rebuilt by world_from_bytes (see copy).
        return (world_from_bytes, (type(self), cy_world_record(self), cy_world_configs(self)))

    def __repr__(self):
        """One line: the class, the name, the radius [km], the mass [kg], and the number of layers."""
        if self._world_ptr.get() == NULL:
            return f"{type(self).__name__}(no world)"
        return (f"{type(self).__name__}({self.name!r}, radius_km={self.radius / 1.0e3:.6g}, "
                f"mass_kg={self.mass:.4e}, num_layers={self.num_layers})")

    def summary(self) -> str:
        """A multi-line description of the world and its layers.

        The first line is the world (class, name, radius, mass, and whether its EOS is solved); then one row per
        layer, inner to outer: its name, radii [km], state, material phases, density [kg m-3] (the range of the solved
        profile, ``-`` before a solve), temperature [K], and shear rheology.

        Returns
        -------
        str
        """
        cdef list lines = [
            f"{type(self).__name__} '{self.name}': radius {self.radius / 1.0e3:.6g} km, mass {self.mass:.4e} kg, "
            f"{self.num_layers} layers, EOS {'solved' if self.eos_solved else 'not solved'}"]
        cdef list rows = [("layer", "radii [km]", "state", "material", "density [kg m-3]", "temperature [K]",
                           "rheology")]
        cdef list phases
        for layer in self._ensure_layer_views():
            material = layer.material
            phases = [] if material is None else [
                slot for slot in ("solid", "liquid") if getattr(material, slot) is not None]
            density_text = "-"
            if self.eos_solved:
                sample_radii = np.linspace(layer.radius_inner, layer.radius_outer, SUMMARY_DENSITY_SAMPLES)
                density = layer.get_density(sample_radii)
                if np.any(np.isfinite(density)):
                    density_text = f"{np.nanmin(density):.0f} to {np.nanmax(density):.0f}"
            rheology = layer.shear_rheology
            rows.append((
                layer.name,
                f"{layer.radius_inner / 1.0e3:.6g} to {layer.radius_outer / 1.0e3:.6g}",
                layer.state,
                "+".join(phases) if phases else "-",
                density_text,
                f"{layer.temperature:.6g}",
                "-" if rheology is None else rheology.model_name))
        cdef list widths = [max(len(row[column]) for row in rows) for column in range(len(rows[0]))]
        for row in rows:
            lines.append("  " + "  ".join(text.ljust(width) for text, width in zip(row, widths)).rstrip())
        return "\n".join(lines)

    def load_binary(self, path, cpp_bool force=False):
        """Load this world's state, its layers included, from a TidalPy binary file.

        The load replaces the world's layers, so layer views taken from this world before it (``world.<name>``,
        ``get_layer``, iteration) no longer refer to a layer of this world and raise if used; take new ones.
        Nothing solved survives the load: run ``solve_eos`` again. The configurations the world was built from
        (:attr:`source_config`, :attr:`portable_config`, and :attr:`built_config`) describe it before the load, so a
        load that succeeds clears them and :meth:`save_to_toml` then writes the loaded world's own
        :meth:`get_config_dict`. A load that fails leaves the world as it was, its layers, views, solved state, and
        configurations included: the file is read into a new world first and reaches this one only once that read
        succeeded.

        Parameters
        ----------
        path : str or os.PathLike
            Source file path.
        force : bool, optional
            Attempt the load even on a schema version mismatch.

        Raises
        ------
        FileNotFoundError
            ``path`` does not exist.
        IOError
            The file holds a record of another class (the message names both, as in "a StarWorld file, not a
            TerrestrialWorld one"), has an incompatible schema version, or is corrupt.
        """
        cdef object view_ref
        cdef Layer view
        StructureBase.load_binary(self, path, force)
        self.source_config   = None
        self.portable_config = None
        self.built_config    = None
        # Every view handed out points at a C++ layer the load replaced: detach them, so a held view raises instead
        # of reading freed memory. The load holds the GIL, so no view is used between the load and this.
        if self._issued_views is not None:
            for view_ref in self._issued_views:
                view = view_ref()
                if view is not None:
                    view._detach()
            self._issued_views = []
        self._layer_views = None
        self._layer_view_by_name = None

    @property
    def radius(self) -> float:
        """World radius [m]."""
        return self._world_ptr.get().get_radius()

    @property
    def mass(self) -> float:
        """World mass [kg]."""
        return self._world_ptr.get().get_mass()

    @property
    def name(self) -> str:
        """World name."""
        return self._world_ptr.get().get_name().decode("utf-8")

    @name.setter
    def name(self, str value):
        self._world_ptr.get().set_name(value.encode("utf-8"))

    @property
    def world_type(self) -> str:
        """World type label."""
        return self._world_ptr.get().get_world_type().decode("utf-8")

    @property
    def albedo(self) -> float:
        """Bond albedo [dimensionless]."""
        return self._world_ptr.get().get_albedo()

    @property
    def emissivity(self) -> float:
        """Surface emissivity [dimensionless]."""
        return self._world_ptr.get().get_emissivity()

    @property
    def obliquity(self) -> float:
        """Axial obliquity [rad]."""
        return self._world_ptr.get().get_obliquity()

    @property
    def spin_frequency(self) -> float:
        """Rotation rate [rad/s]."""
        return self._world_ptr.get().get_spin_frequency()

    def calc_surface_gravity(self) -> float:
        """Surface gravitational acceleration [m/s^2] = G·M/R²."""
        return self._world_ptr.get().calc_surface_gravity()

    def calc_escape_velocity(self) -> float:
        """Escape velocity [m/s] = sqrt(2·G·M/R)."""
        return self._world_ptr.get().calc_escape_velocity()

    def calc_mean_density(self) -> float:
        """Mean density [kg/m^3] = M / V_sphere(R)."""
        return self._world_ptr.get().calc_mean_density()

    def calc_equilibrium_temperature(self, insolation_flux):
        """Radiative-equilibrium temperature [K] for a given insolation flux.

        T_eq = [ (1 − A)·F / (4·ε·σ) ]^(1/4)

        Parameters
        ----------
        insolation_flux : float or np.ndarray
            Incident stellar flux [W/m^2].

        Returns
        -------
        float or np.ndarray
            Equilibrium temperature [K], shaped like ``insolation_flux``; 0.0 for non-positive flux.

        Assumptions
        -----------
        - Fast-rotator, uniform-temperature surface.
        """
        cdef c_BaseWorld* world_ptr = self._world_ptr.get()
        if np.ndim(insolation_flux) == 0:
            return world_ptr.calc_equilibrium_temperature(<double>insolation_flux)
        cdef cnp.ndarray flux_array = np.ascontiguousarray(insolation_flux, dtype=np.float64)
        cdef cnp.ndarray out = np.empty_like(flux_array)
        cdef double[::1] flux_view = flux_array.reshape(-1)
        cdef double[::1] out_view = out.reshape(-1)
        cdef Py_ssize_t flux_i
        for flux_i in range(flux_view.shape[0]):
            out_view[flux_i] = world_ptr.calc_equilibrium_temperature(flux_view[flux_i])
        return out

    def set_spin_frequency(self, double freq):
        """Set the rotation rate [rad/s].

        The stored rate is what a ``System`` reads when it builds the world's tidal state. ``calc_tides`` takes
        its spin rate as an argument, so a tidal result already solved keeps describing the state it was given.
        """
        self._world_ptr.get().set_spin_frequency(freq)

    def set_obliquity(self, double obliq):
        """Set the axial obliquity [rad].

        The stored obliquity is what a ``System`` reads when it builds the world's tidal state. ``calc_tides``
        takes its obliquity as an argument, so a tidal result already solved keeps describing the state it was
        given.
        """
        self._world_ptr.get().set_obliquity(obliq)

    # Global (1D) tidal dissipation (analytic path; common to all world types)
    def set_tide_model(self, tide):
        """Attach a global tide dissipation model.

        A model is shared, not copied (models are not changed in place), so ``tide`` stays usable and can be
        attached to other worlds.

        The analytic models (``cpl``/``ctl``/``ctl_q``) work on any world, layers or not; the ``rheology`` model
        needs layers and a solved EOS.

        Parameters
        ----------
        tide : TideBase, str, dict, or None
            A tide model; a model name (``"fixed_q"``, ``"cpl"``, ...), built by ``make_tide`` with its defaults; or
            a table like a world file's ``[tides]`` table: the model name under ``model`` or ``global_tidal_model``,
            the model's parameters (``fixed_k``, ``fixed_q``, ``fixed_dt_s``; another model's lists are ignored, as
            the world builder ignores them), and optionally any :meth:`set_tide_config` setting, which is applied
            too. ``None`` detaches the model.

        Raises
        ------
        TypeError
            ``tide`` is none of these.
        ValueError
            A table without a model name or with a key that is neither a tide model parameter nor a tide setting,
            an unknown model name, or a bad parameter or setting. Nothing changes then.
        """
        cdef dict settings = {}
        cdef dict table
        cdef str model_name
        cdef frozenset model_keys
        cdef list unknown
        if tide is None or isinstance(tide, TideBase):
            model = tide
        elif isinstance(tide, str):
            model = make_tide(tide)
        elif isinstance(tide, dict):
            table = dict(tide)
            names = [key for key in TIDE_MODEL_NAME_KEYS if key in table]
            if len(names) != 1:
                raise ValueError(
                    f"TidalPy: a tide model table names its model under one of {list(TIDE_MODEL_NAME_KEYS)}; got the "
                    f"keys {sorted(table)}.")
            model_name = str(table.pop(names[0]))
            model_keys = tide_config_keys(model_name)
            parameters = {key: table.pop(key) for key in list(table) if key in model_keys}
            settings = {key: table.pop(key) for key in list(table) if key in TIDE_SETTING_KEYS}
            unknown = [key for key in table if key not in TIDE_CONFIG_KEYS]
            if unknown:
                accepted = sorted(model_keys | TIDE_SETTING_KEYS)
                raise ValueError(
                    f"TidalPy: the tide model table has no key '{unknown[0]}'{did_you_mean(unknown[0], accepted)}; "
                    f"it takes {', '.join(accepted)}, plus the model name.")
            model = make_tide(model_name, parameters)
        else:
            raise TypeError(
                f"TidalPy: a world's tide model is a TideBase, a model name, a tide model table, or None, not "
                f"{type(tide).__name__}.")
        cdef shared_ptr[c_PhysicsBase] handle
        if model is not None:
            handle = (<TideBase>model)._model_sptr
            if handle.get() == NULL:
                raise ValueError("This tide model holds no C++ object.")
        # The settings are checked and applied before the model is attached, so a bad one leaves the world as it was.
        if settings:
            self.set_tide_config(**settings)
        self._world_ptr.get().set_tide_model_handle(handle)

    @property
    def tide_model(self):
        """The attached tide model (a ``TideBase``), or None. Setting it calls :meth:`set_tide_model`."""
        return cy_wrap_model(self._world_ptr.get().share_tide_model())

    @tide_model.setter
    def tide_model(self, value):
        self.set_tide_model(value)

    @property
    def tide_model_set(self) -> bool:
        """Whether a tide dissipation model has been attached."""
        return self._world_ptr.get().get_tide_model_set()

    def set_tide_config(
            self,
            min_degree_l=None,
            max_degree_l=None,
            eccentricity_truncation=None,
            obliquity_truncation=None,
            layer_tidal_heating=None,
            eccentricity_exact_tolerance=None,
            love_method=None,
            love_fixed_q=None,
            love_fixed_dt=None,
            *,
            eccentricity_trunc_lvl=None,
            obliquity_trunc_lvl=None,
            love_fixed_dt_s=None):
        """Change the stored ``[tides]`` truncation/degree configuration and the world's Love-number method.

        Only the arguments given change; every other setting keeps its current value (see
        :meth:`get_tide_config`), so a call can adjust one setting without resetting the rest. Every key
        :meth:`get_tide_config` returns is accepted, so ``world.set_tide_config(**world.get_tide_config())`` restores
        a configuration; ``eccentricity_trunc_lvl``, ``obliquity_trunc_lvl``, and ``love_fixed_dt_s`` (the world file's
        key names) are the same settings as ``eccentricity_truncation``, ``obliquity_truncation``, and
        ``love_fixed_dt``, and a call gives each setting under one name only.

        Parameters
        ----------
        min_degree_l, max_degree_l : int, optional
            Tidal harmonic degree range (2..10).
        eccentricity_truncation, eccentricity_trunc_lvl : int or str, optional
            Eccentricity-function truncation level N: every product of two eccentricity functions is kept through
            e^N. Tabulated levels: ``TidalPy.Tides.eccentricity.ECCENTRICITY_TRUNCATIONS`` (2, 4, 6, 8, 10, 20,
            50), or ``"exact"`` for the functions from the exact orbit (any e < 1).
        eccentricity_exact_tolerance : float, optional
            Heating tail tolerance in (0, 1) that sets the mode range of ``"exact"`` (ignored by the levels).
        obliquity_truncation, obliquity_trunc_lvl : int or str, optional
            Obliquity truncation level N: every product of two obliquity functions is kept through I^N. Tabulated
            levels: ``TidalPy.Tides.obliquity.OBLIQUITY_TRUNCATIONS`` (0 or ``"off"``, 2, 4), or ``"gen"`` for
            the general functions (any obliquity).
        layer_tidal_heating : bool, optional
            Whether ``calc_tides`` also resolves each layer's heating when the Love numbers come from the radial
            solver (a volume integral of the radial solution that costs about as much as the global solve again);
            the other paths share out the heating at no extra cost. Default ``True``.
        love_method : str, optional
            How the world obtains Love numbers when its tide model asks for them (and the default for
            ``solve_love_numbers``): ``'radial_solver'`` (``'shooting'``, ``'rs'``), ``'propagation_matrix'``
            (``'prop_matrix'``, ``'pm'``, ``'prop'``), ``'homogeneous'`` (``'homogen'``), ``'cpl'``, ``'ctl'``,
            or ``'laterally_inhomogeneous'`` (``'3d'``, ``'lat_inhom'``; reserved: raises ``NotImplementedError``).
        love_fixed_q, love_fixed_dt (love_fixed_dt_s) : float, optional
            Quality factor for the ``'cpl'`` method and time lag [s] for the ``'ctl'`` method. A NaN clears the
            value, after which the attached tide model's per-degree fixed Q / time lag is used.

        Raises
        ------
        ValueError
            A setting given under both of its names, or a value out of range.
        """
        cdef tuple spelling_pairs = (
            ("eccentricity_truncation", eccentricity_truncation, "eccentricity_trunc_lvl", eccentricity_trunc_lvl),
            ("obliquity_truncation", obliquity_truncation, "obliquity_trunc_lvl", obliquity_trunc_lvl),
            ("love_fixed_dt", love_fixed_dt, "love_fixed_dt_s", love_fixed_dt_s))
        cdef Py_ssize_t pair_i
        for pair_i in range(len(spelling_pairs)):
            if (spelling_pairs[pair_i][1] is not None) and (spelling_pairs[pair_i][3] is not None):
                raise ValueError(
                    f"TidalPy: set_tide_config got both '{spelling_pairs[pair_i][0]}' and "
                    f"'{spelling_pairs[pair_i][2]}', two names for one setting; give one.")
        if eccentricity_trunc_lvl is not None:
            eccentricity_truncation = eccentricity_trunc_lvl
        if obliquity_trunc_lvl is not None:
            obliquity_truncation = obliquity_trunc_lvl
        if love_fixed_dt_s is not None:
            love_fixed_dt = love_fixed_dt_s
        if eccentricity_truncation is not None:
            eccentricity_truncation = validate_eccentricity_truncation(eccentricity_truncation)
        if eccentricity_exact_tolerance is not None:
            eccentricity_exact_tolerance = validate_eccentricity_exact_tolerance(eccentricity_exact_tolerance)
        if obliquity_truncation is not None:
            obliquity_truncation = validate_obliquity_truncation(obliquity_truncation)
        # Start from the stored configuration so an omitted argument leaves its setting unchanged.
        cdef c_TideConfig cfg = self._world_ptr.get().get_tide_config()
        if min_degree_l is not None:
            cfg.min_degree_l = <int>min_degree_l
        if max_degree_l is not None:
            cfg.max_degree_l = <int>max_degree_l
        if eccentricity_truncation is not None:
            cfg.eccentricity_truncation = <int>eccentricity_truncation
        if eccentricity_exact_tolerance is not None:
            cfg.eccentricity_exact_tolerance = <double>eccentricity_exact_tolerance
        if obliquity_truncation is not None:
            cfg.obliquity_truncation = <int>obliquity_truncation
        if layer_tidal_heating is not None:
            cfg.layer_tidal_heating = <cpp_bool>bool(layer_tidal_heating)
        if love_method is not None:
            cfg.love_method = cy_parse_love_method(str(love_method))
        if love_fixed_q is not None:
            cfg.love_fixed_q = <double>love_fixed_q
        if love_fixed_dt is not None:
            cfg.love_fixed_dt = <double>love_fixed_dt
        self._world_ptr.get().set_tide_config(cfg)

    @contextlib.contextmanager
    def temporary_tide_config(self, **settings):
        """A context manager that applies tide settings for the duration of a ``with`` block.

        On entry it calls :meth:`set_tide_config` with ``settings``; on exit, also when the block raises, it restores
        the configuration :meth:`get_tide_config` reported before. A new tide configuration clears the last tidal
        result, so read the result inside the block (``calc_tides`` returns it).

        Parameters
        ----------
        **settings
            Any :meth:`set_tide_config` keyword.

        Yields
        ------
        BaseWorld
            This world.

        Examples
        --------
        >>> with world.temporary_tide_config(eccentricity_truncation=20):
        ...     result = world.calc_tides()
        """
        cdef dict previous = self.get_tide_config()
        # get_tide_config leaves an unset Love-method Q or time lag out; NaN clears whatever the block set.
        previous.setdefault("love_fixed_q", math.nan)
        previous.setdefault("love_fixed_dt_s", math.nan)
        self.set_tide_config(**settings)
        try:
            yield self
        finally:
            self.set_tide_config(**previous)

    def get_tide_state(self):
        """The orbital state this world's tides are raised in, as the system it belongs to sees it.

        Orbital state never lives on a world: the system a world was added to supplies it, from the world's
        orbit about its tidal host, that host's mass, and the world's own spin and obliquity.

        Returns
        -------
        dict or None
            ``orbital_frequency`` [rad s-1], ``spin_frequency`` [rad s-1], ``eccentricity``, ``obliquity``
            [rad], ``semi_major_axis`` [m], and ``host_mass`` [kg], the arguments of :meth:`calc_tides` in
            its order. ``None`` for a world outside a system, with no tidal host, or with no usable orbit.
        """
        cdef c_TideSolveConfig state
        if not self._world_ptr.get().get_tide_state(state):
            return None
        return {
            "orbital_frequency": state.orbital_frequency,
            "spin_frequency":    state.spin_frequency,
            "eccentricity":      state.eccentricity,
            "obliquity":         state.obliquity,
            "semi_major_axis":   state.semi_major_axis,
            "host_mass":         state.host_mass,
        }

    @property
    def tides_solved(self) -> bool:
        """Whether a successful :meth:`calc_tides` result is held.

        A new tide model or tide config clears it, and so does ``solve_eos``: the result describes the structure it
        was solved with.
        """
        return self._world_ptr.get().get_tides_solved()

    def get_tidal_heating(self) -> float:
        """Total global tidal heating [W] (NaN if unsolved)."""
        return self._world_ptr.get().get_tidal_heating()

    def get_tidal_heat_flux(self) -> float:
        """The tidal heating spread over the surface [W m-2]: :meth:`get_tidal_heating` / (4 pi R^2); NaN before a
        :meth:`calc_tides`."""
        return self._world_ptr.get().get_tidal_heat_flux()

    def get_tidal_potential_derivatives(self) -> dict:
        """The three orbital potential derivatives of the last :meth:`calc_tides` [J kg-1 rad-1]; NaN before one.

        Returns
        -------
        dict
            ``dU_dM`` (by the mean anomaly), ``dU_dw`` (by the argument of pericenter), and ``dU_dO`` (by the
            longitude of the ascending node), the names :meth:`System.calc_world_evolution` uses.
        """
        return {
            "dU_dM": self._world_ptr.get().get_tidal_dU_dM(),
            "dU_dw": self._world_ptr.get().get_tidal_dU_dw(),
            "dU_dO": self._world_ptr.get().get_tidal_dU_dO(),
        }

    def get_tidal_dU_dM_minus_dw(self) -> float:
        """The per-mode sum of ``dUdM - dUdw`` from the last tide solve [J kg-1 rad-1]; NaN before one.

        Pass it as ``dU_dM_minus_dw`` to :meth:`OrbitSolver.calc_de_dt`: at small eccentricity the separate sums of
        :meth:`get_tidal_potential_derivatives` nearly cancel in ``de/dt``, and this sum keeps the rate exact.
        """
        return self._world_ptr.get().get_tidal_dU_dM_minus_dw()

    def get_num_tidal_modes(self) -> int:
        """Number of active (nonzero-frequency) tidal modes summed in the last solve."""
        return self._world_ptr.get().get_num_tidal_modes()

    def get_tidal_love_k(self, degree_l: int, m: int, p: int, q: int) -> complex:
        """Complex potential Love number ``k_l`` for the tidal mode ``(l, m, p, q)``.

        Only populated for the rheology tide model. Returns NaN for the analytic models, which carry no
        displacement Love numbers, or for an inactive mode.
        """
        cdef cpp_complex[double] k = self._world_ptr.get().get_tidal_love_k(
            <int>degree_l, <int>m, <int>p, <int>q)
        return PyComplex_FromDoubles(k.real(), k.imag())

    @staticmethod
    def build(source, overrides=None, force=False, base_dir=None):
        """Build a world from a configuration source (the public builder entry point).

        This is a factory: the concrete subclass returned (``BaseWorld``, ``TerrestrialWorld``,
        ``GasGiantWorld``, or ``StarWorld``) follows the configuration's world ``type``, whichever class ``build`` is
        called on. The normalized configuration is retained on :attr:`source_config` so the world can be written back
        to TOML. A path to a binary file (``save_binary``) is loaded instead (:func:`~TidalPy.Structures.load_world`).

        Parameters
        ----------
        source : str, os.PathLike, or dict
            A bundled world name, a path to a ``.toml`` or binary file, or a configuration dict.
        overrides : dict, optional
            Values merged over the configuration before it is built, table by table, so a nested table only needs the
            keys it changes (``{"tides": {"eccentricity_trunc_lvl": 20}}``). Not accepted with a binary file.
            Default None.
        force : bool, optional
            If True, bypass the schema-version compatibility warning. Default False.
        base_dir : str, optional
            The folder a relative ``data_file`` of a configuration dict is relative to (a system file's folder, for a
            world given inline in it). A file source uses its own folder. Default None: the search of
            ``resolve_data_file`` alone.

        Returns
        -------
        BaseWorld
            The constructed world, with ``source_config`` populated.

        Raises
        ------
        ValueError
            The configuration fails validation, or ``overrides`` is given with a binary file.
        """
        # Deferred imports: the builder helpers import the world subclasses, so
        # importing them at module load would be circular.
        from TidalPy.Structures.configs.world_builder import (
            _construct_owned_world,
            _resolve_source,
            load_world)
        from TidalPy.Structures.configs.toml_loader import load_build_config
        from TidalPy.Structures.configs.worldpack import binary_source_class, resolve_data_file

        # `_resolve_source` hands back whatever the caller gave (a path string, a Path, or an already-parsed
        # mapping), so this one stays `object`; the rest have a single concrete type.
        cdef object resolved = _resolve_source(source)
        if binary_source_class(resolved) is not None:
            if overrides is not None:
                raise ValueError(
                    f"TidalPy: '{resolved}' is a binary file, which overrides cannot change; load it and then change "
                    "the world.")
            return load_world(resolved, force=force)
        cdef dict config = load_build_config(resolved, overrides, force)
        cdef object given_data_file = None
        cdef object data_base_dir
        cdef BaseWorld world
        # Resolve a companion data file (e.g. a PREM profile) relative to the world
        # file's directory so the builder can open it directly.
        if "data_file" in config:
            given_data_file = config["data_file"]
            data_base_dir = os.path.dirname(resolved) if isinstance(resolved, str) else base_dir
            config["data_file"] = resolve_data_file(config["data_file"], data_base_dir)
        # load_toml made this config a private copy, so the world keeps it without copying again.
        world = _construct_owned_world(config)
        if given_data_file is not None and world.portable_config is not None:
            # A saved copy names the file as this one did, not the path it resolved to on this machine.
            world.portable_config["data_file"] = given_data_file
        return world

    @property
    def config(self):
        """The normalized configuration dict the world was built from (None if built directly).

        For a world built from a ``data_file`` this is the expanded form, with the profile's layers; the file
        reference as given is kept on :attr:`portable_config`, which :meth:`save_to_toml` writes instead.
        """
        return self.source_config

    def family_world_type(self) -> str:
        """Builder world ``type`` for this class family, used when the stored label is not a builder type."""
        return "layered"

    def get_builder_world_type(self) -> str:
        """World ``type`` as the TOML builder names it.

        The stored ``world_type`` label is used when it is one of ``BUILDER_WORLD_TYPES``; otherwise the
        class family's default applies (``layered``, ``terrestrial``, ``gasgiant``, or ``star``).
        """
        cdef str stored = self.world_type
        if stored in BUILDER_WORLD_TYPES:
            return stored
        return self.family_world_type()

    def get_tide_config(self) -> dict:
        """Return the stored ``[tides]`` degree and truncation settings under the builder's key names.

        Returns
        -------
        dict
            ``min_degree_l``, ``max_degree_l``, ``eccentricity_trunc_lvl`` (an int, or ``"exact"``),
            ``eccentricity_exact_tolerance``, ``obliquity_trunc_lvl`` (an int, or ``"gen"``),
            ``layer_tidal_heating``, ``love_method``, and
            ``love_fixed_q`` / ``love_fixed_dt_s`` when set.
        """
        cdef c_TideConfig cfg = self._world_ptr.get().get_tide_config()
        cdef dict out = {
            "min_degree_l":                  cfg.min_degree_l,
            "max_degree_l":                  cfg.max_degree_l,
            "eccentricity_trunc_lvl":        eccentricity_truncation_name(cfg.eccentricity_truncation),
            "eccentricity_exact_tolerance":  cfg.eccentricity_exact_tolerance,
            "obliquity_trunc_lvl":           obliquity_truncation_name(cfg.obliquity_truncation),
            "layer_tidal_heating":           bool(cfg.layer_tidal_heating),
            "love_method":                   c_love_method_name_int(cfg.love_method).decode('utf-8'),
        }
        if cfg.love_fixed_q == cfg.love_fixed_q:      # not NaN
            out["love_fixed_q"] = cfg.love_fixed_q
        if cfg.love_fixed_dt == cfg.love_fixed_dt:
            out["love_fixed_dt_s"] = cfg.love_fixed_dt
        return out

    def get_save_config(self, destination_dir=None) -> dict:
        """The configuration :meth:`save_to_toml` writes: the world as it is now.

        A world built from a ``data_file`` and unchanged since its build gives :attr:`portable_config`, the file
        reference and the tables that refined its layers, with the reference rewritten to find the same file from
        ``destination_dir``. Any other world, a data-file world changed after its build included, gives
        :meth:`get_config_dict`, validated against the world schema. The layer masses an EOS solve sets and the
        world's name do not count as changes; the name written is the current one.

        Parameters
        ----------
        destination_dir : str or os.PathLike, optional
            The folder the configuration will be saved into. Default None: the data-file reference as given.

        Returns
        -------
        dict
            A world configuration that builds this world.

        Raises
        ------
        ValueError
            The live configuration fails the world schema.
        """
        from TidalPy.Structures.configs.config_writer import relocated_path, state_changed_since_build
        from TidalPy.Structures.configs.toml_loader import validate_world_config
        cdef dict live = self.get_config_dict()
        cdef dict config
        if (self.portable_config is not None) and not state_changed_since_build(self.built_config, live):
            config = copy.deepcopy(self.portable_config)
            config["name"] = live["name"]
            if (destination_dir is not None) and ("data_file" in (self.source_config or {})):
                config["data_file"] = relocated_path(
                    config["data_file"], self.source_config["data_file"], os.fspath(destination_dir))
            return config
        validate_world_config(live)
        return live

    def save_to_toml(self, file_path, overwrite=True):
        """Write this world's configuration, as it is now, to a TOML file.

        Writes :meth:`get_save_config`: every change made after the build is saved (a new obliquity, a layer's
        model, its temperature). A world built from a ``data_file`` and unchanged since writes the file reference
        instead of the expanded profile, rewritten so it finds the same file from the folder saved into. The file
        starts with a comment header naming the TidalPy, SciPy, and CyRK versions that wrote it.

        Parameters
        ----------
        file_path : str or os.PathLike
            Destination ``.toml`` path.
        overwrite : bool, optional
            Overwrite an existing file. Default True.

        Raises
        ------
        ValueError
            The path does not end in ``.toml``, or the live configuration fails the world schema.
        FileExistsError
            The file exists and ``overwrite`` is False.
        """
        from TidalPy.Structures.configs.config_writer import save_world_to_toml
        cdef str path = os.fspath(file_path)
        cdef dict config = self.get_save_config(os.path.dirname(os.path.abspath(path)))
        return save_world_to_toml(config, path, overwrite=overwrite)

    cdef void _track_view(self, Layer view) except *:
        """Remember a view this world handed out (weakly), so a load that replaces the layers can detach it."""
        if self._issued_views is None:
            self._issued_views = []
        # Drop references to views that no longer exist, so the list stays as long as the live views.
        self._issued_views = [view_ref for view_ref in self._issued_views if view_ref() is not None]
        self._issued_views.append(weakref.ref(view))

    def _layer_moved(self):
        """Called by a layer view after it moved its layer's radii: the solved profile no longer lines up."""
        self._world_ptr.get().update_after_layer_change()

    @property
    def num_layers(self) -> int:
        """Number of layers in the world."""
        return self._world_ptr.get().get_num_layers()

    cdef list _ensure_layer_views(self):
        """Build the per-layer view cache once, lazily, and reuse it thereafter.

        The views are non-owning wrappers onto the world's stable C++ layers; ``add_layer`` invalidates them.
        """
        cdef size_t n, i
        cdef Layer view
        if self._layer_views is None:
            n = self._world_ptr.get().get_num_layers()
            self._layer_views = []
            self._layer_view_by_name = {}
            for i in range(n):
                view = cy_wrap_layer_view(self._world_ptr.get().get_layer(i), self)
                self._track_view(view)
                self._layer_views.append(view)
                self._layer_view_by_name[view.name] = view
        return self._layer_views

    def get_layer(self, layer) -> Layer:
        """Return a wrapper around a layer, given by its index (0 = innermost) or its name.

        The returned object is a non-owning view: the world still owns the C++ layer, so the view exposes the
        matching subclass API but must not outlive the world (it holds a reference to the world to prevent
        that). Views are cached, so repeated access is cheap. Negative indices count from the end. An index out of
        range raises ``IndexError``, and a name the world does not have ``KeyError`` (listing the layer names), as
        every method that takes a layer does.
        """
        cdef list views = self._ensure_layer_views()
        return views[cy_layer_index(self._world_ptr.get(), layer)]

    @property
    def layers(self) -> list:
        """List of (non-owning) layer views, inner to outer (cached). See :meth:`get_layer`."""
        return list(self._ensure_layer_views())

    def __len__(self):
        """Number of layers, so ``len(world)`` works."""
        return self._world_ptr.get().get_num_layers()

    def __bool__(self):
        """A world is always true, so a world with no layers (a star, say) is not taken for a missing one."""
        return True

    def __iter__(self):
        """Iterate the layer views inner-to-outer, so ``for layer in world: ...`` works."""
        return iter(self._ensure_layer_views())

    def __getitem__(self, index):
        """Index, name, or slice the layers: ``world[0]`` or ``world["mantle"]`` (a view), or ``world[:2]`` (a list of
        views).

        Integer indices accept negatives and raise ``IndexError`` out of range, and an unknown name raises
        ``KeyError`` (see :meth:`get_layer`); a slice returns the corresponding list of views.
        """
        if isinstance(index, slice):
            return self._ensure_layer_views()[index]
        return self.get_layer(index)

    def __getattr__(self, name):
        """Resolve ``world.<layer_name>`` to that layer's view (after normal attribute lookup).

        Only consulted when ``name`` is not a real attribute or method, so defined members such as
        ``calc_internal_heating`` always win. Names starting with ``_`` are never treated as layers, so
        internal probes are not intercepted.
        """
        if name.startswith("_") or self._world_ptr.get() == NULL:
            raise AttributeError(name)
        self._ensure_layer_views()
        cdef object view = self._layer_view_by_name.get(name)
        if view is not None:
            return view
        raise AttributeError(
            f"'{type(self).__name__}' object has no attribute or layer named '{name}'")

    def calc_total_mass(self) -> float:
        """Total mass [kg] = sum of all layer masses."""
        return self._world_ptr.get().calc_total_mass()

    def calc_internal_heating(self, double time) -> float:
        """Total internal radiogenic heating [W] at the given time [s].

        Only layers with a radiogenics model contribute.
        """
        return self._world_ptr.get().calc_internal_heating(time)

    def validate_layers(self) -> bool:
        """True if every layer boundary is continuous (innermost starts at 0)."""
        return self._world_ptr.get().validate_layers()

    def solve_eos(
            self,
            *,
            double surface_pressure = 0.0,
            slices_per_layer = None,
            G_to_use = None,
            integration_method = None,
            rtol = None,
            atol = None,
            pressure_tol = None,
            max_iters = None,
            nondimensionalize = None,
            temperature = None,
            solve_temperature = None,
            surface_temperature = None,
            reset_layer_masses = False,
            cpp_bool verbose = False,
            time = None,
            max_thermal_passes = None,
            thermal_tol = None,
            raise_on_fail = False) -> EOSResult:
        """Solve the whole-planet equation of state.

        Integrates gravity, pressure, enclosed mass, and moment of inertia from the planet center to its
        surface with each layer's material (see :attr:`Layer.material`), with the layer's switches, as the local
        density source; a convergence loop on the surface pressure sets the central pressure. On success every
        layer's EOS profile is populated, so :meth:`get_density`, :meth:`get_gravity`, and
        :meth:`get_pressure` work on the world and on the individual layers, and every layer's mass (and so
        its density_bulk) is set to the mass the solved density profile places between its radii.

        Every solver setting left as ``None`` takes the ``[eos_solver]`` value of the TidalPy configuration
        (``TidalPy.config``), the same defaults the standalone ``radial_solver`` uses, unless this world's file
        pinned the key (see :meth:`set_solver_defaults`). Every argument is keyword-only.

        After a successful solve whose mass differs from the world's stated ``mass`` by more than 1 percent, a
        warning (once per world) names both: the Love numbers and tides use the solved structure, while the orbit
        in a ``System``, ``calc_surface_gravity``, and ``calc_mean_density`` use the stated mass.

        Parameters
        ----------
        surface_pressure : float, optional
            Target surface pressure [Pa]. Default 0.0.
        slices_per_layer : int, optional
            Number of radial sample points generated per layer (>= 2).
        G_to_use : float, optional
            Gravitational constant [m^3 kg^-1 s^-2]. ``None`` (default) uses the TidalPy configuration's value.
        integration_method : str, optional
            CyRK integration method: ``'DOP853'``, ``'RK45'``, ``'RK23'``, ``'Tsit5'``, ``'Vern7'``, ``'Vern8'``,
            or the implicit (stiff) methods ``'BDF'``, ``'LSODA'``, ``'Radau'``. The structure ODE is singular at
            the planet's center, where LSODA's startup can fail to take its first step (a clean unsuccessful
            result) while BDF and Radau handle the singular start.
        rtol, atol : float, optional
            Relative and absolute integration tolerances.
        pressure_tol : float, optional
            Convergence tolerance on the surface-pressure mismatch, relative to the central-pressure scale
            (2/3) pi G rho^2 R^2. Keep it above ``rtol``, the integrator's own noise on the surface pressure.
        max_iters : int, optional
            Maximum central-pressure iterations. Hitting the cap sets ``max_iters_hit``; if the surface pressure
            is still off its target by more than ``pressure_tol`` the structure is not hydrostatic, so the solve
            reports ``success = False``, logs a warning, and leaves the world unsolved.
        nondimensionalize : bool, optional
            Integrate in non-dimensional units (the planet radius, its bulk density, and 1/sqrt(pi G rho) as
            the length, density, and time units) so the tolerances mean the same thing for every planet.
            Results are always returned in SI.
        temperature : float, optional
            Temperature [K] given to every layer for this solve in place of the layer's own ``temperature``.
            ``None`` uses each layer's own. It sets the viscosity, melt, and thermal-expansion density of the
            layer's material models.
        solve_temperature : bool, optional
            Carry temperature and heat flow through the structure solve, so each layer's profile follows its
            cooling model (conducting boundary layers around an adiabatic interior for convection) and its
            material models see the local temperature. Without a temperature contrast between the layers or to
            the surface there is no profile to integrate and the solve keeps the four structure variables.
            ``None`` takes the world's ``[eos_solver]`` setting, then the config's.
        surface_temperature : float, optional
            Temperature [K] the outermost layer radiates to. ``None`` leaves no flow through the surface.
        max_thermal_passes : int, optional
            Cap on the passes of a solve that carries temperature, each one structure integration followed by a
            relaxation of the boundary layers, interface temperatures, and heat flows against it. A solve that
            reaches it reports ``thermal_converged = False`` and keeps the profile of its last pass.
        thermal_tol : float, optional
            Largest relative change in the interface temperatures and heat flows between two passes at which the
            passes end, ``thermal_converged`` then being True.
        verbose : bool, optional
            Print solver status messages. Default False.
        time : float, optional
            Time [s] the heat sources of the layers with ``use_heating`` are evaluated at, on the clock their
            radiogenics models share. ``None`` takes each model's own reference time.
        reset_layer_masses : bool, optional
            A layer that holds its mass (``is_volume_fixed = False``) forgets the mass it holds and takes the
            mass its current boundaries hold in this solve. Default False. A layer that holds its mass ends where it
            encloses that mass, which the integration finds itself; the layers above it move with it, keeping their
            volumes, and the world's radius becomes the top of the outermost layer.
        raise_on_fail : bool, optional
            Raise ``SolutionFailedError`` (a ``RuntimeError``) with the solve's message when the solve fails, rather
            than returning a result with ``success = False``. Default False.

        Returns
        -------
        EOSResult
            A ``dict`` whose ``repr`` is a short summary (:class:`EOSResult`) holding ``success``, ``message``,
            ``iterations`` (central-pressure steps), ``structure_integrations`` (integrations of the whole
            structure over every thermal pass, the repeats that keep a converged pass's dense output
            included, the measure of the solve's cost), ``max_iters_hit``, ``pressure_error`` [Pa], the radial
            profile arrays (``radius``, ``gravity``, ``pressure``, ``mass``, ``moi``, ``density``,
            ``temperature``, ``heat_flow``), the scalar results (``surface_gravity``, ``surface_pressure``,
            ``central_pressure``, ``planet_mass``, ``planet_moi``), the iteration report (``thermal_passes``,
            ``thermal_converged``), the solid and liquid zones (``zones``, as :attr:`zones` reports them), and the
            per-layer results (``layer_radius_outer``
            [m], ``layer_temperature`` [K], ``layer_heat_flow_in`` and ``layer_heat_flow_out`` [W],
            ``layer_heating`` [W] (every heat source of the layers with ``use_heating``, with each source's part in
            ``layer_heating_radiogenic``, ``layer_heating_tidal``, and ``layer_heating_prescribed``; zero in a solve
            that carries no temperature), ``layer_node_temperature``,
            ``layer_top_temperature``, and ``layer_base_temperature`` [K] (the two ends of a convecting interior,
            whose top is the layer's own temperature), ``layer_boundary_thickness`` [m], ``layer_rayleigh_number``,
            ``layer_nusselt_number``, where a convecting layer evaluated its viscosity (``layer_reference_pressure``
            [Pa], ``layer_reference_viscosity`` [Pa s], and ``layer_reference_melt_fraction``, at the top of its
            interior and its own temperature; NaN for the other layers), ``layer_magma_ocean`` (a convecting
            interior liquid there, which takes the liquid scaling), ``layer_boundary_fallback`` (a convecting layer
            whose cooling model gave no boundary-layer thickness, so each took the largest share it may), and
            ``layer_in_thermal_network``). A layer's temperature rate and heat capacities are not part of the result:
            :meth:`calc_layer_temperature_rate`, :meth:`calc_layer_thermal_capacity`, and
            :meth:`calc_layer_latent_capacity` compute them from this solve when asked, since their quadratures cost
            more than an isothermal solve.

        Raises
        ------
        ValueError
            If the world has no layers, any layer lacks a material, an unsupported integration
            method is given, or ``slices_per_layer < 2``.
        SolutionFailedError
            With ``raise_on_fail``, when the solve fails.

        Assumptions
        -----------
        - Spherical symmetry; all quantities MKS.
        - Each layer's density comes from its material.
        """
        # The config struct starts from the [eos_solver] section of the TidalPy configuration with the keys this
        # world's file pinned on top (set_solver_defaults); only the arguments given here override it.
        cdef c_WorldEOSSolveConfig cfg = self._world_ptr.get().make_eos_solve_config()
        cfg.surface_pressure = surface_pressure
        # A negative G is the C++ solve's "from the configuration".
        cfg.G_to_use         = -1.0 if G_to_use is None else <double>G_to_use
        cfg.verbose          = <cpp_bool>verbose
        if temperature is not None:
            cfg.temperature = <double>temperature
        if solve_temperature is not None:
            cfg.solve_temperature = <cpp_bool>bool(solve_temperature)
        if surface_temperature is not None:
            cfg.surface_temperature = <double>surface_temperature
        if time is not None:
            cfg.time = <double>time
        cfg.reset_layer_masses = <cpp_bool>bool(reset_layer_masses)
        if slices_per_layer is not None:
            cfg.slices_per_layer = <size_t>int(slices_per_layer)
        if integration_method is not None:
            cfg.integration_method = cy_resolve_integration_method(integration_method)
        if rtol is not None:
            cfg.rtol = <double>rtol
        if atol is not None:
            cfg.atol = <double>atol
        if pressure_tol is not None:
            cfg.pressure_tol = <double>pressure_tol
        if max_iters is not None:
            cfg.max_iters = <size_t>int(max_iters)
        if nondimensionalize is not None:
            cfg.nondimensionalize = <cpp_bool>bool(nondimensionalize)
        if max_thermal_passes is not None:
            cfg.max_thermal_passes = <size_t>int(max_thermal_passes)
        if thermal_tol is not None:
            cfg.thermal_tol = <double>thermal_tol

        # Pure-C++ solve. Input validation throws std::invalid_argument, surfaced here as ValueError via the
        # ``except +`` on the C++ declaration. The result is copied out under the world's call lock in the same call,
        # so another thread's solve_eos cannot replace it before it is read.
        cdef c_WorldEOSReport report
        with nogil:
            report = self._world_ptr.get().solve_eos_report(cfg)

        if report.max_iters_hit and not report.success:
            log_warning(
                f"World '{self.name}' EOS solve stopped at max_iters = {cfg.max_iters} with a surface-pressure "
                f"mismatch above pressure_tol = {cfg.pressure_tol:0.1e}, so the world is left unsolved. Its layers "
                f"may have no hydrostatic structure at this radius and mass; otherwise raise max_iters, or keep "
                f"pressure_tol above the integration rtol ({cfg.rtol:0.1e}).")
        if raise_on_fail and not report.success:
            raise SolutionFailedError(
                f"TidalPy: the EOS solve of world '{self.name}' failed: {report.message.decode('utf-8')}")

        return cy_eos_report_to_dict(report, cy_layer_names(self._world_ptr.get()))

    def _build_eos_result(self):
        """The result (an ``EOSResult``) of the last ``solve_eos``, from a copy of the solution taken under the world's
        call lock."""
        return cy_eos_report_to_dict(self._world_ptr.get().get_eos_report(), cy_layer_names(self._world_ptr.get()))

    @property
    def eos_solved(self) -> bool:
        """True once the world-level EOS solve has populated the layer profiles."""
        return self._world_ptr.get().get_eos_solved()

    @property
    def all_materials_set(self) -> bool:
        """True once every layer has a material."""
        return self._world_ptr.get().get_all_materials_set()

    # Structure and viscoelastic profile queries (delegate to the containing layer).
    #
    # Every getter takes a scalar radius [m] (returning a float or complex) or a NumPy array of radii
    # (returning an array of the same shape), NaN where the EOS is unsolved. Each makes one C++ call for the whole
    # input, which holds the world's call lock throughout, so a read takes turns with solve_eos and the other locked
    # calls on other threads and every value of one call comes from one solve (see cy_eos_field in layers/layer.pyx).
    def _apply_complex(self, radius, double frequency, cpp_bool is_shear):
        # float -> complex; np.ndarray -> complex np.ndarray (same shape), read without the GIL.
        cdef cnp.ndarray in_arr
        cdef cnp.ndarray out_arr
        cdef double[::1] flat_in
        cdef double complex[::1] flat_out
        cdef cpp_complex[double] value
        cdef size_t num_radii
        if isinstance(radius, np.ndarray):
            in_arr  = np.ascontiguousarray(radius, dtype=np.float64)
            out_arr = np.empty_like(in_arr, dtype=np.complex128)
            flat_in = in_arr.reshape(-1)
            flat_out = out_arr.reshape(-1)
            num_radii = <size_t>flat_in.shape[0]
            # double complex and std::complex<double> share one layout.
            if num_radii > 0:
                with nogil:
                    self._world_ptr.get().calc_complex_moduli(
                        is_shear,
                        &flat_in[0],
                        num_radii,
                        frequency,
                        <cpp_complex[double]*><void*>&flat_out[0])
            return out_arr
        if is_shear:
            value = self._world_ptr.get().calc_complex_shear_modulus(<double>radius, frequency)
        else:
            value = self._world_ptr.get().calc_complex_bulk_modulus(<double>radius, frequency)
        return PyComplex_FromDoubles(value.real(), value.imag())

    def get_density(self, radius):
        """Density [kg/m^3] at radius [m] (float or np.ndarray); NaN if unsolved."""
        return cy_eos_field(<const void*>self._world_ptr.get(), cy_world_eos_fields, radius, C_EOS_DENSITY_INDEX)

    def get_gravity(self, radius):
        """Gravitational acceleration [m/s^2] at radius [m] (float or np.ndarray)."""
        return cy_eos_field(<const void*>self._world_ptr.get(), cy_world_eos_fields, radius, C_EOS_GRAVITY_INDEX)

    def get_pressure(self, radius):
        """Pressure [Pa] at radius [m] (float or np.ndarray); NaN if unsolved."""
        return cy_eos_field(<const void*>self._world_ptr.get(), cy_world_eos_fields, radius, C_EOS_PRESSURE_INDEX)

    def get_temperature(self, radius):
        """Temperature [K] at radius [m] (float or np.ndarray) from the solved profile.

        A solve with no temperature contrast reports each layer's own temperature. NaN if unsolved.
        """
        return cy_eos_field(<const void*>self._world_ptr.get(), cy_world_eos_fields, radius, C_EOS_TEMPERATURE_INDEX)

    def get_heat_flow(self, radius):
        """Heat flowing outward through the sphere of radius [m] (float or np.ndarray) \[W\].

        Zero everywhere when the solve carried no temperature. The flow steps across the interior of a
        convecting layer: that difference is the heat the layer stores or releases.
        """
        return cy_eos_field(<const void*>self._world_ptr.get(), cy_world_eos_fields, radius, C_EOS_HEAT_FLOW_INDEX)

    def get_shear_modulus(self, radius):
        """Post-melt static shear modulus [Pa] at radius [m] (float or np.ndarray)."""
        return cy_eos_field(
            <const void*>self._world_ptr.get(), cy_world_eos_fields, radius, C_EOS_SHEAR_MODULUS_INDEX)

    def get_bulk_modulus(self, radius):
        """Post-melt static bulk modulus [Pa] at radius [m] (float or np.ndarray)."""
        return cy_eos_field(
            <const void*>self._world_ptr.get(), cy_world_eos_fields, radius, C_EOS_BULK_MODULUS_INDEX)

    def get_shear_viscosity(self, radius):
        """Post-melt shear viscosity [Pa s] at radius [m] (float or np.ndarray)."""
        return cy_eos_field(
            <const void*>self._world_ptr.get(), cy_world_eos_fields, radius, C_EOS_SHEAR_VISCOSITY_INDEX)

    def get_bulk_viscosity(self, radius):
        """Post-melt bulk viscosity [Pa s] at radius [m] (float or np.ndarray)."""
        return cy_eos_field(
            <const void*>self._world_ptr.get(), cy_world_eos_fields, radius, C_EOS_BULK_VISCOSITY_INDEX)

    def get_melt_fraction(self, radius):
        """Melt fraction at radius [m] (float or np.ndarray): 0.0 where the layer's material does not melt (or its
        layer does not use melting), 1.0 in a liquid-only material."""
        return cy_eos_field(
            <const void*>self._world_ptr.get(), cy_world_eos_fields, radius, C_EOS_MELT_FRACTION_INDEX)

    def calc_complex_shear_modulus(self, radius, double frequency, cpp_bool recalc_eos=False):
        """Complex shear modulus [Pa] at radius [m] (float or np.ndarray) and frequency [rad/s].

        Applies the containing layer's shear rheology to the stored post-melt static
        modulus + viscosity. Solves the EOS first if it has not been solved (or if
        ``recalc_eos``). Without a rheology it is the static modulus as a purely real number.
        """
        self._ensure_solved(recalc_eos)
        return self._apply_complex(radius, frequency, True)

    def calc_complex_bulk_modulus(self, radius, double frequency, cpp_bool recalc_eos=False):
        """Complex bulk modulus [Pa] at radius [m] (float or np.ndarray) and frequency [rad/s].

        Applies the containing layer's bulk rheology to the stored post-melt static
        modulus + viscosity. Solves the EOS first if it has not been solved (or if
        ``recalc_eos``). Without a rheology it is the static modulus as a purely real number.
        """
        self._ensure_solved(recalc_eos)
        return self._apply_complex(radius, frequency, False)

    # Shorthand bundles (one call returns several profiles at once).
    def get_static_viscoelastics(self, radius):
        """``(shear_modulus, shear_viscosity, bulk_modulus, bulk_viscosity)`` (post-melt) at radius.

        Each element is a float (scalar radius) or np.ndarray (array of radii). One evaluation of the solved state
        per radius fills all four.
        """
        return cy_eos_fields(
            <const void*>self._world_ptr.get(), cy_world_eos_fields, radius,
            (C_EOS_SHEAR_MODULUS_INDEX, C_EOS_SHEAR_VISCOSITY_INDEX,
             C_EOS_BULK_MODULUS_INDEX, C_EOS_BULK_VISCOSITY_INDEX))

    def get_state(self, radius):
        """All EOS-related profiles at radius as a dict (float or np.ndarray values), from one evaluation of the
        solved state per radius."""
        values = cy_eos_fields(
            <const void*>self._world_ptr.get(), cy_world_eos_fields, radius,
            (C_EOS_DENSITY_INDEX, C_EOS_GRAVITY_INDEX, C_EOS_PRESSURE_INDEX, C_EOS_SHEAR_MODULUS_INDEX,
             C_EOS_SHEAR_VISCOSITY_INDEX, C_EOS_BULK_MODULUS_INDEX, C_EOS_BULK_VISCOSITY_INDEX,
             C_EOS_MELT_FRACTION_INDEX))
        return dict(zip(("density", "gravity", "pressure", "shear_modulus", "shear_viscosity", "bulk_modulus",
                         "bulk_viscosity", "melt_fraction"), values))

    # calc_* variants: solve the EOS first if it is unsolved (or force_recalc), then read the profile.
    def _ensure_solved(self, cpp_bool force_recalc):
        """Solve the EOS when asked to or when the world is unsolved, raising SolutionFailedError (with the solve's
        message) when that solve fails, rather than letting the calc_* getters return NaN profiles."""
        if force_recalc or not self._world_ptr.get().get_eos_solved():
            self.solve_eos(raise_on_fail=True)

    def calc_density(self, radius, cpp_bool force_recalc=False):
        """Density [kg/m^3]; solves the EOS first if needed."""
        self._ensure_solved(force_recalc)
        return self.get_density(radius)

    def calc_gravity(self, radius, cpp_bool force_recalc=False):
        """Gravitational acceleration [m/s^2]; solves the EOS first if needed."""
        self._ensure_solved(force_recalc)
        return self.get_gravity(radius)

    def calc_pressure(self, radius, cpp_bool force_recalc=False):
        """Pressure [Pa]; solves the EOS first if needed."""
        self._ensure_solved(force_recalc)
        return self.get_pressure(radius)

    def calc_shear_modulus(self, radius, cpp_bool force_recalc=False):
        """Post-melt shear modulus [Pa]; solves the EOS first if needed."""
        self._ensure_solved(force_recalc)
        return self.get_shear_modulus(radius)

    def calc_bulk_modulus(self, radius, cpp_bool force_recalc=False):
        """Post-melt bulk modulus [Pa]; solves the EOS first if needed."""
        self._ensure_solved(force_recalc)
        return self.get_bulk_modulus(radius)

    def calc_shear_viscosity(self, radius, cpp_bool force_recalc=False):
        """Post-melt shear viscosity [Pa s]; solves the EOS first if needed."""
        self._ensure_solved(force_recalc)
        return self.get_shear_viscosity(radius)

    def calc_bulk_viscosity(self, radius, cpp_bool force_recalc=False):
        """Post-melt bulk viscosity [Pa s]; solves the EOS first if needed."""
        self._ensure_solved(force_recalc)
        return self.get_bulk_viscosity(radius)

    def calc_static_viscoelastics(self, radius, cpp_bool force_recalc=False):
        """``get_static_viscoelastics`` but solves the EOS first if needed."""
        self._ensure_solved(force_recalc)
        return self.get_static_viscoelastics(radius)

    def calc_state(self, radius, cpp_bool force_recalc=False):
        """``get_state`` but solves the EOS first if needed."""
        self._ensure_solved(force_recalc)
        return self.get_state(radius)

    @property
    def surface_gravity_eos(self) -> float:
        """Surface gravity [m/s^2] from the last EOS solve, or NaN if not solved."""
        return self._world_ptr.get().get_surface_gravity_eos()

    @property
    def central_pressure(self) -> float:
        """Central pressure [Pa] from the last EOS solve, or NaN if not solved."""
        return self._world_ptr.get().get_central_pressure()

    @property
    def planet_mass_eos(self) -> float:
        """Total planet mass [kg] integrated by the last EOS solve, or NaN if not solved."""
        return self._world_ptr.get().get_planet_mass_eos()

    @property
    def planet_moi_eos(self) -> float:
        """Planet moment of inertia [kg m^2] from the last EOS solve, or NaN if not solved."""
        return self._world_ptr.get().get_planet_moi_eos()

    @property
    def moment_of_inertia_factor(self) -> float:
        """C / (M R^2) [dimensionless] of the solved structure: :attr:`planet_moi_eos` over :attr:`planet_mass_eos`
        times the radius squared. NaN before a successful :meth:`solve_eos`.

        Before a solve, a spin calculation uses the spin model's factor instead (the ``moment_of_inertia_factor`` of
        a world file or of ``[worlds]`` in the TidalPy configuration; see :meth:`get_moment_of_inertia`), which this
        property does not report.
        """
        return self._world_ptr.get().get_moment_of_inertia_factor()

    @property
    def zones(self) -> list:
        """The solid and liquid zones of the last EOS solve, inner to outer, as dicts with ``layer`` (the layer's
        name), ``radius_inner`` and ``radius_outer`` [m], ``mass_inner`` and ``mass_outer`` (the enclosed mass at each
        end) [kg], and ``state`` (``"solid"`` or ``"liquid"``).

        Each layer is one zone unless its material changed state inside it (a layer with ``use_melting`` and
        ``state = "auto"`` whose material melts). The EOS solve finds each boundary as it integrates: where the
        material's post-melt rigidity mu / (rho g R) (the world's stated bulk density, surface gravity, and radius)
        crosses the ``[numerical]`` ``minimum_solid_rigidity`` of the TidalPy configuration (by default 0, so a zone is
        liquid only where the shear modulus vanishes and partial melt stays a weakened solid). A zone thinner than
        ``[numerical] minimum_zone_fraction`` of the world's radius takes the state of its thicker neighbor. The radial
        solver integrates each zone as a layer of its own, with the liquid equations in a liquid zone and the layer's
        ``is_static`` and ``is_incompressible`` choices. Empty before an EOS solve.
        """
        return cy_zones_to_list(self._world_ptr.get().get_zones_copy(), cy_layer_names(self._world_ptr.get()))

    @property
    def molten_regions(self) -> list:
        """The liquid zones of layers that are not liquid throughout (the stretches melting made liquid, see
        :attr:`zones`), as ``(layer_name, radius_inner, radius_outer)`` with radii in [m]. Empty before an EOS solve
        or when nothing is molten.
        """
        cdef vector[c_EOSZone] zones = self._world_ptr.get().get_molten_regions()
        cdef size_t zone_i
        cdef list layer_names = cy_layer_names(self._world_ptr.get())
        regions = []
        for zone_i in range(zones.size()):
            regions.append((layer_names[zones[zone_i].layer_index], zones[zone_i].radius_inner,
                            zones[zone_i].radius_outer))
        return regions

    # Spin dynamics (the Spin model attached to the world; uses the world's EOS moment of inertia)
    def set_spin_model(self, Spin spin not None):
        """Attach a :class:`~TidalPy.Dynamics.Spin` model (its moment-of-inertia factor is the fallback
        when the EOS has not been solved)."""
        self._world_ptr.get().set_spin_model(spin._spin)

    def get_moment_of_inertia(self) -> float:
        """Moment of inertia [kg m^2]: the EOS-solved value when the EOS has been solved, else the spin
        model's ``moment_of_inertia_factor * M R^2`` estimate from the world mass and radius."""
        return self._world_ptr.get().get_moment_of_inertia()

    def calc_spin_derivative(self, double host_mass) -> float:
        """Tidal spin-rate change [rad s-2] ``= M_host * dU/dO / I`` using the world's stored ``dU/dO``
        (from the last :meth:`calc_tides`) and its moment of inertia. Requires a completed tidal solve."""
        return self._world_ptr.get().calc_spin_derivative(host_mass)

    def calc_synchronous_spin(self, double orbital_frequency) -> float:
        """Synchronous spin rate [rad s-1]: equal to the orbital mean motion ``orbital_frequency``."""
        return self._world_ptr.get().calc_synchronous_spin(orbital_frequency)

    def solve_love_numbers(
            self,
            double frequency  = 1.0e-5,
            int degree_l      = 2,
            solve_for         = 'tidal',
            int core_model    = 0,
            starting_method   = None,
            nondimensionalize = None,
            double starting_radius = 0.0,
            start_radius_tol   = None,
            integration_method = None,
            rtol               = None,
            atol               = None,
            scale_rtols        = None,
            max_num_steps      = None,
            expected_size      = None,
            max_ram_MB         = None,
            double max_step    = 0.0,
            cpp_bool verbose   = False,
            cpp_bool warnings  = True,
            love_method        = None,
            fixed_q            = None,
            fixed_dt           = None,
            raise_on_fail      = False,
            cpp_bool love_only = False,
            degree1_frame      = None) -> dict:
        """Solve for whole-planet tidal Love numbers (radial solver, propagation matrix, or analytic methods).

        Requires :meth:`solve_eos` first. Each layer's attached rheology is evaluated at ``frequency`` for the
        complex moduli, then the deformation ODEs are shot from the center to the surface on the solved
        structure to give k, h, l.

        Every solver setting left as ``None`` takes the ``[radial_solver]`` value of the TidalPy configuration
        (``TidalPy.config``), the same defaults the standalone ``radial_solver`` and the world's tidal
        solves use.

        Parameters
        ----------
        frequency : float, optional
            Tidal forcing frequency [rad/s]. Default 1e-5.
        degree_l : int, optional
            Harmonic degree. Default 2. It must be 2 or more, or 1 for a solve for loading alone.
        solve_for : str or sequence of str, optional
            Surface boundary condition: ``'tidal'`` (default; tidal Love numbers k, h, l), ``'loading'``
            (load Love numbers k', h', l'), or ``'free'`` (free-surface response). Same names as the
            standalone ``radial_solver``. A sequence solves for several in one pass, which costs one
            integration rather than one each because the independent solutions do not depend on the
            boundary condition; the Love-number properties then report the first entry.
        core_model : int, optional
            Propagation-matrix core starting condition (0-4). Ignored by the shooting method. Default 0.
        starting_method : str, optional
            Starting conditions of the shooting method: ``'takeuchi'`` (alias ``'ts'``), ``'kamata'``,
            ``'power_series'`` (aliases ``'powerseries'``, ``'ps'``, ``'martens'``), or ``'unity'``. Ignored by the
            other Love methods.
        nondimensionalize : bool, optional
            Non-dimensionalize the problem internally (recommended).
        starting_radius : float, optional
            Minimum radius [m] for the shooting start; 0 selects it automatically. Default 0.
        start_radius_tol : float, optional
            Tolerance of the automatic starting radius, R * tol^(1/l).
        integration_method : str, optional
            CyRK ODE method: ``'DOP853'``, ``'RK45'``, ``'RK23'``, ``'Tsit5'``, ``'Vern7'``, ``'Vern8'``, or the
            implicit (stiff) methods ``'BDF'``, ``'LSODA'``, ``'Radau'``.
        rtol, atol : float, optional
            Relative and absolute ODE tolerances.
        scale_rtols : bool, optional
            Scale tolerances by layer type.
        max_num_steps : int, optional
            Maximum ODE steps.
        expected_size, max_ram_MB : int, optional
            CyRK memory hints.
        max_step : float, optional
            Maximum ODE step size [m]; 0 selects it automatically. Default 0.
        verbose : bool, optional
            Print solver status messages. Default False.
        warnings : bool, optional
            Emit solver warnings, and log a failed solve's message when ``raise_on_fail`` is off. Default True.
        love_method : str, optional
            How the Love numbers are obtained. ``None`` (default) uses the world's configured method
            (``set_tide_config(love_method=...)`` or the ``[tides]`` table). ``'radial_solver'``
            (``'shooting'``, ``'rs'``) integrates the radial ODEs from the center to the surface.
            ``'propagation_matrix'`` (``'prop_matrix'``, ``'pm'``, ``'prop'``) is valid only for a single
            solid, static, incompressible layer; an incompatible world fails the solve gracefully, with
            ``love_success`` False and a non-zero ``love_error_code``. ``'homogeneous'`` (``'homogen'``)
            applies the homogeneous-sphere formulas with the volume-averaged complex shear modulus of the
            tidal layers (``use_tides``), the planet's bulk density, surface gravity, and radius. ``'cpl'``
            and ``'ctl'`` apply those to the static shear modulus with a constant phase lag ``(1 - i/Q)`` or
            time lag ``(1 - i*omega*dt)``. ``'laterally_inhomogeneous'`` (``'3d'``, ``'lat_inhom'``) raises
            ``NotImplementedError``.
        fixed_q, fixed_dt : float, optional
            Quality factor for ``'cpl'`` and time lag [s] for ``'ctl'``. Left unset, the ``[tides]`` config
            values (``set_tide_config(love_fixed_q=..., love_fixed_dt=...)``, the ``love_fixed_q`` and
            ``love_fixed_dt_s`` keys of a world file) are used, then the attached tide
            model's fixed Q or time lag for this degree (``ValueError`` when none is available).
        raise_on_fail : bool, optional
            Raise ``SolutionFailedError`` (a ``RuntimeError``) with the solve's message when the solve fails, rather
            than returning a result with ``success = False``. Default False.
        love_only : bool, optional
            Keep only what the Love numbers need (shooting method): the integration builds no dense output, which
            makes the solve faster, and the radial functions below the surface are then unavailable
            (:meth:`get_love_radial_y` and a released solution's radial functions raise ``ValueError``). The surface
            values (:meth:`get_love_surface_y`) and the Love numbers are unaffected. Default False.
        degree1_frame : str, optional
            Reference frame of degree-1 load Love numbers (Blewitt 2003): ``'CE'`` (center of mass of the solid
            body, k' = 0), ``'CM'`` (center of mass of the body and the load, 1 + k' = 0), ``'CF'`` (center of
            surface figure), ``'CL'`` (center of lateral figure, l' = 0), or ``'CH'`` (center of height figure,
            h' = 0), in any case. The Love numbers and the radial functions are both in this frame. Read only by a
            degree-1 solve for loading; ``CF`` and ``CL`` need l', which a static liquid surface layer does not
            define (error code -16).

        Returns
        -------
        dict
            ``success`` (bool), ``error_code``, ``message``, ``love_method`` (canonical name),
            ``love_number_k``, ``love_number_h``, ``love_number_l`` (complex).

        Raises
        ------
        ValueError
            If the EOS has not been solved or the frequency is too small.
        SolutionFailedError
            With ``raise_on_fail``, when the solve fails.

        Assumptions
        -----------
        - Spherical symmetry; all quantities MKS.
        """
        # The config starts from the world's [tides] Love settings (method, fixed Q, fixed time lag) and the
        # [radial_solver] section of the TidalPy configuration with the keys this world's file pinned on top
        # (set_solver_defaults); only the arguments given here override it.
        cdef c_LoveSolveConfig cfg = self._world_ptr.get().make_love_solve_config()
        cfg.frequency = frequency
        cfg.degree_l  = degree_l
        cy_set_solve_for(&cfg, solve_for)
        if love_method is not None:
            cfg.love_method = cy_parse_love_method(love_method)
        if fixed_q is not None:
            cfg.fixed_q = <double>fixed_q
        if fixed_dt is not None:
            cfg.fixed_dt = <double>fixed_dt
        cfg.core_model      = core_model
        cfg.starting_radius = starting_radius
        cfg.max_step        = max_step
        cfg.verbose         = <cpp_bool>verbose
        cfg.warnings        = <cpp_bool>warnings
        cfg.love_only       = love_only
        cy_apply_love_solve_overrides(
            &cfg, starting_method, degree1_frame, nondimensionalize, start_radius_tol, integration_method, rtol, atol,
            scale_rtols, max_num_steps, expected_size, max_ram_MB)

        with nogil:
            self._world_ptr.get().solve_love_numbers(cfg)

        # The conditioning diagnostic belongs to the radial solvers.
        if warnings and c_love_method_uses_radial_solver_int(cfg.love_method):
            cy_check_surface_solve_conditioning(
                self._world_ptr.get().get_love_surface_amplification(), cfg.rtol,
                self._world_ptr.get().get_love_surface_rcond())
            cy_check_degree1_frame_residual(self._world_ptr.get().get_love_surface_frame_residual())
        result = self._build_love_result()
        if not result["success"]:
            failure = (f"TidalPy: the Love-number solve of world '{self.name}' at frequency {frequency:.6e} rad s-1 "
                       f"failed: {result['message'].strip()}")
            if raise_on_fail:
                raise SolutionFailedError(failure)
            if warnings:
                log_warning(failure)
        return result

    def calc_love_numbers(self, frequencies, **solve_love_numbers_kwargs) -> dict:
        """Love numbers over a set of forcing frequencies, one :meth:`solve_love_numbers` per frequency.

        A solve that fails gives NaN at its frequency and does not stop the sweep; only bad input raises.

        Parameters
        ----------
        frequencies : float or array-like
            Forcing frequencies [rad s-1].
        **solve_love_numbers_kwargs
            Any other :meth:`solve_love_numbers` argument (``degree_l``, ``love_method``, ``rtol``, ...), the same for
            every frequency.

        Returns
        -------
        dict
            ``frequency`` [rad s-1] and, each shaped like it, ``k``, ``h``, ``l`` (complex; NaN where a solve failed,
            and the first entry of a ``solve_for`` sequence), ``success`` (bool), and ``message`` (str objects).

        Raises
        ------
        TypeError
            ``frequency`` or ``raise_on_fail`` given as a keyword.
        ValueError
            Bad input to a solve, such as an unsolved EOS or a frequency that is too small.
        """
        for keyword in ("frequency", "raise_on_fail"):
            if keyword in solve_love_numbers_kwargs:
                raise TypeError(f"TidalPy: calc_love_numbers takes no '{keyword}' argument.")
        cdef cnp.ndarray frequency_array = np.array(frequencies, dtype=np.float64)
        cdef cnp.ndarray flat_frequencies = frequency_array.reshape(-1)
        cdef Py_ssize_t num_frequencies = flat_frequencies.shape[0]
        cdef dict out = {
            "frequency": frequency_array,
            "k": np.full(num_frequencies, np.nan, dtype=np.complex128),
            "h": np.full(num_frequencies, np.nan, dtype=np.complex128),
            "l": np.full(num_frequencies, np.nan, dtype=np.complex128),
            "success": np.zeros(num_frequencies, dtype=np.bool_),
            "message": np.empty(num_frequencies, dtype=object),
        }
        cdef Py_ssize_t frequency_i
        for frequency_i in range(num_frequencies):
            result = self.solve_love_numbers(frequency=flat_frequencies[frequency_i], **solve_love_numbers_kwargs)
            out["success"][frequency_i] = result["success"]
            out["message"][frequency_i] = result["message"]
            if result["success"]:
                out["k"][frequency_i] = result["love_number_k"]
                out["h"][frequency_i] = result["love_number_h"]
                out["l"][frequency_i] = result["love_number_l"]
        for key in ("k", "h", "l", "success", "message"):
            out[key] = out[key].reshape(np.shape(frequency_array))
        return out

    def solve_love_numbers_supplied(
            self,
            double complex[::1] complex_shear_modulus not None,
            double complex[::1] complex_bulk_modulus not None,
            double[::1] radius_array not None,
            double frequency  = 1.0e-5,
            int    degree_l   = 2,
            solve_for         = 'tidal',
            int    core_model = 0,
            starting_method   = None,
            nondimensionalize = None,
            double starting_radius = 0.0,
            start_radius_tol   = None,
            integration_method = None,
            rtol               = None,
            atol               = None,
            scale_rtols        = None,
            max_num_steps      = None,
            expected_size      = None,
            max_ram_MB         = None,
            double max_step    = 0.0,
            cpp_bool verbose   = False,
            cpp_bool warnings  = True,
            str love_method    = 'radial_solver',
            cpp_bool love_only = False,
            degree1_frame      = None) -> dict:
        """Solve Love numbers from externally-supplied complex moduli arrays (instead of layer rheology).

        The supplied shear/bulk moduli [Pa] are defined at ``radius_array`` [m] and are linearly interpolated onto
        the world's internal EOS radius grid. Used by the standalone ``RadialSolver.radial_solver`` API.
        ``solve_eos`` must be called first. Only the radial-solver methods (``love_method``
        ``'radial_solver'`` or ``'propagation_matrix'``) are available here. Solver settings left as ``None``
        take the ``[radial_solver]`` values of the TidalPy configuration, as in :meth:`solve_love_numbers`, which
        describes ``love_only`` and ``degree1_frame``.
        """
        if radius_array.shape[0] == 0:
            raise ValueError("radius_array must not be empty")
        if (complex_shear_modulus.shape[0] != radius_array.shape[0]
                or complex_bulk_modulus.shape[0] != radius_array.shape[0]):
            raise ValueError("complex moduli and radius arrays must have matching length")

        # The world's pinned [radial_solver] keys first, as in solve_love_numbers; the method is always the argument.
        cdef c_LoveSolveConfig cfg = self._world_ptr.get().make_love_solve_config()
        cfg.frequency   = frequency
        cfg.degree_l    = degree_l
        cy_set_solve_for(&cfg, solve_for)
        cfg.love_method = cy_parse_love_method(love_method)
        cfg.core_model  = core_model
        cfg.starting_radius = starting_radius
        cfg.max_step        = max_step
        cfg.verbose         = <cpp_bool>verbose
        cfg.warnings        = <cpp_bool>warnings
        cfg.love_only       = love_only
        cy_apply_love_solve_overrides(
            &cfg, starting_method, degree1_frame, nondimensionalize, start_radius_tol, integration_method, rtol, atol,
            scale_rtols, max_num_steps, expected_size, max_ram_MB)

        cdef size_t n_in = radius_array.shape[0]
        cdef cpp_complex[double]* shear_ptr = <cpp_complex[double]*><void*>&complex_shear_modulus[0]
        cdef cpp_complex[double]* bulk_ptr  = <cpp_complex[double]*><void*>&complex_bulk_modulus[0]
        cdef double* radius_ptr = &radius_array[0]
        with nogil:
            self._world_ptr.get().solve_love_numbers_supplied(
                cfg,
                shear_ptr,
                bulk_ptr,
                radius_ptr,
                n_in)
        if warnings:
            cy_check_surface_solve_conditioning(
                self._world_ptr.get().get_love_surface_amplification(), cfg.rtol,
                self._world_ptr.get().get_love_surface_rcond())
            cy_check_degree1_frame_residual(self._world_ptr.get().get_love_surface_frame_residual())
        return self._build_love_result()

    def release_radial_solution(self):
        """Move the last radial Love solve's storage out of this world into a ``RadialSolverSolution``.

        The solution owns the storage afterwards, so this world's radial cache is emptied and the next
        ``solve_love_numbers`` call rebuilds it. The solution holds a reference back to this world, which is
        what lets its interior getters keep answering: they read the state provider the solve installed, and
        that provider reads this world's solved EOS.

        Returns
        -------
        RadialSolverSolution
            The solved radial functions, Love numbers, and interior of the last radial solve.

        Raises
        ------
        RuntimeError
            If no radial Love solve has been run, or the last one did not use a radial method.
        """
        cdef unique_ptr[c_RadialSolutionStorage] storage_uptr = self._world_ptr.get().release_radial_storage()
        if not storage_uptr:
            raise RuntimeError(
                "No radial solution to release: run solve_love_numbers with the 'radial_solver' or "
                "'propagation_matrix' method first.")
        return RadialSolverSolution._adopt(move(storage_uptr), self)

    def _build_love_result(self):
        """Assemble the Python result dict from the retained C++ Love-number solution."""
        cdef cpp_complex[double] k, h, l
        k = self._world_ptr.get().get_love_number_k(0)
        h = self._world_ptr.get().get_love_number_h(0)
        l = self._world_ptr.get().get_love_number_l(0)
        return {
            'success':       bool(self._world_ptr.get().get_love_solved()),
            'error_code':    self._world_ptr.get().get_love_error_code(),
            'message':       self._world_ptr.get().get_love_message().decode('utf-8'),
            'love_method':   self.love_method,
            'love_number_k': PyComplex_FromDoubles(k.real(), k.imag()),
            'love_number_h': PyComplex_FromDoubles(h.real(), h.imag()),
            'love_number_l': PyComplex_FromDoubles(l.real(), l.imag()),
        }


    @property
    def love_solved(self) -> bool:
        """True while a successful ``solve_love_numbers`` result is held.

        ``solve_eos`` clears it, because Love numbers describe the structure they were solved with; the Love
        number getters return NaN until the next solve.
        """
        return bool(self._world_ptr.get().get_love_solved())

    @property
    def love_success(self) -> bool:
        """True if the last love-number solve reported success."""
        return bool(self._world_ptr.get().get_love_success())

    @property
    def love_error_code(self) -> int:
        """Integer error code from the last love-number solve (-100 if not yet run)."""
        return self._world_ptr.get().get_love_error_code()

    @property
    def love_message(self) -> str:
        """Status message from the last love-number solve."""
        return self._world_ptr.get().get_love_message().decode('utf-8')

    @property
    def love_surface_amplification(self) -> float:
        """Worst-case error amplification of the last solve's surface boundary condition collapse.

        Values near 1 indicate a well conditioned solve; large values mean roundoff and integration error are
        amplified into the Love numbers (the achievable relative accuracy is floor limited to about this value
        times machine epsilon). 0 until a shooting-method solve has run with ``warnings`` enabled; the propagation
        matrix method does not use the shooting surface collapse.
        """
        return self._world_ptr.get().get_love_surface_amplification()

    @property
    def love_surface_rcond(self) -> float:
        """Reciprocal condition number of the last solve's surface boundary condition system.

        Units- and normalization-independent: near 1 is well posed, and a value near machine epsilon means the
        solution constants are undetermined (the solve then fails with error code -13 once it is below
        ``[numerical] minimum_surface_rcond``). A value below the integration rtol draws a conditioning warning.
        NaN before a shooting-method solve and for the other Love methods. The standalone solver reports the same
        number as ``RadialSolverSolution.surface_solve_rcond``.
        """
        return self._world_ptr.get().get_love_surface_rcond()

    @property
    def love_surface_frame_residual(self) -> float:
        """How far the last degree-1 loading solve leaves the surface condition its reference frame replaced.

        Roundoff when every layer is static; about omega^2 R / g when some layers carry inertia and others do not,
        and warned about above 1e-2. NaN for any other solve. The standalone solver reports the same number as
        ``RadialSolverSolution.surface_frame_residual``.
        """
        return self._world_ptr.get().get_love_surface_frame_residual()

    @property
    def love_method(self) -> str:
        """Canonical name of the Love-number method used by the last solve (``'radial_solver'`` before any)."""
        return c_love_method_name_int(self._world_ptr.get().get_love_method_last_int()).decode('utf-8')

    @property
    def love_effective_shear_modulus(self) -> complex:
        """Volume-averaged shear modulus [Pa] of the tidal layers used by the last analytic Love solve.

        NaN after a radial-solver solve. Complex for the ``homogeneous`` method (evaluated at the forcing
        frequency); real (static) for ``cpl`` and ``ctl``.
        """
        cdef cpp_complex[double] v = self._world_ptr.get().get_love_analytic_shear()
        return PyComplex_FromDoubles(v.real(), v.imag())

    @property
    def love_tidal_volume(self) -> float:
        """Volume [m3] of the tidal layers that took part in the last quasi-homogeneous Love solve (NaN otherwise)."""
        return self._world_ptr.get().get_love_analytic_tidal_volume()

    @property
    def love_layer_parts(self) -> list:
        """Each tidal layer's part of the last quasi-homogeneous Love solve (``homogeneous``, ``cpl``, ``ctl``).

        One dict per layer that took part: ``layer`` (its name), ``tidal_scale``, the Love numbers of a homogeneous
        planet made of the layer's averaged material (``love_number_k``, ``love_number_h``, ``love_number_l``), and
        its complex ``shear_modulus`` [Pa] at the solve's frequency. The world's Love numbers are the sum of
        ``tidal_scale`` times these. Empty after a radial-solver solve.
        """
        cdef list parts = []
        # A copy taken under the world's call lock, so another thread's Love solve cannot change it mid-loop.
        cdef vector[c_LayerLove] layer_parts = self._world_ptr.get().get_love_layer_parts()
        cdef size_t i
        for i in range(layer_parts.size()):
            parts.append({
                "layer":         self._world_ptr.get().get_layer(layer_parts[i].layer_index).get_name().decode("utf-8"),
                "tidal_scale":   layer_parts[i].tidal_scale,
                "love_number_k": PyComplex_FromDoubles(layer_parts[i].love.k.real(), layer_parts[i].love.k.imag()),
                "love_number_h": PyComplex_FromDoubles(layer_parts[i].love.h.real(), layer_parts[i].love.h.imag()),
                "love_number_l": PyComplex_FromDoubles(layer_parts[i].love.l.real(), layer_parts[i].love.l.imag()),
                "shear_modulus": PyComplex_FromDoubles(
                    layer_parts[i].shear_modulus.real(), layer_parts[i].shear_modulus.imag()),
            })
        return parts

    @property
    def love_num_ytypes(self) -> int:
        """Number of boundary-condition types solved (typically 1 for tidal-only)."""
        return int(self._world_ptr.get().get_love_num_ytypes())

    @property
    def love_number_k(self) -> complex:
        """Complex potential Love number k2 from the last radial solve (NaN+0j if unsolved)."""
        cdef cpp_complex[double] v = self._world_ptr.get().get_love_number_k(<size_t>0)
        return PyComplex_FromDoubles(v.real(), v.imag())

    @property
    def love_q_k(self) -> float:
        """Quality factor of the last solve's k [dimensionless]: Q = -s |k| / Im(k), with s the sign of Re(k).

        The definition ``RadialSolverSolution.Q_k`` uses (``c_LoveNumbers::get_Q_k``); infinite for a purely elastic
        k, NaN when :attr:`love_number_k` is.
        """
        return self._world_ptr.get().get_love_q_k(<size_t>0)

    @property
    def love_lag_k(self) -> float:
        """Phase lag of the last solve's k [rad]: arctan2(-s Im(k), |Re(k)|), with s the sign of Re(k).

        The definition ``RadialSolverSolution.lag_k`` uses (``c_LoveNumbers::get_lag_k``); 0 for a purely elastic
        k, NaN when :attr:`love_number_k` is.
        """
        return self._world_ptr.get().get_love_lag_k(<size_t>0)

    @property
    def love_number_h(self) -> complex:
        """Complex radial displacement Love number h2 from the last radial solve (NaN+0j if unsolved)."""
        cdef cpp_complex[double] v = self._world_ptr.get().get_love_number_h(<size_t>0)
        return PyComplex_FromDoubles(v.real(), v.imag())

    @property
    def love_number_l(self) -> complex:
        """Complex tangential (Shida) Love number l2 from the last radial solve (NaN+0j if unsolved)."""
        cdef cpp_complex[double] v = self._world_ptr.get().get_love_number_l(<size_t>0)
        return PyComplex_FromDoubles(v.real(), v.imag())

    def get_love_number_k(self, ytype_idx: int = 0) -> complex:
        """Complex k Love number for the given boundary-condition ytype index."""
        cdef cpp_complex[double] v = self._world_ptr.get().get_love_number_k(<size_t>ytype_idx)
        return PyComplex_FromDoubles(v.real(), v.imag())

    def get_love_number_h(self, ytype_idx: int = 0) -> complex:
        """Complex h Love number for the given boundary-condition ytype index."""
        cdef cpp_complex[double] v = self._world_ptr.get().get_love_number_h(<size_t>ytype_idx)
        return PyComplex_FromDoubles(v.real(), v.imag())

    def get_love_number_l(self, ytype_idx: int = 0) -> complex:
        """Complex l (Shida) Love number for the given boundary-condition ytype index."""
        cdef cpp_complex[double] v = self._world_ptr.get().get_love_number_l(<size_t>ytype_idx)
        return PyComplex_FromDoubles(v.real(), v.imag())

    def get_love_radial_y(self, radius, ytype_idx: int = 0, y_idx: int = 0):
        """Radial function y[y_idx + 1] (SI) at ``radius`` from the last radial-solver Love solve.

        The shooting method evaluates its dense per-layer interpolants at the radius; the propagation
        matrix interpolates its grid. NaN if unsolved, after an analytic (homogeneous/cpl/ctl) solve, out
        of range, or below the solver's starting radius. ``y_idx`` 0..5 selects y1..y6.

        Parameters
        ----------
        radius : float or np.ndarray
            Radius [m].
        ytype_idx : int, optional
            The boundary-condition entry of the solve's ``solve_for``. Default 0.
        y_idx : int, optional
            0..5 for y1..y6. Default 0.

        Returns
        -------
        complex or np.ndarray
            A complex for a scalar radius, else a complex array shaped like ``radius``.

        Raises
        ------
        ValueError
            After a ``love_only`` solve, which keeps no radial functions below the surface.
        """
        cdef c_BaseWorld* world_ptr = self._world_ptr.get()
        if world_ptr.get_love_only():
            raise ValueError(
                f"TidalPy: the last Love solve of world '{self.name}' ran with love_only=True, which keeps no radial "
                "functions below the surface; solve again with love_only=False.")
        cdef size_t ytype = <size_t>ytype_idx
        cdef size_t y_index = <size_t>y_idx
        cdef cpp_complex[double] v
        if np.ndim(radius) == 0:
            v = world_ptr.get_radial_solution_y(<double>radius, ytype, y_index)
            return PyComplex_FromDoubles(v.real(), v.imag())
        cdef cnp.ndarray radius_array = np.ascontiguousarray(radius, dtype=np.float64)
        cdef cnp.ndarray out = np.empty(np.shape(radius_array), dtype=np.complex128)
        cdef double[::1] radius_view = radius_array.reshape(-1)
        cdef double complex[::1] out_view = out.reshape(-1)
        cdef Py_ssize_t radius_i
        for radius_i in range(radius_view.shape[0]):
            v = world_ptr.get_radial_solution_y(radius_view[radius_i], ytype, y_index)
            out_view[radius_i] = v.real() + 1j * v.imag()
        return out

    def get_love_surface_y(self, ytype_idx: int, y_idx: int) -> complex:
        """Complex radial y-solution value at the surface for the given ytype and y index."""
        cdef cpp_complex[double] v = self._world_ptr.get().get_love_surface_y(
            <size_t>ytype_idx, <size_t>y_idx)
        return PyComplex_FromDoubles(v.real(), v.imag())

    # Global (1D) tidal dissipation
    def calc_tides(
            self,
            orbital_frequency=None,
            spin_frequency=None,
            eccentricity=None,
            obliquity=None,
            semi_major_axis=None,
            host_mass=None) -> dict:
        """Solve the global tidal dissipation for the given orbital/spin state.

        Requires an attached tide model (:meth:`set_tide_model`). The analytic models (cpl/ctl/ctl_q) collapse from
        their fixed per-degree parameters. The rheology model instead runs the world radial solver at each unique
        tidal frequency to find the global Love numbers, so the EOS must be solved first (:meth:`solve_eos`). The
        result stays on the world for the getters (:meth:`get_tidal_heating`, :meth:`get_tidal_potential_derivatives`,
        :meth:`get_layer_tidal_heating`) until the next solve, and becomes its tidal heat source
        (:attr:`tidal_heat_source`).

        Each orbital-state argument left as ``None`` comes from :meth:`get_tide_state`, the state the system the
        world belongs to gives it (its orbit about its tidal host, the host's mass, and its own spin and obliquity),
        so ``world.calc_tides()`` solves a world of a system in its current state, and
        ``world.calc_tides(eccentricity=0.05)`` changes one value.

        Parameters
        ----------
        orbital_frequency : float, optional
            Orbital mean motion [rad s-1].
        spin_frequency : float, optional
            Spin rate of the deformed body [rad s-1].
        eccentricity : float, optional
            Orbital eccentricity [dimensionless].
        obliquity : float, optional
            Axial tilt [radians].
        semi_major_axis : float, optional
            Orbital semi-major axis [m].
        host_mass : float, optional
            Mass of the tidal host [kg].

        Returns
        -------
        dict
            ``tidal_heating`` [W], the potential derivatives ``dU_dM``, ``dU_dw``, and ``dU_dO`` [J kg-1 rad-1] (see
            :meth:`get_tidal_potential_derivatives`), ``num_tidal_modes``, and ``layer_tidal_heating`` (each layer's
            heating [W] by name, as :meth:`get_layer_tidal_heating` gives it).

        Raises
        ------
        ValueError
            If an orbital-state argument is left out and the world has no tide state to take it from (it is not in a
            system, or has no tidal host), or the state is invalid (an eccentricity outside [0, 1), say).
        RuntimeError
            If no tide model is attached, the rheology model is selected but the EOS has
            not been solved, a radial-solver Love-number solve fails, or the global
            potential solve fails.
        """
        cdef c_TideSolveConfig state = cy_resolve_tide_state(
            self,
            "calc_tides",
            orbital_frequency,
            spin_frequency,
            eccentricity,
            obliquity,
            semi_major_axis,
            host_mass)
        with nogil:
            self._world_ptr.get().calc_tides(state)
        cdef dict result = {"tidal_heating": self._world_ptr.get().get_tidal_heating()}
        result.update(self.get_tidal_potential_derivatives())
        result["num_tidal_modes"] = self._world_ptr.get().get_num_tidal_modes()
        result["layer_tidal_heating"] = self.get_layer_tidal_heating()
        return result

    def get_layer_tidal_heating(self, layer=None):
        """Tidal heating [W] the last ``calc_tides`` put in a layer; NaN before one.

        With a radial-solver Love method it is the volume integral of the radial solution's heating density over the
        layer (NaN when the tides config's ``layer_tidal_heating`` is off); with the ``homogeneous``, ``cpl``, or
        ``ctl`` method, the heating of the layer's own scaled Love numbers; with an analytic tide model, the total
        times the layer's tidal scale.

        Parameters
        ----------
        layer : int or str, optional
            The layer's index (negative counts from the outermost) or name. Left out, every layer's heating.

        Returns
        -------
        float or dict
            The layer's heating [W], or with no ``layer`` a dict of each layer's heating by name, inner to outer.

        Raises
        ------
        IndexError, KeyError
            A layer index out of range, or a name the world does not have (see :meth:`get_layer`).
        """
        cdef c_BaseWorld* world_ptr = self._world_ptr.get()
        cdef list names
        cdef size_t layer_i
        if layer is None:
            names = cy_layer_names(world_ptr)
            return {names[layer_i]: world_ptr.get_layer_tidal_heating(layer_i) for layer_i in range(len(names))}
        return world_ptr.get_layer_tidal_heating(cy_layer_index(world_ptr, layer))

    @property
    def tidal_heat_source(self) -> dict:
        """The tidal heat source: each layer's heating [W] from the last ``calc_tides``.

        Every later :meth:`solve_eos` spreads it over the layers with ``use_heating``, along the radial profile of
        the per-layer heating integral with a radial-solver Love method, else by mass, and
        :meth:`calc_layer_temperature_rate` takes it as soon as ``calc_tides`` has run. So a step of an evolution is
        ``solve_eos``, ``calc_tides``, then the rates, with no second solve. A layer that is not tidal (``use_tides``
        off) takes none. A radial-solver tide with ``layer_tidal_heating`` off gives no per-layer split, so the total
        is spread over the tidal layers by mass. A change to the structure keeps the source; a failed ``calc_tides``,
        :meth:`clear_tidal_heating`, adding a layer, and loading a binary forget it. It is not saved with the world.

        Returns
        -------
        dict
            Heating [W] by layer name; layers ``calc_tides`` gave no heating are left out, and it is empty before
            a ``calc_tides``.
        """
        cdef vector[double] powers = self._world_ptr.get().get_tidal_heat_source()
        cdef list names = cy_layer_names(self._world_ptr.get())
        cdef size_t layer_i
        return {names[layer_i]: powers[layer_i] for layer_i in range(min(powers.size(), len(names)))
                if math.isfinite(powers[layer_i])}

    def clear_tidal_heating(self):
        """Forget the tidal heat source (see :attr:`tidal_heat_source`), so later solves heat no layer tidally."""
        self._world_ptr.get().clear_tidal_heating()

    def set_prescribed_heating(self, layer, power=None, specific_rate=None):
        """Prescribe a layer's internal heating: a ``power`` [W] spread over the layer by mass, or a
        ``specific_rate`` [W kg-1]. With neither, the layer's prescribed heating is cleared.

        It acts in a layer with ``use_heating``, through a :meth:`solve_eos` that carries temperature, beside the
        layer's radiogenics and the tidal heat source. The EOS solve reads it, so the world forgets its solved
        structure. It is saved with the world: :meth:`get_config_dict` writes it in a ``[prescribed_heating]`` table and
        the binary record holds it.

        Parameters
        ----------
        layer : str or int
            The layer's name or index.
        power, specific_rate : float, optional
            One of the two; ``None`` for the other.

        Raises
        ------
        ValueError
            For both given, or a layer the world does not have.
        """
        cdef size_t layer_index = cy_layer_index(self._world_ptr.get(), layer)
        cdef double power_value = math.nan if power is None else <double>power
        cdef double rate_value  = math.nan if specific_rate is None else <double>specific_rate
        self._world_ptr.get().set_prescribed_heating(layer_index, power_value, rate_value)

    @property
    def prescribed_heating(self) -> dict:
        """The prescribed heating by layer name, as ``{"power": watts}`` or ``{"specific_rate": watts_per_kg}``;
        layers without one are left out (see :meth:`set_prescribed_heating`)."""
        cdef list names = cy_layer_names(self._world_ptr.get())
        cdef c_PrescribedLayerHeating prescribed
        cdef size_t layer_i
        out = {}
        for layer_i in range(len(names)):
            prescribed = self._world_ptr.get().get_prescribed_heating(layer_i)
            if math.isfinite(prescribed.specific_rate):
                out[names[layer_i]] = {"specific_rate": prescribed.specific_rate}
            elif math.isfinite(prescribed.power):
                out[names[layer_i]] = {"power": prescribed.power}
        return out

    def calc_layer_temperature_rate(self, layer) -> float:
        """Rate of change of a layer's temperature [K s-1] from the heat entering, leaving, and generated in it.

        (C + C_latent) dT/dt = L_in - L_out + H, with the heat flows of the last :meth:`solve_eos`, H every heat
        source of a layer with ``use_heating`` (the last solve's radiogenic and prescribed heat, and the tidal heat
        of the latest ``calc_tides``; see :attr:`tidal_heat_source`), C the heat its profile stores per kelvin
        (:meth:`calc_layer_thermal_capacity`), and C_latent the latent heat its zone boundaries absorb per kelvin
        (:meth:`calc_layer_latent_capacity`). A step of an evolution is ``solve_eos``, ``calc_tides``, then this rate.

        Parameters
        ----------
        layer : str or int
            The layer's name or index.

        Returns
        -------
        float
            The rate [K s-1]; NaN before a solve or for a layer with no heat capacity.
        """
        cdef size_t layer_index = cy_layer_index(self._world_ptr.get(), layer)
        return self._world_ptr.get().calc_layer_temperature_rate(layer_index)

    def calc_continuation_frequency(self) -> float:
        """The frequency [rad s-1] below which this world's tidal modes are continued linearly to zero.

        A mode of frequency omega takes the Love number at hypot(omega, omega_c), omega_c this frequency, with its
        dissipation scaled by |omega| over that, so the tidal torque passes smoothly through a spin-orbit lock. It is
        ``[numerical] minimum_frequency``, or higher where most of the world's solid would be near-fluid. A solid of
        viscosity eta forced well below its Maxwell rate has a complex modulus near omega eta, and once that falls under
        a near-fluid threshold (1e-6 times rho g R) its Love solve is slow and its result solver noise. The viscosity
        used is the one below which the softer half of the solid volume lies, so a warm ice shell raises it and a thin
        partially molten band does not.

        Returns
        -------
        float
            The continuation frequency [rad s-1] of the last :meth:`solve_eos`; ``minimum_frequency`` before one and for
            an analytic tide model (``fixed_q``, ``fixed_dt``, ...), which has no Love solve.
        """
        return self._world_ptr.get().calc_continuation_frequency()

    def calc_layer_thermal_capacity(self, layer) -> float:
        """The heat a layer's profile stores per kelvin of its temperature [J K-1], from the last :meth:`solve_eos`.

        C = integral of rho c_p S 4 pi r^2 dr over the layer, with rho the solved density, c_p the material's effective
        heat capacity at the solved pressure and temperature (the latent heat of a melting range included), and S how
        far each point moves with the layer's temperature: one through an isothermal layer, T(r) / T along an
        adiabatic interior, and the share of a conducting stretch the layer's temperature moves. Computed on the first
        call after a solve (an adaptive quadrature over the solved profile) and cached until the next.

        Parameters
        ----------
        layer : str or int
            The layer's name or index.

        Returns
        -------
        float
            The capacity [J K-1]; NaN before a successful solve, zero for a layer without a material.

        Assumptions
        -----------
        - The neighbors' interface temperatures hold while the layer's own temperature changes, and a change of the
          layer's temperature moves its profile without moving its boundary layers.
        """
        cdef size_t layer_index = cy_layer_index(self._world_ptr.get(), layer)
        return self._world_ptr.get().calc_layer_thermal_capacity(layer_index)

    def calc_layer_latent_capacity(self, layer) -> float:
        """The latent heat [J K-1] the boundaries between a layer's solid and liquid zones absorb per kelvin of its
        temperature, from the last :meth:`solve_eos` (a Stefan condition).

        Only a material that melts at one temperature has such a boundary; one that melts over a range carries its
        latent heat in :meth:`calc_layer_thermal_capacity` instead. Computed and cached as that method is.

        Parameters
        ----------
        layer : str or int
            The layer's name or index.

        Returns
        -------
        float
            The capacity [J K-1]; zero for a layer without such a boundary, NaN before a successful solve.

        Assumptions
        -----------
        - The neighbors' interface temperatures hold while the layer's own temperature changes, and the melt that
          forms takes the solid's density at the boundary.
        """
        cdef size_t layer_index = cy_layer_index(self._world_ptr.get(), layer)
        return self._world_ptr.get().calc_layer_latent_capacity(layer_index)

    def get_heating(self, radius):
        """Volumetric heating the last solve's heat sources give at a radius, every source summed (radiogenic,
        tidal, and prescribed).

        A solve with ``solve_temperature`` off carries no heat flow, so the heating reported then did not act on it.

        Parameters
        ----------
        radius : float or np.ndarray
            Radius [m].

        Returns
        -------
        float or np.ndarray
            Heating [W m-3]: zero in a layer without ``use_heating``, NaN before a successful solve or outside the
            world.
        """
        cdef double[::1] radius_view
        cdef double[::1] out_view
        cdef Py_ssize_t radius_i
        cdef Py_ssize_t num_radii
        cdef c_BaseWorld* world_ptr = self._world_ptr.get()
        if np.ndim(radius) == 0:
            return world_ptr.get_heating(<double>radius)
        radius_array = np.ascontiguousarray(radius, dtype=np.float64)
        out = np.empty_like(radius_array)
        radius_view = radius_array.reshape(-1)
        out_view = out.reshape(-1)
        num_radii = radius_view.shape[0]
        with nogil:
            for radius_i in range(num_radii):
                out_view[radius_i] = world_ptr.get_heating(radius_view[radius_i])
        return out

    def get_layer_tidal_scale(self, layer) -> float:
        """The tidal scale a layer carries in the quasi-homogeneous Love methods [dimensionless].

        The layer's configured ``tidal_scale``, or its volume over the planet's when none is set; 0 for a layer that
        is not tidal. ``layer`` is its index (negative counts from the outermost) or name; ``IndexError`` or
        ``KeyError`` for one the world does not have (see :meth:`get_layer`).
        """
        return self._world_ptr.get().get_layer_tidal_scale(cy_layer_index(self._world_ptr.get(), layer))

    def get_3d_tidal_heating(
            self,
            orbital_frequency=None,
            spin_frequency=None,
            eccentricity=None,
            obliquity=None,
            semi_major_axis=None,
            host_mass=None,
            radius=None,
            colatitude=None) -> float:
        """Secular (cycle and orbit-averaged) 3D tidal volumetric heating [W m-3] at ``(radius, colatitude)``.

        This is the longitude mean of the time-averaged power density. The active tidal modes come from the
        world's ``[tides]`` truncation config and are merged into coherent waves (modes sharing a real spatial
        function, such as the ``m = 0`` pair at ``+omega`` and ``-omega``, form one wave); the world radial
        response is solved once per ``(l, |omega|)``, and each frequency contributes
        ``(|omega|/2) Im(sigma_c : conj(eps_c))`` of its total complex stress and strain. Cross terms between
        waves at one frequency with different longitude structure make the heating of a synchronously
        rotating body depend on longitude; they average to zero over longitude, so this method returns the
        zonal mean and :meth:`calc_3d_tides` gives the longitude-resolved field. The volume integral over the
        planet equals the 1D :meth:`get_tidal_heating`. Requires the rheology tide model
        (:meth:`set_tide_model`) and a solved EOS (:meth:`solve_eos`). Returns NaN at the center and below the
        solver's starting radius, and 0 in liquid layers.

        Each call solves every radial response again, which costs about as much as the whole of
        :meth:`get_3d_tidal_heating_array` over hundreds of points; for more than one point, use that.

        The six orbital-state arguments are those of :meth:`calc_tides`, each ``None`` taken from
        :meth:`get_tide_state`; ``radius`` [m] and ``colatitude`` [rad] are required (by keyword when the orbital
        state is left out).
        """
        cy_require_arguments("get_3d_tidal_heating", ("radius", "colatitude"), (radius, colatitude))
        cdef double point_radius = <double>radius
        cdef double point_colatitude = <double>colatitude
        cdef c_TideSolveConfig state = cy_resolve_tide_state(
            self,
            "get_3d_tidal_heating",
            orbital_frequency,
            spin_frequency,
            eccentricity,
            obliquity,
            semi_major_axis,
            host_mass)
        cdef double heating
        with nogil:
            heating = self._world_ptr.get().get_3d_tidal_heating(state, point_radius, point_colatitude)
        return heating

    def get_3d_tidal_heating_array(
            self,
            orbital_frequency=None,
            spin_frequency=None,
            eccentricity=None,
            obliquity=None,
            semi_major_axis=None,
            host_mass=None,
            radii=None,
            colatitudes=None,
            int num_threads=0):
        """Longitude-mean secular 3D tidal volumetric heating [W m-3] at ``(radius, colatitude)`` points.

        Batch form of :meth:`get_3d_tidal_heating`: ``radii`` [m] and ``colatitudes`` [rad] are paired, equal-length
        1D arrays (point ``i`` is ``(radii[i], colatitudes[i])``), or one of them is a scalar that pairs with every
        value of the other, and a 1D ``np.ndarray`` of heating, one value per point, is returned. Same physics and
        preconditions as the scalar method, but the world radial response is solved once per unique
        ``(degree l, |omega|)`` and reused across all points, so this is the efficient way to build a zonal-mean
        heating map. The six orbital-state arguments are those of :meth:`calc_tides`, each ``None`` taken from
        :meth:`get_tide_state`.

        ``num_threads`` spreads the radial solves (one per degree and frequency, from ``[numerical]``
        ``love_solve_min_parallel`` of them on) and then the per-point evaluation over that many threads; the result
        is identical for any thread count. The default, 0, uses the logical processors less 4 (at least 1). Pass 1
        inside a process or thread pool that already occupies the machine.
        """
        cy_check_num_threads(num_threads)
        cy_require_arguments("get_3d_tidal_heating_array", ("radii", "colatitudes"), (radii, colatitudes))
        cdef cnp.ndarray radii_arr = np.atleast_1d(np.asarray(radii, dtype=np.float64)).ravel()
        cdef cnp.ndarray colat_arr = np.atleast_1d(np.asarray(colatitudes, dtype=np.float64)).ravel()
        # A scalar (one value) pairs with every point of the other axis.
        if (radii_arr.shape[0] == 1) and (np.ndim(radii) == 0):
            radii_arr = np.full(colat_arr.shape[0], radii_arr[0])
        if (colat_arr.shape[0] == 1) and (np.ndim(colatitudes) == 0):
            colat_arr = np.full(radii_arr.shape[0], colat_arr[0])
        radii_arr = np.ascontiguousarray(radii_arr)
        colat_arr = np.ascontiguousarray(colat_arr)
        if radii_arr.shape[0] != colat_arr.shape[0]:
            raise ValueError("radii and colatitudes must have the same length")

        cdef size_t num_points = <size_t>radii_arr.shape[0]
        cdef cnp.ndarray out_arr = np.empty(num_points, dtype=np.float64)
        if num_points == 0:
            return out_arr

        cdef double[::1] radii_view = radii_arr
        cdef double[::1] colat_view = colat_arr
        cdef double[::1] out_view   = out_arr

        cdef c_TideSolveConfig state = cy_resolve_tide_state(
            self,
            "get_3d_tidal_heating_array",
            orbital_frequency,
            spin_frequency,
            eccentricity,
            obliquity,
            semi_major_axis,
            host_mass)

        with nogil:
            self._world_ptr.get().get_3d_tidal_heating_array(
                state,
                &radii_view[0],
                &colat_view[0],
                num_points,
                &out_view[0],
                num_threads)
        return out_arr

    def calc_3d_displacements(
            self,
            orbital_frequency=None,
            spin_frequency=None,
            eccentricity=None,
            obliquity=None,
            semi_major_axis=None,
            host_mass=None,
            radii=None,
            colatitudes=None,
            longitudes=None,
            times=None,
            int num_threads=0) -> dict:
        """Instantaneous tidal displacements [m] on the grid ``(radius, colatitude, longitude, time)``.

        The active tidal modes are built from the world's ``[tides]`` truncation config, the world radial
        response is solved once per unique ``(l, frequency)``, and each mode's complex displacement
        amplitude ``(y1 U, y3 dU/dtheta, y3 dU/dphi / sin theta)`` at a point is evolved in time as
        ``Re[u e^{i omega t}]`` and summed over the modes (the same phasor convention as the instantaneous
        heating of :meth:`calc_3d_tides`). Requires the rheology tide model and a solved EOS, and a
        radial-solver Love-number method (the analytic methods have no radial functions).

        Parameters
        ----------
        orbital_frequency, spin_frequency, eccentricity, obliquity, semi_major_axis, host_mass : float, optional
            The orbital/spin state, as for :meth:`calc_tides`: each ``None`` comes from :meth:`get_tide_state`.
        radii, colatitudes, longitudes, times : array-like of float
            Grid axes [m], [rad], [rad], [s]; scalars are accepted. Required (by keyword when the orbital state is
            left out).
        num_threads : int, optional
            Threads for the radial solves (one per degree and frequency) and the per-point evaluation.
            Default 0, the logical processors less 4 (at least 1); pass 1 inside a process or thread pool that
            already occupies the machine. The result is identical for any thread count.

        Returns
        -------
        dict
            ``radii``, ``colatitudes``, ``longitudes``, ``times`` (the axes as 1-D arrays) and ``radial``,
            ``polar``, ``azimuthal``: float64 arrays of shape ``(nr, ncolat, nlon, ntime)`` [m]. A radius
            with no depth-resolved solution (the center, below the solver start) is NaN.

        Assumptions
        -----------
        - Linear superposition of the tidal modes; displacements follow the radial functions y1 (radial)
          and y3 (tangential) of each mode's radial solution.
        """
        cy_check_num_threads(num_threads)
        cy_require_arguments(
            "calc_3d_displacements",
            ("radii", "colatitudes", "longitudes", "times"),
            (radii, colatitudes, longitudes, times))
        cdef c_Grid3DAxes axes
        radii_arr, colat_arr, lon_arr, time_arr = cy_grid_axes(radii, colatitudes, longitudes, times, &axes)
        cdef cnp.ndarray out_arr = np.empty(
            (axes.num_radii, axes.num_colatitudes, axes.num_longitudes, axes.num_times, 3), dtype=np.float64)
        cdef double[:, :, :, :, ::1] out_view = out_arr
        cdef c_TideSolveConfig state = cy_resolve_tide_state(
            self,
            "calc_3d_displacements",
            orbital_frequency,
            spin_frequency,
            eccentricity,
            obliquity,
            semi_major_axis,
            host_mass)
        with nogil:
            self._world_ptr.get().get_3d_displacements_grid(
                state,
                axes,
                &out_view[0, 0, 0, 0, 0],
                num_threads)
        return {
            "radii": radii_arr,
            "colatitudes": colat_arr,
            "longitudes": lon_arr,
            "times": time_arr,
            "radial": np.ascontiguousarray(out_arr[..., 0]),
            "polar": np.ascontiguousarray(out_arr[..., 1]),
            "azimuthal": np.ascontiguousarray(out_arr[..., 2]),
        }

    def calc_3d_stress_strain(
            self,
            orbital_frequency=None,
            spin_frequency=None,
            eccentricity=None,
            obliquity=None,
            semi_major_axis=None,
            host_mass=None,
            radii=None,
            colatitudes=None,
            longitudes=None,
            times=None,
            cpp_bool return_stress=True,
            cpp_bool return_strain=True,
            int num_threads=0) -> dict:
        """Instantaneous tidal stress [Pa] and strain on the grid ``(radius, colatitude, longitude, time)``.

        The active tidal modes come from the world's ``[tides]`` truncation config and are merged into
        coherent waves; the radial response is solved once per unique ``(l, frequency)``, and at each point
        every wave's complex stress and strain amplitude is added into the total of its frequency. Each
        component at time ``t`` is the sum over frequencies of ``Re[amplitude e^{i |frequency| t}]``, the
        convention of :meth:`calc_3d_displacements` and of the instantaneous heating of
        :meth:`calc_3d_tides`.

        Parameters
        ----------
        orbital_frequency, spin_frequency, eccentricity, obliquity, semi_major_axis, host_mass : float, optional
            The orbital and spin state, as for :meth:`calc_tides`: each ``None`` comes from :meth:`get_tide_state`.
        radii, colatitudes, longitudes, times : array-like of float
            Grid axes [m], [rad], [rad], [s]; scalars are accepted. Required (by keyword when the orbital state is
            left out).
        return_stress, return_strain : bool, optional
            Which tensors to compute; each takes 48 bytes per grid point and time. Default both.
        num_threads : int, optional
            Threads for the radial solves (one per degree and frequency) and the per-point evaluation.
            Default 0, the logical processors less 4 (at least 1); pass 1 inside a process or thread pool that
            already occupies the machine. The result is identical for any thread count.

        Returns
        -------
        dict
            ``radii``, ``colatitudes``, ``longitudes``, ``times`` (the axes as 1-D arrays), ``components`` (the
            component names in order: ``rr``, ``theta_theta``, ``phi_phi``, ``r_theta``, ``r_phi``,
            ``theta_phi``), and the requested ``stress`` [Pa] and ``strain``: float64 arrays of shape
            ``(nr, ncolat, nlon, ntime, 6)``.

        Raises
        ------
        ValueError
            If an axis is empty, neither tensor is requested, or ``num_threads`` is negative.
        RuntimeError
            If the world has no rheology tide model or no solved EOS, its Love method has no radial functions, or
            a radial solve fails.

        Assumptions
        -----------
        - Linear superposition of the tidal modes and isotropic linear viscoelasticity with the complex moduli at
          each mode's frequency; the strain is the symmetric gradient of :meth:`calc_3d_displacements`.
        - Solid layers only: a point in a liquid layer, or at a radius with no depth-resolved solution (the center,
          below the solver start), is NaN.
        - Only modes with a nonzero forcing frequency are included; the permanent tide is not.
        """
        if not (return_stress or return_strain):
            raise ValueError("At least one of return_stress and return_strain must be True.")
        cy_check_num_threads(num_threads)
        cy_require_arguments(
            "calc_3d_stress_strain",
            ("radii", "colatitudes", "longitudes", "times"),
            (radii, colatitudes, longitudes, times))
        cdef c_Grid3DAxes axes
        radii_arr, colat_arr, lon_arr, time_arr = cy_grid_axes(radii, colatitudes, longitudes, times, &axes)
        cdef tuple tensor_shape = (axes.num_radii, axes.num_colatitudes, axes.num_longitudes, axes.num_times, 6)

        # Output buffers are owned by numpy; a skipped tensor passes a null pointer.
        cdef cnp.ndarray stress_arr = None
        cdef cnp.ndarray strain_arr = None
        cdef double[:, :, :, :, ::1] stress_view
        cdef double[:, :, :, :, ::1] strain_view
        cdef double* stress_ptr = NULL
        cdef double* strain_ptr = NULL
        if return_stress:
            stress_arr = np.empty(tensor_shape, dtype=np.float64)
            stress_view = stress_arr
            stress_ptr = &stress_view[0, 0, 0, 0, 0]
        if return_strain:
            strain_arr = np.empty(tensor_shape, dtype=np.float64)
            strain_view = strain_arr
            strain_ptr = &strain_view[0, 0, 0, 0, 0]

        cdef c_TideSolveConfig state = cy_resolve_tide_state(
            self,
            "calc_3d_stress_strain",
            orbital_frequency,
            spin_frequency,
            eccentricity,
            obliquity,
            semi_major_axis,
            host_mass)
        with nogil:
            self._world_ptr.get().get_3d_stress_strain_grid(
                state,
                axes,
                stress_ptr,
                strain_ptr,
                num_threads)

        cdef dict out = {
            "radii": radii_arr,
            "colatitudes": colat_arr,
            "longitudes": lon_arr,
            "times": time_arr,
            "components": STRESS_STRAIN_COMPONENTS,
        }
        if return_stress:
            out["stress"] = stress_arr
        if return_strain:
            out["strain"] = strain_arr
        return out

    def calc_3d_tides(
            self,
            orbital_frequency=None,
            spin_frequency=None,
            eccentricity=None,
            obliquity=None,
            semi_major_axis=None,
            host_mass=None,
            radii=None,
            colatitudes=None,
            longitudes=None,
            times=None,
            orbit_averaged=True,
            latitude_summed=False,
            longitude_summed=False,
            radial_summed=False,
            latitude_nodes=None,
            longitude_nodes=None,
            radial_slices=None,
            latitude_analytic=True,
            double colatitude_min=0.0,
            double colatitude_max=np.pi,
            int num_threads=0) -> dict:
        """3D tidal heating as a grid over ``(radius, colatitude, longitude[, time])``, optionally reduced.

        With ``orbit_averaged=True`` (default) the quantity is the secular volumetric heating density
        ``h_bar`` [W m-3], the time average of the instantaneous power at each point. It has no time axis and
        depends on longitude wherever waves at one frequency have different longitude structure, as they do
        for a synchronously rotating body. With ``orbit_averaged=False`` it is the instantaneous mechanical
        power density ``sigma_ij(t) * eps_dot_ij(t)`` [W m-3] at each supplied time (a fourth axis), which
        time-averages to ``h_bar`` through e^N: ``h_bar`` cuts every product of two eccentricity functions at the
        truncation level's e^N, while the instantaneous power keeps the partial terms past it.

        Any spatial dimension can be integrated out with ``latitude_summed``, ``longitude_summed``, or
        ``radial_summed``. If any spatial axis is summed, the surviving spatial axes carry their Jacobian
        (``r^2``, ``sin theta``, ``1``) so a plain integral over them recovers the total; if none is summed
        the output is the raw density. The colatitude integral uses an internal Gauss-Legendre grid
        (``latitude_nodes``), the radial integral ``radial_slices`` Gauss-Legendre nodes inside each layer
        (none on a layer boundary), and the longitude integral the analytic ``2*pi`` times the longitude mean
        when averaged or a ``longitude_nodes`` trapezoid when instantaneous. Each of the three left as ``None``
        takes the ``tides_3d_*`` value of the ``[numerical]`` section of the TidalPy configuration (16, 64, and
        16 by default).

        Non-summed spatial axes require the matching ``radii``, ``colatitudes``, or ``longitudes`` array, and
        ``times`` is required when ``orbit_averaged=False``. The returned dict carries the surviving axes plus
        either ``heating`` (ordered radius, colatitude, longitude, time) or, when all three spatial axes are
        summed, ``total`` [W] and ``per_layer`` [W] (innermost first, each an array over time when
        instantaneous). Requires the rheology tide model and a solved EOS. The six orbital-state arguments are those of
        :meth:`calc_tides`, each ``None`` taken from :meth:`get_tide_state`.

        With ``latitude_summed``, ``longitude_summed``, and ``orbit_averaged`` the colatitude integral uses the
        precomputed analytic angular Gram table of the longitude mean (exact, no theta grid);
        ``latitude_analytic=False`` falls back to the Gauss-Legendre quadrature, which agrees to quadrature accuracy,
        and a call that keeps its longitudes always uses the quadrature. A point in a liquid (a liquid layer or a liquid
        zone, see :attr:`zones`) has zero heating. A latitude band can be integrated instead of the full sphere by
        setting ``colatitude_min`` and ``colatitude_max`` [rad] (defaults 0 and pi), so complementary bands add up to
        the full-sphere result; a band
        narrower than the full sphere always uses the quadrature. The band has no effect when colatitude is not summed.

        ``num_threads`` spreads the radial solves (one per degree and frequency, from ``[numerical]``
        ``love_solve_min_parallel`` of them on), then the per-point evaluation over colatitude rows, or for the
        analytic colatitude collapse the radii; the result is identical for any thread count. The default, 0, uses
        the logical processors less 4 (at least 1); pass 1 inside a process or thread pool that already occupies the
        machine.
        """
        cy_check_num_threads(num_threads)
        if not (0.0 <= colatitude_min < colatitude_max <= np.pi + 1.0e-12):
            raise ValueError("colatitude band must satisfy 0 <= colatitude_min < colatitude_max <= pi")
        cdef cpp_bool instantaneous = not orbit_averaged
        cdef c_Heating3DCollapseConfig cfg
        cfg.orbit_averaged   = <cpp_bool>orbit_averaged
        cfg.latitude_summed  = <cpp_bool>latitude_summed
        cfg.longitude_summed = <cpp_bool>longitude_summed
        cfg.radial_summed    = <cpp_bool>radial_summed
        if latitude_nodes is not None:
            cfg.latitude_nodes = <int>int(latitude_nodes)
        if longitude_nodes is not None:
            cfg.longitude_nodes = <int>int(longitude_nodes)
        if radial_slices is not None:
            cfg.radial_slices = <int>int(radial_slices)
        cfg.num_threads      = num_threads
        cfg.latitude_analytic = <cpp_bool>latitude_analytic
        cfg.colatitude_min   = colatitude_min
        # The default (np.pi) can sit a rounding step above the C++ d_PI; clamp so the
        # full-sphere default never trips the C++ band check.
        cfg.colatitude_max   = colatitude_max
        if cfg.colatitude_max > d_PI:
            cfg.colatitude_max = d_PI

        # A summed axis is not supplied: it passes a null pointer and no values. The arrays stay alive to the end.
        cdef const double* radii_ptr = NULL
        cdef const double* colat_ptr = NULL
        cdef const double* lon_ptr = NULL
        cdef const double* time_ptr = NULL
        cdef size_t num_radii = 0
        cdef size_t num_colat = 0
        cdef size_t num_lon = 0
        cdef size_t num_time = 0
        if not radial_summed:
            if radii is None:
                raise ValueError("radii must be provided when radial_summed is False")
            radii_arr = cy_grid_axis(radii, &radii_ptr, &num_radii)
        if not latitude_summed:
            if colatitudes is None:
                raise ValueError("colatitudes must be provided when latitude_summed is False")
            colat_arr = cy_grid_axis(colatitudes, &colat_ptr, &num_colat)
        if not longitude_summed:
            if longitudes is None:
                raise ValueError("longitudes must be provided when longitude_summed is False")
            lon_arr = cy_grid_axis(longitudes, &lon_ptr, &num_lon)
        if instantaneous:
            if times is None:
                raise ValueError("times must be provided when orbit_averaged is False")
            time_arr = cy_grid_axis(times, &time_ptr, &num_time)

        cdef c_TideSolveConfig state = cy_resolve_tide_state(
            self,
            "calc_3d_tides",
            orbital_frequency,
            spin_frequency,
            eccentricity,
            obliquity,
            semi_major_axis,
            host_mass)

        # The layout gives the output shape, so the heating can be written straight into numpy-owned buffers.
        cdef c_Heating3DCollapsed layout
        with nogil:
            layout = self._world_ptr.get().calc_3d_tides_layout(
                radii_ptr,
                num_radii,
                colat_ptr,
                num_colat,
                lon_ptr,
                num_lon,
                time_ptr,
                num_time,
                cfg)

        # Surviving-axis shape (row-major, order radius, colatitude, longitude, time).
        cdef list shape = [<Py_ssize_t>layout.shape[i] for i in range(layout.shape.size())]
        cdef Py_ssize_t num_values = 1
        cdef Py_ssize_t axis_length
        for axis_length in shape:
            num_values *= axis_length
        cdef Py_ssize_t nlayers = <Py_ssize_t>layout.n_layers
        cdef Py_ssize_t ntimes = <Py_ssize_t>layout.n_times
        cdef Py_ssize_t num_layer_totals = nlayers * ntimes if layout.all_spatial_summed else 0
        cdef cnp.ndarray values_arr = np.empty(num_values, dtype=np.float64)
        cdef cnp.ndarray layer_totals_arr = np.empty(num_layer_totals, dtype=np.float64)
        cdef double[::1] values_view
        cdef double[::1] layer_totals_view
        cdef double* values_ptr = NULL
        cdef double* layer_totals_ptr = NULL
        if num_values > 0:
            values_view = values_arr
            values_ptr = &values_view[0]
        if num_layer_totals > 0:
            layer_totals_view = layer_totals_arr
            layer_totals_ptr = &layer_totals_view[0]

        with nogil:
            self._world_ptr.get().calc_3d_tides_into(
                state,
                radii_ptr,
                num_radii,
                colat_ptr,
                num_colat,
                lon_ptr,
                num_lon,
                time_ptr,
                num_time,
                cfg,
                values_ptr,
                layer_totals_ptr)

        cdef dict out = {'radii': cy_vec_to_ndarray(layout.radii),
                         'colatitudes': cy_vec_to_ndarray(layout.colatitudes),
                         'longitudes': cy_vec_to_ndarray(layout.longitudes)}
        if instantaneous:
            out['times'] = cy_vec_to_ndarray(layout.times)

        cdef cnp.ndarray per_layer
        if layout.all_spatial_summed:
            per_layer = layer_totals_arr.reshape(nlayers, ntimes)
            if instantaneous:
                out['total'] = values_arr
                out['per_layer'] = per_layer
            else:
                out['total'] = float(values_arr[0]) if num_values else 0.0
                out['per_layer'] = per_layer[:, 0]
        else:
            out['heating'] = values_arr.reshape(shape) if shape else values_arr
        return out

    def set_solver_defaults(self, eos_solver=None, radial_solver=None):
        """Pin ``[eos_solver]`` and ``[radial_solver]`` settings on this world.

        A world file may carry the two tables so that, with the TidalPy configuration, it reproduces its run on
        another machine. A pinned key replaces the configuration's value for every solve this world runs
        (:meth:`solve_eos`, :meth:`solve_love_numbers`, :meth:`calc_tides`, and the 3D paths); a call's own
        argument still wins over it, and a key left out keeps following the configuration, so a configuration
        changed after the world was built still reaches it. The tables come back from :meth:`get_solver_defaults`
        and :meth:`get_config_dict`.

        Parameters
        ----------
        eos_solver, radial_solver : dict, optional
            Keys of the matching section of ``TidalPy_Configs.toml`` (``EOS_SOLVER_KEYS`` and
            ``RADIAL_SOLVER_KEYS`` of ``TidalPy.Structures.configs``) with their values. A table given replaces
            the one stored, so an empty dict clears it; a table left as ``None`` is untouched.

        Raises
        ------
        ValueError
            For a key the section does not have, a value of the wrong type or out of range, or an unknown
            integration method.
        """
        from TidalPy.Structures.configs.toml_loader import validate_solver_table
        cdef c_EOSSolverOverrides eos
        cdef c_RadialSolverOverrides radial
        cdef str key
        if eos_solver is not None:
            validate_solver_table("eos_solver", eos_solver, "set_solver_defaults")
            for key, value in eos_solver.items():
                eos.set(key.encode("utf-8"), cy_solver_setting_to_double(eos.kind(key.encode("utf-8")), value))
            self._world_ptr.get().set_eos_solver_overrides(eos)
        if radial_solver is not None:
            validate_solver_table("radial_solver", radial_solver, "set_solver_defaults")
            for key, value in radial_solver.items():
                radial.set(key.encode("utf-8"), cy_solver_setting_to_double(radial.kind(key.encode("utf-8")), value))
            self._world_ptr.get().set_radial_solver_overrides(radial)

    def get_solver_defaults(self) -> dict:
        """The ``[eos_solver]`` and ``[radial_solver]`` keys pinned on this world, under the configuration's names.

        Returns
        -------
        dict
            ``eos_solver`` and ``radial_solver`` tables, each present only when it pins a key. Empty when the
            world follows the TidalPy configuration throughout.
        """
        cdef c_EOSSolverOverrides eos = self._world_ptr.get().get_eos_solver_overrides()
        cdef c_RadialSolverOverrides radial = self._world_ptr.get().get_radial_solver_overrides()
        cdef dict out = {}
        cdef dict table = {}
        cdef string key
        for key in c_EOSSolverOverrides.keys():
            if eos.has(key):
                table[key.decode("utf-8")] = cy_solver_setting_from_double(c_EOSSolverOverrides.kind(key), eos.get(key))
        if table:
            out["eos_solver"] = table
        table = {}
        for key in c_RadialSolverOverrides.keys():
            if radial.has(key):
                table[key.decode("utf-8")] = cy_solver_setting_from_double(
                    c_RadialSolverOverrides.kind(key), radial.get(key))
        if table:
            out["radial_solver"] = table
        return out

    cpdef dict get_config_dict(self):
        """Return the world configuration as the TOML builder's world table (MKS).

        The dict validates against the world schema and rebuilds the same world through ``build_world``. It carries
        a ``tides`` table when a tide model is attached (``global_tidal_model`` plus the model's per-degree
        parameters and the stored degree and truncation settings), and a ``layers`` table keyed by layer name when
        the world has layers. Each layer entry is the layer's own ``get_config_dict`` (scalars, the material table,
        attached-model sub-tables) minus the standalone-only keys the builder derives itself (``name``,
        ``radius_inner_m``).

        Returns
        -------
        dict
            Keys: ``schema_version``, ``name``, ``type``, ``radius_m``, ``mass_kg``, ``albedo``, ``emissivity``,
            ``obliquity_rad``, ``spin_frequency_rad_s``, ``moment_of_inertia_factor`` (the attached spin model's, when
            the source configuration gave it or it differs from the ``[worlds]`` default for the world type; the
            solved structure's factor is :attr:`moment_of_inertia_factor`), ``tides`` when set, the ``eos_solver``
            and ``radial_solver`` tables when the world pins any solver key (:meth:`set_solver_defaults`), ``layers``
            when there are any, and ``prescribed_heating`` (by layer name,
            ``power_w`` or ``specific_rate_w_kg``) when a layer has prescribed heating.

        Raises
        ------
        ValueError
            If two layers share a name (the table needs unique keys).
        """
        from TidalPy.Structures.configs.toml_loader import SCHEMA_VERSION, world_type_defaults
        cdef c_BaseWorld* p = self._world_ptr.get()
        cdef dict config = {
            "schema_version":       SCHEMA_VERSION,
            "name":                 p.get_name().decode("utf-8"),
            "type":                 self.get_builder_world_type(),
            "radius_m":             p.get_radius(),
            "mass_kg":              p.get_mass(),
            "albedo":               p.get_albedo(),
            "emissivity":           p.get_emissivity(),
            "obliquity_rad":        p.get_obliquity(),
            "spin_frequency_rad_s": p.get_spin_frequency(),
        }
        cdef dict tides
        if p.get_tide_model_set():
            tides = cy_physics_model_config(<const c_PhysicsBase*>p.get_tide_model())
            tides["global_tidal_model"] = tides.pop("model")
            tides.update(self.get_tide_config())
            config["tides"] = tides
        # The spin model's factor is written when the source configuration gave it or it differs from the [worlds]
        # default this world type would take again; the default itself is the configuration's to supply.
        cdef double moment_of_inertia_factor = p.get_spin_model().get_config().moment_of_inertia_factor
        if (("moment_of_inertia_factor" in (self.source_config or {}))
                or (moment_of_inertia_factor != world_type_defaults(config["type"]).get("moment_of_inertia_factor"))):
            config["moment_of_inertia_factor"] = moment_of_inertia_factor
        config.update(self.get_solver_defaults())
        cdef dict layers = {}
        cdef dict layer_config
        cdef str layer_name
        cdef object view
        for view in self._ensure_layer_views():
            layer_config = view.get_config_dict()
            layer_name = layer_config.pop("name")
            for key in LAYER_STANDALONE_CONFIG_KEYS:
                layer_config.pop(key, None)
            if layer_name in layers:
                raise ValueError(
                    f"Layer name '{layer_name}' is not unique; the world config table needs unique layer names.")
            layers[layer_name] = layer_config
        if layers:
            config["layers"] = layers
        cdef dict prescribed = self.prescribed_heating
        if prescribed:
            config["prescribed_heating"] = {
                name: ({"power_w": heating["power"]} if "power" in heating
                       else {"specific_rate_w_kg": heating["specific_rate"]})
                for name, heating in prescribed.items()}
        return config
