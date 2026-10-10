# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrappers for TidalPy's system class.

A ``System`` links two or more worlds. Each world names its own tidal host (or none) and has a
two-body orbit about it described by a semi-major axis [m] and an eccentricity. The system
owns its worlds through shared pointers, co-owned with the Python world wrappers, so the added
world objects stay usable and are handed straight back by iteration (``for world in system``),
indexing (``system[i]``), and attribute access (``system.<world_name>``).
"""

import copy
import operator
import os
from collections.abc import Mapping
from numbers import Integral

import TidalPy

import numpy as np

from libc.math cimport NAN, isfinite
from libc.string cimport memcpy
from libcpp cimport bool as cpp_bool
from libcpp.vector cimport vector
from libcpp.memory cimport make_unique

from TidalPy.Utilities.logging.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities.classes.classes cimport cy_existing_binary_path
from TidalPy.Structures.worlds.base cimport BaseWorld, cy_set_world_configs, cy_world_configs
from TidalPy.Structures.worlds.terrestrial cimport TerrestrialWorld
from TidalPy.Structures.worlds.gasgiant cimport GasGiantWorld
from TidalPy.Structures.worlds.stellar cimport StarWorld

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())

# Pull in the out-of-line (inline) definition of c_BaseWorld::calc_tides so the system can drive a world's tides
# directly in C++. Its rheology path runs the CyRK-backed radial solver (linked via the CyRK cimport in system.pxd).
cdef extern from "world_tides_.hpp" nogil:
    pass


cdef BaseWorld cy_wrap_world(shared_ptr[c_BaseWorld] ptr):
    """Wrap a C++ world (e.g. one loaded by c_System::read_binary) as the matching Python wrapper.

    Dispatches on the world's concrete type so a terrestrial / gas-giant / star world comes back as its own
    wrapper class (a star with its stellar methods), not a bare BaseWorld.
    """
    cdef int kind = c_world_kind(ptr.get())
    if kind == 1:
        return TerrestrialWorld._wrap(ptr)
    if kind == 2:
        return GasGiantWorld._wrap(ptr)
    if kind == 3:
        return StarWorld._wrap(ptr)
    return BaseWorld._wrap(ptr)


# The name of a member world, by index.
cdef str cy_world_name(c_System* system_ptr, size_t index):
    return system_ptr.get_world(index).get().get_name().decode("utf-8")


cdef dict cy_dissipation_to_dict(c_TidalDissipation dissipation, c_System* system_ptr):
    """Convert a c_TidalDissipation result into a plain Python dict (all values MKS)."""
    return {
        'world_index':       <int>dissipation.world_index,
        'world_name':        cy_world_name(system_ptr, dissipation.world_index),
        'companion_name':    (cy_world_name(system_ptr, dissipation.companion_index)
                              if dissipation.solved else None),
        'solved':            True if dissipation.solved else False,
        'has_tide_model':    True if dissipation.has_tide_model else False,
        'orbital_frequency': dissipation.orbital_frequency,
        'semi_major_axis':   dissipation.semi_major_axis,
        'eccentricity':      dissipation.eccentricity,
        'spin_frequency':    dissipation.spin_frequency,
        'obliquity':         dissipation.obliquity,
        'companion_mass':    dissipation.companion_mass,
        'target_mass':       dissipation.target_mass,
        'tidal_heating':     dissipation.tidal_heating,
        'dU_dM':             dissipation.dU_dM,
        'dU_dw':             dissipation.dU_dw,
        'dU_dO':             dissipation.dU_dO,
        'dU_dM_minus_dw':    dissipation.dU_dM_minus_dw,
        'moment_of_inertia': dissipation.moment_of_inertia,
    }


cdef dict cy_evolution_to_dict(c_WorldEvolution evolution, c_System* system_ptr):
    """Convert a c_WorldEvolution result into a plain Python dict (all values MKS)."""
    return {
        'world_index':       <int>evolution.world_index,
        'world_name':        cy_world_name(system_ptr, evolution.world_index),
        'evolved':           True if evolution.evolved else False,
        'has_tide_model':    True if evolution.has_tide_model else False,
        'orbital_frequency': evolution.orbital_frequency,
        'semi_major_axis':   evolution.semi_major_axis,
        'eccentricity':      evolution.eccentricity,
        'spin_frequency':    evolution.spin_frequency,
        'host_mass':         evolution.host_mass,
        'target_mass':       evolution.target_mass,
        'tidal_heating':     evolution.tidal_heating,
        'dU_dM':             evolution.dU_dM,
        'dU_dw':             evolution.dU_dw,
        'dU_dO':             evolution.dU_dO,
        'dU_dM_minus_dw': evolution.dU_dM_minus_dw,
        'da_dt':             evolution.da_dt,
        'de_dt':             evolution.de_dt,
        'dn_dt':             evolution.dn_dt,
        'dspin_dt':          evolution.dspin_dt,
        'moment_of_inertia': evolution.moment_of_inertia,
        'has_spin':          True if evolution.has_spin else False,
        'dE_orbit_dt':       evolution.dE_orbit_dt,
        'dE_spin_dt':        evolution.dE_spin_dt,
        'energy_residual':   evolution.energy_residual,
    }


cdef dict cy_pair_to_dict(c_PairEvolution pair, c_System* system_ptr):
    """Convert a c_PairEvolution (dual-body) result into a plain Python dict (all values MKS), each body's part keyed
    by its world name."""
    # A world asked for with no tidal host has no partner (the pair names it twice) and no parts.
    cdef cpp_bool has_partner = pair.first_index != pair.second_index
    cdef str first_name = cy_world_name(system_ptr, pair.first_index)
    cdef object second_name = cy_world_name(system_ptr, pair.second_index) if has_partner else None
    cdef dict worlds = {}
    if has_partner:
        worlds[first_name] = cy_evolution_to_dict(pair.first, system_ptr)
        worlds[second_name] = cy_evolution_to_dict(pair.second, system_ptr)
    return {
        'world_names':         (first_name, second_name),
        'evolved':             True if pair.evolved else False,
        'has_tide_model':      True if pair.has_tide_model else False,
        'orbital_frequency':   pair.orbital_frequency,
        'semi_major_axis':     pair.semi_major_axis,
        'eccentricity':        pair.eccentricity,
        'da_dt':               pair.da_dt,
        'de_dt':               pair.de_dt,
        'dn_dt':               pair.dn_dt,
        'tidal_heating_total': pair.tidal_heating_total,
        'dE_orbit_dt':         pair.dE_orbit_dt,
        'dE_spin_dt_total':    pair.dE_spin_dt_total,
        'energy_residual':     pair.energy_residual,
        'worlds':              worlds,
    }


# ======================================================================================================================
# System.evolve results
# ======================================================================================================================
# A copy of a C++ vector of doubles as a numpy array.
cdef object cy_double_array(const vector[double]& values):
    cdef object out = np.empty(values.size(), dtype=np.float64)
    cdef double[::1] out_view = out
    if values.size() > 0:
        memcpy(&out_view[0], values.data(), values.size() * sizeof(double))
    return out


cdef dict cy_segment_to_dict(const c_PairSegment& segment, tuple world_names):
    cdef dict worlds = {}
    cdef size_t b
    for b in range(segment.bodies.size()):
        worlds[world_names[b]] = {
            "reference_ratio": None if segment.bodies[b].rigid else segment.bodies[b].reference_ratio,
            "rigid": bool(segment.bodies[b].rigid),
            "crossing_armed": bool(segment.bodies[b].armed),
        }
    return {
        "time": segment.time,
        "ended": segment.ended.decode("utf-8") if segment.ended.size() > 0 else None,
        "worlds": worlds,
    }


cdef class EvolutionResult:
    """One world's evolution in a pair (:meth:`System.evolve`), at every step the integrator stored: a view of the C++
    record, its arrays copied out on first access.

    Attributes
    ----------
    world_name : str
    time : np.ndarray
        [s], shared with the pair.
    semi_major_axis : np.ndarray
        The shared orbit's semi-major axis [m].
    eccentricity : np.ndarray
        The shared orbit's eccentricity.
    spin_ratio : np.ndarray
        The world's spin frequency over the orbital motion.
    spin_frequency : np.ndarray
        [rad s-1]
    tidal_heating : np.ndarray
        The world's tidal heating [W].
    da_dt, de_dt : np.ndarray
        This world's contribution to the orbit's rates [m s-1], [s-1]; the pair's sums are on the
        :class:`PairedEvolutionResult`.
    dspin_dt : np.ndarray
        The world's spin rate [rad s-2].
    temperature : np.ndarray or None
        Each layer's temperature [K], shape (layers, steps); None when the world's thermal state was not evolved.
    rigid : bool
        The world has no tide model: it raises no tide and its spin frequency stays as it is.
    num_tide_solves, num_eos_solves : int
        This world's tide solves and EOS solves.
    """
    cdef shared_ptr[c_PairEvolutionRecord] _record
    cdef size_t _body
    cdef readonly str world_name
    cdef tuple _world_names
    cdef dict _arrays

    def __init__(self):
        raise TypeError("TidalPy: an EvolutionResult comes from System.evolve.")

    @staticmethod
    cdef EvolutionResult _wrap(shared_ptr[c_PairEvolutionRecord] record, size_t body, tuple world_names):
        cdef EvolutionResult result = EvolutionResult.__new__(EvolutionResult)
        result._record = record
        result._body = body
        result._world_names = world_names
        result.world_name = world_names[body]
        result._arrays = {}
        return result

    def __reduce__(self):
        # Rebuilt from a copy of its pair's arrays (the segments stay with the pair).
        return (evolution_result_from_state, (cy_record_state(self._record.get(), self._world_names, []),
                                              self.world_name))

    cdef object _array(self, str name):
        cdef object array = self._arrays.get(name)
        if array is not None:
            return array
        cdef c_PairEvolutionRecord* record = self._record.get()
        cdef c_BodyEvolutionRecord* body = &record.bodies[self._body]
        if name == "time":
            array = cy_double_array(record.time)
        elif name == "semi_major_axis":
            array = cy_double_array(record.semi_major_axis)
        elif name == "eccentricity":
            array = cy_double_array(record.eccentricity)
        elif name == "spin_ratio":
            array = cy_double_array(body.spin_ratio)
        elif name == "spin_frequency":
            array = cy_double_array(body.spin_frequency)
        elif name == "tidal_heating":
            array = cy_double_array(body.tidal_heating)
        elif name == "da_dt":
            array = cy_double_array(body.da_dt)
        elif name == "de_dt":
            array = cy_double_array(body.de_dt)
        elif name == "dspin_dt":
            array = cy_double_array(body.dspin_dt)
        array.flags.writeable = False   # Shared by every read
        self._arrays[name] = array
        return array

    @property
    def time(self):
        """Time of each stored step [s]."""
        return self._array("time")

    @property
    def semi_major_axis(self):
        """The shared orbit's semi-major axis [m]."""
        return self._array("semi_major_axis")

    @property
    def eccentricity(self):
        """The shared orbit's eccentricity."""
        return self._array("eccentricity")

    @property
    def spin_ratio(self):
        """The world's spin frequency over the orbital motion."""
        return self._array("spin_ratio")

    @property
    def spin_frequency(self):
        """The world's spin frequency [rad s-1]."""
        return self._array("spin_frequency")

    @property
    def tidal_heating(self):
        """The world's tidal heating [W]."""
        return self._array("tidal_heating")

    @property
    def da_dt(self):
        """The orbit's semi-major-axis rate [m s-1]."""
        return self._array("da_dt")

    @property
    def de_dt(self):
        """The orbit's eccentricity rate [s-1]."""
        return self._array("de_dt")

    @property
    def dspin_dt(self):
        """The world's spin rate [rad s-2]."""
        return self._array("dspin_dt")

    @property
    def temperature(self):
        """Each layer's temperature [K], shape (layers, steps); None when the world's thermal state was not evolved."""
        cdef c_BodyEvolutionRecord* body = &self._record.get().bodies[self._body]
        if not body.thermal:
            return None
        cdef object array = self._arrays.get("temperature")
        if array is None:
            array = cy_double_array(body.temperature).reshape(
                (self._record.get().time.size(), body.num_layers)).T.copy()
            array.flags.writeable = False   # Shared by every read
            self._arrays["temperature"] = array
        return array

    @property
    def rigid(self) -> bool:
        """The world has no tide model: it raises no tide and its spin frequency stays as it is."""
        return True if self._record.get().bodies[self._body].rigid else False

    @property
    def num_tide_solves(self) -> int:
        """This world's tide solves."""
        return <long long>self._record.get().bodies[self._body].num_tide_solves

    @property
    def num_eos_solves(self) -> int:
        """This world's EOS solves."""
        return <long long>self._record.get().bodies[self._body].num_eos_solves

    def __repr__(self):
        return (f"EvolutionResult('{self.world_name}', steps={self._record.get().time.size()}, "
                f"rigid={self.rigid}, thermal={self.temperature is not None})")


cdef class PairedEvolutionResult:
    """The evolution of a pair of worlds (:meth:`System.evolve`): a mapping of each world's name to its
    :class:`EvolutionResult`, the world evolved first, then its partner, over one C++ record.

    Attributes
    ----------
    world_names : tuple of str
        The world evolved, then its partner.
    time : np.ndarray
        [s]
    semi_major_axis : np.ndarray
        [m]
    eccentricity : np.ndarray
    da_dt, de_dt, dn_dt : np.ndarray
        The orbit's rates, both worlds' contributions summed [m s-1], [s-1], [rad s-2].
    segments : list of dict
        One entry per integration segment: its start ``time`` [s], ``ended`` (how it ended when it did not reach the
        end: ``"event"`` or the failure; None otherwise), and ``worlds``, each world's ``reference_ratio`` (the
        commensurability its spin offset is measured from; None for a rigid world), ``rigid``, and
        ``crossing_armed`` (whether the segment stops where the spin crosses its reference).
    success : bool
        The run reached the end of its span and left the system at its final state.
    message : str
    elapsed : float
        Wall-clock time [s].
    num_tide_solves, num_eos_solves : int
        Both worlds' tide and EOS solves.
    num_rhs_calls, num_jacobians : int
        The integrator's right-hand-side evaluations and the Jacobians formed (each takes one evaluation per state
        variable).
    """
    cdef shared_ptr[c_PairEvolutionRecord] _record
    cdef readonly tuple world_names
    cdef dict _worlds
    cdef dict _arrays
    cdef list _segments

    def __init__(self):
        raise TypeError("TidalPy: a PairedEvolutionResult comes from System.evolve.")

    @staticmethod
    cdef PairedEvolutionResult _wrap(shared_ptr[c_PairEvolutionRecord] record, tuple world_names):
        cdef PairedEvolutionResult result = PairedEvolutionResult.__new__(PairedEvolutionResult)
        result._record = record
        result.world_names = world_names
        result._worlds = {
            world_names[b]: EvolutionResult._wrap(record, b, world_names) for b in range(len(world_names))}
        result._arrays = {}
        result._segments = None
        return result

    cdef object _array(self, str name):
        cdef object array = self._arrays.get(name)
        if array is not None:
            return array
        cdef c_PairEvolutionRecord* record = self._record.get()
        if name == "time":
            array = cy_double_array(record.time)
        elif name == "semi_major_axis":
            array = cy_double_array(record.semi_major_axis)
        elif name == "eccentricity":
            array = cy_double_array(record.eccentricity)
        elif name == "da_dt":
            array = cy_double_array(record.da_dt)
        elif name == "de_dt":
            array = cy_double_array(record.de_dt)
        elif name == "dn_dt":
            array = cy_double_array(record.dn_dt)
        array.flags.writeable = False   # Shared by every read
        self._arrays[name] = array
        return array

    def __reduce__(self):
        # Rebuilt from a copy of its arrays and segments.
        return (paired_evolution_result_from_state,
                (cy_record_state(self._record.get(), self.world_names, self.segments),))

    # Mapping of world name to EvolutionResult
    def __getitem__(self, world_name):
        return self._worlds[world_name]

    def __iter__(self):
        return iter(self.world_names)

    def __len__(self):
        return len(self.world_names)

    def __contains__(self, world_name):
        return world_name in self._worlds

    def keys(self):
        return self._worlds.keys()

    def values(self):
        return self._worlds.values()

    def items(self):
        return self._worlds.items()

    def get(self, world_name, default=None):
        return self._worlds.get(world_name, default)

    @property
    def worlds(self) -> dict:
        """World name to :class:`EvolutionResult`."""
        return dict(self._worlds)

    @property
    def time(self):
        """Time of each stored step [s]."""
        return self._array("time")

    @property
    def semi_major_axis(self):
        """The shared orbit's semi-major axis [m]."""
        return self._array("semi_major_axis")

    @property
    def eccentricity(self):
        """The shared orbit's eccentricity."""
        return self._array("eccentricity")

    @property
    def da_dt(self):
        """The orbit's semi-major-axis rate [m s-1]."""
        return self._array("da_dt")

    @property
    def de_dt(self):
        """The orbit's eccentricity rate [s-1]."""
        return self._array("de_dt")

    @property
    def dn_dt(self):
        """The orbit's mean-motion rate, both worlds' contributions summed [rad s-2]."""
        return self._array("dn_dt")

    @property
    def segments(self) -> list:
        """One dict per integration segment (see the class)."""
        cdef size_t i
        if self._segments is None:
            self._segments = [cy_segment_to_dict(self._record.get().segments[i], self.world_names)
                              for i in range(self._record.get().segments.size())]
        return self._segments

    @property
    def success(self) -> bool:
        """The run reached the end of its span and left the system at its final state."""
        return True if self._record.get().success else False

    @property
    def message(self) -> str:
        """How the run ended."""
        return self._record.get().message.decode("utf-8")

    @property
    def elapsed(self) -> float:
        """Wall-clock time [s]."""
        return self._record.get().elapsed

    @property
    def num_tide_solves(self) -> int:
        """Both worlds' tide solves."""
        return sum(result.num_tide_solves for result in self._worlds.values())

    @property
    def num_eos_solves(self) -> int:
        """Both worlds' EOS solves."""
        return sum(result.num_eos_solves for result in self._worlds.values())

    @property
    def num_rhs_calls(self) -> int:
        """The integrator's right-hand-side evaluations."""
        return <long long>self._record.get().num_rhs_calls

    @property
    def num_jacobians(self) -> int:
        """The Jacobians formed."""
        return <long long>self._record.get().num_jacobians

    def __repr__(self):
        return (f"PairedEvolutionResult({self.world_names}, steps={self._record.get().time.size()}, "
                f"success={self.success})")


Mapping.register(PairedEvolutionResult)


# A copy of a pair record's arrays, outcome, and segments as plain Python values.
cdef dict cy_record_state(c_PairEvolutionRecord* record, tuple world_names, list segments):
    cdef list bodies = []
    cdef size_t b
    for b in range(record.bodies.size()):
        bodies.append({
            "world_index": <long long>record.bodies[b].world_index,
            "rigid": bool(record.bodies[b].rigid),
            "thermal": bool(record.bodies[b].thermal),
            "num_layers": <long long>record.bodies[b].num_layers,
            "spin_ratio": cy_double_array(record.bodies[b].spin_ratio),
            "spin_frequency": cy_double_array(record.bodies[b].spin_frequency),
            "tidal_heating": cy_double_array(record.bodies[b].tidal_heating),
            "da_dt": cy_double_array(record.bodies[b].da_dt),
            "de_dt": cy_double_array(record.bodies[b].de_dt),
            "dspin_dt": cy_double_array(record.bodies[b].dspin_dt),
            "temperature": cy_double_array(record.bodies[b].temperature),
            "num_tide_solves": <long long>record.bodies[b].num_tide_solves,
            "num_eos_solves": <long long>record.bodies[b].num_eos_solves,
        })
    return {
        "world_names": world_names,
        "time": cy_double_array(record.time),
        "semi_major_axis": cy_double_array(record.semi_major_axis),
        "eccentricity": cy_double_array(record.eccentricity),
        "da_dt": cy_double_array(record.da_dt),
        "de_dt": cy_double_array(record.de_dt),
        "dn_dt": cy_double_array(record.dn_dt),
        "success": bool(record.success),
        "message": record.message.decode("utf-8"),
        "elapsed": record.elapsed,
        "num_rhs_calls": <long long>record.num_rhs_calls,
        "num_jacobians": <long long>record.num_jacobians,
        "segments": copy.deepcopy(segments),
        "bodies": bodies,
    }


# Fills a C++ vector of doubles with a copy of an array.
cdef void cy_fill_vector(vector[double]& values, object array):
    cdef double[::1] array_view = np.ascontiguousarray(array, dtype=np.float64)
    values.resize(array_view.shape[0])
    if array_view.shape[0] > 0:
        memcpy(values.data(), &array_view[0], array_view.shape[0] * sizeof(double))


def paired_evolution_result_from_state(dict state) -> PairedEvolutionResult:
    """A :class:`PairedEvolutionResult` rebuilt from a copy of its arrays and segments (unpickling and copies)."""
    cdef shared_ptr[c_PairEvolutionRecord] record = c_new_pair_evolution_record()
    cdef c_PairEvolutionRecord* record_ptr = record.get()
    cdef size_t b
    cdef dict body_state
    cy_fill_vector(record_ptr.time, state["time"])
    cy_fill_vector(record_ptr.semi_major_axis, state["semi_major_axis"])
    cy_fill_vector(record_ptr.eccentricity, state["eccentricity"])
    cy_fill_vector(record_ptr.da_dt, state["da_dt"])
    cy_fill_vector(record_ptr.de_dt, state["de_dt"])
    cy_fill_vector(record_ptr.dn_dt, state["dn_dt"])
    record_ptr.success = <cpp_bool>bool(state["success"])
    record_ptr.message = state["message"].encode("utf-8")
    record_ptr.elapsed = <double>state["elapsed"]
    record_ptr.num_rhs_calls = <size_t>state["num_rhs_calls"]
    record_ptr.num_jacobians = <size_t>state["num_jacobians"]
    for b in range(len(state["bodies"])):
        body_state = state["bodies"][b]
        record_ptr.bodies[b].world_index = <size_t>body_state["world_index"]
        record_ptr.bodies[b].rigid = <cpp_bool>body_state["rigid"]
        record_ptr.bodies[b].thermal = <cpp_bool>body_state["thermal"]
        record_ptr.bodies[b].num_layers = <size_t>body_state["num_layers"]
        cy_fill_vector(record_ptr.bodies[b].spin_ratio, body_state["spin_ratio"])
        cy_fill_vector(record_ptr.bodies[b].spin_frequency, body_state["spin_frequency"])
        cy_fill_vector(record_ptr.bodies[b].tidal_heating, body_state["tidal_heating"])
        cy_fill_vector(record_ptr.bodies[b].da_dt, body_state["da_dt"])
        cy_fill_vector(record_ptr.bodies[b].de_dt, body_state["de_dt"])
        cy_fill_vector(record_ptr.bodies[b].dspin_dt, body_state["dspin_dt"])
        cy_fill_vector(record_ptr.bodies[b].temperature, body_state["temperature"])
        record_ptr.bodies[b].num_tide_solves = <size_t>body_state["num_tide_solves"]
        record_ptr.bodies[b].num_eos_solves = <size_t>body_state["num_eos_solves"]
    cdef PairedEvolutionResult result = PairedEvolutionResult._wrap(record, tuple(state["world_names"]))
    result._segments = list(state["segments"])
    return result


def evolution_result_from_state(dict state, str world_name) -> EvolutionResult:
    """An :class:`EvolutionResult` rebuilt from a copy of its pair's arrays (unpickling and copies)."""
    return paired_evolution_result_from_state(state)[world_name]


# The configurations a copy or a pickle of a system carries over: its own and each member world's.
cdef dict cy_system_configs(System system):
    return {
        "source_config": system.source_config,
        "source_dir": system.source_dir,
        "world_configs": [cy_world_configs(world) for world in system._world_wrappers],
    }


def system_from_bytes(bytes record, dict configs=None):
    """A system rebuilt from one binary system record held in memory, as :meth:`System.copy` and unpickling make one.

    Parameters
    ----------
    record : bytes
        The system's binary record, the bytes ``save_binary`` writes to a file.
    configs : dict, optional
        ``source_config``, ``source_dir``, and ``world_configs`` (one dict of each member world's configurations, in
        system order) for the new system, each copied; a key left out leaves that configuration ``None``.

    Returns
    -------
    System
        A system with the record's worlds, roles, and orbits, its worlds unsolved.

    Raises
    ------
    IOError
        The record is not a system's, or it is corrupt.
    """
    cdef System system = System()
    try:
        system._system.get().load_binary_bytes(record, b"a System record", False)
    except RuntimeError as exc:
        raise IOError(str(exc)) from exc
    system._rebuild_world_wrappers()
    configs = configs or {}
    system.source_config = copy.deepcopy(configs.get("source_config"))
    system.source_dir = configs.get("source_dir")
    for world, world_configs in zip(system._world_wrappers, configs.get("world_configs") or []):
        cy_set_world_configs(world, world_configs)
    return system


cdef class System:
    """A gravitationally bound set of worlds, each with its own tidal host.

    Parameters
    ----------
    name : str, optional
        Human-readable system name. Default ``""``.

    Notes
    -----
    Add worlds with :meth:`add_world`, naming each one's ``tidal_host`` there or afterwards with
    :meth:`set_tidal_host`. A world carries a two-body orbit about its host (``semi_major_axis`` [m],
    ``eccentricity``); the mean motion is derived from Kepler's third law. Two worlds that host each other
    (the Earth and the Moon) share one orbit, so its elements need to be given on only one of them. A world
    interacts only with its tidal host and the star.
    """

    def __cinit__(self, *args, **kwargs):
        self._world_wrappers = []
        self.source_config = None
        self.source_dir = None

    def __init__(self, str name=""):
        # The owning member is this same type, so make_unique's result moves straight in.
        self._system = make_unique[c_System](<string>name.encode("utf-8"))
        # Wire the inherited TidalPyBaseClass._ptr to the owned c_System so save_binary / load_binary /
        # get_schema_version_str / save_config resolve to it (c_System : c_TidalPyBaseClass).
        self._ptr = self._system.get()

    def __dealloc__(self):
        self._system.reset()
        self._ptr = NULL

    @staticmethod
    def build(source, overrides=None, force=False):
        """Build a system from a configuration source (the public builder entry point).

        Resolves ``source``, validates it, builds each member world, and returns the assembled
        ``System`` with the normalized configuration retained on :attr:`source_config`. A path to a binary file
        (``save_binary``) is loaded instead (:func:`~TidalPy.Structures.load_system`).

        Parameters
        ----------
        source : str, os.PathLike, or dict
            A bundled system name, a path to a ``.toml`` or binary file, or a system configuration dict.
        overrides : dict, optional
            Values merged over the configuration before it is built, table by table, so a nested table only needs the
            keys it changes (``{"worlds": {"earth": {"world": "earth_prem"}}}``). Not accepted with a binary file.
            Default None.
        force : bool, optional
            If True, bypass the schema-version compatibility warning. Default False.

        Returns
        -------
        System
            The constructed system, with ``source_config`` populated.

        Raises
        ------
        ValueError
            The configuration fails validation, or ``overrides`` is given with a binary file.
        """
        # Deferred imports: the builder helpers import the System class, so importing them at module
        # load would be circular.
        from TidalPy.Structures.configs.system_builder import _resolve_source, construct_system, load_system
        from TidalPy.Structures.configs.toml_loader import load_build_config
        from TidalPy.Structures.configs.worldpack import binary_source_class

        # `_resolve_source` accepts a path string, a Path, or an already-parsed mapping, so `object` it is.
        cdef object resolved = _resolve_source(source)
        if binary_source_class(resolved) is not None:
            if overrides is not None:
                raise ValueError(
                    f"TidalPy: '{resolved}' is a binary file, which overrides cannot change; load it and then change "
                    "the system.")
            return load_system(resolved, force=force)
        cdef dict config = load_build_config(resolved, overrides, force)
        # A member's relative world path is relative to the system file.
        base_dir = os.path.dirname(os.path.abspath(resolved)) if isinstance(resolved, str) else None
        return construct_system(config, force=force, base_dir=base_dir)

    def add_world(
            self,
            BaseWorld world not None,
            tidal_host=None,
            cpp_bool is_star=False,
            semi_major_axis=None,
            eccentricity=None,
            cpp_bool synchronous=False,
            stellar_semi_major_axis=None,
            stellar_eccentricity=None):
        """Add a world to the system, returning its index.

        Parameters
        ----------
        world : BaseWorld
            An initialized world (``BaseWorld`` or a subclass: ``TerrestrialWorld``, ``GasGiantWorld``,
            ``StarWorld``). The system co-owns it; the wrapper stays fully usable.
        tidal_host : int or str or BaseWorld, optional
            The world that raises this world's tides, already a member of the system and identified by
            index, name, or the world object. ``None`` leaves the world without a tidal host; name one
            later with :meth:`set_tidal_host`, which is how a host added after the worlds it hosts (or the
            second member of a mutual pair) is named.
        is_star : bool, optional
            If ``True`` the world becomes the system star, the insolation source (the last world added
            as star wins). A star can also be a tidal host.
        semi_major_axis : float, optional
            Two-body semi-major axis about the tidal host [m]. ``None`` leaves it unset.
        eccentricity : float, optional
            Orbital eccentricity about the tidal host. ``None`` (default) leaves it unset: it reads as ``0.0``
            unless an element set describing the same orbit gives it (the partner of a mutual pair, or the orbit
            about the star when the star is the tidal host).
        synchronous : bool, optional
            Set the world's spin frequency to its mean motion about its tidal host
            (:meth:`set_synchronous_rotation`). Needs ``tidal_host`` and ``semi_major_axis``. Default False.
        stellar_semi_major_axis : float, optional
            Semi-major axis about the star [m] (:meth:`set_stellar_semi_major_axis`), for a world whose tidal host
            is not the star. ``None`` (default) leaves it unset.
        stellar_eccentricity : float, optional
            Orbital eccentricity about the star (:meth:`set_stellar_eccentricity`). ``None`` (default) leaves it
            unset.

        Returns
        -------
        int
            The world's index within the system.

        Raises
        ------
        ValueError
            An element is out of range, ``synchronous`` is given without a tidal host and a semi-major axis, or
            stellar elements are given for a world whose tidal host is the star (its orbit about the star is its
            tidal orbit, given by ``semi_major_axis`` and ``eccentricity``). A refused world is not added.
        """
        cdef double a = NAN if semi_major_axis is None else <double>semi_major_axis
        cdef double e = NAN if eccentricity is None else <double>eccentricity
        cdef double stellar_a = NAN if stellar_semi_major_axis is None else <double>stellar_semi_major_axis
        cdef double stellar_e = NAN if stellar_eccentricity is None else <double>stellar_eccentricity
        cdef cpp_bool has_stellar_elements = (stellar_semi_major_axis is not None) or (stellar_eccentricity is not None)
        cdef str world_name = world._world_ptr.get().get_name().decode("utf-8")
        # Everything is checked before the world is added, so a refused world leaves the system as it was.
        cdef Py_ssize_t host_index = -1 if tidal_host is None else self._resolve_index(tidal_host)
        if synchronous and ((host_index < 0) or (semi_major_axis is None)):
            raise ValueError(
                f"TidalPy: world '{world_name}' cannot be added with synchronous rotation: its mean motion needs a "
                "tidal_host and a semi_major_axis about it.")
        if has_stellar_elements and (host_index >= 0) and (host_index == self._system.get().get_star_index()):
            raise ValueError(
                f"TidalPy: world '{world_name}' has the star as its tidal host, so its orbit about the star is its "
                "tidal orbit; give it as semi_major_axis and eccentricity, not the stellar elements.")
        c_check_orbit(stellar_a, stellar_e, world._world_ptr.get().get_name())
        cdef size_t index = self._system.get().add_world(
            world._world_ptr,
            is_star,
            a,
            e)
        self._world_wrappers.append(world)
        if host_index >= 0:
            self._system.get().set_tidal_host(index, <size_t>host_index)
        if stellar_semi_major_axis is not None:
            self._system.get().set_stellar_semi_major_axis(index, stellar_a)
        if stellar_eccentricity is not None:
            self._system.get().set_stellar_eccentricity(index, stellar_e)
        if synchronous:
            self.set_synchronous_rotation(<int>index)
        return <int>index

    def set_synchronous_rotation(self, world) -> float:
        """Set a world's spin frequency to its mean motion about its tidal host, returning it [rad s-1].

        A synchronously rotating world (most large moons) has its spin equal to its orbital mean motion, so its tides
        have no slow forcing term. The spin is set once, from the current orbit (:meth:`calc_orbital_frequency`): a
        later change of the orbit or the masses does not move it, so call this again after one.

        Parameters
        ----------
        world : int or str or BaseWorld
            The world, identified by index, name, or the world object.

        Returns
        -------
        float
            The spin frequency set [rad s-1].

        Raises
        ------
        ValueError
            The world has no tidal host, or no semi-major axis about it, so it has no mean motion.
        """
        cdef size_t index = <size_t>self._resolve_index(world)
        cdef double orbital_frequency = self._system.get().calc_orbital_frequency(index)
        if not isfinite(orbital_frequency):
            raise ValueError(
                f"TidalPy: world '{cy_world_name(self._system.get(), index)}' has no mean motion to rotate "
                "synchronously with: give it a tidal host and a semi-major axis about it first.")
        self._world_wrappers[index].set_spin_frequency(orbital_frequency)
        return orbital_frequency

    @property
    def name(self) -> str:
        """System name."""
        return self._system.get().get_name().decode("utf-8")

    @name.setter
    def name(self, str value):
        self._system.get().set_name(value.encode("utf-8"))

    @property
    def num_worlds(self) -> int:
        """Number of worlds in the system."""
        return <int>self._system.get().get_num_worlds()

    @property
    def worlds(self) -> list:
        """List of the worlds in the system (in the order they were added)."""
        return list(self._world_wrappers)

    # Tidal hosts (one per world, or none; identify a world by index, name, or object)
    def set_tidal_host(self, world, tidal_host):
        """Name the world that raises ``world``'s tides. Both must already be members of the system.

        ``tidal_host=None`` leaves the world without a tidal host. Two worlds may host each other; they then
        share one orbit (see :meth:`is_mutual_pair`).

        Raises
        ------
        ValueError
            If a world is named as its own tidal host.
        """
        cdef size_t index = <size_t>self._resolve_index(world)
        if tidal_host is None:
            self._system.get().clear_tidal_host(index)
        else:
            self._system.get().set_tidal_host(index, <size_t>self._resolve_index(tidal_host))

    def has_tidal_host(self, world) -> bool:
        """Whether ``world`` has a tidal host."""
        return True if self._system.get().has_tidal_host(<size_t>self._resolve_index(world)) else False

    def get_tidal_host_index(self, world) -> int:
        """Index of ``world``'s tidal host, or ``-1`` if it has none."""
        return self._system.get().get_tidal_host_index(<size_t>self._resolve_index(world))

    def get_tidal_host(self, world):
        """The world that raises ``world``'s tides, or ``None`` if it has none."""
        cdef int host_index = self._system.get().get_tidal_host_index(<size_t>self._resolve_index(world))
        if host_index < 0:
            return None
        return self._world_wrappers[host_index]

    def is_mutual_pair(self, world) -> bool:
        """Whether ``world`` and its tidal host host each other.

        The two then share one orbit: one of them may leave its semi-major axis unset and take its
        partner's elements, and when both carry them they must agree.
        """
        return True if self._system.get().is_mutual_pair(<size_t>self._resolve_index(world)) else False

    def set_star(self, world):
        """Designate an already-added world as the star (by index, name, or the world object)."""
        self._system.get().set_star(<size_t>self._resolve_index(world))

    @property
    def has_star(self) -> bool:
        """Whether a star world has been designated."""
        return True if self._system.get().has_star() else False

    @property
    def star_index(self) -> int:
        """Index of the star world, or ``-1`` if none has been set."""
        return self._system.get().get_star_index()

    @property
    def star(self):
        """The star world wrapper, or ``None`` if no star has been set."""
        if not self._system.get().has_star():
            return None
        return self._world_wrappers[self._system.get().get_star_index()]

    def get_star_luminosity(self) -> float:
        """The star's luminosity [W] (NaN if the star world is not a ``StarWorld``).

        Raises ``RuntimeError`` if no star is set.
        """
        return self._system.get().get_star_luminosity()

    # Orbital elements about the tidal host (per orbiting world; identify a world by index, name, or object)
    def set_semi_major_axis(self, world, double semi_major_axis):
        """Set a world's semi-major axis about its tidal host [m]. For a world whose tidal host is the star this
        is also its semi-major axis about the star, and for a mutual pair it is the partner's too (one orbit)."""
        self._system.get().set_semi_major_axis(<size_t>self._resolve_index(world), semi_major_axis)

    def set_eccentricity(self, world, double eccentricity):
        """Set a world's orbital eccentricity about its tidal host. For a world whose tidal host is the star this
        is also its eccentricity about the star, and for a mutual pair it is the partner's too (one orbit)."""
        self._system.get().set_eccentricity(<size_t>self._resolve_index(world), eccentricity)

    def get_semi_major_axis(self, world) -> float:
        """A world's semi-major axis about its tidal host [m]; its partner's, for the member of a mutual pair
        that carries none. Raises ``ValueError`` when the two members of a mutual pair disagree."""
        return self._system.get().get_semi_major_axis(<size_t>self._resolve_index(world))

    def get_eccentricity(self, world) -> float:
        """A world's orbital eccentricity about its tidal host (see :meth:`get_semi_major_axis`)."""
        return self._system.get().get_eccentricity(<size_t>self._resolve_index(world))

    def calc_gravitational_parameter(self, world) -> float:
        """Standard gravitational parameter ``mu = G (M_host + M_world)`` [m^3 s-2].

        NaN for a world with no tidal host.
        """
        return self._system.get().calc_gravitational_parameter(<size_t>self._resolve_index(world))

    def calc_orbital_frequency(self, world) -> float:
        """Mean motion ``n = sqrt(mu / a^3)`` [rad s-1] for a world's two-body orbit about its tidal host.

        Returns NaN for a non-positive semi-major axis or a world with no tidal host.
        """
        return self._system.get().calc_orbital_frequency(<size_t>self._resolve_index(world))

    def calc_semi_major_axis_from_frequency(self, world, double orbital_frequency) -> float:
        """Semi-major axis ``a = (mu / n^2)^(1/3)`` [m] from a mean motion (inverse of the above)."""
        return self._system.get().calc_semi_major_axis_from_frequency(
            <size_t>self._resolve_index(world), orbital_frequency)

    # Orbital elements about the star (the insolation source; may differ from the tidal-host orbit)
    def set_stellar_semi_major_axis(self, world, double semi_major_axis):
        """Set a world's semi-major axis about the star [m]."""
        self._system.get().set_stellar_semi_major_axis(<size_t>self._resolve_index(world), semi_major_axis)

    def set_stellar_eccentricity(self, world, double eccentricity):
        """Set a world's orbital eccentricity about the star."""
        self._system.get().set_stellar_eccentricity(<size_t>self._resolve_index(world), eccentricity)

    def get_stellar_semi_major_axis(self, world) -> float:
        """A world's semi-major axis about the star [m]."""
        return self._system.get().get_stellar_semi_major_axis(<size_t>self._resolve_index(world))

    def get_stellar_eccentricity(self, world) -> float:
        """A world's orbital eccentricity about the star."""
        return self._system.get().get_stellar_eccentricity(<size_t>self._resolve_index(world))

    def calc_stellar_gravitational_parameter(self, world) -> float:
        """Standard gravitational parameter ``mu = G (M_star + M_world)`` [m^3 s-2].

        Raises ``RuntimeError`` if no star is set; returns NaN for the star's own entry.
        """
        return self._system.get().calc_stellar_gravitational_parameter(<size_t>self._resolve_index(world))

    def calc_stellar_orbital_frequency(self, world) -> float:
        """Mean motion ``n = sqrt(mu / a^3)`` [rad s-1] for a world's orbit about the star.

        Returns NaN for a non-positive stellar semi-major axis or the star's own entry.
        """
        return self._system.get().calc_stellar_orbital_frequency(<size_t>self._resolve_index(world))

    def calc_insolation_flux(self, world) -> float:
        """Orbit-averaged incident stellar flux [W m-2] at a world.

        ``F = L_star / (4 pi a^2 sqrt(1-e^2))`` using the world's orbital elements about the star (the
        incident flux before the world's own albedo/emissivity are applied). Raises ``RuntimeError`` if
        no star is set; returns NaN for the star's own entry, an unset stellar semi-major axis, or a
        star with no luminosity.
        """
        return self._system.get().calc_insolation_flux(<size_t>self._resolve_index(world))

    def calc_equilibrium_temperature(self, world) -> float:
        """Surface equilibrium temperature [K] of a world from stellar insolation.

        Gray-body radiative balance using the world's albedo + emissivity,
        ``T = ((1-A) F / (4 eps sigma))^(1/4)``. Returns NaN if the insolation flux is unavailable.
        """
        return self._system.get().calc_equilibrium_temperature(<size_t>self._resolve_index(world))

    def calc_dissipation(self, world) -> dict:
        """One world's tidal dissipation on its orbit about its tidal host, without orbital or spin rates.

        Solves the world's global tides in the current system state, its tidal host raising them: the mean motion
        from Kepler's third law, the spin and obliquity from the world, the semi-major axis and eccentricity from
        its orbit about the host (the orbit both share, for a mutual pair), and the host's mass.
        :meth:`calc_world_evolution` turns this into the world's rates; :meth:`calc_pair_evolution` adds the
        host's dissipation as well.

        Parameters
        ----------
        world : int or str or BaseWorld
            The dissipating world, identified by index, name, or the world object.

        Returns
        -------
        dict
            The world (``world_index``, ``world_name``) and the body raising its tide (``companion_name``), the state
            used (``orbital_frequency``, ``semi_major_axis``, ``eccentricity``, ``spin_frequency``, ``obliquity``,
            ``companion_mass``, ``target_mass``), the tidal outputs (``tidal_heating``, ``dU_dM``, ``dU_dw``,
            ``dU_dO``, ``dU_dM_minus_dw``), and ``moment_of_inertia``, all MKS. ``solved`` is ``False`` for a world
            with no tidal host or no usable orbit about it. ``has_tide_model`` is ``False`` for a rigid world (no
            tide model attached): it raises no tide, so its outputs are zero while ``solved`` stays ``True``. The
            world's per-layer heating is on the world (``get_layer_tidal_heating``).
        """
        cdef size_t index = <size_t>self._resolve_index(world)
        cdef c_TidalDissipation dissipation
        with nogil:
            dissipation = self._system.get().calc_dissipation(index)
        return cy_dissipation_to_dict(dissipation, self._system.get())

    def calc_world_evolution(self, world) -> dict:
        """Evolve one orbiting world for a single tidal solve, returning its rates as a dict.

        Solves the world's global tides in the current system state (:meth:`calc_dissipation`), then turns the
        tidal-potential derivatives into the orbital rates and the world's spin rate. Only this world raises tides;
        its host is treated as a point mass whose state stays as it is, apart from the orbit they share.

        A spin at a spin-orbit equilibrium is not held here: :meth:`evolve` integrates a world through its
        equilibria.

        Parameters
        ----------
        world : int or str or BaseWorld
            The orbiting world, identified by index, name, or the world object.

        Returns
        -------
        dict
            The world (``world_index``, ``world_name``), the orbital and spin state used, the raw tidal outputs
            (``tidal_heating``, ``dU_dM``, ``dU_dw``, ``dU_dO``, ``dU_dM_minus_dw``), the rates (``da_dt``,
            ``de_dt``, ``dn_dt``, ``dspin_dt``), and the energy-balance terms (``dE_orbit_dt``, ``dE_spin_dt``,
            ``energy_residual``), all MKS.
            ``evolved`` is ``False`` for a world with no tidal host or no usable orbit about it.
            ``has_tide_model`` is ``False`` for a rigid world (no tide model attached): it raises no tide, so
            its rates and energy terms are zero while ``evolved`` stays ``True``, and a warning is logged once
            per world.
        """
        cdef size_t index = <size_t>self._resolve_index(world)
        cdef c_WorldEvolution evolution
        with nogil:
            evolution = self._system.get().calc_world_evolution(index)
        return cy_evolution_to_dict(evolution, self._system.get())

    def calc_system_evolution(self) -> list:
        """Evolve every world in the system (single-body dissipation).

        Solves each orbiting world's tides for its current orbit, then returns the rates.

        Returns
        -------
        list of dict
            One :meth:`calc_world_evolution` dict per world, in index order. A world with no tidal host or
            no usable orbit comes back with ``evolved`` set to ``False``, and a rigid world with
            ``has_tide_model`` set to ``False``. The two members of a mutual pair each get an entry, and their
            contributions to the orbit they share add.
        """
        cdef vector[c_WorldEvolution] results
        with nogil:
            results = self._system.get().calc_system_evolution()
        cdef list out = []
        cdef size_t i
        for i in range(results.size()):
            out.append(cy_evolution_to_dict(results[i], self._system.get()))
        return out

    def calc_pair_evolution(self, world, partner=None) -> dict:
        """Evolve two worlds together under dual-body tidal dissipation.

        Both worlds raise a tide on their shared orbit, one being the other's tidal host (most often each other's).
        Each body's tides are solved with the other body as the tide raiser (:meth:`calc_dissipation`); their
        orbital-rate contributions add and each body evolves its own spin. The two are treated alike. A body with
        no tide model is rigid and contributes nothing (with a rigid partner this reduces to
        :meth:`calc_world_evolution`).

        Parameters
        ----------
        world : int or str or BaseWorld
            One world of the pair, identified by index, name, or the world object.
        partner : int or str or BaseWorld, optional
            The other world. None (default) takes ``world``'s tidal host; a world with no tidal host then comes back
            with ``evolved`` set to ``False``.

        Returns
        -------
        dict
            The combined shared-orbit fields (``da_dt``, ``de_dt``, ``dn_dt``, ``tidal_heating_total``,
            ``dE_orbit_dt``, ``dE_spin_dt_total``, ``energy_residual``, ``orbital_frequency``, ``semi_major_axis``,
            ``eccentricity``, ``evolved``, ``has_tide_model``), ``world_names`` (``world``'s name, then
            ``partner``'s), and ``worlds``: each body's full single-body contribution keyed by its world name, each a
            :meth:`calc_world_evolution`-style dict. ``evolved`` is ``False`` for a world with no tidal host
            (``partner`` left out; ``world_names`` then ends in ``None`` and ``worlds`` is empty) or no usable orbit.
            ``has_tide_model`` is ``True`` when at least one body carries a tide model; ``False`` means both are
            rigid, every rate is zero, and a warning is logged once per world. Each body's own flag is in its own
            entry.

        Raises
        ------
        ValueError
            ``partner`` is ``world``, or neither is the other's tidal host (they share no orbit).
        """
        cdef size_t index = <size_t>self._resolve_index(world)
        cdef size_t partner_index = 0
        cdef cpp_bool has_partner = partner is not None
        if has_partner:
            partner_index = <size_t>self._resolve_index(partner)
        cdef c_PairEvolution pair
        with nogil:
            if has_partner:
                pair = self._system.get().calc_pair_evolution(index, partner_index)
            else:
                pair = self._system.get().calc_pair_evolution(index)
        return cy_pair_to_dict(pair, self._system.get())

    def evolve(
            self,
            world,
            time_span,
            *,
            evolve_thermal=None,
            method=None,
            semi_major_axis_rtol=None,
            eccentricity_rtol=None,
            eccentricity_atol=None,
            spin_rtol=None,
            thermal_rtol=None,
            radial_rtol=None,
            radial_atol=None,
            max_wall_time=None):
        """Evolve a world and its tidal host together over a time span.

        The pair is ``world`` and its tidal host (usually each the other's), and the two are treated alike. Integrates
        their shared orbit (a, e), both spins, and, with ``evolve_thermal``, the temperature of every layer of each
        world that has layers, with both bodies raising tides (:meth:`calc_pair_evolution`), all in C++ with CyRK. A
        world with no tide model is rigid: its spin frequency stays as it is. Each spin is integrated as its offset from
        the nearest spin-orbit commensurability j / m, which reaches the tidal mode frequencies exactly. Below a world's
        continuation frequency (:meth:`calc_continuation_frequency
        <TidalPy.Structures.worlds.BaseWorld.calc_continuation_frequency>`) a mode's dissipation falls smoothly and
        linearly to zero, so the tidal torque passes smoothly through each commensurability. A lock is then a stable
        equilibrium of a stiff variable, which the implicit integrator holds with long steps; capture, passage, and
        release need no special handling. The integration restarts where a spin crosses its commensurability, so a lock
        narrower than a step is not stepped over. See the System documentation.

        The orbit, the spins, and the temperatures start from the system's current state, and the system is left at
        the final one. With ``evolve_thermal``, each thermal world's EOS is solved with its temperature profile
        (``solve_temperature``) whenever its temperatures, the time, or its surface temperature change, its surface
        at its insolation temperature from the system's star (no flow through it without a star; add one for a
        thermal run), its heat sources at the integration time, and its layers warming at
        :meth:`calc_layer_temperature_rate <TidalPy.Structures.worlds.BaseWorld.calc_layer_temperature_rate>`;
        without it each world keeps the structure it has (solve its EOS first). A layer that can change state melts
        and freezes inside its own radii.

        Every setting left as None comes from the system file's ``[evolution]`` table, then from the ``[evolution]``
        section of the TidalPy configuration.

        Parameters
        ----------
        world : int or str or BaseWorld
            One world of the pair; it needs a tidal host.
        time_span : tuple of float
            (start, end) [s]; the start is also the heat sources' clock (radiogenic heating counts from the formation
            of the body, so a run of the present day starts near 4.5 Gyr).
        evolve_thermal : bool, optional
            Evolve the layer temperatures of each world with layers.
        method : str, optional
            CyRK's implicit integrator: ``"Radau"`` (the packaged default), ``"BDF"``, or ``"LSODA"``.
        semi_major_axis_rtol, eccentricity_rtol, eccentricity_atol : float, optional
            Tolerances on a / a0 and e.
        spin_rtol : float, optional
            Relative tolerance on each spin's offset from its commensurability (its absolute tolerance follows from
            the world's continuation frequency, ``calc_continuation_frequency``).
        thermal_rtol : float, optional
            Relative tolerance on the layer temperatures.
        radial_rtol, radial_atol : float, optional
            The Love solves' tolerances during the run; the rates must be smooth at the scale of the integrator's
            difference steps. Each world's own settings come back afterwards.
        max_wall_time : float, optional
            Wall-clock cap [s]; a run that reaches it returns with ``success`` False. Infinite: no cap.

        Returns
        -------
        PairedEvolutionResult
            A mapping of each world's name (``world`` first) to its :class:`EvolutionResult` (spin ratio and frequency,
            tidal heating, its contributions to da/dt and de/dt, its dspin/dt, and layer temperatures, at every stored
            step), with the shared time, semi-major axis, eccentricity, and summed rates, the segments, the counts, and
            the outcome.

        Raises
        ------
        ValueError
            The world has no tidal host or no usable orbit about it, a dissipating world's spin is not finite and at
            least zero, the span is not increasing, with ``evolve_thermal`` a world with layers has a layer temperature
            that is not finite and positive, a tolerance is not positive, or the method is not implicit.

        Assumptions
        -----------
        * The worlds have no permanent (triaxial) figure and no rotational flattening; their spins follow the tidal
          torques alone.
        * Layer boundaries are fixed. A thermal world's structure is solved without its tidal heat, which enters only
          its temperature rates.
        * Each world's moment of inertia is constant in the spin equation (no dC/dt term).
        """
        cdef size_t index = <size_t>self._resolve_index(world)
        start, end = time_span
        cdef double t_start = <double>start
        cdef double t_end = <double>end
        given = dict(
            evolve_thermal=evolve_thermal, method=method, semi_major_axis_rtol=semi_major_axis_rtol,
            eccentricity_rtol=eccentricity_rtol, eccentricity_atol=eccentricity_atol, spin_rtol=spin_rtol,
            thermal_rtol=thermal_rtol, radial_rtol=radial_rtol, radial_atol=radial_atol, max_wall_time=max_wall_time)
        options = dict(TidalPy.config["evolution"])
        if self.source_config is not None:
            options.update(self.source_config.get("evolution", {}))
        options.update({key: value for key, value in given.items() if value is not None})
        cdef c_PairEvolveSettings settings
        settings.evolve_thermal       = <cpp_bool>bool(options["evolve_thermal"])
        settings.method               = str(options["method"]).encode("utf-8")
        settings.semi_major_axis_rtol = <double>options["semi_major_axis_rtol"]
        settings.eccentricity_rtol    = <double>options["eccentricity_rtol"]
        settings.eccentricity_atol    = <double>options["eccentricity_atol"]
        settings.spin_rtol            = <double>options["spin_rtol"]
        settings.thermal_rtol         = <double>options["thermal_rtol"]
        settings.radial_rtol          = <double>options["radial_rtol"]
        settings.radial_atol          = <double>options["radial_atol"]
        cdef double wall_cap = <double>float(options["max_wall_time"])
        settings.max_wall_time        = wall_cap if isfinite(wall_cap) else NAN
        cdef shared_ptr[c_PairEvolutionRecord] record
        with nogil:
            record = c_evolve_pair(self._system.get(), index, t_start, t_end, settings)
        cdef size_t b
        cdef tuple world_names = tuple(
            cy_world_name(self._system.get(), record.get().bodies[b].world_index)
            for b in range(record.get().bodies.size()))
        return PairedEvolutionResult._wrap(record, world_names)

    @property
    def config(self):
        """The system configuration dict the system was built from (None if built directly)."""
        return self.source_config

    cpdef dict get_config_dict(self):
        """Return the system's live state as a configuration dict (the ``build_system`` schema).

        Each member world is inlined with its own live configuration (its ``get_config_dict``) under its
        system name, together with its tidal host, its star role, and its orbital elements, so
        ``build_system`` rebuilds the system from it as it stands now and not as it was first described. A world whose
        tidal host is the star carries no stellar elements, since its orbit about the star is its tidal orbit.
        :meth:`get_save_config` builds on it, keeping each unchanged member's original reference.

        Returns
        -------
        dict
            A system configuration dict with ``schema_version``, ``name``, a ``worlds`` table, and the
            ``evolution`` table of the system file it was built from, when it had one.
        """
        from TidalPy.Structures.configs.toml_loader import SCHEMA_VERSION
        cdef c_System* system_ptr = self._system.get()
        cdef int host_index
        cdef int star_index = system_ptr.get_star_index()
        cdef int i
        cdef double a, stellar_a
        cdef dict worlds_table = {}
        # `object`, not `BaseWorld`: `name` and `get_config_dict` are Python-level members defined in
        # base.pyx, so they are not visible through base.pxd and a typed handle could not reach them.
        cdef object world
        cdef dict world_cfg
        cdef dict entry
        for i in range(<int>system_ptr.get_num_worlds()):
            world = self._world_wrappers[i]
            world_cfg = world.get_config_dict()
            entry = {"world": world_cfg}
            host_index = system_ptr.get_tidal_host_index(<size_t>i)
            if host_index >= 0:
                entry["tidal_host"] = self._world_wrappers[host_index].name
            if i == star_index:
                entry["is_star"] = True
            a = system_ptr.get_semi_major_axis(<size_t>i)
            if isfinite(a):
                entry["semi_major_axis_m"] = a
                entry["eccentricity"] = system_ptr.get_eccentricity(<size_t>i)
            # A world whose tidal host is the star has one orbit, written above as its tidal elements; add_world
            # refuses stellar elements for it.
            stellar_a = system_ptr.get_stellar_semi_major_axis(<size_t>i)
            if isfinite(stellar_a) and not system_ptr.is_hosted_by_star(<size_t>i):
                entry["stellar_semi_major_axis_m"] = stellar_a
                entry["stellar_eccentricity"] = system_ptr.get_stellar_eccentricity(<size_t>i)
            worlds_table[world.name] = entry
        config = {"schema_version": SCHEMA_VERSION, "name": self.name, "worlds": worlds_table}
        if (self.source_config is not None) and ("evolution" in self.source_config):
            config["evolution"] = copy.deepcopy(self.source_config["evolution"])
        return config

    def get_save_config(self, destination_dir=None) -> dict:
        """The configuration :meth:`save_to_toml` writes: the system as it is now.

        Every member carries its current tidal host, star role, and orbital elements. A member built from a world
        reference (a bundled name or a file) and unchanged since its build keeps that reference, a relative file path
        rewritten to find the same file from ``destination_dir``; any other member is written inline, as its
        :meth:`~TidalPy.Structures.worlds.base.BaseWorld.get_save_config` gives it, so a change made to a world after
        the build is saved too.

        Parameters
        ----------
        destination_dir : str or os.PathLike, optional
            The folder the configuration will be saved into. Default None: references as given.

        Returns
        -------
        dict
            A system configuration (``name`` and a ``worlds`` table) that builds this system.
        """
        from TidalPy.Structures.configs.config_writer import relocated_path, state_changed_since_build
        cdef dict live = self.get_config_dict()
        cdef dict source_worlds = (self.source_config or {}).get("worlds", {})
        cdef object target_dir = None if destination_dir is None else os.fspath(destination_dir)
        cdef dict worlds_table = {}
        cdef dict entry
        cdef object given
        cdef object world
        # The live table lists the members in system order, so its i-th entry is the i-th wrapper.
        for world, (world_key, live_entry) in zip(self._world_wrappers, live["worlds"].items()):
            entry = {key: value for key, value in live_entry.items() if key != "world"}
            given = source_worlds.get(world_key, {}).get("world")
            if isinstance(given, str) and not state_changed_since_build(world.built_config, live_entry["world"]):
                if ((target_dir is not None) and (self.source_dir is not None) and not os.path.isabs(given)
                        and os.path.isfile(os.path.join(self.source_dir, given))):
                    given = relocated_path(given, os.path.join(self.source_dir, given), target_dir)
                entry["world"] = given
            else:
                entry["world"] = world.get_save_config(target_dir)
            worlds_table[world_key] = entry
        config = {"name": live["name"], "worlds": worlds_table}
        if "evolution" in live:
            config["evolution"] = live["evolution"]
        return config

    def save_to_toml(self, file_path, overwrite=True):
        """Write this system's configuration, as it is now, to a TOML file.

        Writes :meth:`get_save_config`: the current orbits and roles, each member by its original reference when it
        is unchanged since the build, and inline otherwise. Relative references are rewritten to find the same files
        from the folder saved into.

        Parameters
        ----------
        file_path : str or os.PathLike
            Destination ``.toml`` path.
        overwrite : bool, optional
            Overwrite an existing file. Default True.

        Raises
        ------
        ValueError
            The path does not end in ``.toml``, or a member's live configuration fails the world schema.
        FileExistsError
            The file exists and ``overwrite`` is False.
        """
        from TidalPy.Structures.configs.config_writer import save_system_to_toml
        cdef str path = os.fspath(file_path)
        cdef dict config = self.get_save_config(os.path.dirname(os.path.abspath(path)))
        return save_system_to_toml(config, path, overwrite=overwrite)

    cdef void _rebuild_world_wrappers(self):
        """Rebuild the Python world-wrapper list around the C++ system's current worlds.

        Called after load_binary, where c_System::read_binary has replaced the owned worlds with freshly
        deserialized ones; each is re-wrapped as its concrete Python type (co-owning the C++ world).
        """
        cdef c_System* system_ptr = self._system.get()
        cdef size_t num = system_ptr.get_num_worlds()
        cdef size_t i
        self._world_wrappers = []
        for i in range(num):
            self._world_wrappers.append(cy_wrap_world(system_ptr.get_world(i)))

    def load_binary(self, path, cpp_bool force=False):
        """Load this system's state from a TidalPy binary file (overriding the base to rewrap worlds).

        Rebuilds the heterogeneous world list from the stream (each world's concrete type is
        recovered from its record) and the Python wrappers around it. Each world comes back with its tide
        model and ``[tides]`` settings, spin model, pinned solver settings, and (for a star) luminosity
        model; solved state is not saved, so call ``solve_eos`` on each world with layers before evolving.
        :func:`~TidalPy.Structures.load_system` loads a file into a new system.

        Parameters
        ----------
        path : str or os.PathLike
            Source file path.
        force : bool, optional
            If True, attempt to load even on a schema-version mismatch.

        Raises
        ------
        TypeError
            ``path`` is not a string or path-like object.
        FileNotFoundError
            If the file does not exist.
        IOError
            If the file is invalid or the schema version is incompatible, or if its roles or orbits are
            corrupt (a tidal host index that names no other world, an out-of-range star index, a duplicate
            world name, or an orbit that is not bound), or it holds data after the system's record. A failed load
            leaves the system, its worlds, and their wrappers as they were: the file is read into a new system
            first and reaches this one only once that read succeeded.
        """
        cdef str file_path = cy_existing_binary_path(path)
        try:
            self._system.get().load_binary(file_path.encode("utf-8"), force)
        except RuntimeError as exc:
            raise IOError(str(exc)) from exc
        # Only a successful load replaces the worlds, so only then do the wrappers follow the new ones.
        self._rebuild_world_wrappers()
        self.source_config = None
        self.source_dir = None

    def copy(self):
        """An independent system with the same worlds, roles, and orbits.

        The copy goes through the binary record (:meth:`save_binary` without the file), so each world is a copy as
        :meth:`BaseWorld.copy <TidalPy.Structures.worlds.base.BaseWorld.copy>` makes one: same layers, models, and
        settings, unsolved (run ``solve_eos`` on each world with layers again), with its own copies of the
        configurations it was built from. ``copy.copy``, ``copy.deepcopy``, and ``pickle`` use it.

        Returns
        -------
        System
        """
        return system_from_bytes(self._system.get().write_binary_bytes(), cy_system_configs(self))

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def __reduce__(self):
        # Pickled as its binary record and its configurations, rebuilt by system_from_bytes (see copy).
        return (system_from_bytes, (self._system.get().write_binary_bytes(), cy_system_configs(self)))

    def __repr__(self):
        """One line: the class, the name, the member worlds by name, and the star."""
        if self._system.get() == NULL:
            return f"{type(self).__name__}(no system)"
        cdef object star = self.star
        return (f"{type(self).__name__}({self.name!r}, worlds={[world.name for world in self._world_wrappers]!r}, "
                f"star={None if star is None else star.name!r})")

    # World identification: accept an index (int), a world name (str), or the world wrapper object.
    cdef Py_ssize_t _resolve_index(self, object world) except *:
        cdef Py_ssize_t n = len(self._world_wrappers)
        cdef Py_ssize_t i
        cdef int found
        if isinstance(world, bool):
            # A bool is an int to Python, so True would quietly name world 1.
            raise TypeError("a world is named by an int index, a name (str), or the world object, not a bool")
        if isinstance(world, Integral):
            # Any integer type, numpy's included.
            i = <Py_ssize_t>operator.index(world)
            if i < 0:
                i += n
            if i < 0 or i >= n:
                raise IndexError(f"world index {world} out of range (system has {n} worlds)")
            return i
        if isinstance(world, str):
            found = self._system.get().find_world_index(world.encode("utf-8"))
            if found < 0:
                raise KeyError(f"no world named '{world}' in the system")
            return <Py_ssize_t>found
        if isinstance(world, BaseWorld):
            for i in range(n):
                if self._world_wrappers[i] is world:
                    return i
            raise ValueError("world object is not a member of this system")
        raise TypeError(f"world must be an int index, a name (str), or a world object, not {type(world).__name__}")

    def __len__(self):
        """Number of worlds in the system."""
        return len(self._world_wrappers)

    def __iter__(self):
        """Iterate the member worlds."""
        return iter(self._world_wrappers)

    def __getitem__(self, index):
        """Index or slice the worlds: ``system[0]`` (a world) or ``system[:2]`` (a list).

        Integer indices accept negatives; a string returns the world of that name; a slice returns the
        corresponding list of worlds.
        """
        if isinstance(index, slice):
            return self._world_wrappers[index]
        return self._world_wrappers[self._resolve_index(index)]

    def __getattr__(self, name):
        """Resolve ``system.<world_name>`` to that world (after normal attribute lookup).

        Only consulted when ``name`` is not a real attribute/method, so defined members always win.
        Names starting with ``_`` are never treated as worlds.
        """
        if name.startswith("_") or self._system.get() == NULL:
            raise AttributeError(name)
        cdef int found = self._system.get().find_world_index(name.encode("utf-8"))
        if found >= 0:
            return self._world_wrappers[found]
        raise AttributeError(
            f"'{type(self).__name__}' object has no attribute or world named '{name}'")
