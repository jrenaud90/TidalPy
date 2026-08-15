# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
system.pyx
Cython/Python wrapper for TidalPy's system class.

A ``System`` links two or more worlds (a host plus orbiting worlds). Each orbiting world has a two-body
orbit about the host described by a semi-major axis [m] and eccentricity. The system owns its worlds
through shared pointers, co-owned with the Python world wrappers, so the added world objects stay
usable and are handed straight back by iteration (``for world in system``), indexing (``system[i]``),
and attribute access (``system.<world_name>``).

All quantities are MKS: masses [kg], semi-major axes [m], orbital frequencies [rad s-1], gravitational
parameters [m^3 s-2].
"""

from libc.math cimport NAN, isfinite
from libcpp cimport bool as cpp_bool
from libcpp.vector cimport vector

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.structures_x.worlds.base cimport BaseWorld
from TidalPy.structures_x.worlds.layered cimport LayeredWorld
from TidalPy.structures_x.worlds.gasgiant cimport GasGiantWorld
from TidalPy.structures_x.worlds.stellar cimport StarWorld

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())

# Pull in the out-of-line (inline) definitions of the world tidal solves so the system can drive a
# world's tides directly in C++. world_tides_base_.hpp defines the analytic c_BaseWorld::calc_tides;
# world_tides_.hpp defines the layered c_LayeredWorld::calc_tides (rheology + layer heat distribution),
# whose rheology path runs the CyRK-backed radial solver (linked via the CyRK cimport in system.pxd).
cdef extern from "world_tides_base_.hpp" nogil:
    pass
cdef extern from "world_tides_.hpp" nogil:
    pass


cdef BaseWorld _wrap_world(shared_ptr[c_BaseWorld] ptr):
    """Wrap a C++ world (e.g. one loaded by c_System::read_binary) as the matching Python wrapper.

    Dispatches on the world's concrete type so a layered / gas-giant / star world comes back as its own
    wrapper class (with its layered / star methods), not a bare BaseWorld.
    """
    cdef int kind = c_world_kind(ptr.get())
    if kind == 1:
        return LayeredWorld._wrap(ptr)
    if kind == 2:
        return GasGiantWorld._wrap(ptr)
    if kind == 3:
        return StarWorld._wrap(ptr)
    return BaseWorld._wrap(ptr)


cdef dict _evolution_to_dict(c_WorldEvolution evolution):
    """Convert a c_WorldEvolution result into a plain Python dict (all values MKS)."""
    return {
        'world_index':       <int>evolution.world_index,
        'evolved':           True if evolution.evolved else False,
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


cdef dict _pair_to_dict(c_PairEvolution pair):
    """Convert a c_PairEvolution (dual-body) result into a plain Python dict (all values MKS)."""
    return {
        'world_index':         <int>pair.world_index,
        'host_index':          <int>pair.host_index,
        'evolved':             True if pair.evolved else False,
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
        'world':               _evolution_to_dict(pair.world),
        'host':                _evolution_to_dict(pair.host),
    }


# =====================================================================================================================
# System
# =====================================================================================================================
cdef class System:
    """A gravitationally bound set of worlds (a host plus orbiting worlds).

    Parameters
    ----------
    name : str, optional
        Human-readable system name. Default ``""``.

    Notes
    -----
    Add worlds with :meth:`add_world` (one of them flagged ``is_host``). Orbiting worlds carry a
    two-body orbit about the host (``semi_major_axis`` [m], ``eccentricity``); the mean motion is
    derived from Kepler's third law. Worlds interact only with the host, not with one another.
    """

    def __cinit__(self, *args, **kwargs):
        self._world_wrappers = []
        self.source_config = None

    def __init__(self, str name=""):
        self._system.reset(new c_System(name.encode("utf-8")))
        # Wire the inherited TidalPyBaseClass._ptr to the owned c_System so save_binary / load_binary /
        # get_schema_version_str / save_config resolve to it (c_System : c_TidalPyBaseClass).
        self._ptr = self._system.get()

    def __dealloc__(self):
        self._system.reset()
        self._ptr = NULL

    # ------------------------------------------------------------------------------------------------------------------
    # Builder entry point
    # ------------------------------------------------------------------------------------------------------------------
    @staticmethod
    def build(source, force=False):
        """Build a system from a configuration source (the public builder entry point).

        Resolves ``source``, validates it, builds each member world, and returns the assembled
        ``System`` with the normalized configuration retained on :attr:`source_config` (so it can be
        written back to TOML). Mirrors :meth:`BaseWorld.build
        <TidalPy.structures_x.worlds.base.BaseWorld.build>`; the module-level
        :func:`~TidalPy.structures_x.configs.system_builder.build_system` simply calls this method.

        Parameters
        ----------
        source : str or dict
            A bundled system name, a path to a ``.toml`` file, or a system configuration dict.
        force : bool, optional
            If True, bypass the schema-version compatibility warning. Default False.

        Returns
        -------
        System
            The constructed system, with ``source_config`` populated.
        """
        # Deferred imports: the builder helpers import the System class, so importing them at module
        # load would be circular.
        from TidalPy.structures_x.configs.system_builder import _resolve_source, construct_system
        from TidalPy.structures_x.configs.toml_loader import (
            load_toml,
            merge_with_defaults,
            validate_schema_version)

        resolved = _resolve_source(source)
        config = load_toml(resolved)
        config = merge_with_defaults(config)
        validate_schema_version(config, force=force)
        return construct_system(config, force=force)

    # ------------------------------------------------------------------------------------------------------------------
    # Membership
    # ------------------------------------------------------------------------------------------------------------------
    def add_world(
            self,
            BaseWorld world not None,
            cpp_bool is_host=False,
            cpp_bool is_star=False,
            semi_major_axis=None,
            double eccentricity=0.0):
        """Add a world to the system, returning its index.

        Parameters
        ----------
        world : BaseWorld
            An initialized world (``LayeredWorld`` / ``GasGiantWorld`` / ``StarWorld``). The system
            co-owns it; the wrapper stays fully usable.
        is_host : bool, optional
            If ``True`` the world becomes the system tidal host (the last world added as host wins).
        is_star : bool, optional
            If ``True`` the world becomes the system star, the insolation source (the last world added
            as star wins). ``is_host`` and ``is_star`` are independent: a world can be both.
        semi_major_axis : float, optional
            Two-body semi-major axis about the tidal host [m] (orbiting worlds only). ``None`` leaves it
            unset. The orbit about the star is set separately via :meth:`set_stellar_semi_major_axis`.
        eccentricity : float, optional
            Orbital eccentricity about the tidal host. Default ``0.0``.

        Returns
        -------
        int
            The world's index within the system.
        """
        cdef double a = NAN if semi_major_axis is None else <double>semi_major_axis
        cdef size_t index = self._system.get().add_world(world._world_ptr, is_host, is_star, a, eccentricity)
        self._world_wrappers.append(world)
        return <int>index

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

    # ------------------------------------------------------------------------------------------------------------------
    # Host
    # ------------------------------------------------------------------------------------------------------------------
    def set_host(self, world):
        """Designate an already-added world as the host (by index, name, or the world object)."""
        self._system.get().set_host(<size_t>self._resolve_index(world))

    @property
    def has_host(self) -> bool:
        """Whether a host world has been designated."""
        return True if self._system.get().has_host() else False

    @property
    def host_index(self) -> int:
        """Index of the host world, or ``-1`` if none has been set."""
        return self._system.get().get_host_index()

    @property
    def host(self):
        """The host world wrapper, or ``None`` if no host has been set."""
        if not self._system.get().has_host():
            return None
        return self._world_wrappers[self._system.get().get_host_index()]

    # ------------------------------------------------------------------------------------------------------------------
    # Star (the insolation source; may or may not be the tidal host)
    # ------------------------------------------------------------------------------------------------------------------
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

    # ------------------------------------------------------------------------------------------------------------------
    # Orbital elements about the tidal host (per orbiting world; identify a world by index, name, or object)
    # ------------------------------------------------------------------------------------------------------------------
    def set_semi_major_axis(self, world, double semi_major_axis):
        """Set an orbiting world's semi-major axis about the host [m]."""
        self._system.get().set_semi_major_axis(<size_t>self._resolve_index(world), semi_major_axis)

    def set_eccentricity(self, world, double eccentricity):
        """Set an orbiting world's orbital eccentricity about the host."""
        self._system.get().set_eccentricity(<size_t>self._resolve_index(world), eccentricity)

    def get_semi_major_axis(self, world) -> float:
        """An orbiting world's semi-major axis about the host [m]."""
        return self._system.get().get_semi_major_axis(<size_t>self._resolve_index(world))

    def get_eccentricity(self, world) -> float:
        """An orbiting world's orbital eccentricity about the host."""
        return self._system.get().get_eccentricity(<size_t>self._resolve_index(world))

    def calc_gravitational_parameter(self, world) -> float:
        """Standard gravitational parameter ``mu = G (M_host + M_world)`` [m^3 s-2].

        Raises ``RuntimeError`` if no host is set; returns NaN for the host's own entry.
        """
        return self._system.get().calc_gravitational_parameter(<size_t>self._resolve_index(world))

    def calc_orbital_frequency(self, world) -> float:
        """Mean motion ``n = sqrt(mu / a^3)`` [rad s-1] for a world's two-body orbit about the host.

        Returns NaN for a non-positive semi-major axis or the host's own entry.
        """
        return self._system.get().calc_orbital_frequency(<size_t>self._resolve_index(world))

    def calc_semi_major_axis_from_frequency(self, world, double orbital_frequency) -> float:
        """Semi-major axis ``a = (mu / n^2)^(1/3)`` [m] from a mean motion (inverse of the above)."""
        return self._system.get().calc_semi_major_axis_from_frequency(
            <size_t>self._resolve_index(world), orbital_frequency)

    # ------------------------------------------------------------------------------------------------------------------
    # Orbital elements about the star (the insolation source; may differ from the tidal-host orbit)
    # ------------------------------------------------------------------------------------------------------------------
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

    # ------------------------------------------------------------------------------------------------------------------
    # Insolation (stellar irradiation of a world)
    # ------------------------------------------------------------------------------------------------------------------
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

        Grey-body radiative balance using the world's albedo + emissivity,
        ``T = ((1-A) F / (4 eps sigma))^(1/4)``. Returns NaN if the insolation flux is unavailable.
        """
        return self._system.get().calc_equilibrium_temperature(<size_t>self._resolve_index(world))

    # ------------------------------------------------------------------------------------------------------------------
    # Orbital + spin evolution (single-body tidal dissipation)
    # ------------------------------------------------------------------------------------------------------------------
    def calc_world_evolution(self, world) -> dict:
        """Evolve one orbiting world for a single tidal solve, returning its rates as a dict.

        Solves the world's global tides in the current system state (mean motion from Kepler's third law,
        spin + obliquity from the world, eccentricity + semi-major axis from its orbit about the host,
        host mass from the host world), then turns the tidal-potential derivatives into the orbital rates
        and the world's spin rate. Only this world raises tides; the host is treated as a point mass.

        Parameters
        ----------
        world : int or str or BaseWorld
            The orbiting world, identified by index, name, or the world object.

        Returns
        -------
        dict
            Keys (all MKS): 
                - ``world_index``
                - ``evolved`` (``False`` for the host or a world with no usable orbit)
                - ``orbital_frequency``
                - ``semi_major_axis``
                - ``eccentricity``
                - ``spin_frequency``
                - ``host_mass``
                - ``target_mass``
                - ``tidal_heating``
                - ``dU_dM``
                - ``dU_dw``
                - ``dU_dO``
                - ``da_dt``
                - ``de_dt``
                - ``dn_dt``
                - ``dspin_dt``
                - ``moment_of_inertia``
                - ``has_spin``
                - ``dE_orbit_dt``
                - ``dE_spin_dt``
                - ``energy_residual``
        """
        cdef c_WorldEvolution evolution = self._system.get().calc_world_evolution(
            <size_t>self._resolve_index(world))
        return _evolution_to_dict(evolution)

    def calc_system_evolution(self) -> list:
        """Evolve every world in the system (single-body dissipation).

        Solves each orbiting world's tides for its current orbit, then returns the rates.

        Returns
        -------
        list of dict
            One :meth:`calc_world_evolution` dict per world, in index order. The host's own entry and any
            world without a usable orbit come back with ``evolved`` set to ``False``.
        """
        cdef vector[c_WorldEvolution] results = self._system.get().calc_system_evolution()
        cdef list out = []
        cdef size_t i
        for i in range(results.size()):
            out.append(_evolution_to_dict(results[i]))
        return out

    def calc_pair_evolution(self, world) -> dict:
        """Evolve an orbiting world together with its host under dual-body tidal dissipation.

        Both the orbiting world and the tidal host raise a tide on their shared orbit. Each body's tides
        are solved with the other body as the tide raiser; their orbital-rate contributions add and each
        body evolves its own spin. A body with no tide model is rigid and contributes nothing (with a
        rigid host this reduces to :meth:`calc_world_evolution`).

        Parameters
        ----------
        world : int or str or BaseWorld
            The orbiting world, identified by index, name, or the world object.

        Returns
        -------
        dict
            Combined shared-orbit fields (``da_dt``, ``de_dt``, ``dn_dt``, ``tidal_heating_total``,
            ``dE_orbit_dt``, ``dE_spin_dt_total``, ``energy_residual``, plus ``orbital_frequency`` /
            ``semi_major_axis`` / ``eccentricity`` / ``world_index`` / ``host_index`` / ``evolved``), and
            each body's full single-body contribution under keys ``world`` and ``host`` (each a
            :meth:`calc_world_evolution`-style dict). ``evolved`` is ``False`` for the host's own entry, a
            hostless system, or a world with no usable orbit.
        """
        cdef c_PairEvolution pair = self._system.get().calc_pair_evolution(<size_t>self._resolve_index(world))
        return _pair_to_dict(pair)

    # ------------------------------------------------------------------------------------------------------------------
    # Config / serialization
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def config(self):
        """The system configuration dict the system was built from (None if built directly)."""
        return self.source_config

    cpdef dict get_config_dict(self):
        """Return the system's live state as a configuration dict (the ``build_system`` schema).

        Each member world is inlined with its own configuration (its retained ``config`` when present,
        else its world-level ``get_config_dict``) under the world's system name, together with its
        host/star roles and its orbital elements about the tidal host and the star. This is the
        self-contained, live-state expansion; it is used by :meth:`save_to_toml` when no
        ``source_config`` was retained (a directly-built system).

        Returns
        -------
        dict
            A system configuration dict with ``name`` and a ``worlds`` table.
        """
        cdef c_System* system_ptr = self._system.get()
        cdef int host_index = system_ptr.get_host_index()
        cdef int star_index = system_ptr.get_star_index()
        cdef int i
        cdef double a, stellar_a
        worlds_table = {}
        for i in range(<int>system_ptr.get_num_worlds()):
            world = self._world_wrappers[i]
            world_cfg = dict(world.config) if world.config is not None else world.get_config_dict()
            world_cfg["name"] = world.name
            entry = {"world": world_cfg}
            if i == host_index:
                entry["is_host"] = True
            if i == star_index:
                entry["is_star"] = True
            a = system_ptr.get_semi_major_axis(<size_t>i)
            if isfinite(a):
                entry["semi_major_axis_m"] = a
                entry["eccentricity"] = system_ptr.get_eccentricity(<size_t>i)
            stellar_a = system_ptr.get_stellar_semi_major_axis(<size_t>i)
            if isfinite(stellar_a):
                entry["stellar_semi_major_axis_m"] = stellar_a
                entry["stellar_eccentricity"] = system_ptr.get_stellar_eccentricity(<size_t>i)
            worlds_table[world.name] = entry
        return {"name": self.name, "worlds": worlds_table}

    def save_to_toml(self, str file_path, overwrite=True):
        """Write this system's configuration to a TOML file.

        Uses the retained build configuration (:attr:`source_config`) when present for a faithful
        round-trip of the original description (world references intact), otherwise falls back to the
        self-contained live-state expansion :meth:`get_config_dict`.

        Parameters
        ----------
        file_path : str
            Destination ``.toml`` path.
        overwrite : bool, optional
            Overwrite an existing file. Default True.
        """
        from TidalPy.structures_x.configs.config_writer import save_system_to_toml
        config = self.source_config if self.source_config is not None else self.get_config_dict()
        return save_system_to_toml(config, file_path, overwrite=overwrite)

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
            self._world_wrappers.append(_wrap_world(system_ptr.get_world(i)))

    def load_binary(self, str path, bint force=False):
        """Load this system's state from a TidalPy binary file (overriding the base to rewrap worlds).

        Reads the container state and rebuilds the heterogeneous world list from the stream (each world's
        concrete type is recovered from its record), then rebuilds the Python world wrappers around the
        loaded worlds. Physics sub-models a world does not serialize (the star's luminosity model, layer
        EOS data, tide/spin models) are reattached after load, as for a directly-loaded world.

        Parameters
        ----------
        path : str
            Source file path.
        force : bool, optional
            If True, attempt to load even on a schema-version mismatch.

        Raises
        ------
        FileNotFoundError
            If the file does not exist.
        IOError
            If the file is invalid or the schema version is incompatible.
        """
        import os as _os
        if not _os.path.isfile(path):
            raise FileNotFoundError(f"No such file: '{path}'")
        try:
            self._system.get().load_binary(path.encode("utf-8"), force)
        except RuntimeError as exc:
            raise IOError(str(exc)) from exc
        self._rebuild_world_wrappers()
        self.source_config = None

    # ------------------------------------------------------------------------------------------------------------------
    # World identification: accept an index (int), a world name (str), or the world wrapper object.
    # ------------------------------------------------------------------------------------------------------------------
    cdef Py_ssize_t _resolve_index(self, object world) except *:
        cdef Py_ssize_t n = len(self._world_wrappers)
        cdef Py_ssize_t i
        cdef int found
        if isinstance(world, int):
            i = <Py_ssize_t>world
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

    # ------------------------------------------------------------------------------------------------------------------
    # Sequence protocol over the member worlds
    # ------------------------------------------------------------------------------------------------------------------
    def __len__(self):
        """Number of worlds, so ``len(system)`` works."""
        return len(self._world_wrappers)

    def __iter__(self):
        """Iterate the member worlds, so ``for world in system: ...`` works."""
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
