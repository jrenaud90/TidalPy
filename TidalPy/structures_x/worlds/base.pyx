# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrapper for TidalPy's base world class.

BaseWorld holds world-level identity and the orbital and thermal scalars (albedo, emissivity, obliquity, spin
frequency) plus bulk geometry and equilibrium-temperature calculations. Layered worlds and stars subclass it.
"""

from libcpp cimport bool as cpp_bool
from libcpp.utility cimport move
from libcpp.memory cimport make_shared
from libcpp.complex cimport complex as cpp_complex

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport (
    StructureBase,
    c_TidalPyBaseClass,
    c_PhysicsBase,
    cy_physics_model_config,
)
from TidalPy.Tides_x.classes.tide cimport TideBase
from TidalPy.Tides_x.love.love cimport c_parse_love_method_int, c_love_method_name_int
from TidalPy.Tides_x.eccentricity.eccentricity_driver import (
    eccentricity_truncation_name, validate_eccentricity_exact_tolerance, validate_eccentricity_truncation)
from TidalPy.Tides_x.obliquity.obliquity_driver import obliquity_truncation_name, validate_obliquity_truncation

# Pull in the out-of-line definition of c_BaseWorld::calc_tides (the analytic global tidal
# path) plus the heavy global-potential engine it uses, so they compile into this extension.
cdef extern from "world_tides_base_.hpp" nogil:
    pass

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())

# World ``type`` values the TOML builder accepts (kept in step with configs.toml_loader.WORLD_TYPES, which
# cannot be imported here at module load without a circular import).
BUILDER_WORLD_TYPES = ("star", "gasgiant", "terrestrial", "layered")


cdef class BaseWorld(StructureBase):
    """Base world: identity, orbital/thermal scalars, and bulk geometry.

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
        Bond albedo [dimensionless]. Default ``0.3``.
    emissivity : float, optional
        Surface emissivity [dimensionless]. Default ``1.0``.
    obliquity : float, optional
        Axial obliquity [rad]. Default ``0.0``.
    spin_frequency : float, optional
        Rotation rate [rad/s]. Default ``0.0``.

    Assumptions
    -----------
    - Spherically symmetric world.
    """

    def __cinit__(self, *args, **kwargs):
        pass  # unique_ptr<c_BaseWorld> auto-inits to nullptr; _ptr set in __init__

    def __init__(
            self,
            str    name,
            double radius,
            double mass,
            str    world_type = "world",
            double albedo     = 0.3,
            double emissivity = 1.0,
            double obliquity  = 0.0,
            double spin_frequency = 0.0):
        cdef c_WorldConfig config
        config.name           = name.encode("utf-8")
        config.world_type_str = world_type.encode("utf-8")
        config.radius     = radius
        config.mass       = mass
        config.albedo     = albedo
        config.emissivity = emissivity
        config.obliquity  = obliquity
        config.spin_frequency = spin_frequency
        # _world_ptr is a shared_ptr (a System co-owns the world), so build it with make_shared.
        self._world_ptr = make_shared[c_BaseWorld](config)
        self._ptr = <c_TidalPyBaseClass*>self._world_ptr.get()

    def __dealloc__(self):
        self._world_ptr.reset()
        self._ptr = NULL

    @staticmethod
    cdef BaseWorld _wrap(shared_ptr[c_BaseWorld] ptr):
        """Wrap an already-constructed C++ base world (no new C++ object is built)."""
        cdef BaseWorld world = BaseWorld.__new__(BaseWorld)
        world._world_ptr = ptr
        world._ptr = <c_TidalPyBaseClass*>ptr.get()
        return world

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

    def calc_equilibrium_temperature(self, double insolation_flux) -> float:
        """Radiative-equilibrium temperature [K] for a given insolation flux.

        T_eq = [ (1 − A)·F / (4·ε·σ) ]^(1/4)

        Parameters
        ----------
        insolation_flux : float
            Incident stellar flux [W/m^2].

        Returns
        -------
        float
            Equilibrium temperature [K]; 0.0 for non-positive flux.

        Assumptions
        -----------
        - Fast-rotator, uniform-temperature surface.
        """
        return self._world_ptr.get().calc_equilibrium_temperature(insolation_flux)

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
    def set_tide_model(self, TideBase tide not None):
        """Attach a global tide dissipation model (transfers ownership).

        Ownership of the C++ model moves out of ``tide``, which is left an empty shell and must not be reused.

        On a layerless world (e.g. a star) only the analytic models (``cpl``/``ctl``/``ctl_q``)
        are usable; the ``rheology`` model needs the radial solver and a layered world.
        """
        if tide._tide_ptr.get() == NULL:
            raise ValueError("This tide model holds no C++ object (already attached or moved).")
        self._world_ptr.get().set_tide_model(move(tide._tide_ptr))
        tide._ptr = NULL

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
            love_fixed_dt=None):
        """Change the stored ``[tides]`` truncation/degree configuration and the world's Love-number method.

        Only the arguments given change; every other setting keeps its current value (see
        :meth:`get_tide_config`), so a call can adjust one setting without resetting the rest.

        Parameters
        ----------
        min_degree_l, max_degree_l : int, optional
            Tidal harmonic degree range (2..10).
        eccentricity_truncation : int, optional
            Eccentricity-function truncation level N: every product of two eccentricity functions is kept through
            e^N. Tabulated levels: ``TidalPy.Tides_x.eccentricity.ECCENTRICITY_TRUNCATIONS`` (2, 4, 6, 8, 10, 20,
            50), or ``"exact"`` for the functions from the exact orbit (any e < 1).
        eccentricity_exact_tolerance : float, optional
            Heating tail tolerance in (0, 1) that sets the mode range of ``"exact"`` (ignored by the levels).
        obliquity_truncation : int or str, optional
            Obliquity truncation level N: every product of two obliquity functions is kept through I^N. Tabulated
            levels: ``TidalPy.Tides_x.obliquity.OBLIQUITY_TRUNCATIONS`` (0 or ``"off"``, 2, 4), or ``"gen"`` for
            the general functions (any obliquity).
        layer_tidal_heating : bool, optional
            Whether ``calc_tides`` also resolves each layer's heating when the Love numbers come from the radial
            solver (a volume integral of the radial solution that costs about as much as the global solve again);
            the other paths share out the heating at no extra cost. Default ``True``.
        love_method : str, optional
            How the world obtains Love numbers when its tide model asks for them (and the default for
            ``solve_love_numbers``): ``'radial_solver'`` (``'shooting'``, ``'rs'``), ``'propagation_matrix'``
            (``'prop_matrix'``, ``'pm'``, ``'prop'``), ``'homogeneous'`` (``'homogen'``), ``'cpl'``, ``'ctl'``,
            or ``'laterally_inhomogeneous'`` (``'3d'``, ``'lat_inhom'``; reserved, not implemented).
        love_fixed_q, love_fixed_dt : float, optional
            Quality factor for the ``'cpl'`` method and time lag [s] for the ``'ctl'`` method. A NaN clears the
            value, after which the attached tide model's per-degree fixed Q / time lag is used.
        """
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
            cfg.love_method = c_parse_love_method_int(str(love_method).encode('utf-8'))
        if love_fixed_q is not None:
            cfg.love_fixed_q = <double>love_fixed_q
        if love_fixed_dt is not None:
            cfg.love_fixed_dt = <double>love_fixed_dt
        self._world_ptr.get().set_tide_config(cfg)

    def calc_tides(
            self,
            double orbital_frequency,
            double spin_frequency,
            double eccentricity,
            double obliquity,
            double semi_major_axis,
            double host_mass):
        """Solve the global tidal dissipation for the given orbital and spin state.

        Requires an attached tide model (:meth:`set_tide_model`) and populates :attr:`tidal_heating` and the
        three potential derivatives. The base world runs the analytic models (cpl, ctl, ctl_q) only;
        :class:`LayeredWorld` adds the rheology path and per-layer heating.

        Raises
        ------
        RuntimeError
            If no tide model is attached, the rheology model is selected on a non-layered world, or the global
            potential solve fails.
        """
        cdef c_TideSolveConfig state
        state.orbital_frequency = orbital_frequency
        state.spin_frequency    = spin_frequency
        state.eccentricity      = eccentricity
        state.obliquity         = obliquity
        state.semi_major_axis   = semi_major_axis
        state.host_mass         = host_mass
        self._world_ptr.get().calc_tides(state)

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

        A new tide model or tide config clears it, and so does a layered world's ``solve_eos``: the result
        describes the structure it was solved with.
        """
        return self._world_ptr.get().get_tides_solved()

    def get_tidal_heating(self) -> float:
        """Total global tidal heating [W] (NaN if unsolved)."""
        return self._world_ptr.get().get_tidal_heating()

    def get_tidal_potential_derivatives(self) -> tuple:
        """The three orbital potential derivatives ``(dUdM, dUdw, dUdO)`` [J kg-1 rad-1]."""
        return (
            self._world_ptr.get().get_tidal_dU_dM(),
            self._world_ptr.get().get_tidal_dU_dw(),
            self._world_ptr.get().get_tidal_dU_dO(),
        )

    def get_num_tidal_modes(self) -> int:
        """Number of active (nonzero-frequency) tidal modes summed in the last solve."""
        return self._world_ptr.get().get_num_tidal_modes()

    def get_tidal_love_k(self, degree_l: int, m: int, p: int, q: int) -> complex:
        """Complex potential Love number ``k_l`` for the tidal mode ``(l, m, p, q)``.

        Only populated for the rheology tide model (a layered world). Returns NaN for the
        analytic models, which carry no displacement Love numbers, or for an inactive mode.
        """
        cdef cpp_complex[double] k = self._world_ptr.get().get_tidal_love_k(
            <int>degree_l, <int>m, <int>p, <int>q)
        return complex(k.real(), k.imag())

    @staticmethod
    def build(source, force=False):
        """Build a world from a configuration source (the public builder entry point).

        This is a factory: the concrete subclass returned (``LayeredWorld``, ``GasGiantWorld``, or
        ``StarWorld``) follows the configuration's world ``type``, whichever class ``build`` is called on. The
        normalized configuration is retained on :attr:`source_config` so the world can be written back to TOML.

        Parameters
        ----------
        source : str or dict
            A bundled world name, a path to a ``.toml`` file, or a configuration dict.
        force : bool, optional
            If True, bypass the schema-version compatibility warning. Default False.

        Returns
        -------
        BaseWorld
            The constructed world, with ``source_config`` populated.
        """
        # Deferred imports: the builder helpers import the world subclasses, so
        # importing them at module load would be circular.
        import os
        from TidalPy.structures_x.configs.world_builder import (
            _resolve_source,
            construct_world)
        from TidalPy.structures_x.configs.toml_loader import (
            load_toml,
            merge_with_defaults,
            validate_schema_version)
        from TidalPy.structures_x.configs.worldpack import resolve_data_file

        # `_resolve_source` hands back whatever the caller gave (a path string, a Path, or an already-parsed
        # mapping), so this one stays `object`; the rest have a single concrete type.
        cdef object resolved = _resolve_source(source)
        cdef dict config = load_toml(resolved)
        cdef object given_data_file = None
        cdef str base_dir
        cdef BaseWorld world
        # Checked before the defaults fill a missing version in, so a file without one says so.
        validate_schema_version(config, force=force)
        config = merge_with_defaults(config)
        # Resolve a companion data file (e.g. a PREM profile) relative to the world
        # file's directory so construct_world can open it directly.
        if "data_file" in config:
            given_data_file = config["data_file"]
            base_dir = os.path.dirname(resolved) if isinstance(resolved, str) else None
            config["data_file"] = resolve_data_file(config["data_file"], base_dir)
        world = construct_world(config)
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
        class family's default applies (``layered``, ``gasgiant``, or ``star``).
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

    def save_to_toml(self, str file_path, overwrite=True):
        """Write this world's configuration to a TOML file.

        Writes :meth:`get_config_dict`, the world as it is now (a change made after the build, a new obliquity
        say, is saved), validated against the world schema first so the file builds or this raises
        ``ValueError``. A world built from a ``data_file`` writes :attr:`portable_config` instead: the file
        reference as given and the tables that refined its layers, so the saved file builds anywhere the data
        file resolves (changes made to that world after the build are not saved). The file starts with a comment header naming
        the TidalPy, SciPy, and CyRK versions that wrote it.

        Parameters
        ----------
        file_path : str
            Destination ``.toml`` path.
        overwrite : bool, optional
            Overwrite an existing file. Default True.
        """
        from TidalPy.structures_x.configs.config_writer import save_world_to_toml
        cdef dict config
        if self.portable_config is not None:
            config = self.portable_config
        else:
            from TidalPy.structures_x.configs.toml_loader import validate_world_config
            config = self.get_config_dict()
            validate_world_config(config)
        return save_world_to_toml(config, file_path, overwrite=overwrite)

    cpdef dict get_config_dict(self):
        """Return the world configuration as the TOML builder's world table (MKS).

        The dict validates against the world schema and carries a ``tides`` table when a tide model is attached
        (``global_tidal_model`` plus the model's per-degree parameters and the stored degree and truncation
        settings). Subclasses add their layers or stellar values.

        Returns
        -------
        dict
            Keys: ``schema_version``, ``name``, ``type``, ``radius``, ``mass``, ``albedo``, ``emissivity``,
            ``obliquity``, ``spin_frequency``, and ``tides`` when set.
        """
        from TidalPy.structures_x.configs.toml_loader import SCHEMA_VERSION
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
        return config
