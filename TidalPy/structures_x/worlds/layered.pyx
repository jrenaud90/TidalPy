# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
layered.pyx
Cython/Python wrapper for TidalPy's layered world class.

LayeredWorld: a world built from an ordered stack of layers (inner to outer).
Owns its layers, aggregates total mass and internal radiogenic heating, and
validates layer-boundary continuity. Whole-planet EOS / radial solves are added
as methods on this class in later phases.
"""

cimport numpy as cnp
cnp.import_array()

import numpy as np

from libcpp cimport bool as cpp_bool
from libcpp.utility cimport move
from cython.operator cimport dereference as deref

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport c_TidalPyBaseClass
from TidalPy.structures_x.worlds.base cimport BaseWorld, c_BaseWorld, c_WorldConfig
from TidalPy.structures_x.layers.base cimport BaseLayer, c_BaseLayer

# ODEMethod, c_EOSSolution, and c_WorldEOSSolveConfig are provided by layered.pxd.

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())

# =====================================================================================================================
# LayeredWorld
# =====================================================================================================================

cdef class LayeredWorld(BaseWorld):
    """A world built from an ordered (inner-to-outer) stack of layers.

    Parameters
    ----------
    name : str
        Human-readable world name.
    radius_m : float
        World radius [m].
    mass_kg : float
        World mass [kg].
    world_type : str, optional
        Type label. Default ``"terrestrial"``.
    albedo, emissivity, obliquity_rad, spin_frequency_rad_s : float, optional
        See :class:`BaseWorld`.

    Notes
    -----
    Add layers inner-to-outer with :meth:`add_layer`; each layer's inner radius
    must match the previous layer's outer radius (the innermost starts at 0).
    """

    def __cinit__(self, *args, **kwargs):
        self._layered_ptr = NULL

    def __init__(
            self,
            str    name,
            double radius_m,
            double mass_kg,
            str    world_type           = "terrestrial",
            double albedo               = 0.3,
            double emissivity           = 1.0,
            double obliquity_rad        = 0.0,
            double spin_frequency_rad_s = 0.0):
        cdef c_WorldConfig config
        config.name                 = name.encode("utf-8")
        config.world_type_str       = world_type.encode("utf-8")
        config.radius_m             = radius_m
        config.mass_kg              = mass_kg
        config.albedo               = albedo
        config.emissivity           = emissivity
        config.obliquity_rad        = obliquity_rad
        config.spin_frequency_rad_s = spin_frequency_rad_s
        cdef c_LayeredWorld* raw = new c_LayeredWorld(config)
        self._world_ptr.reset(<c_BaseWorld*>raw)
        self._layered_ptr = raw
        self._ptr         = <c_TidalPyBaseClass*>raw

    def __dealloc__(self):
        self._layered_ptr = NULL  # base's unique_ptr owns the C++ object

    # ------------------------------------------------------------------------------------------------------------------
    # Layer management
    # ------------------------------------------------------------------------------------------------------------------
    def add_layer(self, BaseLayer layer not None):
        """Add a layer to the world (inner to outer).

        Ownership of the C++ layer (and its attached physics models) is
        transferred into the world; the passed ``BaseLayer`` becomes an empty,
        non-owning shell and must not be reused.

        Parameters
        ----------
        layer : BaseLayer
            A layer (``BaseLayer``, ``PhysicsLayer``, ``SolidLiquidLayer``, or
            ``GasLayer``). Its inner radius must match the current outermost
            radius (0 for the first layer).

        Raises
        ------
        ValueError
            If ``layer`` has already been added/moved.
        RuntimeError
            If the layer geometry is not continuous with the existing stack.
        """
        if layer._layer_ptr.get() == NULL:
            raise ValueError(
                "This layer holds no C++ object (already added to a world or moved).")
        # Validate continuity before transferring ownership so a rejected layer
        # stays usable (the C++ add_layer would otherwise consume it on throw).
        if not self._layered_ptr.accepts_layer(deref(layer._layer_ptr.get())):
            raise ValueError(
                "Layer geometry is not continuous: the inner radius does not match "
                "the current outermost radius (add layers inner-to-outer, innermost "
                "starting at radius 0).")
        self._layered_ptr.add_layer(move(layer._layer_ptr))

    @property
    def num_layers(self) -> int:
        """Number of layers in the world."""
        return self._layered_ptr.get_num_layers()

    def calc_total_mass(self) -> float:
        """Total mass [kg] = sum of all layer masses."""
        return self._layered_ptr.calc_total_mass()

    def calc_internal_heating(self, double time_s) -> float:
        """Total internal radiogenic heating [W] at the given time [s].

        Only ``SolidLiquidLayer`` layers with an attached radiogenics model
        contribute; all other layers contribute zero.
        """
        return self._layered_ptr.calc_internal_heating(time_s)

    def validate_layers(self) -> bool:
        """True if every layer boundary is continuous (innermost starts at 0)."""
        return self._layered_ptr.validate_layers()

    # ------------------------------------------------------------------------------------------------------------------
    # Equation of state
    # ------------------------------------------------------------------------------------------------------------------
    def solve_eos(
            self,
            double surface_pressure   = 0.0,
            size_t slices_per_layer   = 100,
            double G_to_use           = -1.0,
            str    integration_method = 'DOP853',
            double rtol               = 1.0e-6,
            double atol               = 1.0e-10,
            double pressure_tol       = 1.0e-3,
            size_t max_iters          = 100,
            double temperature        = 0.0,
            bint   verbose            = False) -> dict:
        """Solve the whole-planet equation of state.

        Integrates gravity, pressure, enclosed mass, and moment of inertia
        radially from the planet center to its surface, using each layer's
        attached material EOS model (see :meth:`BaseLayer.set_eos`) as the local
        density source. A convergence loop on the surface pressure determines the
        central pressure. On success, every layer's EOS profile is populated so
        that :meth:`get_density`, :meth:`get_gravity`, and :meth:`get_pressure`
        (on this world or on the individual layers) become available.

        Parameters
        ----------
        surface_pressure : float, optional
            Target surface pressure [Pa]. Default 0.0.
        slices_per_layer : int, optional
            Number of radial sample points generated per layer (>= 2).
            Default 100.
        G_to_use : float, optional
            Gravitational constant [m^3 kg^-1 s^-2]. If negative (default), the
            TidalPy config value is used.
        integration_method : str, optional
            CyRK integration method: ``'DOP853'`` (default), ``'RK45'``, or
            ``'RK23'``.
        rtol, atol : float, optional
            Relative / absolute integration tolerances. Default 1e-6 / 1e-10.
        pressure_tol : float, optional
            Surface-pressure convergence tolerance. Default 1e-3.
        max_iters : int, optional
            Maximum convergence iterations. Default 100.
        temperature : float, optional
            Temperature [K] passed to each EOS model's ``calc_density`` (unused
            by the current isothermal models). Default 0.0.
        verbose : bool, optional
            Print solver status messages. Default False.

        Returns
        -------
        dict
            ``success``, ``message``, ``iterations``, ``pressure_error``, the
            radial profile arrays (``radius``, ``gravity``, ``pressure``,
            ``mass``, ``moi``, ``density``), and the scalar surface/planet
            results (``surface_gravity``, ``surface_pressure``,
            ``central_pressure``, ``planet_mass``, ``planet_moi``).

        Raises
        ------
        ValueError
            If the world has no layers, any layer lacks a material EOS model, an
            unsupported integration method is given, or ``slices_per_layer < 2``.

        Assumptions
        -----------
        - Spherical symmetry; all quantities MKS.
        - Each layer's density comes from its attached material EOS model.

        Notes
        -----
        The whole solve (radius-grid generation, bulk-density estimate, the CyRK
        ODE integration with its surface-pressure loop, and populating each
        layer's profile) runs in pure C++ (``c_LayeredWorld::solve_eos``); this
        wrapper only translates the integration-method string to the CyRK enum and
        builds the Python result dict from the retained C++ solution.
        """
        # Translate the integration-method string to the CyRK enum (string handling
        # stays at the Cython/Python boundary).
        cdef ODEMethod ode_method
        cdef str method_upper = integration_method.upper()
        if method_upper == 'DOP853':
            ode_method = ODEMethod.DOP853
        elif method_upper == 'RK45':
            ode_method = ODEMethod.RK45
        elif method_upper == 'RK23':
            ode_method = ODEMethod.RK23
        else:
            raise ValueError(f"Unsupported integration method: {integration_method}")

        cdef c_WorldEOSSolveConfig cfg
        cfg.surface_pressure   = surface_pressure
        cfg.slices_per_layer   = slices_per_layer
        cfg.G_to_use           = G_to_use
        cfg.integration_method = ode_method
        cfg.rtol               = rtol
        cfg.atol               = atol
        cfg.pressure_tol       = pressure_tol
        cfg.max_iters          = max_iters
        cfg.temperature        = temperature
        cfg.verbose            = <cpp_bool>verbose

        # Pure-C++ solve. Input validation throws std::invalid_argument, surfaced
        # here as ValueError via the ``except +`` on the C++ declaration.
        with nogil:
            self._layered_ptr.solve_eos(cfg)

        return self._build_eos_result()

    def _build_eos_result(self):
        """Assemble the Python result dict from the retained C++ EOS solution."""
        cdef const c_EOSSolution* sol = self._layered_ptr.get_eos_solution()
        cdef size_t n = sol.radius_array_size if sol != NULL else 0
        cdef cnp.ndarray radius_out   = np.empty(n, dtype=np.float64)
        cdef cnp.ndarray gravity_out  = np.empty(n, dtype=np.float64)
        cdef cnp.ndarray pressure_out = np.empty(n, dtype=np.float64)
        cdef cnp.ndarray mass_out     = np.empty(n, dtype=np.float64)
        cdef cnp.ndarray moi_out      = np.empty(n, dtype=np.float64)
        cdef cnp.ndarray density_out  = np.empty(n, dtype=np.float64)
        cdef size_t j
        if sol != NULL and self._layered_ptr.get_eos_solved():
            for j in range(n):
                radius_out[j]   = sol.radius_array_vec[j]
                gravity_out[j]  = sol.gravity_array_vec[j]
                pressure_out[j] = sol.pressure_array_vec[j]
                mass_out[j]     = sol.mass_array_vec[j]
                moi_out[j]      = sol.moi_array_vec[j]
                density_out[j]  = sol.density_array_vec[j]

        return {
            'success':          self._layered_ptr.get_eos_success(),
            'message':          self._layered_ptr.get_eos_message().decode('utf-8'),
            'iterations':       self._layered_ptr.get_eos_iterations(),
            'pressure_error':   self._layered_ptr.get_eos_pressure_error(),
            'radius':           radius_out,
            'gravity':          gravity_out,
            'pressure':         pressure_out,
            'mass':             mass_out,
            'moi':              moi_out,
            'density':          density_out,
            'surface_gravity':  self._layered_ptr.get_surface_gravity_eos(),
            'surface_pressure': self._layered_ptr.get_surface_pressure_eos(),
            'central_pressure': self._layered_ptr.get_central_pressure(),
            'planet_mass':      self._layered_ptr.get_planet_mass_eos(),
            'planet_moi':       self._layered_ptr.get_planet_moi_eos(),
        }

    @property
    def eos_solved(self) -> bool:
        """True once the world-level EOS solve has populated the layer profiles."""
        return self._layered_ptr.get_eos_solved()

    @property
    def all_eos_set(self) -> bool:
        """True once every layer has a material EOS model attached."""
        return self._layered_ptr.get_all_eos_set()

    def get_density(self, double radius_m) -> float:
        """Density [kg/m^3] at radius_m [m] from the solved EOS profile (NaN if unsolved)."""
        return self._layered_ptr.get_density(radius_m)

    def get_gravity(self, double radius_m) -> float:
        """Gravitational acceleration [m/s^2] at radius_m [m] (NaN if EOS unsolved)."""
        return self._layered_ptr.get_gravity(radius_m)

    def get_pressure(self, double radius_m) -> float:
        """Pressure [Pa] at radius_m [m] from the solved EOS profile (NaN if unsolved)."""
        return self._layered_ptr.get_pressure(radius_m)

    @property
    def surface_gravity_eos(self) -> float:
        """Surface gravity [m/s^2] from the last EOS solve, or NaN if not solved."""
        return self._layered_ptr.get_surface_gravity_eos()

    @property
    def central_pressure(self) -> float:
        """Central pressure [Pa] from the last EOS solve, or NaN if not solved."""
        return self._layered_ptr.get_central_pressure()

    @property
    def planet_mass_eos(self) -> float:
        """Total planet mass [kg] integrated by the last EOS solve, or NaN if not solved."""
        return self._layered_ptr.get_planet_mass_eos()

    @property
    def planet_moi_eos(self) -> float:
        """Planet moment of inertia [kg m^2] from the last EOS solve, or NaN if not solved."""
        return self._layered_ptr.get_planet_moi_eos()

    # ------------------------------------------------------------------------------------------------------------------
    # Config
    # ------------------------------------------------------------------------------------------------------------------
    cpdef dict get_config_dict(self):
        """Return the world config plus a list of per-layer base config dicts.

        Returns
        -------
        dict
            All :class:`BaseWorld` keys plus ``num_layers`` and ``layers`` (a
            list of geometry-level dicts for each layer, in index order).
        """
        d = BaseWorld.get_config_dict(self)
        cdef size_t n = self._layered_ptr.get_num_layers()
        cdef size_t i
        cdef c_BaseLayer* layer_ptr
        layers = []
        for i in range(n):
            layer_ptr = self._layered_ptr.get_layer(i)
            layers.append({
                "name":           layer_ptr.get_name().decode("utf-8"),
                "layer_index":    layer_ptr.get_layer_index(),
                "radius_inner_m": layer_ptr.get_radius_inner(),
                "radius_outer_m": layer_ptr.get_radius_outer(),
                "mass_kg":        layer_ptr.get_mass(),
                "material_name":  layer_ptr.get_material_name().decode("utf-8"),
                "is_tidal":       True if layer_ptr.get_is_tidal() else False,
                "tidal_scale":    layer_ptr.get_tidal_scale(),
            })
        d["num_layers"] = n
        d["layers"] = layers
        return d
