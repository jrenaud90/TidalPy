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

from libc.stdint cimport uint32_t
from libcpp cimport bool as cpp_bool
from libcpp.utility cimport move
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport c_TidalPyBaseClass
from TidalPy.structures_x.worlds.base cimport BaseWorld, c_BaseWorld, c_WorldConfig
from TidalPy.structures_x.layers.base cimport BaseLayer, c_BaseLayer, c_tidal_scale_method_name
from TidalPy.structures_x.layers.physics cimport PhysicsLayer, c_PhysicsLayer
from TidalPy.structures_x.layers.solidliquid cimport SolidLiquidLayer, c_SolidLiquidLayer
from TidalPy.structures_x.layers.gas cimport GasLayer, c_GasLayer


# Build the matching layer wrapper as a NON-owning view onto a layer the world owns, dispatched
# by the C++ layer's concrete class id. The view keeps the world alive (see BaseLayer._view).
cdef BaseLayer _wrap_layer_view(c_BaseLayer* ptr, object world):
    cdef uint32_t class_id = ptr.get_layer_class_id()
    if class_id == 101:
        return PhysicsLayer._view(<c_PhysicsLayer*>ptr, world)
    elif class_id == 102:
        return SolidLiquidLayer._view(<c_SolidLiquidLayer*>ptr, world)
    elif class_id == 103:
        return GasLayer._view(<c_GasLayer*>ptr, world)
    return BaseLayer._view(ptr, world)

# Pull in the out-of-line definition of c_LayeredWorld::calc_tides (and the heavy
# global-potential engine it uses) so it compiles into this one extension only.
cdef extern from "world_tides_.hpp" nogil:
    pass


# Copy a C++ vector[double] into a new 1D float64 ndarray.
cdef cnp.ndarray _vec_to_ndarray(const vector[double]& v):
    cdef Py_ssize_t n = <Py_ssize_t>v.size()
    cdef cnp.ndarray out = np.empty(n, dtype=np.float64)
    cdef double[::1] mv
    cdef Py_ssize_t i
    if n > 0:
        mv = out
        for i in range(n):
            mv[i] = v[i]
    return out

# ODEMethod, c_EOSSolution, and c_WorldEOSSolveConfig are provided by layered.pxd.

# Selectors for the vectorized real-valued radius getters (see _eval_real).
cdef enum:
    _KIND_DENSITY        = 0
    _KIND_GRAVITY        = 1
    _KIND_PRESSURE       = 2
    _KIND_SHEAR_MOD      = 3
    _KIND_BULK_MOD       = 4
    _KIND_SHEAR_VISC     = 5
    _KIND_BULK_VISC      = 6
    _KIND_PRE_SHEAR_MOD  = 7
    _KIND_PRE_BULK_MOD   = 8
    _KIND_PRE_SHEAR_VISC = 9
    _KIND_PRE_BULK_VISC  = 10

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
        # The layer set changed; drop the cached views so they rebuild on next access.
        self._layer_views = None
        self._layer_view_by_name = None

    @property
    def num_layers(self) -> int:
        """Number of layers in the world."""
        return self._layered_ptr.get_num_layers()

    cdef list _ensure_layer_views(self):
        """Build the per-layer view cache once (lazily); reuse it thereafter.

        The views are non-owning wrappers onto the world's stable C++ layers, so they are built
        a single time (on first access after the layers are added) and reused, rather than
        re-allocated on every access. ``add_layer`` invalidates the cache.
        """
        cdef size_t n, i
        cdef BaseLayer view
        if self._layer_views is None:
            n = self._layered_ptr.get_num_layers()
            self._layer_views = []
            self._layer_view_by_name = {}
            for i in range(n):
                view = _wrap_layer_view(self._layered_ptr.get_layer(i), self)
                self._layer_views.append(view)
                self._layer_view_by_name[view.name] = view
        return self._layer_views

    def get_layer(self, index: int) -> BaseLayer:
        """Return a wrapper around the layer at ``index`` (0 = innermost).

        The returned object is a **non-owning view**: the world still owns the C++ layer, so
        the view exposes the layer's full Cython API (the matching ``PhysicsLayer`` /
        ``SolidLiquidLayer`` / ``GasLayer`` / ``BaseLayer`` subclass) but must not outlive the
        world (it keeps a reference to the world to prevent that). The views are cached (built
        once), so repeated access is cheap. Negative indices count from the end. Raises
        ``IndexError`` if out of range.
        """
        cdef list views = self._ensure_layer_views()
        cdef Py_ssize_t n = len(views)
        cdef Py_ssize_t i = index
        if i < 0:
            i += n
        if i < 0 or i >= n:
            raise IndexError(f"layer index {index} out of range (world has {n} layers)")
        return views[i]

    @property
    def layers(self) -> list:
        """List of (non-owning) layer views, inner to outer (cached). See :meth:`get_layer`."""
        return list(self._ensure_layer_views())

    def __len__(self):
        """Number of layers, so ``len(world)`` works."""
        return self._layered_ptr.get_num_layers()

    def __iter__(self):
        """Iterate the layer views inner-to-outer, so ``for layer in world: ...`` works."""
        return iter(self._ensure_layer_views())

    def __getitem__(self, index):
        """Index or slice the layers: ``world[0]`` (a view) or ``world[:2]`` (a list of views).

        Integer indices accept negatives and raise ``IndexError`` out of range; a slice returns
        the corresponding list of views.
        """
        if isinstance(index, slice):
            return self._ensure_layer_views()[index]
        return self.get_layer(index)

    def __getattr__(self, name):
        """Resolve ``world.<layer_name>`` to that layer's view (after normal attribute lookup).

        Only consulted when ``name`` is not a real attribute/method, so defined members always
        win. Names starting with ``_`` are never treated as layers (so dunder/internal probes
        are not intercepted). The view is taken from the cache (built once).

        This __getattr__ is only a fallback that is called if the standard getattr fails to find 
        a member attribute. So accessing attributes like `<world>.calc_internal_heating` will
        still work.
        """
        if name.startswith("_") or self._layered_ptr == NULL:
            raise AttributeError(name)
        self._ensure_layer_views()
        view = self._layer_view_by_name.get(name)
        if view is not None:
            return view
        raise AttributeError(
            f"'{type(self).__name__}' object has no attribute or layer named '{name}'")

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

    # ------------------------------------------------------------------------------------------------------------------
    # Structure + viscoelastic profile queries (delegate to the containing layer).
    #
    # Every getter accepts a scalar radius [m] (returns a float / complex) OR a NumPy
    # array of radii (returns an array of the same shape). NaN where the EOS is
    # unsolved, the layer is geometry-only, or no rheology is attached.
    # ------------------------------------------------------------------------------------------------------------------
    cdef double _eval_real(self, int kind, double radius_m) noexcept nogil:
        if   kind == _KIND_DENSITY:        return self._layered_ptr.get_density(radius_m)
        elif kind == _KIND_GRAVITY:        return self._layered_ptr.get_gravity(radius_m)
        elif kind == _KIND_PRESSURE:       return self._layered_ptr.get_pressure(radius_m)
        elif kind == _KIND_SHEAR_MOD:      return self._layered_ptr.get_shear_modulus(radius_m)
        elif kind == _KIND_BULK_MOD:       return self._layered_ptr.get_bulk_modulus(radius_m)
        elif kind == _KIND_SHEAR_VISC:     return self._layered_ptr.get_shear_viscosity(radius_m)
        elif kind == _KIND_BULK_VISC:      return self._layered_ptr.get_bulk_viscosity(radius_m)
        elif kind == _KIND_PRE_SHEAR_MOD:  return self._layered_ptr.get_premelt_shear_modulus(radius_m)
        elif kind == _KIND_PRE_BULK_MOD:   return self._layered_ptr.get_premelt_bulk_modulus(radius_m)
        elif kind == _KIND_PRE_SHEAR_VISC: return self._layered_ptr.get_premelt_shear_viscosity(radius_m)
        elif kind == _KIND_PRE_BULK_VISC:  return self._layered_ptr.get_premelt_bulk_viscosity(radius_m)
        return 0.0

    def _apply_real(self, radius_m, int kind):
        # float -> float; np.ndarray -> np.ndarray (same shape, looped under nogil).
        cdef cnp.ndarray in_arr
        cdef cnp.ndarray out_arr
        cdef double[::1] flat_in
        cdef double[::1] flat_out
        cdef Py_ssize_t i, n
        if isinstance(radius_m, np.ndarray):
            in_arr   = np.ascontiguousarray(radius_m, dtype=np.float64)
            out_arr  = np.empty_like(in_arr)
            flat_in  = in_arr.reshape(-1)
            flat_out = out_arr.reshape(-1)
            n = flat_in.shape[0]
            with nogil:
                for i in range(n):
                    flat_out[i] = self._eval_real(kind, flat_in[i])
            return out_arr
        return self._eval_real(kind, <double>radius_m)

    def _apply_complex(self, radius_m, double frequency_rad_s, cpp_bool is_shear):
        # float -> complex; np.ndarray -> complex np.ndarray (same shape).
        cdef cnp.ndarray in_arr
        cdef cnp.ndarray out_arr
        cdef double[::1] flat_in
        cdef double complex[::1] flat_out
        cdef cpp_complex[double] value
        cdef Py_ssize_t i, n
        if isinstance(radius_m, np.ndarray):
            in_arr   = np.ascontiguousarray(radius_m, dtype=np.float64)
            out_arr  = np.empty_like(in_arr, dtype=np.complex128)
            flat_in  = in_arr.reshape(-1)
            flat_out = out_arr.reshape(-1)
            n = flat_in.shape[0]
            for i in range(n):
                if is_shear:
                    value = self._layered_ptr.calc_complex_shear_modulus(flat_in[i], frequency_rad_s)
                else:
                    value = self._layered_ptr.calc_complex_bulk_modulus(flat_in[i], frequency_rad_s)
                flat_out[i] = value.real() + 1j * value.imag()
            return out_arr
        if is_shear:
            value = self._layered_ptr.calc_complex_shear_modulus(<double>radius_m, frequency_rad_s)
        else:
            value = self._layered_ptr.calc_complex_bulk_modulus(<double>radius_m, frequency_rad_s)
        return complex(value.real(), value.imag())

    def get_density(self, radius_m):
        """Density [kg/m^3] at radius_m [m] (float or np.ndarray); NaN if unsolved."""
        return self._apply_real(radius_m, _KIND_DENSITY)

    def get_gravity(self, radius_m):
        """Gravitational acceleration [m/s^2] at radius_m [m] (float or np.ndarray)."""
        return self._apply_real(radius_m, _KIND_GRAVITY)

    def get_pressure(self, radius_m):
        """Pressure [Pa] at radius_m [m] (float or np.ndarray); NaN if unsolved."""
        return self._apply_real(radius_m, _KIND_PRESSURE)

    def get_shear_modulus(self, radius_m):
        """Post-melt static shear modulus [Pa] at radius_m [m] (float or np.ndarray)."""
        return self._apply_real(radius_m, _KIND_SHEAR_MOD)

    def get_bulk_modulus(self, radius_m):
        """Post-melt static bulk modulus [Pa] at radius_m [m] (float or np.ndarray)."""
        return self._apply_real(radius_m, _KIND_BULK_MOD)

    def get_shear_viscosity(self, radius_m):
        """Post-melt shear viscosity [Pa s] at radius_m [m] (float or np.ndarray)."""
        return self._apply_real(radius_m, _KIND_SHEAR_VISC)

    def get_bulk_viscosity(self, radius_m):
        """Post-melt bulk viscosity [Pa s] at radius_m [m] (float or np.ndarray)."""
        return self._apply_real(radius_m, _KIND_BULK_VISC)

    def get_premelt_shear_modulus(self, radius_m):
        """Pre-melt static shear modulus [Pa] at radius_m [m] (float or np.ndarray)."""
        return self._apply_real(radius_m, _KIND_PRE_SHEAR_MOD)

    def get_premelt_bulk_modulus(self, radius_m):
        """Pre-melt static bulk modulus [Pa] at radius_m [m] (float or np.ndarray)."""
        return self._apply_real(radius_m, _KIND_PRE_BULK_MOD)

    def get_premelt_shear_viscosity(self, radius_m):
        """Pre-melt shear viscosity [Pa s] at radius_m [m] (float or np.ndarray)."""
        return self._apply_real(radius_m, _KIND_PRE_SHEAR_VISC)

    def get_premelt_bulk_viscosity(self, radius_m):
        """Pre-melt bulk viscosity [Pa s] at radius_m [m] (float or np.ndarray)."""
        return self._apply_real(radius_m, _KIND_PRE_BULK_VISC)

    def calc_complex_shear_modulus(self, radius_m, double frequency_rad_s, cpp_bool recalc_eos=False):
        """Complex shear modulus [Pa] at radius_m [m] (float or np.ndarray) and frequency [rad/s].

        Applies the containing layer's shear rheology to the stored post-melt static
        modulus + viscosity. Solves the EOS first if it has not been solved (or if
        ``recalc_eos``). NaN+0j for a geometry-only layer or no rheology.
        """
        self._ensure_solved(recalc_eos)
        return self._apply_complex(radius_m, frequency_rad_s, True)

    def calc_complex_bulk_modulus(self, radius_m, double frequency_rad_s, cpp_bool recalc_eos=False):
        """Complex bulk modulus [Pa] at radius_m [m] (float or np.ndarray) and frequency [rad/s].

        Applies the containing layer's bulk rheology to the stored post-melt static
        modulus + viscosity. Solves the EOS first if it has not been solved (or if
        ``recalc_eos``). NaN+0j for a geometry-only layer or no rheology.
        """
        self._ensure_solved(recalc_eos)
        return self._apply_complex(radius_m, frequency_rad_s, False)

    # ------------------------------------------------------------------------------------------------------------------
    # Shorthand bundles (one call returns several profiles at once).
    # ------------------------------------------------------------------------------------------------------------------
    def get_static_viscoelastics(self, radius_m):
        """``(shear_modulus, shear_viscosity, bulk_modulus, bulk_viscosity)`` (post-melt) at radius_m.

        Each element is a float (scalar radius) or np.ndarray (array of radii).
        """
        return (self.get_shear_modulus(radius_m), self.get_shear_viscosity(radius_m),
                self.get_bulk_modulus(radius_m),  self.get_bulk_viscosity(radius_m))

    def get_state(self, radius_m):
        """All EOS-related profiles at radius_m as a dict (float or np.ndarray values)."""
        return {
            "density":         self.get_density(radius_m),
            "gravity":         self.get_gravity(radius_m),
            "pressure":        self.get_pressure(radius_m),
            "shear_modulus":   self.get_shear_modulus(radius_m),
            "shear_viscosity": self.get_shear_viscosity(radius_m),
            "bulk_modulus":    self.get_bulk_modulus(radius_m),
            "bulk_viscosity":  self.get_bulk_viscosity(radius_m),
        }

    # ------------------------------------------------------------------------------------------------------------------
    # calc_* variants: solve the EOS first if it is unsolved (or force_recalc), then read the profile.
    # ------------------------------------------------------------------------------------------------------------------
    def _ensure_solved(self, cpp_bool force_recalc):
        if force_recalc or not self._layered_ptr.get_eos_solved():
            self.solve_eos()

    def calc_density(self, radius_m, cpp_bool force_recalc=False):
        """Density [kg/m^3]; solves the EOS first if needed."""
        self._ensure_solved(force_recalc)
        return self.get_density(radius_m)

    def calc_gravity(self, radius_m, cpp_bool force_recalc=False):
        """Gravitational acceleration [m/s^2]; solves the EOS first if needed."""
        self._ensure_solved(force_recalc)
        return self.get_gravity(radius_m)

    def calc_pressure(self, radius_m, cpp_bool force_recalc=False):
        """Pressure [Pa]; solves the EOS first if needed."""
        self._ensure_solved(force_recalc)
        return self.get_pressure(radius_m)

    def calc_shear_modulus(self, radius_m, cpp_bool force_recalc=False):
        """Post-melt shear modulus [Pa]; solves the EOS first if needed."""
        self._ensure_solved(force_recalc)
        return self.get_shear_modulus(radius_m)

    def calc_bulk_modulus(self, radius_m, cpp_bool force_recalc=False):
        """Post-melt bulk modulus [Pa]; solves the EOS first if needed."""
        self._ensure_solved(force_recalc)
        return self.get_bulk_modulus(radius_m)

    def calc_shear_viscosity(self, radius_m, cpp_bool force_recalc=False):
        """Post-melt shear viscosity [Pa s]; solves the EOS first if needed."""
        self._ensure_solved(force_recalc)
        return self.get_shear_viscosity(radius_m)

    def calc_bulk_viscosity(self, radius_m, cpp_bool force_recalc=False):
        """Post-melt bulk viscosity [Pa s]; solves the EOS first if needed."""
        self._ensure_solved(force_recalc)
        return self.get_bulk_viscosity(radius_m)

    def calc_static_viscoelastics(self, radius_m, cpp_bool force_recalc=False):
        """``get_static_viscoelastics`` but solves the EOS first if needed."""
        self._ensure_solved(force_recalc)
        return self.get_static_viscoelastics(radius_m)

    def calc_state(self, radius_m, cpp_bool force_recalc=False):
        """``get_state`` but solves the EOS first if needed."""
        self._ensure_solved(force_recalc)
        return self.get_state(radius_m)

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
    # Love-number (radial) solve
    # ------------------------------------------------------------------------------------------------------------------
    def solve_love_numbers(
            self,
            double frequency_rad_s    = 1.0e-5,
            int    degree_l           = 2,
            bint   solve_tidal        = True,
            bint   use_prop_matrix    = False,
            int    core_model         = 0,
            bint   use_kamata         = True,
            bint   nondimensionalize  = True,
            double starting_radius    = 0.0,
            double start_radius_tol   = 1.0e-4,
            str    integration_method = 'DOP853',
            double rtol               = 1.0e-6,
            double atol               = 1.0e-10,
            bint   scale_rtols        = True,
            size_t max_num_steps      = 500000,
            size_t expected_size      = 500,
            size_t max_ram_MB         = 500,
            double max_step           = 0.0,
            bint   verbose            = False,
            bint   warnings           = True,
            double eos_rtol           = 1.0e-6,
            double eos_atol           = 1.0e-10,
            double eos_pressure_tol   = 1.0e-3,
            int    eos_max_iters      = 100) -> dict:
        """Solve for whole-planet tidal Love numbers using the shooting method.

        Requires :meth:`solve_eos` to have been called first.  For each radial
        slice the layer's attached rheology model is evaluated at
        ``frequency_rad_s`` to obtain the complex moduli; the structure ODE is
        re-integrated from those density/modulus profiles and the deformation ODEs
        are shot from the center to the surface to yield k, h, l.

        Parameters
        ----------
        frequency_rad_s : float, optional
            Tidal forcing frequency [rad/s]. Default 1e-5.
        degree_l : int, optional
            Harmonic degree. Default 2.
        solve_tidal : bool, optional
            Use tidal boundary conditions. Default True.
        use_prop_matrix : bool, optional
            Select the radial-solve method: ``False`` (default) uses the shooting
            method; ``True`` uses the propagation-matrix method. The propagation
            matrix is only valid for a single solid, static, incompressible layer;
            an incompatible world fails the solve gracefully (``love_success`` is
            ``False`` with a non-zero ``love_error_code``).
        core_model : int, optional
            Propagation-matrix core starting condition (0-4). Ignored by the
            shooting method. Default 0.
        use_kamata : bool, optional
            Use Kamata starting conditions near the center (shooting method only).
            Default True.
        nondimensionalize : bool, optional
            Non-dimensionalise the problem internally (recommended). Default True.
        starting_radius : float, optional
            Minimum radius [m] for the shooting start.  0 → auto. Default 0.
        start_radius_tol : float, optional
            Tolerance for the starting-radius search. Default 1e-4.
        integration_method : str, optional
            CyRK ODE method: ``'DOP853'`` (default), ``'RK45'``, ``'RK23'``.
        rtol, atol : float, optional
            Relative / absolute ODE tolerances. Default 1e-6 / 1e-10.
        scale_rtols : bool, optional
            Scale tolerances by layer type. Default True.
        max_num_steps : int, optional
            Maximum ODE steps. Default 500000.
        expected_size, max_ram_MB : int, optional
            CyRK memory hints. Default 500 / 500.
        max_step : float, optional
            Maximum ODE step size [m].  0 → auto. Default 0.
        verbose : bool, optional
            Print solver status messages. Default False.
        warnings : bool, optional
            Emit solver warnings. Default True.
        eos_rtol, eos_atol, eos_pressure_tol : float, optional
            Tolerances for the internal EOS re-solve. Default 1e-6 / 1e-10 / 1e-3.
        eos_max_iters : int, optional
            Maximum iterations for the EOS pressure loop. Default 100.

        Returns
        -------
        dict
            ``success`` (bool), ``love_number_k``, ``love_number_h``,
            ``love_number_l`` (complex).

        Raises
        ------
        ValueError
            If the EOS has not been solved or the frequency is too small.

        Assumptions
        -----------
        - Spherical symmetry; all quantities MKS.
        - ``solve_eos`` must be called first.
        """
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

        cdef c_LoveSolveConfig cfg
        cfg.frequency_rad_s   = frequency_rad_s
        cfg.degree_l          = degree_l
        cfg.solve_tidal       = <cpp_bool>solve_tidal
        cfg.use_prop_matrix   = <cpp_bool>use_prop_matrix
        cfg.core_model        = core_model
        cfg.use_kamata        = <cpp_bool>use_kamata
        cfg.nondimensionalize = <cpp_bool>nondimensionalize
        cfg.starting_radius   = starting_radius
        cfg.start_radius_tol  = start_radius_tol
        cfg.integration_method = ode_method
        cfg.rtol              = rtol
        cfg.atol              = atol
        cfg.scale_rtols       = <cpp_bool>scale_rtols
        cfg.max_num_steps     = max_num_steps
        cfg.expected_size     = expected_size
        cfg.max_ram_MB        = max_ram_MB
        cfg.max_step          = max_step
        cfg.verbose           = <cpp_bool>verbose
        cfg.warnings          = <cpp_bool>warnings
        cfg.eos_rtol          = eos_rtol
        cfg.eos_atol          = eos_atol
        cfg.eos_pressure_tol  = eos_pressure_tol
        cfg.eos_max_iters     = eos_max_iters

        with nogil:
            self._layered_ptr.solve_love_numbers(cfg)

        return self._build_love_result()

    def solve_love_numbers_supplied(
            self,
            cnp.ndarray[cnp.complex128_t, ndim=1] complex_shear_modulus not None,
            cnp.ndarray[cnp.complex128_t, ndim=1] complex_bulk_modulus not None,
            cnp.ndarray[cnp.float64_t, ndim=1] radius_array not None,
            double frequency_rad_s    = 1.0e-5,
            int    degree_l           = 2,
            bint   solve_tidal        = True,
            bint   use_prop_matrix    = False,
            int    core_model         = 0,
            bint   use_kamata         = True,
            bint   nondimensionalize  = True,
            double starting_radius    = 0.0,
            double start_radius_tol   = 1.0e-4,
            str    integration_method = 'DOP853',
            double rtol               = 1.0e-6,
            double atol               = 1.0e-10,
            bint   scale_rtols        = True,
            size_t max_num_steps      = 500000,
            size_t expected_size      = 500,
            size_t max_ram_MB         = 500,
            double max_step           = 0.0,
            bint   verbose            = False,
            bint   warnings           = True) -> dict:
        """Solve Love numbers from externally-supplied complex moduli arrays (instead of layer rheology).

        The supplied shear/bulk moduli [Pa] are defined at ``radius_array`` [m] and are linearly interpolated onto
        the world's internal EOS radius grid. Used by the standalone ``RadialSolver_x.radial_solver`` API.
        ``solve_eos`` must be called first.
        """
        if (complex_shear_modulus.shape[0] != radius_array.shape[0]
                or complex_bulk_modulus.shape[0] != radius_array.shape[0]):
            raise ValueError("complex moduli and radius arrays must have matching length")

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

        cdef c_LoveSolveConfig cfg
        cfg.frequency_rad_s    = frequency_rad_s
        cfg.degree_l           = degree_l
        cfg.solve_tidal        = <cpp_bool>solve_tidal
        cfg.use_prop_matrix    = <cpp_bool>use_prop_matrix
        cfg.core_model         = core_model
        cfg.use_kamata         = <cpp_bool>use_kamata
        cfg.nondimensionalize  = <cpp_bool>nondimensionalize
        cfg.starting_radius    = starting_radius
        cfg.start_radius_tol   = start_radius_tol
        cfg.integration_method = ode_method
        cfg.rtol               = rtol
        cfg.atol               = atol
        cfg.scale_rtols        = <cpp_bool>scale_rtols
        cfg.max_num_steps      = max_num_steps
        cfg.expected_size      = expected_size
        cfg.max_ram_MB         = max_ram_MB
        cfg.max_step           = max_step
        cfg.verbose            = <cpp_bool>verbose
        cfg.warnings           = <cpp_bool>warnings

        cdef size_t n_in = radius_array.shape[0]
        cdef cpp_complex[double]* shear_ptr = <cpp_complex[double]*><void*>&complex_shear_modulus[0]
        cdef cpp_complex[double]* bulk_ptr  = <cpp_complex[double]*><void*>&complex_bulk_modulus[0]
        cdef double* radius_ptr = &radius_array[0]
        with nogil:
            self._layered_ptr.solve_love_numbers_supplied(cfg, shear_ptr, bulk_ptr, radius_ptr, n_in)
        return self._build_love_result()

    def _build_love_result(self):
        """Assemble the Python result dict from the retained C++ Love-number solution."""
        cdef cpp_complex[double] k, h, l
        k = self._layered_ptr.get_love_number_k(0)
        h = self._layered_ptr.get_love_number_h(0)
        l = self._layered_ptr.get_love_number_l(0)
        return {
            'success':       bool(self._layered_ptr.get_love_solved()),
            'error_code':    self._layered_ptr.get_love_error_code(),
            'message':       self._layered_ptr.get_love_message().decode('utf-8'),
            'love_number_k': complex(k.real(), k.imag()),
            'love_number_h': complex(h.real(), h.imag()),
            'love_number_l': complex(l.real(), l.imag()),
        }

    # ------------------------------------------------------------------------------------------------------------------
    # Love-number solve state (mirrors EOS getter pattern)
    # ------------------------------------------------------------------------------------------------------------------

    @property
    def love_solved(self) -> bool:
        """True once solve_love_numbers has completed successfully."""
        return bool(self._layered_ptr.get_love_solved())

    @property
    def love_success(self) -> bool:
        """True if the last love-number solve reported success."""
        return bool(self._layered_ptr.get_love_success())

    @property
    def love_error_code(self) -> int:
        """Integer error code from the last love-number solve (-100 if not yet run)."""
        return self._layered_ptr.get_love_error_code()

    @property
    def love_message(self) -> str:
        """Status message from the last love-number solve."""
        return self._layered_ptr.get_love_message().decode('utf-8')

    @property
    def love_num_ytypes(self) -> int:
        """Number of boundary-condition types solved (typically 1 for tidal-only)."""
        return int(self._layered_ptr.get_love_num_ytypes())

    @property
    def love_number_k(self) -> complex:
        """Complex potential Love number k2 from the last radial solve (NaN+0j if unsolved)."""
        cdef cpp_complex[double] v = self._layered_ptr.get_love_number_k(<size_t>0)
        return complex(v.real(), v.imag())

    @property
    def love_number_h(self) -> complex:
        """Complex radial displacement Love number h2 from the last radial solve (NaN+0j if unsolved)."""
        cdef cpp_complex[double] v = self._layered_ptr.get_love_number_h(<size_t>0)
        return complex(v.real(), v.imag())

    @property
    def love_number_l(self) -> complex:
        """Complex tangential (Shida) Love number l2 from the last radial solve (NaN+0j if unsolved)."""
        cdef cpp_complex[double] v = self._layered_ptr.get_love_number_l(<size_t>0)
        return complex(v.real(), v.imag())

    def get_love_number_k(self, ytype_idx: int = 0) -> complex:
        """Complex k Love number for the given boundary-condition ytype index."""
        cdef cpp_complex[double] v = self._layered_ptr.get_love_number_k(<size_t>ytype_idx)
        return complex(v.real(), v.imag())

    def get_love_number_h(self, ytype_idx: int = 0) -> complex:
        """Complex h Love number for the given boundary-condition ytype index."""
        cdef cpp_complex[double] v = self._layered_ptr.get_love_number_h(<size_t>ytype_idx)
        return complex(v.real(), v.imag())

    def get_love_number_l(self, ytype_idx: int = 0) -> complex:
        """Complex l (Shida) Love number for the given boundary-condition ytype index."""
        cdef cpp_complex[double] v = self._layered_ptr.get_love_number_l(<size_t>ytype_idx)
        return complex(v.real(), v.imag())

    def get_love_surface_y(self, ytype_idx: int, y_idx: int) -> complex:
        """Complex radial y-solution value at the surface for the given ytype and y index."""
        cdef cpp_complex[double] v = self._layered_ptr.get_love_surface_y(
            <size_t>ytype_idx, <size_t>y_idx)
        return complex(v.real(), v.imag())

    # ------------------------------------------------------------------------------------------------------------------
    # Global (1D) tidal dissipation
    #
    # The tide model holder, config, and the analytic result accessors are inherited from
    # BaseWorld (set_tide_model / set_tide_config / tide_model_set / tides_solved /
    # get_tidal_heating / get_tidal_potential_derivatives / get_num_tidal_modes /
    # get_tidal_love_k). A layered world overrides calc_tides to add the rheology (radial-
    # solver) path + per-layer heating, and adds the per-layer heating accessor.
    # ------------------------------------------------------------------------------------------------------------------
    def calc_tides(
            self,
            double orbital_frequency,
            double spin_frequency,
            double eccentricity,
            double obliquity,
            double semi_major_axis,
            double host_mass):
        """Solve the global tidal dissipation for the given orbital/spin state.

        Requires an attached tide model (:meth:`set_tide_model`). Populates the world's
        :attr:`tidal_heating`, the three potential derivatives, and per-layer heating.

        The analytic models (cpl/ctl/ctl_q) collapse from their fixed per-degree parameters.
        The rheology model instead runs the world radial solver at each unique tidal
        frequency to find the global Love numbers, so the EOS must be solved first
        (:meth:`solve_eos`).

        Parameters
        ----------
        orbital_frequency : float
            Orbital mean motion [rad s-1].
        spin_frequency : float
            Spin rate of the deformed body [rad s-1].
        eccentricity : float
            Orbital eccentricity [dimensionless].
        obliquity : float
            Axial tilt [radians].
        semi_major_axis : float
            Orbital semi-major axis [m].
        host_mass : float
            Mass of the tidal host [kg].

        Raises
        ------
        RuntimeError
            If no tide model is attached, the rheology model is selected but the EOS has
            not been solved, a radial-solver Love-number solve fails, or the global
            potential solve fails.
        """
        cdef c_TideSolveConfig state
        state.orbital_frequency = orbital_frequency
        state.spin_frequency    = spin_frequency
        state.eccentricity      = eccentricity
        state.obliquity         = obliquity
        state.semi_major_axis   = semi_major_axis
        state.host_mass         = host_mass
        self._layered_ptr.calc_tides(state)

    def get_layer_tidal_heating(self, index: int) -> float:
        """Tidal heating [W] deposited in layer ``index`` (= world heating × its effective scale)."""
        return self._layered_ptr.get_layer_tidal_heating(<size_t>index)

    # ------------------------------------------------------------------------------------------------------------------
    # On-demand 3D tidal stress/strain/heating
    # ------------------------------------------------------------------------------------------------------------------
    def get_3d_tidal_heating(
            self,
            double orbital_frequency,
            double spin_frequency,
            double eccentricity,
            double obliquity,
            double semi_major_axis,
            double host_mass,
            double radius,
            double colatitude) -> float:
        """Secular (cycle/orbit-averaged) 3D tidal volumetric heating [W m-3] at ``(radius, colatitude)``.

        This is the physically time-averaged power density: the active tidal modes are built dynamically
        from the world's ``[tides]`` truncation config (max degree l, eccentricity/obliquity truncation),
        the world radial response is solved once per mode, and each mode contributes
        ``(omega/2) Im(sigma_c : conj(eps_c))`` (a single ``omega/2``, complex amplitudes, no ``abs``),
        summed with sign. It is longitude- and time-independent (the ``e^{i m phi}`` cancels per mode),
        and its volume integral over the planet equals the 1D global tidal heating
        (:meth:`get_tidal_heating`). Requires the rheology tide model (:meth:`set_tide_model`) and a
        solved EOS (:meth:`solve_eos`). Returns NaN in liquid layers / at the center / below the solver's
        starting radius.
        """
        cdef c_TideSolveConfig state
        state.orbital_frequency = orbital_frequency
        state.spin_frequency    = spin_frequency
        state.eccentricity      = eccentricity
        state.obliquity         = obliquity
        state.semi_major_axis   = semi_major_axis
        state.host_mass         = host_mass
        return self._layered_ptr.get_3d_tidal_heating(state, radius, colatitude)

    def get_3d_tidal_heating_array(
            self,
            double orbital_frequency,
            double spin_frequency,
            double eccentricity,
            double obliquity,
            double semi_major_axis,
            double host_mass,
            radii,
            colatitudes):
        """Secular 3D tidal volumetric heating [W m-3] at an array of ``(radius, colatitude)`` points.

        Vectorized batch form of :meth:`get_3d_tidal_heating`: ``radii`` and ``colatitudes`` are paired,
        equal-length 1D arrays (point ``i`` is ``(radii[i], colatitudes[i])``), and a same-shape
        ``np.ndarray`` of heating is returned. Same physics and preconditions as the scalar method, but
        the world radial response is solved once per unique ``(degree l, frequency)`` and reused across
        all points (the solve does not depend on radius or colatitude), so this is the efficient way to
        build a heating map. Entries are NaN for points in a liquid layer / at the center / below the
        solver's starting radius. Requires the rheology tide model and a solved EOS.
        """
        cdef cnp.ndarray radii_arr = np.ascontiguousarray(radii, dtype=np.float64)
        cdef cnp.ndarray colat_arr = np.ascontiguousarray(colatitudes, dtype=np.float64)
        if radii_arr.shape[0] != colat_arr.shape[0]:
            raise ValueError("radii and colatitudes must have the same length")

        cdef size_t num_points = <size_t>radii_arr.shape[0]
        cdef cnp.ndarray out_arr = np.empty(num_points, dtype=np.float64)
        if num_points == 0:
            return out_arr

        cdef double[::1] radii_view = radii_arr
        cdef double[::1] colat_view = colat_arr
        cdef double[::1] out_view   = out_arr

        cdef c_TideSolveConfig state
        state.orbital_frequency = orbital_frequency
        state.spin_frequency    = spin_frequency
        state.eccentricity      = eccentricity
        state.obliquity         = obliquity
        state.semi_major_axis   = semi_major_axis
        state.host_mass         = host_mass

        with nogil:
            self._layered_ptr.get_3d_tidal_heating_array(
                state, &radii_view[0], &colat_view[0], num_points, &out_view[0])
        return out_arr

    def calc_3d_tides(
            self,
            double orbital_frequency,
            double spin_frequency,
            double eccentricity,
            double obliquity,
            double semi_major_axis,
            double host_mass,
            radii=None,
            colatitudes=None,
            longitudes=None,
            times=None,
            orbit_averaged=True,
            latitude_summed=False,
            longitude_summed=False,
            radial_summed=False,
            int latitude_nodes=16,
            int longitude_nodes=64,
            int radial_slices=16) -> dict:
        """3D tidal heating as a full grid over ``(radius, colatitude, longitude[, time])`` or reduced.

        With ``orbit_averaged=True`` (default) the quantity is the secular (cycle-averaged) volumetric
        heating density ``h_bar`` [W m-3], which is longitude- and time-independent (so the longitude
        axis is constant and there is no time axis). With ``orbit_averaged=False`` it is the instantaneous
        mechanical power density ``sigma_ij(t) * eps_dot_ij(t)`` [W m-3] at each supplied time (a 4th
        axis), which orbit-averages back to ``h_bar``.

        Any spatial dimension can be integrated out via ``latitude_summed`` / ``longitude_summed`` /
        ``radial_summed``. Reduction convention: if any spatial axis is summed, the surviving spatial
        axes carry their Jacobian (``r^2``, ``sin theta``, ``1``) so a plain integral over them recovers
        the total; if none is summed the output is the raw density. The colatitude integral uses an
        internal Gauss-Legendre grid (``latitude_nodes``), the radial integral an internal per-layer
        trapezoid (``radial_slices``), and the longitude integral the analytic ``2*pi`` when averaged or
        a ``longitude_nodes`` trapezoid when instantaneous.

        Non-summed spatial axes require the corresponding ``radii`` / ``colatitudes`` / ``longitudes``
        arrays; the ``times`` array is required when ``orbit_averaged=False``. The returned dict carries
        the surviving axes (``radii``, ``colatitudes``, ``longitudes``, and ``times`` when instantaneous)
        plus either ``heating`` (the grid over the surviving axes, in the order radius, colatitude,
        longitude, time) or, when all three spatial axes are summed, ``total`` [W] and ``per_layer`` [W]
        (per-layer totals, innermost first; each an array over time when instantaneous). Requires the
        rheology tide model and a solved EOS.
        """
        cdef cpp_bool instantaneous = not orbit_averaged
        cdef c_Heating3DCollapseConfig cfg
        cfg.orbit_averaged   = <cpp_bool>orbit_averaged
        cfg.latitude_summed  = <cpp_bool>latitude_summed
        cfg.longitude_summed = <cpp_bool>longitude_summed
        cfg.radial_summed    = <cpp_bool>radial_summed
        cfg.latitude_nodes   = latitude_nodes
        cfg.longitude_nodes  = longitude_nodes
        cfg.radial_slices    = radial_slices

        cdef cnp.ndarray radii_arr
        cdef double[::1] radii_view
        cdef const double* radii_ptr = NULL
        cdef size_t num_radii = 0
        if not radial_summed:
            if radii is None:
                raise ValueError("radii must be provided when radial_summed is False")
            radii_arr = np.ascontiguousarray(radii, dtype=np.float64)
            num_radii = <size_t>radii_arr.shape[0]
            if num_radii > 0:
                radii_view = radii_arr
                radii_ptr = &radii_view[0]

        cdef cnp.ndarray colat_arr
        cdef double[::1] colat_view
        cdef const double* colat_ptr = NULL
        cdef size_t num_colat = 0
        if not latitude_summed:
            if colatitudes is None:
                raise ValueError("colatitudes must be provided when latitude_summed is False")
            colat_arr = np.ascontiguousarray(colatitudes, dtype=np.float64)
            num_colat = <size_t>colat_arr.shape[0]
            if num_colat > 0:
                colat_view = colat_arr
                colat_ptr = &colat_view[0]

        cdef cnp.ndarray lon_arr
        cdef double[::1] lon_view
        cdef const double* lon_ptr = NULL
        cdef size_t num_lon = 0
        if not longitude_summed:
            if longitudes is None:
                raise ValueError("longitudes must be provided when longitude_summed is False")
            lon_arr = np.ascontiguousarray(longitudes, dtype=np.float64)
            num_lon = <size_t>lon_arr.shape[0]
            if num_lon > 0:
                lon_view = lon_arr
                lon_ptr = &lon_view[0]

        cdef cnp.ndarray time_arr
        cdef double[::1] time_view
        cdef const double* time_ptr = NULL
        cdef size_t num_time = 0
        if instantaneous:
            if times is None:
                raise ValueError("times must be provided when orbit_averaged is False")
            time_arr = np.ascontiguousarray(times, dtype=np.float64)
            num_time = <size_t>time_arr.shape[0]
            if num_time > 0:
                time_view = time_arr
                time_ptr = &time_view[0]

        cdef c_TideSolveConfig state
        state.orbital_frequency = orbital_frequency
        state.spin_frequency    = spin_frequency
        state.eccentricity      = eccentricity
        state.obliquity         = obliquity
        state.semi_major_axis   = semi_major_axis
        state.host_mass         = host_mass

        cdef c_Heating3DCollapsed res
        with nogil:
            res = self._layered_ptr.calc_3d_tides(
                state, radii_ptr, num_radii, colat_ptr, num_colat,
                lon_ptr, num_lon, time_ptr, num_time, cfg)

        # Surviving-axis shape (row-major, order radius, colatitude, longitude, time).
        cdef list shape = [<Py_ssize_t>res.shape[i] for i in range(res.shape.size())]
        cdef dict out = {'radii': _vec_to_ndarray(res.radii),
                         'colatitudes': _vec_to_ndarray(res.colatitudes),
                         'longitudes': _vec_to_ndarray(res.longitudes)}
        if instantaneous:
            out['times'] = _vec_to_ndarray(res.times)

        cdef Py_ssize_t nlayers = <Py_ssize_t>res.n_layers
        cdef Py_ssize_t ntimes = <Py_ssize_t>res.n_times
        if res.all_spatial_summed:
            total = _vec_to_ndarray(res.values)
            per_layer = _vec_to_ndarray(res.layer_totals).reshape(nlayers, ntimes)
            if not instantaneous:
                total = float(total[0]) if total.size else 0.0
                per_layer = per_layer[:, 0]
            out['total'] = total
            out['per_layer'] = per_layer
        else:
            out['heating'] = _vec_to_ndarray(res.values).reshape(shape) if shape else _vec_to_ndarray(res.values)
        return out

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
        cdef bytes method_bytes
        layers = []
        for i in range(n):
            layer_ptr = self._layered_ptr.get_layer(i)
            method_bytes = c_tidal_scale_method_name(layer_ptr.get_tidal_scale_method())
            layers.append({
                "name":               layer_ptr.get_name().decode("utf-8"),
                "layer_index":        layer_ptr.get_layer_index(),
                "radius_inner_m":     layer_ptr.get_radius_inner(),
                "radius_outer_m":     layer_ptr.get_radius_outer(),
                "mass_kg":            layer_ptr.get_mass(),
                "material_name":      layer_ptr.get_material_name().decode("utf-8"),
                "is_tidal":           True if layer_ptr.get_is_tidal() else False,
                "tidal_scale":        layer_ptr.get_tidal_scale(),
                "tidal_scale_method": method_bytes.decode("utf-8"),
            })
        d["num_layers"] = n
        d["layers"] = layers
        return d
