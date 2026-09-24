# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrapper for TidalPy's layered world class.

LayeredWorld owns an ordered stack of layers (inner to outer), aggregates total mass and internal radiogenic
heating, and validates layer-boundary continuity. The whole-planet EOS and radial (Love number) solves are
methods on this class.
"""

cimport numpy as cnp
cnp.import_array()

import weakref

import numpy as np

from libc.stdint cimport uint32_t
from libc.stdlib cimport malloc, free
from libcpp cimport bool as cpp_bool
from libcpp.utility cimport move
from libcpp.memory cimport make_unique, static_pointer_cast
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address, d_PI, d_NAN
from TidalPy.Utilities_x.classes_x.classes cimport c_TidalPyBaseClass
from TidalPy.structures_x.worlds.base cimport BaseWorld, c_BaseWorld, c_WorldConfig
from TidalPy.structures_x.layers.base cimport (
    BaseLayer, c_BaseLayer, c_layer_class_name, cy_eos_field, cy_eos_fields, C_EOS_DENSITY_INDEX,
    C_EOS_GRAVITY_INDEX, C_EOS_PRESSURE_INDEX, C_EOS_SHEAR_MODULUS_INDEX, C_EOS_SHEAR_VISCOSITY_INDEX,
    C_EOS_BULK_MODULUS_INDEX, C_EOS_BULK_VISCOSITY_INDEX, C_EOS_TEMPERATURE_INDEX, C_EOS_HEAT_FLOW_INDEX,
    C_EOS_MELT_FRACTION_INDEX)
from TidalPy.structures_x.layers.base import LAYER_STANDALONE_CONFIG_KEYS
from TidalPy.structures_x.layers.physics cimport PhysicsLayer, c_PhysicsLayer
from TidalPy.structures_x.layers.solidliquid cimport SolidLiquidLayer, c_SolidLiquidLayer
from TidalPy.structures_x.layers.gas cimport GasLayer, c_GasLayer
from TidalPy.RadialSolver_x.rs_constants cimport C_MAX_NUM_YTYPES
from TidalPy.RadialSolver_x.rs_solution cimport RadialSolverSolution
from TidalPy.RadialSolver_x.rs_solution cimport cy_check_surface_solve_conditioning
from TidalPy.Tides_x.love.love cimport (
    c_parse_love_method_int, c_love_method_name_int, c_love_method_uses_radial_solver_int)
from TidalPy.Utilities_x.logging_x.logger import log_warning


# Build the matching layer wrapper as a non-owning view onto a layer the world owns, dispatched by the C++
# layer's concrete class id. The view keeps the world alive (see BaseLayer._view).
cdef BaseLayer cy_wrap_layer_view(c_BaseLayer* ptr, object world):
    cdef bytes class_name = c_layer_class_name(ptr.get_layer_class_id())
    if class_name == b"physics":
        return PhysicsLayer._view(<c_PhysicsLayer*>ptr, world)
    elif class_name == b"solidliquid":
        return SolidLiquidLayer._view(<c_SolidLiquidLayer*>ptr, world)
    elif class_name == b"gas":
        return GasLayer._view(<c_GasLayer*>ptr, world)
    return BaseLayer._view(ptr, world)

# Pull in the out-of-line definition of c_LayeredWorld::calc_tides (and the heavy
# global-potential engine it uses) so it compiles into this one extension only.
cdef extern from "world_tides_.hpp" nogil:
    pass

# Component order of the stress and strain grids returned by LayeredWorld.calc_3d_stress_strain.
STRESS_STRAIN_COMPONENTS = ("rr", "theta_theta", "phi_phi", "r_theta", "r_phi", "theta_phi")


# Copy the first `n` values of a C++ vector[double] into a typed view, NaN-filling any tail. A solution's
# profile vectors are filled point by point by interpolate_full_planet and can stop short of the radius array,
# so the copy is bounded by both lengths.
cdef void cy_fill_from_vec(double[::1] out, const vector[double]& source, Py_ssize_t n) noexcept nogil:
    cdef Py_ssize_t i
    cdef Py_ssize_t m = <Py_ssize_t>source.size()
    if m > n:
        m = n
    for i in range(m):
        out[i] = source[i]
    for i in range(m, n):
        out[i] = d_NAN


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

# The solve_eos result dict from a report copied under the world's call lock (c_LayeredWorld::get_eos_report).
cdef dict cy_eos_report_to_dict(const c_WorldEOSReport& report):
    cdef size_t j
    cdef size_t num_layers = report.layer_thermal.size()
    cdef list layer_temperature        = []
    cdef list layer_heat_flow_in       = []
    cdef list layer_heat_flow_out      = []
    cdef list layer_heating            = []
    cdef list layer_node_temperature   = []
    cdef list layer_top_temperature    = []
    cdef list layer_boundary_thickness = []
    cdef list layer_rayleigh_number    = []
    cdef list layer_nusselt_number     = []
    cdef list layer_in_thermal_network = []
    for j in range(num_layers):
        layer_temperature.append(report.layer_thermal[j].temperature)
        layer_heat_flow_in.append(report.layer_thermal[j].heat_flow_in)
        layer_heat_flow_out.append(report.layer_thermal[j].heat_flow_out)
        layer_heating.append(report.layer_thermal[j].heating)
        layer_node_temperature.append(report.layer_thermal[j].node_temperature)
        layer_top_temperature.append(report.layer_thermal[j].top_temperature)
        layer_boundary_thickness.append(report.layer_thermal[j].boundary_thickness)
        layer_rayleigh_number.append(report.layer_thermal[j].rayleigh_number)
        layer_nusselt_number.append(report.layer_thermal[j].nusselt_number)
        layer_in_thermal_network.append(bool(report.layer_thermal[j].in_network))

    return {
        'success':          bool(report.success),
        'message':          report.message.decode('utf-8'),
        'iterations':       report.iterations,
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
        'geometry_converged':     bool(report.geometry_converged),
        'layer_radius_outer':     list(report.layer_radius_outer),
        'layer_temperature':      layer_temperature,
        'layer_heat_flow_in':     layer_heat_flow_in,
        'layer_heat_flow_out':    layer_heat_flow_out,
        'layer_heating':          layer_heating,
        'layer_temperature_rate': list(report.layer_temperature_rate),
        'layer_node_temperature':   layer_node_temperature,
        'layer_top_temperature':    layer_top_temperature,
        'layer_boundary_thickness': layer_boundary_thickness,
        'layer_rayleigh_number':    layer_rayleigh_number,
        'layer_nusselt_number':     layer_nusselt_number,
        'layer_in_thermal_network': layer_in_thermal_network,
    }


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
    """Build a :class:`LayeredWorld` whose layers interpolate their own slice of a radial profile.

    The layers are built in C++ (``c_build_world_from_layered_profile``), which is also what the standalone
    ``RadialSolver_x.radial_solver`` calls, so a world built here and the temporary one that API solves are
    the same world. Interface radii appear twice in the profile, once as the top of the lower layer and once
    as the base of the upper one, and each copy belongs to its own layer.

    The world carries its layers and their material EOS models and nothing else: no tide model, no
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
    LayeredWorld

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

    cdef shared_ptr[c_LayeredWorld] world_sptr
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

    return LayeredWorld._wrap(static_pointer_cast[c_BaseWorld, c_LayeredWorld](world_sptr))


cdef int cy_check_num_threads(int num_threads) except -1:
    """Raise ValueError unless ``num_threads`` is at least 1."""
    if num_threads < 1:
        raise ValueError(f"num_threads must be at least 1; got {num_threads}")
    return 0

# Translate an integration-method name to the CyRK enum (string handling stays at the Cython boundary).
# Case-insensitive; covers the explicit Runge-Kutta methods and the implicit, stiff ones.
cdef str cy_integration_method_name(ODEMethod method):
    """The configuration spelling of a CyRK integration method."""
    if method == ODEMethod.DOP853:
        return "DOP853"
    elif method == ODEMethod.RK45:
        return "RK45"
    elif method == ODEMethod.RK23:
        return "RK23"
    elif method == ODEMethod.BDF:
        return "BDF"
    elif method == ODEMethod.LSODA:
        return "LSODA"
    elif method == ODEMethod.RADAU:
        return "Radau"
    raise ValueError(f"Unsupported integration method code: {<int>method}.")


cdef ODEMethod cy_resolve_integration_method(str integration_method) except *:
    cdef str method_upper = integration_method.upper()
    if method_upper == 'DOP853':
        return ODEMethod.DOP853
    elif method_upper == 'RK45':
        return ODEMethod.RK45
    elif method_upper == 'RK23':
        return ODEMethod.RK23
    elif method_upper == 'BDF':
        return ODEMethod.BDF
    elif method_upper == 'LSODA':
        return ODEMethod.LSODA
    elif method_upper == 'RADAU':
        return ODEMethod.RADAU
    raise ValueError(
        f"Unsupported integration method: {integration_method}. "
        "Supported: RK23, RK45, DOP853, BDF, LSODA, Radau.")


cdef void cy_apply_love_solve_overrides(
        c_LoveSolveConfig* cfg,
        object use_kamata,
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
    if use_kamata is not None:
        cfg.use_kamata = <cpp_bool>bool(use_kamata)
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


cdef int cy_resolve_love_method(str love_method) except? -999:
    # Map a Love-number method name (or alias) to its c_LoveMethod index.
    cdef int method = c_parse_love_method_int(love_method.encode('utf-8'))
    if method == 5:
        raise NotImplementedError(
            "The laterally_inhomogeneous Love-number method is reserved for the 3D Love solver and is not "
            "implemented.")
    return method

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
    (<const c_LayeredWorld*>owner).get_eos_fields(field_indices, num_fields, radii, num_radii, values_out)


cdef class LayeredWorld(BaseWorld):
    """A world built from an ordered (inner-to-outer) stack of layers.

    Parameters
    ----------
    name : str
        World name.
    radius : float
        World radius [m].
    mass : float
        World mass [kg].
    world_type : str, optional
        Type label. Default ``"terrestrial"``.
    albedo, emissivity, obliquity, spin_frequency : float, optional
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
            double radius,
            double mass,
            str    world_type = "terrestrial",
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
        # make_unique owns the allocation; ownership then moves into the base-typed member
        # (Cython cannot assign a unique_ptr[Derived] to a unique_ptr[Base] directly).
        cdef unique_ptr[c_LayeredWorld] built = make_unique[c_LayeredWorld](config)
        self._layered_ptr = built.get()
        self._world_ptr.reset(<c_BaseWorld*>built.release())
        self._ptr = <c_TidalPyBaseClass*>self._world_ptr.get()

    def __dealloc__(self):
        self._layered_ptr = NULL  # base's unique_ptr owns the C++ object

    @staticmethod
    cdef LayeredWorld _wrap(shared_ptr[c_BaseWorld] ptr):
        """Wrap an already-constructed C++ layered world (no new C++ object is built)."""
        cdef LayeredWorld world = LayeredWorld.__new__(LayeredWorld)
        world._world_ptr = ptr
        world._ptr = <c_TidalPyBaseClass*>ptr.get()
        world._layered_ptr = <c_LayeredWorld*>ptr.get()
        world._layer_views = None
        world._layer_view_by_name = None
        return world

    def add_layer(self, BaseLayer layer not None):
        """Add a layer to the world (inner to outer).

        Ownership of the C++ layer and its attached physics models moves to the world. ``layer`` stays usable: it
        becomes a non-owning view of the layer the world now holds, the same kind ``world.<layer name>`` returns.

        Parameters
        ----------
        layer : BaseLayer
            A layer (``BaseLayer``, ``PhysicsLayer``, ``SolidLiquidLayer``, or
            ``GasLayer``). Its inner radius must match the current outermost
            radius (0 for the first layer).

        Raises
        ------
        ValueError
            If ``layer`` has already been added or moved, is a view of another world's layer, does not continue
            the existing stack, reaches past the world radius, or shares its name with a layer already added.
        """
        if layer._layer_ptr.get() == NULL:
            raise ValueError(
                "This layer holds no C++ object (already added to a world or moved).")
        if layer._is_view:
            raise ValueError(
                "This layer is a non-owning view into a world-owned layer and cannot be "
                "added to another world. Construct a new layer instead.")
        # Validate continuity before transferring ownership so a rejected layer
        # stays usable (the C++ add_layer would otherwise consume it on throw).
        cdef string rejection = self._layered_ptr.layer_rejection_reason(deref(layer._layer_ptr.get()))
        if rejection.size() > 0:
            raise ValueError(rejection.decode('utf-8'))
        cdef c_BaseLayer* added_layer_ptr = layer._layer_ptr.get()
        self._layered_ptr.add_layer(move(layer._layer_ptr))
        layer._init_view(added_layer_ptr, self)
        self._track_view(layer)
        # The layer set changed; drop the cached views so they rebuild on next access.
        self._layer_views = None
        self._layer_view_by_name = None

    def load_binary(self, str path, cpp_bool force=False):
        """Load this world's state, its layers included, from a TidalPy binary file.

        The load replaces the world's layers, so layer views taken from this world before it (``world.<name>``,
        ``get_layer``, iteration) no longer refer to a layer of this world and must not be used; take new ones.
        Nothing solved survives the load: run ``solve_eos`` again.

        Parameters
        ----------
        path : str
            Source file path.
        force : bool, optional
            Attempt the load even on a schema version mismatch.
        """
        # Every view handed out points at a C++ layer the load replaces: detach them before the layers are freed,
        # so a held view raises instead of reading freed memory.
        cdef object view_ref
        cdef BaseLayer view
        if self._issued_views is not None:
            for view_ref in self._issued_views:
                view = view_ref()
                if view is not None:
                    view._detach()
            self._issued_views = []
        self._layer_views = None
        self._layer_view_by_name = None
        BaseWorld.load_binary(self, path, force)

    cdef void _track_view(self, BaseLayer view) except *:
        """Remember a view this world handed out (weakly), so a load that replaces the layers can detach it."""
        if self._issued_views is None:
            self._issued_views = []
        # Drop references to views that no longer exist, so the list stays as long as the live views.
        self._issued_views = [view_ref for view_ref in self._issued_views if view_ref() is not None]
        self._issued_views.append(weakref.ref(view))

    def _layer_moved(self):
        """Called by a layer view after it moved its layer's radii: the solved profile no longer lines up."""
        self._layered_ptr.update_after_layer_geometry_change()

    @property
    def num_layers(self) -> int:
        """Number of layers in the world."""
        return self._layered_ptr.get_num_layers()

    cdef list _ensure_layer_views(self):
        """Build the per-layer view cache once, lazily, and reuse it thereafter.

        The views are non-owning wrappers onto the world's stable C++ layers; ``add_layer`` invalidates them.
        """
        cdef size_t n, i
        cdef BaseLayer view
        if self._layer_views is None:
            n = self._layered_ptr.get_num_layers()
            self._layer_views = []
            self._layer_view_by_name = {}
            for i in range(n):
                view = cy_wrap_layer_view(self._layered_ptr.get_layer(i), self)
                self._track_view(view)
                self._layer_views.append(view)
                self._layer_view_by_name[view.name] = view
        return self._layer_views

    def get_layer(self, index: int) -> BaseLayer:
        """Return a wrapper around the layer at ``index`` (0 = innermost).

        The returned object is a non-owning view: the world still owns the C++ layer, so the view exposes the
        matching subclass API but must not outlive the world (it holds a reference to the world to prevent
        that). Views are cached, so repeated access is cheap. Negative indices count from the end; an index
        out of range raises ``IndexError``.
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

        Only consulted when ``name`` is not a real attribute or method, so defined members such as
        ``calc_internal_heating`` always win. Names starting with ``_`` are never treated as layers, so
        internal probes are not intercepted.
        """
        if name.startswith("_") or self._layered_ptr == NULL:
            raise AttributeError(name)
        self._ensure_layer_views()
        cdef object view = self._layer_view_by_name.get(name)
        if view is not None:
            return view
        raise AttributeError(
            f"'{type(self).__name__}' object has no attribute or layer named '{name}'")

    def calc_total_mass(self) -> float:
        """Total mass [kg] = sum of all layer masses."""
        return self._layered_ptr.calc_total_mass()

    def calc_internal_heating(self, double time) -> float:
        """Total internal radiogenic heating [W] at the given time [s].

        Only ``SolidLiquidLayer`` layers with an attached radiogenics model
        contribute; all other layers contribute zero.
        """
        return self._layered_ptr.calc_internal_heating(time)

    def validate_layers(self) -> bool:
        """True if every layer boundary is continuous (innermost starts at 0)."""
        return self._layered_ptr.validate_layers()

    def solve_eos(
            self,
            double surface_pressure = 0.0,
            slices_per_layer        = None,
            double G_to_use         = -1.0,
            integration_method      = None,
            rtol                    = None,
            atol                    = None,
            pressure_tol            = None,
            max_iters               = None,
            nondimensionalize       = None,
            temperature             = None,
            solve_temperature       = None,
            surface_temperature     = None,
            reset_layer_masses      = False,
            cpp_bool verbose        = False,
            time                    = None) -> dict:
        """Solve the whole-planet equation of state.

        Integrates gravity, pressure, enclosed mass, and moment of inertia from the planet center to its
        surface with each layer's attached material EOS model (see :meth:`BaseLayer.set_eos`) as the local
        density source; a convergence loop on the surface pressure sets the central pressure. On success every
        layer's EOS profile is populated, so :meth:`get_density`, :meth:`get_gravity`, and
        :meth:`get_pressure` work on the world and on the individual layers, and every layer's mass (and so
        its density_bulk) is set to the mass the solved density profile places between its radii.

        Every solver setting left as ``None`` takes the ``[eos_solver]`` value of the TidalPy configuration
        (``TidalPy.config_x``), the same defaults the standalone ``radial_solver`` uses, unless this world's file
        pinned the key (see :meth:`set_solver_defaults`).

        Parameters
        ----------
        surface_pressure : float, optional
            Target surface pressure [Pa]. Default 0.0.
        slices_per_layer : int, optional
            Number of radial sample points generated per layer (>= 2).
        G_to_use : float, optional
            Gravitational constant [m^3 kg^-1 s^-2]. If negative (default), the TidalPy config value is used.
        integration_method : str, optional
            CyRK integration method: ``'DOP853'``, ``'RK45'``, ``'RK23'``, or the implicit (stiff) methods
            ``'BDF'``, ``'LSODA'``, ``'Radau'``. The structure ODE is singular at the planet's center, where
            LSODA's startup can fail to take its first step (a clean unsuccessful result) while BDF and Radau
            handle the singular start.
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
        verbose : bool, optional
            Print solver status messages. Default False.
        time : float, optional
            Time [s] the heat sources of the layers with ``use_heating`` are evaluated at, on the clock their
            radiogenics models share. ``None`` takes each model's own reference time.
        reset_layer_masses : bool, optional
            A layer that holds its mass (``is_volume_fixed = False``) forgets the mass it holds and takes the
            mass its current boundaries hold in this solve. Default False.

        Returns
        -------
        dict
            ``success``, ``message``, ``iterations``, ``max_iters_hit``, ``pressure_error`` [Pa], the radial
            profile arrays (``radius``, ``gravity``, ``pressure``, ``mass``, ``moi``, ``density``,
            ``temperature``, ``heat_flow``), the scalar results (``surface_gravity``, ``surface_pressure``,
            ``central_pressure``, ``planet_mass``, ``planet_moi``), the iteration report (``thermal_passes``,
            ``thermal_converged``, ``geometry_converged``), and the per-layer results (``layer_radius_outer``
            [m], ``layer_temperature`` [K], ``layer_heat_flow_in`` and ``layer_heat_flow_out`` [W],
            ``layer_heating`` [W], ``layer_temperature_rate`` [K s-1], ``layer_node_temperature`` and
            ``layer_top_temperature`` [K], ``layer_boundary_thickness`` [m], ``layer_rayleigh_number``,
            ``layer_nusselt_number``, and ``layer_in_thermal_network``).

        Raises
        ------
        ValueError
            If the world has no layers, any layer lacks a material EOS model, an unsupported integration
            method is given, or ``slices_per_layer < 2``.

        Assumptions
        -----------
        - Spherical symmetry; all quantities MKS.
        - Each layer's density comes from its attached material EOS model.
        """
        # The config struct starts from the [eos_solver] section of the TidalPy configuration with the keys this
        # world's file pinned on top (set_solver_defaults); only the arguments given here override it.
        cdef c_WorldEOSSolveConfig cfg = self._layered_ptr.make_eos_solve_config()
        cfg.surface_pressure = surface_pressure
        cfg.G_to_use         = G_to_use
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

        # Pure-C++ solve. Input validation throws std::invalid_argument, surfaced here as ValueError via the
        # ``except +`` on the C++ declaration. The result is copied out under the world's call lock in the same call,
        # so another thread's solve_eos cannot replace it before it is read.
        cdef c_WorldEOSReport report
        with nogil:
            report = self._layered_ptr.solve_eos_report(cfg)

        if report.max_iters_hit and not report.success:
            log_warning(
                f"World '{self.name}' EOS solve stopped at max_iters = {cfg.max_iters} with a surface-pressure "
                f"mismatch above pressure_tol = {cfg.pressure_tol:0.1e}, so the world is left unsolved. Its layers "
                f"may have no hydrostatic structure at this radius and mass; otherwise raise max_iters, or keep "
                f"pressure_tol above the integration rtol ({cfg.rtol:0.1e}).")

        return cy_eos_report_to_dict(report)

    def _build_eos_result(self):
        """The result dict of the last ``solve_eos``, from a copy of the solution taken under the world's call lock."""
        return cy_eos_report_to_dict(self._layered_ptr.get_eos_report())

    @property
    def eos_solved(self) -> bool:
        """True once the world-level EOS solve has populated the layer profiles."""
        return self._layered_ptr.get_eos_solved()

    @property
    def all_eos_set(self) -> bool:
        """True once every layer has a material EOS model attached."""
        return self._layered_ptr.get_all_eos_set()

    # Structure and viscoelastic profile queries (delegate to the containing layer).
    #
    # Every getter takes a scalar radius [m] (returning a float or complex) or a NumPy array of radii
    # (returning an array of the same shape). NaN where the EOS is unsolved, the layer is geometry-only, or
    # no rheology is attached. Each makes one C++ call for the whole input, which holds the world's call lock
    # throughout, so a read takes turns with solve_eos and the other locked calls on other threads and every value
    # of one call comes from one solve (see cy_eos_field in layers/base.pyx).
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
                    self._layered_ptr.calc_complex_moduli(
                        is_shear,
                        &flat_in[0],
                        num_radii,
                        frequency,
                        <cpp_complex[double]*><void*>&flat_out[0])
            return out_arr
        if is_shear:
            value = self._layered_ptr.calc_complex_shear_modulus(<double>radius, frequency)
        else:
            value = self._layered_ptr.calc_complex_bulk_modulus(<double>radius, frequency)
        return complex(value.real(), value.imag())

    def get_density(self, radius):
        """Density [kg/m^3] at radius [m] (float or np.ndarray); NaN if unsolved."""
        return cy_eos_field(<const void*>self._layered_ptr, cy_world_eos_fields, radius, C_EOS_DENSITY_INDEX)

    def get_gravity(self, radius):
        """Gravitational acceleration [m/s^2] at radius [m] (float or np.ndarray)."""
        return cy_eos_field(<const void*>self._layered_ptr, cy_world_eos_fields, radius, C_EOS_GRAVITY_INDEX)

    def get_pressure(self, radius):
        """Pressure [Pa] at radius [m] (float or np.ndarray); NaN if unsolved."""
        return cy_eos_field(<const void*>self._layered_ptr, cy_world_eos_fields, radius, C_EOS_PRESSURE_INDEX)

    def get_temperature(self, radius):
        """Temperature [K] at radius [m] (float or np.ndarray) from the solved profile.

        A solve with no temperature contrast reports each layer's own temperature. NaN if unsolved.
        """
        return cy_eos_field(<const void*>self._layered_ptr, cy_world_eos_fields, radius, C_EOS_TEMPERATURE_INDEX)

    def get_heat_flow(self, radius):
        """Heat flowing outward through the sphere of radius [m] (float or np.ndarray) \[W\].

        Zero everywhere when the solve carried no temperature. The flow steps across the interior of a
        convecting layer: that difference is the heat the layer stores or releases.
        """
        return cy_eos_field(<const void*>self._layered_ptr, cy_world_eos_fields, radius, C_EOS_HEAT_FLOW_INDEX)

    def get_shear_modulus(self, radius):
        """Post-melt static shear modulus [Pa] at radius [m] (float or np.ndarray)."""
        return cy_eos_field(
            <const void*>self._layered_ptr, cy_world_eos_fields, radius, C_EOS_SHEAR_MODULUS_INDEX)

    def get_bulk_modulus(self, radius):
        """Post-melt static bulk modulus [Pa] at radius [m] (float or np.ndarray)."""
        return cy_eos_field(
            <const void*>self._layered_ptr, cy_world_eos_fields, radius, C_EOS_BULK_MODULUS_INDEX)

    def get_shear_viscosity(self, radius):
        """Post-melt shear viscosity [Pa s] at radius [m] (float or np.ndarray)."""
        return cy_eos_field(
            <const void*>self._layered_ptr, cy_world_eos_fields, radius, C_EOS_SHEAR_VISCOSITY_INDEX)

    def get_bulk_viscosity(self, radius):
        """Post-melt bulk viscosity [Pa s] at radius [m] (float or np.ndarray)."""
        return cy_eos_field(
            <const void*>self._layered_ptr, cy_world_eos_fields, radius, C_EOS_BULK_VISCOSITY_INDEX)

    def get_melt_fraction(self, radius):
        """Melt fraction at radius [m] (float or np.ndarray); 0.0 where the material has no partial-melt model."""
        return cy_eos_field(
            <const void*>self._layered_ptr, cy_world_eos_fields, radius, C_EOS_MELT_FRACTION_INDEX)

    def calc_complex_shear_modulus(self, radius, double frequency, cpp_bool recalc_eos=False):
        """Complex shear modulus [Pa] at radius [m] (float or np.ndarray) and frequency [rad/s].

        Applies the containing layer's shear rheology to the stored post-melt static
        modulus + viscosity. Solves the EOS first if it has not been solved (or if
        ``recalc_eos``). NaN+0j for a geometry-only layer or no rheology.
        """
        self._ensure_solved(recalc_eos)
        return self._apply_complex(radius, frequency, True)

    def calc_complex_bulk_modulus(self, radius, double frequency, cpp_bool recalc_eos=False):
        """Complex bulk modulus [Pa] at radius [m] (float or np.ndarray) and frequency [rad/s].

        Applies the containing layer's bulk rheology to the stored post-melt static
        modulus + viscosity. Solves the EOS first if it has not been solved (or if
        ``recalc_eos``). NaN+0j for a geometry-only layer or no rheology.
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
            <const void*>self._layered_ptr, cy_world_eos_fields, radius,
            (C_EOS_SHEAR_MODULUS_INDEX, C_EOS_SHEAR_VISCOSITY_INDEX,
             C_EOS_BULK_MODULUS_INDEX, C_EOS_BULK_VISCOSITY_INDEX))

    def get_state(self, radius):
        """All EOS-related profiles at radius as a dict (float or np.ndarray values), from one evaluation of the
        solved state per radius."""
        values = cy_eos_fields(
            <const void*>self._layered_ptr, cy_world_eos_fields, radius,
            (C_EOS_DENSITY_INDEX, C_EOS_GRAVITY_INDEX, C_EOS_PRESSURE_INDEX, C_EOS_SHEAR_MODULUS_INDEX,
             C_EOS_SHEAR_VISCOSITY_INDEX, C_EOS_BULK_MODULUS_INDEX, C_EOS_BULK_VISCOSITY_INDEX,
             C_EOS_MELT_FRACTION_INDEX))
        return dict(zip(("density", "gravity", "pressure", "shear_modulus", "shear_viscosity", "bulk_modulus",
                         "bulk_viscosity", "melt_fraction"), values))

    # calc_* variants: solve the EOS first if it is unsolved (or force_recalc), then read the profile.
    def _ensure_solved(self, cpp_bool force_recalc):
        if force_recalc or not self._layered_ptr.get_eos_solved():
            self.solve_eos()

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

    @property
    def molten_regions(self) -> list:
        """Molten stretches of solid layers from the last EOS solve, as ``(layer_name, radius_inner, radius_outer)``.

        A stretch is molten where the layer's partial-melt model has weakened it past use as a solid: the post-melt
        shear modulus sits at the model's ``liquid_shear`` floor, or its rigidity mu / (rho g R) (planet bulk density,
        surface gravity, and radius) is below the ``[numerical]`` ``minimum_solid_rigidity`` of the TidalPy
        configuration. The radial solver splits the layer at the stretch's edges and solves the stretch as a static
        liquid. Radii in [m]; empty before an EOS solve or when nothing is molten.
        """
        cdef vector[c_RadialSegment] segments = self._layered_ptr.get_molten_regions()
        cdef size_t segment_i
        regions = []
        for segment_i in range(segments.size()):
            regions.append((
                self._layered_ptr.get_layer(segments[segment_i].world_layer).get_name().decode("utf-8"),
                segments[segment_i].radius_inner,
                segments[segment_i].radius_outer))
        return regions

    # Spin dynamics (the Spin model attached to the world; uses the world's EOS moment of inertia)
    def set_spin_model(self, Spin spin not None):
        """Attach a :class:`~TidalPy.dynamics_x.Spin` model (its moment-of-inertia factor is the fallback
        when the EOS has not been solved)."""
        self._layered_ptr.set_spin_model(spin._spin)

    def get_moment_of_inertia(self) -> float:
        """Moment of inertia [kg m^2]: the EOS-solved value when the EOS has been solved, else the spin
        model's ``moment_of_inertia_factor * M R^2`` estimate from the world mass and radius."""
        return self._layered_ptr.get_moment_of_inertia()

    def calc_spin_derivative(self, double host_mass) -> float:
        """Tidal spin-rate change [rad s-2] ``= M_host * dU/dO / I`` using the world's stored ``dU/dO``
        (from the last :meth:`calc_tides`) and its moment of inertia. Requires a completed tidal solve."""
        return self._layered_ptr.calc_spin_derivative(host_mass)

    def calc_synchronous_spin(self, double orbital_frequency) -> float:
        """Synchronous spin rate [rad s-1]: equal to the orbital mean motion ``orbital_frequency``."""
        return self._layered_ptr.calc_synchronous_spin(orbital_frequency)

    def solve_love_numbers(
            self,
            double frequency  = 1.0e-5,
            int degree_l      = 2,
            solve_for         = 'tidal',
            int core_model    = 0,
            use_kamata        = None,
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
            fixed_dt           = None) -> dict:
        """Solve for whole-planet tidal Love numbers (radial solver, propagation matrix, or analytic methods).

        Requires :meth:`solve_eos` first. Each layer's attached rheology is evaluated at ``frequency`` for the
        complex moduli, then the deformation ODEs are shot from the center to the surface on the solved
        structure to give k, h, l.

        Every solver setting left as ``None`` takes the ``[radial_solver]`` value of the TidalPy configuration
        (``TidalPy.config_x``), the same defaults the standalone ``radial_solver`` and the world's tidal
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
        use_kamata : bool, optional
            Use Kamata starting conditions near the center (shooting method only) instead of Takeuchi and
            Saito.
        nondimensionalize : bool, optional
            Non-dimensionalize the problem internally (recommended).
        starting_radius : float, optional
            Minimum radius [m] for the shooting start; 0 selects it automatically. Default 0.
        start_radius_tol : float, optional
            Tolerance of the automatic starting radius, R * tol^(1/l).
        integration_method : str, optional
            CyRK ODE method: ``'DOP853'``, ``'RK45'``, ``'RK23'``, or the implicit (stiff) methods ``'BDF'``,
            ``'LSODA'``, ``'Radau'``.
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
            Emit solver warnings. Default True.
        love_method : str, optional
            How the Love numbers are obtained. ``None`` (default) uses the world's configured method
            (``set_tide_config(love_method=...)`` or the ``[tides]`` table). ``'radial_solver'``
            (``'shooting'``, ``'rs'``) integrates the radial ODEs from the center to the surface.
            ``'propagation_matrix'`` (``'prop_matrix'``, ``'pm'``, ``'prop'``) is valid only for a single
            solid, static, incompressible layer; an incompatible world fails the solve gracefully, with
            ``love_success`` False and a non-zero ``love_error_code``. ``'homogeneous'`` (``'homogen'``)
            applies the homogeneous-sphere formulas with the volume-averaged complex shear modulus of the
            tidal layers (``is_tidal``), the planet's bulk density, surface gravity, and radius. ``'cpl'``
            and ``'ctl'`` apply those to the static shear modulus with a constant phase lag ``(1 - i/Q)`` or
            time lag ``(1 - i*omega*dt)``. ``'laterally_inhomogeneous'`` (``'3d'``, ``'lat_inhom'``) raises
            ``NotImplementedError``.
        fixed_q, fixed_dt : float, optional
            Quality factor for ``'cpl'`` and time lag [s] for ``'ctl'``. Left unset, the ``[tides]`` config
            values (``set_tide_config(love_fixed_q=..., love_fixed_dt=...)``, the ``love_fixed_q`` and
            ``love_fixed_dt_s`` keys of a world file) are used, then the attached tide
            model's fixed Q or time lag for this degree (``ValueError`` when none is available).

        Returns
        -------
        dict
            ``success`` (bool), ``error_code``, ``message``, ``love_method`` (canonical name),
            ``love_number_k``, ``love_number_h``, ``love_number_l`` (complex).

        Raises
        ------
        ValueError
            If the EOS has not been solved or the frequency is too small.

        Assumptions
        -----------
        - Spherical symmetry; all quantities MKS.
        """
        # The config starts from the world's [tides] Love settings (method, fixed Q, fixed time lag) and the
        # [radial_solver] section of the TidalPy configuration with the keys this world's file pinned on top
        # (set_solver_defaults); only the arguments given here override it.
        cdef c_LoveSolveConfig cfg = self._layered_ptr.make_love_solve_config()
        cfg.frequency = frequency
        cfg.degree_l  = degree_l
        cy_set_solve_for(&cfg, solve_for)
        if love_method is not None:
            cfg.love_method = cy_resolve_love_method(love_method)
        if fixed_q is not None:
            cfg.fixed_q = <double>fixed_q
        if fixed_dt is not None:
            cfg.fixed_dt = <double>fixed_dt
        cfg.core_model      = core_model
        cfg.starting_radius = starting_radius
        cfg.max_step        = max_step
        cfg.verbose         = <cpp_bool>verbose
        cfg.warnings        = <cpp_bool>warnings
        cy_apply_love_solve_overrides(
            &cfg, use_kamata, nondimensionalize, start_radius_tol, integration_method, rtol, atol, scale_rtols,
            max_num_steps, expected_size, max_ram_MB)

        with nogil:
            self._layered_ptr.solve_love_numbers(cfg)

        # The conditioning diagnostic belongs to the radial solvers.
        if warnings and c_love_method_uses_radial_solver_int(cfg.love_method):
            cy_check_surface_solve_conditioning(
                self._layered_ptr.get_love_surface_amplification(), cfg.rtol,
                self._layered_ptr.get_love_surface_rcond())
        return self._build_love_result()

    def solve_love_numbers_supplied(
            self,
            double complex[::1] complex_shear_modulus not None,
            double complex[::1] complex_bulk_modulus not None,
            double[::1] radius_array not None,
            double frequency  = 1.0e-5,
            int    degree_l   = 2,
            solve_for         = 'tidal',
            int    core_model = 0,
            use_kamata        = None,
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
            str love_method    = 'radial_solver') -> dict:
        """Solve Love numbers from externally-supplied complex moduli arrays (instead of layer rheology).

        The supplied shear/bulk moduli [Pa] are defined at ``radius_array`` [m] and are linearly interpolated onto
        the world's internal EOS radius grid. Used by the standalone ``RadialSolver_x.radial_solver`` API.
        ``solve_eos`` must be called first. Only the radial-solver methods (``love_method``
        ``'radial_solver'`` or ``'propagation_matrix'``) are available here. Solver settings left as ``None``
        take the ``[radial_solver]`` values of the TidalPy configuration, as in :meth:`solve_love_numbers`.
        """
        if radius_array.shape[0] == 0:
            raise ValueError("radius_array must not be empty")
        if (complex_shear_modulus.shape[0] != radius_array.shape[0]
                or complex_bulk_modulus.shape[0] != radius_array.shape[0]):
            raise ValueError("complex moduli and radius arrays must have matching length")

        cdef c_LoveSolveConfig cfg
        cfg.frequency   = frequency
        cfg.degree_l    = degree_l
        cy_set_solve_for(&cfg, solve_for)
        cfg.love_method = cy_resolve_love_method(love_method)
        cfg.core_model  = core_model
        cfg.starting_radius = starting_radius
        cfg.max_step        = max_step
        cfg.verbose         = <cpp_bool>verbose
        cfg.warnings        = <cpp_bool>warnings
        cy_apply_love_solve_overrides(
            &cfg, use_kamata, nondimensionalize, start_radius_tol, integration_method, rtol, atol, scale_rtols,
            max_num_steps, expected_size, max_ram_MB)

        cdef size_t n_in = radius_array.shape[0]
        cdef cpp_complex[double]* shear_ptr = <cpp_complex[double]*><void*>&complex_shear_modulus[0]
        cdef cpp_complex[double]* bulk_ptr  = <cpp_complex[double]*><void*>&complex_bulk_modulus[0]
        cdef double* radius_ptr = &radius_array[0]
        with nogil:
            self._layered_ptr.solve_love_numbers_supplied(
                cfg,
                shear_ptr,
                bulk_ptr,
                radius_ptr,
                n_in)
        if warnings:
            cy_check_surface_solve_conditioning(
                self._layered_ptr.get_love_surface_amplification(), cfg.rtol,
                self._layered_ptr.get_love_surface_rcond())
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
        cdef unique_ptr[c_RadialSolutionStorage] storage_uptr = self._layered_ptr.release_radial_storage()
        if not storage_uptr:
            raise RuntimeError(
                "No radial solution to release: run solve_love_numbers with the 'radial_solver' or "
                "'propagation_matrix' method first.")
        return RadialSolverSolution._adopt(move(storage_uptr), self)

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
            'love_method':   self.love_method,
            'love_number_k': complex(k.real(), k.imag()),
            'love_number_h': complex(h.real(), h.imag()),
            'love_number_l': complex(l.real(), l.imag()),
        }


    @property
    def love_solved(self) -> bool:
        """True while a successful ``solve_love_numbers`` result is held.

        ``solve_eos`` clears it, because Love numbers describe the structure they were solved with; the Love
        number getters return NaN until the next solve.
        """
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
    def love_surface_amplification(self) -> float:
        """Worst-case error amplification of the last solve's surface boundary condition collapse.

        Values near 1 indicate a well conditioned solve; large values mean roundoff and integration error are
        amplified into the Love numbers (the achievable relative accuracy is floor limited to about this value
        times machine epsilon). 0 until a shooting-method solve has run with ``warnings`` enabled; the propagation
        matrix method does not use the shooting surface collapse.
        """
        return self._layered_ptr.get_love_surface_amplification()

    @property
    def love_surface_rcond(self) -> float:
        """Reciprocal condition number of the last solve's surface boundary condition system.

        Units- and normalization-independent: near 1 is well posed, and a value near machine epsilon means the
        solution constants are undetermined (the solve then fails with error code -13 once it is below
        ``[numerical] minimum_surface_rcond``). A value below the integration rtol draws a conditioning warning.
        NaN before a shooting-method solve and for the other Love methods. The standalone solver reports the same
        number as ``RadialSolverSolution.surface_solve_rcond``.
        """
        return self._layered_ptr.get_love_surface_rcond()

    @property
    def love_method(self) -> str:
        """Canonical name of the Love-number method used by the last solve (``'radial_solver'`` before any)."""
        return c_love_method_name_int(self._layered_ptr.get_love_method_last_int()).decode('utf-8')

    @property
    def love_effective_shear_modulus(self) -> complex:
        """Volume-averaged shear modulus [Pa] of the tidal layers used by the last analytic Love solve.

        NaN after a radial-solver solve. Complex for the ``homogeneous`` method (evaluated at the forcing
        frequency); real (static) for ``cpl`` and ``ctl``.
        """
        cdef cpp_complex[double] v = self._layered_ptr.get_love_analytic_shear()
        return complex(v.real(), v.imag())

    @property
    def love_tidal_volume(self) -> float:
        """Volume [m3] of the tidal layers that took part in the last quasi-homogeneous Love solve (NaN otherwise)."""
        return self._layered_ptr.get_love_analytic_tidal_volume()

    @property
    def love_layer_parts(self) -> list:
        """Each tidal layer's part of the last quasi-homogeneous Love solve (``homogeneous``, ``cpl``, ``ctl``).

        One dict per layer that took part: ``layer`` (its name), ``tidal_scale``, the Love numbers of a homogeneous
        planet made of the layer's averaged material (``love_number_k``, ``love_number_h``, ``love_number_l``), and
        its complex ``shear_modulus`` [Pa] at the solve's frequency. The world's Love numbers are the sum of
        ``tidal_scale`` times these. Empty after a radial-solver solve.
        """
        cdef list parts = []
        cdef const vector[c_LayerLove]* layer_parts = &self._layered_ptr.get_love_layer_parts()
        cdef size_t i
        for i in range(layer_parts.size()):
            parts.append({
                "layer":         self._layered_ptr.get_layer(layer_parts[0][i].layer_index).get_name().decode("utf-8"),
                "tidal_scale":   layer_parts[0][i].tidal_scale,
                "love_number_k": complex(layer_parts[0][i].love.k.real(), layer_parts[0][i].love.k.imag()),
                "love_number_h": complex(layer_parts[0][i].love.h.real(), layer_parts[0][i].love.h.imag()),
                "love_number_l": complex(layer_parts[0][i].love.l.real(), layer_parts[0][i].love.l.imag()),
                "shear_modulus": complex(
                    layer_parts[0][i].shear_modulus.real(), layer_parts[0][i].shear_modulus.imag()),
            })
        return parts

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

    def get_love_radial_y(self, double radius, ytype_idx: int = 0, y_idx: int = 0) -> complex:
        """Radial function y[y_idx + 1] (SI) at ``radius`` from the last radial-solver Love solve.

        The shooting method evaluates its dense per-layer interpolants at the radius; the propagation
        matrix interpolates its grid. NaN if unsolved, after an analytic (homogeneous/cpl/ctl) solve, out
        of range, or below the solver's starting radius. ``y_idx`` 0..5 selects y1..y6.
        """
        cdef cpp_complex[double] v = self._layered_ptr.get_radial_solution_y(
            radius, <size_t>ytype_idx, <size_t>y_idx)
        return complex(v.real(), v.imag())

    def get_love_surface_y(self, ytype_idx: int, y_idx: int) -> complex:
        """Complex radial y-solution value at the surface for the given ytype and y index."""
        cdef cpp_complex[double] v = self._layered_ptr.get_love_surface_y(
            <size_t>ytype_idx, <size_t>y_idx)
        return complex(v.real(), v.imag())

    # Global (1D) tidal dissipation
    #
    # The tide model holder, config, and analytic result accessors are inherited from BaseWorld. A layered
    # world overrides calc_tides to add the rheology (radial-solver) path and the per-layer heating.
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
        with nogil:
            self._layered_ptr.calc_tides(state)

    def get_layer_tidal_heating(self, index: int) -> float:
        """Tidal heating [W] the last ``calc_tides`` put in layer ``index``; NaN before one.

        With a radial-solver Love method it is the volume integral of the radial solution's heating density over the
        layer (NaN when the tides config's ``layer_tidal_heating`` is off); with the ``homogeneous``, ``cpl``, or
        ``ctl`` method, the heating of the layer's own scaled Love numbers; with an analytic tide model, the total
        times the layer's tidal scale.
        """
        return self._layered_ptr.get_layer_tidal_heating(<size_t>index)

    def get_layer_tidal_scale(self, index: int) -> float:
        """The tidal scale layer ``index`` carries in the quasi-homogeneous Love methods [dimensionless].

        The layer's configured ``tidal_scale``, or its volume over the planet's when none is set; 0 for a layer that
        is not tidal.
        """
        return self._layered_ptr.get_layer_tidal_scale(<size_t>index)

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
        """
        cdef c_TideSolveConfig state
        state.orbital_frequency = orbital_frequency
        state.spin_frequency    = spin_frequency
        state.eccentricity      = eccentricity
        state.obliquity         = obliquity
        state.semi_major_axis   = semi_major_axis
        state.host_mass         = host_mass
        cdef double heating
        with nogil:
            heating = self._layered_ptr.get_3d_tidal_heating(state, radius, colatitude)
        return heating

    def get_3d_tidal_heating_array(
            self,
            double orbital_frequency,
            double spin_frequency,
            double eccentricity,
            double obliquity,
            double semi_major_axis,
            double host_mass,
            radii,
            colatitudes,
            int num_threads=1):
        """Longitude-mean secular 3D tidal volumetric heating [W m-3] at ``(radius, colatitude)`` points.

        Batch form of :meth:`get_3d_tidal_heating`: ``radii`` and ``colatitudes`` are paired, equal-length 1D
        arrays (point ``i`` is ``(radii[i], colatitudes[i])``) and a same-shape ``np.ndarray`` of heating is
        returned. Same physics and preconditions as the scalar method, but the world radial response is solved
        once per unique ``(degree l, |omega|)`` and reused across all points, so this is the efficient way to
        build a zonal-mean heating map.

        ``num_threads`` (default 1) spreads the per-point evaluation, which follows the radial solves on the
        calling thread, over that many threads; the result is identical for any thread count. Keep the
        default inside a process pool.
        """
        cy_check_num_threads(num_threads)
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
                state,
                &radii_view[0],
                &colat_view[0],
                num_points,
                &out_view[0],
                num_threads)
        return out_arr

    def calc_3d_displacements(
            self,
            double orbital_frequency,
            double spin_frequency,
            double eccentricity,
            double obliquity,
            double semi_major_axis,
            double host_mass,
            radii,
            colatitudes,
            longitudes,
            times,
            int num_threads=1) -> dict:
        """Instantaneous tidal displacements [m] on the grid ``(radius, colatitude, longitude, time)``.

        The active tidal modes are built from the world's ``[tides]`` truncation config, the world radial
        response is solved once per unique ``(l, frequency)``, and each mode's complex displacement
        amplitude ``(y1 U, y3 dU/dtheta, y3 dU/dphi / sin theta)`` at a point is evolved in time as
        ``Re[u e^{i omega t}]`` and summed over the modes (the same phasor convention as the instantaneous
        heating of :meth:`calc_3d_tides`). Requires the rheology tide model and a solved EOS, and a
        radial-solver Love-number method (the analytic methods have no radial functions).

        Parameters
        ----------
        orbital_frequency, spin_frequency, eccentricity, obliquity, semi_major_axis, host_mass : float
            The orbital/spin state, as for :meth:`calc_tides`.
        radii, colatitudes, longitudes, times : array-like of float
            Grid axes [m], [rad], [rad], [s]; scalars are accepted.
        num_threads : int, optional
            Threads for the per-point evaluation, which follows the radial solves on the calling thread.
            Default 1, which leaves parallelism to the caller, such as a process pool. The result is
            identical for any thread count.

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
        cdef cnp.ndarray radii_arr = np.ascontiguousarray(np.atleast_1d(radii), dtype=np.float64).ravel()
        cdef cnp.ndarray colat_arr = np.ascontiguousarray(np.atleast_1d(colatitudes), dtype=np.float64).ravel()
        cdef cnp.ndarray lon_arr   = np.ascontiguousarray(np.atleast_1d(longitudes), dtype=np.float64).ravel()
        cdef cnp.ndarray time_arr  = np.ascontiguousarray(np.atleast_1d(times), dtype=np.float64).ravel()
        cdef size_t nr = radii_arr.shape[0]
        cdef size_t nth = colat_arr.shape[0]
        cdef size_t nph = lon_arr.shape[0]
        cdef size_t nt = time_arr.shape[0]
        if nr == 0 or nth == 0 or nph == 0 or nt == 0:
            raise ValueError("radii, colatitudes, longitudes, and times must each hold at least one value")
        cdef cnp.ndarray out_arr = np.empty((nr, nth, nph, nt, 3), dtype=np.float64)
        cdef double[::1] radii_view = radii_arr
        cdef double[::1] colat_view = colat_arr
        cdef double[::1] lon_view   = lon_arr
        cdef double[::1] time_view  = time_arr
        cdef double[:, :, :, :, ::1] out_view = out_arr
        cdef c_Grid3DAxes axes
        axes.radii           = &radii_view[0]
        axes.num_radii       = nr
        axes.colatitudes     = &colat_view[0]
        axes.num_colatitudes = nth
        axes.longitudes      = &lon_view[0]
        axes.num_longitudes  = nph
        axes.times           = &time_view[0]
        axes.num_times       = nt
        cdef c_TideSolveConfig state
        state.orbital_frequency = orbital_frequency
        state.spin_frequency    = spin_frequency
        state.eccentricity      = eccentricity
        state.obliquity         = obliquity
        state.semi_major_axis   = semi_major_axis
        state.host_mass         = host_mass
        with nogil:
            self._layered_ptr.get_3d_displacements_grid(
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
            double orbital_frequency,
            double spin_frequency,
            double eccentricity,
            double obliquity,
            double semi_major_axis,
            double host_mass,
            radii,
            colatitudes,
            longitudes,
            times,
            cpp_bool return_stress=True,
            cpp_bool return_strain=True,
            int num_threads=1) -> dict:
        """Instantaneous tidal stress [Pa] and strain on the grid ``(radius, colatitude, longitude, time)``.

        The active tidal modes come from the world's ``[tides]`` truncation config and are merged into
        coherent waves; the radial response is solved once per unique ``(l, frequency)``, and at each point
        every wave's complex stress and strain amplitude is added into the total of its frequency. Each
        component at time ``t`` is the sum over frequencies of ``Re[amplitude e^{i |frequency| t}]``, the
        convention of :meth:`calc_3d_displacements` and of the instantaneous heating of
        :meth:`calc_3d_tides`.

        Parameters
        ----------
        orbital_frequency, spin_frequency, eccentricity, obliquity, semi_major_axis, host_mass : float
            The orbital and spin state, as for :meth:`calc_tides`.
        radii, colatitudes, longitudes, times : array-like of float
            Grid axes [m], [rad], [rad], [s]; scalars are accepted.
        return_stress, return_strain : bool, optional
            Which tensors to compute; each takes 48 bytes per grid point and time. Default both.
        num_threads : int, optional
            Threads for the per-point evaluation, which follows the radial solves on the calling thread.
            Default 1, which leaves parallelism to the caller, such as a process pool. The result is
            identical for any thread count.

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
            If an axis is empty, neither tensor is requested, or ``num_threads`` is below 1.
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
        cdef cnp.ndarray radii_arr = np.ascontiguousarray(np.atleast_1d(radii), dtype=np.float64).ravel()
        cdef cnp.ndarray colat_arr = np.ascontiguousarray(np.atleast_1d(colatitudes), dtype=np.float64).ravel()
        cdef cnp.ndarray lon_arr   = np.ascontiguousarray(np.atleast_1d(longitudes), dtype=np.float64).ravel()
        cdef cnp.ndarray time_arr  = np.ascontiguousarray(np.atleast_1d(times), dtype=np.float64).ravel()
        cdef size_t nr = radii_arr.shape[0]
        cdef size_t nth = colat_arr.shape[0]
        cdef size_t nph = lon_arr.shape[0]
        cdef size_t nt = time_arr.shape[0]
        if nr == 0 or nth == 0 or nph == 0 or nt == 0:
            raise ValueError("radii, colatitudes, longitudes, and times must each hold at least one value")
        cdef double[::1] radii_view = radii_arr
        cdef double[::1] colat_view = colat_arr
        cdef double[::1] lon_view   = lon_arr
        cdef double[::1] time_view  = time_arr
        cdef c_Grid3DAxes axes
        axes.radii           = &radii_view[0]
        axes.num_radii       = nr
        axes.colatitudes     = &colat_view[0]
        axes.num_colatitudes = nth
        axes.longitudes      = &lon_view[0]
        axes.num_longitudes  = nph
        axes.times           = &time_view[0]
        axes.num_times       = nt

        # Output buffers are owned by numpy; a skipped tensor passes a null pointer.
        cdef cnp.ndarray stress_arr = None
        cdef cnp.ndarray strain_arr = None
        cdef double[:, :, :, :, ::1] stress_view
        cdef double[:, :, :, :, ::1] strain_view
        cdef double* stress_ptr = NULL
        cdef double* strain_ptr = NULL
        if return_stress:
            stress_arr = np.empty((nr, nth, nph, nt, 6), dtype=np.float64)
            stress_view = stress_arr
            stress_ptr = &stress_view[0, 0, 0, 0, 0]
        if return_strain:
            strain_arr = np.empty((nr, nth, nph, nt, 6), dtype=np.float64)
            strain_view = strain_arr
            strain_ptr = &strain_view[0, 0, 0, 0, 0]

        cdef c_TideSolveConfig state
        state.orbital_frequency = orbital_frequency
        state.spin_frequency    = spin_frequency
        state.eccentricity      = eccentricity
        state.obliquity         = obliquity
        state.semi_major_axis   = semi_major_axis
        state.host_mass         = host_mass
        with nogil:
            self._layered_ptr.get_3d_stress_strain_grid(
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
            latitude_nodes=None,
            longitude_nodes=None,
            radial_slices=None,
            latitude_analytic=True,
            double colatitude_min=0.0,
            double colatitude_max=np.pi,
            int num_threads=1) -> dict:
        """3D tidal heating as a grid over ``(radius, colatitude, longitude[, time])``, optionally reduced.

        With ``orbit_averaged=True`` (default) the quantity is the secular volumetric heating density
        ``h_bar`` [W m-3], the time average of the instantaneous power at each point. It has no time axis and
        depends on longitude wherever waves at one frequency have different longitude structure, as they do
        for a synchronously rotating body. With ``orbit_averaged=False`` it is the instantaneous mechanical
        power density ``sigma_ij(t) * eps_dot_ij(t)`` [W m-3] at each supplied time (a fourth axis), which
        time-averages to ``h_bar``.

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
        instantaneous). Requires the rheology tide model and a solved EOS.

        With ``latitude_summed``, ``longitude_summed``, and ``orbit_averaged`` the colatitude integral uses the
        precomputed analytic angular Gram table of the longitude mean (exact, no theta grid);
        ``latitude_analytic=False`` falls back to the Gauss-Legendre quadrature, which agrees to quadrature accuracy,
        and a call that keeps its longitudes always uses the quadrature. A point in a liquid (or a molten stretch) has
        zero heating. A latitude band can be integrated instead of the full sphere by setting ``colatitude_min`` and
        ``colatitude_max`` [rad] (defaults 0 and pi), so complementary bands add up to the full-sphere result; a band
        narrower than the full sphere always uses the quadrature. The band has no effect when colatitude is not summed.

        ``num_threads`` (default 1) spreads the per-point evaluation, which follows the radial solves on the
        calling thread, over colatitude rows; the result is identical for any thread count. The analytic
        colatitude collapse has no per-point grid and always runs on one thread. Keep the default inside a
        process pool.
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

        # The layout gives the output shape, so the heating can be written straight into numpy-owned buffers.
        cdef c_Heating3DCollapsed layout
        with nogil:
            layout = self._layered_ptr.calc_3d_tides_layout(
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
            self._layered_ptr.calc_3d_tides_into(
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
            Keys of the matching section of ``TidalPy_Configs_x.toml`` (``EOS_SOLVER_KEYS`` and
            ``RADIAL_SOLVER_KEYS`` of ``TidalPy.structures_x.configs``) with their values. A table given replaces
            the one stored, so an empty dict clears it; a table left as ``None`` is untouched.

        Raises
        ------
        ValueError
            For a key the section does not have, a value of the wrong type or out of range, or an unknown
            integration method.
        """
        from TidalPy.structures_x.configs.toml_loader import validate_solver_table
        cdef c_EOSSolverOverrides eos
        cdef c_RadialSolverOverrides radial
        cdef ODEMethod method
        cdef double number
        cdef size_t count
        cdef cpp_bool flag
        if eos_solver is not None:
            validate_solver_table("eos_solver", eos_solver, "set_solver_defaults")
            if "integration_method" in eos_solver:
                method = cy_resolve_integration_method(str(eos_solver["integration_method"]))
                eos.integration_method = optional[ODEMethod](method)
            if "rtol" in eos_solver:
                number = <double>eos_solver["rtol"]
                eos.rtol = optional[double](number)
            if "atol" in eos_solver:
                number = <double>eos_solver["atol"]
                eos.atol = optional[double](number)
            if "pressure_tol" in eos_solver:
                number = <double>eos_solver["pressure_tol"]
                eos.pressure_tol = optional[double](number)
            if "max_iters" in eos_solver:
                count = <size_t>int(eos_solver["max_iters"])
                eos.max_iters = optional[size_t](count)
            if "slices_per_layer" in eos_solver:
                count = <size_t>int(eos_solver["slices_per_layer"])
                eos.slices_per_layer = optional[size_t](count)
            if "nondimensionalize" in eos_solver:
                flag = <cpp_bool>bool(eos_solver["nondimensionalize"])
                eos.nondimensionalize = optional[cpp_bool](flag)
            if "solve_temperature" in eos_solver:
                flag = <cpp_bool>bool(eos_solver["solve_temperature"])
                eos.solve_temperature = optional[cpp_bool](flag)
            self._layered_ptr.set_eos_solver_overrides(eos)
        if radial_solver is not None:
            validate_solver_table("radial_solver", radial_solver, "set_solver_defaults")
            if "integration_method" in radial_solver:
                method = cy_resolve_integration_method(str(radial_solver["integration_method"]))
                radial.integration_method = optional[ODEMethod](method)
            if "rtol" in radial_solver:
                number = <double>radial_solver["rtol"]
                radial.rtol = optional[double](number)
            if "atol" in radial_solver:
                number = <double>radial_solver["atol"]
                radial.atol = optional[double](number)
            if "use_kamata" in radial_solver:
                flag = <cpp_bool>bool(radial_solver["use_kamata"])
                radial.use_kamata = optional[cpp_bool](flag)
            if "start_radius_tolerance" in radial_solver:
                number = <double>radial_solver["start_radius_tolerance"]
                radial.start_radius_tol = optional[double](number)
            if "scale_rtols" in radial_solver:
                flag = <cpp_bool>bool(radial_solver["scale_rtols"])
                radial.scale_rtols = optional[cpp_bool](flag)
            if "max_num_steps" in radial_solver:
                count = <size_t>int(radial_solver["max_num_steps"])
                radial.max_num_steps = optional[size_t](count)
            if "expected_size" in radial_solver:
                count = <size_t>int(radial_solver["expected_size"])
                radial.expected_size = optional[size_t](count)
            if "max_ram_mb" in radial_solver:
                count = <size_t>int(radial_solver["max_ram_mb"])
                radial.max_ram_MB = optional[size_t](count)
            if "nondimensionalize" in radial_solver:
                flag = <cpp_bool>bool(radial_solver["nondimensionalize"])
                radial.nondimensionalize = optional[cpp_bool](flag)
            self._layered_ptr.set_radial_solver_overrides(radial)

    def get_solver_defaults(self) -> dict:
        """The ``[eos_solver]`` and ``[radial_solver]`` keys pinned on this world, under the configuration's names.

        Returns
        -------
        dict
            ``eos_solver`` and ``radial_solver`` tables, each present only when it pins a key. Empty when the
            world follows the TidalPy configuration throughout.
        """
        cdef c_EOSSolverOverrides eos = self._layered_ptr.get_eos_solver_overrides()
        cdef c_RadialSolverOverrides radial = self._layered_ptr.get_radial_solver_overrides()
        cdef dict out = {}
        cdef dict table = {}
        if eos.integration_method.has_value():
            table["integration_method"] = cy_integration_method_name(eos.integration_method.value())
        if eos.rtol.has_value():
            table["rtol"] = eos.rtol.value()
        if eos.atol.has_value():
            table["atol"] = eos.atol.value()
        if eos.pressure_tol.has_value():
            table["pressure_tol"] = eos.pressure_tol.value()
        if eos.max_iters.has_value():
            table["max_iters"] = <int>eos.max_iters.value()
        if eos.slices_per_layer.has_value():
            table["slices_per_layer"] = <int>eos.slices_per_layer.value()
        if eos.nondimensionalize.has_value():
            table["nondimensionalize"] = bool(eos.nondimensionalize.value())
        if eos.solve_temperature.has_value():
            table["solve_temperature"] = bool(eos.solve_temperature.value())
        if table:
            out["eos_solver"] = table
        table = {}
        if radial.integration_method.has_value():
            table["integration_method"] = cy_integration_method_name(radial.integration_method.value())
        if radial.rtol.has_value():
            table["rtol"] = radial.rtol.value()
        if radial.atol.has_value():
            table["atol"] = radial.atol.value()
        if radial.use_kamata.has_value():
            table["use_kamata"] = bool(radial.use_kamata.value())
        if radial.start_radius_tol.has_value():
            table["start_radius_tolerance"] = radial.start_radius_tol.value()
        if radial.scale_rtols.has_value():
            table["scale_rtols"] = bool(radial.scale_rtols.value())
        if radial.max_num_steps.has_value():
            table["max_num_steps"] = <int>radial.max_num_steps.value()
        if radial.expected_size.has_value():
            table["expected_size"] = <int>radial.expected_size.value()
        if radial.max_ram_MB.has_value():
            table["max_ram_mb"] = <int>radial.max_ram_MB.value()
        if radial.nondimensionalize.has_value():
            table["nondimensionalize"] = bool(radial.nondimensionalize.value())
        if table:
            out["radial_solver"] = table
        return out

    cpdef dict get_config_dict(self):
        """Return the world config with a ``layers`` table keyed by layer name.

        Each entry is the layer's own ``get_config_dict`` (``class``, scalars, attached-model sub-tables)
        minus the standalone-only keys the builder derives itself (``name``, ``radius_inner``, the
        Love-number components), so the result validates against the world schema and rebuilds the same
        structure through ``build_world_from_dict``.

        Returns
        -------
        dict
            All :class:`BaseWorld` keys, ``moment_of_inertia_factor`` (the attached spin model's), ``layers``, and
            the ``eos_solver`` and ``radial_solver`` tables when the world pins any solver key
            (:meth:`set_solver_defaults`).

        Raises
        ------
        ValueError
            If two layers share a name (the table needs unique keys).
        """
        cdef dict config = BaseWorld.get_config_dict(self)
        config["moment_of_inertia_factor"] = self._layered_ptr.get_spin_model().get_config().moment_of_inertia_factor
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
        config["layers"] = layers
        return config
