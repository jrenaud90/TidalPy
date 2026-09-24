# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrapper for TidalPy's base layer class.

BaseLayer is geometry only: inner and outer radii, mass, and material identity. Its EOS profile stays unpopulated
until the world's EOS solve runs or ``update_eos_data`` is called directly.
"""

import os as _os

cimport numpy as cnp
cnp.import_array()

import numpy as np

from libcpp.vector cimport vector
from libcpp cimport bool as cpp_bool
from libcpp.utility cimport move
from libcpp.memory cimport make_unique

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport d_NAN, set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport (
    StructureBase,
    c_TidalPyBaseClass,
    c_PhysicsBase,
    cy_physics_model_config,
)
from TidalPy.Material_x.eos.material_eos cimport MaterialEOSBase, cy_material_config

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())

# Layer config keys that are constructor parameters of a standalone layer but not part of the world
# builder's layer schema: the world writer drops them when nesting a layer under its name.
LAYER_STANDALONE_CONFIG_KEYS = (
    "name",
    "radius_inner_m",
    "love_number_k_re",
    "love_number_k_im",
    "love_number_h_re",
    "love_number_h_im",
    "love_number_l_re",
    "love_number_l_im",
)


# Selectors for the vectorized real-valued radius getters (see _eval_real). Mirrors the world-level
# surface in structures_x/worlds/layered.pyx so layers and worlds share one calling convention.
cdef enum:
    _KIND_DENSITY        = 0
    _KIND_GRAVITY        = 1
    _KIND_PRESSURE       = 2
    _KIND_SHEAR_MOD      = 3
    _KIND_BULK_MOD       = 4
    _KIND_SHEAR_VISC     = 5
    _KIND_BULK_VISC      = 6
    _KIND_MELT_FRACTION  = 11


cdef object cy_eos_fields(const void* owner, cy_eos_state_fn fill, object radius, tuple indices):
    """The dense-layout values at ``indices`` from one evaluation per radius, as a tuple: floats for a scalar
    radius, arrays shaped like ``radius`` for an array."""
    cdef vector[double] state = vector[double](C_EOS_DY_VALUES)
    cdef vector[size_t] field_index
    cdef object index
    for index in indices:
        field_index.push_back(<size_t>index)
    cdef size_t num_fields = field_index.size()
    cdef size_t field_i
    cdef cnp.ndarray in_arr
    cdef cnp.ndarray out_arr
    cdef double[::1] flat_in
    cdef double[:, ::1] flat_out
    cdef Py_ssize_t i, n
    if isinstance(radius, np.ndarray):
        in_arr  = np.ascontiguousarray(radius, dtype=np.float64)
        flat_in = in_arr.reshape(-1)
        n = flat_in.shape[0]
        out_arr = np.empty((num_fields, n), dtype=np.float64)
        flat_out = out_arr
        with nogil:
            for i in range(n):
                fill(owner, flat_in[i], state.data())
                for field_i in range(num_fields):
                    flat_out[field_i, i] = state[field_index[field_i]]
        shape = np.shape(in_arr)
        return tuple([out_arr[field_i].reshape(shape) for field_i in range(num_fields)])
    fill(owner, <double>radius, state.data())
    return tuple([state[field_index[field_i]] for field_i in range(num_fields)])


cdef void _layer_eos_state(const void* owner, double radius, double* y_out) noexcept nogil:
    (<const c_BaseLayer*>owner).get_eos_state(radius, y_out)


cdef class BaseLayer(StructureBase):
    """Geometry base layer: inner and outer radii, mass, and material identity.

    The EOS profile (density, gravity, pressure against radius) starts unpopulated and becomes queryable after the
    world's EOS solve or a direct ``update_eos_data`` call.

    Parameters
    ----------
    name : str
        Layer name (e.g. ``"mantle"``).
    layer_index : int
        Zero-based position in the parent world, innermost layer = 0.
    radius_inner : float
        Inner boundary radius [m].
    radius_outer : float
        Outer boundary radius [m].
    mass : float
        Total layer mass [kg].
    material_name : str, optional
        Material identifier (e.g. ``"perovskite"``). Default ``""``.
    is_volume_fixed : bool, optional
        False lets the layer grow or shrink to hold its mass while the EOS solve redistributes the interior.
        Default ``True``.
    is_tidal : bool, optional
        Whether this layer contributes to tidal dissipation. Default ``True``.
    tidal_scale : float, optional
        The layer's share of the planet in the quasi-homogeneous Love methods (``homogeneous``, ``cpl``, ``ctl``),
        which scale the Im(k) of a homogeneous planet made of the layer's averaged material by it. ``None``
        (default) takes the layer's volume over the planet's.

    Assumptions
    -----------
    - Layer geometry is spherically symmetric.
    - radius_inner <= radius_outer.
    """

    def __cinit__(self, *args, **kwargs):
        # unique_ptr<c_BaseLayer> auto-inits to nullptr; _ptr set in __init__.
        self._is_view   = False
        self._world_ref = None
        self._detached  = False

    cdef void _check_ptr(self) except *:
        if self._detached:
            raise RuntimeError(
                "This layer view no longer refers to a layer: its world loaded a binary file, which replaced the "
                "world's layers. Take a new view from the world (world.<layer name>, get_layer, or iteration).")
        StructureBase._check_ptr(self)

    cdef void _detach(self) noexcept:
        if self._is_view:
            self._layer_ptr.release()
        self._ptr      = NULL
        self._detached = True

    cdef void _notify_world_of_move(self) except *:
        if self._is_view and (self._world_ref is not None):
            self._world_ref._layer_moved()

    def __init__(
            self,
            str    name,
            int    layer_index,
            double radius_inner,
            double radius_outer,
            double mass,
            str    material_name      = "",
            cpp_bool   is_tidal           = True,
            cpp_bool   is_volume_fixed    = True,
            tidal_scale               = None):
        cdef c_BaseLayerConfig config
        config.name         = name.encode("utf-8")
        config.layer_index  = layer_index
        config.radius_inner = radius_inner
        config.radius_outer = radius_outer
        config.mass         = mass
        config.material_name = material_name.encode("utf-8")
        config.is_tidal    = is_tidal
        config.is_volume_fixed = is_volume_fixed
        config.tidal_scale = d_NAN if tidal_scale is None else <double>tidal_scale
        # The owning member is this same type, so make_unique's result moves straight in.
        self._layer_ptr = make_unique[c_BaseLayer](config)
        self._ptr = <c_TidalPyBaseClass*>self._layer_ptr.get()

    def __dealloc__(self):
        if self._is_view:
            # The C++ layer is owned by the world; relinquish without deleting it.
            self._layer_ptr.release()
        else:
            self._layer_ptr.reset()
        self._ptr = NULL

    cdef void _init_view(self, c_BaseLayer* ptr, object world):
        """Set this wrapper up as a non-owning view onto a world-owned C++ layer.

        Keeps a reference to ``world`` so the C++ layer outlives the view; ``__dealloc__`` then releases the
        pointer instead of deleting it.
        """
        self._layer_ptr.reset(ptr)
        self._ptr       = <c_TidalPyBaseClass*>ptr
        self._is_view   = True
        self._world_ref = world

    @staticmethod
    cdef BaseLayer _view(c_BaseLayer* ptr, object world):
        cdef BaseLayer v = BaseLayer.__new__(BaseLayer)
        v._init_view(ptr, world)
        return v

    def load_binary(self, str path, cpp_bool force=False):
        """Load this layer's state from a TidalPy binary file.

        Only a standalone layer can be loaded: a layer view belongs to its world, whose structure a load would
        change behind its back, so load the world instead.

        Parameters
        ----------
        path : str
            Source file path.
        force : bool, optional
            Attempt the load even on a schema version mismatch.
        """
        self._check_ptr()
        if self._is_view:
            raise ValueError(
                f"Layer '{self.name}' belongs to a world and cannot be loaded in place: load the world's binary "
                f"file, or load into a standalone layer.")
        StructureBase.load_binary(self, path, force)

    # _layer_ptr always points at the most-derived C++ layer, so these are safe for subclasses.
    @property
    def radius(self) -> float:
        """Outer radius [m]."""
        self._check_ptr()
        return self._layer_ptr.get().get_radius_outer()

    @property
    def mass(self) -> float:
        """Total layer mass [kg].

        Each successful world EOS solve overwrites it with the mass the solved density profile places between the
        layer's inner and outer radii.
        """
        self._check_ptr()
        return self._layer_ptr.get().get_mass()

    @property
    def name(self) -> str:
        """Layer name."""
        self._check_ptr()
        return self._layer_ptr.get().get_name().decode("utf-8")

    @property
    def layer_index(self) -> int:
        """Zero-based layer index (0 = innermost)."""
        self._check_ptr()
        return self._layer_ptr.get().get_layer_index()

    @property
    def radius_inner(self) -> float:
        """Inner boundary radius [m]."""
        self._check_ptr()
        return self._layer_ptr.get().get_radius_inner()

    @property
    def radius_outer(self) -> float:
        """Outer boundary radius [m]."""
        self._check_ptr()
        return self._layer_ptr.get().get_radius_outer()

    @property
    def thickness(self) -> float:
        """Layer thickness [m] (radius_outer - radius_inner)."""
        self._check_ptr()
        return self._layer_ptr.get().get_thickness()

    @property
    def volume(self) -> float:
        """Layer volume [m^3] (spherical shell)."""
        self._check_ptr()
        return self._layer_ptr.get().get_volume()

    @property
    def density_bulk(self) -> float:
        """Bulk density [kg/m^3] = mass / volume, NaN for a zero-volume layer; follows the EOS-set mass."""
        self._check_ptr()
        return self._layer_ptr.get().get_density_bulk()

    @property
    def surface_area_inner(self) -> float:
        """Inner boundary surface area [m^2]."""
        self._check_ptr()
        return self._layer_ptr.get().get_surface_area_inner()

    @property
    def surface_area_outer(self) -> float:
        """Outer boundary surface area [m^2]."""
        self._check_ptr()
        return self._layer_ptr.get().get_surface_area_outer()

    @property
    def material_name(self) -> str:
        """Material identifier string."""
        self._check_ptr()
        return self._layer_ptr.get().get_material_name().decode("utf-8")

    @property
    def is_volume_fixed(self) -> bool:
        """False if the layer grows or shrinks to hold its mass during an EOS solve."""
        self._check_ptr()
        return self._layer_ptr.get().get_is_volume_fixed()

    @is_volume_fixed.setter
    def is_volume_fixed(self, value: bool):
        self._check_ptr()
        self._layer_ptr.get().set_is_volume_fixed(<cpp_bool>bool(value))

    @property
    def is_tidal(self) -> bool:
        """Whether this layer contributes to tidal dissipation."""
        self._check_ptr()
        return self._layer_ptr.get().get_is_tidal()

    @property
    def tidal_scale(self):
        """The layer's configured tidal scale [dimensionless], or ``None`` when it takes its volume fraction.

        Used only by the quasi-homogeneous Love methods (``homogeneous``, ``cpl``, ``ctl``) and to share out the
        heating of an analytic tide model; the radial solver resolves the layers directly. Settable; ``None``
        returns to the volume fraction. The value in use is ``LayeredWorld.get_layer_tidal_scale``.
        """
        self._check_ptr()
        cdef double value = self._layer_ptr.get().get_tidal_scale()
        return None if value != value else value

    @tidal_scale.setter
    def tidal_scale(self, value):
        self._check_ptr()
        self._layer_ptr.get().set_tidal_scale(d_NAN if value is None else <double>value)

    def get_tidal_heating(self) -> float:
        """Tidal heating [W] deposited in this layer by the world's last tidal solve. NaN before one runs.

        Set by :meth:`LayeredWorld.calc_tides`; how the heating is resolved per layer depends on the world's Love
        method (see the worlds documentation, Tidal Heating of Each Layer).
        """
        self._check_ptr()
        return self._layer_ptr.get().get_tidal_heating()

    @property
    def eos_data_populated(self) -> bool:
        """True after EOS profile data has been populated (world EOS solve or update_eos_data)."""
        self._check_ptr()
        return self._layer_ptr.get().get_eos_data_populated()

    @property
    def eos_set(self) -> bool:
        """True after a material EOS model has been attached via :meth:`set_eos`."""
        self._check_ptr()
        return self._layer_ptr.get().get_eos_set()

    def set_radii(self, double radius_inner, double radius_outer):
        """Move the layer's boundaries [m], keeping every derived geometric quantity in step.

        The world EOS solve calls this itself when a layer holding its mass grows or shrinks. Setting the
        radii by hand on a world's layer leaves the world's own radius and its other layers untouched, so keep
        the stack continuous, and forgets the world's solved profile, which no longer lines up with its layers.
        """
        self._check_ptr()
        self._layer_ptr.get().set_radii(radius_inner, radius_outer)
        self._notify_world_of_move()

    def set_eos(self, MaterialEOSBase eos not None):
        """Attach a material EOS model, the layer's density source for the world-level ``solve_eos``.

        Ownership of the C++ model moves out of ``eos``, which is left an empty shell and must not be reused.

        Parameters
        ----------
        eos : MaterialEOSBase
            A material EOS model, for example ``make_material_eos("constant", ...)``.

        Raises
        ------
        ValueError
            If ``eos`` has already been attached or otherwise moved.
        """
        self._check_ptr()
        if eos._eos_ptr.get() == NULL:
            raise ValueError(
                "This EOS model holds no C++ object (already attached or moved).")
        self._layer_ptr.get().set_eos(move(eos._eos_ptr))
        eos._ptr = NULL

    def update_eos_data(
            self,
            radius,
            density_kgm3,
            gravity_ms2,
            pressure):
        """Populate the EOS profile directly from radius arrays, bypassing the world's EOS solve.

        Parameters
        ----------
        radius : sequence of float
            Radius values [m], sorted ascending, matching the layer bounds.
        density_kgm3 : sequence of float
            Mass density [kg/m^3] at each radius.
        gravity_ms2 : sequence of float
            Gravitational acceleration [m/s^2] at each radius.
        pressure : sequence of float
            Pressure [Pa] at each radius.

        Assumptions
        -----------
        - All four sequences have the same length and radius is strictly ascending.
        """
        self._check_ptr()
        cdef vector[double] r_vec   = radius
        cdef vector[double] rho_vec = density_kgm3
        cdef vector[double] g_vec = gravity_ms2
        cdef vector[double] p_vec = pressure
        cdef c_LayerEOSData eos_data
        eos_data.populate(r_vec, rho_vec, g_vec, p_vec)
        self._layer_ptr.get().update_eos_data(eos_data)

    cdef double _eval_real(self, int kind, double radius) noexcept nogil:
        cdef c_BaseLayer* layer = self._layer_ptr.get()
        if   kind == _KIND_DENSITY:        return layer.get_density(radius)
        elif kind == _KIND_GRAVITY:        return layer.get_gravity(radius)
        elif kind == _KIND_PRESSURE:       return layer.get_pressure(radius)
        elif kind == _KIND_SHEAR_MOD:      return layer.get_shear_modulus(radius)
        elif kind == _KIND_BULK_MOD:       return layer.get_bulk_modulus(radius)
        elif kind == _KIND_SHEAR_VISC:     return layer.get_shear_viscosity(radius)
        elif kind == _KIND_BULK_VISC:      return layer.get_bulk_viscosity(radius)
        elif kind == _KIND_MELT_FRACTION:  return layer.get_melt_fraction(radius)
        return 0.0

    def _apply_real(self, radius, int kind):
        self._check_ptr()
        # float -> float; np.ndarray -> np.ndarray (same shape, looped under nogil).
        cdef cnp.ndarray in_arr
        cdef cnp.ndarray out_arr
        cdef double[::1] flat_in
        cdef double[::1] flat_out
        cdef Py_ssize_t i, n
        if isinstance(radius, np.ndarray):
            in_arr  = np.ascontiguousarray(radius, dtype=np.float64)
            out_arr = np.empty_like(in_arr)
            flat_in = in_arr.reshape(-1)
            flat_out = out_arr.reshape(-1)
            n = flat_in.shape[0]
            with nogil:
                for i in range(n):
                    flat_out[i] = self._eval_real(kind, flat_in[i])
            return out_arr
        return self._eval_real(kind, <double>radius)

    def get_density(self, radius):
        """Density [kg/m^3] at radius [m] (float or np.ndarray); NaN if EOS data not populated."""
        self._check_ptr()
        return self._apply_real(radius, _KIND_DENSITY)

    def get_gravity(self, radius):
        """Gravitational acceleration [m/s^2] at radius [m] (float or np.ndarray); NaN if not populated."""
        self._check_ptr()
        return self._apply_real(radius, _KIND_GRAVITY)

    def get_pressure(self, radius):
        """Pressure [Pa] at radius [m] (float or np.ndarray); NaN if EOS data not populated."""
        self._check_ptr()
        return self._apply_real(radius, _KIND_PRESSURE)

    # Viscoelastic profile (populated by the world EOS solve; NaN before then or on a geometry-only layer)
    @property
    def viscoelastic_populated(self) -> bool:
        """True after the world EOS solve has populated this layer's viscoelastic state."""
        self._check_ptr()
        return self._layer_ptr.get().get_viscoelastic_populated()

    def get_shear_modulus(self, radius):
        """Post-melt static shear modulus [Pa] at radius [m] (float or np.ndarray); NaN if unpopulated."""
        self._check_ptr()
        return self._apply_real(radius, _KIND_SHEAR_MOD)

    def get_bulk_modulus(self, radius):
        """Post-melt static bulk modulus [Pa] at radius [m] (float or np.ndarray); NaN if unpopulated."""
        self._check_ptr()
        return self._apply_real(radius, _KIND_BULK_MOD)

    def get_shear_viscosity(self, radius):
        """Post-melt shear viscosity [Pa s] at radius [m] (float or np.ndarray); NaN if unpopulated."""
        self._check_ptr()
        return self._apply_real(radius, _KIND_SHEAR_VISC)

    def get_bulk_viscosity(self, radius):
        """Post-melt bulk viscosity [Pa s] at radius [m] (float or np.ndarray); NaN if unpopulated."""
        self._check_ptr()
        return self._apply_real(radius, _KIND_BULK_VISC)

    def get_melt_fraction(self, radius):
        """Melt fraction at radius [m] (float or np.ndarray) from the attached partial-melt model.

        0.0 where no partial-melt model is attached; NaN if unpopulated.
        """
        self._check_ptr()
        return self._apply_real(radius, _KIND_MELT_FRACTION)

    # Shorthand bundles (one call returns several profiles at once; mirrors the world-level surface)
    def get_static_viscoelastics(self, radius):
        """``(shear_modulus, shear_viscosity, bulk_modulus, bulk_viscosity)`` (post-melt) at radius [m], each a
        float or np.ndarray. One evaluation of the solved state per radius fills all four.
        """
        self._check_ptr()
        return cy_eos_fields(
            <const void*>self._layer_ptr.get(), _layer_eos_state, radius,
            (C_EOS_SHEAR_MODULUS_INDEX, C_EOS_SHEAR_VISCOSITY_INDEX,
             C_EOS_BULK_MODULUS_INDEX, C_EOS_BULK_VISCOSITY_INDEX))

    def get_state(self, radius):
        """All EOS-related profiles at radius as a dict (float or np.ndarray values), from one evaluation of the
        solved state per radius."""
        self._check_ptr()
        values = cy_eos_fields(
            <const void*>self._layer_ptr.get(), _layer_eos_state, radius,
            (C_EOS_DENSITY_INDEX, C_EOS_GRAVITY_INDEX, C_EOS_PRESSURE_INDEX, C_EOS_SHEAR_MODULUS_INDEX,
             C_EOS_SHEAR_VISCOSITY_INDEX, C_EOS_BULK_MODULUS_INDEX, C_EOS_BULK_VISCOSITY_INDEX,
             C_EOS_MELT_FRACTION_INDEX))
        return dict(zip(("density", "gravity", "pressure", "shear_modulus", "shear_viscosity", "bulk_modulus",
                         "bulk_viscosity", "melt_fraction"), values))

    cpdef dict get_config_dict(self):
        """Return all configuration values as a Python dict (MKS) in the world builder's layer schema.

        Each attached physics model appears as its own sub-table keyed by ``model``. The material ``type`` is
        written as ``"none"`` so that a rebuild does not add the material defaults a typeless layer would take.
        ``name`` and ``radius_inner`` are standalone-layer keys that a world drops when it nests the layer (see
        ``LAYER_STANDALONE_CONFIG_KEYS``).

        Returns
        -------
        dict
            Keys: ``class``, ``type``, ``name``, ``layer_index``, ``radius_inner``, ``radius_outer``, ``mass``,
            ``material_name``, ``is_tidal``, ``is_volume_fixed``, ``tidal_scale`` when one is set, and
            ``material`` when set.
        """
        # Deferred: the configs package imports the layer modules.
        from TidalPy.structures_x.configs.toml_loader import NO_MATERIAL_TYPE

        cdef c_BaseLayer* p = self._layer_ptr.get()
        cdef bytes class_bytes = c_layer_class_name(p.get_layer_class_id())
        cdef double tidal_scale = p.get_tidal_scale()
        cdef dict config = {
            "class":              class_bytes.decode("utf-8"),
            "type":               NO_MATERIAL_TYPE,
            "name":               p.get_name().decode("utf-8"),
            "layer_index":        p.get_layer_index(),
            "radius_inner_m":     p.get_radius_inner(),
            "radius_outer_m":     p.get_radius_outer(),
            "mass_kg":            p.get_mass(),
            "material_name":      p.get_material_name().decode("utf-8"),
            "is_tidal":           bool(p.get_is_tidal()),
            "is_volume_fixed":    bool(p.get_is_volume_fixed()),
        }
        if tidal_scale == tidal_scale:
            config["tidal_scale"] = tidal_scale
        if p.get_eos_set():
            config["material"] = cy_material_config(p.get_eos())
        return config
