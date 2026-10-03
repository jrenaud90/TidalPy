# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrapper for TidalPy's layer class.

A Layer is a shell of one material between two radii: its geometry, the material it is made of, the physics switches
that say how much of that material it uses, its temperature, the radial-solver assumptions, optional rheology
overrides, and optional cooling and radiogenics models. Its profile stays unpopulated until its world's EOS solve runs
or ``update_eos_data`` is called directly.
"""

cimport numpy as cnp
cnp.import_array()

import numpy as np

from libcpp.vector cimport vector
from libcpp.complex cimport complex as cpp_complex
from libcpp cimport bool as cpp_bool
from libcpp.memory cimport make_unique, shared_ptr
from cpython.complex cimport PyComplex_FromDoubles

from TidalPy.Utilities.logging.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport d_NAN, set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities.classes.classes cimport (
    PhysicsBase,
    c_TidalPyBaseClass,
    c_PhysicsBase,
    cy_physics_model_config,
    cy_wrap_model,
)
from TidalPy.Radiogenics.radiogenics import make_radiogenics
from TidalPy.Material import Material, load_material
from TidalPy.Rheology import make_rheology
from TidalPy.Cooling.cooling import make_cooling

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())

# Layer config keys that are constructor parameters of a standalone layer but not part of the world builder's layer
# schema: the world writer drops them when nesting a layer under its name.
LAYER_STANDALONE_CONFIG_KEYS = (
    "name",
    "radius_inner_m",
)

# The material switches, in c_MaterialSwitches order.
LAYER_SWITCHES = ("use_thermal_expansion", "use_melting", "use_pressure_melting", "use_melt_density")


# The profile getters of layers and worlds share these two helpers. Each makes one C++ call for the whole input, which
# holds the owning world's call lock throughout (see c_WorldCallLock), so a read takes turns with a solve_eos on
# another thread and every value of one call comes from one solve. An array is read without the GIL; a scalar keeps
# it, which saves the release on the fast path and is safe because a thread holding the call lock never waits for
# the GIL.
cdef object cy_eos_field(const void* owner, cy_eos_fields_fn fill, object radius, size_t field_index):
    """One dense-layout value at ``field_index``: a float for a scalar radius, an array shaped like ``radius`` for an
    array."""
    cdef cnp.ndarray in_arr
    cdef cnp.ndarray out_arr
    cdef double[::1] flat_in
    cdef double[::1] flat_out
    cdef size_t num_radii
    cdef size_t field = field_index
    cdef double radius_value
    cdef double value = d_NAN
    if isinstance(radius, np.ndarray):
        in_arr   = np.ascontiguousarray(radius, dtype=np.float64)
        out_arr  = np.empty_like(in_arr)
        flat_in  = in_arr.reshape(-1)
        flat_out = out_arr.reshape(-1)
        num_radii = <size_t>flat_in.shape[0]
        if num_radii > 0:
            with nogil:
                fill(owner, &field, 1, &flat_in[0], num_radii, &flat_out[0])
        return out_arr
    radius_value = <double>radius
    fill(owner, &field, 1, &radius_value, 1, &value)
    return value


cdef object cy_eos_fields(const void* owner, cy_eos_fields_fn fill, object radius, tuple indices):
    """The dense-layout values at ``indices`` from one evaluation per radius, as a tuple: floats for a scalar
    radius, arrays shaped like ``radius`` for an array."""
    cdef size_t num_fields = <size_t>len(indices)
    cdef vector[size_t] field_index = vector[size_t](num_fields)
    cdef vector[double] values = vector[double](num_fields)
    cdef size_t field_i
    cdef cnp.ndarray in_arr
    cdef cnp.ndarray out_arr
    cdef double[::1] flat_in
    cdef double[:, ::1] flat_out
    cdef size_t num_radii
    cdef double radius_value
    for field_i in range(num_fields):
        field_index[field_i] = <size_t>indices[field_i]
    if isinstance(radius, np.ndarray):
        in_arr    = np.ascontiguousarray(radius, dtype=np.float64)
        flat_in   = in_arr.reshape(-1)
        num_radii = <size_t>flat_in.shape[0]
        out_arr   = np.empty((num_fields, num_radii), dtype=np.float64)
        if num_radii > 0 and num_fields > 0:
            flat_out = out_arr
            with nogil:
                fill(owner, field_index.data(), num_fields, &flat_in[0], num_radii, &flat_out[0, 0])
        shape = np.shape(in_arr)
        return tuple([out_arr[field_i].reshape(shape) for field_i in range(num_fields)])
    radius_value = <double>radius
    fill(owner, field_index.data(), num_fields, &radius_value, 1, values.data())
    return tuple([values[field_i] for field_i in range(num_fields)])


cdef void cy_layer_eos_fields(
        const void* owner,
        const size_t* field_indices,
        size_t num_fields,
        const double* radii,
        size_t num_radii,
        double* values_out) noexcept nogil:
    (<const c_Layer*>owner).get_eos_fields(field_indices, num_fields, radii, num_radii, values_out)


def _as_material(object material):
    """A material from a Material, a MatPack name, or a material table; None stays None."""
    if material is None or isinstance(material, Material):
        return material
    if isinstance(material, (str, dict)):
        return load_material(material)
    raise TypeError(f"TidalPy: a layer's material is a Material, a MatPack name, or a material table, not "
                    f"{type(material).__name__}.")


def _as_model(object model, str slot, object make_model, str family):
    """A layer's model from a model, a model name, or a config table with a ``model`` key, built by ``make_model``;
    None stays None. ``family`` names the kind of model in the error."""
    if model is None or isinstance(model, PhysicsBase):
        return model
    if isinstance(model, str):
        return make_model(model, {})
    if isinstance(model, dict):
        table = dict(model)
        if "model" not in table:
            raise ValueError(f"TidalPy: the layer's '{slot}' table needs a 'model' key naming its model.")
        return make_model(str(table.pop("model")), table)
    raise TypeError(f"TidalPy: a layer's {slot} is a {family} model, a model name, or a config table, not "
                    f"{type(model).__name__}.")


cdef shared_ptr[c_PhysicsBase] cy_model_handle(object model, str what) except *:
    """The shared C++ model behind a wrapper; empty for None. A wrapper holding no C++ model raises ValueError naming
    ``what``, as a model of another family does in the C++ setter, rather than clearing the slot."""
    cdef shared_ptr[c_PhysicsBase] handle
    if model is not None:
        handle = (<PhysicsBase>model)._model_sptr
        if not handle:
            raise ValueError(f"TidalPy: {what} cannot be a {type(model).__name__}.")
    return handle


cdef class Layer(StructureBase):
    """A spherical shell of one material, as the world's EOS and Love solves see it.

    Parameters
    ----------
    name : str
        Layer name (e.g. ``"mantle"``).
    layer_index : int, optional
        Zero-based position in the parent world, innermost layer = 0. Left out, ``BaseWorld.add_layer`` gives the
        layer its place in the stack; given, ``add_layer`` refuses the layer unless it is that place. A standalone
        layer without one reports 0.
    radius_inner : float, optional
        Inner boundary radius [m]. Left out, ``BaseWorld.add_layer`` starts the layer at the top of the stack (0 for
        the first layer); given, ``add_layer`` refuses the layer unless it matches. A standalone layer without one
        starts at 0.
    radius_outer : float
        Outer boundary radius [m]; required (by keyword when the two arguments before it are left out).
    mass : float, optional
        Layer mass [kg]; each world EOS solve sets it from the solved profile. Default ``0.0``.
    material : Material, str, or dict, optional
        The material: a ``Material``, a MatPack name (``"peridotite"``), or a material table (with or without a
        ``preset``). A world's EOS solve needs every layer to have one.
    use_tides : bool, optional
        Whether the layer dissipates tidal energy. Off, it has no share in the quasi-homogeneous Love methods, and the
        radial solver deforms it with its static moduli and no rheology (purely real moduli), so it adds no heating.
        Default ``True``.
    is_volume_fixed : bool, optional
        False lets the layer grow or shrink to hold its mass while the EOS solve redistributes the interior.
        Default ``True``.
    tidal_scale : float, optional
        The layer's share of the planet in the quasi-homogeneous Love methods; ``None`` (default) takes its volume
        fraction.
    state : str, optional
        How the radial solver treats the layer: ``"auto"`` (default) from its material (liquid for a liquid-only
        material, else solid, split into solid and liquid zones where it melts with ``use_melting``), or
        ``"solid"`` or ``"liquid"``.
    is_static : bool, optional
        Static (no inertia) equations; False takes the dynamic form, which a liquid needs at short periods. Default
        ``True``.
    is_incompressible : bool, optional
        Incompressible equations. Default ``False``.
    temperature : float, optional
        Layer temperature [K]. Default ``0.0``, the cold rigid limit of the viscosity laws.
    use_thermal_expansion, use_melting, use_pressure_melting, use_melt_density : bool, optional
        How much of the material the layer uses; each off, the default, is the simpler and faster case (see
        ``Material.calc_state``).
    use_heating : bool, optional
        Let the world's heat sources act inside the layer during a thermal EOS solve. Default ``False``.
    shear_rheology, bulk_rheology : RheologyBase, str, or dict, optional
        Override the material's default rheology.
    cooling : CoolingBase, str, or dict, optional
        How heat moves through the layer in a thermal solve; without one the layer holds one temperature.
    radiogenics : RadiogenicsBase, str, or dict, optional
        The layer's radiogenic heating; without one the layer has none.

    Assumptions
    -----------
    - Spherically symmetric geometry with radius_inner <= radius_outer.
    """

    def __cinit__(self, *args, **kwargs):
        self._is_view   = False
        self._world_ref = None
        self._detached  = False
        self.p_index_given = True
        self.p_inner_given = True

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
            str name,
            layer_index = None,
            radius_inner = None,
            radius_outer = None,
            double mass = 0.0,
            material = None,
            *,
            use_tides = True,
            is_volume_fixed = True,
            tidal_scale = None,
            str state = "auto",
            is_static = True,
            is_incompressible = False,
            double temperature = 0.0,
            use_thermal_expansion = False,
            use_melting = False,
            use_pressure_melting = False,
            use_melt_density = False,
            use_heating = False,
            shear_rheology = None,
            bulk_rheology = None,
            cooling = None,
            radiogenics = None):
        if radius_outer is None:
            raise TypeError(f"TidalPy: layer '{name}' needs its radius_outer [m].")
        cdef c_LayerConfig config
        config.name              = name.encode("utf-8")
        # An index or inner radius left out holds its standalone default until BaseWorld.add_layer fills it in.
        config.layer_index       = 0 if layer_index is None else <int>layer_index
        config.radius_inner      = 0.0 if radius_inner is None else <double>radius_inner
        config.radius_outer      = <double>radius_outer
        config.mass              = mass
        config.use_tides         = bool(use_tides)
        config.is_volume_fixed   = bool(is_volume_fixed)
        config.tidal_scale       = d_NAN if tidal_scale is None else <double>tidal_scale
        config.state             = c_layer_state_from_name(state.encode("utf-8"))
        config.is_static         = bool(is_static)
        config.is_incompressible = bool(is_incompressible)
        config.temperature       = temperature
        config.switches.use_thermal_expansion = bool(use_thermal_expansion)
        config.switches.use_melting           = bool(use_melting)
        config.switches.use_pressure_melting  = bool(use_pressure_melting)
        config.switches.use_melt_density      = bool(use_melt_density)
        config.use_heating       = bool(use_heating)
        self._layer_ptr = make_unique[c_Layer](config)
        self._ptr = <c_TidalPyBaseClass*>self._layer_ptr.get()
        self.p_index_given = layer_index is not None
        self.p_inner_given = radius_inner is not None
        if material is not None:
            self.material = material
        if shear_rheology is not None:
            self.shear_rheology = shear_rheology
        if bulk_rheology is not None:
            self.bulk_rheology = bulk_rheology
        if cooling is not None:
            self.cooling = cooling
        if radiogenics is not None:
            self.radiogenics = radiogenics

    def __dealloc__(self):
        if self._is_view:
            # The C++ layer is owned by the world; relinquish without deleting it.
            self._layer_ptr.release()
        else:
            self._layer_ptr.reset()
        self._ptr = NULL

    cdef void _init_view(self, c_Layer* ptr, object world):
        """Set this wrapper up as a non-owning view onto a world-owned C++ layer.

        Keeps a reference to ``world`` so the C++ layer outlives the view; ``__dealloc__`` then releases the pointer
        instead of deleting it.
        """
        self._layer_ptr.reset(ptr)
        self._ptr       = <c_TidalPyBaseClass*>ptr
        self._is_view   = True
        self._world_ref = world

    @staticmethod
    cdef Layer _view(c_Layer* ptr, object world):
        cdef Layer view = Layer.__new__(Layer)
        view._init_view(ptr, world)
        return view

    def load_binary(self, path, cpp_bool force=False):
        """Load this layer's state from a TidalPy binary file (a ``str`` or ``os.PathLike`` path).

        Only a standalone layer can be loaded: a layer view belongs to its world, whose structure a load would change
        behind its back, so load the world instead. The loaded layer carries the index and inner radius of the file.
        """
        self._check_ptr()
        if self._is_view:
            raise ValueError(
                f"Layer '{self.name}' belongs to a world and cannot be loaded in place: load the world's binary "
                f"file, or load into a standalone layer.")
        StructureBase.load_binary(self, path, force)
        self.p_index_given = True
        self.p_inner_given = True

    def __repr__(self):
        """One line: the name, index, radii [km], state, and temperature [K]."""
        if self._detached or self._ptr is NULL:
            return "Layer(no layer)"
        cdef c_Layer* layer_ptr = self._layer_ptr.get()
        return (
            f"Layer({self.name!r}, index={layer_ptr.get_layer_index()}, "
            f"radius_inner_km={layer_ptr.get_radius_inner() / 1.0e3:.6g}, "
            f"radius_outer_km={layer_ptr.get_radius_outer() / 1.0e3:.6g}, state={self.state!r}, "
            f"temperature_k={layer_ptr.get_temperature():.6g})")

    # =================================================================================================================
    # Geometry and identification
    # =================================================================================================================
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
    def radius(self) -> float:
        """Outer radius [m]."""
        self._check_ptr()
        return self._layer_ptr.get().get_radius_outer()

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
        """Layer thickness [m]."""
        self._check_ptr()
        return self._layer_ptr.get().get_thickness()

    @property
    def volume(self) -> float:
        """Shell volume [m3]."""
        self._check_ptr()
        return self._layer_ptr.get().get_volume()

    @property
    def mass(self) -> float:
        """Layer mass [kg]; each successful world EOS solve sets it from the solved profile."""
        self._check_ptr()
        return self._layer_ptr.get().get_mass()

    @property
    def density_bulk(self) -> float:
        """Bulk density [kg m-3] = mass / volume; NaN for a zero-volume layer."""
        self._check_ptr()
        return self._layer_ptr.get().get_density_bulk()

    @property
    def surface_area_inner(self) -> float:
        """Inner boundary surface area [m2]."""
        self._check_ptr()
        return self._layer_ptr.get().get_surface_area_inner()

    @property
    def surface_area_outer(self) -> float:
        """Outer boundary surface area [m2]."""
        self._check_ptr()
        return self._layer_ptr.get().get_surface_area_outer()

    @property
    def is_volume_fixed(self) -> bool:
        """False if the layer grows or shrinks to hold its mass during an EOS solve."""
        self._check_ptr()
        return self._layer_ptr.get().get_is_volume_fixed()

    @is_volume_fixed.setter
    def is_volume_fixed(self, value):
        self._check_ptr()
        self._layer_ptr.get().set_is_volume_fixed(bool(value))

    def set_radii(self, double radius_inner, double radius_outer):
        """Move the layer's boundaries [m], keeping every derived geometric quantity in step.

        Setting the radii by hand on a world's layer leaves the world's own radius and its other layers untouched, so
        keep the stack continuous; the world forgets its solved profile, which no longer lines up with its layers.

        Raises
        ------
        ValueError
            Unless the radii are finite with 0 <= radius_inner <= radius_outer; the layer is left as it was.
        """
        self._check_ptr()
        self._layer_ptr.get().set_radii(radius_inner, radius_outer)
        self._notify_world_of_move()

    # =================================================================================================================
    # Tides
    # =================================================================================================================
    @property
    def use_tides(self) -> bool:
        """Whether the layer dissipates tidal energy; off, it responds elastically (see the class docstring)."""
        self._check_ptr()
        return self._layer_ptr.get().get_use_tides()

    @use_tides.setter
    def use_tides(self, value):
        self._check_ptr()
        self._layer_ptr.get().set_use_tides(bool(value))

    @property
    def tidal_scale(self):
        """The layer's configured tidal scale [dimensionless], or ``None`` when it takes its volume fraction.

        Used only by the quasi-homogeneous Love methods (``homogeneous``, ``cpl``, ``ctl``) and to share out the heating
        of an analytic tide model. The value in use is ``BaseWorld.get_layer_tidal_scale``.
        """
        self._check_ptr()
        cdef double value = self._layer_ptr.get().get_tidal_scale()
        return None if value != value else value

    @tidal_scale.setter
    def tidal_scale(self, value):
        self._check_ptr()
        self._layer_ptr.get().set_tidal_scale(d_NAN if value is None else <double>value)

    def get_tidal_heating(self) -> float:
        """Tidal heating [W] put in this layer by the world's last ``calc_tides``; NaN before one runs."""
        self._check_ptr()
        return self._layer_ptr.get().get_tidal_heating()

    # =================================================================================================================
    # Material, switches, and state
    # =================================================================================================================
    @property
    def material(self):
        """The layer's material, or None. Set it to a ``Material``, a MatPack name, or a material table."""
        self._check_ptr()
        return cy_wrap_model(self._layer_ptr.get().share_material_model())

    @material.setter
    def material(self, value):
        self._check_ptr()
        self._layer_ptr.get().set_material_model(cy_model_handle(_as_material(value), "a layer's material"))

    @property
    def material_set(self) -> bool:
        """True once the layer has a material."""
        self._check_ptr()
        return self._layer_ptr.get().get_material_set()

    @property
    def temperature(self) -> float:
        """Layer temperature [K]."""
        self._check_ptr()
        return self._layer_ptr.get().get_temperature()

    @temperature.setter
    def temperature(self, double value):
        self._check_ptr()
        self._layer_ptr.get().set_temperature(value)

    def _get_switch(self, int switch_i):
        self._check_ptr()
        cdef const c_MaterialSwitches* switches = &self._layer_ptr.get().get_switches()
        return (switches.use_thermal_expansion, switches.use_melting, switches.use_pressure_melting,
                switches.use_melt_density)[switch_i]

    def _set_switch(self, int switch_i, value):
        self._check_ptr()
        cdef c_MaterialSwitches switches = self._layer_ptr.get().get_switches()
        if switch_i == 0:
            switches.use_thermal_expansion = bool(value)
        elif switch_i == 1:
            switches.use_melting = bool(value)
        elif switch_i == 2:
            switches.use_pressure_melting = bool(value)
        else:
            switches.use_melt_density = bool(value)
        self._layer_ptr.get().set_switches(switches)

    @property
    def use_thermal_expansion(self) -> bool:
        """Whether the density (and thermal pressure) follow the temperature."""
        return self._get_switch(0)

    @use_thermal_expansion.setter
    def use_thermal_expansion(self, value):
        self._set_switch(0, value)

    @property
    def use_melting(self) -> bool:
        """Whether the material's liquid phase, melt fraction, and weakening take part."""
        return self._get_switch(1)

    @use_melting.setter
    def use_melting(self, value):
        self._set_switch(1, value)

    @property
    def use_pressure_melting(self) -> bool:
        """Whether the solidus and liquidus follow the pressure (else they are read at zero pressure)."""
        return self._get_switch(2)

    @use_pressure_melting.setter
    def use_pressure_melting(self, value):
        self._set_switch(2, value)

    @property
    def use_melt_density(self) -> bool:
        """Whether the density mixes the phases' by melt fraction (else it is the solid's)."""
        return self._get_switch(3)

    @use_melt_density.setter
    def use_melt_density(self, value):
        self._set_switch(3, value)

    @property
    def use_heating(self) -> bool:
        """Whether the world's heat sources act inside this layer during a thermal EOS solve."""
        self._check_ptr()
        return self._layer_ptr.get().get_use_heating()

    @use_heating.setter
    def use_heating(self, value):
        self._check_ptr()
        self._layer_ptr.get().set_use_heating(bool(value))

    @property
    def state(self) -> str:
        """How the radial solver treats the layer: ``"auto"``, ``"solid"``, or ``"liquid"``."""
        self._check_ptr()
        return c_layer_state_name(self._layer_ptr.get().get_state()).decode("utf-8")

    @state.setter
    def state(self, str value):
        self._check_ptr()
        self._layer_ptr.get().set_state(c_layer_state_from_name(value.encode("utf-8")))

    @property
    def is_liquid(self) -> bool:
        """Whether the radial solver treats the whole layer as a liquid."""
        self._check_ptr()
        return self._layer_ptr.get().get_is_liquid()

    @property
    def can_change_state(self) -> bool:
        """Whether part of the layer can turn liquid in a solve: an automatic layer whose material melts, with melting
        on."""
        self._check_ptr()
        return self._layer_ptr.get().get_can_change_state()

    @property
    def is_static(self) -> bool:
        """Static (no inertia) equations; False takes the dynamic form. Read by the radial solver."""
        self._check_ptr()
        return self._layer_ptr.get().get_is_static()

    @is_static.setter
    def is_static(self, value):
        self._check_ptr()
        self._layer_ptr.get().set_is_static(bool(value))

    @property
    def is_incompressible(self) -> bool:
        """Incompressible equations. Read by the radial solver."""
        self._check_ptr()
        return self._layer_ptr.get().get_is_incompressible()

    @is_incompressible.setter
    def is_incompressible(self, value):
        self._check_ptr()
        self._layer_ptr.get().set_is_incompressible(bool(value))

    def calc_state(self, pressure, temperature=None, radius=None) -> dict:
        """The layer's material at a point, with the layer's switches (see ``Material.calc_state``).

        Parameters
        ----------
        pressure : float or np.ndarray
            Pressure [Pa].
        temperature : float or np.ndarray, optional
            Temperature [K]; the layer's own when absent.
        radius : float or np.ndarray, optional
            Radius [m], read only by tabulated laws; the layer's mid-radius when absent.

        Raises
        ------
        ValueError
            The layer has no material.
        """
        self._check_ptr()
        material = self.material
        if material is None:
            raise ValueError(f"TidalPy: layer '{self.name}' has no material.")
        cdef const c_MaterialSwitches* switches = &self._layer_ptr.get().get_switches()
        return material.calc_state(
            pressure,
            self.temperature if temperature is None else temperature,
            self._layer_ptr.get().get_radius_mid() if radius is None else radius,
            use_thermal_expansion=switches.use_thermal_expansion,
            use_melting=switches.use_melting,
            use_pressure_melting=switches.use_pressure_melting,
            use_melt_density=switches.use_melt_density)

    # =================================================================================================================
    # Rheology
    # =================================================================================================================
    @property
    def shear_rheology(self):
        """The shear rheology in effect: the layer's override, else its material's default; None for none. Setting it
        sets the override (a model, a model name, or a config table); None clears the override."""
        self._check_ptr()
        return cy_wrap_model(self._layer_ptr.get().share_shear_rheology_model(False))

    @shear_rheology.setter
    def shear_rheology(self, value):
        self._check_ptr()
        self._layer_ptr.get().set_shear_rheology_model(
            cy_model_handle(_as_model(value, "shear_rheology", make_rheology, "rheology"), "a layer's shear rheology"))

    @property
    def bulk_rheology(self):
        """The bulk rheology in effect (see ``shear_rheology``)."""
        self._check_ptr()
        return cy_wrap_model(self._layer_ptr.get().share_bulk_rheology_model(False))

    @bulk_rheology.setter
    def bulk_rheology(self, value):
        self._check_ptr()
        self._layer_ptr.get().set_bulk_rheology_model(
            cy_model_handle(_as_model(value, "bulk_rheology", make_rheology, "rheology"), "a layer's bulk rheology"))

    def _apply_complex(self, radius, double frequency, cpp_bool is_shear):
        self._check_ptr()
        # Radius-resolved complex modulus: float -> complex; np.ndarray -> complex np.ndarray (same shape).
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
            # One C++ call for the whole array, without the GIL, holding the owning world's call lock throughout.
            # double complex and std::complex<double> share one layout.
            if num_radii > 0:
                with nogil:
                    self._layer_ptr.get().calc_complex_moduli(
                        is_shear,
                        &flat_in[0],
                        num_radii,
                        frequency,
                        <cpp_complex[double]*><void*>&flat_out[0])
            return out_arr
        if is_shear:
            value = self._layer_ptr.get().calc_complex_shear_modulus(<double>radius, frequency)
        else:
            value = self._layer_ptr.get().calc_complex_bulk_modulus(<double>radius, frequency)
        return PyComplex_FromDoubles(value.real(), value.imag())

    def calc_complex_shear_modulus(self, first_arg, frequency=None):
        """Complex shear modulus [Pa]: layer-constant or radius-resolved.

        ``calc_complex_shear_modulus(frequency)`` applies the shear rheology to the layer's material at zero pressure
        and its own temperature. ``calc_complex_shear_modulus(radius, frequency)`` uses the post-melt shear modulus and
        viscosity the world's EOS solve stored at ``radius`` (NaN before one).

        Parameters
        ----------
        first_arg : float or np.ndarray
            Forcing frequency [rad s-1] (one-argument form) or radius [m] (two-argument form).
        frequency : float, optional
            Forcing frequency [rad s-1] for the radius-resolved form.
        """
        self._check_ptr()
        cdef cpp_complex[double] result
        if frequency is None:
            result = self._layer_ptr.get().calc_complex_shear_modulus(<double>first_arg)
            return PyComplex_FromDoubles(result.real(), result.imag())
        return self._apply_complex(first_arg, <double>frequency, True)

    def calc_complex_bulk_modulus(self, first_arg, frequency=None):
        """Complex bulk modulus [Pa] from the adiabatic bulk modulus; the same forms as
        ``calc_complex_shear_modulus``."""
        self._check_ptr()
        cdef cpp_complex[double] result
        if frequency is None:
            result = self._layer_ptr.get().calc_complex_bulk_modulus(<double>first_arg)
            return PyComplex_FromDoubles(result.real(), result.imag())
        return self._apply_complex(first_arg, <double>frequency, False)

    # =================================================================================================================
    # Cooling and radiogenics
    # =================================================================================================================
    @property
    def cooling(self):
        """How heat moves through the layer in a thermal solve, or None, which holds the layer at one temperature. Set
        it to a model, a model name, or a config table; None clears it. The model is shared, not consumed."""
        self._check_ptr()
        return cy_wrap_model(self._layer_ptr.get().share_cooling_model())

    @cooling.setter
    def cooling(self, value):
        self._check_ptr()
        self._layer_ptr.get().set_cooling_model(
            cy_model_handle(_as_model(value, "cooling", make_cooling, "cooling"), "a layer's cooling model"))

    @property
    def cooling_set(self) -> bool:
        """True while a cooling model is attached."""
        self._check_ptr()
        return self._layer_ptr.get().get_cooling_model() != NULL

    @property
    def radiogenics(self):
        """The layer's radiogenic heating, or None, which gives none. Set it to a model, a model name, or a config
        table; None clears it. The model is shared, not consumed."""
        self._check_ptr()
        return cy_wrap_model(self._layer_ptr.get().share_radiogenics_model())

    @radiogenics.setter
    def radiogenics(self, value):
        self._check_ptr()
        self._layer_ptr.get().set_radiogenics_model(cy_model_handle(
            _as_model(value, "radiogenics", make_radiogenics, "radiogenics"), "a layer's radiogenics model"))

    @property
    def radiogenics_set(self) -> bool:
        """True while a radiogenics model is attached."""
        self._check_ptr()
        return self._layer_ptr.get().get_radiogenics_model() != NULL

    def calc_radiogenic_heating(self, double time, double mass) -> float:
        """Radiogenic heating [W] at a time [s] for a mass [kg] from the attached model; 0.0 without one."""
        self._check_ptr()
        return self._layer_ptr.get().calc_radiogenic_heating(time, mass)

    # =================================================================================================================
    # Solved profile
    # =================================================================================================================
    @property
    def eos_data_populated(self) -> bool:
        """True after the profile has been populated (world EOS solve or ``update_eos_data``)."""
        self._check_ptr()
        return self._layer_ptr.get().get_eos_data_populated()

    @property
    def viscoelastic_populated(self) -> bool:
        """True after the world EOS solve has populated this layer's material state."""
        self._check_ptr()
        return self._layer_ptr.get().get_viscoelastic_populated()

    def update_eos_data(self, radius, density_kgm3, gravity_ms2, pressure):
        """Populate the profile directly from arrays (radius [m] ascending, density [kg m-3], gravity [m s-2], pressure
        [Pa], one length), bypassing the world's EOS solve. Such a profile carries no material state."""
        self._check_ptr()
        cdef vector[double] r_vec   = radius
        cdef vector[double] rho_vec = density_kgm3
        cdef vector[double] g_vec   = gravity_ms2
        cdef vector[double] p_vec   = pressure
        cdef c_LayerEOSData eos_data
        eos_data.populate(r_vec, rho_vec, g_vec, p_vec)
        self._layer_ptr.get().update_eos_data(eos_data)

    # A float radius [m] gives a float, an np.ndarray a same-shape array, read in one C++ call under the owning
    # world's call lock (see cy_eos_field).
    def get_density(self, radius):
        """Density [kg m-3] at radius [m]; NaN if the profile is not populated."""
        self._check_ptr()
        return cy_eos_field(<const void*>self._layer_ptr.get(), cy_layer_eos_fields, radius, C_EOS_DENSITY_INDEX)

    def get_gravity(self, radius):
        """Gravitational acceleration [m s-2] at radius [m]; NaN if the profile is not populated."""
        self._check_ptr()
        return cy_eos_field(<const void*>self._layer_ptr.get(), cy_layer_eos_fields, radius, C_EOS_GRAVITY_INDEX)

    def get_pressure(self, radius):
        """Pressure [Pa] at radius [m]; NaN if the profile is not populated."""
        self._check_ptr()
        return cy_eos_field(<const void*>self._layer_ptr.get(), cy_layer_eos_fields, radius, C_EOS_PRESSURE_INDEX)

    def get_shear_modulus(self, radius):
        """Post-melt static shear modulus [Pa] at radius [m]; NaN if the material state is not populated."""
        self._check_ptr()
        return cy_eos_field(
            <const void*>self._layer_ptr.get(), cy_layer_eos_fields, radius, C_EOS_SHEAR_MODULUS_INDEX)

    def get_bulk_modulus(self, radius):
        """Post-melt adiabatic bulk modulus [Pa] at radius [m]; NaN if the material state is not populated."""
        self._check_ptr()
        return cy_eos_field(
            <const void*>self._layer_ptr.get(), cy_layer_eos_fields, radius, C_EOS_BULK_MODULUS_INDEX)

    def get_shear_viscosity(self, radius):
        """Post-melt shear viscosity [Pa s] at radius [m]; NaN if the material state is not populated."""
        self._check_ptr()
        return cy_eos_field(
            <const void*>self._layer_ptr.get(), cy_layer_eos_fields, radius, C_EOS_SHEAR_VISCOSITY_INDEX)

    def get_bulk_viscosity(self, radius):
        """Post-melt bulk viscosity [Pa s] at radius [m]; NaN if the material state is not populated."""
        self._check_ptr()
        return cy_eos_field(
            <const void*>self._layer_ptr.get(), cy_layer_eos_fields, radius, C_EOS_BULK_VISCOSITY_INDEX)

    def get_melt_fraction(self, radius):
        """Melt fraction at radius [m]: 0 for a layer that does not melt, 1 for a liquid-only material; NaN if the
        material state is not populated."""
        self._check_ptr()
        return cy_eos_field(
            <const void*>self._layer_ptr.get(), cy_layer_eos_fields, radius, C_EOS_MELT_FRACTION_INDEX)

    def get_static_viscoelastics(self, radius):
        """``(shear_modulus, shear_viscosity, bulk_modulus, bulk_viscosity)`` (post-melt) at radius [m], each a float or
        np.ndarray, from one evaluation of the solved state per radius."""
        self._check_ptr()
        return cy_eos_fields(
            <const void*>self._layer_ptr.get(), cy_layer_eos_fields, radius,
            (C_EOS_SHEAR_MODULUS_INDEX, C_EOS_SHEAR_VISCOSITY_INDEX,
             C_EOS_BULK_MODULUS_INDEX, C_EOS_BULK_VISCOSITY_INDEX))

    def get_state(self, radius):
        """Every solved profile at radius as a dict (float or np.ndarray values), from one evaluation per radius."""
        self._check_ptr()
        values = cy_eos_fields(
            <const void*>self._layer_ptr.get(), cy_layer_eos_fields, radius,
            (C_EOS_DENSITY_INDEX, C_EOS_GRAVITY_INDEX, C_EOS_PRESSURE_INDEX, C_EOS_SHEAR_MODULUS_INDEX,
             C_EOS_SHEAR_VISCOSITY_INDEX, C_EOS_BULK_MODULUS_INDEX, C_EOS_BULK_VISCOSITY_INDEX,
             C_EOS_MELT_FRACTION_INDEX))
        return dict(zip(("density", "gravity", "pressure", "shear_modulus", "shear_viscosity", "bulk_modulus",
                         "bulk_viscosity", "melt_fraction"), values))

    # =================================================================================================================
    # Configuration
    # =================================================================================================================
    cpdef dict get_config_dict(self):
        """All configuration values as a dict (MKS) in the world builder's layer schema.

        Returns
        -------
        dict
            ``name``, ``layer_index``, ``radius_inner_m``, ``radius_outer_m``, ``mass_kg`` once the layer has a mass
            (one given, or set by a world EOS solve; a mass of 0.0 is the unset default and is left out),
            ``use_tides``, ``is_volume_fixed``, ``tidal_scale`` when one is set, ``state``, ``is_static``,
            ``is_incompressible``,
            ``temperature_k``, the four material switches, ``use_heating``, and the ``material``, ``shear_rheology``
            and ``bulk_rheology`` (the layer's overrides), ``cooling``, and ``radiogenics`` tables when set. ``name``
            and ``radius_inner_m`` are standalone-layer keys that a world drops when it nests the layer (see
            ``LAYER_STANDALONE_CONFIG_KEYS``).
        """
        self._check_ptr()
        cdef c_Layer* p = self._layer_ptr.get()
        cdef double tidal_scale = p.get_tidal_scale()
        cdef const c_MaterialSwitches* switches = &p.get_switches()
        cdef dict config = {
            "name":              p.get_name().decode("utf-8"),
            "layer_index":       p.get_layer_index(),
            "radius_inner_m":    p.get_radius_inner(),
            "radius_outer_m":    p.get_radius_outer(),
        }
        if p.get_mass() != 0.0:
            config["mass_kg"] = p.get_mass()
        config["use_tides"]       = bool(p.get_use_tides())
        config["is_volume_fixed"] = bool(p.get_is_volume_fixed())
        if tidal_scale == tidal_scale:
            config["tidal_scale"] = tidal_scale
        config["state"]                 = c_layer_state_name(p.get_state()).decode("utf-8")
        config["is_static"]             = bool(p.get_is_static())
        config["is_incompressible"]     = bool(p.get_is_incompressible())
        config["temperature_k"]         = p.get_temperature()
        config["use_thermal_expansion"] = bool(switches.use_thermal_expansion)
        config["use_melting"]           = bool(switches.use_melting)
        config["use_pressure_melting"]  = bool(switches.use_pressure_melting)
        config["use_melt_density"]      = bool(switches.use_melt_density)
        config["use_heating"]           = bool(p.get_use_heating())
        material = self.material
        if material is not None:
            config["material"] = material.get_config_dict()
        for key, model in (("shear_rheology", cy_wrap_model(p.share_shear_rheology_model(True))),
                           ("bulk_rheology", cy_wrap_model(p.share_bulk_rheology_model(True)))):
            if model is not None:
                config[key] = model.get_config_dict()
        cdef const c_PhysicsBase* model_ptr
        model_ptr = <const c_PhysicsBase*>p.get_cooling_model()
        if model_ptr != NULL:
            config["cooling"] = cy_physics_model_config(model_ptr)
        model_ptr = <const c_PhysicsBase*>p.get_radiogenics_model()
        if model_ptr != NULL:
            config["radiogenics"] = cy_physics_model_config(model_ptr)
        return config
