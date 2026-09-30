# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrapper for TidalPy's base layer class.

BaseLayer holds a layer's geometry, material, radial-solver flags, temperature, and rheology. Attaching a rheology
gives frequency-dependent complex moduli; without one the static modulus is returned as a real-valued complex number.
Its EOS profile stays unpopulated until the world's EOS solve runs or ``update_eos_data`` is called directly.
"""

import os as _os

cimport numpy as cnp
cnp.import_array()

import numpy as np

from libcpp.vector cimport vector
from libcpp.complex cimport complex as cpp_complex
from libcpp cimport bool as cpp_bool
from libcpp.utility cimport move
from libcpp.memory cimport make_unique

from TidalPy.Utilities.logging.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport d_NAN, set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities.classes.classes cimport (
    StructureBase,
    c_TidalPyBaseClass,
    c_PhysicsBase,
    cy_physics_model_config,
)
from TidalPy.Material.eos.material_eos cimport MaterialEOSBase, cy_material_config
from TidalPy.Rheology.rheology cimport RheologyBase
from TidalPy.Viscosity.viscosity cimport ViscosityBase
from TidalPy.PartialMelt.partial_melt cimport PartialMeltBase

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())

# Layer config keys that are constructor parameters of a standalone layer but not part of the world
# builder's layer schema: the world writer drops them when nesting a layer under its name.
LAYER_STANDALONE_CONFIG_KEYS = (
    "name",
    "radius_inner_m",
)


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


cdef int cy_fill_base_layer_config(
        c_BaseLayerConfig* config,
        str name,
        int layer_index,
        double radius_inner,
        double radius_outer,
        double mass,
        str material_name,
        cpp_bool is_tidal,
        cpp_bool is_volume_fixed,
        object tidal_scale,
        cpp_bool is_solid,
        cpp_bool is_static,
        cpp_bool is_incompressible,
        double temperature,
        cpp_bool use_thermal_eos,
        cpp_bool use_heating) except -1:
    config.name            = name.encode("utf-8")
    config.layer_index     = layer_index
    config.radius_inner    = radius_inner
    config.radius_outer    = radius_outer
    config.mass            = mass
    config.material_name   = material_name.encode("utf-8")
    config.is_tidal        = is_tidal
    config.is_volume_fixed = is_volume_fixed
    config.tidal_scale     = d_NAN if tidal_scale is None else <double>tidal_scale
    config.is_solid          = is_solid
    config.is_static         = is_static
    config.is_incompressible = is_incompressible
    config.temperature       = temperature
    config.use_thermal_eos   = use_thermal_eos
    config.use_heating       = use_heating
    return 0


cdef void cy_layer_eos_fields(
        const void* owner,
        const size_t* field_indices,
        size_t num_fields,
        const double* radii,
        size_t num_radii,
        double* values_out) noexcept nogil:
    (<const c_BaseLayer*>owner).get_eos_fields(field_indices, num_fields, radii, num_radii, values_out)


cdef class BaseLayer(StructureBase):
    """Layer: geometry, material, radial-solver flags, temperature, and optional shear and bulk rheology.

    Rheology models can be attached to give frequency-dependent complex moduli; without one the static modulus is
    returned as a real-valued complex number. The EOS profile (density, gravity, pressure against radius) starts
    unpopulated and becomes queryable after the world's EOS solve or a direct ``update_eos_data`` call.

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
    is_solid : bool, optional
        False marks the layer liquid for the radial Love-number solver. Default ``True``.
    is_static : bool, optional
        Use the static (no inertia) approximation in the radial solver. Default ``True``, so a liquid layer
        is a static liquid unless this is set False.
    is_incompressible : bool, optional
        Use the incompressible approximation in the radial solver. Default ``False``.
    temperature : float, optional
        Layer temperature [K] at which its viscosity and melt models are evaluated. Default ``0.0``, the cold
        rigid limit of the viscosity laws.
    use_thermal_eos : bool, optional
        Pass the temperature to the EOS model, so the density and bulk modulus depend on it. Default ``False``.
    use_heating : bool, optional
        Let the world's heat sources (a solid-liquid layer's radiogenics model among them) act inside the layer
        during a thermal EOS solve. Default ``False``.

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
            str    material_name       = "",
            cpp_bool is_tidal          = True,
            cpp_bool is_volume_fixed   = True,
            tidal_scale                = None,
            cpp_bool is_solid          = True,
            cpp_bool is_static         = True,
            cpp_bool is_incompressible = False,
            double temperature         = 0.0,
            cpp_bool use_thermal_eos   = False,
            cpp_bool use_heating       = False):
        cdef c_BaseLayerConfig config
        cy_fill_base_layer_config(
            &config, name, layer_index, radius_inner, radius_outer, mass, material_name, is_tidal, is_volume_fixed,
            tidal_scale, is_solid, is_static, is_incompressible, temperature, use_thermal_eos, use_heating)
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
        returns to the volume fraction. The value in use is ``BaseWorld.get_layer_tidal_scale``.
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

        Set by :meth:`BaseWorld.calc_tides`; how the heating is resolved per layer depends on the world's Love
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

        Ownership of the C++ model moves out of ``eos``, which is left an empty shell and must not be reused. The
        layer's viscosity and partial-melt models are held by its material, so replacing the material keeps the ones
        attached before (with ``set_shear_viscosity``, ``set_bulk_viscosity``, ``set_partial_melt``, or in the
        material's own config) unless ``eos`` carries its own. A layer in a solved world takes the new material at the
        next ``solve_eos``.

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

    @property
    def shear_modulus_static(self) -> float:
        """Unrelaxed (static) shear modulus [Pa]."""
        self._check_ptr()
        return self._layer_ptr.get().get_shear_modulus_static()

    @property
    def bulk_modulus_static(self) -> float:
        """Unrelaxed (static) bulk modulus [Pa]."""
        self._check_ptr()
        return self._layer_ptr.get().get_bulk_modulus_static()

    @property
    def shear_viscosity_static(self) -> float:
        """Reference dynamic shear viscosity [Pa·s]."""
        self._check_ptr()
        return self._layer_ptr.get().get_shear_viscosity_static()

    @property
    def bulk_viscosity_static(self) -> float:
        """Reference dynamic bulk viscosity [Pa·s]."""
        self._check_ptr()
        return self._layer_ptr.get().get_bulk_viscosity_static()

    @property
    def shear_rheology_set(self) -> bool:
        """True if a shear rheology model has been attached."""
        self._check_ptr()
        return self._layer_ptr.get().get_shear_rheology_set()

    @property
    def bulk_rheology_set(self) -> bool:
        """True if a bulk rheology model has been attached."""
        self._check_ptr()
        return self._layer_ptr.get().get_bulk_rheology_set()

    @property
    def is_solid(self) -> bool:
        """True if this layer is solid, False for liquid. Used by the radial Love-number solver."""
        self._check_ptr()
        return bool(self._layer_ptr.get().get_is_solid())

    @is_solid.setter
    def is_solid(self, value: bool):
        self._check_ptr()
        self._layer_ptr.get().set_is_solid(<cpp_bool>bool(value))

    @property
    def is_static(self) -> bool:
        """True if the static (no-dynamic-terms) approximation is used. Used by the radial solver."""
        self._check_ptr()
        return bool(self._layer_ptr.get().get_is_static())

    @is_static.setter
    def is_static(self, value: bool):
        self._check_ptr()
        self._layer_ptr.get().set_is_static(<cpp_bool>bool(value))

    @property
    def is_incompressible(self) -> bool:
        """True if the incompressible approximation is used. Used by the radial solver."""
        self._check_ptr()
        return bool(self._layer_ptr.get().get_is_incompressible())

    @is_incompressible.setter
    def is_incompressible(self, value: bool):
        self._check_ptr()
        self._layer_ptr.get().set_is_incompressible(<cpp_bool>bool(value))

    @property
    def temperature(self) -> float:
        """Layer temperature [K] at which the viscosity and melt models are evaluated."""
        self._check_ptr()
        return self._layer_ptr.get().get_temperature()

    @temperature.setter
    def temperature(self, double value):
        self._check_ptr()
        self._layer_ptr.get().set_temperature(value)

    @property
    def use_thermal_eos(self) -> bool:
        """True if the EOS model receives the temperature (thermal density and bulk modulus)."""
        self._check_ptr()
        return bool(self._layer_ptr.get().get_use_thermal_eos())

    @use_thermal_eos.setter
    def use_thermal_eos(self, value: bool):
        self._check_ptr()
        self._layer_ptr.get().set_use_thermal_eos(<cpp_bool>bool(value))

    @property
    def use_heating(self) -> bool:
        """True if the world's heat sources act inside this layer during a thermal EOS solve."""
        self._check_ptr()
        return bool(self._layer_ptr.get().get_use_heating())

    @use_heating.setter
    def use_heating(self, value: bool):
        self._check_ptr()
        self._layer_ptr.get().set_use_heating(<cpp_bool>bool(value))

    def set_shear_rheology(self, RheologyBase rheology not None):
        """Attach a rheology model used to compute the complex shear modulus.

        Ownership of the C++ model moves out of ``rheology``, which is left an empty shell and must not be reused.

        Parameters
        ----------
        rheology : RheologyBase
            A rheology model (e.g. ``Maxwell()``, ``make_rheology("andrade")``).

        Raises
        ------
        ValueError
            If ``rheology`` has already been attached or otherwise moved.
        """
        self._check_ptr()
        if rheology._rheology_ptr.get() == NULL:
            raise ValueError(
                "This rheology model holds no C++ object (already attached or moved).")
        self._layer_ptr.get().set_shear_rheology(move(rheology._rheology_ptr))
        rheology._ptr = NULL

    def set_bulk_rheology(self, RheologyBase rheology not None):
        """Attach a rheology model used to compute the complex bulk modulus.

        Ownership of the C++ model moves out of ``rheology``, which is left an empty shell and must not be reused.

        Parameters
        ----------
        rheology : RheologyBase
            A rheology model (e.g. ``Maxwell()``, ``make_rheology("andrade")``).

        Raises
        ------
        ValueError
            If ``rheology`` has already been attached or otherwise moved.
        """
        self._check_ptr()
        if rheology._rheology_ptr.get() == NULL:
            raise ValueError(
                "This rheology model holds no C++ object (already attached or moved).")
        self._layer_ptr.get().set_bulk_rheology(move(rheology._rheology_ptr))
        rheology._ptr = NULL

    @property
    def shear_viscosity_set(self) -> bool:
        """True if a shear viscosity model has been attached."""
        self._check_ptr()
        return self._layer_ptr.get().get_shear_viscosity_set()

    @property
    def bulk_viscosity_set(self) -> bool:
        """True if a bulk viscosity model has been attached."""
        self._check_ptr()
        return self._layer_ptr.get().get_bulk_viscosity_set()

    @property
    def partial_melt_set(self) -> bool:
        """True if a partial-melt model has been attached."""
        self._check_ptr()
        return self._layer_ptr.get().get_partial_melt_set()

    def set_shear_viscosity(self, ViscosityBase viscosity not None):
        """Give the layer's material a viscosity model for its shear viscosity (before partial melt).

        A helper: the material owns this model, so it is handed to the layer's EOS model, which must already be
        attached (``set_eos`` first). Ownership of the C++ model moves out of the argument, which is left an empty
        shell and must not be reused.

        Raises
        ------
        ValueError
            If the model has already been attached or otherwise moved, or the layer has no EOS model yet.
        """
        self._check_ptr()
        if viscosity._visc_ptr.get() == NULL:
            raise ValueError(
                "This viscosity model holds no C++ object (already attached or moved).")
        if not self._layer_ptr.get().get_eos_set():
            raise ValueError(
                f"Attach an EOS model to layer '{self.name}' before giving it a viscosity model: the material owns it.")
        self._layer_ptr.get().set_shear_viscosity(move(viscosity._visc_ptr))
        viscosity._ptr = NULL

    def set_bulk_viscosity(self, ViscosityBase viscosity not None):
        """Give the layer's material a viscosity model for its bulk viscosity (before partial melt).

        A helper: the material owns this model, so it is handed to the layer's EOS model, which must already be
        attached (``set_eos`` first). Ownership of the C++ model moves out of the argument, which is left an empty
        shell and must not be reused.

        Raises
        ------
        ValueError
            If the model has already been attached or otherwise moved, or the layer has no EOS model yet.
        """
        self._check_ptr()
        if viscosity._visc_ptr.get() == NULL:
            raise ValueError(
                "This viscosity model holds no C++ object (already attached or moved).")
        if not self._layer_ptr.get().get_eos_set():
            raise ValueError(
                f"Attach an EOS model to layer '{self.name}' before giving it a viscosity model: the material owns it.")
        self._layer_ptr.get().set_bulk_viscosity(move(viscosity._visc_ptr))
        viscosity._ptr = NULL

    def set_partial_melt(self, PartialMeltBase partial_melt not None):
        """Give the layer's material a partial-melt model that weakens its static moduli and viscosities.

        A helper: the material owns this model, so it is handed to the layer's EOS model, which must already be
        attached (``set_eos`` first). Ownership of the C++ model moves out of the argument, which is left an empty
        shell and must not be reused.

        Raises
        ------
        ValueError
            If the model has already been attached or otherwise moved, or the layer has no EOS model yet.
        """
        self._check_ptr()
        if partial_melt._melt_ptr.get() == NULL:
            raise ValueError(
                "This partial-melt model holds no C++ object (already attached or moved).")
        if not self._layer_ptr.get().get_eos_set():
            raise ValueError(
                f"Attach an EOS model to layer '{self.name}' before giving it a partial-melt model: "
                "the material owns it.")
        self._layer_ptr.get().set_partial_melt(move(partial_melt._melt_ptr))
        partial_melt._ptr = NULL

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
            # One C++ call for the whole array, without the GIL, holding the owning world's call lock throughout
            # so the read takes turns with a solve_eos on another thread. double complex and std::complex<double>
            # share one layout.
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
        return complex(value.real(), value.imag())

    def calc_complex_shear_modulus(self, first_arg, frequency=None):
        """Complex shear modulus [Pa]: layer-constant or radius-resolved.

        ``calc_complex_shear_modulus(frequency)`` applies the shear rheology to the layer-constant static shear
        modulus and viscosity. That static viscosity is the material's (NaN unless it sets one), so a viscous
        rheology then returns NaN: set it explicitly or use the radius-resolved form after the world EOS solve.
        ``calc_complex_shear_modulus(radius, frequency)`` instead uses the post-melt static modulus and viscosity
        stored at ``radius`` by that solve.

        Parameters
        ----------
        first_arg : float or np.ndarray
            Tidal forcing frequency [rad/s] (one-argument form) or query radius [m] (two-argument form).
        frequency : float, optional
            Tidal forcing frequency [rad/s] for the radius-resolved form.

        Returns
        -------
        complex or np.ndarray
            Complex shear modulus [Pa]; a complex ndarray for an array of radii. The radius-resolved form is NaN
            until the world EOS solve populates the layer.

        Assumptions
        -----------
        - Linear viscoelastic response at a single forcing frequency.
        """
        self._check_ptr()
        cdef cpp_complex[double] result
        if frequency is None:
            result = self._layer_ptr.get().calc_complex_shear_modulus(<double>first_arg)
            return complex(result.real(), result.imag())
        return self._apply_complex(first_arg, <double>frequency, True)

    def calc_complex_bulk_modulus(self, first_arg, frequency=None):
        """Complex bulk modulus [Pa]: layer-constant or radius-resolved.

        ``calc_complex_bulk_modulus(frequency)`` applies the bulk rheology to the layer-constant static bulk
        modulus and viscosity. That static viscosity is the material's (NaN unless it sets one), so a viscous
        rheology then returns NaN: set it explicitly or use the radius-resolved form after the world EOS solve.
        ``calc_complex_bulk_modulus(radius, frequency)`` instead uses the post-melt static modulus and viscosity
        stored at ``radius`` by that solve.

        Parameters
        ----------
        first_arg : float or np.ndarray
            Tidal forcing frequency [rad/s] (one-argument form) or query radius [m] (two-argument form).
        frequency : float, optional
            Tidal forcing frequency [rad/s] for the radius-resolved form.

        Returns
        -------
        complex or np.ndarray
            Complex bulk modulus [Pa]; a complex ndarray for an array of radii. The radius-resolved form is NaN
            until the world EOS solve populates the layer.

        Assumptions
        -----------
        - Linear viscoelastic response at a single forcing frequency.
        """
        self._check_ptr()
        cdef cpp_complex[double] result
        if frequency is None:
            result = self._layer_ptr.get().calc_complex_bulk_modulus(<double>first_arg)
            return complex(result.real(), result.imag())
        return self._apply_complex(first_arg, <double>frequency, False)

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

    # The profile getters: a float radius [m] gives a float, an np.ndarray a same-shape array, read in one C++ call
    # under the owning world's call lock (see cy_eos_field).
    def get_density(self, radius):
        """Density [kg/m^3] at radius [m] (float or np.ndarray); NaN if EOS data not populated."""
        self._check_ptr()
        return cy_eos_field(<const void*>self._layer_ptr.get(), cy_layer_eos_fields, radius, C_EOS_DENSITY_INDEX)

    def get_gravity(self, radius):
        """Gravitational acceleration [m/s^2] at radius [m] (float or np.ndarray); NaN if not populated."""
        self._check_ptr()
        return cy_eos_field(<const void*>self._layer_ptr.get(), cy_layer_eos_fields, radius, C_EOS_GRAVITY_INDEX)

    def get_pressure(self, radius):
        """Pressure [Pa] at radius [m] (float or np.ndarray); NaN if EOS data not populated."""
        self._check_ptr()
        return cy_eos_field(<const void*>self._layer_ptr.get(), cy_layer_eos_fields, radius, C_EOS_PRESSURE_INDEX)

    # Viscoelastic profile (populated by the world EOS solve; NaN before then)
    @property
    def viscoelastic_populated(self) -> bool:
        """True after the world EOS solve has populated this layer's viscoelastic state."""
        self._check_ptr()
        return self._layer_ptr.get().get_viscoelastic_populated()

    def get_shear_modulus(self, radius):
        """Post-melt static shear modulus [Pa] at radius [m] (float or np.ndarray); NaN if unpopulated."""
        self._check_ptr()
        return cy_eos_field(
            <const void*>self._layer_ptr.get(), cy_layer_eos_fields, radius, C_EOS_SHEAR_MODULUS_INDEX)

    def get_bulk_modulus(self, radius):
        """Post-melt static bulk modulus [Pa] at radius [m] (float or np.ndarray); NaN if unpopulated."""
        self._check_ptr()
        return cy_eos_field(
            <const void*>self._layer_ptr.get(), cy_layer_eos_fields, radius, C_EOS_BULK_MODULUS_INDEX)

    def get_shear_viscosity(self, radius):
        """Post-melt shear viscosity [Pa s] at radius [m] (float or np.ndarray); NaN if unpopulated."""
        self._check_ptr()
        return cy_eos_field(
            <const void*>self._layer_ptr.get(), cy_layer_eos_fields, radius, C_EOS_SHEAR_VISCOSITY_INDEX)

    def get_bulk_viscosity(self, radius):
        """Post-melt bulk viscosity [Pa s] at radius [m] (float or np.ndarray); NaN if unpopulated."""
        self._check_ptr()
        return cy_eos_field(
            <const void*>self._layer_ptr.get(), cy_layer_eos_fields, radius, C_EOS_BULK_VISCOSITY_INDEX)

    def get_melt_fraction(self, radius):
        """Melt fraction at radius [m] (float or np.ndarray) from the attached partial-melt model.

        0.0 where no partial-melt model is attached; NaN if unpopulated.
        """
        self._check_ptr()
        return cy_eos_field(
            <const void*>self._layer_ptr.get(), cy_layer_eos_fields, radius, C_EOS_MELT_FRACTION_INDEX)

    # Shorthand bundles (one call returns several profiles at once; mirrors the world-level surface)
    def get_static_viscoelastics(self, radius):
        """``(shear_modulus, shear_viscosity, bulk_modulus, bulk_viscosity)`` (post-melt) at radius [m], each a
        float or np.ndarray. One evaluation of the solved state per radius fills all four.
        """
        self._check_ptr()
        return cy_eos_fields(
            <const void*>self._layer_ptr.get(), cy_layer_eos_fields, radius,
            (C_EOS_SHEAR_MODULUS_INDEX, C_EOS_SHEAR_VISCOSITY_INDEX,
             C_EOS_BULK_MODULUS_INDEX, C_EOS_BULK_VISCOSITY_INDEX))

    def get_state(self, radius):
        """All EOS-related profiles at radius as a dict (float or np.ndarray values), from one evaluation of the
        solved state per radius."""
        self._check_ptr()
        values = cy_eos_fields(
            <const void*>self._layer_ptr.get(), cy_layer_eos_fields, radius,
            (C_EOS_DENSITY_INDEX, C_EOS_GRAVITY_INDEX, C_EOS_PRESSURE_INDEX, C_EOS_SHEAR_MODULUS_INDEX,
             C_EOS_SHEAR_VISCOSITY_INDEX, C_EOS_BULK_MODULUS_INDEX, C_EOS_BULK_VISCOSITY_INDEX,
             C_EOS_MELT_FRACTION_INDEX))
        return dict(zip(("density", "gravity", "pressure", "shear_modulus", "shear_viscosity", "bulk_modulus",
                         "bulk_viscosity", "melt_fraction"), values))

    cpdef dict get_config_dict(self):
        """Return all configuration values as a Python dict (MKS) in the world builder's layer schema.

        Each attached physics model appears as its own sub-table keyed by ``model``: the material (the EOS model,
        which carries the static constants, the shear law, and its viscosity and partial-melt models) under
        ``material``, and each rheology under ``shear_rheology`` and ``bulk_rheology``. The material ``type`` is
        written as ``"none"`` so that a rebuild does not add the material defaults a typeless layer would take.
        ``name`` and ``radius_inner`` are standalone-layer keys that a world drops when it nests the layer (see
        ``LAYER_STANDALONE_CONFIG_KEYS``).

        Returns
        -------
        dict
            Keys: ``class``, ``type``, ``name``, ``layer_index``, ``radius_inner``, ``radius_outer``, ``mass``,
            ``material_name``, ``is_tidal``, ``is_volume_fixed``, ``tidal_scale`` when one is set, the radial-solver
            flags ``is_solid``, ``is_static``, and ``is_incompressible``, ``temperature_k``, ``use_thermal_eos``,
            ``use_heating``, and the ``material``, ``shear_rheology``, and ``bulk_rheology`` tables when set.
        """
        # Deferred: the configs package imports the layer modules.
        from TidalPy.Structures.configs.toml_loader import NO_MATERIAL_TYPE

        # The subclasses call this first, so this check covers their typed pointers too.
        self._check_ptr()
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
        config["is_solid"]          = bool(p.get_is_solid())
        config["is_static"]         = bool(p.get_is_static())
        config["is_incompressible"] = bool(p.get_is_incompressible())
        config["temperature_k"]     = p.get_temperature()
        config["use_thermal_eos"]   = bool(p.get_use_thermal_eos())
        config["use_heating"]       = bool(p.get_use_heating())
        if p.get_eos_set():
            config["material"] = cy_material_config(p.get_eos())
        # Attached models, keyed the way the world builder reads them.
        cdef const c_PhysicsBase* model_ptr
        model_ptr = <const c_PhysicsBase*>p.get_shear_rheology_model()
        if model_ptr != NULL:
            config["shear_rheology"] = cy_physics_model_config(model_ptr)
        model_ptr = <const c_PhysicsBase*>p.get_bulk_rheology_model()
        if model_ptr != NULL:
            config["bulk_rheology"] = cy_physics_model_config(model_ptr)
        return config
