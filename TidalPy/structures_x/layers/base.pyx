# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
base.pyx
Cython/Python wrapper for TidalPy's base layer class.

BaseLayer: geometry-only layer storing inner/outer radii, mass, and
material identification. EOS profile (density, gravity, pressure vs. radius)
is unpopulated until the world's EOS solve runs or update_eos_data() is called
directly (e.g., for testing).
"""

import os as _os

from libcpp.vector cimport vector
from libcpp cimport bool as cpp_bool
from libcpp.utility cimport move

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport StructureBase, c_TidalPyBaseClass
from TidalPy.Material_x.eos.material_eos cimport MaterialEOSBase

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


# =====================================================================================================================
# BaseLayer
# =====================================================================================================================

cdef class BaseLayer(StructureBase):
    """Geometry base layer: inner/outer radii, mass, and material identity.

    All spatial data in meters [m] and mass in kilograms [kg] (MKS).
    Derived geometry (thickness, volume, surface areas) is computed at
    construction and accessible as read-only properties.

    The EOS profile (density, gravity, pressure vs. radius) starts unpopulated.
    After the world's EOS solve populates it, or after a direct call to
    ``update_eos_data``, the profile is interpolatable.

    Parameters
    ----------
    name : str
        Human-readable layer name (e.g. ``"mantle"``).
    layer_index : int
        Zero-based position in the parent world, innermost layer = 0.
    radius_inner_m : float
        Inner boundary radius [m].
    radius_outer_m : float
        Outer boundary radius [m].
    mass_kg : float
        Total layer mass [kg].
    material_name : str, optional
        Material identifier string (e.g. ``"perovskite"``). Default ``""``.
    is_tidal : bool, optional
        Whether this layer contributes to tidal dissipation. Default ``True``.
    tidal_scale : float, optional
        Dimensionless scale factor applied to this layer's tidal heating.
        Default ``1.0``.

    Assumptions
    -----------
    - Layer geometry is spherically symmetric.
    - radius_inner_m <= radius_outer_m.
    - All values are in MKS units.
    """

    def __cinit__(self, *args, **kwargs):
        # unique_ptr<c_BaseLayer> auto-inits to nullptr; _ptr set in __init__.
        self._is_view   = False
        self._world_ref = None

    def __init__(
            self,
            str    name,
            int    layer_index,
            double radius_inner_m,
            double radius_outer_m,
            double mass_kg,
            str    material_name      = "",
            cpp_bool   is_tidal           = True,
            double tidal_scale        = 1.0,
            str    tidal_scale_method = "user_provided"):
        cdef c_BaseLayerConfig config
        config.name               = name.encode("utf-8")
        config.layer_index        = layer_index
        config.radius_inner_m     = radius_inner_m
        config.radius_outer_m     = radius_outer_m
        config.mass_kg            = mass_kg
        config.material_name      = material_name.encode("utf-8")
        config.is_tidal           = is_tidal
        config.tidal_scale        = tidal_scale
        config.tidal_scale_method = c_tidal_scale_method_from_name(tidal_scale_method.encode("utf-8"))
        self._layer_ptr.reset(new c_BaseLayer(config))
        self._ptr = <c_TidalPyBaseClass*>self._layer_ptr.get()

    def __dealloc__(self):
        if self._is_view:
            # The C++ layer is owned by the world; relinquish without deleting it.
            self._layer_ptr.release()
        else:
            self._layer_ptr.reset()
        self._ptr = NULL

    # ------------------------------------------------------------------------------------------------------------------
    # Non-owning view construction (used by a world to hand back a wrapper around a layer it owns)
    # ------------------------------------------------------------------------------------------------------------------
    cdef void _init_view(self, c_BaseLayer* ptr, object world):
        """Set this wrapper up as a non-owning view onto a world-owned C++ layer.

        Stores ``ptr`` in ``_layer_ptr`` so all base getters work unchanged, but marks the
        wrapper as a view so ``__dealloc__`` releases (does not delete) it, and keeps a
        reference to the owning ``world`` so the C++ layer outlives the view.
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

    # ------------------------------------------------------------------------------------------------------------------
    # Base class property overrides
    # (_layer_ptr always points to the actual most-derived C++ layer; safe for subclasses)
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def radius(self) -> float:
        """Outer radius [m]."""
        return self._layer_ptr.get().get_radius_outer()

    @property
    def mass(self) -> float:
        """Total layer mass [kg]."""
        return self._layer_ptr.get().get_mass()

    # ------------------------------------------------------------------------------------------------------------------
    # Geometry properties
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def name(self) -> str:
        """Layer name."""
        return self._layer_ptr.get().get_name().decode("utf-8")

    @property
    def layer_index(self) -> int:
        """Zero-based layer index (0 = innermost)."""
        return self._layer_ptr.get().get_layer_index()

    @property
    def radius_inner(self) -> float:
        """Inner boundary radius [m]."""
        return self._layer_ptr.get().get_radius_inner()

    @property
    def radius_outer(self) -> float:
        """Outer boundary radius [m]."""
        return self._layer_ptr.get().get_radius_outer()

    @property
    def thickness(self) -> float:
        """Layer thickness [m] (radius_outer - radius_inner)."""
        return self._layer_ptr.get().get_thickness()

    @property
    def volume(self) -> float:
        """Layer volume [m^3] (spherical shell)."""
        return self._layer_ptr.get().get_volume()

    @property
    def surface_area_inner(self) -> float:
        """Inner boundary surface area [m^2]."""
        return self._layer_ptr.get().get_surface_area_inner()

    @property
    def surface_area_outer(self) -> float:
        """Outer boundary surface area [m^2]."""
        return self._layer_ptr.get().get_surface_area_outer()

    @property
    def material_name(self) -> str:
        """Material identifier string."""
        return self._layer_ptr.get().get_material_name().decode("utf-8")

    @property
    def is_tidal(self) -> bool:
        """Whether this layer contributes to tidal dissipation."""
        return self._layer_ptr.get().get_is_tidal()

    @property
    def tidal_scale(self) -> float:
        """Dimensionless tidal heating scale factor (used when ``tidal_scale_method`` is ``user_provided``)."""
        return self._layer_ptr.get().get_tidal_scale()

    @property
    def tidal_scale_method(self) -> str:
        """How this layer's share of the world tidal heating is determined.

        One of ``"user_provided_scale"`` (use ``tidal_scale``), ``"volume_fraction_scale"``
        (layer volume / planet volume), or ``"tidal_timescale_scale"`` (Maxwell-time bell
        curve; not yet wired). Settable from any alias of those names.
        """
        cdef bytes name_bytes = c_tidal_scale_method_name(
            self._layer_ptr.get().get_tidal_scale_method())
        return name_bytes.decode("utf-8")

    @tidal_scale_method.setter
    def tidal_scale_method(self, str method):
        self._layer_ptr.get().set_tidal_scale_method(
            c_tidal_scale_method_from_name(method.encode("utf-8")))

    def get_tidal_heating(self) -> float:
        """Tidal heating [W] deposited in this layer by the world's last tidal solve.

        Set by :meth:`LayeredWorld.calc_tides` as the world's total tidal heating scaled
        by this layer's contribution. Returns NaN before a tidal solve has run.
        """
        return self._layer_ptr.get().get_tidal_heating()

    # ------------------------------------------------------------------------------------------------------------------
    # EOS profile
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def eos_data_populated(self) -> bool:
        """True after EOS profile data has been populated (world EOS solve or update_eos_data)."""
        return self._layer_ptr.get().get_eos_data_populated()

    @property
    def eos_set(self) -> bool:
        """True after a material EOS model has been attached via :meth:`set_eos`."""
        return self._layer_ptr.get().get_eos_set()

    def set_eos(self, MaterialEOSBase eos not None):
        """Attach a material EOS model (the per-layer density source).

        The model supplies density as a function of pressure and/or radius and is
        consumed by the world-level ``solve_eos`` when integrating the planet's
        radial structure. Ownership of the C++ model is transferred from ``eos``
        into this layer; the passed ``MaterialEOSBase`` becomes an empty,
        non-owning shell and must not be reused.

        Parameters
        ----------
        eos : MaterialEOSBase
            A material EOS model (e.g. ``make_material_eos("constant", ...)``,
            ``BirchMurnaghanEOS(...)``, ``InterpolatedEOS(...)``).

        Raises
        ------
        ValueError
            If ``eos`` has already been attached or otherwise moved.
        """
        if eos._eos_ptr.get() == NULL:
            raise ValueError(
                "This EOS model holds no C++ object (already attached or moved).")
        self._layer_ptr.get().set_eos(move(eos._eos_ptr))

    def update_eos_data(
            self,
            radius_m,
            density_kgm3,
            gravity_ms2,
            pressure_pa):
        """Populate the EOS profile from sorted radius arrays (MKS).

        In normal use this is populated by the world's EOS solve. For unit
        tests or manual construction the arrays can be provided directly.

        Parameters
        ----------
        radius_m : sequence of float
            Radius values [m], sorted ascending, matching the layer bounds.
        density_kgm3 : sequence of float
            Mass density [kg/m^3] at each radius.
        gravity_ms2 : sequence of float
            Gravitational acceleration [m/s^2] at each radius.
        pressure_pa : sequence of float
            Pressure [Pa] at each radius.

        Assumptions
        -----------
        - All four sequences must have the same length.
        - radius_m must be sorted in strictly ascending order.
        - Values are in MKS units.
        """
        cdef vector[double] r_vec   = radius_m
        cdef vector[double] rho_vec = density_kgm3
        cdef vector[double] g_vec   = gravity_ms2
        cdef vector[double] p_vec   = pressure_pa
        cdef c_LayerEOSData eos_data
        eos_data.populate(r_vec, rho_vec, g_vec, p_vec)
        self._layer_ptr.get().update_eos_data(eos_data)

    def get_density(self, double radius_m) -> float:
        """Density at the given radius [kg/m^3]; NaN if EOS data not populated.

        Parameters
        ----------
        radius_m : float
            Query radius [m].

        Assumptions
        -----------
        - Linear interpolation; clamped at layer boundaries.
        """
        return self._layer_ptr.get().get_density(radius_m)

    def get_gravity(self, double radius_m) -> float:
        """Gravitational acceleration at the given radius [m/s^2]; NaN if not populated.

        Parameters
        ----------
        radius_m : float
            Query radius [m].

        Assumptions
        -----------
        - Linear interpolation; clamped at layer boundaries.
        """
        return self._layer_ptr.get().get_gravity(radius_m)

    def get_pressure(self, double radius_m) -> float:
        """Pressure at the given radius [Pa]; NaN if EOS data not populated.

        Parameters
        ----------
        radius_m : float
            Query radius [m].

        Assumptions
        -----------
        - Linear interpolation; clamped at layer boundaries.
        """
        return self._layer_ptr.get().get_pressure(radius_m)

    # ------------------------------------------------------------------------------------------------------------------
    # Viscoelastic profile (populated by the world EOS solve; NaN before then or on a geometry-only layer)
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def viscoelastic_populated(self) -> bool:
        """True after the world EOS solve has populated this layer's viscoelastic state."""
        return self._layer_ptr.get().get_viscoelastic_populated()

    def get_shear_modulus(self, double radius_m) -> float:
        """Post-melt static shear modulus [Pa] at radius_m [m]; NaN if unpopulated."""
        return self._layer_ptr.get().get_shear_modulus(radius_m)

    def get_bulk_modulus(self, double radius_m) -> float:
        """Post-melt static bulk modulus [Pa] at radius_m [m]; NaN if unpopulated."""
        return self._layer_ptr.get().get_bulk_modulus(radius_m)

    def get_shear_viscosity(self, double radius_m) -> float:
        """Post-melt shear viscosity [Pa s] at radius_m [m]; NaN if unpopulated."""
        return self._layer_ptr.get().get_shear_viscosity(radius_m)

    def get_bulk_viscosity(self, double radius_m) -> float:
        """Post-melt bulk viscosity [Pa s] at radius_m [m]; NaN if unpopulated."""
        return self._layer_ptr.get().get_bulk_viscosity(radius_m)

    def get_premelt_shear_modulus(self, double radius_m) -> float:
        """Pre-melt static shear modulus [Pa] at radius_m [m]; NaN if unpopulated."""
        return self._layer_ptr.get().get_premelt_shear_modulus(radius_m)

    def get_premelt_bulk_modulus(self, double radius_m) -> float:
        """Pre-melt static bulk modulus [Pa] at radius_m [m]; NaN if unpopulated."""
        return self._layer_ptr.get().get_premelt_bulk_modulus(radius_m)

    def get_premelt_shear_viscosity(self, double radius_m) -> float:
        """Pre-melt shear viscosity [Pa s] at radius_m [m]; NaN if unpopulated."""
        return self._layer_ptr.get().get_premelt_shear_viscosity(radius_m)

    def get_premelt_bulk_viscosity(self, double radius_m) -> float:
        """Pre-melt bulk viscosity [Pa s] at radius_m [m]; NaN if unpopulated."""
        return self._layer_ptr.get().get_premelt_bulk_viscosity(radius_m)

    # ------------------------------------------------------------------------------------------------------------------
    # Config
    # ------------------------------------------------------------------------------------------------------------------
    cpdef dict get_config_dict(self):
        """Return all configuration values as a Python dict (MKS).

        Returns
        -------
        dict
            Keys: ``name``, ``layer_index``, ``radius_inner_m``,
            ``radius_outer_m``, ``mass_kg``, ``material_name``,
            ``is_tidal``, ``tidal_scale``, ``tidal_scale_method``.
        """
        cdef c_BaseLayer* p = self._layer_ptr.get()
        cdef bytes method_bytes = c_tidal_scale_method_name(p.get_tidal_scale_method())
        return {
            "name":               p.get_name().decode("utf-8"),
            "layer_index":        p.get_layer_index(),
            "radius_inner_m":     p.get_radius_inner(),
            "radius_outer_m":     p.get_radius_outer(),
            "mass_kg":            p.get_mass(),
            "material_name":      p.get_material_name().decode("utf-8"),
            "is_tidal":           bool(p.get_is_tidal()),
            "tidal_scale":        p.get_tidal_scale(),
            "tidal_scale_method": method_bytes.decode("utf-8"),
        }
