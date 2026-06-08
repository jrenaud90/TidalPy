# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
base.pyx
Cython/Python wrapper for TidalPy's base layer class.

BaseLayer: geometry-only layer storing inner/outer radii, mass, and
material identification. EOS profile (density, gravity, pressure vs. radius)
is unpopulated until EOSHandler runs (Phase 8) or update_eos_data() is called
directly (e.g., for testing).
"""

import os as _os

from libcpp.vector cimport vector
from libcpp cimport bool as cpp_bool

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport StructureBase

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
    After EOSHandler solves the Adams-Williamson ODE (Phase 8), or after a
    direct call to ``update_eos_data``, the profile is interpolatable.

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
        # StructureBase.__cinit__ already ran and set _ptr = &_struct.
        # We override _ptr in __init__ once _layer is constructed.
        pass

    def __init__(
            self,
            str    name,
            int    layer_index,
            double radius_inner_m,
            double radius_outer_m,
            double mass_kg,
            str    material_name = "",
            bint   is_tidal      = True,
            double tidal_scale   = 1.0):
        cdef c_BaseLayerConfig config
        config.name           = name.encode("utf-8")
        config.layer_index    = layer_index
        config.radius_inner_m = radius_inner_m
        config.radius_outer_m = radius_outer_m
        config.mass_kg        = mass_kg
        config.material_name  = material_name.encode("utf-8")
        config.is_tidal       = is_tidal
        config.tidal_scale    = tidal_scale
        self._layer = c_BaseLayer(config)
        self._ptr   = &self._layer

    def __dealloc__(self):
        self._ptr = NULL

    # ------------------------------------------------------------------------------------------------------------------
    # Base class property overrides
    # (StructureBase._struct is default-constructed only; use _layer for stored state)
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def radius(self) -> float:
        """Outer radius [m]."""
        return self._layer.get_radius_outer()

    @property
    def mass(self) -> float:
        """Total layer mass [kg]."""
        return self._layer.get_mass()

    # ------------------------------------------------------------------------------------------------------------------
    # Geometry properties
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def name(self) -> str:
        """Layer name."""
        return self._layer.get_name().decode("utf-8")

    @property
    def layer_index(self) -> int:
        """Zero-based layer index (0 = innermost)."""
        return self._layer.get_layer_index()

    @property
    def radius_inner(self) -> float:
        """Inner boundary radius [m]."""
        return self._layer.get_radius_inner()

    @property
    def radius_outer(self) -> float:
        """Outer boundary radius [m]."""
        return self._layer.get_radius_outer()

    @property
    def thickness(self) -> float:
        """Layer thickness [m] (radius_outer - radius_inner)."""
        return self._layer.get_thickness()

    @property
    def volume(self) -> float:
        """Layer volume [m^3] (spherical shell)."""
        return self._layer.get_volume()

    @property
    def surface_area_inner(self) -> float:
        """Inner boundary surface area [m^2]."""
        return self._layer.get_surface_area_inner()

    @property
    def surface_area_outer(self) -> float:
        """Outer boundary surface area [m^2]."""
        return self._layer.get_surface_area_outer()

    @property
    def material_name(self) -> str:
        """Material identifier string."""
        return self._layer.get_material_name().decode("utf-8")

    @property
    def is_tidal(self) -> bool:
        """Whether this layer contributes to tidal dissipation."""
        return self._layer.get_is_tidal()

    @property
    def tidal_scale(self) -> float:
        """Dimensionless tidal heating scale factor."""
        return self._layer.get_tidal_scale()

    # ------------------------------------------------------------------------------------------------------------------
    # EOS profile
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def eos_data_populated(self) -> bool:
        """True after EOS profile data has been populated (EOSHandler or update_eos_data)."""
        return self._layer.get_eos_data_populated()

    def update_eos_data(
            self,
            radius_m,
            density_kgm3,
            gravity_ms2,
            pressure_pa):
        """Populate the EOS profile from sorted radius arrays (MKS).

        In normal use this is called by EOSHandler (Phase 8). For unit tests
        or manual construction the arrays can be provided directly.

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
        self._layer.update_eos_data(eos_data)

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
        return self._layer.get_density(radius_m)

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
        return self._layer.get_gravity(radius_m)

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
        return self._layer.get_pressure(radius_m)

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
            ``is_tidal``, ``tidal_scale``.
        """
        return {
            "name":           self._layer.get_name().decode("utf-8"),
            "layer_index":    self._layer.get_layer_index(),
            "radius_inner_m": self._layer.get_radius_inner(),
            "radius_outer_m": self._layer.get_radius_outer(),
            "mass_kg":        self._layer.get_mass(),
            "material_name":  self._layer.get_material_name().decode("utf-8"),
            "is_tidal":       bool(self._layer.get_is_tidal()),
            "tidal_scale":    self._layer.get_tidal_scale(),
        }
