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
        cdef c_BaseLayer* lp
        layers = []
        for i in range(n):
            lp = self._layered_ptr.get_layer(i)
            layers.append({
                "name":           lp.get_name().decode("utf-8"),
                "layer_index":    lp.get_layer_index(),
                "radius_inner_m": lp.get_radius_inner(),
                "radius_outer_m": lp.get_radius_outer(),
                "mass_kg":        lp.get_mass(),
                "material_name":  lp.get_material_name().decode("utf-8"),
                "is_tidal":       True if lp.get_is_tidal() else False,
                "tidal_scale":    lp.get_tidal_scale(),
            })
        d["num_layers"] = n
        d["layers"] = layers
        return d
