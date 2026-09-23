# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrapper for TidalPy's gas layer class.

GasLayer extends PhysicsLayer with ideal-gas parameters (mean molecular weight, adiabatic index, and a reference
state), stored and serialized for a future gas description; nothing reads them yet, and the layer's density comes
from its material's law. It has no phase changes and no cooling or radiogenics sub-models.
"""

from libcpp.complex cimport complex as cpp_complex
from libcpp cimport bool as cpp_bool
from libcpp.memory cimport unique_ptr, make_unique

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport d_NAN, set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport c_TidalPyBaseClass
from TidalPy.structures_x.layers.base cimport BaseLayer, c_BaseLayer
from TidalPy.structures_x.layers.physics cimport PhysicsLayer, c_PhysicsLayer
from TidalPy.Tides_x.love.love cimport LoveNumbers, c_LoveNumbers

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


cdef class GasLayer(PhysicsLayer):
    """Gas layer: PhysicsLayer plus stored ideal-gas parameters, which nothing reads yet.

    No phase changes, cooling, or radiogenics sub-models are available (use SolidLiquidLayer for those).

    Parameters
    ----------
    name : str
        Layer name.
    layer_index : int
        Zero-based position; innermost layer = 0.
    radius_inner : float
        Inner boundary radius [m].
    radius_outer : float
        Outer boundary radius [m].
    mass : float
        Total layer mass [kg].
    material_name : str, optional
        Material identifier. Default ``""``.
    is_volume_fixed : bool, optional
        False lets the layer grow or shrink to hold its mass during an EOS solve. Default ``True``.
    is_tidal : bool, optional
        Whether this layer contributes to tidal dissipation. Default ``True``.
    tidal_scale : float, optional
        The layer's share of the planet in the quasi-homogeneous Love methods; ``None`` (default) takes
        its volume fraction. See ``BaseLayer``.
    love_number_k : complex, optional
        Potential Love number k (placeholder). Default ``0+0j``.
    love_number_h : complex, optional
        Radial displacement Love number h (placeholder). Default ``0+0j``.
    love_number_l : complex, optional
        Tangential displacement Love number l (placeholder). Default ``0+0j``.
    mean_molecular_weight : float, optional
        Mean molecular weight of the gas [kg/mol]. Default ``2e-3`` (H₂).
    adiabatic_index : float, optional
        Ratio of specific heats γ = c_p/c_v [dimensionless]. Default ``1.4``.
    reference_temperature : float, optional
        Reference temperature [K]. Default ``300.0``.
    reference_density : float, optional
        Reference density [kg/m³]. Default ``1.0``.
    is_solid : bool, optional
        True marks the layer solid for the radial Love-number solver. Default ``False``: a gas carries no
        shear stress, so it is solved as a liquid.
    is_static : bool, optional
        Use the static (no inertia) approximation in the radial solver. Default ``True``.
    is_incompressible : bool, optional
        Use the incompressible approximation in the radial solver. Default ``False``.
    temperature : float, optional
        Layer temperature [K] at which its viscosity and melt models are evaluated. Default ``0.0``, the cold
        rigid limit of the viscosity laws.
    use_thermal_eos : bool, optional
        Pass the temperature to the EOS model, so the density and bulk modulus depend on it. Default ``False``.
    use_heating : bool, optional
        Let the world's heat sources (this layer's radiogenics model among them) act inside the layer during a
        thermal EOS solve. Default ``False``.

    Assumptions
    -----------
    - Spherically symmetric layer geometry.
    - The density, moduli, and viscosities are the material's, as for any physics layer; the ideal-gas parameters
      take no part in any calculation yet.
    """

    def __cinit__(self, *args, **kwargs):
        self._gas_ptr = NULL

    def __init__(
            self,
            str    name,
            int    layer_index,
            double radius_inner,
            double radius_outer,
            double mass,
            str    material_name          = "",
            cpp_bool is_tidal             = True,
            cpp_bool is_volume_fixed      = True,
            tidal_scale                   = None,
            complex love_number_k         = 0+0j,
            complex love_number_h         = 0+0j,
            complex love_number_l         = 0+0j,
            double mean_molecular_weight  = 2.0e-3,
            double adiabatic_index        = 1.4,
            double reference_temperature  = 300.0,
            double reference_density      = 1.0,
            cpp_bool is_solid             = False,
            cpp_bool is_static            = True,
            cpp_bool is_incompressible    = False,
            double temperature            = 0.0,
            cpp_bool use_thermal_eos = False,
            cpp_bool use_heating     = False):
        cdef c_GasConfig config
        config.name                 = name.encode("utf-8")
        config.layer_index          = layer_index
        config.radius_inner         = radius_inner
        config.radius_outer         = radius_outer
        config.mass                 = mass
        config.material_name        = material_name.encode("utf-8")
        config.is_tidal             = is_tidal
        config.is_volume_fixed      = is_volume_fixed
        config.tidal_scale          = d_NAN if tidal_scale is None else <double>tidal_scale
        config.love_numbers = c_LoveNumbers(
            cpp_complex[double](love_number_k.real, love_number_k.imag),
            cpp_complex[double](love_number_h.real, love_number_h.imag),
            cpp_complex[double](love_number_l.real, love_number_l.imag))
        config.is_solid              = is_solid
        config.is_static             = is_static
        config.is_incompressible     = is_incompressible
        config.temperature       = temperature
        config.use_thermal_eos   = use_thermal_eos
        config.use_heating       = use_heating
        config.mean_molecular_weight = mean_molecular_weight
        config.adiabatic_index       = adiabatic_index
        config.reference_temperature = reference_temperature
        config.reference_density     = reference_density
        # make_unique owns the allocation; ownership then moves into the base-typed member
        # (Cython cannot assign a unique_ptr[Derived] to a unique_ptr[Base] directly).
        cdef unique_ptr[c_GasLayer] built = make_unique[c_GasLayer](config)
        self._gas_ptr     = built.get()
        self._physics_ptr = <c_PhysicsLayer*>self._gas_ptr
        self._layer_ptr.reset(<c_BaseLayer*>built.release())
        self._ptr = <c_TidalPyBaseClass*>self._layer_ptr.get()

    def __dealloc__(self):
        self._gas_ptr     = NULL  # base's unique_ptr owns the C++ object
        self._physics_ptr = NULL

    @staticmethod
    cdef GasLayer _view(c_GasLayer* ptr, object world):
        cdef GasLayer v = GasLayer.__new__(GasLayer)
        v._gas_ptr      = ptr
        v._physics_ptr  = <c_PhysicsLayer*>ptr
        v._init_view(<c_BaseLayer*>ptr, world)
        return v

    @property
    def mean_molecular_weight(self) -> float:
        """Mean molecular weight of the gas [kg/mol]."""
        return self._gas_ptr.get_mean_molecular_weight()

    @property
    def adiabatic_index(self) -> float:
        """Ratio of specific heats γ = c_p/c_v [dimensionless]."""
        return self._gas_ptr.get_adiabatic_index()

    @property
    def reference_temperature(self) -> float:
        """Reference temperature [K]."""
        return self._gas_ptr.get_reference_temperature()

    @property
    def reference_density(self) -> float:
        """Reference density [kg/m³]."""
        return self._gas_ptr.get_reference_density()

    cpdef dict get_config_dict(self):
        """Return all configuration values as a Python dict (MKS): the PhysicsLayer keys plus the gas parameters.

        The layer's own ``reference_density`` is left out: a layer file cannot carry it (the layer's density is its
        material's, in the ``material`` table).
        """
        cdef dict d = PhysicsLayer.get_config_dict(self)
        d["mean_molecular_weight_kg_mol"] = self._gas_ptr.get_mean_molecular_weight()
        d["adiabatic_index"]              = self._gas_ptr.get_adiabatic_index()
        d["reference_temperature_k"]      = self._gas_ptr.get_reference_temperature()
        return d
