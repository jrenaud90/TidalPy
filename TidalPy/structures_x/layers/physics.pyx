# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
physics.pyx
Cython/Python wrapper for TidalPy's physics layer class.

PhysicsLayer: extends BaseLayer with static mechanical properties (shear
modulus, bulk modulus, shear and bulk viscosity) and three complex Love
numbers (k, h, l).  Optional rheology objects (Phase 5) enable
frequency-dependent complex moduli; until then the static modulus is
returned as a real-valued complex number.
"""

from libcpp.complex cimport complex as cpp_complex
from libcpp cimport bool as cpp_bool
from libcpp.utility cimport move

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport c_TidalPyBaseClass
from TidalPy.structures_x.layers.base cimport BaseLayer, c_BaseLayer, c_tidal_scale_method_from_name
from TidalPy.Tides_x.love.love cimport LoveNumbers, c_LoveNumbers
from TidalPy.rheology_x.rheology cimport RheologyBase
from TidalPy.viscosity_x.viscosity cimport ViscosityBase
from TidalPy.partial_melt_x.partial_melt cimport PartialMeltBase

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


# =====================================================================================================================
# PhysicsLayer
# =====================================================================================================================

cdef class PhysicsLayer(BaseLayer):
    """Mechanical-properties layer: static shear/bulk modulus, viscosities, Love numbers, and optional rheology.

    Extends BaseLayer with material-mechanical parameters needed for tidal
    calculations.  Rheology objects (Phase 5) can be attached to enable
    frequency-dependent complex moduli; until then the static modulus is
    returned as a real-valued complex number.

    Parameters
    ----------
    name : str
        Human-readable layer name.
    layer_index : int
        Zero-based position; innermost layer = 0.
    radius_inner_m : float
        Inner boundary radius [m].
    radius_outer_m : float
        Outer boundary radius [m].
    mass_kg : float
        Total layer mass [kg].
    material_name : str, optional
        Material identifier (e.g. ``"perovskite"``). Default ``""``.
    is_tidal : bool, optional
        Whether this layer contributes to tidal dissipation. Default ``True``.
    tidal_scale : float, optional
        Dimensionless tidal heating scale. Default ``1.0``.
    shear_modulus_static_pa : float, optional
        Unrelaxed shear modulus [Pa]. Default ``0.0``.
    bulk_modulus_static_pa : float, optional
        Unrelaxed bulk modulus [Pa]. Default ``0.0``.
    shear_viscosity_static_pas : float, optional
        Reference dynamic shear viscosity [Pa·s]. Default ``0.0``.
    bulk_viscosity_static_pas : float, optional
        Reference dynamic bulk viscosity [Pa·s]. Default ``0.0``.
    love_number_k : complex, optional
        Potential Love number (placeholder). Default ``0+0j``.
    love_number_h : complex, optional
        Radial displacement Love number (placeholder). Default ``0+0j``.
    love_number_l : complex, optional
        Tangential displacement Love number (placeholder). Default ``0+0j``.

    Assumptions
    -----------
    - Layer geometry is spherically symmetric.
    - radius_inner_m <= radius_outer_m.
    - All values are in MKS units.
    """

    def __cinit__(self, *args, **kwargs):
        self._physics_ptr = NULL

    def __init__(
            self,
            str    name,
            int    layer_index,
            double radius_inner_m,
            double radius_outer_m,
            double mass_kg,
            str    material_name                = "",
            bint   is_tidal                     = True,
            double tidal_scale                  = 1.0,
            double shear_modulus_static_pa      = 0.0,
            double bulk_modulus_static_pa       = 0.0,
            double shear_viscosity_static_pas   = 0.0,
            double bulk_viscosity_static_pas    = 0.0,
            complex love_number_k               = 0+0j,
            complex love_number_h               = 0+0j,
            complex love_number_l               = 0+0j,
            str    tidal_scale_method           = "user_provided"):
        cdef c_PhysicsConfig config
        config.name                       = name.encode("utf-8")
        config.layer_index                = layer_index
        config.radius_inner_m             = radius_inner_m
        config.radius_outer_m             = radius_outer_m
        config.mass_kg                    = mass_kg
        config.material_name              = material_name.encode("utf-8")
        config.is_tidal                   = is_tidal
        config.tidal_scale                = tidal_scale
        config.tidal_scale_method         = c_tidal_scale_method_from_name(tidal_scale_method.encode("utf-8"))
        config.shear_modulus_static_pa    = shear_modulus_static_pa
        config.bulk_modulus_static_pa     = bulk_modulus_static_pa
        config.shear_viscosity_static_pas = shear_viscosity_static_pas
        config.bulk_viscosity_static_pas  = bulk_viscosity_static_pas
        config.love_numbers = c_LoveNumbers(
            cpp_complex[double](love_number_k.real, love_number_k.imag),
            cpp_complex[double](love_number_h.real, love_number_h.imag),
            cpp_complex[double](love_number_l.real, love_number_l.imag))
        cdef c_PhysicsLayer* raw = new c_PhysicsLayer(config)
        self._layer_ptr.reset(<c_BaseLayer*>raw)
        self._physics_ptr = raw
        self._ptr         = <c_TidalPyBaseClass*>raw

    def __dealloc__(self):
        self._physics_ptr = NULL  # base's unique_ptr owns the C++ object

    @staticmethod
    cdef PhysicsLayer _view(c_PhysicsLayer* ptr, object world):
        cdef PhysicsLayer v = PhysicsLayer.__new__(PhysicsLayer)
        v._physics_ptr = ptr
        v._init_view(<c_BaseLayer*>ptr, world)
        return v

    # ------------------------------------------------------------------------------------------------------------------
    # Static mechanical property properties
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def shear_modulus_static(self) -> float:
        """Unrelaxed (static) shear modulus [Pa]."""
        return self._physics_ptr.get_shear_modulus_static()

    @property
    def bulk_modulus_static(self) -> float:
        """Unrelaxed (static) bulk modulus [Pa]."""
        return self._physics_ptr.get_bulk_modulus_static()

    @property
    def shear_viscosity_static(self) -> float:
        """Reference dynamic shear viscosity [Pa·s]."""
        return self._physics_ptr.get_shear_viscosity_static()

    @property
    def bulk_viscosity_static(self) -> float:
        """Reference dynamic bulk viscosity [Pa·s]."""
        return self._physics_ptr.get_bulk_viscosity_static()

    @property
    def love_numbers(self) -> LoveNumbers:
        """All three complex Love numbers (k, h, l) as a LoveNumbers object."""
        cdef LoveNumbers result = LoveNumbers.__new__(LoveNumbers)
        result._love = self._physics_ptr.get_love_numbers()
        return result

    @property
    def love_number_k(self) -> complex:
        """Complex potential Love number k (placeholder; computed in later phases)."""
        cdef cpp_complex[double] k = self._physics_ptr.get_love_number_k()
        return complex(k.real(), k.imag())

    @property
    def love_number_h(self) -> complex:
        """Complex radial displacement Love number h (placeholder; computed in later phases)."""
        cdef cpp_complex[double] h = self._physics_ptr.get_love_number_h()
        return complex(h.real(), h.imag())

    @property
    def love_number_l(self) -> complex:
        """Complex tangential displacement Love number l (placeholder; computed in later phases)."""
        cdef cpp_complex[double] l = self._physics_ptr.get_love_number_l()
        return complex(l.real(), l.imag())

    @property
    def shear_rheology_set(self) -> bool:
        """True if a shear rheology model has been attached."""
        return self._physics_ptr.get_shear_rheology_set()

    @property
    def bulk_rheology_set(self) -> bool:
        """True if a bulk rheology model has been attached."""
        return self._physics_ptr.get_bulk_rheology_set()

    # ------------------------------------------------------------------------------------------------------------------
    # Radial-solver layer classification flags
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def is_solid(self) -> bool:
        """True if this layer is solid; False for liquid.  Used by the radial Love-number solver."""
        return bool(self._physics_ptr.get_is_solid())

    @is_solid.setter
    def is_solid(self, value: bool):
        self._physics_ptr.set_is_solid(<cpp_bool>bool(value))

    @property
    def is_static(self) -> bool:
        """True if the static (no-dynamic-terms) approximation is used.  Used by the radial solver."""
        return bool(self._physics_ptr.get_is_static())

    @is_static.setter
    def is_static(self, value: bool):
        self._physics_ptr.set_is_static(<cpp_bool>bool(value))

    @property
    def is_incompressible(self) -> bool:
        """True if the incompressible approximation is used.  Used by the radial solver."""
        return bool(self._physics_ptr.get_is_incompressible())

    @is_incompressible.setter
    def is_incompressible(self, value: bool):
        self._physics_ptr.set_is_incompressible(<cpp_bool>bool(value))

    # ------------------------------------------------------------------------------------------------------------------
    # Rheology attachment
    # ------------------------------------------------------------------------------------------------------------------
    def set_shear_rheology(self, RheologyBase rheology not None):
        """Attach a rheology model used to compute the complex shear modulus.

        Ownership of the C++ model is transferred from ``rheology`` into this
        layer; the passed ``RheologyBase`` becomes an empty, non-owning shell and
        must not be reused.

        Parameters
        ----------
        rheology : RheologyBase
            A rheology model (e.g. ``Maxwell()``, ``make_rheology("andrade")``).

        Raises
        ------
        ValueError
            If ``rheology`` has already been attached or otherwise moved.
        """
        if rheology._rheology_ptr.get() == NULL:
            raise ValueError(
                "This rheology model holds no C++ object (already attached or moved).")
        self._physics_ptr.set_shear_rheology(move(rheology._rheology_ptr))

    def set_bulk_rheology(self, RheologyBase rheology not None):
        """Attach a rheology model used to compute the complex bulk modulus.

        Ownership of the C++ model is transferred from ``rheology`` into this
        layer; the passed ``RheologyBase`` becomes an empty, non-owning shell and
        must not be reused.

        Parameters
        ----------
        rheology : RheologyBase
            A rheology model (e.g. ``Maxwell()``, ``make_rheology("andrade")``).

        Raises
        ------
        ValueError
            If ``rheology`` has already been attached or otherwise moved.
        """
        if rheology._rheology_ptr.get() == NULL:
            raise ValueError(
                "This rheology model holds no C++ object (already attached or moved).")
        self._physics_ptr.set_bulk_rheology(move(rheology._rheology_ptr))

    # ------------------------------------------------------------------------------------------------------------------
    # Viscosity + partial-melt attachment
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def shear_viscosity_set(self) -> bool:
        """True if a shear viscosity model has been attached."""
        return self._physics_ptr.get_shear_viscosity_set()

    @property
    def bulk_viscosity_set(self) -> bool:
        """True if a bulk viscosity model has been attached."""
        return self._physics_ptr.get_bulk_viscosity_set()

    @property
    def partial_melt_set(self) -> bool:
        """True if a partial-melt model has been attached."""
        return self._physics_ptr.get_partial_melt_set()

    def set_shear_viscosity(self, ViscosityBase viscosity not None):
        """Attach a viscosity model supplying the pre-melt shear viscosity.

        Ownership of the C++ model transfers into this layer; the passed
        ``ViscosityBase`` becomes an empty, non-owning shell and must not be reused.

        Raises
        ------
        ValueError
            If ``viscosity`` has already been attached or otherwise moved.
        """
        if viscosity._visc_ptr.get() == NULL:
            raise ValueError(
                "This viscosity model holds no C++ object (already attached or moved).")
        self._physics_ptr.set_shear_viscosity(move(viscosity._visc_ptr))

    def set_bulk_viscosity(self, ViscosityBase viscosity not None):
        """Attach a viscosity model supplying the pre-melt bulk viscosity.

        Ownership of the C++ model transfers into this layer; the passed
        ``ViscosityBase`` becomes an empty, non-owning shell and must not be reused.

        Raises
        ------
        ValueError
            If ``viscosity`` has already been attached or otherwise moved.
        """
        if viscosity._visc_ptr.get() == NULL:
            raise ValueError(
                "This viscosity model holds no C++ object (already attached or moved).")
        self._physics_ptr.set_bulk_viscosity(move(viscosity._visc_ptr))

    def set_partial_melt(self, PartialMeltBase partial_melt not None):
        """Attach a partial-melt model that weakens the static moduli and viscosities.

        Ownership of the C++ model transfers into this layer; the passed
        ``PartialMeltBase`` becomes an empty, non-owning shell and must not be reused.

        Raises
        ------
        ValueError
            If ``partial_melt`` has already been attached or otherwise moved.
        """
        if partial_melt._melt_ptr.get() == NULL:
            raise ValueError(
                "This partial-melt model holds no C++ object (already attached or moved).")
        self._physics_ptr.set_partial_melt(move(partial_melt._melt_ptr))

    # ------------------------------------------------------------------------------------------------------------------
    # Calculations
    # ------------------------------------------------------------------------------------------------------------------
    def calc_tidal_susceptibility(self) -> float:
        """Tidal susceptibility (3/2) * r^5 / (G * m^2) [m^3].

        A purely geometrical quantity; Newton's G is taken from the TidalPy
        global config pointer.

        Returns
        -------
        float
            Tidal susceptibility [m^3].  Returns 0.0 when config is not
            initialised or mass is zero.

        Assumptions
        -----------
        - Spherically symmetric layer geometry.
        - G sourced from TidalPy global configuration (set at package init).
        """
        return self._physics_ptr.calc_tidal_susceptibility()

    def calc_complex_shear_modulus(self, double frequency_rad_s) -> complex:
        """Complex shear modulus at the given tidal forcing frequency [Pa].

        When a shear rheology model is attached (Phase 5), the result is
        computed from the model's complex compliance.  Until then the static
        shear modulus is returned as a real-valued complex number.

        Parameters
        ----------
        frequency_rad_s : float
            Tidal forcing frequency [rad/s].

        Returns
        -------
        complex
            Complex shear modulus [Pa].

        Assumptions
        -----------
        - Linear viscoelastic response (single forcing frequency).
        """
        cdef cpp_complex[double] result = self._physics_ptr.calc_complex_shear_modulus(frequency_rad_s)
        return complex(result.real(), result.imag())

    def calc_complex_bulk_modulus(self, double frequency_rad_s) -> complex:
        """Complex bulk modulus at the given tidal forcing frequency [Pa].

        When a bulk rheology model is attached (Phase 5), the result is
        computed from the model's complex compliance.  Until then the static
        bulk modulus is returned as a real-valued complex number.

        Parameters
        ----------
        frequency_rad_s : float
            Tidal forcing frequency [rad/s].

        Returns
        -------
        complex
            Complex bulk modulus [Pa].

        Assumptions
        -----------
        - Linear viscoelastic response (single forcing frequency).
        """
        cdef cpp_complex[double] result = self._physics_ptr.calc_complex_bulk_modulus(frequency_rad_s)
        return complex(result.real(), result.imag())

    # ------------------------------------------------------------------------------------------------------------------
    # Config
    # ------------------------------------------------------------------------------------------------------------------
    cpdef dict get_config_dict(self):
        """Return all configuration values as a Python dict (MKS).

        Returns
        -------
        dict
            Keys: all BaseLayer keys plus ``shear_modulus_static_pa``,
            ``bulk_modulus_static_pa``, ``shear_viscosity_static_pas``,
            ``bulk_viscosity_static_pas``, and Love number components
            ``love_number_k_re``, ``love_number_k_im``, ``love_number_h_re``,
            ``love_number_h_im``, ``love_number_l_re``, ``love_number_l_im``.
        """
        d = BaseLayer.get_config_dict(self)
        d["shear_modulus_static_pa"]      = self._physics_ptr.get_shear_modulus_static()
        d["bulk_modulus_static_pa"]       = self._physics_ptr.get_bulk_modulus_static()
        d["shear_viscosity_static_pas"]   = self._physics_ptr.get_shear_viscosity_static()
        d["bulk_viscosity_static_pas"]    = self._physics_ptr.get_bulk_viscosity_static()
        cdef c_LoveNumbers ln = self._physics_ptr.get_love_numbers()
        d["love_number_k_re"] = ln.k.real()
        d["love_number_k_im"] = ln.k.imag()
        d["love_number_h_re"] = ln.h.real()
        d["love_number_h_im"] = ln.h.imag()
        d["love_number_l_re"] = ln.l.real()
        d["love_number_l_im"] = ln.l.imag()
        return d
