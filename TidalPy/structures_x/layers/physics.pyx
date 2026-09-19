# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrapper for TidalPy's physics layer class.

PhysicsLayer extends BaseLayer with static mechanical properties (shear modulus, bulk modulus, shear and bulk
viscosity) and the three complex Love numbers. Attaching a rheology gives frequency-dependent complex moduli;
without one the static modulus is returned as a real-valued complex number.
"""

cimport numpy as cnp
cnp.import_array()

import numpy as np

from libcpp.complex cimport complex as cpp_complex
from libcpp cimport bool as cpp_bool
from libcpp.utility cimport move
from libcpp.memory cimport make_unique

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport d_NAN, set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport c_TidalPyBaseClass, c_PhysicsBase, cy_physics_model_config
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

    Rheology models can be attached to give frequency-dependent complex moduli; without one the static modulus
    is returned as a real-valued complex number.

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
        Material identifier (e.g. ``"perovskite"``). Default ``""``.
    is_tidal : bool, optional
        Whether this layer contributes to tidal dissipation. Default ``True``.
    tidal_scale : float, optional
        Dimensionless tidal heating scale. Default ``1.0``.
    shear_modulus_static : float, optional
        Unrelaxed shear modulus [Pa]. Default ``0.0``.
    bulk_modulus_static : float, optional
        Unrelaxed bulk modulus [Pa]. Default ``0.0``.
    shear_viscosity_static : float, optional
        Reference dynamic shear viscosity [Pa·s]. Default NaN (unset).
    bulk_viscosity_static : float, optional
        Reference dynamic bulk viscosity [Pa·s]. Default NaN (unset).
    love_number_k : complex, optional
        Potential Love number (placeholder). Default ``0+0j``.
    love_number_h : complex, optional
        Radial displacement Love number (placeholder). Default ``0+0j``.
    love_number_l : complex, optional
        Tangential displacement Love number (placeholder). Default ``0+0j``.
    tidal_scale_method : str, optional
        How the layer's share of the world's tidal heating is set. Default ``"user_provided"``.
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
    shear_modulus_pressure_derivative : float, optional
        Pressure derivative of the static shear modulus [Pa/Pa]. Default ``0.0``.
    shear_modulus_temperature_derivative : float, optional
        Temperature derivative of the static shear modulus [Pa/K]. Default ``0.0``.
    shear_modulus_reference_temperature : float, optional
        Temperature [K] at which ``shear_modulus_static`` applies. ``None`` keeps the default of 300 K.
    use_thermal_eos : bool, optional
        Pass the temperature to the EOS model, so the density and bulk modulus depend on it. Default ``False``.

    Assumptions
    -----------
    - Layer geometry is spherically symmetric.
    - radius_inner <= radius_outer.
    """

    def __cinit__(self, *args, **kwargs):
        self._physics_ptr = NULL

    def __init__(
            self,
            str    name,
            int    layer_index,
            double radius_inner,
            double radius_outer,
            double mass,
            str    material_name        = "",
            cpp_bool is_tidal           = True,
            double tidal_scale          = 1.0,
            double shear_modulus_static = 0.0,
            double bulk_modulus_static  = 0.0,
            double shear_viscosity_static = d_NAN,
            double bulk_viscosity_static = d_NAN,
            complex love_number_k        = 0+0j,
            complex love_number_h        = 0+0j,
            complex love_number_l        = 0+0j,
            str    tidal_scale_method    = "user_provided",
            cpp_bool is_solid            = True,
            cpp_bool is_static           = True,
            cpp_bool is_incompressible   = False,
            double temperature           = 0.0,
            double shear_modulus_pressure_derivative    = 0.0,
            double shear_modulus_temperature_derivative = 0.0,
            shear_modulus_reference_temperature         = None,
            cpp_bool use_thermal_eos     = False):
        cdef c_PhysicsConfig config
        config.name               = name.encode("utf-8")
        config.layer_index        = layer_index
        config.radius_inner       = radius_inner
        config.radius_outer       = radius_outer
        config.mass               = mass
        config.material_name      = material_name.encode("utf-8")
        config.is_tidal           = is_tidal
        config.tidal_scale        = tidal_scale
        config.tidal_scale_method = c_tidal_scale_method_from_name(tidal_scale_method.encode("utf-8"))
        config.shear_modulus_static = shear_modulus_static
        config.bulk_modulus_static  = bulk_modulus_static
        config.shear_viscosity_static = shear_viscosity_static
        config.bulk_viscosity_static  = bulk_viscosity_static
        config.love_numbers = c_LoveNumbers(
            cpp_complex[double](love_number_k.real, love_number_k.imag),
            cpp_complex[double](love_number_h.real, love_number_h.imag),
            cpp_complex[double](love_number_l.real, love_number_l.imag))
        config.is_solid          = is_solid
        config.is_static         = is_static
        config.is_incompressible = is_incompressible
        config.temperature       = temperature
        config.shear_modulus_pressure_derivative    = shear_modulus_pressure_derivative
        config.shear_modulus_temperature_derivative = shear_modulus_temperature_derivative
        # None keeps the C++ default reference temperature.
        if shear_modulus_reference_temperature is not None:
            config.shear_modulus_reference_temperature = <double>shear_modulus_reference_temperature
        config.use_thermal_eos   = use_thermal_eos
        # make_unique owns the allocation; ownership then moves into the base-typed member
        # (Cython cannot assign a unique_ptr[Derived] to a unique_ptr[Base] directly).
        cdef unique_ptr[c_PhysicsLayer] built = make_unique[c_PhysicsLayer](config)
        self._physics_ptr = built.get()
        self._layer_ptr.reset(<c_BaseLayer*>built.release())
        self._ptr = <c_TidalPyBaseClass*>self._layer_ptr.get()

    def __dealloc__(self):
        self._physics_ptr = NULL  # base's unique_ptr owns the C++ object

    @staticmethod
    cdef PhysicsLayer _view(c_PhysicsLayer* ptr, object world):
        cdef PhysicsLayer v = PhysicsLayer.__new__(PhysicsLayer)
        v._physics_ptr = ptr
        v._init_view(<c_BaseLayer*>ptr, world)
        return v

    # ------------------------------------------------------------------------------------------------------------------
    # Static mechanical properties
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
        """Complex potential Love number k (stored value; zero until assigned)."""
        cdef cpp_complex[double] k = self._physics_ptr.get_love_number_k()
        return complex(k.real(), k.imag())

    @property
    def love_number_h(self) -> complex:
        """Complex radial displacement Love number h (stored value; zero until assigned)."""
        cdef cpp_complex[double] h = self._physics_ptr.get_love_number_h()
        return complex(h.real(), h.imag())

    @property
    def love_number_l(self) -> complex:
        """Complex tangential displacement Love number l (stored value; zero until assigned)."""
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
        """True if this layer is solid, False for liquid. Used by the radial Love-number solver."""
        return bool(self._physics_ptr.get_is_solid())

    @is_solid.setter
    def is_solid(self, value: bool):
        self._physics_ptr.set_is_solid(<cpp_bool>bool(value))

    @property
    def is_static(self) -> bool:
        """True if the static (no-dynamic-terms) approximation is used. Used by the radial solver."""
        return bool(self._physics_ptr.get_is_static())

    @is_static.setter
    def is_static(self, value: bool):
        self._physics_ptr.set_is_static(<cpp_bool>bool(value))

    @property
    def is_incompressible(self) -> bool:
        """True if the incompressible approximation is used. Used by the radial solver."""
        return bool(self._physics_ptr.get_is_incompressible())

    @is_incompressible.setter
    def is_incompressible(self, value: bool):
        self._physics_ptr.set_is_incompressible(<cpp_bool>bool(value))

    # ------------------------------------------------------------------------------------------------------------------
    # Material state
    # ------------------------------------------------------------------------------------------------------------------
    @property
    def temperature(self) -> float:
        """Layer temperature [K] at which the viscosity and melt models are evaluated."""
        return self._physics_ptr.get_temperature()

    @temperature.setter
    def temperature(self, double value):
        self._physics_ptr.set_temperature(value)

    @property
    def use_thermal_eos(self) -> bool:
        """True if the EOS model receives the temperature (thermal density and bulk modulus)."""
        return bool(self._physics_ptr.get_use_thermal_eos())

    @use_thermal_eos.setter
    def use_thermal_eos(self, value: bool):
        self._physics_ptr.set_use_thermal_eos(<cpp_bool>bool(value))

    @property
    def shear_modulus_pressure_derivative(self) -> float:
        """Pressure derivative of the static shear modulus [Pa/Pa]."""
        return self._physics_ptr.get_shear_modulus_pressure_derivative()

    @property
    def shear_modulus_temperature_derivative(self) -> float:
        """Temperature derivative of the static shear modulus [Pa/K]."""
        return self._physics_ptr.get_shear_modulus_temperature_derivative()

    @property
    def shear_modulus_reference_temperature(self) -> float:
        """Temperature [K] at which ``shear_modulus_static`` applies."""
        return self._physics_ptr.get_shear_modulus_reference_temperature()

    def calc_material_state(self, double pressure, temperature=None, frequency=None, double radius=0.0) -> dict:
        """Material properties at a pressure [Pa] and temperature [K] from the layer's attached models.

        Evaluates, in order, the EOS density and bulk modulus (athermal unless ``use_thermal_eos``), the static
        moduli (an EOS-provided value takes precedence over the layer's shear law and bulk constant), the
        viscosities (EOS table, else the viscosity model at the temperature and pressure, else the layer
        constant), the partial-melt model, and the rheologies.

        Parameters
        ----------
        pressure : float
            Pressure [Pa].
        temperature : float, optional
            Temperature [K]. ``None`` uses the layer's ``temperature``.
        frequency : float, optional
            Forcing frequency [rad/s] for the complex moduli. ``None`` skips the rheologies, so the complex
            moduli are the post-melt static moduli.
        radius : float, optional
            Radius [m], read only by an interpolated EOS model. Default ``0.0``.

        Returns
        -------
        dict
            ``density`` [kg/m^3], ``melt_fraction``, the ``premelt_`` and post-melt ``shear_modulus``,
            ``bulk_modulus`` [Pa], ``shear_viscosity``, ``bulk_viscosity`` [Pa·s], and ``complex_shear_modulus``,
            ``complex_bulk_modulus`` [Pa].
        """
        cdef double temperature_value = (
            self._physics_ptr.get_temperature() if temperature is None else <double>temperature)
        cdef double frequency_value = d_NAN if frequency is None else <double>frequency
        cdef c_MaterialState state
        self._physics_ptr.calc_material_state(radius, pressure, temperature_value, frequency_value, state)
        return {
            "density":                 state.density,
            "melt_fraction":           state.melt_fraction,
            "premelt_shear_modulus":   state.premelt_shear_modulus,
            "premelt_bulk_modulus":    state.premelt_bulk_modulus,
            "premelt_shear_viscosity": state.premelt_shear_viscosity,
            "premelt_bulk_viscosity":  state.premelt_bulk_viscosity,
            "shear_modulus":           state.shear_modulus,
            "bulk_modulus":            state.bulk_modulus,
            "shear_viscosity":         state.shear_viscosity,
            "bulk_viscosity":          state.bulk_viscosity,
            "complex_shear_modulus":   complex(
                state.complex_shear_modulus.real(), state.complex_shear_modulus.imag()),
            "complex_bulk_modulus":    complex(
                state.complex_bulk_modulus.real(), state.complex_bulk_modulus.imag()),
        }

    # ------------------------------------------------------------------------------------------------------------------
    # Rheology attachment
    # ------------------------------------------------------------------------------------------------------------------
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
        if rheology._rheology_ptr.get() == NULL:
            raise ValueError(
                "This rheology model holds no C++ object (already attached or moved).")
        self._physics_ptr.set_shear_rheology(move(rheology._rheology_ptr))

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

        Ownership of the C++ model moves out of ``viscosity``, which is left an empty shell and must not be reused.

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

        Ownership of the C++ model moves out of ``viscosity``, which is left an empty shell and must not be reused.

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

        Ownership of the C++ model moves out of ``partial_melt``, which is left an empty shell and must not be
        reused.

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
    def _apply_complex(self, radius, double frequency, cpp_bool is_shear):
        # Radius-resolved complex modulus: float -> complex; np.ndarray -> complex np.ndarray (same shape).
        cdef cnp.ndarray in_arr
        cdef cnp.ndarray out_arr
        cdef double[::1] flat_in
        cdef double complex[::1] flat_out
        cdef cpp_complex[double] value
        cdef Py_ssize_t i, n
        if isinstance(radius, np.ndarray):
            in_arr  = np.ascontiguousarray(radius, dtype=np.float64)
            out_arr = np.empty_like(in_arr, dtype=np.complex128)
            flat_in = in_arr.reshape(-1)
            flat_out = out_arr.reshape(-1)
            n = flat_in.shape[0]
            for i in range(n):
                if is_shear:
                    value = self._physics_ptr.calc_complex_shear_modulus(flat_in[i], frequency)
                else:
                    value = self._physics_ptr.calc_complex_bulk_modulus(flat_in[i], frequency)
                flat_out[i] = value.real() + 1j * value.imag()
            return out_arr
        if is_shear:
            value = self._physics_ptr.calc_complex_shear_modulus(<double>radius, frequency)
        else:
            value = self._physics_ptr.calc_complex_bulk_modulus(<double>radius, frequency)
        return complex(value.real(), value.imag())

    def calc_complex_shear_modulus(self, first_arg, frequency=None):
        """Complex shear modulus [Pa]: layer-constant or radius-resolved.

        ``calc_complex_shear_modulus(frequency)`` applies the shear rheology to the layer-constant static shear
        modulus and viscosity. That static viscosity is NaN unless it was given at construction, so a viscous
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
        cdef cpp_complex[double] result
        if frequency is None:
            result = self._physics_ptr.calc_complex_shear_modulus(<double>first_arg)
            return complex(result.real(), result.imag())
        return self._apply_complex(first_arg, <double>frequency, True)

    def calc_complex_bulk_modulus(self, first_arg, frequency=None):
        """Complex bulk modulus [Pa]: layer-constant or radius-resolved.

        ``calc_complex_bulk_modulus(frequency)`` applies the bulk rheology to the layer-constant static bulk
        modulus and viscosity. That static viscosity is NaN unless it was given at construction, so a viscous
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
        cdef cpp_complex[double] result
        if frequency is None:
            result = self._physics_ptr.calc_complex_bulk_modulus(<double>first_arg)
            return complex(result.real(), result.imag())
        return self._apply_complex(first_arg, <double>frequency, False)

    # ------------------------------------------------------------------------------------------------------------------
    # Config
    # ------------------------------------------------------------------------------------------------------------------
    cpdef dict get_config_dict(self):
        """Return all configuration values as a Python dict (MKS).

        Returns
        -------
        dict
            The BaseLayer keys plus ``shear_modulus_static``, ``bulk_modulus_static``,
            ``shear_viscosity_static``, ``bulk_viscosity_static``, the radial-solver flags ``is_solid``,
            ``is_static``, and ``is_incompressible``, the material-state keys (``temperature_k``, the three
            shear-law keys, ``use_thermal_eos``), the six Love number components, and one sub-table per
            attached model (``shear_rheology``, ``bulk_rheology``, ``shear_viscosity``, ``bulk_viscosity``,
            ``partial_melt``).
        """
        d = BaseLayer.get_config_dict(self)
        d["shear_modulus_static_pa"]      = self._physics_ptr.get_shear_modulus_static()
        d["bulk_modulus_static_pa"]       = self._physics_ptr.get_bulk_modulus_static()
        d["shear_viscosity_static_pas"]   = self._physics_ptr.get_shear_viscosity_static()
        d["bulk_viscosity_static_pas"]    = self._physics_ptr.get_bulk_viscosity_static()
        d["is_solid"]          = bool(self._physics_ptr.get_is_solid())
        d["is_static"]         = bool(self._physics_ptr.get_is_static())
        d["is_incompressible"] = bool(self._physics_ptr.get_is_incompressible())
        d["temperature_k"]     = self._physics_ptr.get_temperature()
        d["shear_modulus_pressure_derivative"] = self._physics_ptr.get_shear_modulus_pressure_derivative()
        d["shear_modulus_temperature_derivative_pa_k"] = (
            self._physics_ptr.get_shear_modulus_temperature_derivative())
        d["shear_modulus_reference_temperature_k"] = (
            self._physics_ptr.get_shear_modulus_reference_temperature())
        d["use_thermal_eos"]   = bool(self._physics_ptr.get_use_thermal_eos())
        cdef c_LoveNumbers ln = self._physics_ptr.get_love_numbers()
        d["love_number_k_re"] = ln.k.real()
        d["love_number_k_im"] = ln.k.imag()
        d["love_number_h_re"] = ln.h.real()
        d["love_number_h_im"] = ln.h.imag()
        d["love_number_l_re"] = ln.l.real()
        d["love_number_l_im"] = ln.l.imag()
        # Attached models, keyed the way the world builder reads them.
        cdef const c_PhysicsBase* model_ptr
        model_ptr = <const c_PhysicsBase*>self._physics_ptr.get_shear_rheology_model()
        if model_ptr != NULL:
            d["shear_rheology"] = cy_physics_model_config(model_ptr)
        model_ptr = <const c_PhysicsBase*>self._physics_ptr.get_bulk_rheology_model()
        if model_ptr != NULL:
            d["bulk_rheology"] = cy_physics_model_config(model_ptr)
        model_ptr = <const c_PhysicsBase*>self._physics_ptr.get_shear_viscosity_model()
        if model_ptr != NULL:
            d["shear_viscosity"] = cy_physics_model_config(model_ptr)
        model_ptr = <const c_PhysicsBase*>self._physics_ptr.get_bulk_viscosity_model()
        if model_ptr != NULL:
            d["bulk_viscosity"] = cy_physics_model_config(model_ptr)
        model_ptr = <const c_PhysicsBase*>self._physics_ptr.get_partial_melt_model()
        if model_ptr != NULL:
            d["partial_melt"] = cy_physics_model_config(model_ptr)
        return d
