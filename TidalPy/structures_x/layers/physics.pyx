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
from TidalPy.structures_x.layers.base cimport BaseLayer, c_BaseLayer
from TidalPy.Tides_x.love.love cimport LoveNumbers, c_LoveNumbers
from TidalPy.rheology_x.rheology cimport RheologyBase
from TidalPy.viscosity_x.viscosity cimport ViscosityBase
from TidalPy.partial_melt_x.partial_melt cimport PartialMeltBase

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


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
    is_volume_fixed : bool, optional
        False lets the layer grow or shrink to hold its mass during an EOS solve. Default ``True``.
    is_tidal : bool, optional
        Whether this layer contributes to tidal dissipation. Default ``True``.
    tidal_scale : float, optional
        The layer's share of the planet in the quasi-homogeneous Love methods; ``None`` (default) takes
        its volume fraction. See ``BaseLayer``.
    love_number_k : complex, optional
        Potential Love number (placeholder). Default ``0+0j``.
    love_number_h : complex, optional
        Radial displacement Love number (placeholder). Default ``0+0j``.
    love_number_l : complex, optional
        Tangential displacement Love number (placeholder). Default ``0+0j``.
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
        Let the world's heat sources (this layer's radiogenics model among them) act inside the layer during a
        thermal EOS solve. Default ``False``.

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
            str    material_name          = "",
            cpp_bool is_tidal             = True,
            cpp_bool is_volume_fixed      = True,
            tidal_scale                   = None,
            complex love_number_k         = 0+0j,
            complex love_number_h         = 0+0j,
            complex love_number_l         = 0+0j,
            cpp_bool is_solid             = True,
            cpp_bool is_static            = True,
            cpp_bool is_incompressible    = False,
            double temperature            = 0.0,
            cpp_bool use_thermal_eos = False,
            cpp_bool use_heating     = False):
        cdef c_PhysicsConfig config
        config.name               = name.encode("utf-8")
        config.layer_index        = layer_index
        config.radius_inner       = radius_inner
        config.radius_outer       = radius_outer
        config.mass               = mass
        config.material_name      = material_name.encode("utf-8")
        config.is_tidal           = is_tidal
        config.is_volume_fixed    = is_volume_fixed
        config.tidal_scale        = d_NAN if tidal_scale is None else <double>tidal_scale
        config.love_numbers = c_LoveNumbers(
            cpp_complex[double](love_number_k.real, love_number_k.imag),
            cpp_complex[double](love_number_h.real, love_number_h.imag),
            cpp_complex[double](love_number_l.real, love_number_l.imag))
        config.is_solid          = is_solid
        config.is_static         = is_static
        config.is_incompressible = is_incompressible
        config.temperature       = temperature
        config.use_thermal_eos   = use_thermal_eos
        config.use_heating       = use_heating
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
    def use_heating(self) -> bool:
        """True if the world's heat sources act inside this layer during a thermal EOS solve."""
        return bool(self._physics_ptr.get_use_heating())

    @use_heating.setter
    def use_heating(self, value: bool):
        self._physics_ptr.set_use_heating(<cpp_bool>bool(value))

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
        if rheology._rheology_ptr.get() == NULL:
            raise ValueError(
                "This rheology model holds no C++ object (already attached or moved).")
        self._physics_ptr.set_bulk_rheology(move(rheology._rheology_ptr))
        rheology._ptr = NULL

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
        """Give the layer's material a viscosity model for its shear viscosity (before partial melt).

        A helper: the material owns this model, so it is handed to the layer's EOS model, which must already be
        attached (``set_eos`` first). Ownership of the C++ model moves out of the argument, which is left an empty
        shell and must not be reused.

        Raises
        ------
        ValueError
            If the model has already been attached or otherwise moved, or the layer has no EOS model yet.
        """
        if viscosity._visc_ptr.get() == NULL:
            raise ValueError(
                "This viscosity model holds no C++ object (already attached or moved).")
        if not self._layer_ptr.get().get_eos_set():
            raise ValueError(
                f"Attach an EOS model to layer '{self.name}' before giving it a viscosity model: the material owns it.")
        self._physics_ptr.set_shear_viscosity(move(viscosity._visc_ptr))
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
        if viscosity._visc_ptr.get() == NULL:
            raise ValueError(
                "This viscosity model holds no C++ object (already attached or moved).")
        if not self._layer_ptr.get().get_eos_set():
            raise ValueError(
                f"Attach an EOS model to layer '{self.name}' before giving it a viscosity model: the material owns it.")
        self._physics_ptr.set_bulk_viscosity(move(viscosity._visc_ptr))
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
        if partial_melt._melt_ptr.get() == NULL:
            raise ValueError(
                "This partial-melt model holds no C++ object (already attached or moved).")
        if not self._layer_ptr.get().get_eos_set():
            raise ValueError(
                f"Attach an EOS model to layer '{self.name}' before giving it a partial-melt model: "
                "the material owns it.")
        self._physics_ptr.set_partial_melt(move(partial_melt._melt_ptr))
        partial_melt._ptr = NULL

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
            # Every name in this loop is a C type, so it runs without the interpreter, as the matching
            # per-radius loops in worlds/layered.pyx and layers/base.pyx already do.
            with nogil:
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

    cpdef dict get_config_dict(self):
        """Return all configuration values as a Python dict (MKS).

        Returns
        -------
        dict
            The BaseLayer keys (with the ``material`` table of the attached EOS model, which carries the static
            constants, the shear law, and its viscosity and partial-melt models) plus the radial-solver flags
            ``is_solid``, ``is_static``, and ``is_incompressible``, ``temperature_k``, ``use_thermal_eos``,
            ``use_heating``, the six
            Love number components, and a sub-table for each attached rheology (``shear_rheology``,
            ``bulk_rheology``).
        """
        cdef dict d = BaseLayer.get_config_dict(self)
        d["is_solid"]          = bool(self._physics_ptr.get_is_solid())
        d["is_static"]         = bool(self._physics_ptr.get_is_static())
        d["is_incompressible"] = bool(self._physics_ptr.get_is_incompressible())
        d["temperature_k"]     = self._physics_ptr.get_temperature()
        d["use_thermal_eos"]   = bool(self._physics_ptr.get_use_thermal_eos())
        d["use_heating"]       = bool(self._physics_ptr.get_use_heating())
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
        return d
