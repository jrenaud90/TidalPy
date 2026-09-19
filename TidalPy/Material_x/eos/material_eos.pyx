# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrappers for the material EOS models.

A model returns a material's density [kg/m^3] from the local pressure [Pa] (analytic models) or radius [m]
(interpolated model) and supplies a layer's density source for the whole-planet EOS solve. Models:
ConstantDensityEOS ("constant"), BirchMurnaghanEOS ("bm"), VinetEOS ("vinet"), InterpolatedEOS ("interp").

Every model takes a thermal expansivity and a reference temperature. Birch-Murnaghan and Vinet add the thermal
pressure alpha0 K0 (T - T_ref) to their pressure law; the constant and interpolated models scale their density by
exp(-alpha0 (T - T_ref)). A zero expansivity (the default) or no temperature gives the athermal EOS.
"""

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector
from libcpp.utility cimport move

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport d_NAN, set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_TidalPyBaseClass
from TidalPy.Utilities_x.classes_x.classes import check_config_keys

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


# =====================================================================================================================
# Free analytic pressure laws (Python-accessible for cross-checks)
# =====================================================================================================================
def birch_murnaghan_pressure(double eta, double reference_bulk_modulus,
                             double bulk_modulus_derivative=4.0) -> float:
    """3rd-order Birch-Murnaghan pressure [Pa] at compression eta = rho/rho0."""
    return eos_bm_pressure(eta, reference_bulk_modulus, bulk_modulus_derivative)


def vinet_pressure(double eta, double reference_bulk_modulus,
                   double bulk_modulus_derivative=4.0) -> float:
    """Vinet pressure [Pa] at compression eta = rho/rho0."""
    return eos_vinet_pressure(eta, reference_bulk_modulus, bulk_modulus_derivative)


# =====================================================================================================================
# MaterialEOSBase
# =====================================================================================================================
cdef class MaterialEOSBase(PhysicsBase):
    """Abstract base for material EOS models. Instantiate a concrete subclass."""

    def __cinit__(self, *args, **kwargs):
        pass  # unique_ptr<c_MaterialEOSBase> auto-inits to nullptr

    def __init__(self, *args, **kwargs):
        raise TypeError(
            "MaterialEOSBase is abstract; instantiate a concrete model "
            "(ConstantDensityEOS, BirchMurnaghanEOS, VinetEOS, InterpolatedEOS).")

    def __dealloc__(self):
        self._eos_ptr.reset()
        self._ptr = NULL

    def calc_density(self, double pressure, temperature=None, double radius=0.0) -> float:
        """Density [kg/m^3] at pressure [Pa] (analytic models) or radius [m] (interpolated model).

        ``temperature`` [K] enters through the model's thermal expansivity; ``None`` gives the athermal density.
        """
        cdef double temperature_value = d_NAN if temperature is None else <double>temperature
        return self._eos_ptr.get().calc_density(pressure, temperature_value, radius)

    def calc_bulk_modulus(self, double pressure, temperature=None, double radius=0.0) -> float:
        """Isothermal bulk modulus K = rho dP/drho [Pa] at pressure [Pa] and temperature [K].

        Birch-Murnaghan and Vinet return the modulus of their pressure law at the solved compression. The other
        models return their bulk table at ``radius`` [m], or NaN when they have none.
        """
        cdef double temperature_value = d_NAN if temperature is None else <double>temperature
        return self._eos_ptr.get().calc_bulk_modulus(pressure, temperature_value, radius)

    @property
    def thermal_expansion(self) -> float:
        """Thermal expansivity alpha0 [1/K]; zero is the athermal EOS."""
        return self._eos_ptr.get().get_thermal_expansion()

    @property
    def reference_temperature(self) -> float:
        """Temperature [K] at which the reference density (and bulk modulus) apply."""
        return self._eos_ptr.get().get_reference_temperature()

    def calc_static_shear_modulus(self, double radius) -> float:
        """Static shear modulus [Pa] at radius [m]; NaN unless the model carries a shear table."""
        return self._eos_ptr.get().calc_static_shear_modulus(radius)

    def calc_static_bulk_modulus(self, double radius) -> float:
        """Static bulk modulus [Pa] at radius [m]; NaN if the model has no bulk table."""
        return self._eos_ptr.get().calc_static_bulk_modulus(radius)

    def calc_shear_viscosity(self, double radius) -> float:
        """Shear viscosity [Pa s] at radius [m]; NaN if the model has no shear-viscosity table."""
        return self._eos_ptr.get().calc_shear_viscosity(radius)

    def calc_bulk_viscosity(self, double radius) -> float:
        """Bulk viscosity [Pa s] at radius [m]; NaN if the model has no bulk-viscosity table."""
        return self._eos_ptr.get().calc_bulk_viscosity(radius)


# =====================================================================================================================
# EOS Models
# =====================================================================================================================
cdef class ConstantDensityEOS(MaterialEOSBase):
    """Incompressible (uniform) density EOS."""

    def __cinit__(self, *args, **kwargs):
        self._constant_ptr = NULL

    def __init__(self, double reference_density=3500.0, double thermal_expansion=0.0, reference_temperature=None):
        # None keeps the C++ default reference temperature.
        cdef c_MaterialEOSConfig config
        config.reference_density = reference_density
        config.thermal_expansion = thermal_expansion
        if reference_temperature is not None:
            config.reference_temperature = <double>reference_temperature
        cdef unique_ptr[c_MaterialEOSBase] ptr = c_find_material_eos(
            c_MaterialEOSModel.Constant, config)
        self._constant_ptr = <c_ConstantDensityEOS*>ptr.get()
        self._eos_ptr      = move(ptr)
        self._ptr          = <c_TidalPyBaseClass*>self._eos_ptr.get()

    def __dealloc__(self):
        self._constant_ptr = NULL

    @property
    def reference_density(self) -> float:
        """Reference (uniform) density [kg/m^3]."""
        return self._constant_ptr.get_reference_density()


cdef class BirchMurnaghanEOS(MaterialEOSBase):
    """3rd-order Birch-Murnaghan EOS; density from pressure."""

    def __cinit__(self, *args, **kwargs):
        self._bm_ptr = NULL

    def __init__(
            self,
            double reference_density=3500.0,
            double reference_bulk_modulus=1.0e11,
            double bulk_modulus_derivative=4.0,
            invert_rtol=None,
            invert_max_iters=None,
            double thermal_expansion=0.0,
            reference_temperature=None):
        # None keeps the C++ default inversion settings and reference temperature.
        cdef c_MaterialEOSConfig config
        config.reference_density   = reference_density
        config.reference_bulk_modulus = reference_bulk_modulus
        config.bulk_modulus_derivative   = bulk_modulus_derivative
        config.thermal_expansion = thermal_expansion
        if reference_temperature is not None:
            config.reference_temperature = <double>reference_temperature
        if invert_rtol is not None:
            config.invert_rtol = <double>invert_rtol
        if invert_max_iters is not None:
            config.invert_max_iters = <int>invert_max_iters
        cdef unique_ptr[c_MaterialEOSBase] ptr = c_find_material_eos(
            c_MaterialEOSModel.BirchMurnaghan, config)
        self._bm_ptr  = <c_BirchMurnaghanEOS*>ptr.get()
        self._eos_ptr = move(ptr)
        self._ptr     = <c_TidalPyBaseClass*>self._eos_ptr.get()

    def __dealloc__(self):
        self._bm_ptr = NULL

    @property
    def reference_density(self) -> float:
        """Reference density rho0 [kg/m^3]."""
        return self._bm_ptr.get_reference_density()

    @property
    def reference_bulk_modulus(self) -> float:
        """Reference bulk modulus K0 [Pa]."""
        return self._bm_ptr.get_reference_bulk_modulus()

    @property
    def bulk_modulus_derivative(self) -> float:
        """Pressure derivative of the bulk modulus K0' [dimensionless]."""
        return self._bm_ptr.get_bulk_modulus_derivative()

    @property
    def invert_rtol(self) -> float:
        """Relative convergence tolerance for the density-from-pressure inversion."""
        return self._bm_ptr.get_invert_rtol()

    @property
    def invert_max_iters(self) -> int:
        """Hard iteration cap (termination safeguard) for the inversion."""
        return self._bm_ptr.get_invert_max_iters()


cdef class VinetEOS(MaterialEOSBase):
    """Vinet (universal) EOS; density from pressure."""

    def __cinit__(self, *args, **kwargs):
        self._vinet_ptr = NULL

    def __init__(
            self,
            double reference_density=3500.0,
            double reference_bulk_modulus=1.0e11,
            double bulk_modulus_derivative=4.0,
            invert_rtol=None,
            invert_max_iters=None,
            double thermal_expansion=0.0,
            reference_temperature=None):
        # None keeps the C++ default inversion settings and reference temperature.
        cdef c_MaterialEOSConfig config
        config.reference_density   = reference_density
        config.reference_bulk_modulus = reference_bulk_modulus
        config.bulk_modulus_derivative   = bulk_modulus_derivative
        config.thermal_expansion = thermal_expansion
        if reference_temperature is not None:
            config.reference_temperature = <double>reference_temperature
        if invert_rtol is not None:
            config.invert_rtol = <double>invert_rtol
        if invert_max_iters is not None:
            config.invert_max_iters = <int>invert_max_iters
        cdef unique_ptr[c_MaterialEOSBase] ptr = c_find_material_eos(
            c_MaterialEOSModel.Vinet, config)
        self._vinet_ptr = <c_VinetEOS*>ptr.get()
        self._eos_ptr   = move(ptr)
        self._ptr       = <c_TidalPyBaseClass*>self._eos_ptr.get()

    def __dealloc__(self):
        self._vinet_ptr = NULL

    @property
    def reference_density(self) -> float:
        """Reference density rho0 [kg/m^3]."""
        return self._vinet_ptr.get_reference_density()

    @property
    def reference_bulk_modulus(self) -> float:
        """Reference bulk modulus K0 [Pa]."""
        return self._vinet_ptr.get_reference_bulk_modulus()

    @property
    def bulk_modulus_derivative(self) -> float:
        """Pressure derivative of the bulk modulus K0' [dimensionless]."""
        return self._vinet_ptr.get_bulk_modulus_derivative()

    @property
    def invert_rtol(self) -> float:
        """Relative convergence tolerance for the density-from-pressure inversion."""
        return self._vinet_ptr.get_invert_rtol()

    @property
    def invert_max_iters(self) -> int:
        """Hard iteration cap (termination safeguard) for the inversion."""
        return self._vinet_ptr.get_invert_max_iters()


cdef class InterpolatedEOS(MaterialEOSBase):
    """density(radius) lookup table (e.g. PREM-style profiles)."""

    def __cinit__(self, *args, **kwargs):
        self._interp_ptr = NULL

    def __init__(
            self,
            radius,
            density,
            shear_modulus=None,
            bulk_modulus=None,
            shear_viscosity=None,
            bulk_viscosity=None,
            double thermal_expansion=0.0,
            reference_temperature=None):
        cdef c_MaterialEOSConfig config
        config.thermal_expansion = thermal_expansion
        if reference_temperature is not None:
            config.reference_temperature = <double>reference_temperature
        config.radius      = <vector[double]>radius
        config.density = <vector[double]>density
        if config.radius.size() != config.density.size():
            raise ValueError("radius and density must have the same length.")
        if shear_modulus is not None:
            config.shear_modulus = <vector[double]>shear_modulus
            if config.shear_modulus.size() != config.radius.size():
                raise ValueError("shear_modulus must match radius in length.")
        if bulk_modulus is not None:
            config.bulk_modulus = <vector[double]>bulk_modulus
            if config.bulk_modulus.size() != config.radius.size():
                raise ValueError("bulk_modulus must match radius in length.")
        if shear_viscosity is not None:
            config.shear_viscosity = <vector[double]>shear_viscosity
            if config.shear_viscosity.size() != config.radius.size():
                raise ValueError("shear_viscosity must match radius in length.")
        if bulk_viscosity is not None:
            config.bulk_viscosity = <vector[double]>bulk_viscosity
            if config.bulk_viscosity.size() != config.radius.size():
                raise ValueError("bulk_viscosity must match radius in length.")
        cdef unique_ptr[c_MaterialEOSBase] ptr = c_find_material_eos(
            c_MaterialEOSModel.Interpolated, config)
        self._interp_ptr = <c_InterpolatedEOS*>ptr.get()
        self._eos_ptr    = move(ptr)
        self._ptr        = <c_TidalPyBaseClass*>self._eos_ptr.get()

    def __dealloc__(self):
        self._interp_ptr = NULL

    @property
    def num_points(self) -> int:
        """Number of (radius, density) table points."""
        return self._interp_ptr.get_num_points()


# =====================================================================================================================
# Factory
# =====================================================================================================================
# Every config key some material EOS model reads; make_material_eos rejects anything else.
MATERIAL_EOS_CONFIG_KEYS = frozenset({
    "reference_density_kg_m3", "reference_bulk_modulus_pa", "bulk_modulus_derivative", "invert_rtol",
    "invert_max_iters", "radius_m", "density_kg_m3", "shear_modulus_pa", "bulk_modulus_pa",
    "shear_viscosity_pas", "bulk_viscosity_pas", "thermal_expansion_1_k", "reference_temperature_k"})


def make_material_eos(str model_name, dict config=None) -> MaterialEOSBase:
    """Build a material EOS model by name, returning the matching rich subclass.

    Parameters
    ----------
    model_name : str
        ``"constant"``, ``"bm"``/``"birch_murnaghan"``, ``"vinet"``, or ``"interpolate"`` (case-insensitive).
    config : dict, optional
        ``reference_density_kg_m3``, ``reference_bulk_modulus_pa``, ``bulk_modulus_derivative``, ``invert_rtol``,
        ``invert_max_iters`` (analytic models); ``radius_m`` and ``density_kg_m3`` sequences plus the optional
        ``shear_modulus_pa``, ``bulk_modulus_pa``, ``shear_viscosity_pas``, ``bulk_viscosity_pas`` tables
        (interpolated model); ``thermal_expansion_1_k`` and ``reference_temperature_k`` (every model).

    Returns
    -------
    MaterialEOSBase
        The concrete model subclass.

    Raises
    ------
    ValueError
        Unknown model name, or a ``config`` key that no material EOS model reads.
    """
    check_config_keys(config, MATERIAL_EOS_CONFIG_KEYS, "material EOS")
    if config is None:
        config = {}
    # Only the supplied keys override the C++ struct defaults.
    cdef c_MaterialEOSConfig cfg
    if "reference_density_kg_m3" in config:
        cfg.reference_density = config["reference_density_kg_m3"]
    if "reference_bulk_modulus_pa" in config:
        cfg.reference_bulk_modulus = config["reference_bulk_modulus_pa"]
    if "bulk_modulus_derivative" in config:
        cfg.bulk_modulus_derivative = config["bulk_modulus_derivative"]
    if "invert_rtol" in config:
        cfg.invert_rtol = config["invert_rtol"]
    if "invert_max_iters" in config:
        cfg.invert_max_iters = config["invert_max_iters"]
    if "thermal_expansion_1_k" in config:
        cfg.thermal_expansion = config["thermal_expansion_1_k"]
    if "reference_temperature_k" in config:
        cfg.reference_temperature = config["reference_temperature_k"]
    if "radius_m" in config:
        cfg.radius = <vector[double]>config["radius_m"]
    if "density_kg_m3" in config:
        cfg.density = <vector[double]>config["density_kg_m3"]
    if "shear_modulus_pa" in config:
        cfg.shear_modulus = <vector[double]>config["shear_modulus_pa"]
    if "bulk_modulus_pa" in config:
        cfg.bulk_modulus = <vector[double]>config["bulk_modulus_pa"]
    if "shear_viscosity_pas" in config:
        cfg.shear_viscosity = <vector[double]>config["shear_viscosity_pas"]
    if "bulk_viscosity_pas" in config:
        cfg.bulk_viscosity = <vector[double]>config["bulk_viscosity_pas"]

    cdef c_MaterialEOSModel model = c_material_eos_model_from_name(model_name.encode("utf-8"))
    cdef unique_ptr[c_MaterialEOSBase] ptr = c_find_material_eos(model, cfg)

    cdef ConstantDensityEOS constant_eos
    cdef BirchMurnaghanEOS  bm_eos
    cdef VinetEOS           vinet_eos
    cdef InterpolatedEOS    interp_eos
    if model == c_MaterialEOSModel.Constant:
        constant_eos = ConstantDensityEOS.__new__(ConstantDensityEOS)
        constant_eos._constant_ptr = <c_ConstantDensityEOS*>ptr.get()
        constant_eos._eos_ptr = move(ptr)
        constant_eos._ptr = <c_TidalPyBaseClass*>constant_eos._eos_ptr.get()
        return constant_eos
    elif model == c_MaterialEOSModel.BirchMurnaghan:
        bm_eos = BirchMurnaghanEOS.__new__(BirchMurnaghanEOS)
        bm_eos._bm_ptr = <c_BirchMurnaghanEOS*>ptr.get()
        bm_eos._eos_ptr = move(ptr)
        bm_eos._ptr = <c_TidalPyBaseClass*>bm_eos._eos_ptr.get()
        return bm_eos
    elif model == c_MaterialEOSModel.Vinet:
        vinet_eos = VinetEOS.__new__(VinetEOS)
        vinet_eos._vinet_ptr = <c_VinetEOS*>ptr.get()
        vinet_eos._eos_ptr = move(ptr)
        vinet_eos._ptr = <c_TidalPyBaseClass*>vinet_eos._eos_ptr.get()
        return vinet_eos
    else:
        interp_eos = InterpolatedEOS.__new__(InterpolatedEOS)
        interp_eos._interp_ptr = <c_InterpolatedEOS*>ptr.get()
        interp_eos._eos_ptr = move(ptr)
        interp_eos._ptr = <c_TidalPyBaseClass*>interp_eos._eos_ptr.get()
        return interp_eos
