# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
material_eos.pyx
Cython/Python wrapper for TidalPy's material EOS model hierarchy.

A material EOS model returns a material's density [kg/m^3] from the local pressure
[Pa] (analytic models) or radius [m] (interpolated model). Models are attached to
layers and supply the per-layer density source for the whole-planet EOS solve.

Models:
- ConstantDensityEOS (alias "constant") — incompressible.
- BirchMurnaghanEOS  (alias "bm")       — 3rd-order Birch-Murnaghan.
- VinetEOS           (alias "vinet")    — Vinet/UBER EOS.
- InterpolatedEOS    (alias "interp")   — density(radius) lookup table.
"""

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector
from libcpp.utility cimport move

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_TidalPyBaseClass

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

    def calc_density(self, double pressure, double temperature=0.0,
                     double radius=0.0) -> float:
        """Material density [kg/m^3] at the given pressure [Pa] (and radius [m]).

        Analytic models use pressure; the interpolated model uses radius.
        Temperature is accepted for API uniformity (currently unused; isothermal).
        """
        return self._eos_ptr.get().calc_density(pressure, temperature, radius)

    def calc_static_shear_modulus(self, double radius) -> float:
        """Static shear modulus [Pa] at radius [m]; NaN if the model has no shear table.

        Only the interpolated model carries radius-varying viscoelastic tables; the
        analytic models return NaN (the world solve then uses the layer constant).
        """
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

    def __init__(self, double reference_density=3500.0):
        cdef c_MaterialEOSConfig config
        config.reference_density = reference_density
        # Build through the C++ factory (make_unique) and adopt ownership; no raw new/delete.
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

    def __init__(self, double reference_density=3500.0,
                 double reference_bulk_modulus=1.0e11,
                 double bulk_modulus_derivative=4.0,
                 invert_rtol=None, invert_max_iters=None):
        # A default-constructed config carries the C++ default inversion settings;
        # only override them when the caller explicitly supplies a value.
        cdef c_MaterialEOSConfig config
        config.reference_density   = reference_density
        config.reference_bulk_modulus = reference_bulk_modulus
        config.bulk_modulus_derivative   = bulk_modulus_derivative
        if invert_rtol is not None:
            config.invert_rtol = <double>invert_rtol
        if invert_max_iters is not None:
            config.invert_max_iters = <int>invert_max_iters
        # Build through the C++ factory (make_unique) and adopt ownership; no raw new/delete.
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

    def __init__(self, double reference_density=3500.0,
                 double reference_bulk_modulus=1.0e11,
                 double bulk_modulus_derivative=4.0,
                 invert_rtol=None, invert_max_iters=None):
        # A default-constructed config carries the C++ default inversion settings;
        # only override them when the caller explicitly supplies a value.
        cdef c_MaterialEOSConfig config
        config.reference_density   = reference_density
        config.reference_bulk_modulus = reference_bulk_modulus
        config.bulk_modulus_derivative   = bulk_modulus_derivative
        if invert_rtol is not None:
            config.invert_rtol = <double>invert_rtol
        if invert_max_iters is not None:
            config.invert_max_iters = <int>invert_max_iters
        # Build through the C++ factory (make_unique) and adopt ownership; no raw new/delete.
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
            bulk_viscosity=None):
        cdef c_MaterialEOSConfig config
        config.radius      = <vector[double]>radius
        config.density = <vector[double]>density
        if config.radius.size() != config.density.size():
            raise ValueError("radius and density must have the same length.")
        # Optional radius-varying viscoelastic tables; each must match radius length.
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
        # Build through the C++ factory (make_unique) and adopt ownership; no raw new/delete.
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
def make_material_eos(str model_name, dict config=None) -> MaterialEOSBase:
    """Build a material EOS model by name, returning the matching rich subclass.

    Parameters
    ----------
    model_name : str
        One of ``"constant"``, ``"bm"``/``"birch_murnaghan"``, ``"vinet"``,
        ``"interpolate"`` (case-insensitive; aliases accepted).
    config : dict, optional
        Model parameters: ``reference_density``, ``reference_bulk_modulus``,
        ``bulk_modulus_derivative`` (analytic models); ``radius`` and
        ``density`` sequences (interpolated model).

    Returns
    -------
    MaterialEOSBase
        The concrete model subclass.

    Raises
    ------
    ValueError
        If the model name is unknown.
    """
    if config is None:
        config = {}
    # A default-constructed config carries the C++ defaults; only override the
    # fields the caller actually supplies (keeps defaults in one place: the C++ struct).
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
    if "radius_m" in config:
        cfg.radius = <vector[double]>config["radius_m"]
    if "density_kg_m3" in config:
        cfg.density = <vector[double]>config["density_kg_m3"]
    # Optional radius-varying viscoelastic tables for the interpolated model.
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
