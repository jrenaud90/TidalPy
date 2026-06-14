# distutils: language = c++
"""
material_eos.pxd
Cython declarations for TidalPy's material EOS model hierarchy.

Exports the C++ EOS models, the combined config struct, the enum factory, and the
Python wrapper classes so other extensions (layers, worlds) can cimport and build
or attach EOS models at C speed.
"""

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector

from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_PhysicsBase


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "material_eos_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_MaterialEOSBase(c_PhysicsBase):
        double calc_density(double pressure_pa, double temperature_k, double radius_m) const

    cdef cppclass c_MaterialEOSConfig:
        double reference_density_kg_m3
        double reference_bulk_modulus_pa
        double bulk_modulus_derivative
        double invert_rtol
        int    invert_max_iters
        vector[double] radius_m
        vector[double] density_kg_m3

    cdef cppclass c_ConstantDensityEOS(c_MaterialEOSBase):
        c_ConstantDensityEOS() except +
        c_ConstantDensityEOS(const c_MaterialEOSConfig& cfg) except +
        double get_reference_density() const

    cdef cppclass c_BirchMurnaghanEOS(c_MaterialEOSBase):
        c_BirchMurnaghanEOS() except +
        c_BirchMurnaghanEOS(const c_MaterialEOSConfig& cfg) except +
        double get_reference_density()       const
        double get_reference_bulk_modulus()  const
        double get_bulk_modulus_derivative() const
        double get_invert_rtol()             const
        int    get_invert_max_iters()        const

    cdef cppclass c_VinetEOS(c_MaterialEOSBase):
        c_VinetEOS() except +
        c_VinetEOS(const c_MaterialEOSConfig& cfg) except +
        double get_reference_density()       const
        double get_reference_bulk_modulus()  const
        double get_bulk_modulus_derivative() const
        double get_invert_rtol()             const
        int    get_invert_max_iters()        const

    cdef cppclass c_InterpolatedEOS(c_MaterialEOSBase):
        c_InterpolatedEOS() except +
        c_InterpolatedEOS(const c_MaterialEOSConfig& cfg) except +
        size_t get_num_points() const

    # Free analytic pressure laws (for tests / cross-checks).
    double eos_bm_pressure(double eta, double K0, double K0_prime)
    double eos_vinet_pressure(double eta, double K0, double K0_prime)

    cdef enum class c_MaterialEOSModel:
        Constant
        BirchMurnaghan
        Vinet
        Interpolated

    c_MaterialEOSModel c_material_eos_model_from_name(const string& model_name) except +
    unique_ptr[c_MaterialEOSBase] c_find_material_eos(
        c_MaterialEOSModel model, const c_MaterialEOSConfig& cfg) except +


# =====================================================================================================================
# Cython wrapper class declarations
# =====================================================================================================================
cdef class MaterialEOSBase(PhysicsBase):
    cdef unique_ptr[c_MaterialEOSBase] _eos_ptr   # owns the most-derived C++ model object
    cpdef dict get_config_dict(self)


cdef class ConstantDensityEOS(MaterialEOSBase):
    cdef c_ConstantDensityEOS* _constant_ptr   # non-owning; ownership via MaterialEOSBase._eos_ptr
    cpdef dict get_config_dict(self)


cdef class BirchMurnaghanEOS(MaterialEOSBase):
    cdef c_BirchMurnaghanEOS* _bm_ptr          # non-owning; ownership via MaterialEOSBase._eos_ptr
    cpdef dict get_config_dict(self)


cdef class VinetEOS(MaterialEOSBase):
    cdef c_VinetEOS* _vinet_ptr                # non-owning; ownership via MaterialEOSBase._eos_ptr
    cpdef dict get_config_dict(self)


cdef class InterpolatedEOS(MaterialEOSBase):
    cdef c_InterpolatedEOS* _interp_ptr        # non-owning; ownership via MaterialEOSBase._eos_ptr
    cpdef dict get_config_dict(self)
