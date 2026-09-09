# distutils: language = c++
"""
rheology.pxd
Cython declarations for TidalPy's rheology model hierarchy.

Exports the seven concrete C++ rheology models, the combined config struct, and
the Python wrapper classes so other extensions can cimport and build/attach
rheology models at C speed.

Usage::

    from TidalPy.rheology_x.rheology cimport (
        RheologyBase, Maxwell, c_RheologyBase, c_Maxwell, c_RheologyConfig)
"""

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector
from libcpp.complex cimport complex as cpp_complex

from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_PhysicsBase


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "rheology_base_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_RheologyBase(c_PhysicsBase):
        cpp_complex[double] calc_complex_modulus(
            double modulus_pa, double viscosity_pas, double frequency_rad_s) const
        void calc_complex_modulus_vectorize_modulus(
            const vector[double]& modulus_pa,
            const vector[double]& viscosity_pas,
            double frequency_rad_s,
            vector[cpp_complex[double]]& out_complex_modulus) except +
        void calc_complex_modulus_vectorize_frequency(
            double modulus_pa,
            double viscosity_pas,
            const vector[double]& frequency_rad_s,
            vector[cpp_complex[double]]& out_complex_modulus) except +
        void calc_complex_modulus_vectorize_all(
            const vector[double]& modulus_pa,
            const vector[double]& viscosity_pas,
            const vector[double]& frequency_rad_s,
            vector[cpp_complex[double]]& out_complex_modulus) except +


cdef extern from "rheology_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_RheologyConfig:
        double alpha
        double zeta
        double voigt_modulus_frac
        double voigt_viscosity_frac

    cdef cppclass c_Elastic(c_RheologyBase):
        c_Elastic() except +
        c_Elastic(const c_RheologyConfig& cfg) except +

    cdef cppclass c_Viscous(c_RheologyBase):
        c_Viscous() except +
        c_Viscous(const c_RheologyConfig& cfg) except +

    cdef cppclass c_Maxwell(c_RheologyBase):
        c_Maxwell() except +
        c_Maxwell(const c_RheologyConfig& cfg) except +

    cdef cppclass c_Voigt(c_RheologyBase):
        c_Voigt() except +
        c_Voigt(const c_RheologyConfig& cfg) except +
        double get_voigt_modulus_frac()   const
        double get_voigt_viscosity_frac() const

    cdef cppclass c_Burgers(c_RheologyBase):
        c_Burgers() except +
        c_Burgers(const c_RheologyConfig& cfg) except +
        double get_voigt_modulus_frac()   const
        double get_voigt_viscosity_frac() const

    cdef cppclass c_Andrade(c_RheologyBase):
        c_Andrade() except +
        c_Andrade(const c_RheologyConfig& cfg) except +
        double get_alpha() const
        double get_zeta()  const

    cdef cppclass c_Sundberg(c_RheologyBase):
        c_Sundberg() except +
        c_Sundberg(const c_RheologyConfig& cfg) except +
        double get_alpha()                const
        double get_zeta()                 const
        double get_voigt_modulus_frac()   const
        double get_voigt_viscosity_frac() const

    # Enum naming each concrete rheology model.
    cdef enum class c_RheologyModel:
        Elastic
        Viscous
        Voigt
        Maxwell
        Burgers
        Andrade
        Sundberg

    # Map a name/alias to the enum (raises ValueError on unknown name).
    c_RheologyModel c_rheology_model_from_name(const string& model_name) except +

    # Build the model named by the enum; returns an owning unique_ptr.
    unique_ptr[c_RheologyBase] c_find_rheology(
        c_RheologyModel model, const c_RheologyConfig& cfg) except +


# =====================================================================================================================
# Cython wrapper class declarations
# =====================================================================================================================
cdef class RheologyBase(PhysicsBase):
    cdef unique_ptr[c_RheologyBase] _rheology_ptr   # owns the most-derived C++ model object
    cpdef dict get_config_dict(self)


cdef class Elastic(RheologyBase):
    pass


cdef class Viscous(RheologyBase):
    pass


cdef class Maxwell(RheologyBase):
    pass


cdef class Voigt(RheologyBase):
    cdef c_Voigt* _voigt_ptr   # non-owning; ownership via RheologyBase._rheology_ptr
    cpdef dict get_config_dict(self)


cdef class Burgers(RheologyBase):
    cdef c_Burgers* _burgers_ptr   # non-owning; ownership via RheologyBase._rheology_ptr
    cpdef dict get_config_dict(self)


cdef class Andrade(RheologyBase):
    cdef c_Andrade* _andrade_ptr   # non-owning; ownership via RheologyBase._rheology_ptr
    cpdef dict get_config_dict(self)


cdef class Sundberg(RheologyBase):
    cdef c_Sundberg* _sundberg_ptr   # non-owning; ownership via RheologyBase._rheology_ptr
    cpdef dict get_config_dict(self)
