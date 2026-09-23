# distutils: language = c++

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector
from libcpp.complex cimport complex as cpp_complex

from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_PhysicsBase


cdef extern from "rheology_base_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_RheologyBase(c_PhysicsBase):
        cpp_complex[double] calc_complex_modulus(
            double modulus, double viscosity, double frequency) const
        void calc_complex_modulus_vectorize_modulus(
            const vector[double]& modulus,
            const vector[double]& viscosity,
            double frequency,
            vector[cpp_complex[double]]& out_complex_modulus) except +
        void calc_complex_modulus_vectorize_frequency(
            double modulus,
            double viscosity,
            const vector[double]& frequency,
            vector[cpp_complex[double]]& out_complex_modulus) except +
        void calc_complex_modulus_vectorize_all(
            const vector[double]& modulus,
            const vector[double]& viscosity,
            const vector[double]& frequency,
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

    cdef enum class c_RheologyModel:
        Elastic
        Viscous
        Voigt
        Maxwell
        Burgers
        Andrade
        Sundberg

    # Raises ValueError on an unknown name.
    c_RheologyModel c_rheology_model_from_name(const string& model_name) except +

    unique_ptr[c_RheologyBase] c_find_rheology(
        c_RheologyModel model, const c_RheologyConfig& cfg) except +


cdef class RheologyBase(PhysicsBase):
    cdef unique_ptr[c_RheologyBase] _rheology_ptr   # owns the most-derived C++ model object


cdef class Elastic(RheologyBase):
    pass


cdef class Viscous(RheologyBase):
    pass


cdef class Maxwell(RheologyBase):
    pass


cdef class Voigt(RheologyBase):
    cdef c_Voigt* _voigt_ptr   # non-owning; ownership via RheologyBase._rheology_ptr


cdef class Burgers(RheologyBase):
    cdef c_Burgers* _burgers_ptr   # non-owning; ownership via RheologyBase._rheology_ptr


cdef class Andrade(RheologyBase):
    cdef c_Andrade* _andrade_ptr   # non-owning; ownership via RheologyBase._rheology_ptr


cdef class Sundberg(RheologyBase):
    cdef c_Sundberg* _sundberg_ptr   # non-owning; ownership via RheologyBase._rheology_ptr
