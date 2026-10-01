# distutils: language = c++

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector
from libcpp.complex cimport complex as cpp_complex

from TidalPy.Utilities.classes.classes cimport PhysicsBase, c_PhysicsBase, c_ParamMap


cdef extern from "rheology_base_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_RheologyBase(c_PhysicsBase):
        cpp_complex[double] calc_complex_modulus(
            double modulus, double viscosity, double frequency) const
        void calc_complex_modulus_vectorize(
            const vector[double]& modulus,
            const vector[double]& viscosity,
            const vector[double]& frequency,
            vector[cpp_complex[double]]& out_complex_modulus) except +


cdef extern from "rheology_.hpp" namespace "tidalpy" nogil:

    unique_ptr[c_RheologyBase] c_find_rheology(const string& model_name, const c_ParamMap& params) except +
    string c_rheology_canonical_name(const string& model_name) except +
    vector[string] c_rheology_model_names() except +

    # A copy of a model as its family type, for the layer setters that take ownership.
    unique_ptr[c_RheologyBase] c_clone_rheology "tidalpy::c_clone_as<tidalpy::c_RheologyBase>"(
        const c_RheologyBase& model) except +


cdef class RheologyBase(PhysicsBase):
    cdef c_RheologyBase* _rheology(self) except NULL
