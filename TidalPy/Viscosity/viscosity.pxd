# distutils: language = c++

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector

from TidalPy.Utilities.classes.classes cimport PhysicsBase, c_PhysicsBase, c_ParamMap, c_ThermoPoint


cdef extern from "viscosity_base_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_ViscosityBase(c_PhysicsBase):
        double calc_viscosity(const c_ThermoPoint& point) const
        void calc_viscosity_vectorize(
            const vector[double]& temperature,
            const vector[double]& pressure,
            const vector[double]& radius,
            vector[double]& out_viscosity) except +


cdef extern from "viscosity_.hpp" namespace "tidalpy" nogil:

    unique_ptr[c_ViscosityBase] c_find_viscosity(const string& model_name, const c_ParamMap& params) except +
    string c_viscosity_canonical_name(const string& model_name) except +
    vector[string] c_viscosity_model_names() except +

    # A copy of a model as its family type, for the layer and material setters that take ownership.
    unique_ptr[c_ViscosityBase] c_clone_viscosity "tidalpy::c_clone_as<tidalpy::c_ViscosityBase>"(
        const c_ViscosityBase& model) except +


cdef class ViscosityBase(PhysicsBase):
    cdef c_ViscosityBase* _viscosity(self) except NULL
