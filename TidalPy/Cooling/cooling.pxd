# distutils: language = c++

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector

from TidalPy.Utilities.classes.classes cimport PhysicsBase, c_PhysicsBase, c_ParamMap


cdef extern from "cooling_base_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_CoolingInputs:
        c_CoolingInputs() except +
        double delta_temp
        double thickness
        double gravity
        double density
        double viscosity
        double thermal_conductivity
        double thermal_diffusivity
        double thermal_expansion
        cpp_bool liquid

    cdef cppclass c_CoolingResult:
        c_CoolingResult() except +
        double cooling_flux
        double blt
        double rayleigh_number
        double nusselt_number

    cdef cppclass c_CoolingBase(c_PhysicsBase):
        c_CoolingResult calc_cooling(const c_CoolingInputs& inputs) except +
        void calc_cooling_vectorize(
            const vector[double]& delta_temp,
            const vector[double]& viscosity,
            const c_CoolingInputs& base_inputs,
            size_t num_points,
            double* out_cooling_flux,
            double* out_blt,
            double* out_rayleigh,
            double* out_nusselt) except +


cdef extern from "cooling_.hpp" namespace "tidalpy" nogil:

    unique_ptr[c_CoolingBase] c_find_cooling(const string& model_name, const c_ParamMap& params) except +
    string c_cooling_canonical_name(const string& model_name) except +
    vector[string] c_cooling_model_names() except +


cdef class CoolingResult:
    cdef public object cooling_flux              # [W/m^2]
    cdef public object boundary_layer_thickness  # [m]
    cdef public object rayleigh                  # [dimensionless]
    cdef public object nusselt                   # [dimensionless]


cdef class CoolingBase(PhysicsBase):
    cdef c_CoolingBase* _cooling(self) except NULL
