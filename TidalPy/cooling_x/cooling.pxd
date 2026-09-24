# distutils: language = c++

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector

from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_PhysicsBase


cdef extern from "cooling_base_.hpp" namespace "tidalpy" nogil:

    cdef enum class c_CoolingModel:
        Off
        Convection
        Conduction

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

    cdef cppclass c_CoolingResult:
        c_CoolingResult() except +
        double cooling_flux
        double blt
        double rayleigh_number
        double nusselt_number

    cdef cppclass c_CoolingBase(c_PhysicsBase):
        c_CoolingResult calc_cooling(const c_CoolingInputs& inputs) const
        void calc_cooling_vectorize(
            const vector[double]& delta_temp,
            const vector[double]& viscosity,
            const c_CoolingInputs& base_inputs,
            vector[c_CoolingResult]& out_results) except +


cdef extern from "cooling_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_CoolingConfig:
        double convection_alpha
        double convection_beta
        double critical_rayleigh

    cdef cppclass c_OffCooling(c_CoolingBase):
        c_OffCooling() except +
        c_OffCooling(const c_CoolingConfig& cfg) except +

    cdef cppclass c_ConductiveCooling(c_CoolingBase):
        c_ConductiveCooling() except +
        c_ConductiveCooling(const c_CoolingConfig& cfg) except +

    cdef cppclass c_ConvectiveCooling(c_CoolingBase):
        c_ConvectiveCooling() except +
        c_ConvectiveCooling(const c_CoolingConfig& cfg) except +
        double get_convection_alpha()  const
        double get_convection_beta()   const
        double get_critical_rayleigh() const


    # Raises ValueError on an unknown name.
    c_CoolingModel c_cooling_model_from_name(const string& model_name) except +

    unique_ptr[c_CoolingBase] c_find_cooling(
        c_CoolingModel model, const c_CoolingConfig& cfg) except +


cdef class CoolingResult:
    cdef public object cooling_flux              # [W/m^2]
    cdef public object boundary_layer_thickness  # [m]
    cdef public object rayleigh                  # [dimensionless]
    cdef public object nusselt                   # [dimensionless]


cdef class CoolingBase(PhysicsBase):
    cdef unique_ptr[c_CoolingBase] _cooling_ptr   # owns the most-derived C++ model object
    cdef void _adopt(self, unique_ptr[c_CoolingBase]& model) noexcept


cdef class OffCooling(CoolingBase):
    pass


cdef class ConductiveCooling(CoolingBase):
    pass


cdef class ConvectiveCooling(CoolingBase):
    pass
