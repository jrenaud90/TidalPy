# distutils: language = c++
"""
cooling.pxd
Cython declarations for TidalPy's cooling model hierarchy.

Exports the three C++ cooling models, the input/result/config structs, and the
Python wrapper classes so other extensions can cimport and build/attach cooling
models at C speed.

Usage::

    from TidalPy.cooling_x.cooling cimport (
        CoolingBase, ConvectiveCooling,
        c_CoolingBase, c_ConvectiveCooling, c_CoolingInputs, c_CoolingResult)
"""

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector

from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_PhysicsBase


# =====================================================================================================================
# C++ class declarations
# =====================================================================================================================
cdef extern from "cooling_base_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_CoolingInputs:
        c_CoolingInputs() except +
        double delta_temp_k
        double thickness_m
        double gravity_m_s2
        double density_kg_m3
        double viscosity_pas
        double thermal_conductivity_w_mk
        double thermal_diffusivity_m2_s
        double thermal_expansion_1_k

    cdef cppclass c_CoolingResult:
        c_CoolingResult() except +
        double cooling_flux_w_m2
        double blt_m
        double rayleigh_number
        double nusselt_number

    cdef cppclass c_CoolingBase(c_PhysicsBase):
        c_CoolingResult calc_cooling(const c_CoolingInputs& inputs) const
        void calc_cooling_vectorize_temperature(
            const vector[double]& delta_temp_k,
            const c_CoolingInputs& base_inputs,
            vector[c_CoolingResult]& out_results) except +
        void calc_cooling_vectorize_viscosity(
            const vector[double]& viscosity_pas,
            const c_CoolingInputs& base_inputs,
            vector[c_CoolingResult]& out_results) except +
        void calc_cooling_vectorize_all(
            const vector[double]& delta_temp_k,
            const vector[double]& viscosity_pas,
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

    # Enum naming each cooling model.
    cdef enum class c_CoolingModel:
        Off
        Convection
        Conduction

    # Map a name/alias to the enum (raises ValueError on unknown name).
    c_CoolingModel c_cooling_model_from_name(const string& model_name) except +

    # Build the model named by the enum; returns an owning unique_ptr.
    unique_ptr[c_CoolingBase] c_find_cooling(
        c_CoolingModel model, const c_CoolingConfig& cfg) except +


# =====================================================================================================================
# Cython wrapper class declarations
# =====================================================================================================================
cdef class CoolingResult:
    cdef public object cooling_flux              # [W/m^2]
    cdef public object boundary_layer_thickness  # [m]
    cdef public object rayleigh                  # [dimensionless]
    cdef public object nusselt                   # [dimensionless]


cdef class CoolingBase(PhysicsBase):
    cdef unique_ptr[c_CoolingBase] _cooling_ptr   # owns the most-derived C++ model object
    cpdef dict get_config_dict(self)


cdef class OffCooling(CoolingBase):
    pass


cdef class ConductiveCooling(CoolingBase):
    pass


cdef class ConvectiveCooling(CoolingBase):
    cdef c_ConvectiveCooling* _convective_ptr   # non-owning; ownership via CoolingBase._cooling_ptr
    cpdef dict get_config_dict(self)
