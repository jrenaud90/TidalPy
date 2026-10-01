# distutils: language = c++

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string
from libcpp.memory cimport shared_ptr, unique_ptr
from libcpp.vector cimport vector

from TidalPy.Utilities.classes.classes cimport PhysicsBase, c_PhysicsBase, c_ParamMap, c_ThermoPoint


cdef extern from "material_.hpp" namespace "tidalpy" nogil:

    cdef cppclass c_MaterialSwitches:
        cpp_bool use_thermal_expansion
        cpp_bool use_melting
        cpp_bool use_pressure_melting
        cpp_bool use_melt_density

    cdef enum class c_MaterialPhase:
        Solid
        Partial
        Liquid

    cdef cppclass c_PhaseState:
        double density
        double bulk_modulus
        double adiabatic_bulk_modulus
        double thermal_expansion
        double heat_capacity
        double thermal_conductivity
        double shear_modulus
        double shear_viscosity
        double bulk_viscosity

    cdef cppclass c_MaterialState:
        c_MaterialPhase phase
        double density
        double bulk_modulus
        double adiabatic_bulk_modulus
        double thermal_expansion
        double heat_capacity
        double thermal_conductivity
        double shear_modulus
        double shear_viscosity
        double bulk_viscosity
        double melt_fraction
        double solidus
        double liquidus

    cdef cppclass c_PhaseComponents:
        void set(const string& slot, const shared_ptr[c_PhysicsBase]& model) except +
        shared_ptr[c_PhysicsBase] get(const string& slot) except +

    cdef cppclass c_Phase(c_PhysicsBase):
        const c_PhaseComponents& get_components() const
        void calc_phase_state(const c_ThermoPoint& point, cpp_bool thermal, c_PhaseState& out) const

    cdef cppclass c_MaterialComponents:
        void set(const string& slot, const shared_ptr[c_PhysicsBase]& model) except +
        shared_ptr[c_PhysicsBase] get(const string& slot) except +

    cdef cppclass c_Material(c_PhysicsBase):
        const c_MaterialComponents& get_components() const
        cpp_bool get_can_melt() const
        void calc_melting_range(
            double pressure, const c_MaterialSwitches& switches, double& solidus, double& liquidus) const
        void calc_state(const c_ThermoPoint& point, const c_MaterialSwitches& switches, c_MaterialState& out) const
        double calc_rigidity_margin(
            const c_ThermoPoint& point, const c_MaterialSwitches& switches, double minimum_shear_modulus) const
        void calc_state_vectorize(
            const vector[double]& pressure,
            const vector[double]& temperature,
            const vector[double]& radius,
            const c_MaterialSwitches& switches,
            vector[c_MaterialState]& out_states) except +

    unique_ptr[c_Phase] c_make_phase(const c_ParamMap& params, const c_PhaseComponents& components) except +
    unique_ptr[c_Material] c_make_material(const c_ParamMap& params, const c_MaterialComponents& components) except +


cdef class Phase(PhysicsBase):
    cdef c_Phase* _phase(self) except NULL


cdef class Material(PhysicsBase):
    cdef c_Material* _material(self) except NULL
