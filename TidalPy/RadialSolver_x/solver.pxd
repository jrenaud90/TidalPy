from libcpp cimport bool as cpp_bool
from libcpp.string cimport string as cpp_string
from libcpp.vector cimport vector
from libcpp.complex cimport complex as cpp_complex

from CyRK cimport ODEMethod
from TidalPy.RadialSolver_x.rs_solution cimport c_RadialSolutionStorage

cdef extern from "solver_.hpp":
    void c_validate_and_prep_radial_inputs(
        size_t total_slices,
        const double* radius_array,
        const double* density_array,
        double frequency,
        size_t num_layers,
        const vector[cpp_string]& layer_types,
        const cpp_bool* is_static_bylayer,
        const cpp_bool* is_incompressible_bylayer,
        const double* upper_radius_bylayer_array,
        cpp_bool use_prop_matrix,
        double starting_radius,
        const vector[cpp_string]& solve_for,
        const cpp_string& integration_method,
        const vector[cpp_string]& eos_method_bylayer,
        const cpp_string& eos_integration_method,
        cpp_bool warnings,
        int* layer_types_out_ptr,
        int* bc_models_out_ptr,
        size_t& num_bc_models_out,
        ODEMethod& integration_method_out,
        vector[int]& eos_integration_method_int_bylayer_out,
        ODEMethod& eos_integration_method_out
    ) except +
