from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr
from libcpp cimport bool as cpp_bool

from CyRK cimport PreEvalFunc, ODEMethod

from TidalPy.Material_x.eos.ode cimport c_EOS_ODEInput
from TidalPy.Material_x.eos.eos_solution cimport c_EOSSolution


cdef extern from "solver_.hpp" nogil:

    void c_solve_eos(
        c_EOSSolution* eos_solution_ptr,
        vector[PreEvalFunc] eos_function_bylayer_ptr_vec,
        vector[c_EOS_ODEInput] eos_input_bylayer_vec,
        double planet_bulk_density,
        double surface_pressure,
        double G_to_use,
        ODEMethod integration_method,
        double rtol,
        double atol,
        double pressure_tol,
        size_t max_iters,
        cpp_bool verbose
        ) noexcept
