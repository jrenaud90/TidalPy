cimport numpy as cnp
cnp.import_array()

from libcpp cimport bool as cpp_bool
from libcpp.string cimport string as cpp_string
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr
from libcpp.complex cimport complex as cpp_complex

from TidalPy.Utilities_x.dimensions.nondimensional cimport c_NonDimensionalScales
from TidalPy.Material_x.eos.eos_solution cimport c_EOSSolution
from TidalPy.RadialSolver_x.love cimport c_LoveNumbers

cdef extern from "rs_solution_.hpp" nogil:

    cdef cppclass c_RadialSolutionStorage:

        c_RadialSolutionStorage()
        c_RadialSolutionStorage(
            size_t num_ytypes,
            double* upper_radius_bylayer_ptr,
            size_t num_layers,
            double* radius_array_ptr,
            size_t size_radius_array,
            int degree_l)

        cpp_bool success
        int error_code
        int degree_l
        cpp_string message
        size_t num_ytypes
        size_t num_slices
        size_t num_layers
        size_t total_size
        unique_ptr[c_EOSSolution] eos_solution_uptr
        vector[double] full_solution_vec
        vector[c_LoveNumbers] complex_love_vec
        vector[int] p_bc_models
        vector[size_t] shooting_method_steps_taken_vec
        double surface_amplification
        double p_love_frequency_si

        c_EOSSolution* get_eos_solution_ptr()
        void change_radius_array(
            double* new_radius_array_ptr,
            size_t new_size_radius_array,
            cpp_bool array_changed)
        cpp_bool get_radial_solution(
            double radius_si,
            size_t ytype_i,
            cpp_complex[double]* out6)
        void get_radial_solution_array(
            const double* radii_si,
            size_t n,
            size_t ytype_i,
            cpp_complex[double]* out)
        cpp_bool get_surface_y(size_t ytype_i, cpp_complex[double]* out6)
        cpp_bool get_eos_si(double radius_si, double* out)
        void get_complex_moduli_si(double radius_si, cpp_complex[double]& shear_out, cpp_complex[double]& bulk_out)
        void find_love()
        void dimensionalize_data(
            c_NonDimensionalScales* nondim_scales,
            cpp_bool redimensionalize)


cdef class RadialSolverSolution:

    # Size and state information
    cdef size_t radius_array_size
    cdef public size_t num_ytypes
    cdef public size_t num_layers
    cdef cpp_bool ytype_names_set
    cdef char* ytypes[5]

    # Main storage container
    cdef unique_ptr[c_RadialSolutionStorage] solution_storage_uptr
    cdef c_RadialSolutionStorage* solution_storage_ptr

    # The world this solution was released from, if any.
    cdef object p_source_world

    # Result pointers and data
    cdef cnp.ndarray full_solution_arr

    # EOS solution arrays

    # Shooting method diagnostics
    cdef cnp.ndarray shooting_method_steps_taken_array
    cdef cnp.ndarray eos_steps_taken_array

    cdef void finalize_python_storage(self) noexcept

    # Complex shear (which = 0) or bulk (1) modulus over an array of radii, swept in C.
    cdef object _complex_moduli_sweep(self, object radius, size_t which)

    # Adopt a storage released by a world, instead of building one (see the .pyx).
    @staticmethod
    cdef RadialSolverSolution _adopt(
        unique_ptr[c_RadialSolutionStorage] storage_uptr,
        object source_world)

    cdef void set_model_names(
        self,
        int* bc_models_ptr) noexcept nogil

    cdef void change_radius_array(
        self,
        double* new_radius_array_ptr,
        size_t new_size_radius_array,
        cpp_bool array_changed = *) noexcept


# Warn when the surface boundary condition solve is poorly conditioned. Returns True when it warned. Holds
# the GIL: it formats and logs a message.
cdef bint cy_check_surface_solve_conditioning(
    double surface_amplification,
    double integration_rtol) except *
