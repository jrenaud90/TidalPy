# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libcpp.memory cimport make_unique
from libcpp.string cimport string as cpp_string

from TidalPy.RadialSolver_x.constants cimport C_MAX_NUM_Y
from TidalPy.constants cimport d_PI

cimport numpy as cnp
import numpy as np
cnp.import_array()


cdef class RadialSolverSolution:

    def __init__(
            self,
            size_t num_ytypes,
            double[::1] upper_radius_bylayer_view,
            double[::1] radius_array_view,
            int degree_l
            ):
        # Build pointers
        cdef double* upper_radius_bylayer_ptr = &upper_radius_bylayer_view[0]
        self.num_layers                       = upper_radius_bylayer_view.size
        cdef double* radius_array_ptr         = &radius_array_view[0]
        cdef size_t radius_array_size         = radius_array_view.size

        # Set state information
        self.ytype_names_set = False
        self.num_ytypes      = num_ytypes

        # Create C++ storage instance
        self.solution_storage_uptr = make_unique[c_RadialSolutionStorage](
            self.num_ytypes,
            upper_radius_bylayer_ptr,
            self.num_layers,
            radius_array_ptr,
            radius_array_size,
            degree_l)
        self.solution_storage_ptr = self.solution_storage_uptr.get()

        if not self.solution_storage_ptr:
            raise RuntimeError("c_RadialSolutionStorage extension class could not be initialized.")

        # Finish initialization with provided array
        self.change_radius_array(radius_array_ptr, radius_array_size, array_changed=False)

    def __dealloc__(self):
        self.solution_storage_uptr.reset()
        self.solution_storage_ptr = NULL

    cdef void set_model_names(self, int* bc_models_ptr) noexcept nogil:
        cdef size_t ytype_i
        cdef int bc_model
        for ytype_i in range(self.num_ytypes):
            bc_model = bc_models_ptr[ytype_i]
            if bc_model == 0:
                self.ytypes[ytype_i] = "free"
            elif bc_model == 1:
                self.ytypes[ytype_i] = "tidal"
            elif bc_model == 2:
                self.ytypes[ytype_i] = "loading"
            else:
                self.solution_storage_ptr.error_code = -2
                self.solution_storage_ptr.message = cpp_string(b"ArgumentException:: Unknown boundary condition provided")
        self.ytype_names_set = True

    cdef void change_radius_array(
            self,
            double* new_radius_array_ptr,
            size_t new_size_radius_array,
            cpp_bool array_changed = True) noexcept:

        self.radius_array_size = new_size_radius_array

        if array_changed:
            self.solution_storage_ptr.change_radius_array(new_radius_array_ptr, new_size_radius_array, array_changed)

        cdef cnp.npy_intp[2] full_solution_shape   = [self.radius_array_size, self.num_ytypes * C_MAX_NUM_Y]
        cdef cnp.npy_intp* full_solution_shape_ptr = &full_solution_shape[0]
        cdef cnp.npy_intp full_solution_shape_ndim = 2

        cdef cnp.npy_intp[1] love_shape   = [self.num_ytypes * 3]
        cdef cnp.npy_intp* love_shape_ptr = &love_shape[0]
        cdef cnp.npy_intp love_shape_ndim = 1

        cdef cnp.npy_intp[1] eos_float_shape     = [self.radius_array_size]
        cdef cnp.npy_intp* eos_float_shape_ptr   = &eos_float_shape[0]
        cdef cnp.npy_intp[1] eos_complex_shape   = [self.radius_array_size]
        cdef cnp.npy_intp* eos_complex_shape_ptr = &eos_complex_shape[0]
        cdef cnp.npy_intp eos_ndim               = 1

        cdef c_EOSSolution* eos_solution_ptr = self.solution_storage_ptr.get_eos_solution_ptr()

        if not self.solution_storage_ptr:
            raise RuntimeError("RadialSolverSolution:: c_RadialSolutionStorage is not initialized.")
        else:
            self.full_solution_arr = cnp.PyArray_SimpleNewFromData(
                full_solution_shape_ndim,
                full_solution_shape_ptr,
                cnp.NPY_COMPLEX128,
                <double complex*>self.solution_storage_ptr.full_solution_vec.data())

            self.complex_love_arr = cnp.PyArray_SimpleNewFromData(
                love_shape_ndim,
                love_shape_ptr,
                cnp.NPY_COMPLEX128,
                <double complex*>self.solution_storage_ptr.complex_love_vec.data())

            if not eos_solution_ptr:
                raise RuntimeError("RadialSolverSolution:: c_EOSSolution is not initialized.")
            else:
                self.radius_array_cnp = cnp.PyArray_SimpleNewFromData(
                    eos_ndim, eos_float_shape_ptr, cnp.NPY_FLOAT64,
                    eos_solution_ptr.radius_array_vec.data())

                self.gravity_array_cnp = cnp.PyArray_SimpleNewFromData(
                    eos_ndim, eos_float_shape_ptr, cnp.NPY_FLOAT64,
                    eos_solution_ptr.gravity_array_vec.data())

                self.pressure_array_cnp = cnp.PyArray_SimpleNewFromData(
                    eos_ndim, eos_float_shape_ptr, cnp.NPY_FLOAT64,
                    eos_solution_ptr.pressure_array_vec.data())

                self.mass_array_cnp = cnp.PyArray_SimpleNewFromData(
                    eos_ndim, eos_float_shape_ptr, cnp.NPY_FLOAT64,
                    eos_solution_ptr.mass_array_vec.data())

                self.moi_array_cnp = cnp.PyArray_SimpleNewFromData(
                    eos_ndim, eos_float_shape_ptr, cnp.NPY_FLOAT64,
                    eos_solution_ptr.moi_array_vec.data())

                self.density_array_cnp = cnp.PyArray_SimpleNewFromData(
                    eos_ndim, eos_float_shape_ptr, cnp.NPY_FLOAT64,
                    eos_solution_ptr.density_array_vec.data())

                self.shear_modulus_array_cnp = cnp.PyArray_SimpleNewFromData(
                    eos_ndim, eos_complex_shape_ptr, cnp.NPY_COMPLEX128,
                    <double complex*>eos_solution_ptr.complex_shear_array_vec.data())

                self.bulk_modulus_array_cnp = cnp.PyArray_SimpleNewFromData(
                    eos_ndim, eos_complex_shape_ptr, cnp.NPY_COMPLEX128,
                    <double complex*>eos_solution_ptr.complex_bulk_array_vec.data())

    cdef void finalize_python_storage(self) noexcept:

        cdef cnp.npy_intp[2] eos_steps_taken_shape   = [self.solution_storage_ptr.get_eos_solution_ptr().num_cyolver_calls / self.num_layers, self.num_layers]
        cdef cnp.npy_intp* eos_steps_taken_shape_ptr = &eos_steps_taken_shape[0]
        cdef cnp.npy_intp eos_steps_taken_ndim       = 2
        self.eos_steps_taken_array = cnp.PyArray_SimpleNewFromData(
            eos_steps_taken_ndim,
            eos_steps_taken_shape_ptr,
            cnp.NPY_UINT64,
            self.solution_storage_ptr.get_eos_solution_ptr().steps_taken_vec.data())

        cdef cnp.npy_intp[2] steps_taken_shape   = [self.num_layers, 3]
        cdef cnp.npy_intp* steps_taken_shape_ptr = &steps_taken_shape[0]
        cdef cnp.npy_intp steps_taken_ndim       = 2
        self.shooting_method_steps_taken_array = cnp.PyArray_SimpleNewFromData(
            steps_taken_ndim,
            steps_taken_shape_ptr,
            cnp.NPY_UINT64,
            self.solution_storage_ptr.shooting_method_steps_taken_vec.data())

    # Properties
    @property
    def error_code(self):
        return self.solution_storage_ptr.error_code

    @property
    def message(self):
        return self.solution_storage_ptr.message.decode('utf-8')

    @property
    def success(self):
        return self.solution_storage_ptr.success

    @property
    def eos_error_code(self):
        return self.solution_storage_ptr.get_eos_solution_ptr().error_code

    @property
    def eos_message(self):
        return self.solution_storage_ptr.get_eos_solution_ptr().message.decode('utf-8')

    @property
    def eos_success(self):
        return self.solution_storage_ptr.get_eos_solution_ptr().success

    @property
    def eos_pressure_error(self):
        return self.solution_storage_ptr.get_eos_solution_ptr().pressure_error

    @property
    def eos_iterations(self):
        return self.solution_storage_ptr.get_eos_solution_ptr().iterations

    @property
    def eos_steps_taken(self):
        return np.copy(self.eos_steps_taken_array)

    @property
    def radius_array(self):
        return np.copy(self.radius_array_cnp)

    @property
    def gravity_array(self):
        return np.copy(self.gravity_array_cnp)

    @property
    def pressure_array(self):
        return np.copy(self.pressure_array_cnp)

    @property
    def mass_array(self):
        return np.copy(self.mass_array_cnp)

    @property
    def moi_array(self):
        return np.copy(self.moi_array_cnp)

    @property
    def density_array(self):
        return np.copy(self.density_array_cnp)

    @property
    def shear_modulus_array(self):
        return np.copy(self.shear_modulus_array_cnp)

    @property
    def bulk_modulus_array(self):
        return np.copy(self.bulk_modulus_array_cnp)

    @property
    def radius(self):
        return self.solution_storage_ptr.get_eos_solution_ptr().radius

    @property
    def volume(self):
        return (4.0 / 3.0) * d_PI * self.radius**3

    @property
    def mass(self):
        return self.solution_storage_ptr.get_eos_solution_ptr().mass

    @property
    def moi(self):
        return self.solution_storage_ptr.get_eos_solution_ptr().moi

    @property
    def moi_factor(self):
        cdef double ideal_moi = (2.0 / 5.0) * self.mass * self.radius**2
        return self.moi / ideal_moi

    @property
    def density_bulk(self):
        return self.mass / self.volume

    @property
    def central_pressure(self):
        return self.solution_storage_ptr.get_eos_solution_ptr().central_pressure

    @property
    def surface_pressure(self):
        return self.solution_storage_ptr.get_eos_solution_ptr().surface_pressure

    @property
    def surface_gravity(self):
        return self.solution_storage_ptr.get_eos_solution_ptr().surface_gravity

    @property
    def layer_upper_radius_array(self):
        cdef cnp.ndarray[cnp.float64_t, ndim=1] upper_radius_array = np.empty(self.num_layers, dtype=np.float64)
        cdef c_EOSSolution* eos_solution_ptr = self.solution_storage_ptr.get_eos_solution_ptr()
        cdef size_t layer_i
        for layer_i in range(self.num_layers):
            upper_radius_array[layer_i] = eos_solution_ptr.upper_radius_bylayer_vec[layer_i]
        return upper_radius_array

    @property
    def degree_l(self):
        return self.solution_storage_ptr.degree_l

    @property
    def result(self):
        if self.success and (self.error_code == 0):
            return np.copy(self.full_solution_arr).T
        else:
            return None

    @property
    def love(self):
        if self.success and (self.error_code == 0):
            return np.copy(self.complex_love_arr).reshape(self.num_ytypes, 3)
        else:
            return np.nan * np.ones((self.num_ytypes, 3), dtype=np.complex128)

    @property
    def k(self):
        if self.success and (self.error_code == 0):
            if self.num_ytypes == 1:
                return self.complex_love_arr[0]
            else:
                return np.copy(self.complex_love_arr[0::3])
        else:
            if self.num_ytypes == 1:
                return np.nan
            else:
                return np.nan * np.ones(self.num_ytypes, dtype=np.float64)

    @property
    def h(self):
        if self.success and (self.error_code == 0):
            if self.num_ytypes == 1:
                return self.complex_love_arr[1]
            else:
                return np.copy(self.complex_love_arr[1::3])
        else:
            if self.num_ytypes == 1:
                return np.nan
            else:
                return np.nan * np.ones(self.num_ytypes, dtype=np.float64)

    @property
    def l(self):
        if self.success and (self.error_code == 0):
            if self.num_ytypes == 1:
                return self.complex_love_arr[2]
            else:
                return np.copy(self.complex_love_arr[2::3])
        else:
            if self.num_ytypes == 1:
                return np.nan
            else:
                return np.nan * np.ones(self.num_ytypes, dtype=np.float64)

    @property
    def Q(self):
        if self.success and (self.error_code == 0):
            return np.abs(self.k) / -np.imag(self.k)
        else:
            if self.num_ytypes == 1:
                return np.nan
            else:
                return np.nan * np.ones(self.num_ytypes, dtype=np.float64)

    @property
    def lag(self):
        if self.success and (self.error_code == 0):
            return np.arctan(-np.imag(self.k) / np.real(self.k))
        else:
            if self.num_ytypes == 1:
                return np.nan
            else:
                return np.nan * np.ones(self.num_ytypes, dtype=np.float64)

    @property
    def steps_taken(self):
        return np.copy(self.shooting_method_steps_taken_array)

    def get_result_by_ytype_name(self, str ytype_name):
        cdef size_t ytype_i
        cdef size_t requested_sol_num = 0
        cdef cpp_bool found = False
        cdef str sol_test_name
        if self.ytype_names_set and self.success and (self.error_code == 0):
            for ytype_i in range(self.num_ytypes):
                sol_test_name = str(self.ytypes[ytype_i], 'UTF-8')
                if sol_test_name == ytype_name.lower():
                    requested_sol_num = ytype_i
                    found = True
                    break
            if not found:
                raise ValueError('Unknown solution type requested.')

            return np.copy(self.result[C_MAX_NUM_Y * (requested_sol_num): C_MAX_NUM_Y * (requested_sol_num + 1)])
        else:
            return None

    def __len__(self):
        return <Py_ssize_t>self.num_ytypes

    def __getitem__(self, str ytype_name):
        return self.get_result_by_ytype_name(ytype_name)
