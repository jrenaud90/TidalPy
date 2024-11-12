# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libc.stdlib cimport exit, EXIT_FAILURE
from libc.stdio cimport printf

import numpy as np
cimport numpy as np

from TidalPy.RadialSolver.constants cimport MAX_NUM_Y


cdef class RadialSolverSolution:

    def __init__(
            self,
            size_t num_slices,
            char num_ytypes
            ):

        # Set state information
        self.ytype_names_set = False
        self.num_slices = num_slices
        self.num_ytypes = num_ytypes

        # Create C++ instance
        self.solution_storage_ptr = new RadialSolutionStorageCC(self.num_slices, self.num_ytypes)

        # The RadialSolutionStorage class has full control of memory. This class simply wraps it.
        cdef np.npy_intp[2] shape   = [num_slices, num_ytypes * MAX_NUM_Y]
        cdef np.npy_intp* shape_ptr = &shape[0]
        cdef np.npy_intp ndim       = 2

        # `solution_storage_ptr.full_solution_ptr` is a double pointer but it is really storing double complex data
        # ordered by y0_real, y0_imag, y1_real, y1_image, ... so we can safely convert it to a complex128 np.ndarray
        if not (self.solution_storage_ptr.full_solution_ptr is NULL):
            self.full_solution_arr = np.PyArray_SimpleNewFromData(
                ndim,
                shape_ptr,
                np.NPY_COMPLEX128,
                self.solution_storage_ptr.full_solution_ptr)
        
        # Initialize Love numbers (3 love numbers for each requested y-type)
        # Love numbers are stored (k, h, l)_ytype0, (k, h, l)_ytype1, (k, h, l)_ytype2, ...
        # If there is only 1 ytype then return a 1-D array, otherwise return a 2D one where the first index is by y-type
        if num_ytypes == 1:
            shape_ptr[0] = 3
            shape_ptr[1] = 0
            ndim = 1
        else:
            shape_ptr[0] = num_ytypes
            shape_ptr[1] = 3
            ndim = 2

        # Same note as above, `solution_storage_ptr.complex_love_ptr` is a double pointer that we are converting to a
        # complex128 np.ndarray.
        if not (self.solution_storage_ptr.complex_love_ptr is NULL):
            self.complex_love_arr = np.PyArray_SimpleNewFromData(
                ndim,
                shape_ptr,
                np.NPY_COMPLEX128,
                self.solution_storage_ptr.complex_love_ptr)
        
        # Make arrays for all of the equation of state variables in a similar manner to the above.
        cdef np.npy_intp[1] eos_float_shape     = [num_slices]
        cdef np.npy_intp* eos_float_shape_ptr   = &eos_float_shape[0]
        cdef np.npy_intp[1] eos_complex_shape   = [2 * num_slices]
        cdef np.npy_intp* eos_complex_shape_ptr = &eos_complex_shape[0]
        cdef np.npy_intp eos_ndim                = 1

        if not (self.solution_storage_ptr.eos_properties_ptr is NULL):
            if not (self.solution_storage_ptr.gravity_ptr is NULL):
                self.gravity_array = np.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_float_shape_ptr,
                    np.NPY_FLOAT64,
                    self.solution_storage_ptr.gravity_ptr)

            if not (self.solution_storage_ptr.pressure_ptr is NULL):
                self.pressure_array = np.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_float_shape_ptr,
                    np.NPY_FLOAT64,
                    self.solution_storage_ptr.pressure_ptr)
            
            if not (self.solution_storage_ptr.density_ptr is NULL):
                self.density_array = np.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_float_shape_ptr,
                    np.NPY_FLOAT64,
                    self.solution_storage_ptr.density_ptr)
            
            if not (self.solution_storage_ptr.shear_mod_ptr is NULL):
                self.shear_modulus_array = np.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_complex_shape_ptr,
                    np.NPY_COMPLEX128,
                    self.solution_storage_ptr.shear_mod_ptr)
            
            if not (self.solution_storage_ptr.bulk_mod_ptr is NULL):
                self.bulk_modulus_array = np.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_complex_shape_ptr,
                    np.NPY_COMPLEX128,
                    self.solution_storage_ptr.bulk_mod_ptr)

    cdef void set_model_names(self, int* bc_models_ptr) noexcept nogil:
        # Unfortunately this must be done outside of __init__ because the argument is a pointer and python interpretor
        # is used during __init__ and it does not like pointers.

        # Find solution types
        cdef char ytype_i
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
                printf("AttributeError:: Unknown boundary condition")
                exit(EXIT_FAILURE)
        self.ytype_names_set = True

    def __dealloc__(self):

        # Release the heap allocated storage
        del self.solution_storage_ptr

    @property
    def message(self):
        """ Return solver's message """
        return str(self.solution_storage_ptr.message_ptr, 'UTF-8')
    
    @property
    def success(self):
        """ Return if the solver was successful message """
        return self.solution_storage_ptr.success

    @property
    def eos_message(self):
        """ Return solver's equation of state message """
        return str(self.solution_storage_ptr.eos_message_ptr, 'UTF-8')
    
    @property
    def eos_success(self):
        """ Return if the solver's equation of state sub-solver was successful """
        return self.solution_storage_ptr.eos_success

    @property
    def result(self):
        """ Return result array. """

        if self.solution_storage_ptr.success:
            # TODO: Optimize solution storage so that transpose is not required?
            return self.full_solution_arr.T
        else:
            return None

    @property
    def love(self):
        """ Return all complex love numbers. """
        if self.solution_storage_ptr.success:
            return self.complex_love_arr
        else:
            return None

    @property
    def k(self):
        """ Tidal Love number k. """
        if self.success:
            if self.num_ytypes == 1:
                # 1D Love
                return self.complex_love_arr[0]
            else:
                # 2D Slice skipping other Love Numbers.
                return self.complex_love_arr[0::3]
        else:
            return None

    @property
    def h(self):
        """ Tidal Love number h. """
        if self.success:
            if self.num_ytypes == 1:
                # 1D Love
                return self.complex_love_arr[1]
            else:
                # 2D Slice skipping other Love Numbers.
                return self.complex_love_arr[1::3]
        else:
            return None
    
    @property
    def l(self):
        """ Tidal Shida number l. """
        if self.success:
            if self.num_ytypes == 1:
                # 1D Love
                return self.complex_love_arr[2]
            else:
                # 2D Slice skipping other Love Numbers.
                return self.complex_love_arr[2::3]
        else:
            return None

    def __len__(self):
        """Return number of solution types."""
        return <Py_ssize_t>self.num_ytypes
    
    def __getitem__(self, str ytype_name):
        """Get a specific solution type array."""
        
        cdef char ytype_i
        cdef char requested_sol_num = 0
        cdef cpp_bool found = False
        cdef str sol_test_name
        if self.ytype_names_set and self.solution_storage_ptr.success:
            for ytype_i in range(self.num_ytypes):
                sol_test_name = str(self.ytypes[ytype_i], 'UTF-8')
                if sol_test_name == ytype_name:
                    requested_sol_num = ytype_i
                    found = True
                    break
            if not found:
                raise AttributeError('Unknown solution type requested. Key must match names passed to radial_solver "solve_for" argument and be lower case.')
            
            # Slice the result and return only the requested solution type.
            return self.result[MAX_NUM_Y * (requested_sol_num): MAX_NUM_Y * (requested_sol_num + 1)]
        else:
            return None


cdef size_t cf_find_num_shooting_solutions(
        int layer_type,
        bint is_static,
        bint is_incompressible
        ) noexcept nogil:
    """ Determine number of solutions required for layer based on assumptions.
    
    Parameters
    ----------
    layer_type : int
        - 0: Layer is solid 
        - 1: Layer is liquid
    is_static : bool
        Use static (True) or dynamic (False) assumption.
    is_incompressible : bool
        Use incompressible (True) or compressible (False) assumption.

    Returns
    -------
    num_sols : int
        Number of solutions required for layer.

    """

    # Initialize
    cdef Py_ssize_t num_sols
    num_sols = 0

    if (layer_type == 0):
        # Solid
        if is_static:
            if is_incompressible:
                # TODO: Confirm
                num_sols = 3
            else:
                num_sols = 3
        else:
            # Dynamic
            if is_incompressible:
                # TODO: Confirm
                num_sols = 3
            else:
                num_sols = 3
    else:
        # Liquid
        if is_static:
            if is_incompressible:
                # TODO: Confirm
                num_sols = 1
            else:
                num_sols = 1
        else:
            # Dynamic
            if is_incompressible:
                # TODO: Confirm
                num_sols = 2
            else:
                num_sols = 2
    return num_sols


def find_num_shooting_solutions(
        int layer_type,
        bint is_static,
        bint is_incompressible
        ):
    
    return cf_find_num_shooting_solutions(layer_type, is_static, is_incompressible)
