# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libc.stdlib cimport exit, EXIT_FAILURE
from libc.stdio cimport printf
from libcpp.memory cimport make_shared

import numpy as np
np.import_array()

from TidalPy.RadialSolver.constants cimport MAX_NUM_Y


cdef class RadialSolverSolution:

    def __init__(
            self,
            char num_ytypes,
            double* upper_radius_bylayer_ptr,
            const size_t num_layers,
            double* radius_array_ptr,
            const size_t radius_array_size
            ):

        printf("DEBUG-\t RadialSolverSolution Point 1\n")
        # Set state information
        self.ytype_names_set = False
        self.num_ytypes      = num_ytypes

        # Create C++ storage instance
        printf("DEBUG-\t RadialSolverSolution Point 2\n")
        self.solution_storage_sptr = make_shared[RadialSolutionStorageCC](
            num_ytypes,
            upper_radius_bylayer_ptr,
            num_layers,
            radius_array_ptr,
            radius_array_size)

        if not self.solution_storage_sptr.get():
            # C++ solution storage could not be initialized.
            raise RuntimeError("RadialSolutionStorageCC extension class could not be initialized.")
        
        # Finish initialization with provided array
        self.change_radius_array(radius_array_ptr, radius_array_size, array_changed=False)

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
                self.solution_storage_sptr.get().error_code = -2
                self.solution_storage_sptr.get().set_message("AttributeError:: Unknown boundary condition provided")
        self.ytype_names_set = True

    cdef change_radius_array(
            self,
            double* radius_array_ptr,
            const size_t radius_array_size,
            cpp_bool array_changed=True) noexcept:
        """ Change the radius array used by RadialSolver and Equation of State solvers. """

        # Update attributes in this class
        self.num_slices = radius_array_size

        if array_changed:
            # Radius array held by the underlying c++ classes is no longer correct.
            # Update it using its methods.
            self.solution_storage_sptr.get().change_radius_array(radius_array_ptr, radius_array_size, array_changed)
        
        # Now that the underlying storage has been updated we can build our numpy array wrappers.
    
        # The RadialSolutionStorageCC class has full control of memory. This class simply wraps it.
        # The shape needs to be twice what we expect because the underlying C++ class only works with double arrays 
        # but at this level we want arrays to be double complex (the 2x factor is built into MAX_NUM_Y_REAL)
        printf("DEBUG-\t RadialSolverSolution Point 3\n")
        cdef np.npy_intp[2] full_solution_shape   = [self.num_slices, self.num_ytypes * MAX_NUM_Y_REAL]
        cdef np.npy_intp* full_solution_shape_ptr = &full_solution_shape[0]
        cdef np.npy_intp full_solution_shape_ndim = 2

        # Initialize Love numbers (3 love numbers for each requested y-type)
        # Love numbers are stored (k, h, l)_ytype0, (k, h, l)_ytype1, (k, h, l)_ytype2, ...
        cdef np.npy_intp[2] love_shape   = [3, 0]
        cdef np.npy_intp* love_shape_ptr = &love_shape[0]
        cdef np.npy_intp love_shape_ndim = 1
        # If there is only 1 ytype then return a 1-D array, otherwise return a 2D one where the first index is by y-type
        if self.num_ytypes == 1:
            love_shape_ptr[0] = 3
            love_shape_ptr[1] = 0
            love_shape_ndim   = 1
        else:
            love_shape_ptr[0] = self.num_ytypes
            love_shape_ptr[1] = 3
            love_shape_ndim   = 2
        
        # Make numpy arrays that wrap all of the equation of state class vectors in a similar manner to the above.
        cdef np.npy_intp[1] eos_float_shape     = [self.num_slices]
        cdef np.npy_intp* eos_float_shape_ptr   = &eos_float_shape[0]
        cdef np.npy_intp[1] eos_complex_shape   = [2 * self.num_slices]
        cdef np.npy_intp* eos_complex_shape_ptr = &eos_complex_shape[0]
        cdef np.npy_intp eos_ndim               = 1

        # `solution_storage_sptr.full_solution_vec` is a double vector but it is really storing double complex data
        # ordered by y0_real, y0_imag, y1_real, y1_image, ... so we can safely convert it to a complex128 np.ndarray
        printf("DEBUG-\t RadialSolverSolution Point 4\n")
        if not self.solution_storage_sptr.get():
            raise RuntimeError("RadialSolutionStorageCC extension class could not be initialized.")
        else:
            printf("DEBUG-\t RadialSolverSolution Point 4a\n")
            self.full_solution_arr = np.PyArray_SimpleNewFromData(
                full_solution_shape_ndim,
                full_solution_shape_ptr,
                np.NPY_COMPLEX128,
                &self.solution_storage_sptr.get().full_solution_vec[0])

            # Same note as above, `solution_storage_sptr.complex_love_vec` is a double vector that we are converting to a
            # complex128 np.ndarray.
            printf("DEBUG-\t RadialSolverSolution Point 4b\n")
            self.complex_love_arr = np.PyArray_SimpleNewFromData(
                love_shape_ndim,
                love_shape_ptr,
                np.NPY_COMPLEX128,
                &self.solution_storage_sptr.get().complex_love_vec[0])

            if not self.solution_storage_sptr.get().eos_solution_sptr.get():
                raise RuntimeError("EOSSolutionCC extension class could not be initialized.")
            else:
                printf("DEBUG-\t RadialSolverSolution Point 4c\n")
                self.gravity_array = np.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_float_shape_ptr,
                    np.NPY_FLOAT64,
                    &self.solution_storage_sptr.get().eos_solution_sptr.get().gravity_array_vec[0])

                self.pressure_array = np.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_float_shape_ptr,
                    np.NPY_FLOAT64,
                    &self.solution_storage_sptr.get().eos_solution_sptr.get().pressure_array_vec[0])
                
                self.mass_array = np.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_float_shape_ptr,
                    np.NPY_FLOAT64,
                    &self.solution_storage_sptr.get().eos_solution_sptr.get().mass_array_vec[0])
                
                self.moi_array = np.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_float_shape_ptr,
                    np.NPY_FLOAT64,
                    &self.solution_storage_sptr.get().eos_solution_sptr.get().moi_array_vec[0])

                self.density_array = np.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_float_shape_ptr,
                    np.NPY_FLOAT64,
                    &self.solution_storage_sptr.get().eos_solution_sptr.get().density_array_vec[0])

                printf("DEBUG-\t RadialSolverSolution Point 4d\n")
                # These arrays are converted to complex128
                self.shear_modulus_array = np.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_complex_shape_ptr,
                    np.NPY_COMPLEX128,
                    &self.solution_storage_sptr.get().eos_solution_sptr.get().complex_shear_array_vec[0])
                
                self.bulk_modulus_array = np.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_complex_shape_ptr,
                    np.NPY_COMPLEX128,
                    &self.solution_storage_sptr.get().eos_solution_sptr.get().complex_bulk_array_vec[0])

    def __dealloc__(self):

        # Release the heap allocated storage
        del self.solution_storage_ptr

    @property
    def error_code(self):
        """ Return solution storage's error code """
        return self.solution_storage_sptr.get().error_code

    @property
    def message(self):
        """ Return solver's message """
        return str(self.solution_storage_sptr.get().message_ptr, 'UTF-8')
    
    @property
    def success(self):
        """ Return if the solver was successful message """
        return self.solution_storage_sptr.get().success

    @property
    def eos_error_code(self):
        """ Return solver's equation of state message """
        return self.solution_storage_sptr.get().eos_solution_sptr.get().error_code

    @property
    def eos_message(self):
        """ Return solver's equation of state message """
        return str(self.solution_storage_sptr.get().eos_solution_sptr.get().message_ptr, 'UTF-8')
    
    @property
    def eos_success(self):
        """ Return if the solver's equation of state sub-solver was successful """
        return self.solution_storage_sptr.get().eos_solution_sptr.get().success

    @property
    def result(self):
        """ Return result array. """

        if self.success and (self.error_code == 0):
            # TODO: Optimize solution storage so that transpose is not required?
            return self.full_solution_arr.T
        else:
            return None

    @property
    def love(self):
        """ Return all complex love numbers. """
        if self.success and (self.error_code == 0):
            return self.complex_love_arr
        else:
            return None

    @property
    def k(self):
        """ Tidal Love number k. """
        if self.success and (self.error_code == 0):
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
        if self.success and (self.error_code == 0):
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
        if self.success and (self.error_code == 0):
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
        if self.ytype_names_set and self.success and (self.error_code == 0):
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
