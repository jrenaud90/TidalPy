# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libc.stdio cimport printf
from libcpp.memory cimport make_shared

from TidalPy.RadialSolver.constants cimport MAX_NUM_Y
cimport numpy as cnp
cnp.import_array()

cdef class RadialSolverSolution:

    def __init__(
            self,
            size_t num_ytypes,
            double[::1] upper_radius_bylayer_view,
            double[::1] radius_array_view
            ):

        printf("DEBUG-\t RadialSolverSolution Point 1\n")
        # Build pointers
        cdef double* upper_radius_bylayer_ptr = &upper_radius_bylayer_view[0]
        cdef size_t num_layers                = upper_radius_bylayer_view.size
        cdef double* radius_array_ptr         = &radius_array_view[0]
        cdef size_t radius_array_size         = radius_array_view.size

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
                self.solution_storage_sptr.get().error_code = -2
                self.solution_storage_sptr.get().set_message("AttributeError:: Unknown boundary condition provided")
        self.ytype_names_set = True

    cdef change_radius_array(
            self,
            double* new_radius_array_ptr,
            size_t new_size_radius_array,
            cpp_bool array_changed = True) noexcept:
        """ Change the radius array used by RadialSolver and Equation of State solvers. """

        # Update attributes in this class
        self.radius_array_size = new_size_radius_array

        if array_changed:
            # Radius array held by the underlying C++ classes is no longer correct.
            # Update it using its methods.
            self.solution_storage_sptr.get().change_radius_array(new_radius_array_ptr, new_size_radius_array, array_changed)
        
        # Now that the underlying storage has been updated we can build our numpy array wrappers.
    
        # The RadialSolutionStorageCC class has full control of memory. This class simply wraps it.
        # The shape needs to be twice what we expect because the underlying C++ class only works with double arrays 
        # but at this level we want arrays to be double complex (the 2x factor is built into MAX_NUM_Y_REAL)
        printf("DEBUG-\t RadialSolverSolution Point 3\n")
        cdef cnp.npy_intp[2] full_solution_shape   = [self.radius_array_size, self.num_ytypes * MAX_NUM_Y]
        cdef cnp.npy_intp* full_solution_shape_ptr = &full_solution_shape[0]
        cdef cnp.npy_intp full_solution_shape_ndim = 2

        # Initialize Love numbers (3 love numbers for each requested y-type)
        # Love numbers are stored (k, h, l)_ytype0, (k, h, l)_ytype1, (k, h, l)_ytype2, ...
        cdef cnp.npy_intp[2] love_shape   = [3, 0]
        cdef cnp.npy_intp* love_shape_ptr = &love_shape[0]
        cdef cnp.npy_intp love_shape_ndim = 1
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
        cdef cnp.npy_intp[1] eos_float_shape     = [self.radius_array_size]
        cdef cnp.npy_intp* eos_float_shape_ptr   = &eos_float_shape[0]
        cdef cnp.npy_intp[1] eos_complex_shape   = [self.radius_array_size]
        cdef cnp.npy_intp* eos_complex_shape_ptr = &eos_complex_shape[0]
        cdef cnp.npy_intp eos_ndim               = 1

        # `solution_storage_sptr.full_solution_vec` is a double vector but it is really storing double complex data
        # ordered by y0_real, y0_imag, y1_real, y1_image, ... so we can safely convert it to a complex128 np.ndarray
        cdef EOSSolutionCC* eos_solution_ptr = self.solution_storage_sptr.get().get_eos_solution_ptr()

        printf("DEBUG-\t RadialSolverSolution Point 4\n")
        if not self.solution_storage_sptr.get():
            raise RuntimeError("RadialSolverSolution (PyClass):: RadialSolutionStorageCC extension class is not initialized.")
        else:
            printf("DEBUG-\t RadialSolverSolution Point 4a\n")
            self.full_solution_arr = cnp.PyArray_SimpleNewFromData(
                full_solution_shape_ndim,
                full_solution_shape_ptr,
                cnp.NPY_COMPLEX128,
                <double complex*>self.solution_storage_sptr.get().full_solution_vec.data())

            # Same note as above, `solution_storage_sptr.complex_love_vec` is a double vector that we are converting to a
            # complex128 np.ndarray.
            printf("DEBUG-\t RadialSolverSolution Point 4b\n")
            self.complex_love_arr = cnp.PyArray_SimpleNewFromData(
                love_shape_ndim,
                love_shape_ptr,
                cnp.NPY_COMPLEX128,
                <double complex*>self.solution_storage_sptr.get().complex_love_vec.data())

            if not eos_solution_ptr:
                raise RuntimeError("RadialSolverSolution (PyClass):: EOSSolutionCC extension class is not initialized.")
            else:
                printf("DEBUG-\t RadialSolverSolution Point 4c\n")
                self.gravity_array = cnp.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_float_shape_ptr,
                    cnp.NPY_FLOAT64,
                    eos_solution_ptr.gravity_array_vec.data())

                self.pressure_array = cnp.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_float_shape_ptr,
                    cnp.NPY_FLOAT64,
                    eos_solution_ptr.pressure_array_vec.data())
                
                self.mass_array = cnp.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_float_shape_ptr,
                    cnp.NPY_FLOAT64,
                    eos_solution_ptr.mass_array_vec.data())
                
                self.moi_array = cnp.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_float_shape_ptr,
                    cnp.NPY_FLOAT64,
                    eos_solution_ptr.moi_array_vec.data())

                self.density_array = cnp.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_float_shape_ptr,
                    cnp.NPY_FLOAT64,
                    eos_solution_ptr.density_array_vec.data())

                printf("DEBUG-\t RadialSolverSolution Point 4d\n")
                # These arrays are converted to complex128
                self.shear_modulus_array = cnp.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_complex_shape_ptr,
                    cnp.NPY_COMPLEX128,
                    <double complex*>eos_solution_ptr.complex_shear_array_vec.data())
                
                self.bulk_modulus_array = cnp.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_complex_shape_ptr,
                    cnp.NPY_COMPLEX128,
                    <double complex*>eos_solution_ptr.complex_bulk_array_vec.data())

    def __dealloc__(self):

        # Release the heap allocated storage
        self.solution_storage_sptr.reset()
    
    # Radial solver storage, feedback properties
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

    # EOS class properties
    @property
    def eos_error_code(self):
        """ Return solver's equation of state message """
        return self.solution_storage_sptr.get().get_eos_solution_ptr().error_code

    @property
    def eos_message(self):
        """ Return solver's equation of state message """
        return str(self.solution_storage_sptr.get().get_eos_solution_ptr().message_ptr, 'UTF-8')
    
    @property
    def eos_success(self):
        """ Return if the solver's equation of state sub-solver was successful """
        return self.solution_storage_sptr.get().get_eos_solution_ptr().success
    
    @property
    def eos_pressure_error(self):
        """ Return the surface pressure error found by the equation of state sub-solver """
        return self.solution_storage_sptr.get().get_eos_solution_ptr().pressure_error
    
    @property
    def eos_iterations(self):
        """ Return the number of iterations performed by the EOS sub-solver to find convergence on surface pressure """
        return self.solution_storage_sptr.get().get_eos_solution_ptr().iterations

    @property
    def radius(self):
        """ Return's the planet's radius, set by user """
        return self.solution_storage_sptr.get().get_eos_solution_ptr().radius

    @property
    def mass(self):
        """ Return's the total planet mass, found by the EOS solver """
        return self.solution_storage_sptr.get().get_eos_solution_ptr().mass
    
    @property
    def moi(self):
        """ Return's the planet's real moment of inertia, found by the EOS solver """
        return self.solution_storage_sptr.get().get_eos_solution_ptr().moi
    
    @property
    def moi_factor(self):
        """ Return's the planet's moment of inertia factor, found by the EOS solver """
        cdef double ideal_moi = (2. / 5.) * self.mass * self.radius**2
        return self.moi / ideal_moi
    
    @property
    def central_pressure(self):
        """ Return's the pressure at the planet's center, found by the EOS solver """
        return self.solution_storage_sptr.get().get_eos_solution_ptr().central_pressure
    
    @property
    def surface_pressure(self):
        """ Return's the pressure at the planet's surface, found by the EOS solver """
        return self.solution_storage_sptr.get().get_eos_solution_ptr().surface_pressure
    
    @property
    def surface_gravity(self):
        """ Return's the acceleration due to gravity at the planet's surface, found by the EOS solver """
        return self.solution_storage_sptr.get().get_eos_solution_ptr().surface_gravity

    # RadialSolver result properties
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
