# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libcpp.memory cimport make_shared

from TidalPy.RadialSolver.constants cimport MAX_NUM_Y
from TidalPy.constants cimport d_PI_DBL

cimport numpy as cnp
import numpy as np
cnp.import_array()

from TidalPy.exceptions import ArgumentException
from TidalPy.logger import get_logger

log = get_logger("TidalPy")


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
        self.solution_storage_sptr = make_shared[RadialSolutionStorageCC](
            self.num_ytypes,
            upper_radius_bylayer_ptr,
            self.num_layers,
            radius_array_ptr,
            radius_array_size,
            degree_l)
        self.solution_storage_ptr = self.solution_storage_sptr.get()

        if not self.solution_storage_ptr:
            # C++ solution storage could not be initialized.
            raise RuntimeError("RadialSolutionStorageCC extension class could not be initialized.")

        # Finish initialization with provided array
        self.change_radius_array(radius_array_ptr, radius_array_size, array_changed=False)
    
    def __dealloc__(self):
        self.solution_storage_sptr.reset()
        self.solution_storage_ptr = NULL

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
                self.solution_storage_ptr.error_code = -2
                self.solution_storage_ptr.set_message("ArgumentException:: Unknown boundary condition provided")
        self.ytype_names_set = True

    cdef void change_radius_array(
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
            self.solution_storage_ptr.change_radius_array(new_radius_array_ptr, new_size_radius_array, array_changed)
        
        # Now that the underlying storage has been updated we can build our numpy array wrappers.
    
        # The RadialSolutionStorageCC class has full control of memory. This class simply wraps it.
        # The shape needs to be twice what we expect because the underlying C++ class only works with double arrays 
        # but at this level we want arrays to be double complex (the 2x factor is built into MAX_NUM_Y_REAL)
        cdef cnp.npy_intp[2] full_solution_shape   = [self.radius_array_size, self.num_ytypes * MAX_NUM_Y]
        cdef cnp.npy_intp* full_solution_shape_ptr = &full_solution_shape[0]
        cdef cnp.npy_intp full_solution_shape_ndim = 2

        # Initialize Love numbers (3 love numbers for each requested y-type)
        # Love numbers are stored (k, h, l)_ytype0, (k, h, l)_ytype1, (k, h, l)_ytype2, ...
        cdef cnp.npy_intp[1] love_shape   = [self.num_ytypes * 3]
        cdef cnp.npy_intp* love_shape_ptr = &love_shape[0]
        cdef cnp.npy_intp love_shape_ndim = 1
        
        # Make numpy arrays that wrap all of the equation of state class vectors in a similar manner to the above.
        cdef cnp.npy_intp[1] eos_float_shape     = [self.radius_array_size]
        cdef cnp.npy_intp* eos_float_shape_ptr   = &eos_float_shape[0]
        cdef cnp.npy_intp[1] eos_complex_shape   = [self.radius_array_size]
        cdef cnp.npy_intp* eos_complex_shape_ptr = &eos_complex_shape[0]
        cdef cnp.npy_intp eos_ndim               = 1

        # `solution_storage_sptr.full_solution_vec` is a double vector but it is really storing double complex data
        # ordered by y0_real, y0_imag, y1_real, y1_image, ... so we can safely convert it to a complex128 np.ndarray
        cdef EOSSolutionCC* eos_solution_ptr = self.solution_storage_ptr.get_eos_solution_ptr()

        if not self.solution_storage_ptr:
            raise RuntimeError("RadialSolverSolution (PyClass):: RadialSolutionStorageCC extension class is not initialized.")
        else:
            self.full_solution_arr = cnp.PyArray_SimpleNewFromData(
                full_solution_shape_ndim,
                full_solution_shape_ptr,
                cnp.NPY_COMPLEX128,
                <double complex*>self.solution_storage_ptr.full_solution_vec.data())

            # Same note as above, `solution_storage_sptr.complex_love_vec` is a double vector that we are converting to a
            # complex128 np.ndarray.
            self.complex_love_arr = cnp.PyArray_SimpleNewFromData(
                love_shape_ndim,
                love_shape_ptr,
                cnp.NPY_COMPLEX128,
                <double complex*>self.solution_storage_ptr.complex_love_vec.data())

            if not eos_solution_ptr:
                raise RuntimeError("RadialSolverSolution (PyClass):: EOSSolutionCC extension class is not initialized.")
            else:
                self.radius_array_cnp = cnp.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_float_shape_ptr,
                    cnp.NPY_FLOAT64,
                    eos_solution_ptr.radius_array_vec.data())

                self.gravity_array_cnp = cnp.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_float_shape_ptr,
                    cnp.NPY_FLOAT64,
                    eos_solution_ptr.gravity_array_vec.data())

                self.pressure_array_cnp = cnp.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_float_shape_ptr,
                    cnp.NPY_FLOAT64,
                    eos_solution_ptr.pressure_array_vec.data())
                
                self.mass_array_cnp = cnp.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_float_shape_ptr,
                    cnp.NPY_FLOAT64,
                    eos_solution_ptr.mass_array_vec.data())
                
                self.moi_array_cnp = cnp.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_float_shape_ptr,
                    cnp.NPY_FLOAT64,
                    eos_solution_ptr.moi_array_vec.data())

                self.density_array_cnp = cnp.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_float_shape_ptr,
                    cnp.NPY_FLOAT64,
                    eos_solution_ptr.density_array_vec.data())

                # These arrays are converted to complex128
                self.shear_modulus_array_cnp = cnp.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_complex_shape_ptr,
                    cnp.NPY_COMPLEX128,
                    <double complex*>eos_solution_ptr.complex_shear_array_vec.data())
                
                self.bulk_modulus_array_cnp = cnp.PyArray_SimpleNewFromData(
                    eos_ndim,
                    eos_complex_shape_ptr,
                    cnp.NPY_COMPLEX128,
                    <double complex*>eos_solution_ptr.complex_bulk_array_vec.data())
    
    cdef void finalize_python_storage(self) noexcept:

        # Setup EOS diagnostic arrays
        cdef cnp.npy_intp[2] eos_steps_taken_shape   = [self.solution_storage_ptr.get_eos_solution_ptr().num_cyolver_calls / self.num_layers, self.num_layers]
        cdef cnp.npy_intp* eos_steps_taken_shape_ptr = &eos_steps_taken_shape[0]
        cdef cnp.npy_intp eos_steps_taken_ndim       = 2
        self.eos_steps_taken_array = cnp.PyArray_SimpleNewFromData(
            eos_steps_taken_ndim,
            eos_steps_taken_shape_ptr,
            # The below should be "cnp.NPY_UINTP" to ensure that the pointer is for a size_t*. 
            # But it looks like numpy does not have an enum for this, even though their documentation says it does:
            # https://numpy.org/devdocs/reference/c-api/dtype.html#c.NPY_TYPES.NPY_UINTP
            # Track question related to this here: https://stackoverflow.com/questions/79231405/no-enum-for-numpy-uintp
            cnp.NPY_UINT64,
            self.solution_storage_ptr.get_eos_solution_ptr().steps_taken_vec.data())

        # Setup shooting method diagnostic arrays
        cdef cnp.npy_intp[2] steps_taken_shape   = [self.num_layers, 3]
        cdef cnp.npy_intp* steps_taken_shape_ptr = &steps_taken_shape[0]
        cdef cnp.npy_intp steps_taken_ndim       = 2
        self.shooting_method_steps_taken_array = cnp.PyArray_SimpleNewFromData(
            steps_taken_ndim,
            steps_taken_shape_ptr,
            # The below should be "cnp.NPY_UINTP" to ensure that the pointer is for a size_t*. 
            # But it looks like numpy does not have an enum for this, even though their documentation says it does:
            # https://numpy.org/devdocs/reference/c-api/dtype.html#c.NPY_TYPES.NPY_UINTP
            # Track question related to this here: https://stackoverflow.com/questions/79231405/no-enum-for-numpy-uintp
            cnp.NPY_UINT64,
            self.solution_storage_ptr.shooting_method_steps_taken_vec.data())
    
    def eos_call(self, double radius):
        
        cdef EOSSolutionCC* eos_solution_ptr = self.solution_storage_ptr.get_eos_solution_ptr()

        cdef int layer_index = -1
        cdef size_t layer_i
        cdef double layer_r = 0.0
        cdef double last_layer_r = 0.0

        for layer_i in range(eos_solution_ptr.upper_radius_bylayer_vec.size()):
            layer_r = eos_solution_ptr.upper_radius_bylayer_vec[layer_i]
            if last_layer_r <= radius < layer_r:
                layer_index = <int>layer_i
                break
            last_layer_r = layer_r
            
        if layer_index < 0:
            raise ValueError("Could not find correct layer for provided radius.")

        cdef cnp.ndarray[cnp.float64_t, ndim=1] eos_interp = np.empty(9, dtype=np.float64, order='C')
        cdef double[::1] eos_interp_view = eos_interp
        cdef double* eos_interp_ptr      = &eos_interp_view[0]

        eos_solution_ptr.call(<size_t>layer_index, radius, eos_interp_ptr)
        return eos_interp

    def plot_ys(self):
        cdef list result_list
        cdef list radius_list
        cdef list labels
        cdef size_t ytype_i
        cdef size_t stored_ytypes = 0
        cdef str ytype_name

        if self.success:
            from TidalPy.utilities.graphics.multilayer import yplot
            if self.num_ytypes <= 0:
                raise AttributeError("`RadialSolverSolution` can not plot ys because number of ytypes is less than 1.")
            elif self.num_ytypes == 1:
                if self.result is not None:
                    return yplot(self.result, self.radius_array)
                else:
                    raise AttributeError("`RadialSolverSolution` can not plot ys because result is None (perhaps failed solution?).")
            else:
                result_list = list()
                radius_list = list()
                labels      = list()
                for ytype_i in range(self.num_ytypes):
                    ytype_name = str(self.ytypes[ytype_i], 'UTF-8')
                    if self.get_result_by_ytype_name(ytype_name) is not None:
                        result_list.append(self.get_result_by_ytype_name(ytype_name))
                        radius_list.append(self.radius_array)
                        labels.append(ytype_name.title())
                        stored_ytypes += 1
                if stored_ytypes > 1:
                    return yplot(result_list, radius_list, labels=labels)
                else:
                    raise AttributeError("`RadialSolverSolution` can not plot ys because result is None (perhaps failed solution?).")
        else:
            raise AttributeError("`RadialSolverSolution` can not plot ys because result was not successful.")
        
    def plot_interior(self):
        if self.eos_success:
            from TidalPy.utilities.graphics.planet_plot import planet_plot

            return planet_plot(
                self.radius_array,
                self.gravity_array,
                self.pressure_array,
                self.density_array,
                None,
                self.shear_modulus_array,
                self.bulk_modulus_array,
                self.radius,
                self.density_bulk,
                planet_name=None,
                use_scatter=False,
                depth_plot=False,
                auto_show=True,
                annotate=True)
    
    def print_diagnostics(self, cpp_bool print_diagnostics = True, cpp_bool log_diagnostics = False):
        cdef str log_message = ""
        log_message += "\n\tEquation of State Solver:"
        log_message += f"\n\t\tSuccess:           {self.eos_success}"
        log_message += f"\n\t\tError code:        {self.eos_error_code}"
        log_message += f"\n\t\tMessage:           {self.eos_message}"
        if self.eos_success:
            log_message += f"\n\t\tIterations:        {self.eos_iterations}"
            log_message += f"\n\t\tPressure Error:    {self.eos_pressure_error:0.3e}"
            log_message += f"\n\t\tCentral Pressure:  {self.central_pressure:0.3e}"
            log_message += f"\n\t\tMass:              {self.eos_pressure_error:0.3e}"
            log_message += f"\n\t\tMOI (factor):      {self.moi:0.3e} ({self.moi_factor:0.3f})"
            log_message += f"\n\t\tSurface gravity:   {self.surface_gravity:0.3e}\n"
        log_message += "\n\tRadial Solver Results:"
        log_message += f"\n\t\tSuccess:     {self.success}"
        log_message += f"\n\t\tError code:  {self.error_code}"
        log_message += f"\n\t\tMessage:     {self.message}"
        log_message += f"\n\t\tSteps Taken (per sub-solution):"
        cdef size_t layer_i
        for layer_i in range(self.num_layers):
            log_message += f"\n\t\t\tLayer {layer_i} = {self.steps_taken[layer_i]}"
        if self.success:
            log_message += f"\n\t\tk_{self.degree_l} = {self.k}"
            log_message += f"\n\t\th_{self.degree_l} = {self.h}"
            log_message += f"\n\t\tl_{self.degree_l} = {self.l}"
        if print_diagnostics:
            print(log_message)
            return None
        if log_diagnostics:
            log.info(log_message)
            return None
        if not print_diagnostics and not log_diagnostics:
            return log_message

    def __dealloc__(self):

        # Release the heap allocated storage
        self.solution_storage_sptr.reset()
    
    # Radial solver storage, feedback properties
    @property
    def error_code(self):
        """ Return solution storage's error code """
        return self.solution_storage_ptr.error_code

    @property
    def message(self):
        """ Return solver's message """
        return str(self.solution_storage_ptr.message_ptr, 'UTF-8')
    
    @property
    def success(self):
        """ Return if the solver was successful message """
        return self.solution_storage_ptr.success

    # EOS class properties
    @property
    def eos_error_code(self):
        """ Return solver's equation of state message """
        return self.solution_storage_ptr.get_eos_solution_ptr().error_code

    @property
    def eos_message(self):
        """ Return solver's equation of state message """
        return str(self.solution_storage_ptr.get_eos_solution_ptr().message_ptr, 'UTF-8')
    
    @property
    def eos_success(self):
        """ Return if the solver's equation of state sub-solver was successful """
        return self.solution_storage_ptr.get_eos_solution_ptr().success
    
    @property
    def eos_pressure_error(self):
        """ Return the surface pressure error found by the equation of state sub-solver """
        return self.solution_storage_ptr.get_eos_solution_ptr().pressure_error
    
    @property
    def eos_iterations(self):
        """ Return the number of iterations performed by the EOS sub-solver to find convergence on surface pressure """
        return self.solution_storage_ptr.get_eos_solution_ptr().iterations
    
    @property
    def eos_steps_taken(self):
        """ Return the number of integration steps required by the last iteration of the EOS solver """
        return np.copy(self.eos_steps_taken_array)
    
    # EOS Solution arrays
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

    # EOS Scalar parameters
    @property
    def radius(self):
        """ Return's the planet's radius, set by user """
        return self.solution_storage_ptr.get_eos_solution_ptr().radius
    
    @property
    def volume(self):
        """ Return's the planet's volume, calculated by its radius """
        return (4.0 / 3.0) * d_PI_DBL * self.radius**3

    @property
    def mass(self):
        """ Return's the total planet mass, found by the EOS solver """
        return self.solution_storage_ptr.get_eos_solution_ptr().mass
    
    @property
    def moi(self):
        """ Return's the planet's real moment of inertia, found by the EOS solver """
        return self.solution_storage_ptr.get_eos_solution_ptr().moi
    
    @property
    def moi_factor(self):
        """ Return's the planet's moment of inertia factor, found by the EOS solver """
        cdef double ideal_moi = (2. / 5.) * self.mass * self.radius**2
        return self.moi / ideal_moi
    
    @property
    def density_bulk(self):
        """ Return's the planet's bulk density """
        return self.mass / self.volume
    
    @property
    def central_pressure(self):
        """ Return's the pressure at the planet's center, found by the EOS solver """
        return self.solution_storage_ptr.get_eos_solution_ptr().central_pressure
    
    @property
    def surface_pressure(self):
        """ Return's the pressure at the planet's surface, found by the EOS solver """
        return self.solution_storage_ptr.get_eos_solution_ptr().surface_pressure
    
    @property
    def surface_gravity(self):
        """ Return's the acceleration due to gravity at the planet's surface, found by the EOS solver """
        return self.solution_storage_ptr.get_eos_solution_ptr().surface_gravity
    
    @property
    def layer_upper_radius_array(self):
        """ Return's the upper radius at each layer in the planet. """
        cdef cnp.ndarray[cnp.float64_t, ndim=1] upper_radius_array = np.empty(self.num_layers, dtype=np.float64)
        cdef EOSSolutionCC* eos_solution_ptr = self.solution_storage_ptr.get_eos_solution_ptr()
        cdef size_t layer_i
        for layer_i in range(self.num_layers):
            upper_radius_array[layer_i] = eos_solution_ptr.upper_radius_bylayer_vec[layer_i]
        return upper_radius_array

    # RadialSolver result properties
    @property
    def degree_l(self):
        """ Return the harmonic degree used to perform the calculation. """
        return self.solution_storage_ptr.degree_l

    @property
    def result(self):
        """ Return result array. """

        if self.success and (self.error_code == 0):
            # TODO: Optimize solution storage so that transpose is not required?
            return np.copy(self.full_solution_arr).T
        else:
            return None

    @property
    def love(self):
        """ Return all complex love numbers. """
        if self.success and (self.error_code == 0):
            return np.copy(self.complex_love_arr).reshape(self.num_ytypes, 3)
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
                return np.copy(self.complex_love_arr[0::3])
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
                return np.copy(self.complex_love_arr[1::3])
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
                return np.copy(self.complex_love_arr[2::3])
        else:
            return None
    
    @property
    def steps_taken(self):
        """ Number of integration steps by layer and solution. """
        return np.copy(self.shooting_method_steps_taken_array)

    def get_result_by_ytype_name(self, str ytype_name):
        """Get a specific solution type array."""
        cdef char ytype_i
        cdef char requested_sol_num = 0
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
                raise ArgumentException('Unknown solution type requested. Key must match names passed to radial_solver "solve_for" argument and be lower case.')
            
            # Slice the result and return only the requested solution type.
            return np.copy(self.result[MAX_NUM_Y * (requested_sol_num): MAX_NUM_Y * (requested_sol_num + 1)])
        else:
            return None

    def __len__(self):
        """Return number of solution types."""
        return <Py_ssize_t>self.num_ytypes
    
    def __getitem__(self, str ytype_name):
        """Wrapper for `get_result_by_ytype_name`."""
        return self.get_result_by_ytype_name(ytype_name)
