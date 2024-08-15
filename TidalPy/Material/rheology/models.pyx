# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
from libc.stdlib cimport exit, EXIT_FAILURE
from libc.stdio cimport printf
from libc.string cimport strcmp

from libcpp cimport bool as cpp_bool

cdef dict known_rhologies = {
    'elastic': (ELASTIC_RHEOLOGY_INT, ELASTIC_NUM_PARAMETERS),
    'newton': (NEWTON_RHEOLOGY_INT, NEWTON_NUM_PARAMETERS),
    'maxwell': (MAXWELL_RHEOLOGY_INT, MAXWELL_NUM_PARAMETERS),
    'voigt': (VOIGT_RHEOLOGY_INT, VOIGT_NUM_PARAMETERS),
    'burgers': (BURGERS_RHEOLOGY_INT, BURGERS_NUM_PARAMETERS),
    'andrade': (ANDRADE_RHEOLOGY_INT, ANDRADE_NUM_PARAMETERS),
    'sundberg': (SUNDBERG_RHEOLOGY_INT, SUNDBERG_NUM_PARAMETERS),

    # Alias other names
    'sundberg-cooper': (SUNDBERG_RHEOLOGY_INT, SUNDBERG_NUM_PARAMETERS),
    'sundbergcooper': (SUNDBERG_RHEOLOGY_INT, SUNDBERG_NUM_PARAMETERS),
}

cdef class RheologyModel:

    def __init__(self,
            str rheology_name,
            double[::1] rheology_parameters = None):
        
        # Clean up rheology name and find the correct C++ rheology class
        self.rheology_name = rheology_name.lower()
        if self.rheology_name not in known_rhologies:
            raise AttributeError(f"Incorrect or unsupported rheological model requested: {rheology_name}.")
        
        cdef (char, char) rheo_data  = known_rhologies[self.rheology_name]
        self.rheology_int        = rheo_data[0]
        self.required_num_params = rheo_data[1]

        # Make a heap allocated instance of the rheology class
        if self.rheology_int == ELASTIC_RHEOLOGY_INT:
            self.rheology_class_ptr = new Elastic()
        elif self.rheology_int == NEWTON_RHEOLOGY_INT:
            self.rheology_class_ptr = new Newton()
        elif self.rheology_int == MAXWELL_RHEOLOGY_INT:
            self.rheology_class_ptr = new Maxwell()
        elif self.rheology_int == VOIGT_RHEOLOGY_INT:
            self.rheology_class_ptr = new Voigt()
        elif self.rheology_int == BURGERS_RHEOLOGY_INT:
            self.rheology_class_ptr = new Burgers()
        elif self.rheology_int == ANDRADE_RHEOLOGY_INT:
            self.rheology_class_ptr = new Andrade()
        elif self.rheology_int == SUNDBERG_RHEOLOGY_INT:
            self.rheology_class_ptr = new SundbergCooper()
        
        # Update parameters if user provided them
        self.update_parameters(rheology_parameters)
    
    def __dealloc__(self):
        del self.rheology_class_ptr
    
    def update_parameters(self,
            double[::1] new_rheology_parameters):

         # Parse additional parameters
        cdef double* parameter_storage_ptr = NULL
        cdef cpp_bool user_provided_params = False
        cdef int num_provided_params = 0
        if new_rheology_parameters is not None:
            user_provided_params = True
            num_provided_params = len(new_rheology_parameters)

        if user_provided_params:
             # Check if user provided parameters are correct
            if num_provided_params != self.required_num_params:
                raise AttributeError(f"Incorrect number of provided parameters for rheology. Provided {num_provided_params} when {self.rheology_name} requires {self.required_num_params}.")

            # Pass a pointer to the new parameters to the rheology C++ class. This will create a copy of the data stored in `new_rheology_parameters`.
            parameter_storage_ptr = &new_rheology_parameters[0]
            self.rheology_class_ptr.update_parameters(parameter_storage_ptr)
    
    def call(self,
            double frequency,
            double modulus,
            double viscosity):

        cdef double complex result

        # The C++ class only works with doubles so we need to create a pointer to a double[2] array that stores the double complex result.
        cdef double* result_dbl_ptr = <double*>&result[0]

        self.rheology_class_ptr.call(result_dbl_ptr, frequency, modulus, viscosity)
    
    def vectorize_frequency(self,
            double complex[::1] output_array,
            double[::1] frequency_array,
            double modulus,
            double viscosity):
        
        cdef size_t output_size = len(output_array)
        if len(frequency_array) != output_size:
            raise AttributeError("Frequency array and output array must be the same size in rheology `vectorize_frequency` method.")

        cdef double* output_ptr    = <double*>&output_array[0]
        cdef double* frequency_ptr = &frequency_array[0]

        self.rheology_class_ptr.vectorize_frequency(output_ptr, output_size, frequency_ptr, modulus, viscosity)

    def vectorize_modulus_viscosity(self,
            double complex[::1] output_array,
            double frequency,
            double[::1] modulus_array,
            double[::1] viscosity_array):
        
        cdef size_t output_size = len(output_array)
        if len(modulus_array) != output_size:
            raise AttributeError("Modulus array and output array must be the same size in rheology `vectorize_modulus_viscosity` method.")
        if len(viscosity_array) != output_size:
            raise AttributeError("Viscosity array and output array must be the same size in rheology `vectorize_modulus_viscosity` method.")

        cdef double* output_ptr    = <double*>&output_array[0]
        cdef double* modulus_ptr   = &modulus_array[0]
        cdef double* viscosity_ptr = &viscosity_array[0]

        self.rheology_class_ptr.vectorize_frequency(output_ptr, output_size, frequency, modulus_ptr, viscosity_ptr)
    
    def vectorize_all(self,
            double complex[::1] output_array,
            double[::1] frequency_array,
            double[::1] modulus_array,
            double[::1] viscosity_array):
        
        cdef size_t output_size = len(output_array)
        if len(frequency_array) != output_size:
            raise AttributeError("Frequency array and output array must be the same size in rheology `vectorize_all` method.")
        if len(modulus_array) != output_size:
            raise AttributeError("Modulus array and output array must be the same size in rheology `vectorize_all` method.")
        if len(viscosity_array) != output_size:
            raise AttributeError("Viscosity array and output array must be the same size in rheology `vectorize_all` method.")

        cdef double* output_ptr    = <double*>&output_array[0]
        cdef double* frequency_ptr = &frequency_array[0]
        cdef double* modulus_ptr   = &modulus_array[0]
        cdef double* viscosity_ptr = &viscosity_array[0]

        self.rheology_class_ptr.vectorize_all(output_ptr, output_size, frequency_ptr, modulus_ptr, viscosity_ptr)

    def __call__(self, frequency, modulus, viscosity):
        
        cdef cpp_bool freq_is_array = False
        cdef cpp_bool mod_is_array = False
        cdef cpp_bool visc_is_array = False
        if type(frequency) == np.ndarray:
            freq_is_array = True
        if type(modulus) == np.ndarray:
            mod_is_array = True
        if type(viscosity) == np.ndarray:
            visc_is_array = True
        
        # Deal with all float case
        if (not freq_is_array) and (not mod_is_array) and (not visc_is_array):
            return self.call(frequency, modulus, viscosity)
        
        # Deal with only frequency as an array
        cdef np.ndarray output_arr
        cdef double complex[::1] output_view
        if freq_is_array and (not mod_is_array) and (not visc_is_array):
            assert frequency.ndim == 1
            output_arr = np.empty(frequency.size, dtype=np.complex128, order="C")
            output_view = output_arr
            self.vectorize_frequency(output_view, frequency, modulus, viscosity)
            return output_arr
        
        # Deal with only viscosity and modulus as an array
        cdef double[::1] mod_array
        cdef double[::1] visc_array
        if (not freq_is_array) and (mod_is_array or visc_is_array):

            # First ensure that modulus and viscosity are both arrays
            if not mod_is_array:
                assert viscosity.ndim == 1
                visc_array = viscosity
                mod_array = modulus * np.ones(viscosity.size, dtype=np.float64, order='C')
            elif not visc_is_array:
                assert modulus.ndim == 1
                mod_array = modulus
                visc_array = viscosity * np.ones(modulus.size, dtype=np.float64, order='C')
            
            output_arr = np.empty(len(modulus), dtype=np.complex128, order='C')
            output_view = output_arr
            self.vectorize_modulus_viscosity(output_view, frequency, mod_array, visc_array)
            return output_arr
        
        # Handle case where all are passed in as arrays
        cdef size_t freq_size, mod_size, visc_size
        cdef size_t i
        cdef int freq_ndim, mod_ndim, visc_ndim
        if freq_is_array and mod_is_array and visc_is_array:
            
            freq_size = frequency.size
            freq_ndim = frequency.ndim
            mod_size  = modulus.size
            mod_ndim  = modulus.ndim
            visc_size = viscosity.size
            visc_ndim = viscosity.ndim

            # Case that they are all 1D and the same size
            if freq_ndim == 1 and mod_ndim == 1 and visc_ndim == 1:
                if (freq_size != mod_size) or (freq_size != visc_size):
                    raise AttributeError("Unsupported array sizes encountered in rheology call.")
                output_arr = np.empty(freq_size, dtype=np.complex128, order='C')
                output_view = output_arr
                self.vectorize_all(output_view, frequency, modulus, viscosity)
                return output_arr
            
            # Case where inputs are higher dimensions. Note in this case they must all be the same dimension.
            elif (freq_ndim == 2) and (freq_ndim == mod_ndim) and (freq_ndim == visc_ndim):
                if (freq_size != mod_size) or (freq_size != visc_size):
                    raise AttributeError("Unsupported array sizes encountered in rheology call.")
                
                output_arr = np.empty(frequency.shape, dtype=np.complex128, order='C')
                # Step through each dimesnion and call vectorized function
                for i in range(frequency.shape[0]):
                    output_view = output_arr[i]
                    self.vectorize_all(output_view, frequency[i], modulus[i], viscosity[i])

            else:
                raise AttributeError("Unsupported number of dimensions encoutnered in input arrays to rheology call.")


cdef shared_ptr[RheologyModelBaseCC] cf_find_and_build_rheology(
        const char rheology_int,
        double* rheology_parameters) noexcept nogil:
    
    cdef shared_ptr[RheologyModelBaseCC] rheology_class_ptr

    if rheology_int == ELASTIC_RHEOLOGY_INT:
        rheology_class_ptr = make_shared[Elastic](rheology_parameters)
    elif rheology_int == NEWTON_RHEOLOGY_INT:
        rheology_class_ptr = make_shared[Newton](rheology_parameters)
    elif rheology_int == MAXWELL_RHEOLOGY_INT:
        rheology_class_ptr = make_shared[Maxwell](rheology_parameters)
    elif rheology_int == VOIGT_RHEOLOGY_INT:
        rheology_class_ptr = make_shared[Voigt](rheology_parameters)
    elif rheology_int == BURGERS_RHEOLOGY_INT:
        rheology_class_ptr = make_shared[Burgers](rheology_parameters)
    elif rheology_int == ANDRADE_RHEOLOGY_INT:
        rheology_class_ptr = make_shared[Andrade](rheology_parameters)
    elif rheology_int == SUNDBERG_RHEOLOGY_INT:
        rheology_class_ptr = make_shared[SundbergCooper](rheology_parameters)
    else:
        printf("Unknown Rheology model integer: %d.", rheology_int)
        exit(EXIT_FAILURE)
    
    return rheology_class_ptr


cdef char cf_get_rheology_int(
        const char* rheology_name) noexcept nogil:
    # Note that the input rheology name is assumed to be lower case.
    cdef char ouput = -1

    if strcmp(rheology_name, "elastic") == 0:
        ouput = ELASTIC_RHEOLOGY_INT
    elif strcmp(rheology_name, "newton") == 0:
        ouput = NEWTON_RHEOLOGY_INT
    elif strcmp(rheology_name, "maxwell") == 0:
        ouput = MAXWELL_RHEOLOGY_INT
    elif strcmp(rheology_name, "voigt") == 0:
        ouput = VOIGT_RHEOLOGY_INT
    elif strcmp(rheology_name, "burgers") == 0:
        ouput = BURGERS_RHEOLOGY_INT
    elif strcmp(rheology_name, "andrade") == 0:
        ouput = ANDRADE_RHEOLOGY_INT
    elif strcmp(rheology_name, "sundberg") == 0:
        ouput = SUNDBERG_RHEOLOGY_INT
    
     # Other alias
    elif strcmp(rheology_name, "sundberg-cooper") == 0:
        ouput = SUNDBERG_RHEOLOGY_INT
    elif strcmp(rheology_name, "sundbergcooper") == 0:
        ouput = SUNDBERG_RHEOLOGY_INT
    
    # Unknown / Not Implemented
    else:
        printf("Unknown or not implemented rheology %s", rheology_name)
        exit(EXIT_FAILURE)

    return ouput


def find_rheology(str rheology_name, double[::1] rheology_parameters = None):

    return RheologyModel(rheology_name, rheology_parameters)
