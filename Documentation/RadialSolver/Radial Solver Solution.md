# Radial Solver Solution
`TidalPy.RadialSolver.radial_solver` functions stores the viscoelastic-gravitational solution, results of solving a 
planet's equation of state, and other parameters and meta data in a cythonized python class `TidalPy.RadialSolver.rs_solution.RadialSolverSolution`. This document details how to access this API from both Python and Cython.

## Python API
The example below outlines which parameters are availble via Python from the `RadialSolverSolution` class instance.

```python

radial_solver_solution = radial_solver(...)  # See inputs from the previous section

# Check is the solution was successful.
# See the following sections for information on why a solution may not be successful.
radial_solver_solution.success  # == (True / False)

# Check any integration messages (e.g., hints about why the solver failed)
radial_solver_solution.message

# There is a helpful way to directly get the result you are looking for if you solved for multiple Love types.
radial_solver_solution['tidal']    # Tidal gravitational-viscoelastic results
radial_solver_solution['loading']  # Loading gravitational-viscoelastic results

# Or you can get all requested complex Love numbers
# This is a np.ndarray (complex) with the shape of [num_solve_for, 3]
# These are returned in a 2D array structured like [num_solution_types : love_number] where the first column is the
# solution type. If you are only solving for tidal Love numbers then there will only be one solution in the "0" location.
# If you are solving for multiple solutions then they will be stored in order (e.g., if you solve for tidal and load then
# tidal will be at 0 and load at 1). The second column is used to access the three Love numbers:
# - index 0: k Love number.
# - index 1: h Love number.
# - index 2: l Shida number.
# For example, if you only set `solve_for=('tidal',)` then `radial_solver_solution.love[0, 1]` == tidal k.
radial_solver_solution.love

# There are also helpful shortcuts for each individual Love number.
# These are each np.ndarrays but of only size 1 if `num_solve_for==1`
radial_solver_solution.k
radial_solver_solution.h
radial_solver_solution.l

# Additional physical properties of the planet can be accessed. Many of these are solved for during the radial_solver
# process via the equation of state solver. 
radial_solver_solution.radius            # Planet's radius (user provided)
radial_solver_solution.mass              # Planet's mass (calculated)
radial_solver_solution.moi               # Planet's moment of inertia (calculated)
radial_solver_solution.moi_factor        # Planet's moment of inertia factor (calculated as moi / ideal sphere moi)
radial_solver_solution.central_pressure  # Pressure at the center of the planet (calculated)
radial_solver_solution.surface_pressure  # Pressure at the surface of the planet (calculated)
radial_solver_solution.surface_gravity   # Acceleration due to gravity at the surface (calculated)

# Numpy arrays of these properties defined throughout the planet are also available. These are interpolated at the user
# provided radial steps (defined by the provided radius array).
radial_solver_solution.gravity_array
radial_solver_solution.pressure_array
radial_solver_solution.mass_array
radial_solver_solution.moi_array
radial_solver_solution.density_array
radial_solver_solution.shear_modulus_array
radial_solver_solution.bulk_modulus_array

```
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

    # Radial solver storage, physical properties
    @property
    def planet_radius(self):
        return self.solution_storage_sptr.get().radius_array_ptr[self.radius_array_size - 1]

    # EOS class properties
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
    def mass(self):
        """ Return's the total planet mass, found by the EOS solver """
        return self.solution_storage_sptr.get().eos_solution_sptr.get().mass
    
    @property
    def moi(self):
        """ Return's the planet's real moment of inertia, found by the EOS solver """
        return self.solution_storage_sptr.get().eos_solution_sptr.get().moi
    
    @property
    def moi_factor(self):
        """ Return's the planet's moment of inertia factor, found by the EOS solver """
        cdef double ideal_moi = (2. / 5.) * self.mass * self.planet_radius**2
        return self.moi / ideal_moi
    
    @property
    def central_pressure(self):
        """ Return's the pressure at the planet's center, found by the EOS solver """
        return self.solution_storage_sptr.get().eos_solution_sptr.get().central_pressure
    
    @property
    def surface_pressure(self):
        """ Return's the pressure at the planet's surface, found by the EOS solver """
        return self.solution_storage_sptr.get().eos_solution_sptr.get().surface_pressure
    
    @property
    def surface_gravity(self):
        """ Return's the acceleration due to gravity at the planet's surface, found by the EOS solver """
        return self.solution_storage_sptr.get().eos_solution_sptr.get().surface_gravity

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