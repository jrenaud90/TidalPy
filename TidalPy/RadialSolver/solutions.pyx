# distutils: language = c
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libc.math cimport NAN
from libc.stdlib cimport exit, EXIT_FAILURE
from libc.stdio cimport printf

import numpy as np
cimport numpy as np

from CyRK.utils.utils cimport allocate_mem, reallocate_mem, free_mem
from TidalPy.RadialSolver.constants cimport MAX_NUM_Y, MAX_NUM_Y_REAL, MAX_NUM_SOL


cdef class RadialSolverSolution():

    def __init__(
            self,
            size_t num_slices,
            size_t num_ytypes
            ):

        # loop indicies
        cdef size_t i
        cdef size_t love_array_size

        # Initialize pointers
        self.full_solution_ptr = NULL
        self.complex_love_ptr = NULL

        # Initialize status
        self._message = 'RadialSolverSolution has not had its status set.'
        self.success = False

        # Store number of solution types
        self.num_ytypes = num_ytypes
        self.ytypes_set = False
        
        # Store size information
        self.num_slices = num_slices
        self.total_size = MAX_NUM_Y * self.num_slices * self.num_ytypes

        # Have the radial solver take control of the full solution memory
        self.full_solution_ptr = <double complex*> allocate_mem(
            self.total_size * sizeof(double complex),
            'full_solution_ptr (RadialSolverSolution; init)'
            )
        if not (self.full_solution_ptr is NULL):
            self.full_solution_view = <double complex[:self.total_size]> self.full_solution_ptr
        
        # Set all values of the solution array to NANs. This ensures that if there is a problem with the solver then
        #  the solutions will be defined (but nan).
        for i in range(self.total_size):
            self.full_solution_ptr[i] = NAN
        
        # Initialize Love numbers (3 love numbers for each requested y-type)
        # Love numbers are stored (k, h, l)_ytype0, (k, h, l)_ytype1, (k, h, l)_ytype2, ...
        love_array_size = 3 * self.num_ytypes
        self.complex_love_ptr = <double complex*> allocate_mem(
            love_array_size * sizeof(double complex),
            'complex_love_ptr (RadialSolverSolution; init)'
            )
        if not (self.complex_love_ptr is NULL):
            self.complex_love_view = <double complex[:love_array_size]> self.complex_love_ptr
        for i in range(love_array_size):
            self.complex_love_ptr[i] = NAN

    cdef void set_models(self, int* bc_models_ptr) noexcept nogil:
        # Unfortunately this must be done outside of __init__ because the argument is a pointer and python interpretor
        # is used during __init__ and it does not like pointers.

        # Find solution types
        cdef int bc_model
        for i in range(self.num_ytypes):
            bc_model = bc_models_ptr[i]
            if bc_model == 0:
                self.ytypes[i] = "free"
            elif bc_model == 1:
                self.ytypes[i] = "tidal"
            elif bc_model == 2:
                self.ytypes[i] = "loading"
            else:
                printf("AttributeError:: Unknown boundary condition")
                exit(EXIT_FAILURE)
        self.ytypes_set = True

    cdef void set_message(self, str new_message) noexcept:
        self._message_bytes = new_message.encode('utf-8') + b'\x00'
        self._message = self._message_bytes

    @property
    def message(self):
        """ Return solver's message """
        return str(self._message, 'UTF-8')
    
    @property
    def result(self):
        """ Return result array. """

        if self.success and self.ytypes_set:
            return np.ascontiguousarray(
                self.full_solution_view,
                dtype=np.complex128
                ).reshape((self.num_slices, self.num_ytypes * MAX_NUM_Y)).T
        else:
            return None
        
    @property
    def love(self):
        """ Return all complex love numbers. """
        if self.success and self.ytypes_set:
            return np.ascontiguousarray(
                self.complex_love_view,
                dtype=np.complex128
            ).reshape((self.num_ytypes, 3))
        else:
            return None

    @property
    def k(self):
        """ Tidal Love number k. """
        if self.success and self.ytypes_set:
            return np.ascontiguousarray(
                self.complex_love_view[0::3],
                dtype=np.complex128
            )
        else:
            return None

    @property
    def h(self):
        """ Tidal Love number h. """
        if self.success and self.ytypes_set:
            return np.ascontiguousarray(
                self.complex_love_view[1::3],
                dtype=np.complex128
            )
        else:
            return None
    
    @property
    def l(self):
        """ Tidal Shida number l. """
        if self.success and self.ytypes_set:
            return np.ascontiguousarray(
                self.complex_love_view[2::3],
                dtype=np.complex128
            )
        else:
            return None
    
    def __len__(self):
        """Return number of solution types."""
        return <Py_ssize_t>self.num_ytypes
    
    def __getitem__(self, str ytype_name):
        """Get a specific solution type array."""
        
        cdef size_t i
        cdef size_t requested_sol_num = 0
        cdef bint found = False
        cdef str sol_test_name
        for i in range(self.num_ytypes):
            sol_test_name = str(self.ytypes[i], 'UTF-8')
            if sol_test_name == ytype_name:
                requested_sol_num = i
                found = True
                break
        if not found:
            raise AttributeError('Unknown solution type requested. Key must match names passed to radial_solver "solve_for" argument and be lower case.')
        
        # Slice the result and return only the requested solution type.
        if self.success and self.ytypes_set:
            return self.result[MAX_NUM_Y * (requested_sol_num): MAX_NUM_Y * (requested_sol_num + 1)]
        else:
            return None

    def __dealloc__(self):

        # The RadialSolverSolution class has full control of the solution so it is responsible for releasing its memory.
        if not (self.full_solution_ptr is NULL):
            free_mem(self.full_solution_ptr)
        if not (self.complex_love_ptr is NULL):
            free_mem(self.complex_love_ptr)


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
