# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

# The headers this module compiles read the shared TidalPy configuration, so its pointer is wired here as well.
from TidalPy.constants cimport get_shared_config_address, set_tidalpy_config_ptr
set_tidalpy_config_ptr(get_shared_config_address())

# The matrix method's entry point is the top-level solver (solver.pyx), which calls
# c_matrix_propagate directly in C++; this module only re-exports the declaration for testing.

from TidalPy.RadialSolver_x.matrix cimport c_matrix_propagate
