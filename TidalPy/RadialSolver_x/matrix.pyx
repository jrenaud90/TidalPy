# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

# The matrix method's entry point is the top-level solver (solver.pyx), which calls
# c_matrix_propagate directly in C++; this module only re-exports the declaration for testing.

from libcpp cimport bool as cpp_bool

from TidalPy.RadialSolver_x.matrix cimport c_matrix_propagate
