# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

# Minimal Cython wrapper for c_matrix_propagate.
# The primary entry point for the matrix method is through the top-level solver (solver.pyx),
# which calls c_matrix_propagate directly in C++. This file exists for potential direct testing.

from libcpp cimport bool as cpp_bool

from TidalPy.RadialSolver_x.matrix cimport c_matrix_propagate
