// shooting_.hpp - Shooting method declarations
// Ported from TidalPy/RadialSolver/shooting.pyx
//
// NOTE: The actual shooting solver (cf_shooting_solver) lives in shooting.pyx
// because it depends on CyRK's Cython-only API (baseline_cysolve_ivp_noreturn).
// c_find_num_shooting_solutions is already defined in derivatives/odes_.hpp.
//
// This header exists as a placeholder for any future C++ helper functions
// that may be extracted from the shooting solver.
#pragma once

#include "constants_.hpp"
