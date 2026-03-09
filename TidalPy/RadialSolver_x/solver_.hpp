// solver_.hpp - Top-level radial solver declarations
// Ported from TidalPy/RadialSolver/solver.pyx
//
// NOTE: The actual solver (cf_radial_solver) and Python wrapper (radial_solver)
// live in solver.pyx because they depend on:
//   - CyRK's Cython-only API (ODEMethod, PreEvalFunc)
//   - Material_x EOS solver (c_solve_eos via Cython)
//   - The shooting method (cf_shooting_solver in shooting.pyx)
//
// This header exists as a placeholder for any future C++ helper functions.
#pragma once

#include "constants_.hpp"
#include "rs_solution_.hpp"
