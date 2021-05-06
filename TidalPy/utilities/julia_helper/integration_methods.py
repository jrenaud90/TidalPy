from diffeqpy import ode, de

# Read more about Julia's ode solvers here: https://diffeq.sciml.ai/dev/solvers/ode_solve/


non_stiff_solvers = {
    # The canonical Runge-Kutta Order 4 method. Uses a defect control for adaptive stepping using maximum error over the whole interval.
    'rk4'        : ode.RK4,
    # Tsitouras 5/4 Runge-Kutta method. (free 4th order interpolant).
    'tsit5'      : ode.Tsit5,
    # Tsitouras' Runge-Kutta-Oliver 6 stage 5th order method. This method is robust on problems which have a singularity at t=0.
    'rko65'      : ode.RKO65,
    # Tsitouras-Papakostas 8/7 Runge-Kutta method.
    'tsitpap8'   : ode.TsitPap8,
    # Feagin's 10th-order Runge-Kutta method.
    'feagin10'   : ode.Feagin10,
    # Feagin's 12th-order Runge-Kutta method.
    'feagin12'   : ode.Feagin12,
    # Feagin's 14th-order Runge-Kutta method.
    'feagin14'   : ode.Feagin14,
    # Additionally, the following algorithms have a lazy interpolant:
    # BS5 - Bogacki-Shampine 3/2 Runge-Kutta method. (lazy 5th order interpolant).
    'bs3'        : ode.BS3,
    # BS5 - Bogacki-Shampine 5/4 Runge-Kutta method. (lazy 5th order interpolant).
    'bs5'        : ode.BS5,
    # Verner's "Most Efficient" 6/5 Runge-Kutta method. (lazy 6th order interpolant).
    'vern6'      : ode.Vern6,
    # Verner's "Most Efficient" 7/6 Runge-Kutta method. (lazy 7th order interpolant).
    'vern7'      : ode.Vern7,
    # Verner's "Most Efficient" 8/7 Runge-Kutta method. (lazy 8th order interpolant).
    'vern8'      : ode.Vern8,
    # Verner's "Most Efficient" 9/8 Runge-Kutta method. (lazy 9th order interpolant).
    'vern9'      : ode.Vern9,
    # A 5 parallel, 2 processor explicit Runge-Kutta method of 5th order.
    # These methods utilize multithreading on the f calls to parallelize the problem. This requires that simultaneous calls to f are thread-safe.
    'kuttaprk2p5': ode.KuttaPRK2p5
}

stiff_solvers = {
    'rosenbrock23': de.Rosenbrock23,
    # The ODEInterface algorithms are the classic Fortran algorithms. While the non-stiff algorithms are superseded by the more featured and higher performance Julia implementations from OrdinaryDiffEq.jl, the stiff solvers such as radau are some of the most efficient methods available (but are restricted for use on arrays of Float64).
    # Rosenbrock 4(3) method.
    'rodas4'      : de.Rodas4,
    'rodas5'      : de.Rodas5
}

def get_julia_solver(solver_name: str):
    try:
        solver = non_stiff_solvers[solver_name.lower()]
        ode_system = ode
    except KeyError:
        try:
            solver = stiff_solvers[solver_name.lower()]
            ode_system = de
        except KeyError:
            raise KeyError(f'Unknown Julia Integration Model: {solver_name}.')

    return ode_system, solver
