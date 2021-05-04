from typing import Union

import numpy as np

from ....utilities.types import ComplexArray, FloatArray
from .helper import build_static_solid_solver, build_dynamic_solid_solver
from ....tides.multilayer.numerical_int import solid_dynamic_guess, solid_static_guess

def calculate_homogen_solid(radius: np.ndarray, shear_modulus: np.ndarray, bulk_modulus: np.ndarray,
                            density: np.ndarray, gravity: np.ndarray, frequency: float,
                            order_l: int = 2, use_static: bool = False, use_ts74: bool = False):

    # Initial (base) guess will be for a solid layer
    if use_static:
        initial_values = solid_static_guess(radius[0], shear_modulus[0], bulk_modulus[0], density[0], order_l=order_l)
    else:
        initial_values = solid_dynamic_guess(radius[0], shear_modulus[0], bulk_modulus[0], density[0], order_l=order_l)

    # Find the differential equation
    if use_static:
        radial_derivative = \
            build_static_solid_solver(radius, shear_modulus, bulk_modulus, density, gravity, order_l=order_l)
    else:
        radial_derivative = \
            build_dynamic_solid_solver(radius, shear_modulus, bulk_modulus, density, gravity, frequency,
                                       order_l=order_l)

        # Make sure the initial value remains complex.
        initial_value = np.asarray(initial_value, dtype=np.complex128)
        solution = solve_ivp(deriv_liq, span_core, initial_value,
                             method=integration_method, rtol=integration_tol,
                             t_eval=radius_array[radius_array <= R_core])
        if solution.status == 0:
            print(f'Solution {sn} solved.')
            solution_rs_core.append(solution.t)
            solution_ys_core.append(solution.y)
        else:
            print(solution.message)


