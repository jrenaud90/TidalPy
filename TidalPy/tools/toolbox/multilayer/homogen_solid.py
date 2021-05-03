from typing import Union

import numpy as np

from ....utilities.types import ComplexArray, FloatArray
from .helper import build_static_solid_solver, build_dynamic_solid_solver
from ....tides.multilayer.numerical_int import solid_dynamic_guess, solid_static_guess, \
    liquid_dynamic_guess, liquid_static_guess

CmplxFltArray = Union[ComplexArray, FloatArray]


def calculate_homogen_solid(radius: FloatArray, shear_modulus: CmplxFltArray, bulk_modulus: CmplxFltArray,
                            density: FloatArray, gravity: FloatArray, frequency: float,
                            order_l: int = 2, use_static: bool = False, use_ts74: bool = False):

    if use_static:
        radial_derivatives_func = \
            build_static_solid_solver(radius, shear_modulus, bulk_modulus, density, gravity, order_l=order_l)
        if use_ts74:
            initial_guess = solid_guess_takeuchi_static(radius, shear_modulus, bulk_modulus, density, order_l=order_l)
        else:
            initial_guess = solid_guess_kamata_static(radius, shear_modulus, bulk_modulus, density, order_l=order_l)
    else:
        radial_derivatives_func = \
            build_dynamic_solid_solver(radius, shear_modulus, bulk_modulus, density, gravity, frequency, order_l=order_l)
        if use_ts74:
            initial_guess = \
                solid_guess_takeuchi_dynamic(radius, shear_modulus, bulk_modulus, density, frequency, order_l=order_l)
        else:
            initial_guess = \
                solid_guess_kamata_dynamic(radius, shear_modulus, bulk_modulus, density, frequency, order_l=order_l)




