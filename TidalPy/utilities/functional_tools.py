from inspect import isfunction

from numba.targets.registry import CPUTarget, CPUDispatcher
from numba.array_analysis import MAP_TYPES
from numba import njit as real_njit

from ..performance import njit as tpy_njit


def is_function(potential_func, ignore_njit: bool = True) -> bool:
    """ Checks if a function is a python or numba function """

    if ignore_njit and (potential_func is real_njit or potential_func is tpy_njit):
        return False

    func_check = False
    if isfunction(potential_func):
        func_check = True
    if isinstance(potential_func, CPUDispatcher):
        func_check = True
    if isinstance(potential_func, CPUTarget):
        func_check = True
    if type(potential_func) in MAP_TYPES:
        func_check = True

    return func_check