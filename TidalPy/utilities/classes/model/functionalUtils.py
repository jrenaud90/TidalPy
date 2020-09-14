from inspect import isfunction
from typing import List, Tuple, Callable

from numba import njit as real_njit
try:
    from numba.array_analysis import MAP_TYPES
    from numba.targets.registry import CPUTarget, CPUDispatcher
except (AttributeError, ImportError, ModuleNotFoundError):
    # At some point after numba v0.45 numba experienced a big refactoring
    from numba.parfors.array_analysis import MAP_TYPES
    from numba.core.registry import CPUTarget, CPUDispatcher

from ...performance.numba import njit as tpy_njit
from ....exceptions import ModelException


def is_function(potential_func: Callable, ignore_njit: bool = True) -> bool:
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


def parse_model_docstring(model_func: Callable) -> Tuple[Tuple[str, ...], Tuple[str, ...]]:
    """ Parses a function's docstring looking for TidalPy-specific information regarding the function's input.

    Constant Argument Signature:
        "!TPY_args const: arg1, arg2, ..."

    Live Argument Signature:
        "!TPY_args live: self.<reference_1>.<reference_2>.__etc__.arg1, self.<reference_1>.<reference_2>.__etc__.arg2, ...

    Parameters
    ----------
    model_func : Callable
        Function to parse

    Returns
    -------
    const_args : List[str]
        List of constant arguments (if any) declared in the function's docstrings
    live_args : List[str]
        List of live arguments (if any) declared in the function's docstrings
    """

    live_indicator = '!tpy_args live:'
    const_indicator = '!tpy_args const:'
    func_name = model_func.__name__

    const_args = None
    live_args = None

    if not (model_func.__doc__ is None or model_func.__doc__ == ''):

        for line in model_func.__doc__.split('\n'):
            lowered_line = line.lower()

            # Look for live arguments
            if live_indicator in lowered_line:
                if live_args is not None:
                    raise ModelException(f'Multiple live arg lines found for {func_name}.')

                starting_index = lowered_line.index(live_indicator)
                rest_of_line = line[starting_index + len(live_indicator):].strip()
                live_args = rest_of_line.split(',')
                live_args = [arg.strip() for arg in live_args]
                for arg_i, arg in enumerate(live_args):
                    if 'self.' not in arg:
                        raise ModelException("Live args must be prepended with 'self.'.")

            # Look for constant arguments
            if const_indicator in lowered_line:
                if const_args is not None:
                    raise ModelException(f'Multiple constant arg lines found for {func_name}.')

                starting_index = lowered_line.index(const_indicator)
                rest_of_line = line[starting_index + len(const_indicator):].strip()
                const_args = rest_of_line.split(',')
                const_args = [arg.strip() for arg in const_args]

    # Convert output to tuples
    if const_args is None:
        const_args = tuple()
    else:
        const_args = tuple(const_args)
    if live_args is None:
        live_args = tuple()
    else:
        live_args = tuple(live_args)

    return const_args, live_args
