import copy
from types import ModuleType
from typing import Callable, Tuple, Dict, Union, Any
from warnings import warn

from .functionalUtils import is_function, parse_model_docstring
from ..config.dictionary_utils import nested_get


def find_all_models(module: ModuleType, ignore_functional_types: tuple = tuple()) -> \
        Tuple[Dict[str, Callable], Dict[str, Tuple[str, ...]], Dict[str, Tuple[str, ...]]]:
    """ Find all functions within the provided module

    Provides the user with information regarding the functions inputs.

    Parameters
    ----------
    module : ModuleType
        A python module full of function definitions.
    ignore_functional_types: tuple = (,)
        A tuple of additional types to ignore when looking for functions.

    Returns
    -------
    models : Dict[str, Callable]
        Dictionary of function names: function.
    model_const_args : Dict[str, Tuple[str, ...]]
        Dictionary of function names: tuple of names of constant arguments for that function.
    model_live_args : Dict[str, Tuple[str, ...]]
        Dictionary of function names: tuple of names of live arguments for that function.
    """

    models = dict()  # type: Dict['str', Callable]
    model_const_args = dict()  # type: Dict[str, Tuple[str, ...]]
    model_live_args = dict()  # type: Dict[str, Tuple[str, ...]]

    for item_name, item in module.__dict__.items():

        skip = False
        if type(item) in ignore_functional_types:
            # Function was in the additional types to be ignored.
            skip = True
        else:
            for ignore_func in ignore_functional_types:
                if item is ignore_func:
                    skip = True
                    break
        if skip:
            continue

        if is_function(item):
            const_args, live_args = parse_model_docstring(item)
            models[item_name] = item
            model_const_args[item_name] = const_args
            model_live_args[item_name] = live_args

    return models, model_const_args, model_live_args


def build_model_default_inputs(const_arg_dict: Dict[str, Tuple[str, ...]], dict_of_defaults: dict,
                               inner_keys: Union[Tuple[str, ...], str] = None) -> Dict[str, Dict[str, Tuple[Any, ...]]]:
    """ Builds a dictionary of default input parameters using a constant argument dictionary and a dictionary of
    default values.

    Parameters
    ----------
    const_arg_dict : Dict[str, Tuple[str, ...]]
        Constant argument dictionary is assumed to have been built by model_utils.find_all_models.
    dict_of_defaults : dict
        Dictionary of default parameters, assumed to have an outer key (generally used for layer types).
    inner_keys : Union[Tuple[str, ...], str] (optional)
        Optional list of keys to use after the outer key to find the correct parameter sub dict.


    Returns
    -------
    default_args_byfunc : Dict[str, Dict[str, Tuple[Any, ...]]]
        Dictionary of default constant arguments, broken up by function names and then by layer types.
    """

    default_args_byfunc = dict()

    # Make a copy of the default dict to ensure that no values are changed or that any mutable types are stored that
    #    could be later changed.
    dict_of_defaults = copy.deepcopy(dict_of_defaults)

    for func_name, const_arg_names in const_arg_dict.items():
        default_args_byfunc[func_name] = dict()

        # Assume outer key is used to pick out layer types or similar
        for outer_key, inner_dict in dict_of_defaults.items():

            if inner_keys is not None:
                # If parameters are stored in a subdict, then we need to find locate that subdict.
                inner_dict = nested_get(inner_dict, nested_keys=inner_keys, raiseon_nolocate=True)

            default_args = list()
            force_none = False
            for const_arg_name in const_arg_names:
                if const_arg_name in inner_dict:
                    default_args.append(inner_dict[const_arg_name])
                else:
                    # What to do if a required constant argument is not in the default dictionary.
                    # For now we will raise a warning.
                    warn(f'Can not find the required constant parameter {const_arg_name} (used in function {func_name})'
                         f' for layer type {outer_key}.')
                    # Then just put make this entire functions default args == None.
                    force_none = True

            if force_none:
                default_args_byfunc[func_name][outer_key] = None
            else:
                default_args_byfunc[func_name][outer_key] = tuple(default_args)

    return default_args_byfunc