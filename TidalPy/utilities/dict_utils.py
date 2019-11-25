from typing import List, Tuple, Union, Any


PossibleKeyType = Union[str, int, float]


def nested_get(input_dict: dict,
               nested_keys: Union[PossibleKeyType, List[PossibleKeyType, ...], Tuple[PossibleKeyType, ...]],
               default: Any = None,
               raiseon_nolocate: bool = False) -> Any:
    """ Returns a value from a series of nested dictionaries given a list of keys

    Parameters
    ----------
    input_dict : dict
        Dictionary containing nested dictionary(ies)
    nested_keys : Union[PossibleKeyType, List[PossibleKeyType, ...], Tuple[PossibleKeyType, ...]]
        List of keys to search the dictionary with.
    default : Any = None
        Default value if key can't be found
    raiseon_nolocate : bool = False
        If key can't be found and this is True then a KeyError exception will be raised.

    Returns
    -------
    value : Any
        Value stored in the nested dictionaries
    """
    # TODO: When numba .44 comes out test this with @njit to speed it up. That version of numba should support
    #     dict types
    #     --Update: The new version of numba is out, but you need to specifiy the types of the keys and items.
    #         So I don't think njiting is going to work except for very specfic cases.

    if type(nested_keys) not in [str, int, float, list, tuple]:
        raise TypeError

    internal_dict_value = input_dict
    if type(nested_keys) not in [tuple, list]:
        nested_keys = [nested_keys]

    for key in nested_keys:
        if type(key) not in [str, int, float]:
            # Items in nested_key can not be lists or Tuples
            raise TypeError

        if key not in internal_dict_value:
            if raiseon_nolocate:
                raise KeyError
            else:
                return default
        else:
            internal_dict_value = internal_dict_value[key]

    return internal_dict_value
