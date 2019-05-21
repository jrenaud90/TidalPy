from typing import Union, List, Tuple

PossibleKeyType = Union[str, int, float]

def nested_get(input_dict: dict, nested_keys: Union[PossibleKeyType, List[PossibleKeyType], Tuple[PossibleKeyType]],
               default = None,
               raiseon_nolocate: bool = False):
    """ Returns a value from a series of nested dictionaries given a list of keys

    :param input_dict: <dict> Dictionary containing nested dictionary(ies)
    :param nested_keys: <list or tuple> List of keys to search the dictionary with.
    :param default:    Default value if key can't be found
    :param raiseon_nolocate: <bool> If key can't be found and this is True then a KeyError exception will be raised.
    :return:
    """
    # TODO: When numba .44 comes out test this with @njit to speed it up. That version of numba should support dict types

    if type(nested_keys) not in [str, int, float, list, tuple]:
        raise TypeError

    internal_dict_value = input_dict
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