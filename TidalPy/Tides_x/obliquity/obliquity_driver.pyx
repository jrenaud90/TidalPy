# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from TidalPy.Utilities_x.lookups cimport IntMap1, IntMap3, c_Key2, c_Key1, c_IntMap

def obliquity_func(
        double obliquity,
        int degree_l,
        object truncation = 'gen'):
    """Obliquity functions F_lmp(I) of one degree at one obliquity.

    Parameters
    ----------
    obliquity : float
        Obliquity I [rad].
    degree_l : int
        Tidal harmonic degree, 2 to 10.
    truncation : int or str, optional
        ``0`` or ``'off'`` (I = 0), ``1`` or ``2`` (every term through I^1 or I^2), or ``10`` or ``'gen'``
        (the exact form; default).

    Returns
    -------
    result_by_lmp : IntMap3
        Non-zero F_lmp keyed by (l, m, p).
    results_by_lm : dict
        The same values as an ``IntMap1`` of F by (p,) for each (l, m).
    """
    if isinstance(truncation, str):
        if truncation.lower() in ('gen', 'general'):
            truncation = 10
        elif truncation.lower() in ('off', '0', 'none'):
            truncation = 0
        elif truncation == '1':
            truncation = 1
        elif truncation == '2':
            truncation = 2
        else:
            raise NotImplementedError(
                "Unsupported truncation provided for obliquity function. "
                "Options are: ('general', 'off', '1', '2'). Use 'general' if unsure.")
    elif isinstance(truncation, int):
        if truncation not in (0, 1, 2, 10):
            raise NotImplementedError(
                "Unsupported truncation provided for obliquity function. "
                "Options are: ('general', 'off', '1', '2'). Use 'general' if unsure.")
    else:
        raise TypeError("Unexpected type found for `truncation`.")

    if degree_l not in (2, 3, 4, 5, 6, 7, 8, 9, 10):
        raise NotImplementedError(
            f"Degree l = {degree_l} is not currently supported for obliquity function calculations. "
            "Supported degrees: l = 2 through 10.")
    
    cdef int error_code = 0
    cdef ObliquityFuncOutput result_pair = c_obliquity_func(
        &error_code,
        obliquity,
        degree_l,
        truncation)
    
    if error_code != 0:
        if error_code == -1:
            raise NotImplementedError(
                "Obliquity function error code -1: Unsupported / Not implemented truncation provided.")
        elif error_code == -2:
            raise NotImplementedError(
                "Obliquity function error code -2: Unsupported / Not implemented degree l provided.")
        else:
            raise RuntimeError(f"Unknown obliquity function error code: {error_code}.")

    # Convert output to a python safe structure.
    cdef IntMap3 result_by_lmp = IntMap3()
    result_by_lmp.intmap_cinst = result_pair.first
    # The Python-accessible `IntMap` does not support non-numeric keys, so the results by (l, m) go into a
    # plain dict of inner IntMap1's.
    cdef dict results_by_lm = dict()
    cdef size_t i
    cdef pair[c_Key2, c_IntMap[c_Key1, double]] c_key_value
    cdef IntMap1 tmp_map
    for i in range(result_pair.second.size()):
        c_key_value = result_pair.second.data[i]
        tmp_map = IntMap1()
        tmp_map.intmap_cinst = c_key_value.second
        results_by_lm[(c_key_value.first.a, c_key_value.first.b)] = tmp_map

    return result_by_lmp, results_by_lm
