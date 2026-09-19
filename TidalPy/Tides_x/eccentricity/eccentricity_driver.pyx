# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from TidalPy.Utilities_x.lookups cimport IntMap1, IntMap3, c_Key2, c_Key1, c_IntMap

def eccentricity_func(
        double eccentricity,
        int degree_l,
        object truncation = 3):
    """Eccentricity functions G_lpq(e) of one degree at one eccentricity.

    Parameters
    ----------
    eccentricity : float
        Orbital eccentricity, 0 <= e < 1.
    degree_l : int
        Tidal harmonic degree, 2 to 10.
    truncation : int or str, optional
        Truncation level n: every term of G_lpq through e^n is kept. Tabulated at 1, 2, 3, 4, 5, 10, 15, 20.

    Returns
    -------
    result_by_lpq : IntMap3
        Non-zero G_lpq keyed by (l, p, q).
    results_by_lp : dict
        The same values as an ``IntMap1`` of G by (q,) for each (l, p).
    """

    # Numeric strings such as "10" are accepted.
    if isinstance(truncation, str):
        try:
            truncation = int(truncation)
        except ValueError:
            raise NotImplementedError(
                "Unsupported truncation provided for eccentricity function. "
                "Tabulated levels: 1, 2, 3, 4, 5, 10, 15, 20.")
    elif isinstance(truncation, int):
        pass
    else:
        raise TypeError("Unexpected type found for `truncation`.")
    if truncation not in (1, 2, 3, 4, 5, 10, 15, 20):
            raise NotImplementedError(
                "Unsupported truncation provided for eccentricity function. "
                "Tabulated levels: 1, 2, 3, 4, 5, 10, 15, 20.")

    if degree_l not in (2, 3, 4, 5, 6, 7, 8, 9, 10):
        raise NotImplementedError(
            f"Degree l = {degree_l} is not currently supported for eccentricity function calculations. "
            "Supported degrees: l = 2 through 10.")
    
    cdef int error_code = 0
    cdef EccentricityFuncOutput result_pair = c_eccentricity_func(
        &error_code,
        eccentricity,
        degree_l,
        truncation)
    
    if error_code != 0:
        if error_code == -1:
            raise NotImplementedError("Eccentricity function error code -1: Unsupported / Not implemented truncation provided.")
        elif error_code == -2:
            raise NotImplementedError("Eccentricity function error code -2: Unsupported / Not implemented degree l provided.")
        else:
            raise RuntimeError(f"Unknown eccentricity function error code: {error_code}.")

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
