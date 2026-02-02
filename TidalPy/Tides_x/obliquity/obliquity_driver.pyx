# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from TidalPy.utilities.lookups cimport IntMap1, IntMap3, c_Key2, c_Key1, c_IntMap

def obliquity_func(
        double obliquity,
        int degree_l,
        object truncation = 'gen'):
    
    # Clean up input and check for issues.
    if isinstance(truncation, str):
        if truncation.lower() in ('gen', 'general'):
            truncation = 10
        elif truncation.lower() in ('off', '0', 'none'):
            truncation = 0
        elif truncation == '2':
            truncation = 2
        elif truncation == '4':
            truncation = 4
        else:
            raise NotImplementedError("Unsupported truncation provided for obliquity function. Options are: ('general', 'off', '2', '4'). Use 'general' if unsure.")
    elif isinstance(truncation, int):
        if truncation not in (0, 2, 4, 10):
            raise NotImplementedError("Unsupported truncation provided for obliquity function. Options are: ('general', 'off', '2', '4'). Use 'general' if unsure.")
    else:
        raise TypeError("Unexpected type found for `truncation`.")
    
    if degree_l not in (2,):
        raise NotImplementedError(f"Degree l = {degree_l} is not currently supported for obliquity function calculations. Options are l = 2")
    
    # Call c function to get result.
    cdef int error_code = 0
    cdef ObliquityFuncOutput result_pair = c_obliquity_func(
        &error_code,
        obliquity,
        degree_l,
        truncation)
    
    # Check for issues during calculation.
    if error_code != 0:
        if error_code == -1:
            raise NotImplementedError("Obliquity function error code -1: Unsupported / Not implemented truncation provided.")
        elif error_code == -2:
            raise NotImplementedError("Obliquity function error code -2: Unsupported / Not implemented degree l provided.")
        else:
            raise RuntimeError(f"Unknown obliquity function error code: {error_code}.")

    # Convert output to a python safe structure.
    cdef IntMap3 result_by_lmp = IntMap3()
    result_by_lmp.intmap_cinst = result_pair.first
    # For the results by l, m - the Python-accessible `IntMap` does not (currently) support non-numeric types.
    # Instead we will just make a regular python dictionary and store the inner IntMap1's in it.
    cdef dict results_by_lm = dict()
    cdef size_t i
    cdef pair[c_Key2, c_IntMap[c_Key1, double]] c_key_value
    cdef IntMap1 tmp_map
    for i in range(result_pair.second.size()):
        c_key_value = result_pair.second.data[i]
        # Store the IntMap as a Python tuple in the dict.
        tmp_map = IntMap1()
        tmp_map.intmap_cinst = c_key_value.second
        results_by_lm[(c_key_value.first.a, c_key_value.first.b)] = tmp_map

    return result_by_lmp, results_by_lm
