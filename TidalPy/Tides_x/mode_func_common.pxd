# Shared by the eccentricity and obliquity drivers: the lookup maps their C++ builders return and the Python objects
# made from them, the builders' error codes, the TidalPy configuration reads, and the truncation parsing. Every
# function here is inline, so no extension module backs this file.

from libcpp.pair cimport pair
from libcpp.utility cimport move

from TidalPy.Utilities_x.lookups cimport IntMap1, IntMap3, c_IntMap, c_Key1, c_Key2, c_Key3


# Declared in mode_func_common_.hpp, which each driver's own common header includes.
cdef extern from * nogil:
    ctypedef pair[c_IntMap[c_Key3, double], c_IntMap[c_Key2, c_IntMap[c_Key1, double]]] c_ModeFuncOutput


cdef inline tuple cy_mode_func_output(c_ModeFuncOutput& result_pair):
    """Move the C++ maps into an IntMap3 by the full key and a dict of IntMap1 by the last key part for each leading
    (l, row) pair; the maps are left empty."""
    cdef IntMap3 result_by_key3 = IntMap3()
    result_by_key3.intmap_cinst = move(result_pair.first)
    # The Python-accessible `IntMap` does not support non-numeric keys, so the results by (l, row) go into a plain dict
    # of inner IntMap1's.
    cdef dict results_by_key2 = dict()
    cdef size_t i
    cdef pair[c_Key2, c_IntMap[c_Key1, double]]* entry_ptr
    cdef IntMap1 inner_map
    for i in range(result_pair.second.size()):
        entry_ptr = &result_pair.second.data[i]
        inner_map = IntMap1()
        inner_map.intmap_cinst = move(entry_ptr.second)
        results_by_key2[(entry_ptr.first.a, entry_ptr.first.b)] = inner_map
    return result_by_key3, results_by_key2


cdef inline void cy_check_error(int error_code, str kind) except *:
    """Raise for a non-zero error code of an eccentricity or obliquity builder; `kind` names which."""
    if error_code == -1:
        raise NotImplementedError(f"{kind} function error code -1: the truncation level is not tabulated.")
    elif error_code == -2:
        raise NotImplementedError(
            f"{kind} function error code -2: the degree l is not supported (l = 2 through 10).")
    elif error_code != 0:
        raise RuntimeError(f"Unknown {kind.lower()} function error code: {error_code}.")


cdef inline object cy_config_value(str section, str key, object fallback):
    """A ``[section]`` value of the TidalPy configuration, or `fallback` when it (or the configuration) is absent."""
    import TidalPy
    return ((getattr(TidalPy, "config_x", None) or {}).get(section, {}) or {}).get(key, fallback)


cdef inline object cy_validate_truncation(object truncation, str kind, tuple levels, dict names, str names_hint):
    """The integer truncation level of `truncation` given directly: a name of `names` (matched case-insensitively), an
    int, an integral float, or a numeric string, which must be one of `levels` or a code of `names`.

    `kind` ("eccentricity" or "obliquity") and `names_hint` (the names as listed after the levels) go into the messages.

    Raises
    ------
    TypeError
        For a bool, or a value that is not an integer level.
    NotImplementedError
        For a string that is neither a name nor a number, or a level that is not tabulated.
    """
    if isinstance(truncation, bool):
        raise TypeError(f"An {kind} truncation is an integer level or a name, not a bool.")
    if isinstance(truncation, str):
        text = truncation.strip().lower()
        if text in names:
            return names[text]
        try:
            truncation = int(text)
        except ValueError:
            raise NotImplementedError(
                f"{kind.capitalize()} truncation {truncation!r} is not tabulated. Tabulated levels: {levels}"
                f"{names_hint}.")
    elif isinstance(truncation, float) and truncation.is_integer():
        truncation = int(truncation)
    try:
        level = int(truncation)
    except (TypeError, ValueError):
        raise TypeError(f"Unexpected {kind} truncation {truncation!r}.")
    if level != truncation:
        raise TypeError(f"Unexpected {kind} truncation {truncation!r}.")
    if (level not in levels) and (level not in names.values()):
        raise NotImplementedError(
            f"{kind.capitalize()} truncation {level} is not tabulated. Tabulated levels: {levels}{names_hint}.")
    return level
