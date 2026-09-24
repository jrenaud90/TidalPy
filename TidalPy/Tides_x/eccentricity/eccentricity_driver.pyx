# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libcpp.pair cimport pair

from TidalPy.Utilities_x.lookups cimport IntMap1, IntMap3, c_Key2, c_Key1, c_IntMap
from TidalPy.Tides_x.eccentricity.eccentricity_common cimport (
    EccentricityFuncOutput, C_ECCENTRICITY_TRUNCATIONS, C_NUM_ECCENTRICITY_TRUNCATIONS)

import warnings

import TidalPy

# The tabulated eccentricity truncation levels (the C++ list): level N keeps every product of two eccentricity
# functions through e^N.
cdef list cy_truncations = []
cdef int cy_level_i
for cy_level_i in range(C_NUM_ECCENTRICITY_TRUNCATIONS):
    cy_truncations.append(C_ECCENTRICITY_TRUNCATIONS[cy_level_i])
ECCENTRICITY_TRUNCATIONS = tuple(cy_truncations)


# Untabulated configured levels already warned about, so a stale configuration file warns once per session per level.
_WARNED_PROMOTIONS = set()


def promote_eccentricity_truncation(object level, set warned_levels=None, object warn=None) -> int:
    """Resolve a configured eccentricity truncation to a tabulated level.

    A configured level that is not tabulated (a configuration file written before the levels changed, say) is promoted
    to the next tabulated level with a once-per-session warning, so the file keeps working while the accuracy never
    silently decreases. A level the caller passes directly is checked strictly instead
    (``validate_eccentricity_truncation``).

    Parameters
    ----------
    level : int
        The configured level.
    warned_levels : set, optional
        Levels already warned about; the module's own set when None.
    warn : bool, optional
        Whether to warn; the ``[warnings]`` ``truncation_promotion`` switch of the TidalPy configuration when None.

    Raises
    ------
    ValueError
        For a level above every tabulated one.
    """
    level = int(level)
    if level in ECCENTRICITY_TRUNCATIONS:
        return level
    if warned_levels is None:
        warned_levels = _WARNED_PROMOTIONS
    if warn is None:
        warn = bool(((getattr(TidalPy, "config_x", None) or {}).get("warnings", {}) or {}).get(
            "truncation_promotion", True))
    for supported in ECCENTRICITY_TRUNCATIONS:
        if supported > level:
            if warn and level not in warned_levels:
                warned_levels.add(level)
                warnings.warn(
                    f"Eccentricity truncation {level} is not tabulated; using {supported} instead. "
                    f"Supported levels: {ECCENTRICITY_TRUNCATIONS}.")
            return supported
    raise ValueError(
        f"Eccentricity truncation {level} is not supported. Supported levels: {ECCENTRICITY_TRUNCATIONS}.")


def validate_eccentricity_truncation(object truncation=None) -> int:
    """Return an eccentricity truncation level as a tabulated int.

    Parameters
    ----------
    truncation : int, str, or None, optional
        A level (an int or a numeric string such as ``"10"``). None takes the ``[tides]`` ``eccentricity_trunc_lvl``
        of the TidalPy configuration, promoted to the next tabulated level when the configuration holds an untabulated
        one (``promote_eccentricity_truncation``).

    Raises
    ------
    TypeError
        For a value that is not an integer level.
    NotImplementedError
        For a level passed directly that is not tabulated (see ``ECCENTRICITY_TRUNCATIONS``).
    """
    if truncation is None:
        return promote_eccentricity_truncation(
            ((getattr(TidalPy, "config_x", None) or {}).get("tides", {}) or {}).get("eccentricity_trunc_lvl", 10))
    if isinstance(truncation, bool):
        raise TypeError("An eccentricity truncation is an integer level, not a bool.")
    if isinstance(truncation, str):
        try:
            truncation = int(truncation)
        except ValueError:
            raise NotImplementedError(
                f"Eccentricity truncation {truncation!r} is not tabulated. "
                f"Tabulated levels: {ECCENTRICITY_TRUNCATIONS}.")
    elif isinstance(truncation, float) and truncation.is_integer():
        truncation = int(truncation)
    try:
        level = int(truncation)
    except (TypeError, ValueError):
        raise TypeError(f"Unexpected eccentricity truncation {truncation!r}.")
    if level != truncation:
        raise TypeError(f"Unexpected eccentricity truncation {truncation!r}.")
    if level not in ECCENTRICITY_TRUNCATIONS:
        raise NotImplementedError(
            f"Eccentricity truncation {level} is not tabulated. Tabulated levels: {ECCENTRICITY_TRUNCATIONS}.")
    return level


cdef tuple cy_eccentricity_output(EccentricityFuncOutput& result_pair):
    """Convert the C++ maps into an IntMap3 by (l, p, q) and a dict of IntMap1 by (q,) for each (l, p)."""
    cdef IntMap3 result_by_lpq = IntMap3()
    result_by_lpq.intmap_cinst = result_pair.first
    # The Python-accessible `IntMap` does not support non-numeric keys, so the results by (l, p) go into a
    # plain dict of inner IntMap1's.
    cdef dict results_by_lp = dict()
    cdef size_t i
    cdef pair[c_Key2, c_IntMap[c_Key1, double]] c_key_value
    cdef IntMap1 tmp_map
    for i in range(result_pair.second.size()):
        c_key_value = result_pair.second.data[i]
        tmp_map = IntMap1()
        tmp_map.intmap_cinst = c_key_value.second
        results_by_lp[(c_key_value.first.a, c_key_value.first.b)] = tmp_map
    return result_by_lpq, results_by_lp


cdef void cy_check_degree(int degree_l) except *:
    if degree_l not in (2, 3, 4, 5, 6, 7, 8, 9, 10):
        raise NotImplementedError(
            f"Degree l = {degree_l} is not currently supported for eccentricity function calculations. "
            "Supported degrees: l = 2 through 10.")


cdef void cy_check_error(int error_code) except *:
    if error_code == -1:
        raise NotImplementedError("Eccentricity function error code -1: the truncation level is not tabulated.")
    elif error_code == -2:
        raise NotImplementedError("Eccentricity function error code -2: the degree l is not supported.")
    elif error_code != 0:
        raise RuntimeError(f"Unknown eccentricity function error code: {error_code}.")


def eccentricity_func(
        double eccentricity,
        int degree_l,
        object truncation = None):
    """Eccentricity functions G_lpq(e) of one degree at one eccentricity, unsquared.

    Truncation level N keeps every mode with |q| <= N, each G_lpq through e^N, so a tidal potential built from them is
    complete through e^N. The k = l - 2p + q = 0 modes are exact.

    Parameters
    ----------
    eccentricity : float
        Orbital eccentricity, 0 <= e < 1.
    degree_l : int
        Tidal harmonic degree, 2 to 10.
    truncation : int or str, optional
        Truncation level N, one of ``ECCENTRICITY_TRUNCATIONS``. None takes the ``[tides]`` ``eccentricity_trunc_lvl``
        of the TidalPy configuration.

    Returns
    -------
    result_by_lpq : IntMap3
        Non-zero G_lpq keyed by (l, p, q).
    results_by_lp : dict
        The same values as an ``IntMap1`` of G by (q,) for each (l, p).
    """
    cdef int level = validate_eccentricity_truncation(truncation)
    cy_check_degree(degree_l)
    cdef int error_code = 0
    cdef EccentricityFuncOutput result_pair = c_eccentricity_func(&error_code, eccentricity, degree_l, level)
    cy_check_error(error_code)
    return cy_eccentricity_output(result_pair)


def eccentricity_squared_func(
        double eccentricity,
        int degree_l,
        object truncation = None):
    """Squared eccentricity functions G_lpq(e)^2 of one degree at one eccentricity, cut at e^N.

    These are what the global (1D) heating uses: every mode with |q| <= N / 2, each square through e^N, so the sum
    over modes is the Taylor series of the heating through e^N. A cut square can be negative for the highest-|q|
    modes at large e. The k = l - 2p + q = 0 modes are exact.

    Parameters
    ----------
    eccentricity : float
        Orbital eccentricity, 0 <= e < 1.
    degree_l : int
        Tidal harmonic degree, 2 to 10.
    truncation : int or str, optional
        Truncation level N, one of ``ECCENTRICITY_TRUNCATIONS``. None takes the ``[tides]`` ``eccentricity_trunc_lvl``
        of the TidalPy configuration.

    Returns
    -------
    result_by_lpq : IntMap3
        Non-zero cut G_lpq^2 keyed by (l, p, q).
    results_by_lp : dict
        The same values as an ``IntMap1`` by (q,) for each (l, p).
    """
    cdef int level = validate_eccentricity_truncation(truncation)
    cy_check_degree(degree_l)
    cdef int error_code = 0
    cdef EccentricityFuncOutput result_pair = c_eccentricity_squared_func(&error_code, eccentricity, degree_l, level)
    cy_check_error(error_code)
    return cy_eccentricity_output(result_pair)
