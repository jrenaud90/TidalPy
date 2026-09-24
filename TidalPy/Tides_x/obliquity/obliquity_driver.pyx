# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libcpp.pair cimport pair

from TidalPy.Utilities_x.lookups cimport IntMap1, IntMap3, c_Key2, c_Key1, c_IntMap
from TidalPy.Tides_x.obliquity.obliquity_common cimport (
    ObliquityFuncOutput, C_OBLIQUITY_TRUNCATIONS, C_NUM_OBLIQUITY_TRUNCATIONS, C_OBLIQUITY_OFF, C_OBLIQUITY_GENERAL,
    c_obliquity_accuracy_limit, c_recommend_obliquity_truncation)

import warnings

import TidalPy

# The tabulated obliquity truncation levels (the C++ list): level N keeps every product of two obliquity functions
# through I^N; level 0 is the obliquity off.
cdef list cy_truncations = []
cdef int cy_level_i
for cy_level_i in range(C_NUM_OBLIQUITY_TRUNCATIONS):
    cy_truncations.append(C_OBLIQUITY_TRUNCATIONS[cy_level_i])
OBLIQUITY_TRUNCATIONS = tuple(cy_truncations)

# Obliquity off (``"off"``), and the truncation code of the general (exact) obliquity functions (``"gen"`` or
# ``"general"`` wherever a truncation is given by name).
OBLIQUITY_OFF = C_OBLIQUITY_OFF
OBLIQUITY_GENERAL = C_OBLIQUITY_GENERAL

cdef dict cy_names = {"off": C_OBLIQUITY_OFF, "none": C_OBLIQUITY_OFF, "gen": C_OBLIQUITY_GENERAL,
                      "general": C_OBLIQUITY_GENERAL}


# Untabulated configured levels already warned about, so a stale configuration file warns once per session per level.
_WARNED_PROMOTIONS = set()


def promote_obliquity_truncation(object level, set warned_levels=None, object warn=None) -> int:
    """Resolve a configured obliquity truncation to a tabulated level.

    A configured level that is not tabulated (a configuration file written before the levels changed, say) is promoted
    to the next tabulated level, and anything past the highest to the general functions, with a once-per-session
    warning, so the file keeps working while the accuracy never silently decreases. The old levels map as 1 -> 2 and
    10 (the old general code) -> general. A level the caller passes directly is checked strictly instead
    (``validate_obliquity_truncation``).

    Parameters
    ----------
    level : int or str
        The configured level, or ``"off"`` / ``"gen"``.
    warned_levels : set, optional
        Levels already warned about; the module's own set when None.
    warn : bool, optional
        Whether to warn; the ``[warnings]`` ``truncation_promotion`` switch of the TidalPy configuration when None.

    Raises
    ------
    ValueError
        For a negative level or an unknown name.
    """
    if isinstance(level, str):
        text = level.strip().lower()
        if text in cy_names:
            return cy_names[text]
        try:
            level = int(text)
        except ValueError:
            raise ValueError(
                f"Obliquity truncation {level!r} is not supported. Supported levels: {OBLIQUITY_TRUNCATIONS}, 'off', "
                "or 'gen'.")
    level = int(level)
    if level in OBLIQUITY_TRUNCATIONS or level == C_OBLIQUITY_GENERAL:
        return level
    if level < 0:
        raise ValueError(
            f"Obliquity truncation {level} is not supported. Supported levels: {OBLIQUITY_TRUNCATIONS}, 'off', or "
            "'gen'.")
    promoted = C_OBLIQUITY_GENERAL
    for supported in OBLIQUITY_TRUNCATIONS:
        if supported > level:
            promoted = supported
            break
    if warned_levels is None:
        warned_levels = _WARNED_PROMOTIONS
    if warn is None:
        warn = bool(((getattr(TidalPy, "config_x", None) or {}).get("warnings", {}) or {}).get(
            "truncation_promotion", True))
    if warn and level not in warned_levels:
        warned_levels.add(level)
        warnings.warn(
            f"Obliquity truncation {level} is not tabulated; using {obliquity_truncation_name(promoted)!r} instead. "
            f"Supported levels: {OBLIQUITY_TRUNCATIONS}, 'off', or 'gen'.")
    return promoted


def validate_obliquity_truncation(object truncation=None) -> int:
    """Return an obliquity truncation as a tabulated int, or ``OBLIQUITY_GENERAL``.

    Parameters
    ----------
    truncation : int, str, or None, optional
        A level (an int or a numeric string such as ``"2"``), ``"off"`` (level 0), or ``"gen"`` / ``"general"``
        (``OBLIQUITY_GENERAL``). None takes the ``[tides]`` ``obliquity_trunc_lvl`` of the TidalPy configuration,
        promoted to the next tabulated level when the configuration holds an untabulated one
        (``promote_obliquity_truncation``).

    Raises
    ------
    TypeError
        For a value that is not an integer level or a name.
    NotImplementedError
        For a level passed directly that is not tabulated (see ``OBLIQUITY_TRUNCATIONS``).
    """
    if truncation is None:
        return promote_obliquity_truncation(
            ((getattr(TidalPy, "config_x", None) or {}).get("tides", {}) or {}).get("obliquity_trunc_lvl", "off"))
    if isinstance(truncation, bool):
        raise TypeError("An obliquity truncation is an integer level or a name, not a bool.")
    if isinstance(truncation, str):
        text = truncation.strip().lower()
        if text in cy_names:
            return cy_names[text]
        try:
            truncation = int(text)
        except ValueError:
            raise NotImplementedError(
                f"Obliquity truncation {truncation!r} is not tabulated. Tabulated levels: {OBLIQUITY_TRUNCATIONS}, "
                "'off', or 'gen'.")
    elif isinstance(truncation, float) and truncation.is_integer():
        truncation = int(truncation)
    try:
        level = int(truncation)
    except (TypeError, ValueError):
        raise TypeError(f"Unexpected obliquity truncation {truncation!r}.")
    if level != truncation:
        raise TypeError(f"Unexpected obliquity truncation {truncation!r}.")
    if level == C_OBLIQUITY_GENERAL:
        return level
    if level not in OBLIQUITY_TRUNCATIONS:
        raise NotImplementedError(
            f"Obliquity truncation {level} is not tabulated. Tabulated levels: {OBLIQUITY_TRUNCATIONS} (0 is off), or "
            "'gen' for the general functions.")
    return level


def obliquity_truncation_name(int truncation):
    """The level as written in a configuration: the int, or ``"gen"`` for ``OBLIQUITY_GENERAL``."""
    return "gen" if truncation == C_OBLIQUITY_GENERAL else truncation


def obliquity_accuracy_limit(object truncation, double tolerance=0.01, int max_degree_l=2) -> float:
    """The largest obliquity [rad] at which a level's heating stays within ``tolerance`` of the general value.

    Measured against the general functions, worst case over constant-phase-lag, constant-time-lag, and Maxwell tides at
    spin rates of 0.5, 1, and 2.3 times the mean motion and eccentricities of 0 and 0.2, for degree 2
    (``max_degree_l == 2``) or degree 3 (any higher ``max_degree_l``). The tabulated tolerances are 1e-8, 1e-6, 1e-4,
    1e-3, 1e-2, and 1e-1; a tolerance in between uses the next smaller one. Zero for ``"off"``, which holds only at
    zero obliquity, and NaN for ``"gen"``, which has no such limit.
    """
    cdef int level = validate_obliquity_truncation(truncation)
    return c_obliquity_accuracy_limit(level, tolerance, max_degree_l)


def recommend_obliquity_truncation(double obliquity, double tolerance=0.01, int max_degree_l=2):
    """The lowest obliquity truncation level that keeps the heating within ``tolerance`` at ``obliquity``.

    Parameters
    ----------
    obliquity : float
        Obliquity [rad].
    tolerance : float, optional
        Largest acceptable relative error in the heating (default 1%).
    max_degree_l : int, optional
        Highest tidal degree of the solve; degree 3 and higher lose accuracy at a lower obliquity than degree 2.

    Returns
    -------
    int or str
        0 (off) at zero obliquity, a tabulated level, or ``"gen"`` when no level holds the tolerance there (past about
        0.47 rad, 27 degrees, at the default 1%). The limits barely depend on the spin rate or the eccentricity.

    Raises
    ------
    ValueError
        For a tolerance that is not positive.
    """
    if not (tolerance > 0.0):
        raise ValueError(f"The tolerance must be positive; got {tolerance}.")
    return obliquity_truncation_name(c_recommend_obliquity_truncation(obliquity, tolerance, max_degree_l))


cdef tuple cy_obliquity_output(ObliquityFuncOutput& result_pair):
    """Convert the C++ maps into an IntMap3 by (l, m, p) and a dict of IntMap1 by (p,) for each (l, m)."""
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


cdef void cy_check_degree(int degree_l) except *:
    if degree_l not in (2, 3, 4, 5, 6, 7, 8, 9, 10):
        raise NotImplementedError(
            f"Degree l = {degree_l} is not currently supported for obliquity function calculations. "
            "Supported degrees: l = 2 through 10.")


cdef void cy_check_error(int error_code) except *:
    if error_code == -1:
        raise NotImplementedError("Obliquity function error code -1: the truncation level is not tabulated.")
    elif error_code == -2:
        raise NotImplementedError("Obliquity function error code -2: the degree l is not supported.")
    elif error_code != 0:
        raise RuntimeError(f"Unknown obliquity function error code: {error_code}.")


def obliquity_func(
        double obliquity,
        int degree_l,
        object truncation = 'gen'):
    """Obliquity functions F_lmp(I) of one degree at one obliquity, unsquared.

    Truncation level N keeps every function whose Taylor series starts at or below I^N, each through I^N, so a tidal
    potential built from them is complete through I^N. ``"gen"`` is exact at any obliquity.

    Parameters
    ----------
    obliquity : float
        Obliquity I [rad].
    degree_l : int
        Tidal harmonic degree, 2 to 10.
    truncation : int or str, optional
        ``"off"`` or 0 (I = 0), 2 or 4 (level N), or ``"gen"`` / ``"general"`` (the exact functions; default). None
        takes the ``[tides]`` ``obliquity_trunc_lvl`` of the TidalPy configuration.

    Returns
    -------
    result_by_lmp : IntMap3
        Non-zero F_lmp keyed by (l, m, p).
    results_by_lm : dict
        The same values as an ``IntMap1`` of F by (p,) for each (l, m).
    """
    cdef int level = validate_obliquity_truncation(truncation)
    cy_check_degree(degree_l)
    cdef int error_code = 0
    cdef ObliquityFuncOutput result_pair = c_obliquity_func(&error_code, obliquity, degree_l, level)
    cy_check_error(error_code)
    return cy_obliquity_output(result_pair)


def obliquity_squared_func(
        double obliquity,
        int degree_l,
        object truncation = 'gen'):
    """Squared obliquity functions F_lmp(I)^2 of one degree at one obliquity, cut at I^N.

    These are what the global (1D) heating uses: every function starting at or below I^(N / 2), each square through
    I^N, so the sum over modes is the Taylor series of the heating through I^N. ``"gen"`` gives the plain squares of
    the exact functions.

    Parameters
    ----------
    obliquity : float
        Obliquity I [rad].
    degree_l : int
        Tidal harmonic degree, 2 to 10.
    truncation : int or str, optional
        As in ``obliquity_func``.

    Returns
    -------
    result_by_lmp : IntMap3
        Non-zero cut F_lmp^2 keyed by (l, m, p).
    results_by_lm : dict
        The same values as an ``IntMap1`` by (p,) for each (l, m).
    """
    cdef int level = validate_obliquity_truncation(truncation)
    cy_check_degree(degree_l)
    cdef int error_code = 0
    cdef ObliquityFuncOutput result_pair = c_obliquity_squared_func(&error_code, obliquity, degree_l, level)
    cy_check_error(error_code)
    return cy_obliquity_output(result_pair)
