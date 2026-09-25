# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from TidalPy.Tides_x.mode_func_common cimport (
    c_ModeFuncOutput, cy_check_error, cy_config_value, cy_mode_func_output, cy_validate_truncation)
from TidalPy.Tides_x.eccentricity.eccentricity_common cimport (
    C_ECCENTRICITY_TRUNCATIONS, C_NUM_ECCENTRICITY_TRUNCATIONS, C_ECCENTRICITY_EXACT, C_ECCENTRICITY_EXACT_TOLERANCE,
    c_eccentricity_accuracy_limit, c_recommend_eccentricity_truncation)

import warnings

# The tabulated eccentricity truncation levels (the C++ list): level N keeps every product of two eccentricity
# functions through e^N.
cdef list cy_truncations = []
cdef int cy_level_i
for cy_level_i in range(C_NUM_ECCENTRICITY_TRUNCATIONS):
    cy_truncations.append(C_ECCENTRICITY_TRUNCATIONS[cy_level_i])
ECCENTRICITY_TRUNCATIONS = tuple(cy_truncations)

# The truncation code of the exact eccentricity functions (Hansen coefficients from the exact Kepler orbit, the modes
# chosen by the exact tolerance); ``"exact"`` wherever a truncation is given by name.
ECCENTRICITY_EXACT = C_ECCENTRICITY_EXACT

cdef dict cy_names = {"exact": C_ECCENTRICITY_EXACT}


def validate_eccentricity_exact_tolerance(object tolerance=None) -> float:
    """Return the heating tail tolerance of the exact eccentricity functions, in (0, 1).

    None takes the ``[tides]`` ``eccentricity_exact_tolerance`` of the TidalPy configuration.

    Raises
    ------
    ValueError
        For a tolerance outside (0, 1).
    """
    if tolerance is None:
        tolerance = cy_config_value("tides", "eccentricity_exact_tolerance", C_ECCENTRICITY_EXACT_TOLERANCE)
    if isinstance(tolerance, bool):
        raise ValueError("The exact eccentricity tolerance is a number in (0, 1), not a bool.")
    value = float(tolerance)
    if not (0.0 < value < 1.0):
        raise ValueError(f"The exact eccentricity tolerance must be in (0, 1); got {tolerance!r}.")
    return value


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
    if isinstance(level, str) and level.strip().lower() == "exact":
        return C_ECCENTRICITY_EXACT
    level = int(level)
    if level in ECCENTRICITY_TRUNCATIONS or level == C_ECCENTRICITY_EXACT:
        return level
    if warned_levels is None:
        warned_levels = _WARNED_PROMOTIONS
    if warn is None:
        warn = bool(cy_config_value("warnings", "truncation_promotion", True))
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
        A level (an int or a numeric string such as ``"10"``) or ``"exact"`` (``ECCENTRICITY_EXACT``). None takes the
        ``[tides]`` ``eccentricity_trunc_lvl`` of the TidalPy configuration, promoted to the next tabulated level when
        the configuration holds an untabulated one (``promote_eccentricity_truncation``).

    Raises
    ------
    TypeError
        For a value that is not an integer level.
    NotImplementedError
        For a level passed directly that is not tabulated (see ``ECCENTRICITY_TRUNCATIONS``).
    """
    if truncation is None:
        return promote_eccentricity_truncation(cy_config_value("tides", "eccentricity_trunc_lvl", 10))
    return cy_validate_truncation(truncation, "eccentricity", ECCENTRICITY_TRUNCATIONS, cy_names, ", or 'exact'")


def eccentricity_truncation_name(int truncation):
    """The level as written in a configuration: the int, or ``"exact"`` for ``ECCENTRICITY_EXACT``."""
    return "exact" if truncation == C_ECCENTRICITY_EXACT else truncation


def eccentricity_accuracy_limit(object truncation, double tolerance=0.01, int max_degree_l=2) -> float:
    """The largest eccentricity at which a truncation level's heating stays within ``tolerance`` of the exact value.

    Measured against the exact heating, worst case over constant-phase-lag, constant-time-lag, and Maxwell tides at spin
    rates of 0.5, 1, and 2.3 times the mean motion, for degree 2 (``max_degree_l == 2``) or degree 3 (any higher
    ``max_degree_l``). The tabulated tolerances are 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, and 1e-1; a tolerance in between
    uses the next smaller one. NaN for ``"exact"``, which has no such limit.
    """
    cdef int level = validate_eccentricity_truncation(truncation)
    return c_eccentricity_accuracy_limit(level, tolerance, max_degree_l)


def recommend_eccentricity_truncation(double eccentricity, double tolerance=0.01, int max_degree_l=2):
    """The lowest eccentricity truncation level that keeps the heating within ``tolerance`` at ``eccentricity``.

    Parameters
    ----------
    eccentricity : float
        Orbital eccentricity, 0 <= e < 1.
    tolerance : float, optional
        Largest acceptable relative error in the heating (default 1%).
    max_degree_l : int, optional
        Highest tidal degree of the solve; degree 3 and higher lose accuracy at a lower eccentricity than degree 2.

    Returns
    -------
    int or str
        A tabulated level, or ``"exact"`` when no level holds the tolerance there (eccentricities past about 0.75, or a
        tolerance tighter than the tables were measured to). The limits hold for any spin rate measured; they barely
        depend on it.

    Raises
    ------
    ValueError
        For an eccentricity outside [0, 1) or a tolerance that is not positive.
    """
    if not (0.0 <= eccentricity < 1.0):
        raise ValueError(f"The eccentricity must be in [0, 1); got {eccentricity}.")
    if not (tolerance > 0.0):
        raise ValueError(f"The tolerance must be positive; got {tolerance}.")
    return eccentricity_truncation_name(c_recommend_eccentricity_truncation(eccentricity, tolerance, max_degree_l))


def eccentricity_func(
        double eccentricity,
        int degree_l,
        object truncation = None,
        object exact_tolerance = None):
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
        Truncation level N, one of ``ECCENTRICITY_TRUNCATIONS``, or ``"exact"`` for the functions from the exact orbit.
        None takes the ``[tides]`` ``eccentricity_trunc_lvl`` of the TidalPy configuration.
    exact_tolerance : float, optional
        For ``"exact"``: the modes kept are those whose q^2-weighted squares leave a tail below this fraction of the
        total (so it bounds the relative error of the synchronous constant-time-lag heating). None takes the
        ``[tides]`` ``eccentricity_exact_tolerance`` (1e-4 by default).

    Returns
    -------
    result_by_lpq : IntMap3
        Non-zero G_lpq keyed by (l, p, q).
    results_by_lp : dict
        The same values as an ``IntMap1`` of G by (q,) for each (l, p).
    """
    cdef int level = validate_eccentricity_truncation(truncation)
    cdef double tolerance = validate_eccentricity_exact_tolerance(exact_tolerance)
    cdef int error_code = 0
    cdef c_ModeFuncOutput result_pair = c_eccentricity_func(&error_code, eccentricity, degree_l, level, tolerance)
    cy_check_error(error_code, "Eccentricity")
    return cy_mode_func_output(result_pair)


def eccentricity_squared_func(
        double eccentricity,
        int degree_l,
        object truncation = None,
        object exact_tolerance = None):
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
        Truncation level N, one of ``ECCENTRICITY_TRUNCATIONS``, or ``"exact"`` for the plain squares of the functions
        from the exact orbit. None takes the ``[tides]`` ``eccentricity_trunc_lvl`` of the TidalPy configuration.
    exact_tolerance : float, optional
        For ``"exact"``, as in ``eccentricity_func``.

    Returns
    -------
    result_by_lpq : IntMap3
        Non-zero cut G_lpq^2 keyed by (l, p, q).
    results_by_lp : dict
        The same values as an ``IntMap1`` by (q,) for each (l, p).
    """
    cdef int level = validate_eccentricity_truncation(truncation)
    cdef double tolerance = validate_eccentricity_exact_tolerance(exact_tolerance)
    cdef int error_code = 0
    cdef c_ModeFuncOutput result_pair = c_eccentricity_squared_func(
        &error_code, eccentricity, degree_l, level, tolerance)
    cy_check_error(error_code, "Eccentricity")
    return cy_mode_func_output(result_pair)
