# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Standalone global (1D) tidal-mode collapse.

``collapse_global_tides`` runs the global-potential engine for an orbital and spin state, then collapses
the per-mode potential terms with an analytic tide model's -Im[k_l] to give the global tidal heating and
the three orbital potential derivatives. Only the analytic models are supported here: the rheology model
needs per-mode Love numbers from the radial solver and is driven by the world's ``calc_tides``.
"""

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr

from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.Tides_x.classes.tide import TIDE_CONFIG_KEYS, _same_model
from TidalPy.Utilities_x.classes_x.classes import check_config_keys, factory_defaults
import TidalPy
from TidalPy.Tides_x.eccentricity.eccentricity_driver import (
    validate_eccentricity_exact_tolerance, validate_eccentricity_truncation)
from TidalPy.Tides_x.classes.tide cimport (
    c_TideBase, c_TideModel, c_TideModelConfig, c_tide_model_from_name, c_find_tide,
)

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


cdef extern from "global_.hpp" nogil:

    cdef cppclass c_GlobalPotentialStorage:
        int error_code
        int working_on_l

    c_GlobalPotentialStorage c_global_potential(
        double planet_radius,
        double semi_major_axis,
        double orbital_frequency,
        double spin_frequency,
        double obliquity,
        double eccentricity,
        double host_mass,
        double G_to_use,
        int min_degree_l,
        int max_degree_l,
        int obliquity_truncation,
        int eccentricity_truncation,
        double eccentricity_exact_tolerance) except +


cdef extern from "tide_collapse_.hpp" nogil:

    cdef cppclass c_GlobalTideResult:
        double tidal_heating
        double dU_dM
        double dU_dw
        double dU_dO
        double dU_dM_minus_dw
        int num_modes
        int error_code

    # The C++ function takes a third, defaulted pointer to per-mode radial-solver Love numbers. Only a
    # world's rheology path passes it, and that call is made in C++, so it is left out here and the
    # analytic collapse below always takes the null default.
    c_GlobalTideResult c_collapse_global_tides(
        const c_GlobalPotentialStorage& potential,
        const c_TideBase& tide_model) except +


cdef int cy_resolve_obliquity_truncation(object obliquity_truncation) except? -999:
    """Normalize an obliquity truncation ('gen', 'off', or an int) to the C++ integer."""
    cdef int value = 0
    cdef str text
    if isinstance(obliquity_truncation, str):
        text = obliquity_truncation.lower()
        if text in ("gen", "general"):
            value = 10
        elif text in ("off",):
            value = 0
        else:
            try:
                value = int(obliquity_truncation)
            except ValueError:
                raise ValueError("Unexpected obliquity truncation encountered.")
    elif isinstance(obliquity_truncation, bool):
        raise ValueError("An obliquity truncation is 'off', 'gen', or an integer level, not a bool.")
    elif isinstance(obliquity_truncation, (int, float)) and float(obliquity_truncation).is_integer():
        # A whole-valued float (2.0) is the level it names, rather than falling through to 'off'.
        value = int(obliquity_truncation)
    else:
        raise ValueError(f"Unexpected obliquity truncation {obliquity_truncation!r}.")
    if value not in (0, 1, 2, 10):
        raise NotImplementedError(
            f"Obliquity truncation {value} is not tabulated. "
            "Supported levels: 0 ('off'), 1, 2, 10 ('gen', fully general).")
    return value


cdef c_TideModelConfig cy_build_tide_config(dict config) except *:
    """Build a c_TideModelConfig from the optional per-degree list keys, indexed from l = 2."""
    cdef c_TideModelConfig cfg
    if config is None:
        return cfg
    if "fixed_k" in config:
        for value in config["fixed_k"]:
            cfg.fixed_k.push_back(<double>value)
    if "fixed_q" in config:
        for value in config["fixed_q"]:
            cfg.fixed_q.push_back(<double>value)
    if "fixed_dt_s" in config:
        for value in config["fixed_dt_s"]:
            cfg.fixed_dt.push_back(<double>value)
    return cfg


def collapse_global_tides(
        double planet_radius,
        double orbital_frequency,
        double spin_frequency,
        double eccentricity,
        double obliquity,
        double semi_major_axis,
        double host_mass,
        double G_to_use,
        str tide_model,
        dict tide_config=None,
        int min_degree_l=2,
        int max_degree_l=2,
        object eccentricity_truncation=None,
        object eccentricity_exact_tolerance=None,
        object obliquity_truncation=None) -> dict:
    """Collapse the global tidal modes into heating and orbital potential derivatives.

    The body radius comes first, then the orbital state in the same order as the world's
    ``calc_tides`` (orbital frequency, spin frequency, eccentricity, obliquity, semi-major axis,
    host mass), then Newton's constant.

    Parameters
    ----------
    planet_radius : float
        Radius of the tidally deformed body [m].
    orbital_frequency : float
        Orbital mean motion [rad s-1].
    spin_frequency : float
        Spin rate of the deformed body [rad s-1].
    eccentricity : float
        Orbital eccentricity [dimensionless].
    obliquity : float
        Axial tilt [radians].
    semi_major_axis : float
        Orbital semi-major axis [m].
    host_mass : float
        Mass of the tidal host [kg].
    G_to_use : float
        Gravitational constant to use [m3 kg-1 s-2].
    tide_model : str
        Analytic tide model name: 
            - ``"cpl"``/``"fixed_q"``
            - ``"ctl"``/``"fixed_dt"``
            - ``"ctl_q"``/``"fixed_dt_q"``
        The ``"rheology"`` model is not supported here (use the world's ``calc_tides``).
    tide_config : dict, optional
        Per-degree model parameters (``fixed_k``, ``fixed_q``, ``fixed_dt_s`` [s] lists indexed
        from degree l = 2). A key left out takes the ``[tides]`` default of the TidalPy configuration, as
        ``make_tide`` does; any other key raises ``ValueError``.
    min_degree_l, max_degree_l : int
        Tidal harmonic degree range (2..10).
    eccentricity_truncation : int, optional
        Eccentricity truncation level N, one of ``ECCENTRICITY_TRUNCATIONS`` (every product of two eccentricity
        functions through e^N), or ``"exact"``. None takes the ``[tides]`` ``eccentricity_trunc_lvl`` of the TidalPy
        configuration.
    eccentricity_exact_tolerance : float, optional
        Heating tail tolerance of the ``"exact"`` functions; None takes the ``[tides]`` value.
    obliquity_truncation : str or int, optional
        Obliquity truncation: ``"off"`` (0), 1, 2, or ``"gen"``/``"general"`` (10). None takes the ``[tides]``
        ``obliquity_trunc_lvl`` of the TidalPy configuration (``"off"`` by default, which ignores the obliquity).

    Returns
    -------
    dict
        ``{"tidal_heating": W, "dUdM": ..., "dUdw": ..., "dUdO": ..., "num_modes": int}``.

    Raises
    ------
    ValueError
        If the tide model name or a config key is unknown, the degree range is out of order, the eccentricity is
        outside [0, 1), or the semi-major axis is not positive.
    NotImplementedError
        If the rheology model is requested, or a truncation/degree is unsupported.
    """
    if obliquity_truncation is None:
        # The same default a built world takes: the [tides] obliquity_trunc_lvl of TidalPy.config_x.
        obliquity_truncation = ((getattr(TidalPy, "config_x", None) or {}).get("tides", {}) or {}).get(
            "obliquity_trunc_lvl", "off")
    cdef int i_obliquity_truncation = cy_resolve_obliquity_truncation(obliquity_truncation)

    if not (2 <= min_degree_l <= max_degree_l <= 10):
        raise ValueError(
            f"The degree range must satisfy 2 <= min_degree_l <= max_degree_l <= 10; got {min_degree_l} to "
            f"{max_degree_l}.")
    if not (0.0 <= eccentricity < 1.0):
        raise ValueError(f"The eccentricity must be in [0, 1); got {eccentricity}.")
    if not (semi_major_axis > 0.0):
        raise ValueError(f"The semi-major axis must be positive; got {semi_major_axis} m.")
    # None takes the [tides] eccentricity_trunc_lvl of the TidalPy configuration, as a built world does.
    cdef int i_eccentricity_truncation = validate_eccentricity_truncation(eccentricity_truncation)
    cdef double eccentricity_tolerance = validate_eccentricity_exact_tolerance(eccentricity_exact_tolerance)

    # The same defaults and key check as make_tide: an absent config takes the [tides] defaults, and an unknown or
    # misspelled key raises instead of silently leaving a list empty (which would give no heating).
    if tide_config is not None:
        check_config_keys(tide_config, TIDE_CONFIG_KEYS, "tide")
    tide_config = {**factory_defaults("tides", TIDE_CONFIG_KEYS, tide_model, _same_model), **(tide_config or {})}
    cdef c_TideModelConfig cfg = cy_build_tide_config(tide_config)
    cdef c_TideModel model_enum = c_tide_model_from_name(tide_model.encode("utf-8"))
    if model_enum == c_TideModel.Rheology:
        raise NotImplementedError(
            "collapse_global_tides supports the analytic tide models only "
            "(cpl/fixed_q, ctl/fixed_dt, ctl_q/fixed_dt_q). The rheology model needs the "
            "radial solver; use the world's calc_tides method.")
    cdef unique_ptr[c_TideBase] tide_ptr = c_find_tide(model_enum, cfg)

    cdef c_GlobalPotentialStorage potential = c_global_potential(
        planet_radius,
        semi_major_axis,
        orbital_frequency,
        spin_frequency,
        obliquity,
        eccentricity,
        host_mass,
        G_to_use,
        min_degree_l,
        max_degree_l,
        i_obliquity_truncation,
        i_eccentricity_truncation,
        eccentricity_tolerance)

    if potential.error_code != 0:
        raise RuntimeError(
            f"Global potential failed with error code {potential.error_code} "
            f"(working on degree l={potential.working_on_l}).")

    cdef c_GlobalTideResult result = c_collapse_global_tides(potential, tide_ptr.get()[0])

    return {
        "tidal_heating": result.tidal_heating,
        "dUdM": result.dU_dM,
        "dUdw": result.dU_dw,
        "dUdO": result.dU_dO,
        "dUdM_minus_dw": result.dU_dM_minus_dw,
        "num_modes": result.num_modes,
    }
