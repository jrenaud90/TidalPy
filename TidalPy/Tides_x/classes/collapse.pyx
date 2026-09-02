# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
collapse.pyx
Standalone global (1D) tidal-mode collapse.

``collapse_global_tides`` runs the global-potential engine (the eccentricity/obliquity
functions + tidal potential of Renaud et al. 2021) for a given orbital/spin state, then
collapses the per-mode potential terms with an analytic tide model's dissipation
multiplier -Im[k_l] to give the world's global tidal heating and the three orbital
potential derivatives (dU/dM, dU/dw, dU/dO).

This standalone path supports the analytic models only ("cpl", "ctl", "ctl_q").
The rheology model needs per-mode Love numbers from the radial
solver and is driven by the world's ``calc_tides`` method instead.
"""

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr

from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.Tides_x.classes.tide cimport (
    c_TideBase, c_TideModel, c_TideModelConfig, c_tide_model_from_name, c_find_tide,
)

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


# =====================================================================================================================
# C++ declarations
# =====================================================================================================================
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
        int eccentricity_truncation) except +


cdef extern from "tide_collapse_.hpp" nogil:

    cdef cppclass c_GlobalTideResult:
        double tidal_heating
        double dU_dM
        double dU_dw
        double dU_dO
        int num_modes
        int error_code

    c_GlobalTideResult c_collapse_global_tides(
        const c_GlobalPotentialStorage& potential,
        const c_TideBase& tide_model) except +


# =====================================================================================================================
# Helpers
# =====================================================================================================================
cdef int _resolve_obliquity_truncation(object obliquity_truncation) except? -999:
    """Normalize the obliquity truncation (string 'gen'/'off' or int) to the C++ integer."""
    cdef int value = 0
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
    elif isinstance(obliquity_truncation, int):
        value = obliquity_truncation
    if value not in (0, 2, 4, 10):
        raise NotImplementedError(
            f"Obliquity truncation {value} is not tabulated. "
            "Supported levels: 0 ('off'), 2, 4, 10 ('gen', fully general).")
    return value


cdef c_TideModelConfig _build_tide_config(dict config) except *:
    """Build a c_TideModelConfig from optional per-degree list keys (indexed from l=2)."""
    cdef c_TideModelConfig cfg
    if config is None:
        return cfg
    if "fixed_k" in config:
        for value in config["fixed_k"]:
            cfg.fixed_k.push_back(<double>value)
    if "fixed_q" in config:
        for value in config["fixed_q"]:
            cfg.fixed_q.push_back(<double>value)
    if "fixed_dt" in config:
        for value in config["fixed_dt"]:
            cfg.fixed_dt.push_back(<double>value)
    return cfg


# =====================================================================================================================
# Public API
# =====================================================================================================================
def collapse_global_tides(
        double planet_radius,
        double semi_major_axis,
        double orbital_frequency,
        double spin_frequency,
        double obliquity,
        double eccentricity,
        double host_mass,
        double G_to_use,
        str tide_model,
        dict tide_config=None,
        int min_degree_l=2,
        int max_degree_l=2,
        object obliquity_truncation='gen',
        int eccentricity_truncation=3) -> dict:
    """Collapse the global tidal modes into heating and orbital potential derivatives.

    Parameters
    ----------
    planet_radius : float
        Radius of the tidally deformed body [m].
    semi_major_axis : float
        Orbital semi-major axis [m].
    orbital_frequency : float
        Orbital mean motion [rad s-1].
    spin_frequency : float
        Spin rate of the deformed body [rad s-1].
    obliquity : float
        Axial tilt [radians].
    eccentricity : float
        Orbital eccentricity [dimensionless].
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
        Per-degree model parameters (``fixed_k``, ``fixed_q``, ``fixed_dt`` lists indexed
        from degree l = 2).
    min_degree_l, max_degree_l : int
        Tidal harmonic degree range (2..10).
    obliquity_truncation : str or int
        Obliquity truncation: ``"off"`` (0), 2, 4, or ``"gen"``/``"general"`` (10).
    eccentricity_truncation : int
        Eccentricity truncation level. Tabulated levels: 1..5, 10, 15, 20.

    Returns
    -------
    dict
        ``{"tidal_heating": W, "dUdM": ..., "dUdw": ..., "dUdO": ..., "num_modes": int}``.

    Raises
    ------
    ValueError
        If the tide model name is unknown.
    NotImplementedError
        If the rheology model is requested, or a truncation/degree is unsupported.
    """
    cdef int i_obliquity_truncation = _resolve_obliquity_truncation(obliquity_truncation)

    if eccentricity_truncation not in (1, 2, 3, 4, 5, 10, 15, 20):
        raise NotImplementedError(
            f'Eccentricity truncation {eccentricity_truncation} is not tabulated. '
            'Supported levels: 1, 2, 3, 4, 5, 10, 15, 20.')

    # Build the tide model (analytic only).
    cdef c_TideModelConfig cfg = _build_tide_config(tide_config)
    cdef c_TideModel model_enum = c_tide_model_from_name(tide_model.encode("utf-8"))
    if model_enum == c_TideModel.Rheology:
        raise NotImplementedError(
            "collapse_global_tides supports the analytic tide models only "
            "(cpl/fixed_q, ctl/fixed_dt, ctl_q/fixed_dt_q). The rheology model needs the "
            "radial solver; use the world's calc_tides method.")
    cdef unique_ptr[c_TideBase] tide_ptr = c_find_tide(model_enum, cfg)

    # Run the global-potential engine.
    cdef c_GlobalPotentialStorage potential = c_global_potential(
        planet_radius, semi_major_axis, orbital_frequency, spin_frequency,
        obliquity, eccentricity, host_mass, G_to_use,
        min_degree_l, max_degree_l, i_obliquity_truncation, eccentricity_truncation)

    if potential.error_code != 0:
        raise RuntimeError(
            f"Global potential failed with error code {potential.error_code} "
            f"(working on degree l={potential.working_on_l}).")

    # Collapse with the tide model's dissipation multiplier.
    cdef c_GlobalTideResult result = c_collapse_global_tides(potential, tide_ptr.get()[0])

    return {
        "tidal_heating": result.tidal_heating,
        "dUdM": result.dU_dM,
        "dUdw": result.dU_dw,
        "dUdO": result.dU_dO,
        "num_modes": result.num_modes,
    }
