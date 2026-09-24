# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from TidalPy.constants cimport tidalpy_config_ptr, get_shared_config_address, set_tidalpy_config_ptr
set_tidalpy_config_ptr(get_shared_config_address())
from TidalPy.Tides_x.potential.potential_common cimport ModeMap, UniqueFrequencyMap
from TidalPy.Tides_x.potential.potential_common import ModeMap, UniqueFrequencyMap

from TidalPy.Tides_x.eccentricity.eccentricity_driver import (
    validate_eccentricity_exact_tolerance, validate_eccentricity_truncation)
from TidalPy.Tides_x.obliquity.obliquity_driver import validate_obliquity_truncation

def global_potential(
        double planet_radius,
        double orbital_frequency,
        double spin_frequency,
        double eccentricity,
        double obliquity,
        double semi_major_axis,
        double host_mass,
        double G_to_use,
        int min_degree_l=2,
        int max_degree_l=2,
        object eccentricity_truncation=None,
        object obliquity_truncation=None,
        object eccentricity_exact_tolerance=None
    ):
    """Build the global (1D) tidal potential mode tables for one orbital state.

    Arguments follow the world's ``calc_tides`` order after the body radius [m]; frequencies in rad s-1,
    angles in radians, MKS throughout. ``min_degree_l`` and ``max_degree_l`` must satisfy
    2 <= min <= max <= 10. ``obliquity_truncation`` is ``'off'`` (0), 2 or 4 (every product of two obliquity
    functions through I^N), or ``'gen'`` (the general functions); None takes the ``[tides]`` ``obliquity_trunc_lvl``
    of the TidalPy configuration (``'off'`` by default, which ignores the obliquity). ``eccentricity_truncation``
    is a level of ``ECCENTRICITY_TRUNCATIONS`` (every product of two eccentricity functions through e^N) or
    ``"exact"``; None takes the ``[tides]`` ``eccentricity_trunc_lvl``.
    ``eccentricity_exact_tolerance`` sets the mode range of ``"exact"``; None takes the ``[tides]`` value.

    Returns
    -------
    tuple
        ``(mode_map, unique_freq_index_map, unique_freq_list, potential_dict)``; ``potential_dict`` maps
        ``(l, m, p, q)`` to ``(dU_dM, dU_dw, dU_dO, E_dot)``.
    """
    if not (2 <= min_degree_l <= max_degree_l <= 10):
        raise ValueError(
            f"The degree range must satisfy 2 <= min_degree_l <= max_degree_l <= 10; got {min_degree_l} to "
            f"{max_degree_l}.")
    # None takes the [tides] obliquity_trunc_lvl of the TidalPy configuration, as a built world does.
    cdef int i_obliquity_truncation = validate_obliquity_truncation(obliquity_truncation)
    # None takes the [tides] eccentricity_trunc_lvl of the TidalPy configuration, as a built world does.
    cdef int i_eccentricity_truncation = validate_eccentricity_truncation(eccentricity_truncation)
    cdef double eccentricity_tolerance = validate_eccentricity_exact_tolerance(eccentricity_exact_tolerance)

    cdef c_GlobalPotentialStorage c_result = c_global_potential(
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
        eccentricity_tolerance
    )

    if c_result.error_code != 0:
        if c_result.error_code == -20:
            raise NotImplementedError(
                f"Global potential error code -20: Could not find l,m coefficient "
                f"(working on degree l={c_result.working_on_l}). Perhaps unsupported degree l.")
        elif c_result.error_code == -1:
            raise NotImplementedError(
                f"Global potential error code -1: Unsupported obliquity truncation "
                f"(working on degree l={c_result.working_on_l}).")
        elif c_result.error_code == -2:
            raise NotImplementedError(
                f"Global potential error code -2: Unsupported degree l={c_result.working_on_l} "
                f"(the eccentricity and obliquity functions are tabulated for l = 2 to 10).")
        else:
            raise RuntimeError(
                f"Unknown global potential error code: {c_result.error_code} "
                f"(working on degree l={c_result.working_on_l}).")

    # Convert C++ results to Python-accessible objects.
    cdef ModeMap mode_map = ModeMap()
    mode_map._cinst = c_result.mode_map

    cdef UniqueFrequencyMap unique_freq_index_map = UniqueFrequencyMap()
    unique_freq_index_map._cinst = c_result.unique_freq_index_map

    # Unique frequency map (vector of c_FrequencyStorage -> list of (frequency, num_instances))
    cdef list unique_freq_list = []
    cdef size_t i
    for i in range(c_result.unique_freq_map.size()):
        unique_freq_list.append(
            (c_result.unique_freq_map[i].frequency, c_result.unique_freq_map[i].num_instances)
        )

    # Potential map (c_IntMap<c_Key4, c_GlobalPotentialResultAtMode> -> dict)
    cdef dict potential_dict = dict()
    cdef c_Key4 pkey
    # One index per mode; the four derivatives below came from four separate lookups.
    # A pointer, not a copy: the result struct has no default constructor.
    cdef const c_GlobalPotentialResultAtMode* mode_result_ptr = NULL
    for i in range(c_result.potential_map.size()):
        pkey = c_result.potential_map.data[i].first
        mode_result_ptr = &c_result.potential_map.data[i].second
        potential_dict[(pkey.a, pkey.b, pkey.c, pkey.d)] = (
            mode_result_ptr.dU_dM,
            mode_result_ptr.dU_dw,
            mode_result_ptr.dU_dO,
            mode_result_ptr.E_dot
        )

    return mode_map, unique_freq_index_map, unique_freq_list, potential_dict
