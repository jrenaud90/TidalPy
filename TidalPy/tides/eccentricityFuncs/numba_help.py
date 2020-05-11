from . import EccenOutput
from .orderl2 import eccentricity_funcs_trunc2 as eccentricity_funcs_l2_trunc2, \
    eccentricity_funcs_trunc4 as eccentricity_funcs_l2_trunc4, \
    eccentricity_funcs_trunc6 as eccentricity_funcs_l2_trunc6, \
    eccentricity_funcs_trunc8 as eccentricity_funcs_l2_trunc8, \
    eccentricity_funcs_trunc10 as eccentricity_funcs_l2_trunc10, \
    eccentricity_funcs_trunc12 as eccentricity_funcs_l2_trunc12, \
    eccentricity_funcs_trunc14 as eccentricity_funcs_l2_trunc14, \
    eccentricity_funcs_trunc16 as eccentricity_funcs_l2_trunc16, \
    eccentricity_funcs_trunc18 as eccentricity_funcs_l2_trunc18, \
    eccentricity_funcs_trunc20 as eccentricity_funcs_l2_trunc20
from .orderl3 import eccentricity_funcs_trunc2 as eccentricity_funcs_l3_trunc2, \
    eccentricity_funcs_trunc4 as eccentricity_funcs_l3_trunc4, \
    eccentricity_funcs_trunc6 as eccentricity_funcs_l3_trunc6, \
    eccentricity_funcs_trunc8 as eccentricity_funcs_l3_trunc8, \
    eccentricity_funcs_trunc10 as eccentricity_funcs_l3_trunc10, \
    eccentricity_funcs_trunc12 as eccentricity_funcs_l3_trunc12, \
    eccentricity_funcs_trunc14 as eccentricity_funcs_l3_trunc14, \
    eccentricity_funcs_trunc16 as eccentricity_funcs_l3_trunc16, \
    eccentricity_funcs_trunc18 as eccentricity_funcs_l3_trunc18, \
    eccentricity_funcs_trunc20 as eccentricity_funcs_l3_trunc20
from .orderl4 import eccentricity_funcs_trunc2 as eccentricity_funcs_l4_trunc2, \
    eccentricity_funcs_trunc4 as eccentricity_funcs_l4_trunc4, \
    eccentricity_funcs_trunc6 as eccentricity_funcs_l4_trunc6, \
    eccentricity_funcs_trunc8 as eccentricity_funcs_l4_trunc8, \
    eccentricity_funcs_trunc10 as eccentricity_funcs_l4_trunc10, \
    eccentricity_funcs_trunc12 as eccentricity_funcs_l4_trunc12, \
    eccentricity_funcs_trunc14 as eccentricity_funcs_l4_trunc14, \
    eccentricity_funcs_trunc16 as eccentricity_funcs_l4_trunc16, \
    eccentricity_funcs_trunc18 as eccentricity_funcs_l4_trunc18, \
    eccentricity_funcs_trunc20 as eccentricity_funcs_l4_trunc20
from .orderl5 import eccentricity_funcs_trunc2 as eccentricity_funcs_l5_trunc2, \
    eccentricity_funcs_trunc4 as eccentricity_funcs_l5_trunc4, \
    eccentricity_funcs_trunc6 as eccentricity_funcs_l5_trunc6, \
    eccentricity_funcs_trunc8 as eccentricity_funcs_l5_trunc8, \
    eccentricity_funcs_trunc10 as eccentricity_funcs_l5_trunc10, \
    eccentricity_funcs_trunc12 as eccentricity_funcs_l5_trunc12, \
    eccentricity_funcs_trunc14 as eccentricity_funcs_l5_trunc14, \
    eccentricity_funcs_trunc16 as eccentricity_funcs_l5_trunc16, \
    eccentricity_funcs_trunc18 as eccentricity_funcs_l5_trunc18, \
    eccentricity_funcs_trunc20 as eccentricity_funcs_l5_trunc20
from .orderl6 import eccentricity_funcs_trunc2 as eccentricity_funcs_l6_trunc2, \
    eccentricity_funcs_trunc4 as eccentricity_funcs_l6_trunc4, \
    eccentricity_funcs_trunc6 as eccentricity_funcs_l6_trunc6, \
    eccentricity_funcs_trunc8 as eccentricity_funcs_l6_trunc8, \
    eccentricity_funcs_trunc10 as eccentricity_funcs_l6_trunc10, \
    eccentricity_funcs_trunc12 as eccentricity_funcs_l6_trunc12, \
    eccentricity_funcs_trunc14 as eccentricity_funcs_l6_trunc14, \
    eccentricity_funcs_trunc16 as eccentricity_funcs_l6_trunc16, \
    eccentricity_funcs_trunc18 as eccentricity_funcs_l6_trunc18, \
    eccentricity_funcs_trunc20 as eccentricity_funcs_l6_trunc20
from .orderl7 import eccentricity_funcs_trunc2 as eccentricity_funcs_l7_trunc2, \
    eccentricity_funcs_trunc4 as eccentricity_funcs_l7_trunc4, \
    eccentricity_funcs_trunc6 as eccentricity_funcs_l7_trunc6, \
    eccentricity_funcs_trunc8 as eccentricity_funcs_l7_trunc8, \
    eccentricity_funcs_trunc10 as eccentricity_funcs_l7_trunc10, \
    eccentricity_funcs_trunc12 as eccentricity_funcs_l7_trunc12, \
    eccentricity_funcs_trunc14 as eccentricity_funcs_l7_trunc14, \
    eccentricity_funcs_trunc16 as eccentricity_funcs_l7_trunc16, \
    eccentricity_funcs_trunc18 as eccentricity_funcs_l7_trunc18, \
    eccentricity_funcs_trunc20 as eccentricity_funcs_l7_trunc20
from ...utilities.performance import njit
from ...utilities.types import FloatArray
from ...configurations import cache_numba


@njit(cache=cache_numba)
def eccentricity_func_trunc_2(eccentricity: FloatArray, tidal_order: int = 2) -> EccenOutput:
    """ Calculates eccentricity functions (squared) for a given tidal order, at truncation level 2.

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    eccentricity_result : FloatArray
        Eccentricity function G^2_lpq truncated to e^2
    """

    if tidal_order == 2:
        result = eccentricity_funcs_l2_trunc2(eccentricity)
    elif tidal_order == 3:
        result = eccentricity_funcs_l3_trunc2(eccentricity)
    elif tidal_order == 4:
        result = eccentricity_funcs_l4_trunc2(eccentricity)
    elif tidal_order == 5:
        result = eccentricity_funcs_l5_trunc2(eccentricity)
    elif tidal_order == 6:
        result = eccentricity_funcs_l6_trunc2(eccentricity)
    elif tidal_order == 7:
        result = eccentricity_funcs_l7_trunc2(eccentricity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Tidal order not supported for level 2 truncation.')

    return result


@njit(cache=cache_numba)
def eccentricity_func_trunc_4(eccentricity: FloatArray, tidal_order: int = 2) -> EccenOutput:
    """ Calculates eccentricity functions (squared) for a given tidal order, at truncation level 4.

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    eccentricity_result : FloatArray
        Eccentricity function G^2_lpq truncated to e^4
    """

    if tidal_order == 2:
        result = eccentricity_funcs_l2_trunc4(eccentricity)
    elif tidal_order == 3:
        result = eccentricity_funcs_l3_trunc4(eccentricity)
    elif tidal_order == 4:
        result = eccentricity_funcs_l4_trunc4(eccentricity)
    elif tidal_order == 5:
        result = eccentricity_funcs_l5_trunc4(eccentricity)
    elif tidal_order == 6:
        result = eccentricity_funcs_l6_trunc4(eccentricity)
    elif tidal_order == 7:
        result = eccentricity_funcs_l7_trunc4(eccentricity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Tidal order not supported for level 4 truncation.')

    return result

@njit(cache=cache_numba)
def eccentricity_func_trunc_6(eccentricity: FloatArray, tidal_order: int = 2) -> EccenOutput:
    """ Calculates eccentricity functions (squared) for a given tidal order, at truncation level 6.

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    eccentricity_result : FloatArray
        Eccentricity function G^2_lpq truncated to e^6
    """

    if tidal_order == 2:
        result = eccentricity_funcs_l2_trunc6(eccentricity)
    elif tidal_order == 3:
        result = eccentricity_funcs_l3_trunc6(eccentricity)
    elif tidal_order == 4:
        result = eccentricity_funcs_l4_trunc6(eccentricity)
    elif tidal_order == 5:
        result = eccentricity_funcs_l5_trunc6(eccentricity)
    elif tidal_order == 6:
        result = eccentricity_funcs_l6_trunc6(eccentricity)
    elif tidal_order == 7:
        result = eccentricity_funcs_l7_trunc6(eccentricity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Tidal order not supported for level 6 truncation.')

    return result

@njit(cache=cache_numba)
def eccentricity_func_trunc_8(eccentricity: FloatArray, tidal_order: int = 2) -> EccenOutput:
    """ Calculates eccentricity functions (squared) for a given tidal order, at truncation level 8.

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    eccentricity_result : FloatArray
        Eccentricity function G^2_lpq truncated to e^8
    """

    if tidal_order == 2:
        result = eccentricity_funcs_l2_trunc8(eccentricity)
    elif tidal_order == 3:
        result = eccentricity_funcs_l3_trunc8(eccentricity)
    elif tidal_order == 4:
        result = eccentricity_funcs_l4_trunc8(eccentricity)
    elif tidal_order == 5:
        result = eccentricity_funcs_l5_trunc8(eccentricity)
    elif tidal_order == 6:
        result = eccentricity_funcs_l6_trunc8(eccentricity)
    elif tidal_order == 7:
        result = eccentricity_funcs_l7_trunc8(eccentricity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Tidal order not supported for level 8 truncation.')

    return result

@njit(cache=cache_numba)
def eccentricity_func_trunc_10(eccentricity: FloatArray, tidal_order: int = 2) -> EccenOutput:
    """ Calculates eccentricity functions (squared) for a given tidal order, at truncation level 10.

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    eccentricity_result : FloatArray
        Eccentricity function G^2_lpq truncated to e^10
    """

    if tidal_order == 2:
        result = eccentricity_funcs_l2_trunc10(eccentricity)
    elif tidal_order == 3:
        result = eccentricity_funcs_l3_trunc10(eccentricity)
    elif tidal_order == 4:
        result = eccentricity_funcs_l4_trunc10(eccentricity)
    elif tidal_order == 5:
        result = eccentricity_funcs_l5_trunc10(eccentricity)
    elif tidal_order == 6:
        result = eccentricity_funcs_l6_trunc10(eccentricity)
    elif tidal_order == 7:
        result = eccentricity_funcs_l7_trunc10(eccentricity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Tidal order not supported for level 10 truncation.')

    return result

@njit(cache=cache_numba)
def eccentricity_func_trunc_12(eccentricity: FloatArray, tidal_order: int = 2) -> EccenOutput:
    """ Calculates eccentricity functions (squared) for a given tidal order, at truncation level 12.

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    eccentricity_result : FloatArray
        Eccentricity function G^2_lpq truncated to e^12
    """

    if tidal_order == 2:
        result = eccentricity_funcs_l2_trunc12(eccentricity)
    elif tidal_order == 3:
        result = eccentricity_funcs_l3_trunc12(eccentricity)
    elif tidal_order == 4:
        result = eccentricity_funcs_l4_trunc12(eccentricity)
    elif tidal_order == 5:
        result = eccentricity_funcs_l5_trunc12(eccentricity)
    elif tidal_order == 6:
        result = eccentricity_funcs_l6_trunc12(eccentricity)
    elif tidal_order == 7:
        result = eccentricity_funcs_l7_trunc12(eccentricity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Tidal order not supported for level 12 truncation.')

    return result

@njit(cache=cache_numba)
def eccentricity_func_trunc_14(eccentricity: FloatArray, tidal_order: int = 2) -> EccenOutput:
    """ Calculates eccentricity functions (squared) for a given tidal order, at truncation level 14.

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    eccentricity_result : FloatArray
        Eccentricity function G^2_lpq truncated to e^12
    """

    if tidal_order == 2:
        result = eccentricity_funcs_l2_trunc14(eccentricity)
    elif tidal_order == 3:
        result = eccentricity_funcs_l3_trunc14(eccentricity)
    elif tidal_order == 4:
        result = eccentricity_funcs_l4_trunc14(eccentricity)
    elif tidal_order == 5:
        result = eccentricity_funcs_l5_trunc14(eccentricity)
    elif tidal_order == 6:
        result = eccentricity_funcs_l6_trunc14(eccentricity)
    elif tidal_order == 7:
        result = eccentricity_funcs_l7_trunc14(eccentricity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Tidal order not supported for level 14 truncation.')

    return result

@njit(cache=cache_numba)
def eccentricity_func_trunc_16(eccentricity: FloatArray, tidal_order: int = 2) -> EccenOutput:
    """ Calculates eccentricity functions (squared) for a given tidal order, at truncation level 16.

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    eccentricity_result : FloatArray
        Eccentricity function G^2_lpq truncated to e^16
    """

    if tidal_order == 2:
        result = eccentricity_funcs_l2_trunc16(eccentricity)
    elif tidal_order == 3:
        result = eccentricity_funcs_l3_trunc16(eccentricity)
    elif tidal_order == 4:
        result = eccentricity_funcs_l4_trunc16(eccentricity)
    elif tidal_order == 5:
        result = eccentricity_funcs_l5_trunc16(eccentricity)
    elif tidal_order == 6:
        result = eccentricity_funcs_l6_trunc16(eccentricity)
    elif tidal_order == 7:
        result = eccentricity_funcs_l7_trunc16(eccentricity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Tidal order not supported for level 16 truncation.')

    return result

@njit(cache=cache_numba)
def eccentricity_func_trunc_18(eccentricity: FloatArray, tidal_order: int = 2) -> EccenOutput:
    """ Calculates eccentricity functions (squared) for a given tidal order, at truncation level 18.

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    eccentricity_result : FloatArray
        Eccentricity function G^2_lpq truncated to e^18
    """

    if tidal_order == 2:
        result = eccentricity_funcs_l2_trunc18(eccentricity)
    elif tidal_order == 3:
        result = eccentricity_funcs_l3_trunc18(eccentricity)
    elif tidal_order == 4:
        result = eccentricity_funcs_l4_trunc18(eccentricity)
    elif tidal_order == 5:
        result = eccentricity_funcs_l5_trunc18(eccentricity)
    elif tidal_order == 6:
        result = eccentricity_funcs_l6_trunc18(eccentricity)
    elif tidal_order == 7:
        result = eccentricity_funcs_l7_trunc18(eccentricity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Tidal order not supported for level 18 truncation.')

    return result

@njit(cache=cache_numba)
def eccentricity_func_trunc_20(eccentricity: FloatArray, tidal_order: int = 2) -> EccenOutput:
    """ Calculates eccentricity functions (squared) for a given tidal order, at truncation level 20.

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    eccentricity_result : FloatArray
        Eccentricity function G^2_lpq truncated to e^20
    """

    if tidal_order == 2:
        result = eccentricity_funcs_l2_trunc20(eccentricity)
    elif tidal_order == 3:
        result = eccentricity_funcs_l3_trunc20(eccentricity)
    elif tidal_order == 4:
        result = eccentricity_funcs_l4_trunc20(eccentricity)
    elif tidal_order == 5:
        result = eccentricity_funcs_l5_trunc20(eccentricity)
    elif tidal_order == 6:
        result = eccentricity_funcs_l6_trunc20(eccentricity)
    elif tidal_order == 7:
        result = eccentricity_funcs_l7_trunc20(eccentricity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Tidal order not supported for level 20 truncation.')

    return result

eccentricity_func_by_truncation_lvl = {
    2: eccentricity_func_trunc_2,
    4: eccentricity_func_trunc_4,
    6: eccentricity_func_trunc_6,
    8: eccentricity_func_trunc_8,
    10: eccentricity_func_trunc_10,
    12: eccentricity_func_trunc_12,
    14: eccentricity_func_trunc_14,
    16: eccentricity_func_trunc_16,
    18: eccentricity_func_trunc_18,
    20: eccentricity_func_trunc_20,
}

# The _fast versions below are restricted to only order 2 and 3. This reduces JIT compile times. The above general
#    functions must compile all orders (2 -- max) in order to tun. In contrast these fast functions only need to compile
#    l=2 and l=3.
MAX_L_FOR_FAST_SUPPORT = 3

@njit(cache=cache_numba)
def eccentricity_func_trunc_2_fast(eccentricity: FloatArray, tidal_order: int = 2) -> EccenOutput:
    """ Calculates eccentricity functions (squared) for a given tidal order (max=3), at truncation level 2.

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    eccentricity_result : FloatArray
        Eccentricity function G^2_lpq truncated to e^2
    """

    if tidal_order == 2:
        result = eccentricity_funcs_l2_trunc2(eccentricity)
    elif tidal_order == 3:
        result = eccentricity_funcs_l3_trunc2(eccentricity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Fast eccentricity function only supports l=2 and l=3. Use general form instead.')

    return result


@njit(cache=cache_numba)
def eccentricity_func_trunc_4_fast(eccentricity: FloatArray, tidal_order: int = 2) -> EccenOutput:
    """ Calculates eccentricity functions (squared) for a given tidal order (max=3), at truncation level 4.

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    eccentricity_result : FloatArray
        Eccentricity function G^2_lpq truncated to e^4
    """

    if tidal_order == 2:
        result = eccentricity_funcs_l2_trunc4(eccentricity)
    elif tidal_order == 3:
        result = eccentricity_funcs_l3_trunc4(eccentricity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Fast eccentricity function only supports l=2 and l=3. Use general form instead.')

    return result

@njit(cache=cache_numba)
def eccentricity_func_trunc_6_fast(eccentricity: FloatArray, tidal_order: int = 2) -> EccenOutput:
    """ Calculates eccentricity functions (squared) for a given tidal order (max=3), at truncation level 6.

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    eccentricity_result : FloatArray
        Eccentricity function G^2_lpq truncated to e^6
    """

    if tidal_order == 2:
        result = eccentricity_funcs_l2_trunc6(eccentricity)
    elif tidal_order == 3:
        result = eccentricity_funcs_l3_trunc6(eccentricity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Fast eccentricity function only supports l=2 and l=3. Use general form instead.')

    return result

@njit(cache=cache_numba)
def eccentricity_func_trunc_8_fast(eccentricity: FloatArray, tidal_order: int = 2) -> EccenOutput:
    """ Calculates eccentricity functions (squared) for a given tidal order (max=3), at truncation level 8.

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    eccentricity_result : FloatArray
        Eccentricity function G^2_lpq truncated to e^8
    """

    if tidal_order == 2:
        result = eccentricity_funcs_l2_trunc8(eccentricity)
    elif tidal_order == 3:
        result = eccentricity_funcs_l3_trunc8(eccentricity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Fast eccentricity function only supports l=2 and l=3. Use general form instead.')

    return result

@njit(cache=cache_numba)
def eccentricity_func_trunc_10_fast(eccentricity: FloatArray, tidal_order: int = 2) -> EccenOutput:
    """ Calculates eccentricity functions (squared) for a given tidal order (max=3), at truncation level 10.

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    eccentricity_result : FloatArray
        Eccentricity function G^2_lpq truncated to e^10
    """

    if tidal_order == 2:
        result = eccentricity_funcs_l2_trunc10(eccentricity)
    elif tidal_order == 3:
        result = eccentricity_funcs_l3_trunc10(eccentricity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Fast eccentricity function only supports l=2 and l=3. Use general form instead.')

    return result

@njit(cache=cache_numba)
def eccentricity_func_trunc_12_fast(eccentricity: FloatArray, tidal_order: int = 2) -> EccenOutput:
    """ Calculates eccentricity functions (squared) for a given tidal order (max=3), at truncation level 12.

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    eccentricity_result : FloatArray
        Eccentricity function G^2_lpq truncated to e^12
    """

    if tidal_order == 2:
        result = eccentricity_funcs_l2_trunc12(eccentricity)
    elif tidal_order == 3:
        result = eccentricity_funcs_l3_trunc12(eccentricity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Fast eccentricity function only supports l=2 and l=3. Use general form instead.')

    return result

@njit(cache=cache_numba)
def eccentricity_func_trunc_14_fast(eccentricity: FloatArray, tidal_order: int = 2) -> EccenOutput:
    """ Calculates eccentricity functions (squared) for a given tidal order (max=3), at truncation level 14.

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    eccentricity_result : FloatArray
        Eccentricity function G^2_lpq truncated to e^12
    """

    if tidal_order == 2:
        result = eccentricity_funcs_l2_trunc14(eccentricity)
    elif tidal_order == 3:
        result = eccentricity_funcs_l3_trunc14(eccentricity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Fast eccentricity function only supports l=2 and l=3. Use general form instead.')

    return result

@njit(cache=cache_numba)
def eccentricity_func_trunc_16_fast(eccentricity: FloatArray, tidal_order: int = 2) -> EccenOutput:
    """ Calculates eccentricity functions (squared) for a given tidal order (max=3), at truncation level 16.

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    eccentricity_result : FloatArray
        Eccentricity function G^2_lpq truncated to e^16
    """

    if tidal_order == 2:
        result = eccentricity_funcs_l2_trunc16(eccentricity)
    elif tidal_order == 3:
        result = eccentricity_funcs_l3_trunc16(eccentricity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Fast eccentricity function only supports l=2 and l=3. Use general form instead.')

    return result

@njit(cache=cache_numba)
def eccentricity_func_trunc_18_fast(eccentricity: FloatArray, tidal_order: int = 2) -> EccenOutput:
    """ Calculates eccentricity functions (squared) for a given tidal order (max=3), at truncation level 18.

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    eccentricity_result : FloatArray
        Eccentricity function G^2_lpq truncated to e^18
    """

    if tidal_order == 2:
        result = eccentricity_funcs_l2_trunc18(eccentricity)
    elif tidal_order == 3:
        result = eccentricity_funcs_l3_trunc18(eccentricity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Fast eccentricity function only supports l=2 and l=3. Use general form instead.')

    return result

@njit(cache=cache_numba)
def eccentricity_func_trunc_20_fast(eccentricity: FloatArray, tidal_order: int = 2) -> EccenOutput:
    """ Calculates eccentricity functions (squared) for a given tidal order (max=3), at truncation level 20.

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    eccentricity_result : FloatArray
        Eccentricity function G^2_lpq truncated to e^20
    """

    if tidal_order == 2:
        result = eccentricity_funcs_l2_trunc20(eccentricity)
    elif tidal_order == 3:
        result = eccentricity_funcs_l3_trunc20(eccentricity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Fast eccentricity function only supports l=2 and l=3. Use general form instead.')

    return result

eccentricity_func_by_truncation_lvl_fast = {
    2: eccentricity_func_trunc_2_fast,
    4: eccentricity_func_trunc_4_fast,
    6: eccentricity_func_trunc_6_fast,
    8: eccentricity_func_trunc_8_fast,
    10: eccentricity_func_trunc_10_fast,
    12: eccentricity_func_trunc_12_fast,
    14: eccentricity_func_trunc_14_fast,
    16: eccentricity_func_trunc_16_fast,
    18: eccentricity_func_trunc_18_fast,
    20: eccentricity_func_trunc_20_fast,
}
