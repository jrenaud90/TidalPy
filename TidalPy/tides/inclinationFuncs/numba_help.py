from . import InclinOutput
from .orderl2 import calc_inclination as calc_inclin_l2, calc_inclination_off as calc_inclin_l2_off
from .orderl3 import calc_inclination as calc_inclin_l3, calc_inclination_off as calc_inclin_l3_off
from .orderl4 import calc_inclination as calc_inclin_l4, calc_inclination_off as calc_inclin_l4_off
from .orderl5 import calc_inclination as calc_inclin_l5, calc_inclination_off as calc_inclin_l5_off
from .orderl6 import calc_inclination as calc_inclin_l6, calc_inclination_off as calc_inclin_l6_off
from .orderl7 import calc_inclination as calc_inclin_l7, calc_inclination_off as calc_inclin_l7_off
from ...utilities.performance import njit
from ...utilities.types import FloatArray
from ...configurations import cache_numba


@njit(cache=cache_numba)
def inclination_func_off(obliquity: FloatArray, tidal_order: int = 2) -> InclinOutput:
    """ Calculate F^2_lmp (assuming I=0) for provided order l.

    Parameters
    ----------
    obliquity : FloatArray
        Planet's obliquity (axial tilt) relative to the orbital plane [rads]
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    result : InclinOutput
        F^2_lmp for provided eccentricity and tidal order l.
    """

    if tidal_order == 2:
        result = calc_inclin_l2_off(obliquity)
    elif tidal_order == 3:
        result = calc_inclin_l3_off(obliquity)
    elif tidal_order == 4:
        result = calc_inclin_l4_off(obliquity)
    elif tidal_order == 5:
        result = calc_inclin_l5_off(obliquity)
    elif tidal_order == 6:
        result = calc_inclin_l6_off(obliquity)
    elif tidal_order == 7:
        result = calc_inclin_l7_off(obliquity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Tidal order is not yet supported.')

    return result

@njit(cache=cache_numba)
def inclination_func_on(obliquity: FloatArray, tidal_order: int = 2) -> InclinOutput:
    """ Calculate F^2_lmp for provided order l.

    Parameters
    ----------
    obliquity : FloatArray
        Planet's obliquity (axial tilt) relative to the orbital plane [rads]
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    result : InclinOutput
        F^2_lmp for provided eccentricity and tidal order l.
    """

    if tidal_order == 2:
        result = calc_inclin_l2(obliquity)
    elif tidal_order == 3:
        result = calc_inclin_l3(obliquity)
    elif tidal_order == 4:
        result = calc_inclin_l4(obliquity)
    elif tidal_order == 5:
        result = calc_inclin_l5(obliquity)
    elif tidal_order == 6:
        result = calc_inclin_l6(obliquity)
    elif tidal_order == 7:
        result = calc_inclin_l7(obliquity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Tidal order is not yet supported.')

    return result

inclination_funcs = {
    True: inclination_func_on,
    False: inclination_func_off
}

# The _fast versions below are restricted to only order 2 and 3. This reduces JIT compile times. The above general
#    functions must compile all orders (2 -- max) in order to tun. In contrast these fast functions only need to compile
#    l=2 and l=3.

MAX_L_FOR_FAST_SUPPORT = 3

@njit(cache=cache_numba)
def inclination_func_off_fast(obliquity: FloatArray, tidal_order: int = 2) -> InclinOutput:
    """ Calculate F^2_lmp (assuming I=0) for provided order l. Max l supported = 3

    Parameters
    ----------
    obliquity : FloatArray
        Planet's obliquity (axial tilt) relative to the orbital plane [rads]
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    result : InclinOutput
        F^2_lmp for provided eccentricity and tidal order l.
    """

    if tidal_order == 2:
        result = calc_inclin_l2_off(obliquity)
    elif tidal_order == 3:
        result = calc_inclin_l3_off(obliquity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Fast inclination function only supports l=2 and l=3. Use general form instead.')

    return result

@njit(cache=cache_numba)
def inclination_func_on_fast(obliquity: FloatArray, tidal_order: int = 2) -> InclinOutput:
    """ Calculate F^2_lmp for provided order l. Max l supported = 3

    Parameters
    ----------
    obliquity : FloatArray
        Planet's obliquity (axial tilt) relative to the orbital plane [rads]
    tidal_order : int = 2
        Tidal harmonic order

    Returns
    -------
    result : InclinOutput
        F^2_lmp for provided eccentricity and tidal order l.
    """

    if tidal_order == 2:
        result = calc_inclin_l2(obliquity)
    elif tidal_order == 3:
        result = calc_inclin_l3(obliquity)
    else:
        # Numba does not support exceptions other than "Exception", can't use TidalPy exceptions
        raise Exception('Fast inclination function only supports l=2 and l=3. Use general form instead.')

    return result


inclination_funcs_fast = {
    True:  inclination_func_on_fast,
    False: inclination_func_off_fast
}