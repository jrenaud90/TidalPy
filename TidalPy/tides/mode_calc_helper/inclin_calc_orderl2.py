from typing import Dict

from ..inclination_funcs import orderl2, InclinOutput
from ...utilities.performance.numba import njit
from ...utilities.types import FloatArray


@njit(cacheable=True)
def inclination_off_maxl_2(obliquity: FloatArray) -> Dict[int, InclinOutput]:
    """ Calculates inclination functions (squared) for a given maximum tidal order (going through each l) - Off Mode

    Obliquity is assumed to be zero.

    Max Supported l = 2

    Parameters
    ----------
    obliquity : FloatArray
        Planet's obliquity [radians]
    Returns
    -------
    result_by_orderl : Dict[int, InclinOutput]
        Inclination function L^2_lmp truncated. Stored by order-l.
    """

    result_by_orderl = {
        2: orderl2.calc_inclination_off(obliquity)
    }

    return result_by_orderl


@njit(cacheable=True)
def inclination_on_maxl_2(obliquity: FloatArray) -> Dict[int, InclinOutput]:
    """ Calculates inclination functions (squared) for a given maximum tidal order (going through each l) - On Mode

    Obliquity can be arbitrary.

    Max Supported l = 2

    Parameters
    ----------
    obliquity : FloatArray
        Planet's obliquity [radians]
    Returns
    -------
    result_by_orderl : Dict[int, InclinOutput]
        Inclination function L^2_lmp truncated. Stored by order-l.
    """

    result_by_orderl = {
        2: orderl2.calc_inclination(obliquity)
    }

    return result_by_orderl