from typing import Dict

from ..eccentricity_funcs import orderl2, EccenOutput
from ...utilities.performance.numba import njit
from ...utilities.types import FloatArray


@njit(cacheable=True)
def eccentricity_truncation_2_maxl_2(eccentricity: FloatArray) -> Dict[int, EccenOutput]:
    """ Calculates eccentricity functions (squared) for a given maximum tidal order (going through each l)

    Truncation level = 2
    Max Supported l = 2

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    result_by_orderl : Dict[int, EccenOutput]
        Eccentricity function G^2_lpq truncated. Stored by order-l.
    """

    result_by_orderl = {
        2: orderl2.eccentricity_funcs_trunc2(eccentricity)
    }

    return result_by_orderl

@njit(cacheable=True)
def eccentricity_truncation_4_maxl_2(eccentricity: FloatArray) -> Dict[int, EccenOutput]:
    """ Calculates eccentricity functions (squared) for a given maximum tidal order (going through each l)

    Truncation level = 4
    Max Supported l = 2

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    result_by_orderl : Dict[int, EccenOutput]
        Eccentricity function G^2_lpq truncated. Stored by order-l.
    """

    result_by_orderl = {
        2: orderl2.eccentricity_funcs_trunc4(eccentricity)
    }

    return result_by_orderl

@njit(cacheable=True)
def eccentricity_truncation_6_maxl_2(eccentricity: FloatArray) -> Dict[int, EccenOutput]:
    """ Calculates eccentricity functions (squared) for a given maximum tidal order (going through each l)

    Truncation level = 6
    Max Supported l = 2

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    result_by_orderl : Dict[int, EccenOutput]
        Eccentricity function G^2_lpq truncated. Stored by order-l.
    """

    result_by_orderl = {
        2: orderl2.eccentricity_funcs_trunc6(eccentricity)
    }

    return result_by_orderl

@njit(cacheable=True)
def eccentricity_truncation_8_maxl_2(eccentricity: FloatArray) -> Dict[int, EccenOutput]:
    """ Calculates eccentricity functions (squared) for a given maximum tidal order (going through each l)

    Truncation level = 8
    Max Supported l = 2

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    result_by_orderl : Dict[int, EccenOutput]
        Eccentricity function G^2_lpq truncated. Stored by order-l.
    """

    result_by_orderl = {
        2: orderl2.eccentricity_funcs_trunc8(eccentricity)
    }

    return result_by_orderl

@njit(cacheable=True)
def eccentricity_truncation_10_maxl_2(eccentricity: FloatArray) -> Dict[int, EccenOutput]:
    """ Calculates eccentricity functions (squared) for a given maximum tidal order (going through each l)

    Truncation level = 10
    Max Supported l = 2

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    result_by_orderl : Dict[int, EccenOutput]
        Eccentricity function G^2_lpq truncated. Stored by order-l.
    """

    result_by_orderl = {
        2: orderl2.eccentricity_funcs_trunc10(eccentricity)
    }

    return result_by_orderl

@njit(cacheable=True)
def eccentricity_truncation_12_maxl_2(eccentricity: FloatArray) -> Dict[int, EccenOutput]:
    """ Calculates eccentricity functions (squared) for a given maximum tidal order (going through each l)

    Truncation level = 12
    Max Supported l = 2

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    result_by_orderl : Dict[int, EccenOutput]
        Eccentricity function G^2_lpq truncated. Stored by order-l.
    """

    result_by_orderl = {
        2: orderl2.eccentricity_funcs_trunc12(eccentricity)
    }

    return result_by_orderl

@njit(cacheable=True)
def eccentricity_truncation_14_maxl_2(eccentricity: FloatArray) -> Dict[int, EccenOutput]:
    """ Calculates eccentricity functions (squared) for a given maximum tidal order (going through each l)

    Truncation level = 14
    Max Supported l = 2

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    result_by_orderl : Dict[int, EccenOutput]
        Eccentricity function G^2_lpq truncated. Stored by order-l.
    """

    result_by_orderl = {
        2: orderl2.eccentricity_funcs_trunc14(eccentricity)
    }

    return result_by_orderl

@njit(cacheable=True)
def eccentricity_truncation_16_maxl_2(eccentricity: FloatArray) -> Dict[int, EccenOutput]:
    """ Calculates eccentricity functions (squared) for a given maximum tidal order (going through each l)

    Truncation level = 16
    Max Supported l = 2

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    result_by_orderl : Dict[int, EccenOutput]
        Eccentricity function G^2_lpq truncated. Stored by order-l.
    """

    result_by_orderl = {
        2: orderl2.eccentricity_funcs_trunc16(eccentricity)
    }

    return result_by_orderl

@njit(cacheable=True)
def eccentricity_truncation_18_maxl_2(eccentricity: FloatArray) -> Dict[int, EccenOutput]:
    """ Calculates eccentricity functions (squared) for a given maximum tidal order (going through each l)

    Truncation level = 18
    Max Supported l = 2

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    result_by_orderl : Dict[int, EccenOutput]
        Eccentricity function G^2_lpq truncated. Stored by order-l.
    """

    result_by_orderl = {
        2: orderl2.eccentricity_funcs_trunc18(eccentricity)
    }

    return result_by_orderl

@njit(cacheable=True)
def eccentricity_truncation_20_maxl_2(eccentricity: FloatArray) -> Dict[int, EccenOutput]:
    """ Calculates eccentricity functions (squared) for a given maximum tidal order (going through each l)

    Truncation level = 20
    Max Supported l = 2

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    result_by_orderl : Dict[int, EccenOutput]
        Eccentricity function G^2_lpq truncated. Stored by order-l.
    """

    result_by_orderl = {
        2: orderl2.eccentricity_funcs_trunc20(eccentricity)
    }

    return result_by_orderl

@njit(cacheable=True)
def eccentricity_truncation_22_maxl_2(eccentricity: FloatArray) -> Dict[int, EccenOutput]:
    """ Calculates eccentricity functions (squared) for a given maximum tidal order (going through each l)

    Truncation level = 22
    Max Supported l = 2

    Parameters
    ----------
    eccentricity : FloatArray
        Orbital Eccentricity
    Returns
    -------
    result_by_orderl : Dict[int, EccenOutput]
        Eccentricity function G^2_lpq truncated. Stored by order-l.
    """

    result_by_orderl = {
        2: orderl2.eccentricity_funcs_trunc22(eccentricity)
    }

    return result_by_orderl