from ...utilities.performance import njit
from ...utilities.types import FloatArray


@njit(cacheable=True)
def linear_dt(frequency: FloatArray, fixed_dt: float):
    """ Estimates dissipative term of the Love number assuming a function that is inversely linear with frequency.

    Parameters
    ----------
    frequency : FloatArray
        Frequency (the absolute value of the tidal modes)
    fixed_dt : float
        Inverse proportionality coefficient [s]

    Returns
    -------
    effective_q : FloatArray
        The effective Q for the world
    """

    effective_q = frequency * fixed_dt

    return effective_q

@njit(cacheable=True)
def linear_dt_with_q(frequency: FloatArray, fixed_dt: float, fixed_q: float):
    """ Estimates dissipative term of the Love number assuming a function that is linear with frequency. The fixed Q
        acts as an additional inverse proportionality constant.

    Parameters
    ----------
    frequency : FloatArray
        Frequency (the absolute value of the tidal modes)
    fixed_dt : float
        Inverse proportionality coefficient [s]
    fixed_q : float
        Additional inverse proportionality coefficient

    Returns
    -------
    effective_q : FloatArray
        The effective Q for the world
    """

    effective_q = frequency * fixed_dt / fixed_q

    return effective_q