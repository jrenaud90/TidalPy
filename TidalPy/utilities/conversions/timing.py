from typing import Tuple

from TidalPy.utilities.performance.numba import njit


@njit(cacheable=True)
def convert_to_hms(seconds: float) -> Tuple[int, int, int, float]:
    """ Convert seconds to a tuple of days, hours, minutes, seconds

    Parameters
    ----------
    seconds : float
        Time in seconds

    Returns
    -------
    days : int
        Days
    hours : int
        Hours
    minutes : int
        Minutes
    seconds : float
        Remaining seconds after conversion
    """

    days = int(seconds / (24. * 3600.))
    seconds = seconds % (24. * 3600.)
    hours = int(seconds / 3600.)
    seconds = seconds % 3600.
    minutes = int(seconds / 60.)
    seconds = seconds % 60.

    return days, hours, minutes, seconds
