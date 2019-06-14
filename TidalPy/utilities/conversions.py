from ..performance import njit

@njit
def m2Au(meters):
    """ Convert Meters to Astronomical Units

    :param meters: <FloatArray> distance in meters
    :return:       <FloatArray> distance in Au
    """

    return meters / 1.496e11

@njit
def Au2m(astronomical_units):
    """ Convert Astronomical Units to Meters

    :param astronomical_units: <FloatArray> distance in meters
    :return:       <FloatArray> distance in Au
    """

    return astronomical_units * 1.496e11