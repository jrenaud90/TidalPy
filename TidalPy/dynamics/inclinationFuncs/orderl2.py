import numpy as np
from ...performance import njit

sin = np.sin
cos = np.cos

@njit
def calc_inclination_off(inclination):
    """ Calculate inclination functions squared for inclination == 0

    Parameters
    ----------
    inclination

    Returns
    -------
    output : dict[(l, m, p)] = F^2
    """

    # Output is structured (p0_(m0_result, m1_result, ...), p1_(...), ...)
    output = (
        (
            np.zeros_like(inclination),
            np.zeros_like(inclination),
            9. * np.ones_like(inclination)
        ),
        (
            (1. / 4.) * np.ones_like(inclination),
            np.zeros_like(inclination),
            np.zeros_like(inclination)
        ),
        (
            np.zeros_like(inclination),
            np.zeros_like(inclination),
            np.zeros_like(inclination)
        )
    )

    return output

@njit
def calc_inclination(inclination):
    """ Calculate inclination functions squared for a given inclination

    Parameters
    ----------
    inclination

    Returns
    -------
    output : Tuple[Tuple[float, float, float], ...]

    """

    # Optimizations
    sin_i_4 = sin(inclination)**4
    i_half = inclination / 2.
    i_double = 2. * inclination

    # Output is structured (p0_(m0_result, m1_result, ...), p1_(...), ...)
    output = (
        (
            (9. / 64.) * sin_i_4,
            9. * (sin(i_half)**2) * (cos(i_half)**6),
            9. * (cos(i_half)**8)
        ),
        (
            (3. * sin(inclination)**2 - 2.)**2,
            (9. / 16.) * (sin(i_double)**2),
            (9. / 4.) * sin_i_4
        ),
        (
            (9. / 64.)*sin_i_4,
            9. * (sin(i_half)**6) * (cos(i_half)**2),
            9. * (sin(i_half)**8)
        )
    )

    return output