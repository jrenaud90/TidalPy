from numba import njit
from ..types import FloatArray

@njit
def complex_love(complex_compliance: FloatArray, shear_modulus: FloatArray, eff_rigidity: FloatArray):
    """ Calculates the 2nd order complex Love number

    :param complex_compliance: <FloatArray> Complex compliance (rheology based) [Pa-1]
    :param shear_modulus:      <FloatArray> Temperature modulated rigidity [Pa]
    :param eff_rigidity:       <FloatArray> 2nd order effective rigidity
    :return:                   <FloatArray> Complex Love Number
    """

    return (3. / 2.) * (1. + eff_rigidity / (shear_modulus * complex_compliance))**(-1)


@njit
def complex_love_general(complex_compliance: FloatArray, shear_modulus: FloatArray, eff_rigidity_general: FloatArray,
                         order_l: int = 2):
    """ Calculates the l-th order complex Love number

    :param complex_compliance:      <FloatArray> Complex compliance (rheology based) [Pa-1]
    :param shear_modulus:           <FloatArray> Temperature modulated rigidity [Pa]
    :param eff_rigidity_general:    <FloatArray> l-th order effective rigidity
    :param order_l:                 <int> (optional) Outer-most fourier summation index
    :return:                        <FloatArray> Complex Love Number
    """

    return (3. / (2. * (order_l - 1.))) * (1. + eff_rigidity_general / (shear_modulus * complex_compliance))**(-1)


# The functions below could, in theory, be cached but with shear_modulus potentially being an array lru_cache does not
#    work. Testing showed that a modified lru_cache is slower than just using njit (or the combination of the two).
#    Also, the modified cache would cause nasty bugs if an array were changed.
@njit
def effective_rigidity(shear_modulus: FloatArray, gravity: float, radius: float, density: float):
    """ Calculates the 2nd order effective rigidity

    :param shear_modulus: <FloatArray> Temperature modulated rigidity
    :param gravity:       <float> Surface gravity [m s-2]
    :param radius:        <float> Surface radius [m]
    :param density:       <float> Bulk density [kg m-3]
    :return:              <FloatArray> 2nd order Effective Rigidity
    """

    return (19. / 2.) * shear_modulus / (gravity * radius * density)

@njit
def effective_rigidity_general(shear_modulus: FloatArray, gravity: float, radius: float, density: float,
                               order_l: int = 2):
    """ Calculates the l-th order effective rigidity

    :param shear_modulus: <FloatArray> Temperature modulated rigidity
    :param gravity:       <float> Surface gravity [m s-2]
    :param radius:        <float> Surface radius [m]
    :param density:       <float> Bulk density [kg m-3]
    :param order_l:       <int> (optional) Outer-most fourier summation index
    :return:              <FloatArray> 2nd order Effective Rigidity
    """

    return (2. * order_l**2 + 4. * order_l + 3. / order_l) * shear_modulus / (gravity * radius * density)