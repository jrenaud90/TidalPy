import numpy as np

from ..performance import njit
from ..types import FloatArray


# @njit
def complex_love(complex_compliance: FloatArray, shear_modulus: FloatArray, eff_rigidity: FloatArray) -> FloatArray:
    """ Calculates the 2nd order complex Love number

    :param complex_compliance: <FloatArray> Complex compliance (rheology based) [Pa-1]
    :param shear_modulus:      <FloatArray> Temperature modulated rigidity [Pa]
    :param eff_rigidity:       <FloatArray> 2nd order effective rigidity
    :return:                   <FloatArray> Complex Love Number
    """

    real_j = np.real(complex_compliance)
    imag_j = np.imag(complex_compliance)
    real_j2 = real_j**2
    imag_j2 = imag_j**2
    common_factor = (3. / 2.) * ((real_j + eff_rigidity / shear_modulus)**2 + imag_j2)**-1
    real_love = (real_j2 + imag_j2 + real_j * eff_rigidity / shear_modulus) * common_factor
    imag_love = (imag_j * eff_rigidity / shear_modulus) * common_factor
    complex_love = real_love + 1.0j * imag_love

    return complex_love


@njit
def complex_love_general(complex_compliance: FloatArray, shear_modulus: FloatArray, eff_rigidity_general: FloatArray,
                         order_l: int = 2) -> FloatArray:
    """ Calculates the l-th order complex Love number

    :param complex_compliance:      <FloatArray> Complex compliance (rheology based) [Pa-1]
    :param shear_modulus:           <FloatArray> Temperature modulated rigidity [Pa]
    :param eff_rigidity_general:    <FloatArray> l-th order effective rigidity
    :param order_l:                 <int> (optional) Outer-most fourier summation index
    :return:                        <FloatArray> Complex Love Number
    """

    real_j = np.real(complex_compliance)
    imag_j = np.imag(complex_compliance)
    real_j2 = real_j**2
    imag_j2 = imag_j**2
    common_factor = (3. / (order_l - 1.)) * ((real_j + eff_rigidity_general / shear_modulus)**2 + imag_j2)**-1
    real_love = (real_j2 + imag_j2 + real_j * eff_rigidity_general / shear_modulus) * common_factor
    imag_love = (imag_j * eff_rigidity_general / shear_modulus) * common_factor
    complex_love = real_love + 1.0j * imag_love

    return complex_love


# The functions below could, in theory, be cached but with shear_modulus potentially being an array lru_cache does not
#    work. Testing showed that a modified lru_cache is slower than just using njit (or the combination of the two).
#    Also, the modified cache would cause nasty bugs if an array were changed.
@njit
def effective_rigidity(shear_modulus: FloatArray, gravity: float, radius: float, density: float) -> FloatArray:
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
                               order_l: int = 2) -> FloatArray:
    """ Calculates the l-th order effective rigidity

    :param shear_modulus: <FloatArray> Temperature modulated rigidity
    :param gravity:       <float> Surface gravity [m s-2]
    :param radius:        <float> Surface radius [m]
    :param density:       <float> Bulk density [kg m-3]
    :param order_l:       <int> (optional) Outer-most fourier summation index
    :return:              <FloatArray> 2nd order Effective Rigidity
    """

    return (2. * order_l**2 + 4. * order_l + 3. / order_l) * shear_modulus / (gravity * radius * density)
