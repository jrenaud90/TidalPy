from typing import TYPE_CHECKING

from TidalPy.utilities.performance.numba import njit

if TYPE_CHECKING:
    from TidalPy.utilities.types import FloatArray


@njit(cacheable=True)
def complex_love(complex_compliance: 'FloatArray', shear_modulus: 'FloatArray',
                 eff_rigidity: 'FloatArray') -> 'FloatArray':
    """ Calculates the 2nd order complex Love number

    Parameters
    ----------
    complex_compliance : FloatArray
        Complex compliance (rheology based) [Pa-1]
    shear_modulus : FloatArray
        Temperature modulated rigidity [Pa]
    eff_rigidity : FloatArray
        2nd order effective rigidity

    Returns
    -------
    cmplx_love : FloatArray
        2nd order complex Love number
    """

    rheology_factor = complex_compliance * shear_modulus
    cmplx_love = (3. / 2.) * (1. / (1. + (eff_rigidity / rheology_factor)))

    return cmplx_love


@njit(cacheable=True)
def complex_love_general(
    complex_compliance: 'FloatArray', shear_modulus: 'FloatArray', eff_rigidity_general: 'FloatArray',
    order_l: int = 2
    ) -> 'FloatArray':
    """ Calculates the l-th order complex Love number

    Parameters
    ----------
    complex_compliance : FloatArray
        Complex compliance (rheology based) [Pa-1]
    shear_modulus : FloatArray
        Temperature modulated rigidity [Pa]
    eff_rigidity_general : FloatArray
        l-th order effective rigidity
    order_l : int
        Outermost Fourier summation integer

    Returns
    -------
    cmplx_love : FloatArray
        l-th order complex Love number
    """

    rheology_factor = complex_compliance * shear_modulus
    cmplx_love = (3. / (2. * (order_l - 1))) * (1. / (1. + (eff_rigidity_general / rheology_factor)))

    return cmplx_love


@njit(cacheable=True)
def static_love(eff_rigidity: 'FloatArray') -> 'FloatArray':
    """ Calculate the static (non-complex) 2nd order Love number

    Parameters
    ----------
    eff_rigidity : FloatArray
        2nd order effective rigidity

    Returns
    -------
    static_love : FloatArray
        2nd order complex Love number
    """

    static_love_ = (3. / 2.) * (1. / (1. + eff_rigidity))
    return static_love_


@njit(cacheable=True)
def static_love_general(eff_rigidity_general: 'FloatArray', order_l: int = 2) -> 'FloatArray':
    """ Calculate the static (non-complex) tidal Love number k.

    Parameters
    ----------
    eff_rigidity_general : FloatArray
        l-th order effective rigidity
    order_l : int
        Outermost Fourier summation integer

    Returns
    -------
    static_love : FloatArray
        l-th order complex Love number
    """

    static_love = (3. / (2. * (order_l - 1))) * (1. / (1. + eff_rigidity_general))
    return static_love


@njit(cacheable=True)
def effective_rigidity(shear_modulus: 'FloatArray', gravity: float, radius: float, density: float) -> 'FloatArray':
    """ Calculates the 2nd order effective rigidity

    Parameters
    ----------
    shear_modulus : FloatArray
        Temperature modulated rigidity [Pa]
    gravity : float
        Surface (of planet or layer) gravity [m s-2]
    radius : float
        Surface (of planet or layer) radius [m]
    density : float
        Bulk density (of planet or layer) [kg m-3]

    Returns
    -------
    eff_rigid : FloatArray
        2nd order Effective Rigidity
    """

    eff_rigid = (19. / 2.) * shear_modulus / (gravity * radius * density)

    return eff_rigid


@njit(cacheable=True)
def effective_rigidity_general(
    shear_modulus: 'FloatArray', gravity: float, radius: float, density: float,
    order_l: int = 2
    ) -> 'FloatArray':
    """ Calculates the l-th order effective rigidity

    Parameters
    ----------
    shear_modulus : FloatArray
        Temperature modulated rigidity [Pa]
    gravity : float
        Surface (of planet or layer) gravity [m s-2]
    radius : float
        Surface (of planet or layer) radius [m]
    density : float
        Bulk density (of planet or layer) [kg m-3]
    order_l : int
        Outermost Fourier summation integer

    Returns
    -------
    eff_rigid : FloatArray
        l-th order Effective Rigidity
    """

    eff_rigid = (2. * order_l**2 + 4. * order_l + 3. / order_l) * shear_modulus / (gravity * radius * density)

    return eff_rigid
