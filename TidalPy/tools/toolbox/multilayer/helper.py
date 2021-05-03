import numpy as np

from ....utilities.performance import njit
from ....tides.multilayer.numerical_int import radial_derivatives_liquid_dynamic, radial_derivatives_liquid_static, \
    radial_derivatives_solid_dynamic, radial_derivatives_solid_static, RadialFuncSolidStaticType, \
    RadialFuncSolidDynamicType, RadialFuncLiquidDynamicType, RadialFuncLiquidStaticType


def build_dynamic_solid_solver(radii: np.ndarray, shear_moduli: np.ndarray, bulk_moduli: np.ndarray,
                               densities: np.ndarray, gravities: np.ndarray, frequency: float,
                               order_l: int = 2) -> callable:
    """ Build a njit-safe dynamic radial derivative function for solid layers.

    Parameters
    ----------
    radii : np.ndarray
        Array of radii for the interior of a planet or layer. The bottom most layer should not be equal to zero [m]
    shear_moduli : np.ndarray
        Shear modulus at each `radii` [Pa] (can be complex for shear dissipation)
    bulk_moduli : np.ndarray
        Bulk modulus at each `radii` [Pa] (can be complex for bulk dissipation)
    densities : np.ndarray
        Density at each `radii` [kg m-3]
    gravities : np.ndarray
        Acceleration due to gravity calculated at each `radii` [m s-2]
    frequency : float
        Forcing frequency [rad s-1]
    order_l : int = 2
        Tidal harmonic order

    Returns
    -------
    solid_dynamic_derivatives : callable
        Callable (njit-safe) function for calculating the radial derivatives for a solid layer using the dynamic
        assumption.
    """

    @njit(cacheable=False)
    def solid_dynamic_derivatives(radius: float, y_vector: RadialFuncSolidDynamicType) -> RadialFuncSolidDynamicType:
        shear_modulus = np.interp(radius, radii, shear_moduli)
        bulk_modulus = np.interp(radius, radii, bulk_moduli)
        density = np.interp(radius, radii, densities)
        gravity = np.interp(radius, radii, gravities)

        y_derivatives = radial_derivatives_solid_dynamic(radius, y_vector,
                                                         shear_modulus, bulk_modulus, density, gravity,
                                                         frequency, order_l=order_l)
        return y_derivatives

    return solid_dynamic_derivatives


def build_dynamic_liquid_solver(radii: np.ndarray, bulk_moduli: np.ndarray,
                                densities: np.ndarray, gravities: np.ndarray, frequency: float,
                                order_l: int = 2) -> callable:
    """ Build a njit-safe dynamic radial derivative function for liquid layers.

        Parameters
        ----------
        radii : np.ndarray
            Array of radii for the interior of a planet or layer. The bottom most layer should not be equal to zero [m]
        bulk_moduli : np.ndarray
            Bulk modulus at each `radii` [Pa] (can be complex for bulk dissipation)
        densities : np.ndarray
            Density at each `radii` [kg m-3]
        gravities : np.ndarray
            Acceleration due to gravity calculated at each `radii` [m s-2]
        frequency : float
            Forcing frequency [rad s-1]
        order_l : int = 2
            Tidal harmonic order

        Returns
        -------
        liquid_dynamic_derivatives : callable
            Callable (njit-safe) function for calculating the radial derivatives for a liquid layer using the dynamic
            assumption.
        """

    @njit(cacheable=False)
    def liquid_dynamic_derivatives(radius: float, y_vector: RadialFuncLiquidDynamicType) -> RadialFuncLiquidDynamicType:
        bulk_modulus = np.interp(radius, radii, bulk_moduli)
        density = np.interp(radius, radii, densities)
        gravity = np.interp(radius, radii, gravities)

        y_derivatives = radial_derivatives_liquid_dynamic(radius, y_vector,
                                                          bulk_modulus, density, gravity,
                                                          frequency, order_l=order_l)
        return y_derivatives

    return liquid_dynamic_derivatives


def build_static_solid_solver(radii: np.ndarray, shear_moduli: np.ndarray, bulk_moduli: np.ndarray,
                              densities: np.ndarray, gravities: np.ndarray,
                              order_l: int = 2) -> callable:
    """ Build a njit-safe static radial derivative function for solid layers.

    Parameters
    ----------
    radii : np.ndarray
        Array of radii for the interior of a planet or layer. The bottom most layer should not be equal to zero [m]
    bulk_moduli : np.ndarray
        Bulk modulus at each `radii` [Pa] (can be complex for bulk dissipation)
    densities : np.ndarray
        Density at each `radii` [kg m-3]
    gravities : np.ndarray
        Acceleration due to gravity calculated at each `radii` [m s-2]
    order_l : int = 2
        Tidal harmonic order

    Returns
    -------
    solid_static_derivatives : callable
        Callable (njit-safe) function for calculating the radial derivatives for a solid layer using the static
        assumption.
    """

    @njit(cacheable=False)
    def solid_static_derivatives(radius: float, y_vector: RadialFuncSolidStaticType) -> RadialFuncSolidStaticType:
        shear_modulus = np.interp(radius, radii, shear_moduli)
        bulk_modulus = np.interp(radius, radii, bulk_moduli)
        density = np.interp(radius, radii, densities)
        gravity = np.interp(radius, radii, gravities)

        y_derivatives = radial_derivatives_solid_static(radius, y_vector,
                                                        shear_modulus, bulk_modulus, density, gravity,
                                                        order_l=order_l)
        return y_derivatives

    return solid_static_derivatives


def build_static_liquid_solver(radii: np.ndarray, bulk_moduli: np.ndarray,
                               densities: np.ndarray, gravities: np.ndarray,
                               order_l: int = 2) -> callable:
    """ Build a njit-safe static radial derivative function for liquid layers.

    Parameters
    ----------
    radii : np.ndarray
        Array of radii for the interior of a planet or layer. The bottom most layer should not be equal to zero [m]
    bulk_moduli : np.ndarray
        Bulk modulus at each `radii` [Pa] (can be complex for bulk dissipation)
    densities : np.ndarray
        Density at each `radii` [kg m-3]
    gravities : np.ndarray
        Acceleration due to gravity calculated at each `radii` [m s-2]
    order_l : int = 2
        Tidal harmonic order

    Returns
    -------
    liquid_static_derivatives : callable
        Callable (njit-safe) function for calculating the radial derivatives for a liquid layer using the static
        assumption.
    """

    @njit(cacheable=False)
    def liquid_static_derivatives(radius: float, y_vector: RadialFuncLiquidStaticType) -> RadialFuncLiquidStaticType:
        bulk_modulus = np.interp(radius, radii, bulk_moduli)
        density = np.interp(radius, radii, densities)
        gravity = np.interp(radius, radii, gravities)

        y_derivatives = radial_derivatives_liquid_static(radius, y_vector,
                                                         bulk_modulus, density, gravity,
                                                         order_l=order_l)
        return y_derivatives

    return liquid_static_derivatives
