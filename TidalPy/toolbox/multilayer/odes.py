import numpy as np

from ...constants import G
from ...tides.multilayer.numerical_int import radial_derivatives_liquid_dynamic, radial_derivatives_liquid_static, \
    radial_derivatives_solid_dynamic, radial_derivatives_solid_static
from ...utilities.performance import njit
from ...utilities.types import FloatArray


@njit(cacheable=True)
def dynamic_solid_ode(radius: FloatArray, y_vector: np.ndarray,
                      radii: np.ndarray, shear_moduli: np.ndarray, bulk_moduli: np.ndarray,
                      densities: np.ndarray, gravities: np.ndarray, frequency: float,
                      order_l: int = 2, G_to_use: float = G) -> callable:
    """ A njit-safe radial derivative function for static, solid layers.

    Parameters
    ----------
    radius : FloatArray
        Requested radius to at which the radial derivatives are calculated [m]
    y_vector : np.ndarray
        The radial functions at (or near) the provided radius. Used to estimate the derivatives.
    radii : np.ndarray
        Array of radii for the interior of a planet or layer where this radial derivative function is valid.
        The bottom most radii of a planet should not be equal to zero [m]
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
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.

    Returns
    -------
    solid_dynamic_derivatives : np.ndarray
        The radial derivatives estimated using the provided radius for a solid layer under the dynamic assumption.
    """

    # These physical parameters are generally given in a coarser resolution than what is required during integration.
    #    Therefore, we interpolate their values at the integrator's requested `radius`.
    shear_modulus = np.interp(radius, radii, shear_moduli)
    bulk_modulus = np.interp(radius, radii, bulk_moduli)
    density = np.interp(radius, radii, densities)
    gravity = np.interp(radius, radii, gravities)

    y_derivatives = \
        radial_derivatives_solid_dynamic(
                radius, y_vector, shear_modulus, bulk_modulus, density, gravity, frequency,
                order_l=order_l, G_to_use=G_to_use)

    return y_derivatives


@njit(cacheable=True)
def static_solid_ode(radius: FloatArray, y_vector: np.ndarray,
                     radii: np.ndarray, shear_moduli: np.ndarray, bulk_moduli: np.ndarray,
                     densities: np.ndarray, gravities: np.ndarray,
                     order_l: int = 2, G_to_use: float = G) -> callable:
    """ A njit-safe radial derivative function for static, solid layers.

    Parameters
    ----------
    radius : FloatArray
        Requested radius to at which the radial derivatives are calculated [m]
    y_vector : np.ndarray
        The radial functions at (or near) the provided radius. Used to estimate the derivatives.
    radii : np.ndarray
        Array of radii for the interior of a planet or layer where this radial derivative function is valid.
        The bottom most radii of a planet should not be equal to zero [m]
    shear_moduli : np.ndarray
        Shear modulus at each `radii` [Pa] (can be complex for shear dissipation)
    bulk_moduli : np.ndarray
        Bulk modulus at each `radii` [Pa] (can be complex for bulk dissipation)
    densities : np.ndarray
        Density at each `radii` [kg m-3]
    gravities : np.ndarray
        Acceleration due to gravity calculated at each `radii` [m s-2]
    order_l : int = 2
        Tidal harmonic order
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.

    Returns
    -------
    solid_static_derivatives : np.ndarray
        The radial derivatives estimated using the provided radius for a solid layer under the static assumption.
    """

    # These physical parameters are generally given in a coarser resolution than what is required during integration.
    #    Therefore, we interpolate their values at the integrator's requested `radius`.
    shear_modulus = np.interp(radius, radii, shear_moduli)
    bulk_modulus = np.interp(radius, radii, bulk_moduli)
    density = np.interp(radius, radii, densities)
    gravity = np.interp(radius, radii, gravities)

    y_derivatives = \
        radial_derivatives_solid_static(
                radius, y_vector, shear_modulus, bulk_modulus, density, gravity,
                order_l=order_l, G_to_use=G_to_use)

    return y_derivatives


@njit(cacheable=True)
def dynamic_liquid_ode(radius: FloatArray, y_vector: np.ndarray,
                       radii: np.ndarray, bulk_moduli: np.ndarray,
                       densities: np.ndarray, gravities: np.ndarray, frequency: float,
                       order_l: int = 2, G_to_use: float = G) -> callable:
    """ A njit-safe radial derivative function for dynamic, liquid layers.

    Parameters
    ----------
    radius : FloatArray
        Requested radius to at which the radial derivatives are calculated [m]
    y_vector : np.ndarray
        The radial functions at (or near) the provided radius. Used to estimate the derivatives.
    radii : np.ndarray
        Array of radii for the interior of a planet or layer where this radial derivative function is valid.
        The bottom most radii of a planet should not be equal to zero [m]
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
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.

    Returns
    -------
    solid_dynamic_derivatives : np.ndarray
        The radial derivatives estimated using the provided radius for a liquid layer under the dynamic assumption.
    """

    # These physical parameters are generally given in a coarser resolution than what is required during integration.
    #    Therefore, we interpolate their values at the integrator's requested `radius`.
    bulk_modulus = np.interp(radius, radii, bulk_moduli)
    density = np.interp(radius, radii, densities)
    gravity = np.interp(radius, radii, gravities)

    y_derivatives = radial_derivatives_liquid_dynamic(
            radius, y_vector, bulk_modulus, density, gravity, frequency,
            order_l=order_l, G_to_use=G_to_use)

    return y_derivatives


@njit(cacheable=True)
def static_liquid_ode(radius: FloatArray, y_vector: np.ndarray,
                      radii: np.ndarray, densities: np.ndarray, gravities: np.ndarray,
                      order_l: int = 2, G_to_use: float = G) -> callable:
    """ A njit-safe radial derivative function for static, liquid layers.

    Parameters
    ----------
    radius : FloatArray
        Requested radius to at which the radial derivatives are calculated [m]
    y_vector : np.ndarray
        The radial functions at (or near) the provided radius. Used to estimate the derivatives.
    radii : np.ndarray
        Array of radii for the interior of a planet or layer where this radial derivative function is valid.
        The bottom most radii of a planet should not be equal to zero [m]
    densities : np.ndarray
        Density at each `radii` [kg m-3]
    gravities : np.ndarray
        Acceleration due to gravity calculated at each `radii` [m s-2]
    order_l : int = 2
        Tidal harmonic order
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.

    Returns
    -------
    solid_dynamic_derivatives : np.ndarray
        The radial derivatives estimated using the provided radius for a liquid layer under the static assumption.
    """

    # These physical parameters are generally given in a coarser resolution than what is required during integration.
    #    Therefore, we interpolate their values at the integrator's requested `radius`.
    density = np.interp(radius, radii, densities)
    gravity = np.interp(radius, radii, gravities)

    y_derivatives = radial_derivatives_liquid_static(
            radius, y_vector, density, gravity,
            order_l=order_l, G_to_use=G_to_use)

    return y_derivatives
