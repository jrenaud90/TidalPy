""" This module provides the interior differential equations for the numerical shooting method of solving the
viscoelastic-gravitational problem within a planet.

The functions here wrap those found in the radial_derivatives_dynamic.py and radial_derivatives_static.py modules to
allow for a coarser-than-required interior model to be used. For example: the user might provide a planet's interior
with radius = np.linspace(0., 6.e6, 10) with shear_modulus and bulk_modulus defined at those steps. However, the
radial solution integrator may need more steps, or steps inbetween the ones provided.

These functions use a linear interpolation to fill in the missing steps.

"""

import numpy as np

from TidalPy.constants import G
from TidalPy.utilities.performance import njit

from . import (
    radial_derivatives_liquid_dynamic, radial_derivatives_liquid_static,
    radial_derivatives_solid_dynamic, radial_derivatives_solid_static,
    radial_derivatives_liquid_dynamic_incomp, radial_derivatives_liquid_static_incomp,
    radial_derivatives_solid_dynamic_incomp, radial_derivatives_solid_static_incomp
)


@njit(cacheable=True)
def dynamic_solid_ode(
    radius: float,
    y_vector: np.ndarray,
    radius_array: np.ndarray, shear_modulus_array: np.ndarray, bulk_modulus_array: np.ndarray,
    density_array: np.ndarray, gravity_array: np.ndarray, frequency: float,
    order_l: int = 2, G_to_use: float = G, incompressible: bool = False
    ) -> np.ndarray:
    """ A njit-safe radial derivative function for static, solid layers.

    Parameters
    ----------
    radius : float
        Requested radius to at which the radial derivatives are calculated [m]
    y_vector : np.ndarray
        The radial functions at (or near) the provided radius. Used to estimate the derivatives.
    radius_array : np.ndarray
        Array of radius_array for the interior of a planet or layer where this radial derivative function is valid.
        The bottom most radius_array of a planet should not be equal to zero [m]
    shear_modulus_array : np.ndarray
        Shear modulus at each `radius_array` [Pa] (can be complex for shear dissipation)
    bulk_modulus_array : np.ndarray
        Bulk modulus at each `radius_array` [Pa] (can be complex for bulk dissipation)
    density_array : np.ndarray
        Density at each `radius_array` [kg m-3]
    gravity_array : np.ndarray
        Acceleration due to gravity calculated at each `radius_array` [m s-2]
    frequency : float
        Forcing frequency [rad s-1]
    order_l : int = 2
        Tidal harmonic order
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.
    incompressible : bool = False
        If `True`, the incompressible assumption will be used.

    Returns
    -------
    solid_dyn_derivatives : np.ndarray
        The radial derivatives estimated using the provided radius for a solid layer under the dynamic assumption.
    """

    # These physical parameters are generally given in a coarser resolution than what is required during integration.
    #    Therefore, we interpolate their values at the integrator's requested `radius`.
    shear_modulus = np.interp(radius, radius_array, shear_modulus_array)
    density = np.interp(radius, radius_array, density_array)
    gravity = np.interp(radius, radius_array, gravity_array)

    if incompressible:
        solid_dyn_derivatives = \
            radial_derivatives_solid_dynamic_incomp(
                radius, y_vector, shear_modulus, density, gravity, frequency,
                order_l=order_l, G_to_use=G_to_use
                )
    else:
        bulk_modulus = np.interp(radius, radius_array, bulk_modulus_array)
        solid_dyn_derivatives = \
            radial_derivatives_solid_dynamic(
                radius, y_vector, shear_modulus, bulk_modulus, density, gravity, frequency,
                order_l=order_l, G_to_use=G_to_use
                )

    return solid_dyn_derivatives


@njit(cacheable=True)
def static_solid_ode(
    radius: float, y_vector: np.ndarray,
    radius_array: np.ndarray, shear_modulus_array: np.ndarray, bulk_modulus_array: np.ndarray,
    density_array: np.ndarray, gravity_array: np.ndarray,
    order_l: int = 2, G_to_use: float = G, incompressible: bool = False
    ) -> np.ndarray:
    """ A njit-safe radial derivative function for static, solid layers.

    Parameters
    ----------
    radius : float
        Requested radius to at which the radial derivatives are calculated [m]
    y_vector : np.ndarray
        The radial functions at (or near) the provided radius. Used to estimate the derivatives.
    radius_array : np.ndarray
        Array of radius_array for the interior of a planet or layer where this radial derivative function is valid.
        The bottom most radius_array of a planet should not be equal to zero [m]
    shear_modulus_array : np.ndarray
        Shear modulus at each `radius_array` [Pa] (can be complex for shear dissipation)
    bulk_modulus_array : np.ndarray
        Bulk modulus at each `radius_array` [Pa] (can be complex for bulk dissipation)
    density_array : np.ndarray
        Density at each `radius_array` [kg m-3]
    gravity_array : np.ndarray
        Acceleration due to gravity calculated at each `radius_array` [m s-2]
    order_l : int = 2
        Tidal harmonic order
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.
    incompressible : bool = False
        If `True`, the incompressible assumption will be used.

    Returns
    -------
    solid_static_derivatives : np.ndarray
        The radial derivatives estimated using the provided radius for a solid layer under the static assumption.
    """

    # These physical parameters are generally given in a coarser resolution than what is required during integration.
    #    Therefore, we interpolate their values at the integrator's requested `radius`.
    shear_modulus = np.interp(radius, radius_array, shear_modulus_array)
    density = np.interp(radius, radius_array, density_array)
    gravity = np.interp(radius, radius_array, gravity_array)

    if incompressible:
        solid_static_derivatives = \
            radial_derivatives_solid_static_incomp(
                radius, y_vector, shear_modulus, density, gravity,
                order_l=order_l, G_to_use=G_to_use
                )
    else:
        bulk_modulus = np.interp(radius, radius_array, bulk_modulus_array)
        solid_static_derivatives = \
            radial_derivatives_solid_static(
                radius, y_vector, shear_modulus, bulk_modulus, density, gravity,
                order_l=order_l, G_to_use=G_to_use
                )

    return solid_static_derivatives


@njit(cacheable=True)
def dynamic_liquid_ode(
    radius: float, y_vector: np.ndarray,
    radius_array: np.ndarray, bulk_modulus_array: np.ndarray,
    density_array: np.ndarray, gravity_array: np.ndarray, frequency: float,
    order_l: int = 2, G_to_use: float = G, incompressible: bool = False
    ) -> np.ndarray:
    """ A njit-safe radial derivative function for dynamic, liquid layers.

    Parameters
    ----------
    radius : float
        Requested radius to at which the radial derivatives are calculated [m]
    y_vector : np.ndarray
        The radial functions at (or near) the provided radius. Used to estimate the derivatives.
    radius_array : np.ndarray
        Array of radius_array for the interior of a planet or layer where this radial derivative function is valid.
        The bottom most radius_array of a planet should not be equal to zero [m]
    bulk_modulus_array : np.ndarray
        Bulk modulus at each `radius_array` [Pa] (can be complex for bulk dissipation)
    density_array : np.ndarray
        Density at each `radius_array` [kg m-3]
    gravity_array : np.ndarray
        Acceleration due to gravity calculated at each `radius_array` [m s-2]
    frequency : float
        Forcing frequency [rad s-1]
    order_l : int = 2
        Tidal harmonic order
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.
    incompressible : bool = False
        If `True`, the incompressible assumption will be used.

    Returns
    -------
    liquid_dyn_derivatives : np.ndarray
        The radial derivatives estimated using the provided radius for a liquid layer under the dynamic assumption.
    """

    # These physical parameters are generally given in a coarser resolution than what is required during integration.
    #    Therefore, we interpolate their values at the integrator's requested `radius`.
    density = np.interp(radius, radius_array, density_array)
    gravity = np.interp(radius, radius_array, gravity_array)

    if incompressible:
        liquid_dyn_derivatives = radial_derivatives_liquid_dynamic_incomp(
            radius, y_vector, density, gravity, frequency,
            order_l=order_l, G_to_use=G_to_use
            )
    else:
        bulk_modulus = np.interp(radius, radius_array, bulk_modulus_array)
        liquid_dyn_derivatives = radial_derivatives_liquid_dynamic(
            radius, y_vector, bulk_modulus, density, gravity, frequency,
            order_l=order_l, G_to_use=G_to_use
            )

    return liquid_dyn_derivatives


@njit(cacheable=True)
def static_liquid_ode(
    radius: float, y_vector: np.ndarray,
    radius_array: np.ndarray, density_array: np.ndarray, gravity_array: np.ndarray,
    order_l: int = 2, G_to_use: float = G, incompressible: bool = False
    ) -> np.ndarray:
    """ A njit-safe radial derivative function for static, liquid layers.

    Parameters
    ----------
    radius : float
        Requested radius to at which the radial derivatives are calculated [m]
    y_vector : np.ndarray
        The radial functions at (or near) the provided radius. Used to estimate the derivatives.
    radius_array : np.ndarray
        Array of radius_array for the interior of a planet or layer where this radial derivative function is valid.
        The bottom most radius_array of a planet should not be equal to zero [m]
    density_array : np.ndarray
        Density at each `radius_array` [kg m-3]
    gravity_array : np.ndarray
        Acceleration due to gravity calculated at each `radius_array` [m s-2]
    order_l : int = 2
        Tidal harmonic order
    G_to_use : float = G
        Gravitational constant. Provide a non-dimensional version if the rest of the inputs are non-dimensional.
    incompressible : bool = False
        If `True`, the incompressible assumption will be used.

    Returns
    -------
    liquid_static_derivatives : np.ndarray
        The radial derivatives estimated using the provided radius for a liquid layer under the static assumption.
    """

    # These physical parameters are generally given in a coarser resolution than what is required during integration.
    #    Therefore, we interpolate their values at the integrator's requested `radius`.
    density = np.interp(radius, radius_array, density_array)
    gravity = np.interp(radius, radius_array, gravity_array)

    if incompressible:
        liquid_static_derivatives = radial_derivatives_liquid_static_incomp(
            radius, y_vector, density, gravity,
            order_l=order_l, G_to_use=G_to_use
            )
    else:
        liquid_static_derivatives = radial_derivatives_liquid_static(
            radius, y_vector, density, gravity,
            order_l=order_l, G_to_use=G_to_use
            )

    return liquid_static_derivatives
