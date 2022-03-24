""" This module provides the interior differential equations for the numerical shooting method of solving the
viscoelastic-gravitational problem within a planet.

The functions here wrap those found in the radial_derivatives_dynamic.py and radial_derivatives_static.py modules to
allow for a coarser-than-required interior model to be used. For example: the user might provide a planet's interior
with radius = np.linspace(0., 6.e6, 10) with shear_modulus and bulk_modulus defined at those steps. However, the
radial solution integrator may need more steps, or steps inbetween the ones provided.

These functions use a linear interpolation to fill in the missing steps.

"""

from typing import Tuple

import numpy as np

from . import (radial_derivatives_liquid_dynamic, radial_derivatives_liquid_static,
               radial_derivatives_solid_dynamic, radial_derivatives_solid_static)
from .....constants import G
from .....utilities.performance import njit
from .....utilities.types import FloatArray, NumArray


@njit(cacheable=True)
def dynamic_solid_ode(
    radius: FloatArray,
    y_vector: Tuple[NumArray, NumArray, NumArray, NumArray, NumArray, NumArray],
    radii: np.ndarray, shear_moduli: np.ndarray, bulk_moduli: np.ndarray,
    densities: np.ndarray, gravities: np.ndarray, frequency: float,
    order_l: int = 2, G_to_use: float = G
    ) -> Tuple[NumArray, NumArray, NumArray, NumArray, NumArray, NumArray]:
    """ A njit-safe radial derivative function for static, solid layers.

    Parameters
    ----------
    radius : FloatArray
        Requested radius to at which the radial derivatives are calculated [m]
    y_vector : Tuple[NumArray, NumArray, NumArray, NumArray, NumArray, NumArray]
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
    solid_dyn_derivatives : Tuple[NumArray, NumArray, NumArray, NumArray, NumArray, NumArray]
        The radial derivatives estimated using the provided radius for a solid layer under the dynamic assumption.
    """

    # These physical parameters are generally given in a coarser resolution than what is required during integration.
    #    Therefore, we interpolate their values at the integrator's requested `radius`.
    shear_modulus = np.interp(radius, radii, shear_moduli)
    bulk_modulus = np.interp(radius, radii, bulk_moduli)
    density = np.interp(radius, radii, densities)
    gravity = np.interp(radius, radii, gravities)

    solid_dyn_derivatives = \
        radial_derivatives_solid_dynamic(
            radius, y_vector, shear_modulus, bulk_modulus, density, gravity, frequency,
            order_l=order_l, G_to_use=G_to_use
            )

    return solid_dyn_derivatives


@njit(cacheable=True)
def static_solid_ode(
    radius: FloatArray, y_vector: Tuple[NumArray, NumArray, NumArray, NumArray, NumArray, NumArray],
    radii: np.ndarray, shear_moduli: np.ndarray, bulk_moduli: np.ndarray,
    densities: np.ndarray, gravities: np.ndarray,
    order_l: int = 2, G_to_use: float = G
    ) -> Tuple[NumArray, NumArray, NumArray, NumArray, NumArray, NumArray]:
    """ A njit-safe radial derivative function for static, solid layers.

    Parameters
    ----------
    radius : FloatArray
        Requested radius to at which the radial derivatives are calculated [m]
    y_vector : Tuple[NumArray, NumArray, NumArray, NumArray, NumArray, NumArray]
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
    solid_static_derivatives : Tuple[NumArray, NumArray, NumArray, NumArray, NumArray, NumArray]
        The radial derivatives estimated using the provided radius for a solid layer under the static assumption.
    """

    # These physical parameters are generally given in a coarser resolution than what is required during integration.
    #    Therefore, we interpolate their values at the integrator's requested `radius`.
    shear_modulus = np.interp(radius, radii, shear_moduli)
    bulk_modulus = np.interp(radius, radii, bulk_moduli)
    density = np.interp(radius, radii, densities)
    gravity = np.interp(radius, radii, gravities)

    solid_static_derivatives = \
        radial_derivatives_solid_static(
            radius, y_vector, shear_modulus, bulk_modulus, density, gravity,
            order_l=order_l, G_to_use=G_to_use
            )

    return solid_static_derivatives


@njit(cacheable=True)
def dynamic_liquid_ode(
    radius: FloatArray, y_vector: Tuple[NumArray, NumArray, NumArray, NumArray],
    radii: np.ndarray, bulk_moduli: np.ndarray,
    densities: np.ndarray, gravities: np.ndarray, frequency: float,
    order_l: int = 2, G_to_use: float = G
    ) -> callable:
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
    liquid_dyn_derivatives : Tuple[NumArray, NumArray, NumArray, NumArray]
        The radial derivatives estimated using the provided radius for a liquid layer under the dynamic assumption.
    """

    # These physical parameters are generally given in a coarser resolution than what is required during integration.
    #    Therefore, we interpolate their values at the integrator's requested `radius`.
    bulk_modulus = np.interp(radius, radii, bulk_moduli)
    density = np.interp(radius, radii, densities)
    gravity = np.interp(radius, radii, gravities)

    liquid_dyn_derivatives = radial_derivatives_liquid_dynamic(
        radius, y_vector, bulk_modulus, density, gravity, frequency,
        order_l=order_l, G_to_use=G_to_use
        )

    return liquid_dyn_derivatives


@njit(cacheable=True)
def static_liquid_ode(
    radius: FloatArray, y_vector: Tuple[NumArray, NumArray],
    radii: np.ndarray, densities: np.ndarray, gravities: np.ndarray,
    order_l: int = 2, G_to_use: float = G
    ) -> Tuple[NumArray, NumArray]:
    """ A njit-safe radial derivative function for static, liquid layers.

    Parameters
    ----------
    radius : FloatArray
        Requested radius to at which the radial derivatives are calculated [m]
    y_vector : Tuple[NumArray, NumArray]
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
    liquid_static_derivatives : Tuple[NumArray, NumArray]
        The radial derivatives estimated using the provided radius for a liquid layer under the static assumption.
    """

    # These physical parameters are generally given in a coarser resolution than what is required during integration.
    #    Therefore, we interpolate their values at the integrator's requested `radius`.
    density = np.interp(radius, radii, densities)
    gravity = np.interp(radius, radii, gravities)

    liquid_static_derivatives = radial_derivatives_liquid_static(
        radius, y_vector, density, gravity,
        order_l=order_l, G_to_use=G_to_use
        )

    return liquid_static_derivatives
