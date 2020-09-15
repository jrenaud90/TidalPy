from typing import Union

import numpy as np

from ...exceptions import ParameterMissingError
from ...utilities.types import NoneType


def find_geometry_from_config(config: dict, layer_index: int, is_top_layer: bool,
                              world_radius: float, world_mass: float,
                              layer_below_radius: Union[float, NoneType] = None):
    """ Parse a configuration dictionary for geometry information

    Parameters
    ----------
    config : dict
        Object configuration dictionary.
    layer_index : int
        Index of where the layer is inside a world.
    is_top_layer : float
        If `True`, this layer is the top-most layer.
    world_radius : float
        World's radius [m]
    world_mass : float
        World's mass [kg]
    layer_below_radius : Union[float, NoneType] = None
        The layer below this one's radius [m]

    Returns
    -------
    radius : float
        Object's radius [m]
    thickness : float
        Object's thickness [m]
    volume : float
        Object's volume [m3]
    mass : float
        Object's mass [kg]
    density : float
        Object's density [kg m-3]
    """

    # Physical and geometric properties set up on initial call.
    radius = config.get('radius', None)
    mass = config.get('mass', None)
    thickness = config.get('thickness', None)
    density = config.get('density', None)
    if density is None:
        density = config.get('density_bulk', None)
    mass_frac = config.get('mass_frac', None)

    # Try to determine the layer's physical geometry
    geo_fail = False
    if radius is None:
        if thickness is None:
            if is_top_layer:
                if layer_below_radius is None:
                    geo_fail = True
                else:
                    radius = world_radius
                    thickness = radius - layer_below_radius
            else:
                geo_fail = True
        else:
            if layer_index == 0:
                radius = thickness
            else:
                if layer_below_radius is None:
                    geo_fail = True
                else:
                    radius = layer_below_radius + thickness
    else:
        if thickness is None:
            if layer_index == 0:
                thickness = radius
            else:
                if layer_below_radius is None:
                    geo_fail = True
                else:
                    thickness = radius - layer_below_radius
    if geo_fail:
        raise ParameterMissingError(f'Not enough information provided to determine geometry for layer.')
    volume = (4. / 3.) * np.pi * (radius**3 - (radius - thickness)**3)

    # Try to determine the layer's mass
    mass_fail = False
    if mass is None:
        if density is None:
            if mass_frac is None:
                mass_fail = True
            else:
                # Mass fraction provided.
                if world_mass is None:
                    mass_fail = True
                else:
                    mass = world_mass * mass_frac
        else:
            mass = density * volume
    if mass_fail:
        raise ParameterMissingError(f'Not enough information provided to calculate mass for layer.')

    return radius, thickness, volume, mass, density