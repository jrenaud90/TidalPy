from TidalPy.exceptions import ParameterMissingError, UnusualRealValueError, BadValueError
from TidalPy.types import float_eps

from .material.helper import find_material
from .. import debug_mode
from ..configurations import save_burnman_planets_to_run_folder, save_burnman_planets_to_tidalpy_dir

import burnman
import numpy as np

default_layer_params = {
    'material_source': None,
    'slices': 40,
    'temperature_mode': 'adiabatic',
    'fixed_temperature': None,
    'top_temperature': None,
    'interp_lookup_method': 'mid'
}

def build_layer(layer_name: str, layer_config: dict):
    """ Builds a BurnMan layer from a provided configuration dictionary.

    """

    # Load in defaults
    layer_config = {**default_layer_params, **layer_config}

    # Check for missing parameters
    for param_name in ['radii', 'material', 'material_source', 'temperature_mode']:
        if param_name not in layer_config:
            raise ParameterMissingError(f'BurnMan Layer requires parameter: {param_name}')

    # Find BurnMan material based on material and source name
    material_class = find_material(layer_config['material'], layer_config['material_source'])

    # Initialize material
    material = material_class()

    # Make layer and load in material
    radii = layer_config['radii']
    layer = burnman.Layer(name=layer_name, radii=radii, verbose=debug_mode)
    layer.set_material(material)

    # Load in temperature model
    if layer_config['temperature_mode'] == 'adiabatic':
        if layer_config.get('temperature_top', None) is None:
            layer.set_temperature_mode('adiabatic')
        else:
            layer.set_temperature_mode('adiabatic', temperature_top=layer_config['temperature_top'])
    elif layer_config['temperature_mode'] == 'user-defined':
        layer.set_temperature_mode('user-defined', layer_config['temperature_fixed'] * np.ones_like(radii))
    else:
        layer.set_temperature_mode(layer_config['temperature_mode'])

    return layer


def build_planet(planet_config: dict):
    """ Builds a BurnMan planet from a provided configuration dictionary.

    NOTE: This is a slow process and should not be repeated unless physical parameters of the planet are changing.
    """

    # Store Layer information
    try:
        layers = planet_config['layers']
    except KeyError:
        raise ParameterMissingError('BurnMan Planet must have at least one layer.')

    # Set defaults if none are provided
    last_layer_radius = 0.
    for layer_i, (layer_name, layer_config_user) in enumerate(layers.items()):
        # Check if layer has must-have parameters
        for param_name in ['type', 'material', 'radius']:
            if param_name not in layer_config_user:
                raise ParameterMissingError(f'BurnMan Layer requires parameter: {param_name}')

        # Load in defaults for non-required parameters if not present
        layer_config = {**default_layer_params, **layer_config_user}

        # Load in calculated parameters
        if 'thickness' not in layer_config:
            if layer_i == 0:
                layer_config['thickness'] = layer_config['radius']
            else:
                layer_config['thickness'] = layer_config['radius'] - last_layer_radius
        layer_config['radius_upper'] = layer_config['radius']
        layer_config['radius_lower'] = layer_config['radius_upper'] - layer_config['thickness']
        layer_config['radii'] = np.linspace(layer_config['radius_lower'], layer_config['radius_upper'], layer_config['slices'])

        # Check for physical consistency
        if layer_config['thickness'] <= float_eps:
            raise BadValueError('Layer thickness must be non-zero')
        if layer_i > 0:
            if layer_config['radius_lower'] != last_layer_radius:
                raise BadValueError('Layers must stack on top of one another')

        last_layer_radius = layer_config['radius']
        # Put the new config back into the layers storage
        layers[layer_name] = layer_config

    # Build BurnMan layers
    burnman_layers_byname = {layer_name: build_layer(layer_name, layer_data) for layer_name, layer_data in
                             layers.items()}
    burnman_layers = [burnman_layers_byname[layer_name] for layer_name in layers]

    # Build BurnMan Planet
    burnman_planet = burnman.Planet(planet_config['name'], burnman_layers, verbose=debug_mode)
    # Note: This section is slow!!
    burnman_planet.make()

    return burnman_layers, burnman_planet
