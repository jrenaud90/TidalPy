from typing import Tuple

import burnman
import numpy as np

from .defaults import default_burnman_layer_params
from .material.helper import find_material
from .. import configurations, debug_mode, log
from ..exceptions import BadValueError, ParameterMissingError, ParameterValueError
from ..utilities.types import float_eps

BURNMAN_VERBOSE = debug_mode and not configurations['force_burnman_quiet']


def build_layer(layer_name: str, layer_config: dict) -> burnman.Layer:
    """ Build a Burnman layer

    Parameters
    ----------
    layer_name : str
        Name of the layer
    layer_config : dict
        Dictionary of layer configurations used to initialize the layer

    Returns
    -------
    layer : burnman.Layer
        An initialized burnman layer
    """

    # Load in defaults
    layer_config = {**default_burnman_layer_params, **layer_config}

    # Check for missing parameters
    for param_name in ['radii', 'material', 'material_source', 'temperature_mode']:
        if param_name not in layer_config:
            raise ParameterMissingError(f'BurnMan Layer requires parameter: {param_name}')

    material = layer_config['material']
    material_source = layer_config['material_source']

    # Determine if the material is a composite or not
    composite = False
    if type(layer_config['material']) in [list, tuple]:
        log.debug(f'Composite material encountered for layer: {layer_name}.')
        composite = True

    if composite:
        try:
            material_fractions = layer_config['material_fractions']
        except KeyError:
            raise ParameterMissingError(f"Multiple materials provided but no material_fraction's were given.")
        if type(material_fractions) not in [list, tuple]:
            raise ParameterValueError('Material fractions must be provided as a list or tuple.')
        if type(material_source) not in [list, tuple]:
            raise ParameterValueError('Multiple materials provided with only one source.')
        if len(material) != len(material_source):
            raise ParameterValueError('Number of composite materials is not the same as number of material sources.')
        if len(material) != len(material_fractions):
            raise ParameterValueError('Number of composite materials is not the same as number of material fractions.')

        if type(material_fractions) == tuple:
            # Burnman only supports lists for this
            material_fractions = list(material_fractions)

        material_instances = list()
        for mat_name, mat_source in zip(material, material_source):
            # Find BurnMan material based on material and source name
            material_class = find_material(mat_name, mat_source)
            material_instances.append(material_class())

        # Build composite material
        init_material = burnman.Composite(material_instances, material_fractions)

    else:
        # Find BurnMan material based on material and source name
        material_class = find_material(layer_config['material'], layer_config['material_source'])

        # Initialize material
        init_material = material_class()

    # Make layer and load in material
    radii = layer_config['radii']
    layer = burnman.Layer(name=layer_name, radii=radii, verbose=BURNMAN_VERBOSE)
    layer.set_material(init_material)

    # Don't keep radii in the config, otherwise it will save to the json file as a large list
    del layer_config['radii']

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


def build_burnman_world(planet_config: dict, verbose: bool = False) -> Tuple[burnman.Planet, Tuple[burnman.Layer, ...]]:
    """ Builds a BurnMan planet from a provided configuration dictionary.

    Parameters
    ----------
    planet_config : dict
        Dictionary of planet configurations used to initialize layers and the planet.
    verbose : bool = False
        If True, BurnMan will print more information to console.

    Returns
    -------
    burnman_world : burnman.Planet
        An initialized and built burnman planet class
    burnman_layers : Tuple[burnman.Layer, ...]
        An ordered list of initialized burnman layers
    """

    log.debug(f"Building BurnMan world: {planet_config['name']}")

    # Store Layer information
    try:
        layers = planet_config['layers']
    except KeyError:
        raise ParameterMissingError('BurnMan Planet must have at least one layer.')

    # Set defaults if none are provided
    last_layer_radius = 0.
    for layer_i, (layer_name, layer_config_user) in enumerate(layers.items()):

        log.debug(f"Initializing layer: {layer_name}")

        # Check if layer has must-have parameters
        for param_name in ['type', 'material', 'radius']:
            if param_name not in layer_config_user:
                raise ParameterMissingError(f'BurnMan Layer requires parameter: {param_name}')

        # Load in defaults for non-required parameters if not present
        layer_config = {**default_burnman_layer_params, **layer_config_user}

        # Load in calculated parameters
        if 'thickness' not in layer_config:
            if layer_i == 0:
                layer_config['thickness'] = layer_config['radius']
            else:
                layer_config['thickness'] = layer_config['radius'] - last_layer_radius
        layer_config['radius_upper'] = layer_config['radius']
        layer_config['radius_lower'] = layer_config['radius_upper'] - layer_config['thickness']
        layer_config['radii'] = np.linspace(
            layer_config['radius_lower'], layer_config['radius_upper'],
            layer_config['slices'], endpoint=True
            )

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
    log.debug(f"Building BurnMan Layers")
    burnman_layers_byname = {layer_name: build_layer(layer_name, layer_data) for layer_name, layer_data in
                             layers.items()}
    burnman_layers = [burnman_layers_byname[layer_name] for layer_name in layers]

    # Build BurnMan Planet
    log.debug(f"Initializing BurnMan Planet for {planet_config['name']}")
    verbose = BURNMAN_VERBOSE or verbose
    burnman_world = burnman.Planet(planet_config['name'], burnman_layers, verbose=verbose)
    # Note: This section can be slow!
    burnman_world.make()
    log.debug('BurnMan world construction completed')

    # Change layer list to tuple to emphasize that these should not change or mutate
    burnman_layers = tuple(burnman_layers)

    return burnman_world, burnman_layers
