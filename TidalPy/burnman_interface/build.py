import burnman
import numpy as np

from .material.helper import find_material
from .. import debug_mode
from ..exceptions import BadValueError, ParameterValueError, ParameterMissingError
from ..initialize import log
from ..types import float_eps
from ..configurations import force_burnman_quiet

burnman_verbose = debug_mode and not force_burnman_quiet

default_layer_params = {
    'material_source'     : None,
    'slices'              : 40,
    'temperature_mode'    : 'adiabatic',
    'fixed_temperature'   : None,
    'top_temperature'     : None,
    'interp_lookup_method': 'mid'
}


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
    layer_config = {**default_layer_params, **layer_config}

    # Check for missing parameters
    for param_name in ['radii', 'material', 'material_source', 'temperature_mode']:
        if param_name not in layer_config:
            raise ParameterMissingError(f'BurnMan Layer requires parameter: {param_name}')

    material = layer_config['material']
    material_source = layer_config['material_source']

    # Determine if the material is a composite or not
    composite = False
    if type(layer_config['material']) in [list, tuple]:
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
            # Burnman supports lists for this
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
    layer = burnman.Layer(name=layer_name, radii=radii, verbose=burnman_verbose)
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


def build_planet(planet_config: dict):
    """ Builds a BurnMan planet from a provided configuration dictionary.

    NOTE: Building a planet can be a slow process as it requires iteration between the various layers' EOS and the
     pressure. Therefore, it is recommended you make use of dilled or pickled planets to speed up initial computation
     time.

    Parameters
    ----------
    planet_config : dict
        Dictionary of planet configurations used to initialize layers and the planet

    Returns
    -------
    burnman_layers : list
        An ordered list of initialized burnman layers
    burnman_planet : burnman.Planet
        An initialized and built burnman planet class
    """

    log(f"Building planet: {planet_config['name']}", level='debug')

    # Store Layer information
    try:
        layers = planet_config['layers']
    except KeyError:
        raise ParameterMissingError('BurnMan Planet must have at least one layer.')

    # Set defaults if none are provided
    last_layer_radius = 0.
    for layer_i, (layer_name, layer_config_user) in enumerate(layers.items()):

        log(f"Initializing layer: {layer_name}", level='debug')

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
        layer_config['radii'] = np.linspace(layer_config['radius_lower'], layer_config['radius_upper'],
                                            layer_config['slices'])

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
    log(f"Building BurnMan Layers", level='debug')
    burnman_layers_byname = {layer_name: build_layer(layer_name, layer_data) for layer_name, layer_data in
                             layers.items()}
    burnman_layers = [burnman_layers_byname[layer_name] for layer_name in layers]

    # Build BurnMan Planet
    log(f"Initializing BurnMan Planet for {planet_config['name']}", level='debug')
    burnman_planet = burnman.Planet(planet_config['name'], burnman_layers, verbose=burnman_verbose)
    # Note: This section is slow!!
    log(f"Building BurnMan Planet for {planet_config['name']}", level='debug')
    burnman_planet.make()
    log('Planet construction completed', level='debug')

    return burnman_layers, burnman_planet
