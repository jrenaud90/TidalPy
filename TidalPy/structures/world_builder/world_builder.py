import copy
from typing import TextIO, Union

import toml

from TidalPy.paths import get_worlds_dir
from TidalPy.exceptions import (MissingArgumentError, NotYetImplementedError, TidalPyWorldError, UnknownWorld,
                                UnknownWorldType)
from TidalPy.utilities.dictionary_utils import nested_merge
from TidalPy.structures.world_builder.config_handler import clean_world_config, get_world_configs
from TidalPy.structures.world_types import world_types

from TidalPy.logger import get_logger
log = get_logger("TidalPy")


def build_world(world_name: str, world_config: Union[dict, TextIO] = None):
    """ Build a TidalPy world based on a pre-built config or a user-provided configuration dictionary

    Note: 'world' is used for any object: stars, gas giants, earth-like world_types, moons, and even Pluto!

    TidalPy ships with a directory of world configurations. By default these can be found in:
        <TidalPy install directory>/TidalPy/structures/world_types/

    These are TOML files that contain all the information that TidalPy needs to build a new world. Use the pre-built
        configuration files as templates to create new ones (which should be saved in the same directory).

    Alternatively, you can make a python dictionary object that contains the same information and pass it to this
        function using the world_config argument.

    Parameters
    ----------
    world_name : str
        World's name. This is used to search the pre-built configuration files.
    world_config : dict = None
        Alternatively, the user can provide a python dictionary with all of the required information.
        It is recommended that a pre-built configuration is used as a template for this dictionary.

    Returns
    -------
    world :
        The initialized TidalPy world object.
    """

    log.debug(f'Preparing to build world: {world_name}.')

    # If world_config is a file then load it through toml and get a dict
    if world_config is not None:
        if type(world_config) != dict:
            log.debug(f'Converting user provided planet configuration file to dictionary for {world_name}.')
            world_config = toml.load(world_config)

        # Make a copy of the dict so any subsequent changes do not affect the original dictionary
        world_config = copy.deepcopy(world_config)
        log.debug(f'User provided world configurations for {world_name}.')

    if world_config is None:
        log.debug(
            f'No manual configuration dictionary was provided for {world_name}, '
            f'attempting to locate saved configuration file.'
            )

        known_worlds = get_world_configs()

        if world_name in known_worlds:
            world_config = known_worlds[world_name]
        elif world_name.lower() in known_worlds:
            world_config = known_worlds[world_name.lower()]
        else:
            raise UnknownWorld(
                f'The user provided world name, {world_name}, can not be found in the directory of pre-built '
                f'world configs. Please add a new config to this directory or provide a manual world '
                f'configuration dictionary. Pre-built world configs can be found in:\n{get_worlds_dir()}'
                )

        log.debug(f'World configuration dictionary found for {world_name}.')

    # Determine world's type
    world_type = world_config['type'].lower()

    if world_type == 'burnman':
        log.debug(f'BurnMan world type detected.')
        from TidalPy.Extending.burnman import build_burnman_world, BurnManWorld
        
        log.debug('Attempting to build the BurnMan class for the world.')
        burnman_world, burnman_layers = build_burnman_world(world_config)
        log.debug(f'Burnman world building completed!')

        log.debug('Installing Burnman world into TidalPy world class.')
        world = BurnManWorld(world_config, burnman_world, burnman_layers, name=world_name)

    elif world_type not in world_types:
        raise UnknownWorldType(f'The world type, {world_type}, for {world_name} is unknown or not yet implemented.')
    
    else:
        # Get the world class
        WorldClass = world_types[world_type]
        log.debug(f'{WorldClass.world_class} world type detected.')
        world = WorldClass(world_config, world_name)

    if world is None:
        raise TidalPyWorldError('World building failed.')

    log.debug('World building was successful!')

    # TODO: Add back in this killing feature?
    # if configurations['exit_planets']:
    #     atexit.register(planet.kill_world)

    return world


def build_from_world(old_world, new_config: dict, new_name: str = None):
    """ Constructs a new world based on a previously built one.

    The old_world should have characteristics close to the one you desire. For example, if you wanted to make an Io-like
        world except that it has an ice layer overtop of it, then choose Europa as a starting point since it would have
        a similar layer structure and size to the desired world.

    If instead, you wish to simply scale up or down a world (turn Earth into a super-Earth or Io into a super-io,
        then look into using the scale_from_world function instead.

    See Also
    --------
    'world_builder.py'.scale_from_world

    Parameters
    ----------
    old_world :
        Already initialized world object.
    new_config : dict
        Planet config for the changes to the new planet. This will override configs found in old_planet.
    new_name : str
        Optional name provided to the new planet.

    Returns
    -------
    new_world :
        The newly configured and initialized planet object.
    """

    # Make a deep copy of the original planet's config
    log.debug(f'Building a new world based off of {old_world.name}.')
    old_config = old_world.config

    # Delete items that are most likely going to change
    old_config_copy = clean_world_config(old_config, make_copy=True)

    # Combine new and old dictionaries allowing the new dict to over write the old
    combo_dict = nested_merge(old_config_copy, new_config, make_copies=True)

    # Make any additional changes to the configs before planet build
    variant = False
    if new_name is None:
        if 'name' in new_config:
            # There is a name in the new configuration file. Is it the same as the previous world's name?
            new_name = new_config['name']
            if new_name == old_config_copy['name']:
                variant = True
        else:
            new_name = combo_dict['name']
            variant = True
    else:
        # User provided a new name as an argument, use it over anything else
        # Check if it matches the old world's name.
        if new_name == old_config_copy['name']:
            variant = True

    if variant:
        # Come up with a new name so that it is clear that this new planet was built from a different one
        if '_variant' in new_name:
            # Already was a variant - add a number
            world_name_pre = new_name.split('_variant')[0]
            i = 2
            while True:
                new_variant_name = f'{world_name_pre}_variant_{i}'
                if new_variant_name != new_name:
                    break
            new_name = new_variant_name
        else:
            new_name = f'{new_name}_variant'

    combo_dict['name'] = new_name

    # Build the world
    new_world = build_world(new_name, combo_dict)
    log.debug(f'{old_world.name} world built.')

    return new_world


def scale_from_world(old_world, new_name: str = None, mass_scale: float = None, radius_scale: float = None):
    """ Constructs a new planet that is a scaled up/down version of an old planet

    The old_planet should not be too different from your goal planet in terms of mass, radius, layer structure, and
        composition. For example, if you are trying to make a slightly larger/smaller Earth, don't use Io or Europa as
        for old_planet. Likewise, if you want to make a 2x mass Pluto, don't use the Earth as your base.

    Making sure that the layer structure of old_world is the same as the desired planet is more important than
        ensuring that old_world's radius/mass is close to the new world_types. So, if you want to make a super-sized Europa
        then it is actually better to start from Ganymede as it will have a high-pressure ice layer that Europa does
        not, but the new larger super-Europa might.

    It is recommended to call the new planet's .paint() method before using it in for any analysis to ensure its
        internal structure matches your expectations.

    Parameters
    ----------
    old_world :
        World to base the scale off of. This should be an already initialized TidalPy world object
    new_name : str = None
        New name to call the planet
    mass_scale : FloatNone = None
        Providing this will treat mass fractions of the planet's layer's as constants.
        Radii will be estimated from density
    radius_scale : FloatNone = None
        Providing this will treat volume fractions of the planet's layer's as constants.
        Masses will be calculated from EOS

    Returns
    -------
    new_world :
        The newly scaled planet.
    """

    if mass_scale is None:
        if radius_scale is None:
            raise MissingArgumentError(
                'radius_scale must be provided to the scale_from_world function '
                '(until mass_scale argument is implemented).'
                )
    else:
        # TODO: This may require some sort of iterative method w/ Burnman
        raise NotYetImplementedError('Right now planets can only be scaled by their radius')

    log.debug(f'Creating a new world by scaling {old_world.name}, with radius scale of {radius_scale}.')

    # Make a copy of the old world's configuration file
    scaled_config = clean_world_config(old_world.config, make_copy=True)
    old_name = scaled_config['name']

    # Change layer radius
    scaled_config['radius'] = radius_scale * old_world.config['radius']
    prev_layer_radius = 0.
    for layer_name, layer_dict in old_world.config['layers'].items():

        old_radius = layer_dict['radius']
        scaled_config['layers'][layer_name]['radius'] = radius_scale * old_radius
        scaled_config['layers'][layer_name]['radius_inner'] = prev_layer_radius

        # Use this layer's upper radius as the next layer's lower radius
        prev_layer_radius = scaled_config['layers'][layer_name]['radius']

        # Update other items
        scaled_config['layers'][layer_name]['thickness'] = scaled_config['layers'][layer_name]['radius'] - \
                                                           scaled_config['layers'][layer_name]['radius_inner']

    # Give it an identifiable name
    if new_name is None:
        if radius_scale >= 1:
            new_name = f'super-{old_name}'
        else:
            new_name = f'mini-{old_name}'

    # Make new planet based on the scaled dictionary
    new_planet = build_from_world(old_world, new_config=scaled_config, new_name=new_name)

    return new_planet
