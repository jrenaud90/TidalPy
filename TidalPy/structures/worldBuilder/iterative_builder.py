""" Planet Iteration Package

    The goal of this package is to provide the user an easy, but less accurate, way to estiamte the interior structure
        of planets. It was designed with exoplanets in mind where some measurements (moment of inertia, etc) are
        completely unknown.

    TODO: This is *not* the fastest nor most accurate way to do this. It was a quick and dirty way to get some estimates
     In the future a more robust method should be developed. It should also be able to handle MOI, etc. for a more
     accurate answer. Ideally it would also be able to change the number/composition of layers to better fit the data.
     That latter part is probably a much bigger task

    TODO: This could be sped up by only calling on the burnman builder (instead of having to initialize all of TidalPy
     at each iteration step).

"""
import copy

import numpy as np

from . import build_from_world
from ... import log
from ...exceptions import IncompatibleModelError


def world_iterative_builder(base_world, goal_radius, goal_mass, ice_mass_frac: float = None,
                            initial_core_density: float = 10000., initial_mantle_density: float = 4000.,
                            ice_density: float = 920.,
                            tolerance = 0.05):
    """ Takes a baseline planet that is close to ideal and iterates on the constructor until there is convergence in
    the mass and radius

    Parameters
    ----------
    base_world :
        Instantiated world object that acts as the starting point of the iteration
    goal_radius : float
        The final desired radius
    goal_mass : float
        The final desired mass
    ice_mass_frac : float
        The mass fraction of ice (or other envelope, just change the density) on the planet.
    initial_core_density : float
        Initial (average) core density used in the first iteration. The closer this is to the final value, the faster
        the iteration.
    initial_mantle_density : float
        Initial (average) mantle density used in the first iteration. The closer this is to the final value, the faster
        the iteration.
    ice_density : float
        Density of the ice layer
    tolerance : float
        The tolerance on final vs. goal mass and radius.

    Returns
    -------
    iterated_world :
        Final, instantiated, planet which the iterated mass and radius values
    """

    log.debug(f'Building a planet based on an iteration on planet: {base_world.name}')
    log.warning(f'The world_iterative_builder is very much a work in progress. '
                f'Double check that results match expectations.')

    # This is a two-layer iterator so a ice_mass_fraction is required. If none is provided then we will use whatever
    #    value was initially used in base_world.
    num_layers = len(base_world.layers)

    if num_layers < 2 or num_layers > 3:
        log.error(f'world_iterative_builder requires at least 2, and at most 3, layers in the '
                  f'base_world, {base_world.name}.')
        raise IncompatibleModelError('world_iterative_builder requires at least 2, and at most 3, layers in the base_world.')
    if num_layers == 2:
        # No ice layer is used
        ice_mass_frac = None
    if num_layers == 3 and ice_mass_frac is None:
        ice_mass_frac = base_world.layers[-1].mass / base_world.mass

    # This iterator always assumes that the ice shell is the same mass.
    if ice_mass_frac is None:
        ice_mass = 0.
    else:
        ice_mass = goal_mass * ice_mass_frac
    ice_volume = ice_mass / ice_density
    # Assume that the ice shell is the top-most layer
    if ice_volume == 0.:
        ice_radius_upper = None
        ice_radius_lower = None
    else:
        ice_radius_upper = goal_radius
        ice_radius_lower = (ice_radius_upper**3 - (3. * ice_volume / (4. * np.pi)))**(1/3)

    # Make a copy of the base planet's dictionary
    new_config = copy.deepcopy(base_world.config)
    new_config['radius'] = goal_radius
    bulk_density = goal_mass / ((4. / 3.) * np.pi * goal_radius**3)

    # Ge the names of the various layers
    core_layer_name = base_world.layers[0].name
    mantle_layer_name = base_world.layers[1].name
    if ice_mass_frac is not None:
        ice_layer_name = base_world.layers[2].name
    else:
        ice_layer_name = None

    new_config['layers'][core_layer_name]['radius_lower'] = 0.
    # Add the ice and mantle layer information now
    if ice_radius_upper is not None:
        new_config['layers'][ice_layer_name]['radius'] = ice_radius_upper
        new_config['layers'][ice_layer_name]['radius_upper'] = ice_radius_upper
        new_config['layers'][ice_layer_name]['radius_lower'] = ice_radius_lower

        mantle_radius = ice_radius_lower
        new_config['layers'][mantle_layer_name]['radius'] = mantle_radius
        new_config['layers'][mantle_layer_name]['radius_upper'] = mantle_radius
    else:
        # No Ice Shell
        ice_density = 0.
        mantle_radius = goal_radius
        new_config['layers'][mantle_layer_name]['radius'] = goal_radius
        new_config['layers'][mantle_layer_name]['radius_upper'] = goal_radius

    # Now use a two layer system
    core_density = initial_core_density
    mantle_density = initial_mantle_density
    iterated_world = base_world
    iteration_step = 0
    while abs(goal_mass - iterated_world.mass) / goal_mass > tolerance:
        iteration_step += 1
        log.debug(f'\tIteration Step {iteration_step}')
        if iteration_step > 40:
            log.error(f'world_iterative_builder could not find convergence working from the '
                      f'base_world, {base_world.name}.')
            raise StopIteration('world_iterative_builder could not find convergence.')

        # Pull out information from the burnman build
        mantle_density = iterated_world.layers[1].density_bulk
        core_density = iterated_world.layers[0].density_bulk

        core_radius = ((goal_radius**3 * (bulk_density - ice_density) +
                        mantle_radius**3 * (ice_density - mantle_density)) / (core_density - mantle_density))**(1/3)

        # Update config with new values
        new_config['layers'][core_layer_name]['radius'] = core_radius
        new_config['layers'][core_layer_name]['radius_upper'] = core_radius

        new_config['layers'][mantle_layer_name]['radius_lower'] = core_radius

        # Build new planet from the updated config
        iterated_world = build_from_world(iterated_world, new_config, new_name=base_world.name)

    return iterated_world







