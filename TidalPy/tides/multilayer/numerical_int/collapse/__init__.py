from typing import Tuple

import numpy as np

from .solid_solid_liquid_solid import collapse_ssls_dynamic_liq, collapse_ssls_static_liq
from .homogen_solid import collapse_homogen_solid
from .liquid_solid import collapse_ls_dynamic_liq, collapse_ls_static_liq
from .solid_liquid_solid import collapse_sls_dynamic_liq, collapse_sls_static_liq

# Known models are stored by the model name and then by if the liquid layer is static or not.
known_collapse_models = {
    # The homogen solid model is repeated because there is no liquid layer.
    'homogeneous_solid'       : {
        True : collapse_homogen_solid,
        False: collapse_homogen_solid
        },
    'liquid_solid'            : {
        True : collapse_ls_static_liq,
        False: collapse_ls_dynamic_liq
        },
    'solid_liquid_solid'      : {
        True : collapse_sls_static_liq,
        False: collapse_sls_dynamic_liq
        },
    'solid_solid_liquid_solid': {
        True : collapse_ssls_static_liq,
        False: collapse_ssls_dynamic_liq
        }
    }

aliased_names = {
    'homogeneous_solid'       : 'homogeneous_solid',
    'homogen_solid'           : 'homogeneous_solid',
    'liquid_solid'            : 'liquid_solid',
    'ls'                      : 'liquid_solid',
    'solid_liquid_solid'      : 'solid_liquid_solid',
    'sls'                     : 'solid_liquid_solid',
    'solid_solid_liquid_solid': 'solid_solid_liquid_solid',
    'ssls'                    : 'solid_solid_liquid_solid',
    'homogen-solid'           : 'homogeneous_solid',
    'liquid-solid'            : 'liquid_solid',
    'solid-liquid-solid'      : 'solid_liquid_solid',
    'solid-solid-liquid-solid': 'solid_solid_liquid_solid',
    }


def find_collapse_func(
    model: str,
    surface_boundary_condition: np.ndarray,
    is_liquid_static: bool = None,
    radius_array: np.ndarray = None,
    gravity_array: np.ndarray = None,
    density_array: np.ndarray = None,
    liquid_layer_indices: np.ndarray = None,
    frequency: float = None
    ) -> Tuple[callable, tuple]:
    # Clean up input
    model = model.lower()

    # Check for aliases
    model = aliased_names[model]

    if model in ['homogeneous_solid']:
        # No liquid model required. We can return right away.
        return known_collapse_models[model][True], (surface_boundary_condition,)
    elif model in ['liquid_solid']:
        if is_liquid_static is None:
            raise ValueError(
                'Collapse model contains a liquid layer but additional liquid layer parameters were'
                'not provided.'
                )
        elif is_liquid_static:
            # The static liquid case for a liquid-solid model requires no other input.
            return known_collapse_models[model][is_liquid_static], (surface_boundary_condition,)

    # There should be a liquid layer. Make sure all the information was provided.
    if is_liquid_static is None or radius_array is None or gravity_array is None or density_array is None or \
            liquid_layer_indices is None or frequency is None:
        raise ValueError(
            'Collapse model contains a liquid layer but additional liquid layer parameters were'
            'not provided.'
            )
    liquid_gravity = gravity_array[liquid_layer_indices]
    liquid_density = density_array[liquid_layer_indices]

    if is_liquid_static:
        # Static liquid layers only additionally require gravity and density.
        collapse_inputs = (surface_boundary_condition, liquid_gravity, liquid_density)
    else:

        liquid_radii = radius_array[liquid_layer_indices]
        collapse_inputs = (surface_boundary_condition, liquid_gravity, liquid_density, liquid_radii, frequency)

    return known_collapse_models[model][is_liquid_static], collapse_inputs
