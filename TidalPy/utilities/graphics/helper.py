""" Various helper functions used by other graphics functions and modules. """
import os.path
from typing import Union

import matplotlib.pyplot as plt
from matplotlib.colors import Colormap

from .cmaps import KNOWN_CMAPS


def get_cmap(cmap: Union[str, Colormap]):
    """ Looks at TidalPy and Matplotlib's colormap cataloger for requested cmap name.

    Parameters
    ----------
    cmap : Union[str, Colormap]

    Returns
    -------
    cmap : Union[str, Colormap]

    """

    # Find color map
    if type(cmap) is str:
        # Name of a cmap is provided. Look at TidalPy's or Matplotlib's catalog to see if we know it.
        if cmap in KNOWN_CMAPS:
            cmap = KNOWN_CMAPS[cmap]
        elif cmap.lower() in KNOWN_CMAPS:
            cmap = KNOWN_CMAPS[cmap]
        elif cmap in plt.colormaps():
            pass
        elif cmap.lower() in plt.colormaps():
            pass
        else:
            raise KeyError(f'Unknown color map name provided: {cmap}.')
    else:
        # Assume that the user is attempting to pass a properly formatted color map.
        # Check that it is the correct instance.
        if not isinstance(cmap, Colormap):
            raise TypeError('Unexpected type found for user provided cmap.')

    return cmap
