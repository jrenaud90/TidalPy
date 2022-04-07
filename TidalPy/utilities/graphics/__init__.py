# Cartopy is usually not installed during testing which is fine. But this module will throw an error when tests are run
#    and Cartopy is not installed. So check if it is installed before bringing these packages up.
import importlib.util

from .grid_plot import GridPlot
from .planet_plot import geotherm_plot

if (spec := importlib.util.find_spec('cartopy')) is not None:
    from .global_map import projection_map