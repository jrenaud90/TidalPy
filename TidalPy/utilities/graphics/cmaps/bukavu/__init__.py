import numpy as np

from .bukavu import LinearSegmentedColormap, cm_data, bukavu_map

bukavu_map_r = LinearSegmentedColormap.from_list('bukavu_r', cm_data[::-1])


_bukavu_bottom = bukavu_map(np.linspace(0., 0.5, 128)[::-1])
_bukavu_top = bukavu_map(np.linspace(0.5, 1, 128))

# combine them and build a new colormap
_colors = np.vstack((_bukavu_bottom, _bukavu_top))
bukavu_continuous_map = LinearSegmentedColormap.from_list('bukavu_cont', _colors)
bukavu_continuous_map_r = LinearSegmentedColormap.from_list('bukavu_cont_r', _colors[::-1])