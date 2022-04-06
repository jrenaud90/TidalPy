from .tofino import LinearSegmentedColormap, cm_data, tofino_map

tofino_map_r = LinearSegmentedColormap.from_list('tofino_r', cm_data[::-1])