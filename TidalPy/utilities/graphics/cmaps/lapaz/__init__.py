from .lapaz import LinearSegmentedColormap, cm_data, lapaz_map

lapaz_map_r = LinearSegmentedColormap.from_list('lapaz_r', cm_data[::-1])
