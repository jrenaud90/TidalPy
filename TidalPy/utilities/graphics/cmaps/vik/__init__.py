from .vik import LinearSegmentedColormap, cm_data, vik_map

vik_map_r = LinearSegmentedColormap.from_list('vik_r', cm_data[::-1])