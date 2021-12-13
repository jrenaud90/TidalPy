from .lajolla import LinearSegmentedColormap, cm_data, lajolla_map

lajolla_map_r = LinearSegmentedColormap.from_list('lajolla_r', cm_data[::-1])
