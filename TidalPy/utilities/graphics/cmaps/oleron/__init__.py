from .oleron import LinearSegmentedColormap, cm_data, oleron_map

oleron_map_r = LinearSegmentedColormap.from_list('oleron_r', cm_data[::-1])
