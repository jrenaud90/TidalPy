from .cork import LinearSegmentedColormap, cm_data, cork_map

cork_map_r = LinearSegmentedColormap.from_list('cork_r', cm_data[::-1])
