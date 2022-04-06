""" TidalPy utilizes the open source scientific color maps created by Fabio Crameri and Grace Shephard.

If you use these color maps in your publications please consider citing their great work.
URL: https://zenodo.org/record/5501399
DOI: 10.5281/zenodo.5501399

"""

from .vik import vik_map, vik_map_r
from .lapaz import lapaz_map, lapaz_map_r
from .tofino import tofino_map, tofino_map_r
from .lajolla import lajolla_map, lajolla_map_r
from .cork import cork_map, cork_map_r
from .oleron import oleron_map_r, oleron_map
from .bukavu import bukavu_map, bukavu_map_r, bukavu_continuous_map, bukavu_continuous_map_r

KNOWN_CMAPS = {
    'vik': vik_map,
    'vik_r': vik_map_r,
    'lapaz': lapaz_map,
    'lapaz_r': lapaz_map_r,
    'tofino': tofino_map,
    'tofino_r': tofino_map_r,
    'lajolla': lajolla_map,
    'lajolla_r': lajolla_map_r,
    'cork': cork_map,
    'cork_r': cork_map_r,
    'oleron': oleron_map,
    'oleron_r': oleron_map_r,
    'bukavu': bukavu_map,
    'bukavu_r': bukavu_map_r,
    'bukavu_continuous': bukavu_continuous_map,
    'bukavu_continuous_r': bukavu_continuous_map_r,
    }