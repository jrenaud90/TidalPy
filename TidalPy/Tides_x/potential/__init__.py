import importlib as _importlib

from TidalPy.Tides_x.potential.potential_common import ModeMap, UniqueFrequencyMap, test_mode_map
from TidalPy.Tides_x.potential.tidal_potential import (
    TidalPotentialBase, SyncLowEPotential, NSRModesPotential, make_tidal_potential)

# `global` is a Python reserved keyword, so we use importlib to import it.
_global_mod = _importlib.import_module('TidalPy.Tides_x.potential.global')
global_potential = _global_mod.global_potential
