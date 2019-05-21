from ..utilities.search import ModelSearcher
from . import viscosity

# Pre-melt Strength using the same parameters as rheology
from ..rheology.rheology import rheology_param_defaults

find_viscosity = ModelSearcher(viscosity, rheology_param_defaults)