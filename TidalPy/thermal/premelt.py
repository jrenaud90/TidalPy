from ..utilities.model import ModelSearcher
from . import viscosity_models

# Pre-melt Strength using the same parameters as rheology
from ..rheology.rheology import rheology_param_defaults

find_viscosity = ModelSearcher(viscosity_models, rheology_param_defaults)