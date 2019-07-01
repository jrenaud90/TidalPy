from . import viscosity_models
# Pre-melt Strength using the same parameters as rheology
from ..rheology.rheology import rheology_param_defaults
from ..utilities.model import ModelSearcher


find_viscosity = ModelSearcher(viscosity_models, rheology_param_defaults)
