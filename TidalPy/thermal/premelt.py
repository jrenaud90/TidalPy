# Pre-melt Strength using the same parameters as rheology
from TidalPy.tides.tides import rheology_param_defaults
from . import viscosity_models
from ..utilities.model import ModelSearcher


find_viscosity = ModelSearcher(viscosity_models, rheology_param_defaults)
