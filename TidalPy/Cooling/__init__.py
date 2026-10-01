"""C++ cooling (heat-transport) models and their name-based factory."""

from TidalPy.Cooling.cooling import (
    CoolingResult,
    CoolingBase,
    OffCooling,
    ConvectiveCooling,
    ConductiveCooling,
    make_cooling,
    cooling_model_names,
    cooling_config_keys,
    cooling_off,
    convective,
    conductive,
)

__all__ = [
    "CoolingResult",
    "CoolingBase",
    "OffCooling",
    "ConvectiveCooling",
    "ConductiveCooling",
    "make_cooling",
    "cooling_model_names",
    "cooling_config_keys",
    # Direct functions; each accepts floats or ndarrays.
    "cooling_off",
    "convective",
    "conductive",
]
