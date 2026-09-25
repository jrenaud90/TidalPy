"""C++ cooling (heat-transport) models and their name-based factory."""

from TidalPy.cooling_x.cooling import (
    CoolingResult,
    CoolingBase,
    OffCooling,
    ConvectiveCooling,
    ConductiveCooling,
    make_cooling,
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
    # Direct functions; each accepts floats or ndarrays.
    "cooling_off",
    "convective",
    "conductive",
]
