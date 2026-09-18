"""C++ cooling (heat-transport) models and their name-based factory.

``OffCooling`` (alias "none"), ``ConvectiveCooling`` (parameterized boundary-layer convection), and
``ConductiveCooling`` map a layer's thermal state to a ``CoolingResult`` (heat flux [W/m^2],
boundary-layer thickness [m], Rayleigh and Nusselt numbers) through ``calc_cooling``.
"""

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
