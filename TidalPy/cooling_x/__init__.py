"""TidalPy cooling_x — C++ cooling (heat-transport) model hierarchy.

Exposes the three cooling models and a name-based factory:

- ``OffCooling``         (alias ``"none"``)  — cooling disabled (zero flux).
- ``ConvectiveCooling``  (alias ``"convective"``) — parameterized boundary-layer convection.
- ``ConductiveCooling``  (alias ``"conductive"``) — conduction across the layer.

Each model returns a ``CoolingResult`` (heat flux [W/m^2], boundary-layer
thickness [m], Rayleigh and Nusselt numbers) via ``calc_cooling``.
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
    # Result container
    "CoolingResult",
    # Model classes
    "CoolingBase",
    "OffCooling",
    "ConvectiveCooling",
    "ConductiveCooling",
    # Factory
    "make_cooling",
    # Direct convenience functions (float or np.ndarray inputs)
    "cooling_off",
    "convective",
    "conductive",
]
