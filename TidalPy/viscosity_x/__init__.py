"""TidalPy viscosity_x — C++ viscosity model hierarchy.

Exposes the three viscosity models and a name-based factory:

- ``ArrheniusViscosity``  (alias ``"arr"``)   — Arrhenius flow law.
- ``ReferenceViscosity``  (alias ``"ref"``)   — relative-activation law.
- ``ConstantViscosity``   (alias ``"const"``) — temperature/pressure independent.

Each model returns the dynamic viscosity [Pa·s] as a function of temperature [K]
and pressure [Pa] via ``calc_viscosity``. This is the pre-melt ("solid")
viscosity that the partial-melt step weakens in the love-number pipeline.
"""

from TidalPy.viscosity_x.viscosity import (
    ViscosityBase,
    ConstantViscosity,
    ReferenceViscosity,
    ArrheniusViscosity,
    make_viscosity,
)

__all__ = [
    "ViscosityBase",
    "ConstantViscosity",
    "ReferenceViscosity",
    "ArrheniusViscosity",
    "make_viscosity",
]
