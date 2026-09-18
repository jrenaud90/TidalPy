"""C++ viscosity models and their name-based factory.

``ArrheniusViscosity`` (alias "arr"), ``ReferenceViscosity`` (alias "ref"), and
``ConstantViscosity`` (alias "const") return the dynamic viscosity [Pa s] at a temperature [K] and
pressure [Pa] through ``calc_viscosity``. This is the pre-melt (solid) viscosity that the
partial-melt step weakens in the Love number pipeline.
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
