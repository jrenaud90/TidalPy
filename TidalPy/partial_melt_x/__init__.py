"""C++ partial-melt (melt-weakening) models and their name-based factory.

Each maps a material's pre-melt viscosity and shear modulus, plus its temperature, to the post-melt values.
Frequency independent, so these feed the downstream rheology (complex modulus) step.
"""

from TidalPy.partial_melt_x.partial_melt import (
    PartialMeltBase,
    OffPartialMelt,
    SpohnPartialMelt,
    HenningPartialMelt,
    make_partial_melt,
)

__all__ = [
    "PartialMeltBase",
    "OffPartialMelt",
    "SpohnPartialMelt",
    "HenningPartialMelt",
    "make_partial_melt",
]
