"""TidalPy partial_melt_x — C++ partial-melt (melt-weakening) model hierarchy.

Exposes the three partial-melt models and a name-based factory:

- ``OffPartialMelt``     (alias ``"none"``)    — no melt weakening (returns pre-melt).
- ``SpohnPartialMelt``   (alias ``"fischer"``) — Fischer & Spohn (1990) temperature law.
- ``HenningPartialMelt``                       — Henning (2009/2010) three-regime law.

Each model maps a material's pre-melt viscosity + shear modulus (and temperature)
to the post-melt viscosity + shear modulus, and reports the volumetric melt
fraction. These quantities are frequency-independent and feed the downstream
rheology (complex modulus) step in the whole-planet love-number pipeline.
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
