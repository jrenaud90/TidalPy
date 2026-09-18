"""TidalPy partial_melt_x: the C++ partial-melt (melt-weakening) model hierarchy.

Each model maps a material's pre-melt viscosity and shear modulus, plus its temperature, to the post-melt
values and reports the volumetric melt fraction. These are frequency-independent and feed the downstream
rheology (complex modulus) step. Models: ``OffPartialMelt`` (alias ``"none"``), ``SpohnPartialMelt``
(alias ``"fischer"``), and ``HenningPartialMelt``, plus the ``make_partial_melt`` factory.
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
