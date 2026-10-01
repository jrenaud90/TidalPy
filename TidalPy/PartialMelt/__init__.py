"""Melting laws: melting curves (solidus and liquidus), melt weakening of the shear modulus and viscosity, and the
optional bulk-modulus and bulk-viscosity mixing of melt (``melting``); and the partial-melt models that the world
pipeline still reads (``partial_melt``).
"""

from TidalPy.PartialMelt.partial_melt import (
    PartialMeltBase,
    OffPartialMelt,
    SpohnPartialMelt,
    HenningPartialMelt,
    make_partial_melt,
)
from TidalPy.PartialMelt.melting import (
    MeltingCurveBase,
    ConstantMeltingCurve,
    SimonGlatzelCurve,
    SimonGlatzel2Curve,
    InterpolatedMeltingCurve,
    make_melting_curve,
    melting_curve_model_names,
    MeltWeakeningBase,
    NoMeltWeakening,
    SpohnMeltWeakening,
    HenningMeltWeakening,
    make_melt_weakening,
    melt_weakening_model_names,
    BulkModulusMixingBase,
    HashinShtrikmanMixing,
    make_bulk_modulus_mixing,
    BulkViscosityMixingBase,
    CompactionViscosity,
    make_bulk_viscosity_mixing,
)

__all__ = [
    "PartialMeltBase",
    "OffPartialMelt",
    "SpohnPartialMelt",
    "HenningPartialMelt",
    "make_partial_melt",
    "MeltingCurveBase",
    "ConstantMeltingCurve",
    "SimonGlatzelCurve",
    "SimonGlatzel2Curve",
    "InterpolatedMeltingCurve",
    "make_melting_curve",
    "melting_curve_model_names",
    "MeltWeakeningBase",
    "NoMeltWeakening",
    "SpohnMeltWeakening",
    "HenningMeltWeakening",
    "make_melt_weakening",
    "melt_weakening_model_names",
    "BulkModulusMixingBase",
    "HashinShtrikmanMixing",
    "make_bulk_modulus_mixing",
    "BulkViscosityMixingBase",
    "CompactionViscosity",
    "make_bulk_viscosity_mixing",
]
