"""Melting laws: melting curves (solidus and liquidus), melt weakening of the shear modulus and viscosity, and the
optional bulk-modulus and bulk-viscosity mixing of melt (``melting``). A material (``TidalPy.Material``) combines them
with its solid and liquid phases.
"""

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
