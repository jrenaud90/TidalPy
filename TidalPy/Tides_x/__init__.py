"""Global (1D) tide models, Love-number helpers, tidal potentials, and the point-wise 3D stress, strain,
heating, and displacement kernels.

The world's ``calc_tides`` and ``calc_3d_tides`` drive these internally; this package exposes the pieces for
direct use. Each name is also importable from its own subpackage.
"""
from TidalPy.Tides_x.classes import (
    TideBase,
    RheologyTide,
    FixedQTide,
    FixedLagTide,
    CTLQTide,
    make_tide,
    collapse_global_tides,
)
from TidalPy.Tides_x.love import (
    LoveNumbers,
    apply_fixed_dt,
    apply_fixed_q,
    calc_effective_rigidity,
    calc_homogeneous_love_numbers,
    love_method_name,
)
from TidalPy.Tides_x.potential import (
    ModeMap,
    UniqueFrequencyMap,
    tidal_potential_3d_modes,
    global_potential,
)
from TidalPy.Tides_x.multilayer import (
    angular_gram,
    displacement_point,
    strain_stress_heating_point,
    volumetric_heating,
)
from TidalPy.Tides_x.eccentricity import eccentricity_func
from TidalPy.Tides_x.obliquity import obliquity_func

__all__ = [
    # Global tide models
    "TideBase",
    "RheologyTide",
    "FixedQTide",
    "FixedLagTide",
    "CTLQTide",
    "make_tide",
    "collapse_global_tides",
    # Love numbers
    "LoveNumbers",
    "apply_fixed_dt",
    "apply_fixed_q",
    "calc_effective_rigidity",
    "calc_homogeneous_love_numbers",
    "love_method_name",
    # Tidal potentials
    "ModeMap",
    "UniqueFrequencyMap",
    "tidal_potential_3d_modes",
    "global_potential",
    # Point-wise 3D kernels
    "angular_gram",
    "displacement_point",
    "strain_stress_heating_point",
    "volumetric_heating",
    # Eccentricity and obliquity functions
    "eccentricity_func",
    "obliquity_func",
]
