"""TidalPy radiogenics_x — C++ radiogenics model hierarchy.

Exposes the three radiogenic-heating models and a name-based factory:

- ``OffRadiogenics``     (alias ``"none"``)     — radiogenics disabled.
- ``IsotopeRadiogenics``                        — sum of decaying isotopes.
- ``FixedRadiogenics``   (alias ``"constant"``) — lumped rate with optional decay.

Each model computes the total radiogenic heating [W] via ``calc_heating``.
"""

from TidalPy.radiogenics_x.radiogenics import (
    RadiogenicsBase,
    OffRadiogenics,
    IsotopeRadiogenics,
    FixedRadiogenics,
    make_radiogenics,
    available_isotope_datasets,
    isotope_dataset,
    off,
    isotope,
    fixed,
)

__all__ = [
    # Model classes
    "RadiogenicsBase",
    "OffRadiogenics",
    "IsotopeRadiogenics",
    "FixedRadiogenics",
    # Factory
    "make_radiogenics",
    # Built-in literature isotope datasets
    "available_isotope_datasets",
    "isotope_dataset",
    # Direct heating convenience functions (float or np.ndarray inputs)
    "off",
    "isotope",
    "fixed",
]
