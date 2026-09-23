"""C++ radiogenic heating models and their name-based factory."""

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
