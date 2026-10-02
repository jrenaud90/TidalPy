"""C++ radiogenic heating models and their name-based factory."""

from TidalPy.Radiogenics.radiogenics import (
    RadiogenicsBase,
    OffRadiogenics,
    IsotopeRadiogenics,
    FixedRadiogenics,
    make_radiogenics,
    radiogenics_model_names,
    radiogenics_config_keys,
    available_isotope_datasets,
    isotope_dataset,
    isotope_dataset_parameters,
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
    "radiogenics_model_names",
    "radiogenics_config_keys",
    # Literature isotope datasets
    "available_isotope_datasets",
    "isotope_dataset",
    "isotope_dataset_parameters",
    # Direct heating convenience functions (float or np.ndarray inputs)
    "off",
    "isotope",
    "fixed",
]
