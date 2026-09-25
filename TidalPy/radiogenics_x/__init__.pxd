# distutils: language = c++

from TidalPy.radiogenics_x.radiogenics cimport (
    RadiogenicsBase,
    OffRadiogenics,
    IsotopeRadiogenics,
    FixedRadiogenics,
    c_RadiogenicsBase,
    c_RadiogenicsConfig,
    c_RadiogenicsModel,
    c_Isotope,
    c_IsotopeDataset,
    c_OffRadiogenics,
    c_IsotopeRadiogenics,
    c_FixedRadiogenics,
    c_find_radiogenics,
    c_radiogenics_model_from_name,
    c_get_isotope_dataset,
    c_isotope_dataset_names,
)
