# distutils: language = c++

from TidalPy.Radiogenics.radiogenics cimport (
    RadiogenicsBase,
    OffRadiogenics,
    IsotopeRadiogenics,
    FixedRadiogenics,
    c_RadiogenicsBase,
    c_Isotope,
    c_IsotopeDataset,
    c_IsotopeRadiogenics,
    c_find_radiogenics,
    c_make_isotope_radiogenics,
    c_radiogenics_canonical_name,
    c_radiogenics_model_names,
    c_get_isotope_dataset,
    c_isotope_dataset_names,
)
