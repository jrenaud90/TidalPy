# distutils: language = c++

from TidalPy.Cooling.cooling cimport (
    CoolingResult,
    CoolingBase,
    c_CoolingBase,
    c_CoolingInputs,
    c_CoolingResult,
    c_find_cooling,
    c_cooling_canonical_name,
    c_cooling_model_names,
)
