# distutils: language = c++

from TidalPy.cooling_x.cooling cimport (
    CoolingResult,
    CoolingBase,
    OffCooling,
    ConvectiveCooling,
    ConductiveCooling,
    c_CoolingBase,
    c_CoolingConfig,
    c_CoolingInputs,
    c_CoolingResult,
    c_CoolingModel,
    c_OffCooling,
    c_ConvectiveCooling,
    c_ConductiveCooling,
    c_find_cooling,
    c_cooling_model_from_name,
)
