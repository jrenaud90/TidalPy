# distutils: language = c++

from TidalPy.stellar_x.luminosity cimport (
    LuminosityBase,
    FixedLuminosity,
    MassToLuminosity,
    PowerLawLuminosity,
    c_LuminosityBase,
    c_LuminosityConfig,
    c_LuminosityModel,
    c_FixedLuminosity,
    c_MassToLuminosity,
    c_PowerLawLuminosity,
    c_find_luminosity,
    c_luminosity_model_from_name,
)
