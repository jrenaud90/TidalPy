# distutils: language = c++

from TidalPy.Stellar.luminosity cimport (
    LuminosityBase,
    FixedLuminosity,
    MassToLuminosity,
    PowerLawLuminosity,
    c_LuminosityBase,
    c_find_luminosity,
    c_luminosity_canonical_name,
    c_luminosity_model_names,
)
