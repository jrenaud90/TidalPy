"""C++ stellar luminosity models and their name-based factory."""

from TidalPy.stellar_x.luminosity import (
    LuminosityBase,
    FixedLuminosity,
    MassToLuminosity,
    PowerLawLuminosity,
    make_luminosity,
    fixed,
    mass_to_luminosity,
    power_law,
)

__all__ = [
    # Model classes
    "LuminosityBase",
    "FixedLuminosity",
    "MassToLuminosity",
    "PowerLawLuminosity",
    # Factory
    "make_luminosity",
    # Direct luminosity convenience functions (float or np.ndarray mass)
    "fixed",
    "mass_to_luminosity",
    "power_law",
]
